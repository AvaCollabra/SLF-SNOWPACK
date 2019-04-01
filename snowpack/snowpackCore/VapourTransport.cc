/*
 *  SNOWPACK stand-alone
 *
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
*/
/*  This file is part of Snowpack.
	Snowpack is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	Snowpack is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with Snowpack.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <snowpack/snowpackCore/VapourTransport.h>

#include <assert.h>
#include <sstream>
#include <errno.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/
/// @brief Define the assembly macro
void VapourTransport::EL_INCID(const int &e, int Ie[]) {
	Ie[0] = e;
	Ie[1] = e+1;
}

/// @brief Define the node to element temperature macro
void VapourTransport::EL_TEMP( const int Ie[], double Te0[], double Tei[], const std::vector<NodeData> &T0, const double Ti[] ) {
	Te0[ 0 ] = T0[ Ie[ 0 ] ].rhov;
	Te0[ 1 ] = T0[ Ie[ 1 ] ].rhov;
	Tei[ 0 ] = Ti[ Ie[ 0 ] ];
	Tei[ 1 ] = Ti[ Ie[ 1 ] ];
}

/// @brief Element right-hand side macro
void VapourTransport::EL_RGT_ASSEM(double F[], const int Ie[], const double Fe[]) {
	F[Ie[0]] += Fe[0];
	F[Ie[1]] += Fe[1];
}

/************************************************************
 * non-static section                                       *
 ************************************************************/
VapourTransport::VapourTransport(const SnowpackConfig& cfg)
               : WaterTransport(cfg), RichardsEquationSolver1d(cfg, false), variant(),
                 iwatertransportmodel_snow(BUCKET), iwatertransportmodel_soil(BUCKET), watertransportmodel_snow("BUCKET"), watertransportmodel_soil("BUCKET"),
                 sn_dt(IOUtils::nodata),
                 hoar_thresh_rh(IOUtils::nodata), hoar_thresh_vw(IOUtils::nodata), hoar_thresh_ta(IOUtils::nodata),
                 /*hoar_density_buried(IOUtils::nodata), hoar_density_surf(IOUtils::nodata), hoar_min_size_buried(IOUtils::nodata),
                 minimum_l_element(IOUtils::nodata),*/ useSoilLayers(false), water_layer(false), enable_vapour_transport(false)
{
	cfg.getValue("VARIANT", "SnowpackAdvanced", variant);

	// Defines whether soil layers are used
	cfg.getValue("SNP_SOIL", "Snowpack", useSoilLayers);

	//To build a thin top rain-water layer over a thin top ice layer, rocks, roads etc.
	cfg.getValue("WATER_LAYER", "SnowpackAdvanced", water_layer);

	/**
	 * @brief No surface hoar will form for rH above threshold (1)
	 * - Original calibration with the 98/99 data set: 0.9
	 * - r141: HOAR_THRESH_RH set to 0.9
	 * - r719: HOAR_THRESH_RH set to 0.97
	 */
	cfg.getValue("HOAR_THRESH_RH", "SnowpackAdvanced", hoar_thresh_rh);

	/**
	 * @brief No surface hoar will form at wind speeds above threshold (m s-1)
	 * - Original calibration with the 98/99 data set: 3.5
	 * - r141: HOAR_THRESH_VW set to 3.0
	 * - r242: HOAR_THRESH_VW set to 3.5
	 */
	cfg.getValue("HOAR_THRESH_VW", "SnowpackAdvanced", hoar_thresh_vw);

	/**
	 * @brief No surface hoar will form at air temperatures above threshold (m s-1)
	 * - Originaly, using THRESH_RAIN
	 * - r787: HOAR_THRESH_TA set to 1.2
	 */
	cfg.getValue("HOAR_THRESH_TA", "SnowpackAdvanced", hoar_thresh_ta);

	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	const double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	sn_dt = M_TO_S(calculation_step_length);

	//Enable vapour transport
	cfg.getValue("ENABLE_VAPOUR_TRANSPORT", "SnowpackAdvanced", enable_vapour_transport);

	//Water transport model snow
	cfg.getValue("WATERTRANSPORTMODEL_SNOW", "SnowpackAdvanced", watertransportmodel_snow);
	iwatertransportmodel_snow=UNDEFINED;
	if (watertransportmodel_snow=="BUCKET") {
		iwatertransportmodel_snow=BUCKET;
	} else if (watertransportmodel_snow=="NIED") {
		iwatertransportmodel_snow=NIED;
	} else if (watertransportmodel_snow=="RICHARDSEQUATION") {
		iwatertransportmodel_snow=RICHARDSEQUATION;
	}

	//Water transport model soil
	cfg.getValue("WATERTRANSPORTMODEL_SOIL", "SnowpackAdvanced", watertransportmodel_soil);
	iwatertransportmodel_soil=UNDEFINED;
	if (watertransportmodel_soil=="BUCKET") {
		iwatertransportmodel_soil=BUCKET;
	} else if (watertransportmodel_soil=="NIED") {
		iwatertransportmodel_soil=NIED;
	} else if (watertransportmodel_soil=="RICHARDSEQUATION") {
		iwatertransportmodel_soil=RICHARDSEQUATION;
	}
}

void VapourTransport::compTransportMass(const CurrentMeteo& Mdata, double& ql,
                                       SnowStation& Xdata, SurfaceFluxes& Sdata, const double& surfaceVaporPressure)
{

	// First, consider no soil with no snow on the ground
	if (!useSoilLayers && Xdata.getNumberOfNodes() == Xdata.SoilNode+1) {
		return;
	}

    try {
	    LayerToLayer(Mdata, Xdata, Sdata, ql, surfaceVaporPressure);
	    WaterTransport::adjustDensity(Xdata);
	    //WaterTransport::mergingElements(Xdata, Sdata);
    } catch(const exception&)
    {
	    prn_msg( __FILE__, __LINE__, "err", Mdata.date, "Error in transportVapourMass()");
	    throw;
    }


}

void VapourTransport::LayerToLayer(const CurrentMeteo& Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata, double& ql, const double& surfaceVaporPressure)
 {
	const size_t nN = Xdata.getNumberOfNodes();
	size_t nE = nN-1;
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;
	size_t e = nE;
	std::vector<double> deltaM(nE, 0.);//Calculate the limited layer mass change
	std::vector<double> totalMassChange(nE, 0.);// store the total mass change
	std::vector<double> vaporFluxDiv(nE, 0.);// store the vapor flux divergence

	std::vector<double> factor_(nE, 1.);// this is for source term in vapor transport equation
	for(size_t i=0; i<Xdata.SoilNode; i++)
	{
		factor_[i]=0.;
	}

	std::vector<double> D(nE, 0.);
	for (size_t i=0; i<=nE-1; i++)
	{
		double aaa = EMS[i].theta[AIR];
		double nnn = 1.- EMS[i].theta[SOIL];
		double D_vapSoil = Constants::diffusion_coefficient_in_air * pow(aaa,10./3.)/nnn/nnn; // based on jury1983
		D[i] = factor_[i]*Constants::diffusion_coefficient_in_snow + (1.0-factor_[i])*D_vapSoil;
	}

	// For snow, update theta_s (i.e., pore space)
	if (nE > Xdata.SoilNode) {
		for (size_t el = Xdata.SoilNode; el < nE; el++) {
			EMS[el].VG.theta_s = (1. - EMS[el].theta[ICE])*(Constants::density_ice/Constants::density_water);	// TODO: link to van Genuchten parameterisation
		}
	}

	double sumMassChange = 0.0;
	double sumMassChange2 = 0.0;
	double sumMassChange3=0.0;
	double inletVapourMass = std::abs(ql)/Constants::lh_sublimation*sn_dt;
	//double ql_bcup = ql;
	bool forcingMassBalance= false;

	if (enable_vapour_transport) {

        // compute vapor density profile
        compDensityProfile(Mdata, Xdata, true, ql, surfaceVaporPressure);
       	ql = 0.; //Now that we used the remaining ql, put it to 0.

		//calculation of vapor diffusion flux
		for (size_t i = 0; i < nE; i++) {
			//double aaa = EMS[i].theta[AIR];
			//double nnn = 1.- EMS[i].theta[SOIL];
			//double D_vapSoil = Constants::diffusion_coefficient_in_air * pow(aaa,10./3.)/nnn/nnn; // Jafari added, based on jury1983
			//double D = factor_[i]*Constants::diffusion_coefficient_in_snow + (1.0-factor_[i])*D_vapSoil;
			//double topSaturatedVapor=Atmosphere::waterVaporDensity(NDS[i+1].T, Atmosphere::vaporSaturationPressure(NDS[i+1].T));
			//double botSaturatedVapor=Atmosphere::waterVaporDensity(NDS[i].T, Atmosphere::vaporSaturationPressure(NDS[i].T));
			//EMS[i].vapTrans_fluxDiff =-D*(topSaturatedVapor-botSaturatedVapor)/EMS[i].L;
			EMS[i].vapTrans_fluxDiff =-D[i]*(NDS[i+1].rhov-NDS[i].rhov)/EMS[i].L;
		}

		// scaling to force the mass to be conserved
		if(forcingMassBalance)
		{
			for (size_t i=0; i<=nE-1; i++)
			{
				double saturationVapor = Atmosphere::waterVaporDensity(EMS[i].Te, Atmosphere::vaporSaturationPressure(EMS[i].Te));
				totalMassChange[i] =(EMS[i].rhov -saturationVapor)*EMS[i].theta[AIR]*EMS[i].L; //total mass change, (kg m-2 )
				sumMassChange=sumMassChange+std::abs(totalMassChange[i]);
			}
			for (size_t i=0; i<=nE-1; i++)
			{
				totalMassChange[i]=(totalMassChange[i]/sumMassChange)*inletVapourMass;
				sumMassChange2=sumMassChange2+totalMassChange[i];
			}
		}else
		{
			for (size_t i=0; i<=nE-1; i++)
			{
				double saturationVapor = Atmosphere::waterVaporDensity(EMS[i].Te, Atmosphere::vaporSaturationPressure(EMS[i].Te));
				totalMassChange[i] =(EMS[i].rhov -saturationVapor)*EMS[i].theta[AIR]*EMS[i].L; //total mass change, (kg m-2 )
			}
		}

		// consider the mass change due to vapour transport in snow/soil
		while (e-- > 0) {
			const double massPhaseChange = totalMassChange[e];

			double dM = 0.;	//mass change induced by vapor flux (kg m-2)

			// Now, the mass change is limited by:
			// - we cannot remove more WATER and ICE than available
			// - we cannot add more WATER and ICE than pore space available
			if ( EMS[e].theta[SOIL] < Constants::eps ) {// there is no soil in element to keep element not to merge
			    dM = std::max(  -((EMS[e].theta[WATER] - EMS[e].VG.theta_r * (1. + Constants::eps)) * Constants::density_water * EMS[e].L + (EMS[e].theta[ICE] - Snowpack::min_ice_content) * Constants::density_ice * EMS[e].L)  ,
				                std::min(  (EMS[e].theta[AIR] * Constants::density_ice * EMS[e].L), massPhaseChange  )
				     ); // mass change due to difference in water vapor flux (kg m-2), at most can fill the pore space.
			} else {

				dM = std::max(  -((EMS[e].theta[WATER] - EMS[e].VG.theta_r * (1. + Constants::eps)) * Constants::density_water * EMS[e].L + EMS[e].theta[ICE] * Constants::density_ice * EMS[e].L)  ,
				                std::min(  (EMS[e].theta[AIR] * Constants::density_ice * EMS[e].L), massPhaseChange  )
				     ); // mass change due to difference in water vapor flux (kg m-2), at most can fill the pore space.

			}

			// If there is no pore space, or, in fact, only so much pore space to accomodate the larger volume occupied by ice when all water freezes,
			// we inhibit vapour flux. This is necessary to maintain saturated conditions when present, and this is in turn necessary for the stability in the Richards equation solver.
			if(EMS[e].theta[AIR] < EMS[e].theta[WATER]*(Constants::density_water/Constants::density_ice - 1.) + Constants::eps) {
				dM = 0.;
			}

			deltaM[e] += dM;
		}
	} else {
		compSurfaceSublimation(Mdata, ql, Xdata, Sdata);
		// Only deal with the remaining ql (i.e., latent heat exchange at the surface)
		/*
		const double topFlux = -ql / Constants::lh_sublimation;										//top layer flux (kg m-2 s-1)
		deltaM[nE-1] += std::max(-EMS[nE-1].theta[ICE] * (Constants::density_ice * EMS[nE-1].L), -(topFlux * sn_dt));
		// HACK: note that if we cannot satisfy the ql at this point, we overestimated the latent heat from soil.
		// We will not get mass from deeper layers, as to do that, one should work with enable_vapour_transport == true.
		Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] -= topFlux * sn_dt;
		*/
	}

	double dHoar = 0.;

	for (size_t i=0; i<=nE-1; i++)
	{
		sumMassChange3=sumMassChange3+deltaM[i];
	}

	for (e = 0; e < nE; e++) {
		EMS[e].M += deltaM[e];
		assert(EMS[e].M >= (-Constants::eps2)); //mass must be positive

		EMS[e].Qmm = 0.0;

		if (deltaM[e] < 0.) {
			// Mass loss: apply mass change first to water, then to ice, based on energy considerations
			// We can only do this partitioning here in this "simple" way, without checking if the mass is available, because we already limited dM above, based on available ICE + WATER.
			const double dTh_water = std::max( (EMS[e].VG.theta_r * (1. + Constants::eps) - EMS[e].theta[WATER])  ,  deltaM[e] / (Constants::density_water * EMS[e].L) );
			const double dTh_ice = ( deltaM[e] - (dTh_water * Constants::density_water * EMS[e].L) ) / (Constants::density_ice * EMS[e].L);
			EMS[e].theta[WATER] += dTh_water;
			EMS[e].theta[ICE] += dTh_ice;

			Sdata.mass[SurfaceFluxes::MS_EVAPORATION] += dTh_water * Constants::density_water * EMS[e].L;
			Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += dTh_ice * Constants::density_ice * EMS[e].L;
			EMS[e].M += dTh_water * Constants::density_water * EMS[e].L+dTh_ice * Constants::density_ice * EMS[e].L;
			assert(EMS[e].M >= (-Constants::eps2)); //mass must be positive

			EMS[e].Qmm += (dTh_water*Constants::density_water*Constants::lh_vaporization +
							dTh_ice*Constants::density_ice*Constants::lh_sublimation)/sn_dt;//[w/m^3]

			// If present at surface, surface hoar is sublimated away
			if (e == nE-1 && deltaM[e]<0) {
				dHoar = std::max(-NDS[nN-1].hoar, deltaM[e]);
			}
		} else {		// Mass gain: add water in case temperature at or above melting point, ice otherwise
			if (EMS[e].Te >= EMS[e].meltfreeze_tk) {
				EMS[e].theta[WATER] += deltaM[e] / (Constants::density_water * EMS[e].L);
				EMS[e].Qmm += (deltaM[e]*Constants::lh_vaporization)/sn_dt/EMS[e].L;//  [w/m^3]
				Sdata.mass[SurfaceFluxes::MS_EVAPORATION] += deltaM[e]; //
			} else {
				EMS[e].theta[ICE] += deltaM[e] / (Constants::density_ice * EMS[e].L);
				EMS[e].Qmm += (deltaM[e]*Constants::lh_sublimation)/sn_dt/EMS[e].L;// [w/m^3]
				Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += deltaM[e]; //
			}
			EMS[e].M += deltaM[e];
		}

    	EMS[e].theta[AIR] = std::max(1. - EMS[e].theta[WATER] - EMS[e].theta[WATER_PREF] - EMS[e].theta[ICE] - EMS[e].theta[SOIL],0.);

	    if(std::fabs(EMS[e].theta[AIR])<1.e-15)
    	{
       		EMS[e].theta[AIR]=0;
    	}

		EMS[e].Rho = (EMS[e].theta[ICE] * Constants::density_ice)
			      + ((EMS[e].theta[WATER] + EMS[e].theta[WATER_PREF]) * Constants::density_water)
			      + (EMS[e].theta[SOIL] * EMS[e].soil[SOIL_RHO]);
		assert(EMS[e].Rho > 0 || EMS[e].Rho==IOUtils::nodata); //density must be positive

		if (!(EMS[e].Rho > Constants::eps && EMS[e].theta[AIR] >= 0.)) {
				prn_msg(__FILE__, __LINE__, "err", Date(),
				    "Volume contents: e=%d nE=%d rho=%lf ice=%lf wat=%lf wat_pref=%lf soil=%lf air=%le", e, nE, EMS[e].Rho, EMS[e].theta[ICE],
                        EMS[e].theta[WATER], EMS[e].theta[WATER_PREF], EMS[e].theta[SOIL], EMS[e].theta[AIR]);
				throw IOException("Cannot evaluate mass balance in vapour transport LayerToLayer routine", AT);
		}

		//some useful output in case of vapor transport
		EMS[e].vapTrans_snowDenChangeRate = deltaM[e]/sn_dt/EMS[e].L;
		EMS[e].vapTrans_cumulativeDenChange += deltaM[e]/EMS[e].L;
	}

	Sdata.hoar += dHoar;
	NDS[nN-1].hoar += dHoar;
	if (NDS[nN-1].hoar < 0.) {
		NDS[nN-1].hoar = 0.;
	}
}

void VapourTransport::compSurfaceSublimation(const CurrentMeteo& Mdata, double& ql, SnowStation& Xdata, SurfaceFluxes& Sdata)
{
	double dL = 0., dM = 0.;     // Length and mass changes
	double M = 0.;               // Initial mass and volumetric content (water or ice)
	double dHoar = 0.;           // Actual change in hoar mass
	double cH_old;               // Temporary variable to hold height of snow

	const size_t nN = Xdata.getNumberOfNodes();
	size_t nE = nN-1;
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;
	const double Tss = NDS[nE].T; // Surface Temperature

	/*
	 * If ql > 0:
	 * Surface hoar is formed when surface temperature is below freezing.
	 * If no surface hoar can be formed, ql is kept and is used as boundary condition
	 * when calculating vapour flux.
	 * If there are elements and ql < 0:
	 * If ql is large enough to remove full surface elements, remove them.
	 * left over ql is used as boundary condition when calculating vapour flux.
	 *
	 * In both cases: add/subtract mass to MS_SUBLIMATION
	 */
	if (ql > Constants::eps2) { // Add Mass
		const double meltfreeze_tk = (Xdata.getNumberOfElements()>0)? Xdata.Edata[Xdata.getNumberOfElements()-1].meltfreeze_tk : Constants::meltfreeze_tk;
		if (Tss < meltfreeze_tk) { // Add Ice
			dM = ql*sn_dt/Constants::lh_sublimation;
			//if rh is very close to 1, vw too high or ta too high, surface hoar is destroyed and should not be formed
			if (!((Mdata.rh > hoar_thresh_rh) || (Mdata.vw > hoar_thresh_vw) || (Mdata.ta >= IOUtils::C_TO_K(hoar_thresh_ta)))) {
				// Under these conditions, form surface hoar
				ql = 0.;
				Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += dM;
				dHoar = dM;

				// In this case adjust properties of element, keeping snow density constant
				const double L_top = EMS[nE-1].L;
				const double theta_i0 = EMS[nE-1].theta[ICE];
				dL = dM/(EMS[nE-1].Rho); // length change
				if (nE == Xdata.SoilNode) {
					dL = 0.;
					dM = std::min(dM,EMS[nE-1].theta[AIR]*(Constants::density_ice*EMS[nE-1].L));
				}
				NDS[nE].z += dL + NDS[nE].u; NDS[nE].u = 0.0;
				EMS[nE-1].L0 = EMS[nE-1].L = L_top + dL;
				EMS[nE-1].E = EMS[nE-1].Eps = EMS[nE-1].dEps = EMS[nE-1].Eps_e = EMS[nE-1].Eps_v = EMS[nE-1].S = 0.0;
				EMS[nE-1].theta[ICE] *= L_top/EMS[nE-1].L;
				EMS[nE-1].theta[ICE] += dM/(Constants::density_ice*EMS[nE-1].L);
				EMS[nE-1].theta[ICE] = std::max(0., std::min(1., EMS[nE-1].theta[ICE]));
				EMS[nE-1].theta[WATER] *= L_top/EMS[nE-1].L;
				EMS[nE-1].theta[WATER] = std::max(0., std::min(1., EMS[nE-1].theta[WATER]));
				EMS[nE-1].theta[WATER_PREF] *= L_top/EMS[nE-1].L;
				EMS[nE-1].theta[WATER_PREF] = std::max(0., std::min(1., EMS[nE-1].theta[WATER_PREF]));

				for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[nE-1].conc[ICE][ii] *= L_top*theta_i0/(EMS[nE-1].theta[ICE]*EMS[nE-1].L);
				}

				EMS[nE-1].M += dM;
				assert(EMS[nE-1].M >= (-Constants::eps2)); //mass must be positive

				// Update remaining volumetric contents and density
				EMS[nE-1].theta[AIR] = std::max(0., 1.0 - EMS[nE-1].theta[WATER] - EMS[nE-1].theta[WATER_PREF] - EMS[nE-1].theta[ICE] - EMS[nE-1].theta[SOIL]);
				EMS[nE-1].updDensity();
			}
		}
	} else if ((ql < (-Constants::eps2)) && (nE > 0)) {
		// If ql < 0, SUBLIMATE mass off
		std::vector<double> M_Solutes(Xdata.number_of_solutes, 0.); // Mass of solutes from disappearing phases
		size_t e = nE;
		while ((e > 0) && (ql < (-Constants::eps2))) {  // While energy is available
			e--;
			/*
			* Determine the amount of potential sublimation and collect some variables
			* that will be continuously used: L0 and M
			*/
			const double L0 = EMS[e].L;
			const double theta_i0 = EMS[e].theta[ICE];
			M = theta_i0*Constants::density_ice*L0;
			dM = ql*sn_dt/Constants::lh_sublimation;
			if (-dM > M) {
				// Only if mass change is sufficient to remove the full element
				dM = -M;
				// Add solutes to Storage
				for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					M_Solutes[ii] += EMS[e].conc[ICE][ii]*theta_i0*L0;
				}
				EMS[e].theta[ICE] = 0.;
				dL = 0.;

				EMS[e].M += dM;
				Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += dM;
				ql -= dM*Constants::lh_sublimation/sn_dt;     // Update the energy used

				// If present at surface, surface hoar is sublimated away
				if (e == nE-1) {
					dHoar = std::max(-NDS[nN-1].hoar, dM);
				}

				// Update remaining volumetric contents and density
				EMS[e].theta[AIR] = std::max(0., 1.0 - EMS[e].theta[WATER] - EMS[e].theta[WATER_PREF] - EMS[e].theta[ICE] - EMS[e].theta[SOIL]);
				EMS[e].updDensity();
				// Merge the element if it is a snow layer. This will take care of possible left over liquid water (will be put one layer down)
				// Keep layer if it is a soil layer inside the snowpack (for example with snow farming)
				if(e>=Xdata.SoilNode) {
					if(EMS[e].theta[SOIL]<Constants::eps) {
						if (e>0) SnowStation::mergeElements(EMS[e-1], EMS[e], false, true);
						// Now reduce the number of elements by one.
						nE--;
					}
					//In case e==Xdata.SoilNode, we removed the last snow element and we should break out of the loop.
					if(e==Xdata.SoilNode) break;
				}
			} else {
				// Not enough energy anymore to remove complete element, so we should break out of the loop.
				break;
			}

			//check that thetas and densities are consistent
			assert(EMS[e].theta[SOIL] >= (-Constants::eps2) && EMS[e].theta[SOIL] <= (1.+Constants::eps2));
			assert(EMS[e].theta[ICE] >= (-Constants::eps2) && EMS[e].theta[ICE]<=(1.+Constants::eps2));
			assert(EMS[e].theta[WATER] >= (-Constants::eps2) && EMS[e].theta[WATER]<=(1.+Constants::eps2));
			assert(EMS[e].theta[WATER_PREF] >= (-Constants::eps2) && EMS[e].theta[WATER_PREF]<=(1.+Constants::eps2));
			assert(EMS[e].theta[AIR] >= (-Constants::eps2) && EMS[e].theta[AIR]<=(1.+Constants::eps2));
			assert(EMS[e].Rho >= (-Constants::eps2) || EMS[e].Rho==IOUtils::nodata); //we want positive density
		}

		// Now take care of left over solute mass.
		if (nE == Xdata.SoilNode) { // Add Solute Mass to Runoff TODO HACK CHECK
			for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
				Sdata.load[ii] += M_Solutes[ii]/S_TO_H(sn_dt);
			}
		} else { // Add Solute Mass to Element below
			if (EMS[e].theta[WATER] > 0.) {
				for(size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[e].conc[WATER][ii] += M_Solutes[ii]/EMS[e].theta[WATER]/EMS[e].L;
				}
			} else if (EMS[e].theta[ICE] > 0.) {
				for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[e].conc[ICE][ii] += M_Solutes[ii]/EMS[e].theta[ICE]/EMS[e].L;
				}
			} else {
				for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[e].conc[SOIL][ii] += M_Solutes[ii]/EMS[e].theta[SOIL]/EMS[e].L;
				}
			}
		}
		Xdata.reduceNumberOfElements(nE);
	}

	// HACK: this code is under verification. The comment reads "surface hoar *is* destroyed, but the next line says surface hoar *may be* destroyed, depending on the sign of the latent heat flux.
	// If the code is correct, we can delete this part, if the comment is correct, we should modify the code to read: hoar = -NDS[nE].hoar;
	// Check for surface hoar destruction or formation (once upon a time ml_sn_SurfaceHoar)
	/*if ((Mdata.rh > hoar_thresh_rh) || (Mdata.vw > hoar_thresh_vw) || (Mdata.ta >= IOUtils::C_TO_K(hoar_thresh_ta))) {
		//if rh is very close to 1, vw too high or ta too high, surface hoar is destroyed
		hoar = std::min(hoar, 0.);
	}*/

	Sdata.hoar += dHoar;
	NDS[nN-1].hoar += dHoar;
	if (NDS[nN-1].hoar < 0.) {
		NDS[nN-1].hoar = 0.;
	}

	// Surface hoar cannot exist when the top element is wet
	if (nE > 0) {
		const double theta_r=((iwatertransportmodel_snow==RICHARDSEQUATION && nE-1>=Xdata.SoilNode) || (iwatertransportmodel_soil==RICHARDSEQUATION && nE-1<Xdata.SoilNode)) ? (PhaseChange::RE_theta_r) : (PhaseChange::theta_r);
		if (Xdata.Edata[nE-1].theta[WATER] > theta_r) {
			NDS[nE].hoar = 0.;
		}
	}

	// At the end also update the overall height
	cH_old = Xdata.cH;
	Xdata.cH = NDS[Xdata.getNumberOfNodes()-1].z + NDS[Xdata.getNumberOfNodes()-1].u;
	if (Xdata.mH!=Constants::undefined) Xdata.mH -= std::min(Xdata.mH - Xdata.Ground, (cH_old - Xdata.cH));	// TODO/HACK: why is this correction for Xdata.mH necessary?
}

/**
 * @brief Computes the vapor density which are given by the following formula: \n
 * @par
 *            drho_v            d^2rho_v
 *   		  ----     -  D   * ------   =  Q
 *             dt                dz^2
 * \n
 * with initial and Dirichlet and/or Neumann boundary conditions.
 *
 * rho_v(z,t) = vapor density (kg m-3);
 * rho = snow density (kg m-3); c = specific heat capacity (J K-1 kg-1); k = heat conductivity (W K-1 m-1); t = time (s);
 * \n
 * Note:  The equations are solved with a fully implicit time-integration scheme and the
 * system of finite element matrices are solved using a sparse matrix solver.
 * @param Xdata Snow profile data
 * @param[in] Mdata Meteorological forcing
 * @param Bdata Boundary conditions
 * @param[in] ThrowAtNoConvergence	If true, throw exception when temperature equation does not converge; if false, function will return false after non convergence and true otherwise.
 * @return true when temperature equation converged, false if it did not.
 */
bool VapourTransport::compDensityProfile(const CurrentMeteo& Mdata, SnowStation& Xdata, const bool& ThrowAtNoConvergence, double& ql, const double& surfaceVaporPressure)
{
	bool topDirichletBCtype= false;
	bool bottomDirichletBCtype= false;
	//ChangingBoundaryToDirichlet: //label for goto to change B.C. to Dirichlet

	int Ie[N_OF_INCIDENCES];                     // Element incidences
	double T0[N_OF_INCIDENCES];                  // Element nodal vapor density at t0
	double TN[N_OF_INCIDENCES];                  // Iterated element nodal vapor densities
	double Se[N_OF_INCIDENCES][N_OF_INCIDENCES]; // Element stiffnes matrix
	double Fe[N_OF_INCIDENCES];                  // Element right hand side vector

	double *U=NULL, *dU=NULL, *ddU=NULL;         // Solution vectors

	// Dereference the pointers
	void *Kt_vapor = Xdata.Kt_vapor; //I think we neeed to define Xdata.Kt_vapor for vapor densiry profile....
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;

	const size_t nN = Xdata.getNumberOfNodes();
	const size_t nE = Xdata.getNumberOfElements();

	std::vector<double> factor_(nE, 1.);// this is for source term in vapor transport equation
	for(size_t i=0; i<Xdata.SoilNode; i++)
	{
		factor_[i]=0.;
	}
	for(size_t i=0; i<nN; i++)// update the vapor saturarion density according to the current temperature
	{
		double sVapor = Atmosphere::waterVaporDensity(NDS[i].T, Atmosphere::vaporSaturationPressure(NDS[i].T));
		NDS[i].rhov=sVapor;
	}

	if (Kt_vapor != NULL)
		ds_Solve(ReleaseMatrixData, (SD_MATRIX_DATA*)Kt_vapor, 0);
	ds_Initialize(static_cast<int>(nN), (SD_MATRIX_DATA**)&Kt_vapor);
	/*
	 * Define the structure of the matrix, i.e. its connectivity. For each element
	 * we compute the element incidences and pass the incidences to the solver.
	 * The solver assumes that the element incidences build a crique, i.e. the
	 * equations specified by the incidence set are all connected to each other.
	 * Initialize element data.
	*/
	for (int e = 0; e < static_cast<int>(nE); e++) {
		int Nodes[2] = {e, e+1};
		ds_DefineConnectivity( (SD_MATRIX_DATA*)Kt_vapor, 2, Nodes , 1, 0 );
	}

	/*
	 * Perform the symbolic factorization. By specifying the element incidences, we
	 * have simply declared which coefficients of the global matrix are not zero.
	 * However, when we factorize the matrix in a LU form there is some fill-in.
	 * Coefficients that were zero prior to start the factorization process will
	 * have a value different from zero thereafter. At this step the solver compute
	 * exactly how many memory is required to solve the problem and allocate this
	 * memory in order to store the numerical matrix. Then reallocate all the
	 * solution vectors.
	*/
	ds_Solve(SymbolicFactorize, (SD_MATRIX_DATA*)Kt_vapor, 0);

	// Make sure that these vectors are always available for use ....
	errno=0;
	U=(double *) realloc(U, nN*sizeof(double));
	if (errno != 0 || U==NULL) {
		free(U);
		prn_msg(__FILE__, __LINE__, "err", Date(), "%s (allocating  solution vector U)", strerror(errno));
		throw IOException("Runtime error in compDensityProfile", AT);
	}
	dU=(double *) realloc(dU, nN*sizeof(double));
	if (errno != 0 || dU==NULL) {
		free(U); free(dU);
		prn_msg(__FILE__, __LINE__, "err", Date(), "%s (allocating  solution vector dU)", strerror(errno));
		throw IOException("Runtime error in compDensityProfile", AT);
	}
	ddU=(double *) realloc(ddU, nN*sizeof(double));
	if (errno != 0 || ddU==NULL) {
		free(U); free(dU); free(ddU);
		prn_msg(__FILE__, __LINE__, "err", Date(), "%s (allocating  solution vector ddU)", strerror(errno));
		throw IOException("Runtime error in compDensityProfile", AT);
	}

	// Make sure that the global data structures know where the pointers are for the next integration step after the reallocation ....
	Xdata.Kt_vapor = Kt_vapor;

	// Copy Temperature at time0 into First Iteration
	for (size_t n = 0; n < nN; n++) {
		U[n] = NDS[n].rhov;
		dU[n] = 0.0;
		ddU[n] = 0.0;
		/*if (!(U[n] > 0. && U[n] < 10.)) {
		if (0) {
			prn_msg(__FILE__, __LINE__, "err", Mdata.date, "vapor density out of bound at beginning of iteration!");
			prn_msg(__FILE__, __LINE__, "msg", Date(), "At node n=%d (nN=%d, SoilNode=%d): T=%.2lf", n, nN, Xdata.SoilNode, U[n]);

			free(U); free(dU); free(ddU);
			throw IOException("Runtime error in compDensityProfile", AT);
		}*/
	}

	// Set the iteration counters, as well as the phase change boolean values
	unsigned int iteration = 0;   // iteration counter (not really required)
	bool NotConverged = true;     // true if iteration did not converge
	// Set the default solution routine convergence parameters
	unsigned int MaxItnTemp = 10000; // maximum 40 iterations for temperature field
	double ControlTemp = 1.e-12;    // solution convergence to within 0.01 degC

	// The temperature equation was found to show slow convergence with only 1 or 2 elements left.
	// Likely, the reason is that the LW-radiation is only approximated as linear, but in reality it is not. When only 1 or 2 elements
	// are left, their temperature gets very sensitive to energy input and during the iterations, the temperature gets out of the
	// validity range for the linearization. Therefore, we increase the MaxItnTemp for these cases:
	if (nN==3) MaxItnTemp = 10000;
	if (nN==2) MaxItnTemp = 10000;

	// IMPLICIT INTEGRATION LOOP
	bool TempEqConverged = true;	// Return value of this function compTemperatureProfile(...)
	do {
		iteration++;
		// Reset the matrix data and zero out all the increment vectors
		ds_Solve(ResetMatrixData, (SD_MATRIX_DATA*)Kt_vapor, 0);
		for (size_t n = 0; n < nN; n++) {
			ddU[n] = dU[n];
			dU[n] = 0.0;
		}

		// Assemble matrix
		for(int e = static_cast<int>(nE); e -->0; ) {
			EL_INCID( e, Ie );
			EL_TEMP( Ie, T0, TN, NDS, U );
			// Now fill the matrix
			if (!sn_ElementKtMatrix(EMS[e], sn_dt, T0, Se, Fe, factor_, e)) {
				prn_msg(__FILE__, __LINE__, "msg+", Mdata.date, "Error in sn_ElementKtMatrix @ element %d:", e);
				for (size_t n = 0; n < nN; n++)
					fprintf(stdout, "U[%u]=%g K\n", (unsigned int)n, U[n]);
				free(U); free(dU); free(ddU);
				throw IOException("Runtime error in compDensityProfile", AT);
			}
			ds_AssembleMatrix( (SD_MATRIX_DATA*)Kt_vapor, 2, Ie, 2,  (double*) Se );
			EL_RGT_ASSEM( dU, Ie, Fe );
		}

		// the upper B.C.
		if(topDirichletBCtype){
			double p_vapor = surfaceVaporPressure; //Mdata.rh * Atmosphere::vaporSaturationPressure( Mdata.ta );
			NDS[nE].rhov = Atmosphere::waterVaporDensity(NDS[nE].T, p_vapor); //Jafari added vapor density of the new node
			double Big = Constants::big;	// big number for DIRICHLET boundary conditions)
			// Dirichlet BC at surface: prescribed temperature value
			// NOTE Insert Big at this location to hold the temperature constant at the prescribed value.
			Ie[0] = static_cast<int>(nE);
			ds_AssembleMatrix((SD_MATRIX_DATA*) Kt_vapor, 1, Ie, 1, &Big);
		}
		if(!topDirichletBCtype){
			EL_INCID(static_cast<int>(nE-1), Ie);
			EL_TEMP(Ie, T0, TN, NDS, U);
			double D_drhov_dn = ql/(0.5*Constants::lh_sublimation+0.5*Constants::lh_vaporization);
			neumannBoundaryConditions(Se, Fe, D_drhov_dn);
			ds_AssembleMatrix( (SD_MATRIX_DATA*)Kt_vapor, 2, Ie, 2,  (double*) Se );
			EL_RGT_ASSEM( dU, Ie, Fe );
		}
        // the lower B.C.
		if(bottomDirichletBCtype){
			double Big = Constants::big;	// big number for DIRICHLET boundary conditions)
			double elementSaturationVaporDensity=Atmosphere::waterVaporDensity(EMS[0].Te, Atmosphere::vaporSaturationPressure(EMS[0].Te));
			NDS[0].rhov=elementSaturationVaporDensity;
			// Bottom node
			// Dirichlet BC at bottom: prescribed temperature value
			// NOTE Insert Big at this location to hold the temperature constant at the prescribed value.
			Ie[0] = 0;
			ds_AssembleMatrix((SD_MATRIX_DATA*) Kt_vapor, 1, Ie, 1, &Big);
		}
		if(!bottomDirichletBCtype){
			EL_INCID(0, Ie);
			EL_TEMP(Ie, T0, TN, NDS, U);
			double zeroFlux =0.0;
			neumannBoundaryConditions(Se, Fe, zeroFlux);
			ds_AssembleMatrix( (SD_MATRIX_DATA*)Kt_vapor, 2, Ie, 2,  (double*) Se );
			EL_RGT_ASSEM( dU, Ie, Fe );
		}

		/*
		 * Solve the linear system of equation. The te_F vector is used first as right-
		 * hand-side vector for the linear system. The solver stores in this vector
		 * the solution of the system of equations, the new temperature.
		 * It will throw an exception whenever the linear solver failed
		 */
		if (!ds_Solve(ComputeSolution, (SD_MATRIX_DATA*) Kt_vapor, dU)) {
			  prn_msg(__FILE__, __LINE__, "err", Mdata.date,
			  "Linear solver failed to solve for dU on the %d-th iteration.",
			  iteration);
			  throw IOException("Runtime error in compDensityProfile", AT);
		}
		// Update the solution vectors and check for convergence
		for (size_t n = 0; n < nN; n++)
			ddU[n] = dU[n] - ddU[n];
		double MaxTDiff = fabs(ddU[0]);  // maximum temperature difference for convergence
		for (size_t n = 1; n < nN; n++) {
			const double TDiff = fabs(ddU[n]); // temperature difference for convergence check
			if (TDiff > MaxTDiff)
				MaxTDiff = TDiff;
		}
		/*
		 * This is to increase the number of iterations for a phase changing uppermost element.
		 * Otherwise we would violate energy conservation because the implicit scheme
		 * leads to an adaptation of surface fluxes to the iterated surface temperature.
		 * In reality, the surface temperature will not change and therefore the fluxes
		 * must be constant. This means that the fluxes must be treated explicitely
		 * (see neumannBoundaryConditions)
		 */
//		if (U[nE] + ddU[nE] > EMS[nE-1].melting_tk || EMS[nE-1].theta[WATER] > 0.) {
//			ControlTemp = 0.007;
//			MaxItnTemp = std::max(MaxItnTemp, (unsigned)200); // NOTE originally 100;
//		}
		NotConverged = (MaxTDiff > ControlTemp);

		if (iteration > MaxItnTemp) {
			if (ThrowAtNoConvergence) {
				prn_msg(__FILE__, __LINE__, "err", Mdata.date,
				        "Temperature did not converge (azi=%.0lf, slope=%.0lf)!",
				        Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
				prn_msg(__FILE__, __LINE__, "msg", Date(),
				        "%d iterations > MaxItnTemp=%d; ControlTemp=%.4lf; nN=%d; sn_dt=%f",
				        iteration, MaxItnTemp, ControlTemp, nN, sn_dt);
				for (size_t n = 0; n < nN; n++) {
					if (n > 0)
						prn_msg(__FILE__, __LINE__, "msg-", Date(),
						        "U[%03d]:%6.1lf K, ddU:%8.4lf K;  NDS.T(t-1)=%6.1lf K; EMS[n-1].th_w(t-1)=%.5lf",
						        n, U[n], ddU[n], NDS[n].T, EMS[n-1].theta[WATER]);
					else
						prn_msg(__FILE__, __LINE__, "msg-", Date(),
						        "U[%03d]:%6.1lf K, ddU:%8.4lf K;  NDS.T(t-1)=%6.1lf K;",
						        n, U[n], ddU[n], NDS[n].T);
				}
				free(U); free(dU); free(ddU);
				throw IOException("Runtime error in compDensityProfile", AT);
			} else {
				TempEqConverged = false;	// Set return value of function
				NotConverged = false;		// Ensure we leave the do...while loop
			}
		}
		for (size_t n = 0; n < nN; n++) {
			U[n] += ddU[ n ];
		}
	} while ( NotConverged ); // end Convergence Loop

	if (TempEqConverged) {
		//bool prn_date = true;
		for (size_t n = 0; n < nN; n++) {
			if (U[n] > 0.) {
				NDS[n].rhov = U[n];
				//NDS[n].rhov = (0.999999*U[n]+0.000001*NDS[n].rhov);
			} else {
				//prn_msg(__FILE__, __LINE__, "wrn", Mdata.date,
				//        "negative vapor density, switching to Drichlet B.C. for top surface rho= %lf, nN=%d, index=%d", U[n], nN, n);
				//free(U); free(dU); free(ddU);
				//topDirichletBCtype=true;
				//goto ChangingBoundaryToDirichlet;
				NDS[n].rhov=U[n];
				//NDS[n].rhov = 0.;
			}
		}
		for (size_t e = 0; e < nE; e++) {
			EMS[e].rhov = (NDS[e].rhov + NDS[e+1].rhov) / 2.0;
		}
	}

	free(U); free(dU); free(ddU);

	return TempEqConverged;
}

/**
 * @brief Computes the element stiffness matrix and right hand side vector for a fully implicit time integration scheme \n
 * The matrices that must be evaluated are : \n
 * - [Se] = [Ce]/dt + [Ke] :  [Ce] is the element heat capacity matrix and [Ke] is the
 *                            element conductivity matrix.
 * - {Fe} = {Q} - [Ke]*{T0} : {Fe} is the element right-hand side vector containing the element heat flux. \n
 *                            {Q} arises from shortwave radiation -- the other heat fluxes are treated separately
 *                            when determining the meteo parameters.
 * @param Edata
 * @param dt Calculation time step length (s)
 * @param dvdz Wind pumping velocity gradient (s-1)
 * @param T0 Initial nodal temperatures (K)
 * @param Se Element heat capacitity (stiffness) matrix
 * @param Fe Element right hand side vector
 * @param VaporEnhance Vapor transport enhancement factor
 * @return false on error, true if no error occurred
 */
bool VapourTransport::sn_ElementKtMatrix(ElementData &Edata, double dt, double T0[ N_OF_INCIDENCES ], double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ], double Fe[ N_OF_INCIDENCES ], const std::vector<double> factor_, const int index)
{
	if (Edata.L < 0.0) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Negative length L=%e", Edata.L);
		return false;
	}

	//reset Fe_0 and Fe_1
	Fe[0] = Fe[1] = 0.0;

	double aaa = Edata.theta[AIR];
	double nnn = 1.- Edata.theta[SOIL];
	double D_vapSoil = Constants::diffusion_coefficient_in_air * pow(aaa,10./3.)/nnn/nnn; // based on jury1983
	double D = factor_[index]*Constants::diffusion_coefficient_in_snow + (1.0-factor_[index])*D_vapSoil;
	D = D/Edata.L;   //Conductivity. Divide by the length to save from doing it during the matrix operations

	// Evaluate the stiffness matrix
	Se[0][0] = Se[1][1] = D;
	Se[0][1] = Se[1][0] = -D;

	Fe[1] += 0.0; // this is the body source term

	// Add the implicit time integration term to the right hand side
	Fe[0] -= (Se[0][0] * T0[0] + Se[0][1] * T0[1]);
	Fe[1] -= (Se[1][0] * T0[0] + Se[1][1] * T0[1]);

	const double c = 1. * Edata.L * 1. / (6. * dt);
	Se[0][0] += 2. * c;
	Se[1][1] += 2. * c;
	Se[0][1] += c;
	Se[1][0] += c;

    double phaseChangeTerm = 0.0;
	Fe[0] += 0.5 * phaseChangeTerm * Edata.L;
	Fe[1] += 0.5 * phaseChangeTerm * Edata.L;

	return true;
}

/**
 * @brief Imposes Neumann boundary conditions at the surface. \n
 * The fluxes are first updated in updateBoundHeatFluxes. \n
 * This splitting ensures that these initial fluxes are used for each iteration in the
 * semi-explicit solution while they are not altered by the implicit solution.
 * Note that long wave radiation heat transfer is treated as a convection boundary condition
 * for the implicit solution, similar to the sensible heat exchange: \n
 *            lw_net = delta*(e*Ta - T_iter) \n
 * where delta is a function of both Ta and T_iter and is computed by SnLaws::compLWRadCoefficient
 * @param Mdata
 * @param Bdata
 * @param Xdata
 * @param T_snow Initial (snow) surface temperature (K)
 * @param T_iter Iterated (snow) surface temperature (K)
 * @param Se Element stiffness matrix
 * @param Fe Element right hand side vector
 */
void VapourTransport::neumannBoundaryConditions(double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ],
                                         double Fe[ N_OF_INCIDENCES ],
                                         double& X)
{
	// First zero out the interiour node contribution
	Se[0][0] = Se[0][1] = Se[1][0] = Se[1][1] = Fe[0] = Fe[1] = 0.0;
	Fe[1] += X;
}

/*
 * End of VapourTransport.cc
 */
