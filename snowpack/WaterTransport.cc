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

#include <snowpack/WaterTransport.h>
#include <snowpack/Snowpack.h>

#include <assert.h>

using namespace std;
using namespace mio;

WaterTransport::WaterTransport(const mio::Config& cfg)
                : variant(), thresh_rain(IOUtils::nodata), sn_dt(IOUtils::nodata),
                  hoar_thresh_rh(IOUtils::nodata), hoar_thresh_vw(IOUtils::nodata),
                  hoar_density_buried(IOUtils::nodata), hoar_density_surf(IOUtils::nodata), hoar_min_size_buried(IOUtils::nodata),
                  minimum_l_element(IOUtils::nodata), useSoilLayers(false), water_layer(false), jam(false)
{
	cfg.getValue("VARIANT", "SnowpackAdvanced", variant);

	// Defines whether soil layers are used
	cfg.getValue("SNP_SOIL", "Snowpack", useSoilLayers);

	//To build a thin top rain-water layer over a thin top ice layer, rocks, roads etc.
	cfg.getValue("WATER_LAYER", "SnowpackAdvanced", water_layer);

	//Rain only for air temperatures warmer than threshold (degC)
	cfg.getValue("THRESH_RAIN", "SnowpackAdvanced", thresh_rain);

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

	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	sn_dt = M_TO_S(calculation_step_length);

	//To build up a water table over impermeable layers
	cfg.getValue("JAM", "SnowpackAdvanced", jam);

	// Density of BURIED surface hoar (kg m-3), default: 125./ Antarctica: 200.
	cfg.getValue("HOAR_DENSITY_BURIED", "SnowpackAdvanced", hoar_density_buried);

	//Minimum surface hoar size to be buried (mm). Increased by 50% for Dirichlet bc.
	cfg.getValue("HOAR_MIN_SIZE_BURIED", "SnowpackAdvanced", hoar_min_size_buried);

	//Density of surface hoar (-> hoar index of surface node) (kg m-3)
	cfg.getValue("HOAR_DENSITY_SURF", "SnowpackAdvanced", hoar_density_surf);

	//Minimum element length (m)
	cfg.getValue("MINIMUM_L_ELEMENT", "SnowpackAdvanced", minimum_l_element);
}

/**
 * @brief This part of the code is EXTREMELY IMPORTANT -- especially for predicting SURFACE HOAR and BURIED DEPTH HOAR layers \n
 * The total latent heat flux is predicted.  If positive (and above a certain cutoff level) then there
 * is a good possibility that SURFACE HOAR crystal have grown.  Of course, if negative
 * then we are also loosing mass from the surface.  These are seperate routines since
 * they might want to be changed or updated in future. \n
 * Just before his wedding, when Michael was implementing solute transport as initiated by
 * Peter Waldner, he realized that sublimation and evaporation was not possible from blank
 * soil layers. So he changed the routine on 29 of April 2002. \n
 * Another very important case of WATER movement through the snowpack is the SUBLIMATION of
 * VAPOR; this piece of code was taken from phase change and placed here because there is the
 * good possibility that the an ELEMENT might be SUBLIMATED away. \n
 * TODO Revise description!
 * @param *Xdata
 * @param ql Latent heat flux (W m-2)
 * @param *Sdata
 * @param *Mdata
 */
void WaterTransport::compSurfaceSublimation(const CurrentMeteo& Mdata, double ql, SnowStation& Xdata,
                                            SurfaceFluxes& Sdata)
{
	double dL=0., dM=0.;     // Length and mass chamges
	double M=0.;             // Initial mass and volumetric content (water or ice)
	double Tss;              // Surface Temperature
	double hoar=0.0;         // Actual change in hoar mass
	double cH_old;           // Temporary variable to hold height of snow

	size_t nN = Xdata.getNumberOfNodes();
	size_t nE = nN-1;
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;
	Tss = NDS[nE].T;

	/*
	 * If there are elements and ql > 0:
	 * update densities and volumetric contents (ELEMENT data),
	 * add/subtract mass to MS_SUBLIMATION and/or MS_EVAPORATION,
	 * potential surface hoar formation will be tested at the end of this routine (NODAL data);
	*/
	if (ql > Constants::eps2) { // Add Mass
		const double melting_tk = (Xdata.getNumberOfElements()>0)? Xdata.Edata[Xdata.getNumberOfElements()-1].melting_tk : Constants::melting_tk;
		if (Tss < melting_tk) { // Add Ice
			dM = ql*sn_dt/Constants::lh_sublimation;
			Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += dM;
			hoar = dM;

			// In this case adjust properties of element, keeping snow density constant
			const double L_top = EMS[nE-1].L;
			const double theta_i0 = EMS[nE-1].theta[ICE];
			dL = dM/(EMS[nE-1].Rho); // length change
			if (nE == Xdata.SoilNode) {
				dL = 0.;
				dM = MIN(dM,EMS[nE-1].theta[AIR]*(Constants::density_ice*EMS[nE-1].L));
			}
			NDS[nE].z += dL + NDS[nE].u; NDS[nE].u = 0.0;
			EMS[nE-1].L0 = EMS[nE-1].L = L_top + dL;
			EMS[nE-1].E = EMS[nE-1].dE = EMS[nE-1].Ee = EMS[nE-1].Ev = EMS[nE-1].S = 0.0;
			EMS[nE-1].theta[ICE] *= L_top/EMS[nE-1].L;
			EMS[nE-1].theta[ICE] += dM/(Constants::density_ice*EMS[nE-1].L);
			EMS[nE-1].theta[WATER] *= L_top/EMS[nE-1].L;

			for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
				EMS[nE-1].conc[ICE][ii] *= L_top*theta_i0/(EMS[nE-1].theta[ICE]*EMS[nE-1].L);
			}
		} else {
			// Add Water
			const double theta_w0 = EMS[nE-1].theta[WATER];
			dM = ql*sn_dt/Constants::lh_vaporization;
			Sdata.mass[SurfaceFluxes::MS_EVAPORATION] += dM;
			if (nE == Xdata.SoilNode) {
				dM = MIN(dM,EMS[nE-1].theta[AIR]*(Constants::density_water*EMS[nE-1].L));
			}
			EMS[nE-1].theta[WATER] += dM/(Constants::density_water*EMS[nE-1].L);

			for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
				EMS[nE-1].conc[WATER][ii] *= theta_w0/EMS[nE-1].theta[WATER];
			}
		}
		EMS[nE-1].M += dM;
		assert(EMS[nE-1].M>=0.); //mass must be positive

		// Update remaining volumetric contents and density
		EMS[nE-1].theta[AIR] = MAX(0., 1.0 - EMS[nE-1].theta[WATER] - EMS[nE-1].theta[ICE] - EMS[nE-1].theta[SOIL]);
		EMS[nE-1].Rho = (EMS[nE-1].theta[ICE] * Constants::density_ice)
		                     + (EMS[nE-1].theta[WATER] * Constants::density_water)
		                         + (EMS[nE-1].theta[SOIL] * EMS[nE-1].soil[SOIL_RHO]);
	} else if (ql < -Constants::eps2 && nE > 0) {
		// If  there is water in some form and ql < 0, SUBLIMATE and/or EVAPORATE some mass off
		std::vector<double> M_Solutes(Xdata.number_of_solutes, 0.); // Mass of solutes from disappearing phases
		size_t e = nE;
		while ((e > 0) && (ql < -Constants::eps2)) {  // While energy is available
			e--;
			/*
			* Determine the amount of potential sublimation/evaporation and collect some variables
			* that will be continuously used: L0 and M
			*  - NOTE: if water is present, evaporate first, then sublimate ice matrix.
			*          Otherwise sublimate ice matrix only.
			*/
			const double L0 = EMS[e].L;
			// If there is water ...
			if (EMS[e].theta[WATER] > 0.) {
				const double theta_w0 = EMS[e].theta[WATER];
				dM = ql*sn_dt/Constants::lh_vaporization;
				M = theta_w0*Constants::density_water*L0;
				// Check that you only take the available mass of water
				if (-dM >= M) {
					dM = -M;
					// Add solutes to Storage
					for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
						M_Solutes[ii] += EMS[e].conc[WATER][ii]*theta_w0*L0;
					}
					EMS[e].theta[WATER] = 0.0;
				} else {
					EMS[e].theta[WATER] += dM/(Constants::density_water*L0);
					for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
						EMS[e].conc[WATER][ii] *= theta_w0/EMS[e].theta[WATER];
					}
				}
				EMS[e].M += dM;
				assert(EMS[e].M>=0.); //mass must be positive
				Sdata.mass[SurfaceFluxes::MS_EVAPORATION] += dM;
				ql -= dM*Constants::lh_vaporization/sn_dt; // Update the energy used
			}
			if (ql < -Constants::eps2) {
				// If there is no water or if there was not enough water ...
				const double theta_i0 = EMS[e].theta[ICE];
				M = theta_i0*Constants::density_ice*L0;
				dM = ql*sn_dt/Constants::lh_sublimation;
				if (-dM > M) {
					dM = -M;
					// Add solutes to Storage
					for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
						M_Solutes[ii] += EMS[e].conc[ICE][ii]*theta_i0*L0;
					}
					EMS[e].theta[ICE]=0.0; dL = 0.;
				} else {
					dL = dM/(EMS[e].Rho);
					if (e < Xdata.SoilNode) {
						dL = 0.;
					}
					NDS[e+1].z += dL; EMS[e].L0 = EMS[e].L = L0 + dL;
					NDS[e+1].z += NDS[e+1].u; NDS[e+1].u = 0.0;

					EMS[e].E = EMS[e].dE = EMS[e].Ee = EMS[e].Ev = EMS[e].S = 0.0;
					EMS[e].theta[ICE] *= L0/EMS[e].L;
					EMS[e].theta[ICE] += dM/(Constants::density_ice*EMS[e].L);
					EMS[e].theta[WATER] *= L0/EMS[e].L;
					for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
						EMS[e].conc[ICE][ii] *= L0*theta_i0/(EMS[e].theta[ICE]*EMS[e].L);
					}
				}
				EMS[e].M += dM;
				 //if we remove the whole mass, we might have some small inconcistencies between mass and theta[ICE]*density*L -> negative
				//but the whole element will be removed anyway when getting out of here
				assert(EMS[e].M>=(-Constants::eps2));
				Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += dM;
				ql -= dM*Constants::lh_sublimation/sn_dt;     // Update the energy used

				// If present at surface, surface hoar is sublimated away
				if (e == nE-1) {
					hoar = dM;
				}
			}
			// Update remaining volumetric contents and density
			EMS[e].theta[AIR] = MAX(0., 1.0 - EMS[e].theta[WATER] - EMS[e].theta[ICE] - EMS[e].theta[SOIL]);
			EMS[e].Rho = (EMS[e].theta[ICE] * Constants::density_ice) + (EMS[e].theta[WATER] * Constants::density_water) + (EMS[e].theta[SOIL] * EMS[e].soil[SOIL_RHO]);
			assert(EMS[e].Rho>=0.); //we want positive density
		}
		// Now take care of left over solute mass.
		if (e == 0) { // Add Solute Mass to Runoff TODO HACK CHECK
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
	}
	// Check for surface hoar destruction or formation (once upon a time ml_sn_SurfaceHoar)
	if ((Mdata.rh > hoar_thresh_rh) || (Mdata.vw > hoar_thresh_vw) || (Mdata.ta >= C_TO_K(thresh_rain))) {
		hoar = MIN(hoar,0.);
	}
	Sdata.hoar += hoar;
	NDS[nN-1].hoar += hoar;
	if (NDS[nN-1].hoar < 0.) {
		NDS[nN-1].hoar = 0.;
	}
	for (size_t e = 0; e<nE-1; e++) {
		if (Xdata.Edata[e].theta[WATER] > 0.) {
			NDS[e+1].hoar = 0.;
		}
	}
	// At the end also update the overall height
	cH_old = Xdata.cH;
	Xdata.cH = NDS[Xdata.getNumberOfNodes()-1].z + NDS[Xdata.getNumberOfNodes()-1].u;
	Xdata.mH -= (cH_old - Xdata.cH);
}

/**
 * @brief Remove snow elements that either have no more ice or are too thin \n
 * -# Starting from the surface, remove or join as many snow elements as you can; check also the lowest snow element
 * -# Remove or join:
 * 	- Remove watery elements by adding their contents to the lower element \n
 *	  => Use sn_ReduceNumberElements() to get rid of unwanted elements
 * 	- Join too thin elements with their lower neighbour; keep buried SH and tagged layers longer on \n
 * 	  but enforce join for uncovered but previously buried SH \n
 * 	  => Use SnowStation::mergeElements() to compute the properties of the lower element
 * - NOTE
 * 	- Water will be transported AFTER elements have been removed.
 * 	- If WATER_LAYER is set with soil, make sure that you keep a potential wet water layer over soil or ice
 * 	- Reset ground surface temperature if no snow is left and there is no soil
 * @param *Xdata
 * @param *Sdata
 * @param ta Air temperature (K)
 */
void WaterTransport::removeElements(SnowStation& Xdata, SurfaceFluxes& Sdata)
{
	bool enforce_join = false;  // To enforce merging in special cases

	size_t nN = Xdata.getNumberOfNodes();
	size_t rnN = nN, nE = nN-1, rnE = nN-1;
	vector<ElementData>& EMS = Xdata.Edata;

	if ((nN == Xdata.SoilNode+1)
	        || (water_layer && useSoilLayers && (nN == Xdata.SoilNode+2)
	                && (EMS[nE-1].theta[ICE] < Snowpack::min_ice_content)
	                    && (EMS[nE-1].theta[WATER] > SnowStation::thresh_moist_snow) && (EMS[nE-1].L > 0.0001))) {
		return;
	}

	size_t eUpper = nE; // Index of the upper element, the properties of which will be transferred to the lower adjacent one
	while (eUpper-- > Xdata.SoilNode) {
		enforce_join = true;
		if ((EMS[eUpper].L < minimum_l_element) || (EMS[eUpper].mk%100 == 3)) {
			if ((EMS[eUpper].mk >= 100) && (EMS[eUpper].L >= 0.5 * minimum_l_element)) {
				enforce_join = false;
			}
			if (EMS[eUpper].mk%100 == 3) {
				if ((eUpper < nE-1) && (EMS[eUpper].L >= MM_TO_M(0.75*hoar_min_size_buried
				                                     * (hoar_density_surf/hoar_density_buried)))
				                    && (EMS[eUpper].Rho <= 300.)) {
					enforce_join = false;
				} else {
					enforce_join = true;
				}
			}
		} else {
			enforce_join = false;
		}
		if (((EMS[eUpper].theta[ICE] < Snowpack::min_ice_content) || enforce_join)
		       && (EMS[eUpper].theta[SOIL] < Constants::eps2)
		           && (EMS[eUpper].mk % 100 != 9)) {  // no PLASTIC or WATER_LAYER please
			if (eUpper > Xdata.SoilNode) { //If we have snow elements below to join with
				SnowStation::mergeElements(EMS[eUpper-1], EMS[eUpper], enforce_join);
			} else {
				enforce_join=false;
				if (eUpper==Xdata.SoilNode && Xdata.SoilNode>0.) {
					// In case of soil and removal of first snow element above soil:
					SnowStation::mergeElements(EMS[eUpper-1], EMS[eUpper], enforce_join);
				}
				// route liquid water and solute load to runoff
				Sdata.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF] += EMS[eUpper].M;
				if (Xdata.SoilNode == 0) { // In case of no soil
					Sdata.mass[SurfaceFluxes::MS_SOIL_RUNOFF] += EMS[eUpper].M;
					for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
						Sdata.load[ii] += (EMS[eUpper].L*EMS[eUpper].theta[WATER]*EMS[eUpper].conc[WATER][ii]
						    + EMS[eUpper].L*EMS[eUpper].theta[ICE] * EMS[eUpper].conc[ICE][ii]) / S_TO_H(sn_dt);
					}
				}
			}
			rnE--;
			rnN--;
			EMS[eUpper].Rho = Constants::undefined;
			if (!enforce_join)
				EMS[eUpper].L *= -1.;
			if ((eUpper < nE-1) && (EMS[eUpper+1].Rho < 0.) && (EMS[eUpper+1].L > 0.)) {
				EMS[eUpper+1].L *= -1.;
			}
		}
	}
	if (rnE < nE) {
		Xdata.reduceNumberOfElements(rnE);
		if (!useSoilLayers && (rnE == Xdata.SoilNode)) {
			Xdata.Ndata[Xdata.SoilNode].T = MIN(Constants::melting_tk, Xdata.Ndata[Xdata.SoilNode].T);
		}
	}

	if(rnE>=Xdata.SoilNode) {
		Xdata.ColdContent = 0.;
		for (size_t e=Xdata.SoilNode; e<rnE; e++) {
			Xdata.ColdContent += EMS[e].coldContent();
		}
	}
}

/**
 * @brief Surface sublimation and melt artificially create surface elements that have a much
 * too low density and this needs to be corrected. \n
 * TODO Check description!
 * @param Xdata
 */
void WaterTransport::adjustDensity(SnowStation& Xdata)
{
	size_t nN = Xdata.getNumberOfNodes();
	if (nN == Xdata.SoilNode + 1) return;

	size_t nE = nN-1;
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;

	size_t e = nE;
	while ((e-- > Xdata.SoilNode) && (EMS[e].theta[SOIL] < Constants::eps2) && (EMS[e].mk%100 > 9)
	           && (EMS[e].theta[WATER] < 0.3) && (EMS[e].theta[ICE] < Snowpack::min_ice_content)
	               && (EMS[e].L > minimum_l_element)) {
		const double L = EMS[e].L;
		double dL=0.; // Element length change
		// For water layer go to water density
		if (water_layer && (EMS[e].theta[WATER] < 0.95) && (EMS[e].theta[ICE] < Snowpack::min_ice_content)
			&& (e > Xdata.SoilNode) && ((EMS[e-1].theta[SOIL] > 0.95) || (EMS[e-1].theta[ICE] > 0.95))) {
			dL = -L * (1. - EMS[e].theta[WATER]);
			EMS[e].theta[WATER] = 1.;
			EMS[e].theta[ICE] = 0.;
		} else if (variant == "JAPAN") {
			//NIED (H. Hirashima) //Fz: Please check this adaptation still works as you want it to work!
			double multif = 0.05 / EMS[e].theta[ICE];
			dL  = -L*((multif-1)/multif);
			EMS[e].theta[WATER] *= multif;
			EMS[e].theta[ICE]   *= multif;
		} else {
			dL  = -L / 3.; //TODO check whether this approach is correct even though it is a "very old SNOWPACK approach"
			EMS[e].theta[WATER] *= 1.5;
			EMS[e].theta[ICE]   *= 1.5;
		}

		for (size_t eAbove = e; eAbove < nE; e++) {
			NDS[eAbove+1].z += dL + NDS[eAbove+1].u;
			NDS[eAbove+1].u = 0.0;
			EMS[eAbove].E = EMS[eAbove].dE= EMS[eAbove].Ee=EMS[eAbove].Ev=0.0;
		}
		EMS[e].L0 = EMS[e].L = L + dL;
		EMS[e].theta[AIR] = 1.0 - EMS[e].theta[WATER] - EMS[e].theta[ICE];
		EMS[e].Rho = (EMS[e].theta[ICE]*Constants::density_ice) + (EMS[e].theta[WATER]*Constants::density_water);
		if (!(EMS[e].Rho > Constants::min_rho && EMS[e].Rho <= Constants::max_rho)) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "Volume contents: e:%d nE:%d rho:%lf ice:%lf wat:%lf air:%le",
			        e, nE, EMS[e].Rho, EMS[e].theta[ICE], EMS[e].theta[WATER], EMS[e].theta[AIR]);
			throw IOException("Cannot evaluate mass balance in adjust density routine", AT);
		}
	}
	double cH_old = Xdata.cH;
	Xdata.cH = NDS[Xdata.getNumberOfNodes()-1].z + NDS[Xdata.getNumberOfNodes()-1].u;
	Xdata.mH -= (cH_old - Xdata.cH);
}

/**
 * @brief Default version of WaterTransport \n
 * Now that the snowpack has been updated, you must move the water down and compute the amount
 * of water being released from the snowpack AFTER having gone through removeElements() \n
 * TODO Revise description!
 * @param *Xdata
 * @param *Sdata
 * @param *Mdata
 */
void WaterTransport::transportWater(const CurrentMeteo& Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata)
{
	double Wres;          // Residual water content depending on snow or soil element
	double Store;         // Depth of liquid precipitation ready to infiltrate snow and/or soil (m)
	double z_water;       // Position of upper node of top water-film layer (m)

	size_t nN = Xdata.getNumberOfNodes();
	size_t nE = nN-1;
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;

	// First, consider no soil with no snow on the ground
	if (!useSoilLayers && nN == 1) {
		return;
	} else { // add rainfall to snow/soil pack
		if ((Mdata.hnw > 0.) && (Mdata.ta >= C_TO_K(thresh_rain))) {
			Store = Mdata.hnw / Constants::density_water;
			// Now find out whether you are on an impermeable surface and want to create a water layer ...
			if (water_layer && (Store > 0.)
			        && ((useSoilLayers && (nE == Xdata.SoilNode)
			                && (EMS[nE-1].theta[SOIL] > 0.95)) || ((nE-1 > 0) && (EMS[nE-2].theta[ICE] > 0.95)))) {
				nE++;
				nN++;
				Xdata.ErosionLevel = nE-1;
				Xdata.resize(nE);

				// Temperature of the uppermost node
				NDS[nN-1].T = Mdata.ta;
				// The new nodal position
				z_water = MIN(Store, MAX(0.001, 0.01 * cos(DEG_TO_RAD(Xdata.meta.getSlopeAngle()))));
				NDS[nN-1].z = NDS[nN-2].z + NDS[nN-2].u + z_water;
				Store -= z_water;
				// Fill the element data
				EMS[nE-1].depositionDate = Mdata.date;

				EMS[nE-1].Te = Mdata.ta;
				EMS[nE-1].L0 = EMS[nE-1].L = z_water;
				EMS[nE-1].Rho = Constants::density_water;
				EMS[nE-1].M = EMS[nE-1].L0 * EMS[nE-1].Rho;
				assert(EMS[nE-1].M>=0.); //mass must be positive
				EMS[nE-1].theta[WATER] = 1.0;
				EMS[nE-1].mk = 19;
				//NOTE all other microstructure parameters should better be set to Constants::undefined but ...
				EMS[nE-1].N3 = 1.;
				EMS[nE-1].dd = 0.;
				EMS[nE-1].sp = 1.;
				EMS[nE-1].rg = 1.0;
				EMS[nE-1].rb = 0.5;
				Xdata.cH = Xdata.mH = NDS[nN-1].z + NDS[nN-1].u;
			} else if (water_layer && (Store > 0.)
			               && ((useSoilLayers && (nE == Xdata.SoilNode+1) && (EMS[nE-2].theta[SOIL] > 0.95))
			                       || ((nE > 1) && (EMS[nE-2].theta[ICE] > 0.95)))) {
				// Put rain water in existing wet layer
				z_water = MIN(Store, MAX(0.0, (0.01 * cos(DEG_TO_RAD(Xdata.meta.getSlopeAngle())) - EMS[nE-1].L)));
				NDS[nN-1].z += z_water;
				Store -= z_water;
				EMS[nE-1].L0 = EMS[nE-1].L = (NDS[nN-1].z + NDS[nN-1].u) - (NDS[nN-2].z + NDS[nN-2].u);
				EMS[nE-1].M = EMS[nE-1].L0 * EMS[nE-1].Rho;
				assert(EMS[nE-1].M>=0.); //mass must be positive
				Xdata.cH = Xdata.mH = NDS[nN-1].z + NDS[nN-1].u;
			}

			//Put rain water in the layers, starting from the top element.
			size_t e = nE;
			while (Store > 0.0 && e > 0.) {
				e--;	//This construct with e=nE and e-- is to circumvent troubles with the unsigned ints.
				const double L = EMS[e].L;
				double dThetaW;
				if(e>0) {	//Check how much space is left to put rain water in.
					dThetaW = MIN(Constants::density_ice/Constants::density_water*EMS[e].theta[AIR],Store/L);
				} else {	//For the lowest layer, put all rain water in, to prevent loosing the mass.
					dThetaW = Store / L;
				}
				Store -= dThetaW*L;
				for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[e].conc[WATER][ii] = (L * dThetaW * Mdata.conc[ii]
					                                  + L * EMS[e].theta[WATER] * EMS[e].conc[WATER][ii])
					                              / (L * (EMS[e].theta[WATER] + dThetaW));
				}
				EMS[e].theta[WATER] += dThetaW;
				EMS[e].theta[AIR] -= dThetaW;
				EMS[e].M += dThetaW * L * Constants::density_water;
				assert(EMS[e].M>=0.); //mass must be positive
				// Update snowpack runoff with rain infiltrating into soil (equal to Store when e == Xdata.SoilNode)
				if (e == Xdata.SoilNode) {
					Sdata.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF] += Store * Constants::density_water;
				}
				// Update soil runoff with rain (equal to Store when e == 0)
				if (e == 0) {
					Sdata.mass[SurfaceFluxes::MS_SOIL_RUNOFF] += Store * Constants::density_water;
				}
			}
			Sdata.mass[SurfaceFluxes::MS_RAIN] += Mdata.hnw;
		}
	}

	// Preferential flow system: excess water that cannot be retained in lower element is stored in
	//   excess_water and moved down the domain; NOTE units: [m^3/m^2]
	double excess_water=0.;

	// Now move water as needed, starting from the top element ...
	for (size_t eUpper = nE-1, eLower = nE-2; eUpper >= 1; eUpper--, eLower-- ) {
		// Determine the additional storage capacity due to refreezing
		const double dth_w = EMS[eUpper].c[TEMPERATURE] * EMS[eUpper].Rho / Constants::lh_fusion / Constants::density_water
		                         * MAX(0., EMS[eUpper].melting_tk-EMS[eUpper].Te);
		if ((eUpper == nE-1) && (EMS[eLower].theta[AIR] <= 0.05) && water_layer) {
			// allow for a water table in the last layer above road/rock
			Wres = Constants::density_ice/Constants::density_water
			           * (1. - EMS[eUpper].theta[ICE] - EMS[eUpper].theta[SOIL] - 0.05);
		} else if (EMS[eUpper].theta[SOIL] < Constants::eps2) {
			Wres = MIN((1. - EMS[eUpper].theta[ICE]) * Constants::density_ice / Constants::density_water,
			           EMS[eUpper].res_wat_cont + dth_w);
		} else { // treat soil separately
			Wres = MIN(Constants::density_ice / Constants::density_water
			               * (1. - EMS[eUpper].theta[ICE] - EMS[eUpper].theta[SOIL]),
			           EMS[eUpper].soilFieldCapacity() + dth_w);
		}
		Wres = MAX (0., Wres);

		const double W_upper = EMS[eUpper].theta[WATER];
		if (eUpper == nE-1 && (W_upper > 0.0 && W_upper <= Wres)) {
			// In that case you need to update the volumetric air content and the density of the top element
			// as it may have caught some rain! Only top element should be considered, as when rain would have
			// infiltrated lower elements as well, W_upper>Wres.
			EMS[eUpper].theta[AIR] = MAX(0., 1. - EMS[eUpper].theta[WATER] - EMS[eUpper].theta[ICE] - EMS[eUpper].theta[SOIL]);
			EMS[eUpper].Rho = (EMS[eUpper].theta[ICE] * Constants::density_ice)
			                   + (EMS[eUpper].theta[WATER] * Constants::density_water)
			                       + (EMS[eUpper].theta[SOIL] * EMS[eUpper].soil[SOIL_RHO]);
			assert(EMS[eUpper].Rho>=0.); //we want positive density
			if ( EMS[eUpper].theta[SOIL] < Constants::eps2 ) {
				if ( !(EMS[eUpper].Rho > Constants::min_rho && EMS[eUpper].Rho <= Constants::max_rho) ) {
					prn_msg(__FILE__, __LINE__, "err", Mdata.date,
					        "Volume contents: e:%d nE:%d rho:%lf ice:%lf wat:%lf air:%le",
					        eUpper, nE, EMS[eUpper].Rho, EMS[eUpper].theta[ICE], EMS[eUpper].theta[WATER], EMS[eUpper].theta[AIR]);
					throw IOException("Cannot transfer water within the snowpack in transportWater()", AT);
				}
			}
		}

		if (W_upper > Wres || excess_water > 0.) {
			// Then water is being transferred between elements
			const double L_upper = EMS[eUpper].L;
			const double L_lower = EMS[eLower].L;
			const double W_lower = EMS[eLower].theta[WATER];
			double dThetaW_upper = W_upper - Wres;

			// dThetaW_lower is determined by also taking excess_water into account. Maybe excess_water can be stored in this layer.
			dThetaW_upper = MAX(0, dThetaW_upper);
			if (dThetaW_upper > 0. || excess_water > 0.) {
				// dThetaW_lower is determined by also taking excess_water into account. Maybe excess_water can be stored in this layer.
				double dThetaW_lower = dThetaW_upper*(L_upper/L_lower)+(excess_water/L_lower);
				// Now check whether there is enough air left - in case of ice, rock or heavy
				// soil you might not be able to move the water or/and water may refreeze and expand.
				// Specifically, you might want to create a water table over ice or frozen soil
				if ((dThetaW_lower + W_lower) > (Constants::density_ice / Constants::density_water
				        * (1. - EMS[eLower].theta[ICE] - EMS[eLower].theta[SOIL]))) {
					// Deal with excess water ... Look how much you can leave in the lower element eLower.
					// If you have too much water even for the lower element (more melt or rain per time
					// step than can be kept in this element), water is transferred to excess_water.
					// excess_water moves the water downward, trying to insert the water in lower elements.
					const double i_dThetaW_lower=dThetaW_lower;		//Make backup
					dThetaW_lower = MAX(0., (Constants::density_ice/Constants::density_water
					                        * (1. - EMS[eLower].theta[ICE] - EMS[eLower].theta[SOIL]) - W_lower));
					if (jam) {
						// In case of jam, the change in water content in the upper layer is the maximum possible
						// change in water content in the lower layer minus excess_water. This means excess_water
						// has the right of way to fill up the lower element. If the lower element cannot contain
						// all excess_water, it is transferred to the upper element. TODO: actually, if we fill
						// the upper layer completely and we still have excess_water left, we should put it in
						// the layer above the upper layer and so forth.
						dThetaW_upper = dThetaW_lower*L_lower/L_upper - excess_water/L_upper;
						if ((W_upper - dThetaW_upper) > (Constants::density_ice/Constants::density_water
						                           * (1. - EMS[eUpper].theta[ICE] - EMS[eUpper].theta[SOIL]))) {
							const double i_dThetaW_upper = dThetaW_upper;	//Make backup
							dThetaW_upper = W_upper - Constants::density_ice/Constants::density_water
							                    * (1. - EMS[eUpper].theta[ICE] - EMS[eUpper].theta[SOIL]);
							excess_water = (i_dThetaW_upper-dThetaW_upper)*L_upper;
						} else {
							excess_water = 0.;
						}
						if (EMS[eLower].theta[SOIL] < Constants::eps2) {
							Sdata.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF] += excess_water*Constants::density_water;
						}
						Sdata.mass[SurfaceFluxes::MS_SOIL_RUNOFF] += excess_water*Constants::density_water;
						// Take care of Solutes
						for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
							Sdata.load[ii] += (EMS[eLower].conc[WATER][ii] * excess_water
						                           * Constants::density_water/S_TO_H(sn_dt));
						}
						// We have moved the excess_water out of the snowpack.
						excess_water=0.;
					} else {
						// All the water that could not be stored in eLower is considered excess_water.
						excess_water = (i_dThetaW_lower-dThetaW_lower)*L_lower;
					}
				} else {
					// else EMS[eLower] can contain all water, so we have no excess_water anymore.
					excess_water=0.;
				}

				// Water movement from element eUpper to element eLower: move solutes also
				for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[eLower].conc[WATER][ii] = (W_lower * EMS[eLower].conc[WATER][ii] + dThetaW_lower * EMS[eUpper].conc[WATER][ii])
				                                  / (W_lower+dThetaW_lower);
				}

				// update volumetric contents, masses and density
				EMS[eUpper].theta[WATER]=W_upper-dThetaW_upper;
				EMS[eLower].theta[WATER]=W_lower+dThetaW_lower;
				EMS[eUpper].theta[AIR] = 1. - EMS[eUpper].theta[WATER] - EMS[eUpper].theta[ICE] - EMS[eUpper].theta[SOIL];
				EMS[eLower].theta[AIR] = 1. - EMS[eLower].theta[WATER] - EMS[eLower].theta[ICE] - EMS[eLower].theta[SOIL];
				EMS[eUpper].M -= L_upper * Constants::density_water * dThetaW_upper;
				assert(EMS[eUpper].M>=0.); //mass must be positive
				EMS[eLower].M += L_lower * Constants::density_water * dThetaW_lower;
				assert(EMS[eLower].M>=0.); //mass must be positive
				EMS[eUpper].Rho = (EMS[eUpper].theta[ICE] * Constants::density_ice)
				                  + (EMS[eUpper].theta[WATER] * Constants::density_water)
				                      + (EMS[eUpper].theta[SOIL] * EMS[eUpper].soil[SOIL_RHO]);
				assert(EMS[eUpper].Rho>=0.); //we want positive density
				EMS[eLower].Rho = (EMS[eLower].theta[ICE] * Constants::density_ice)
				                  + (EMS[eLower].theta[WATER] * Constants::density_water)
				                      + (EMS[eLower].theta[SOIL] * EMS[eLower].soil[SOIL_RHO]);
				assert(EMS[eLower].Rho>=0.); //we want positive density
				if (EMS[eUpper].theta[SOIL] < Constants::eps2) {
					if (!(EMS[eUpper].theta[AIR] >= -Constants::eps)) {
						prn_msg(__FILE__, __LINE__, "err", Mdata.date,
						        "Volume contents: e:%d nE:%d rho:%lf ice:%lf wat:%lf air:%le",
						        eUpper, nE, EMS[eUpper].Rho, EMS[eUpper].theta[ICE], EMS[eUpper].theta[WATER], EMS[eUpper].theta[AIR]);
						throw IOException("Cannot transfer water within the snowpack in transportWater()", AT);
					}
				}
				// Update snowpack runoff with soil. Note: in case of no soil layers, or lowest soil element: the runoff for the lowest element is updated outside the loop.
				if (useSoilLayers && eUpper == Xdata.SoilNode) {
					Sdata.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF] += L_lower * Constants::density_water * dThetaW_lower + excess_water * Constants::density_water;
				}
			} // end positive water movement
		}  // end if( W_upper > Wres )
	}  // end FOR loop over the number of elements

	// The TOP element is very important because it is always losing mass--the strain state
	// is becoming more and more deformed.  Update the strain state of the top element.
	const size_t eTop = nE-1;
	EMS[eTop].L0 = EMS[eTop].L;
	NDS[nN-1].z += NDS[nN-1].u; NDS[nN-1].u = 0.0;
	NDS[nN-2].z += NDS[nN-2].u; NDS[nN-2].u = 0.0;
	EMS[eTop].E = EMS[eTop].dE = EMS[eTop].Ee = EMS[eTop].Ev = EMS[eTop].S = 0.0;

	// RUNOFF at bottom of either snowpack or soil
	// Determine the additional storage capacity due to refreezing
	const double dth_w = EMS[0].c[TEMPERATURE] * EMS[0].Rho / Constants::lh_fusion / Constants::density_water
	                         * MAX(0., EMS[0].melting_tk-EMS[0].Te);
	if (EMS[0].theta[SOIL] < Constants::eps2) {
		Wres = MIN((1. - EMS[0].theta[ICE]) * Constants::density_ice / Constants::density_water,
		           EMS[0].res_wat_cont + dth_w);
	} else { // treat soil separately
		Wres = MIN(Constants::density_ice/Constants::density_water*(1. - EMS[0].theta[ICE] - EMS[0].theta[SOIL]),
		       EMS[0].soilFieldCapacity() + dth_w);
	}
	Wres = MAX (0., Wres);

	const double W0 = EMS[0].theta[WATER];
	if ((W0 > Wres) // NOTE: if water_layer is set, do not drain water element on top of soil
	        && !(water_layer && (EMS[0].theta[ICE] < Snowpack::min_ice_content)
	                 && (EMS[0].theta[SOIL] < Constants::eps2))) {
		const double dM = EMS[0].L * Constants::density_water * (W0 - Wres);
		EMS[0].M -= dM;
		assert(EMS[0].M>=0.); //mass must be positive
		EMS[0].theta[WATER] = Wres;
		EMS[0].theta[AIR] = 1. - EMS[0].theta[WATER] - EMS[0].theta[ICE] - EMS[0].theta[SOIL];
		EMS[0].Rho = (EMS[0].theta[ICE] * Constants::density_ice)
		                 + (EMS[0].theta[WATER] * Constants::density_water)
		                     + (EMS[0].theta[SOIL] * EMS[0].soil[SOIL_RHO]);
		assert(EMS[0].Rho>=0.); //we want positive density
		// Note that remaining excess_water should also be routed to MS_SOIL_RUNOFF and MS_SNOWPACK_RUNOFF
		if (EMS[0].theta[SOIL] < Constants::eps2) {
			Sdata.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF] += dM + (excess_water * Constants::density_water);
		}
		Sdata.mass[SurfaceFluxes::MS_SOIL_RUNOFF] += dM + (excess_water * Constants::density_water);
		for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
			Sdata.load[ii] +=  (EMS[0].conc[WATER][ii] * dM / S_TO_H(sn_dt));
		}
	}
}

/**
 * @brief The mass transport procedure is called from sn_Snowpack -- AFTER calling
 * the NEWSNOW (sn_SnowFall) or SNOWDRIFT (sn_SnowDrift) modules but BEFORE
 * calling the TEMPERATURE (sn_SnowTemperature), PHASECHANGE (pc_PhaseChange)
 * or CREEP (sn_SnowCreep) routines. \n
 * The mass transport routines were inserted at this location since they can set the NEWMESH
 * variable which means the FEM data structures must be reallocated before solving the instationary heat equations(?) \n
 * These routines are responsible for MOVING MASS (water) in, out and through the
 * snowpack.  They are subsequently responsible for WATER TRANSPORT and SURFACE
 * SUBLIMATION.  Since surface sublimation does not change the FE data structure
 * it is treated FIRST. \n
 * The phase change routines will increment the volumetric water content of the
 * elements, then the WATER TRANSPORT routines will move the excess water from
 * element "e" to element "e-1". \n
 * NOTES:
 * -#  The water will only be moved if it is above the residual water content
 * -#   The water will only be moved if there is enough VOID SPACE in the element
 *      receiving the water
 * -#   Water in the last element will be DISCHARGED from the snowpack.  This
 *      amount of water is termed the MELTWATER RUNOFF
 * -#   It would be very possible to make the RESIDUAL_WATER_CONTENT a function
 *      of the MICRO-properties of the snow.  This is simply a few lines of code. \n
 *      This was done on 3 Dec 2006 -> ElementData::snowResidualWaterContent().
 * -#   IMPORTANT: the top surface element can be removed if the VOLUMETRIC ICE
 *      content is 0; that is, when the element contains only water and voids. \n
 * The routines were changed dramatically by Perry on June 3rd 1998 after Perry
 * and Michael worked the entirity of June 2nd together.  We were running
 * the model operationally in Davos and discovered that the code was bombing
 * during heavy melt periods.  Michael identified the problem that the elements
 * were being REMOVED after the water was TRANSPORTED.  Meaning that water contents
 * greater than the 0.03 were being picked up and if the melting was strong
 * enough, then negative volumetric AIR contents were predicted.  The solution
 * to the problem is simple enough: FIRST remove the elements, then caculate
 * the WATER TRANSPORT.
 * @param Xdata
 * @param ql Latent heat flux (W m-2)
 * @param Sdata
 * @param Mdata
 */
void WaterTransport::compTransportMass(const CurrentMeteo& Mdata, const double& ql,
                                       SnowStation& Xdata, SurfaceFluxes& Sdata)
{
	// First, consider no soil with no snow on the ground and deal with possible rain water
	if (!useSoilLayers && (Xdata.getNumberOfNodes() == Xdata.SoilNode+1)) {
		if (Mdata.ta >= C_TO_K(thresh_rain)) {
			Sdata.mass[SurfaceFluxes::MS_RAIN] += Mdata.hnw;
			Sdata.mass[SurfaceFluxes::MS_SOIL_RUNOFF] += Mdata.hnw;
			for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
				Sdata.load[ii] += Mdata.conc[ii] * Mdata.hnw /*/ S_TO_H(sn_dt)*/;
			}
		}
		return;
	}

	compSurfaceSublimation(Mdata, ql, Xdata, Sdata);
	removeElements(Xdata, Sdata);

	try {
		adjustDensity(Xdata);
		transportWater(Mdata, Xdata, Sdata);
	} catch(const exception&){
		prn_msg( __FILE__, __LINE__, "err", Mdata.date, "Error in transportMass()");
		throw;
	}
}

/*
 * End of WaterTransport.cc
 */
