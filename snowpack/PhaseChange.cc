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

#include <snowpack/PhaseChange.h>
#include <snowpack/Snowpack.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/

/*
 * Residual Water Content,  for now we say  0.0
 * - This constant is used in the PHASE CHANGE ROUTINES; moreover we will freeze
 *   ALL water in the snowpack.  This ensures that we can have DRY snow. (Perryanic comment!)
 */
const double PhaseChange::theta_r = 0.0;
const double PhaseChange::RE_theta_r = 1E-5/1000.;	//It is recommended that this value is REQUIRED_ACCURACY_THETA/1000. (see ReSolver1d.cc)
const double PhaseChange::RE_theta_threshold = 1E-5; 	//It is recommended that this value is REQUIRED_ACCURACY_THETA

//Saturated Water Content, for now we say 1.0
const double PhaseChange::theta_s = 1.0;

/************************************************************
 * non-static section                                       *
 ************************************************************/

PhaseChange::PhaseChange(const SnowpackConfig& cfg)
             : iwatertransportmodel_snow(BUCKET), iwatertransportmodel_soil(BUCKET),
	       watertransportmodel_snow("BUCKET"), watertransportmodel_soil("BUCKET"),
	       sn_dt(0.), cold_content_in(IOUtils::nodata), cold_content_out(IOUtils::nodata)
{
	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	sn_dt = M_TO_S(calculation_step_length);

	//Water transport model snow
	cfg.getValue("WATERTRANSPORTMODEL_SNOW", "SnowpackAdvanced", watertransportmodel_snow);
	if (watertransportmodel_snow=="BUCKET") {
		iwatertransportmodel_snow=BUCKET;
	} else if (watertransportmodel_snow=="NIED") {
		iwatertransportmodel_snow=NIED;
	} else if (watertransportmodel_snow=="RICHARDSEQUATION") {
		iwatertransportmodel_snow=RICHARDSEQUATION;
	}

	//Water transport model soil
	cfg.getValue("WATERTRANSPORTMODEL_SOIL", "SnowpackAdvanced", watertransportmodel_soil);
	if (watertransportmodel_soil=="BUCKET") {
		iwatertransportmodel_soil=BUCKET;
	} else if (watertransportmodel_soil=="NIED") {
		iwatertransportmodel_soil=NIED;
	} else if (watertransportmodel_soil=="RICHARDSEQUATION") {
		iwatertransportmodel_soil=RICHARDSEQUATION;
	}
}

/**
 * @brief Subsurface Melting:
 * -# Consistency check.
 * -# Temperature ice and water check: If the  element temperature is below melting
 * temperature, or if the volumetric ice content is nonpositive, or if the water content
 * equals or exceeds PhaseChange::theta_s, no melting occurs.
 * If the check is ok,  then the  difference dT between element and melting temperature
 * is computed.
 * -# The coefficient A is computed.
 * -# Analogous to Subsurface Freezing ...
 * -# Check if there is enough ice to be melted: If it is not the case, the values have to
 * be corrected.
 * -# If there is  TOO MUCH  water being melted  then d_th(w) ( = PhaseChange::theta_s - th(w) ) as welL
 * as the changes d_th(i) and dT are corrected.
 * -# Computation of contents.
 * -# Characterization of element ( Liquid, Dry, Void ).
 * -# Computation of densities and water generation rate.
 * @param *Edata
 * @param dt Time step (s)
 * @param *ql_Rest Latent heat flux balance (J m-2)
 */
void PhaseChange::compSubSurfaceMelt(ElementData& Edata, const unsigned int nSolutes, const double& dt,
                                     double& ql_Rest, const mio::Date& date_in)
{
	const double T_melt=Edata.melting_tk;		// Retrieve melting temperature from ElementData

	if(!Edata.checkVolContent()) prn_msg(__FILE__, __LINE__, "wrn", Date(), "wrong volumetric content");
	/*
	 * Now see if any melting is going on -- this implies that (1) the temperature of the element
	 * is above the melting temperature (2) there is something to melt and (3) there is enough
	 * room to place the meltwater ...
	*/
	if (((Edata.Te < T_melt) && (ql_Rest < Constants::eps2))
	        || (Edata.theta[ICE] <= 0.0) || (Edata.theta[WATER] >= PhaseChange::theta_s)) {
		return;
	} else {
		double dT = T_melt - Edata.Te; // Edata.melting_tk - Te > 0
		// Now we take into account that there might be some extra energy that could not
		// be used by the element above because of complete melting
		dT -= ql_Rest / (Edata.c[TEMPERATURE] * Edata.Rho * Edata.L);
		if (dT > 0.) {
			return;
		}
		// Determine the DECREASE in ice content and the INCREASE of water content
		// Adapt A to compute mass changes
		const double A = (Edata.c[TEMPERATURE] * Edata.Rho) / ( Constants::density_ice * Constants::lh_fusion );
		double dth_i = A * dT; // change in volumetric ice content
		double dth_w = - (Constants::density_ice / Constants::density_water) * dth_i; // change in volumetric water content
		// It could happen that there is enough energy available to melt more ice than is present.
		// You can only melt so much ice as is there ....
		if ( (Edata.theta[ICE] + dth_i) < 0.0 ) {
			dth_i = - Edata.theta[ICE];
			dth_w = - (Constants::density_ice / Constants::density_water) * dth_i;
			dT = dth_i / A;
		}
		// It could also be that you are trying to produce more water than is allowed.
		if ( (Edata.theta[WATER] + dth_w) > PhaseChange::theta_s ) {
			dth_w = PhaseChange::theta_s - Edata.theta[WATER];
			dth_i = - (Constants::density_water / Constants::density_ice) * dth_w;
			dT = dth_i / A;
		}
		// Reset element properties
		Edata.Te += dT;
		if (Edata.Te <= T_melt) {
			ql_Rest = 0.0;
			Edata.Te = T_melt;
		} else {
			ql_Rest = Edata.c[TEMPERATURE] * Edata.Rho * Edata.L * (Edata.Te - T_melt);
			Edata.Te = T_melt;
		}
		Edata.Qmf += (dth_i * Constants::density_ice * Constants::lh_fusion) / dt; // (W m-3)
		Edata.dth_w = dth_w; // (1)
		for (unsigned int ii = 0; ii < nSolutes; ii++) {
			if( dth_w > 0. ) {
				Edata.conc[WATER][ii] = (Edata.theta[WATER] * Edata.conc[WATER][ii]
				    + dth_w * Edata.conc[ICE][ii]) / (Edata.theta[WATER] + dth_w);
			}
		}
		Edata.theta[ICE] += dth_i;
		Edata.theta[WATER] += dth_w;
		Edata.theta[AIR] = MAX (0.0, 1.0 - Edata.theta[ICE] - Edata.theta[WATER] - Edata.theta[SOIL]);
		// State when you have solid element
		if ( Edata.theta[AIR] <= 0.0 ) {
			Edata.theta[AIR] = 0.0;
		}
		// State when the ice content has disappeared (PERMAFROST)
		if ( Edata.theta[ICE] <= 0.0 ) {
			Edata.theta[ICE] = 0.0;
		}
		// State when the water content has disappeared (PERMAFROST)
		if ( Edata.theta[WATER] <= 0.0 ) {
			Edata.theta[WATER] = 0.0;
		}
		// State when the element is wet (PERMAFROST)
		if ( Edata.theta[WATER] >= 1.0 ) {
			Edata.theta[WATER] = 1.0;
		}

		// Make sure the sum of all volumetric contents is near 1 (Can make a 1% error)
		if (!Edata.checkVolContent()) {
			prn_msg(__FILE__, __LINE__, "err", date_in, "Sum theta[I,W,A,S] > 1");
			prn_msg(__FILE__, __LINE__, "msg-", Date(),
			        "Ice: %f, Water: %f, Air: %f Soil: %f",
			        Edata.theta[ICE], Edata.theta[WATER], Edata.theta[AIR], Edata.theta[SOIL]);
			throw IOException("In compSubSurfaceMelt!", AT);
		}
		Edata.Rho = Constants::density_ice * Edata.theta[ICE]
		                + (Constants::density_water * Edata.theta[WATER] )
		                    + (Edata.theta[SOIL] * Edata.soil[SOIL_RHO]);
		Edata.heatCapacity();
	}
}

/**
 * @brief Subsurface Freezing:
 * -# Check whether the conditions
 *    Sum(th(x)) = 1 [sum of volumetric contents soil, ice, water, and void] \n
 *    and \n
 *    rho <= th(s) * rho(s) + ( 1 - th(s) ) * rho(i) = rhomax \n
 *    with th(s), rho(s) volumetric content and density of soil, rho element density
 *    and rhomax maximal possible element density are fulfilled;
 *    if not return with an error code (violation  of the consistency requirements).
 * -# Temperature and water check: If the  element temperature is equal to or above
 *    freezing temperature, or if the volumetric water content is less or equal to the
 *    residual water content, no subsurface  freezing takes place. \n
 *    If the temperature and water check is ok, then the difference dT between actual
 *    element and freezing temperature is evaluated.
 * -# The coefficient A is computed.
 * -# The changes d th(i) and d th(w) in the volumetric ice and water content are computed.
 * -# If there is  NOT  enough water to perform the computed d th(w) then
 *    d_th(w) ( = th(w) ) as well as the changes d_th(i) and dT are corrected.
 * -# Check whether the element is frozen solid and edit a message.  Computation of the
 *    the corrected changes and contents is then performed:
 *    -# volumetric contents ( soil, ice, water, air )
 *    -# dry density
 *    -# element density
 *    -# fusion power.
 *    -# water loss rate
 * -# Again the consistency requirement of (1) is checked.
 * @param *Edata
 * @param dt Time step (s)
 */
void PhaseChange::compSubSurfaceFrze(ElementData& Edata, const unsigned int nSolutes, const double& dt,
                                     const mio::Date& date_in)
{
	const double T_freeze=Edata.freezing_tk;	// Retrieve melting temperature from ElementData

	if(!Edata.checkVolContent()) prn_msg(__FILE__, __LINE__, "wrn", Date(), "wrong volumetric content");
	/*
	 * Freezing within the snowpack can occur if (1) the temperature of the element is below freezing
	 * and if water is present to be refrozen
	*/
	const double cmp_theta_r=((iwatertransportmodel_snow==RICHARDSEQUATION && Edata.theta[SOIL]<Constants::eps) || (iwatertransportmodel_soil==RICHARDSEQUATION && Edata.theta[SOIL]>Constants::eps)) ? (PhaseChange::RE_theta_r) : (PhaseChange::theta_r);
	if ((Edata.Te >= T_freeze) || (Edata.theta[WATER] <= cmp_theta_r)) {
		return;
	} else {
		double dT = T_freeze - Edata.Te;
		// Adapt A to compute mass changes
		double A = (Edata.c[TEMPERATURE] * Edata.Rho) / ( Constants::density_ice * Constants::lh_fusion );
		// Compute the change in volumetric ice and water contents
		double dth_i = A * dT;
		double dth_w = - (Constants::density_ice / Constants::density_water) * dth_i;
		// Make sure that there is enough water to refreeze
		if ((Edata.theta[WATER] + dth_w) < cmp_theta_r) {
			dth_w = - fabs( Edata.theta[WATER] - cmp_theta_r );
			dth_i = - (Constants::density_water / Constants::density_ice) * dth_w;
			dT = dth_i / A;
		}
		// See if the element is pure ICE
		if ((Edata.theta[ICE] + cmp_theta_r + dth_i + Edata.theta[SOIL]) >= 1.0) {
			dth_w = - fabs( Edata.theta[WATER] - cmp_theta_r );
			dth_i = - (Constants::density_water / Constants::density_ice) * dth_w;
			Edata.theta[WATER] = cmp_theta_r;
			Edata.theta[ICE] = 1.0 - Edata.theta[SOIL] - Edata.theta[WATER];
			Edata.theta[AIR] = 0.0;
		} else {
			// Concentration of solutes
			for (unsigned int ii = 0; ii < nSolutes; ii++) {
				if (dth_i > 0.) {
					Edata.conc[ICE][ii] = (Edata.theta[ICE] * Edata.conc[ICE][ii] +
					                           dth_i * Edata.conc[WATER][ii]) / ( Edata.theta[ICE] + dth_i);
				}
			}
			Edata.theta[ICE] += dth_i;
			Edata.theta[WATER] += dth_w;
			Edata.theta[AIR] = MAX(0., 1.0 - Edata.theta[ICE] - Edata.theta[WATER] - Edata.theta[SOIL]);
		}
		// State when the element is wet (PERMAFROST)
		if (Edata.theta[WATER] >= 1.0) {
			prn_msg(__FILE__, __LINE__, "msg+", Date(), "Wet Element! (dth_w=%e) (compSubSurfaceFrze)", dth_w);
			Edata.theta[WATER] = 1.0;
		}
		// Make sure the sum of all volumetric contents is near 1 (Can make a 1% error)
		if (!Edata.checkVolContent()) {
			prn_msg(__FILE__, __LINE__, "err", date_in, "Sum theta[I,W,A,S] > 1");
			prn_msg(__FILE__, __LINE__, "msg-", Date(),
			        "Ice: %f, Water: %f, Air: %f Soil: %f",
			        Edata.theta[ICE], Edata.theta[WATER], Edata.theta[AIR], Edata.theta[SOIL]);
			throw IOException("In compSubSurfaceFrze!", AT);
		}
		Edata.Rho = Constants::density_ice * Edata.theta[ICE] +
		                (Constants::density_water * Edata.theta[WATER]) +
		                    (Edata.theta[SOIL] * Edata.soil[SOIL_RHO]);
		Edata.heatCapacity();
		// Compute the volumetric refreezing power
		Edata.Qmf += (dth_i * Constants::density_ice * Constants::lh_fusion) / dt; // (W m-3)
		Edata.dth_w = dth_w;
		Edata.Te += dT;
	}
}




void PhaseChange::initialize(SnowStation& Xdata)
{
	// Initialize PhaseChange: execute this function before doing any call to PhaseChange::compPhaseChange for the current time step, to reset the energy balance values.
	size_t e, nE;
	ElementData* EMS;
	nE = Xdata.getNumberOfElements(); EMS = &Xdata.Edata[0];

	// Initialize and Determine Energy Content
	for (e = 0; e < nE; e++) {
		EMS[e].dth_w = EMS[e].Qmf = 0.;
	}

	// Get cold content
	cold_content_in=Xdata.ColdContent;

	// Reset meltFreezeEnergy and dIntEnergy
	Xdata.meltFreezeEnergy=0.;
	Xdata.dIntEnergy=0.;

	return;
}



void PhaseChange::finalize(const SurfaceFluxes& Sdata, SnowStation& Xdata, const mio::Date& date_in)
{
	// After all PhaseChange::compPhaseChange calls for the current time step, execute this function to finalize temperature structure calculations and energy balance values.
	size_t e, nE;
	double sum_Qmf=0.;
	cold_content_out=0.;

	ElementData* EMS;
	bool prn_CK = false;
	nE = Xdata.getNumberOfElements(); EMS = &Xdata.Edata[0]; vector<NodeData>& NDS = Xdata.Ndata;

	try {
		// In the final step compute temperature and temperature gradient, check both density and mass balance
		for (e = 0; e < nE; e++) {
			//Check Nodal temperatures
			double thresh_th_w;
			// In soils, some water may still be liquid below freezing
			if (EMS[e].theta[SOIL] > Constants::eps2) {
				thresh_th_w = PhaseChange::theta_r;
				thresh_th_w = (iwatertransportmodel_soil==RICHARDSEQUATION) ? 1. : (thresh_th_w);
			} else {
				thresh_th_w = 0.;
				thresh_th_w = (iwatertransportmodel_snow==RICHARDSEQUATION) ? (PhaseChange::RE_theta_r) : (thresh_th_w);
			}
			// Now, if you have water in the element -- the nodal temperatures must be 0 degC!!!
			// Element temperature will become zero (if it is not already) in the next step, where we calculate element temperature
			// as average over nodal temperatures.
			if ((EMS[e].theta[WATER] > thresh_th_w) &&
				((EMS[e].theta[ICE] > Snowpack::min_ice_content) || (EMS[e].Te <= EMS[e].melting_tk))) {
				NDS[e+1].T = EMS[e].melting_tk;
				if (e > 0) // NOTE Bottom soil node temperature cannot be changed
					NDS[e].T = EMS[e].melting_tk;
			}

			//Restructure temperature arrays
			EMS[e].gradT = (NDS[e+1].T - NDS[e].T) / EMS[e].L;
			EMS[e].Te = (NDS[e].T + NDS[e+1].T) / 2.0;
			if (((EMS[e].Te - EMS[e].melting_tk) > 0.2) && (e > Xdata.SoilNode))
				prn_msg(__FILE__, __LINE__, "wrn", date_in,
				        "Snow temperature Te=%f K is above melting point in element %d (nE=%d; T0=%f K, T1=%f K)",
				        EMS[e].Te, e, nE, NDS[e].T, NDS[e+1].T);
			if (EMS[e].theta[SOIL] < Constants::eps2) {
				if (!(EMS[e].Rho > Constants::eps && EMS[e].Rho <= Constants::max_rho)) {
					prn_msg(__FILE__, __LINE__, "err", date_in, "Phase Change End: rho_snow[%d]=%f", e, EMS[e].Rho);
					throw IOException("Run-time error in compPhaseChange()", AT);
				}
			}
			cold_content_out += EMS[e].c[TEMPERATURE] * EMS[e].Rho * (EMS[e].Te - EMS[e].melting_tk) * EMS[e].L;
			sum_Qmf += EMS[e].Qmf * EMS[e].L;
		}
		if (prn_CK && (sum_Qmf > 0.)) {
			prn_msg(__FILE__, __LINE__, "msg+", date_in, "Checking energy balance  (W/m2):");
			prn_msg(__FILE__, __LINE__, "msg+", date_in, " E1: %f   E0: %f  E1-E0: %f  sum_Qmf: %f  Surface EB : %f",
			           (cold_content_out) / sn_dt, (cold_content_in) / sn_dt,
			               (cold_content_out - cold_content_in) / sn_dt, sum_Qmf,
			                   Sdata.qs + Sdata.ql + Sdata.lw_net + Sdata.qr + Sdata.qw);
		}
	} catch (const exception& ) {
		throw;
	}

	return;
}



/**
 * @brief Driving routine for subsurface melting and refreezing as well as surface melting.
 * The basic equation in both subsurface processes is: \n
 * d_th(i) = A * dt \n
 * with th(i) volumetric ice content (1), c_p(T) specific heat capacity of ice (J kg-1 K-1),
 * Q_f the freezing / melting energy (J kg-1), T the absolute temperature (K),
 * and the coefficient: \n
 * A = c_p(T) * th(i) * Q_f \n
 * ql_Rest is the Energy that is transferred from the upper element to the lower one (J m-2)
 * in case of complete melting of the former
 * @param Xdata
 * @param date_in is the current date
 * @param verbose print detailed warnings for various situations? (default=true)
 */
void PhaseChange::compPhaseChange(SnowStation& Xdata, const mio::Date& date_in, const bool& verbose)
{
	size_t e, nE;
	double ql_Rest;
	ElementData* EMS;
	nE = Xdata.getNumberOfElements(); EMS = &Xdata.Edata[0]; vector<NodeData>& NDS = Xdata.Ndata;

	try {
		// In the first step:
		// 1) check the density
		// 2) set flags for SubSurfaceMelt and SubSurfaceFrze

		//Reset flags
		Xdata.SubSurfaceMelt = false;
		Xdata.SubSurfaceFrze = false;

		for (e = 0; e < nE; e++) {
			if (EMS[e].theta[SOIL] == 0.0) {
				if (verbose && !(EMS[e].Rho > 0. && EMS[e].Rho <= Constants::max_rho)) {
					prn_msg(__FILE__, __LINE__, "wrn", date_in, "Phase Change Begin: rho[%d]=%f", e, EMS[e].Rho);
				}
			}
			// and make sure the sum of all volumetric contents is near 1 (Can make a 1% error)
			if (verbose && !EMS[e].checkVolContent()) {
				prn_msg(__FILE__, __LINE__, "msg+", date_in,
				        "Phase Change Begin: Element=%d, nE=%d  ICE %f, Water %f, Air %f Soil %f",
				        e, nE, EMS[e].theta[ICE], EMS[e].theta[WATER], EMS[e].theta[AIR], EMS[e].theta[SOIL]);
			}

			const double cmp_theta_r=((iwatertransportmodel_snow==RICHARDSEQUATION && EMS[e].theta[SOIL]<Constants::eps) || (iwatertransportmodel_soil==RICHARDSEQUATION && EMS[e].theta[SOIL]>Constants::eps)) ? (PhaseChange::RE_theta_r) : (PhaseChange::theta_r);
			// Set flags
			if (EMS[e].Te > EMS[e].melting_tk && EMS[e].theta[ICE] > Constants::eps2)
				Xdata.SubSurfaceMelt = true;
			if (EMS[e].Te < EMS[e].freezing_tk && EMS[e].theta[WATER] > cmp_theta_r + Constants::eps2)
				Xdata.SubSurfaceFrze = true;
		}

		if (Xdata.SubSurfaceMelt) {
			ql_Rest = 0.;
			e = nE;
			while (e > 0) {
				e--;
				const double i_Te = EMS[e].Te;
				try {
					if(!(iwatertransportmodel_soil==RICHARDSEQUATION && e<Xdata.SoilNode)) {
						compSubSurfaceMelt(EMS[e], Xdata.number_of_solutes, sn_dt, ql_Rest, date_in);
					}
				} catch(...) {
					prn_msg(__FILE__, __LINE__, "msg-", Date(), "in compSubSurfaceMelt at element %d of %d", e, nE);
					throw IOException("Run-time error in compPhaseChange()", AT);
				}
				// See whether hoar has melted
				if (EMS[e].theta[WATER] > Constants::eps) {
					NDS[e].hoar = 0.;
					NDS[e+1].hoar = 0.;
				}
				// Make sure all nodal temperatures are consistent with the temperature change of the element
				if( (e < Xdata.SoilNode && iwatertransportmodel_soil!=RICHARDSEQUATION) || (e >= Xdata.SoilNode) ) {
					NDS[e].T += (EMS[e].Te - i_Te);
					NDS[e+1].T += (EMS[e].Te - i_Te);
					if(EMS[e].theta[ICE]>0.) {	// If ice is present, nodal temperature cannot exceed melting_tk
						NDS[e].T = MIN(NDS[e].T, EMS[e].melting_tk);
						NDS[e+1].T = MIN(NDS[e+1].T, EMS[e].melting_tk);
					}
				}
			}
		}
		// TODO If WATER_LAYER && ql_rest > 0, consider evaporating water left in the last element above soil!

		// Subsurface refreezing check
		if (Xdata.SubSurfaceFrze) {
			for (e = 0; e < nE; e++) {
				const double i_Te = EMS[e].Te;
				try {
					if(!(iwatertransportmodel_soil==RICHARDSEQUATION && e<Xdata.SoilNode)) {	// For Richards Equation, phase changes in soil are taken care of in WaterTransport
						compSubSurfaceFrze(EMS[e], Xdata.number_of_solutes, sn_dt, date_in);
					}
				} catch(...) {
					prn_msg(__FILE__, __LINE__, "msg-", date_in, "SubSurfaceFrze at element %d of %d", e, nE);
					throw IOException("Run-time error in compPhaseChange()", AT);
				}
				// Make sure all nodal temperatures are consistent with the temperature change of the element
				if( (e < Xdata.SoilNode && iwatertransportmodel_soil!=RICHARDSEQUATION) || (e >= Xdata.SoilNode) ) {
					NDS[e].T += (EMS[e].Te - i_Te);
					NDS[e+1].T += (EMS[e].Te - i_Te);
					if(EMS[e].theta[ICE]>0.) {	// If ice is present, nodal temperature cannot exceed melting_tk
  						NDS[e].T = MIN(NDS[e].T, EMS[e].melting_tk);
						NDS[e+1].T = MIN(NDS[e+1].T, EMS[e].melting_tk);
					}
				}
			}
		}
	} catch (const exception& ) {
		throw;
	}
}
