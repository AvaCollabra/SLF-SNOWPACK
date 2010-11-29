/*
 *  SNOWPACK VERSION 9.x
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

//Saturated Water Content, for now we say  1.0
const double PhaseChange::theta_s = 1.0; 

/************************************************************
 * non-static section                                       *
 ************************************************************/

PhaseChange::PhaseChange(const mio::Config& i_cfg) : cfg(i_cfg) 
{
	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Parameters");
	sn_dt = M_TO_S(calculation_step_length);
}

/**
 * @brief Subsurface Melting:
 * -# Consistency check.
 * -# Temperature ice and water check: If the  element temperature is below melting
 * temperature, or if the volumetric ice content is nonpositive, or if the water content
 * equals or exceeds PhaseChange::theta_s, no melting occurs.
 * If the check is ok,  then the  difference dT between element and melting temperature
 * is calculated.
 * -# The coefficient A is calculated.
 * -# Analogous to Subsurface Freezing ...
 * -# Check if there is enough ice to be melted: If it is not the case, the values have to
 * be corrected.
 * -# If there is  TOO MUCH  water being melted  then d_th(w) ( = PhaseChange::theta_s - th(w) ) as welL
 * as the changes d_th(i) and dT are corrected.
 * -# Calculation of contents.
 * -# Characterization of element ( Liquid, Dry, Void ).
 * -# Calculation of densities and water generation rate.
 * @param *Edata
 * @param dt Time step (s)
 * @param *ql_Rest Latent heat flux balance (W m-2)
 */
void PhaseChange::calcSubSurfaceMelt(SN_ELEM_DATA& Edata, const double& dt, double& ql_Rest)
{
	double sum;   // sum of components
	double dT;    // Constants::freezing_tk - Te > 0
	double A;     // coefficient A ( see notes below )
	double dth_i; // change in volumetric ice content
	double dth_w; // change in volumetric water content
	int i;

	// First make sure the sum of all volumetric contents is near 1 (Can make a 1% error)
	for (sum = 0.0, i = 0; i < N_COMPONENTS; i++) {
		sum += Edata.theta[i];
	}
	if (sum <= 0.99 || sum >= 1.01) {
		prn_msg(__FILE__, __LINE__, "wrn", -1., "SUM of volumetric contents = %1.4f", sum);
	}
	// The volumetric contents must remain positive ....
	if ( (Edata.theta[SOIL] < -0.0001) || (Edata.theta[ICE] < -0.01) || (Edata.theta[WATER] < -0.01) || (Edata.theta[AIR] < -0.01) ) {
		prn_msg(__FILE__, __LINE__, "wrn", -1., "Negative volumetric content! [I/W/A/S] [%e/%e/%e/%e]", Edata.theta[ICE], Edata.theta[WATER], Edata.theta[AIR], Edata.theta[SOIL]);
	}
	/*
	 * Now see if any melting is going on -- this implies that (1) the temperature of the element
	 * is above the melting temperature (2) there is something to melt and (3) there is enough
	 * room to place the meltwater ...
	*/
	if (((Edata.Te < Constants::melting_tk) && (ql_Rest == 0.)) 
		|| (Edata.theta[ICE] <= 0.0) || (Edata.theta[WATER] >= PhaseChange::theta_s)) {
		return;
	} else {
		dT = Constants::melting_tk - Edata.Te;
		// Now we take into account that there might be some extra energy that could not be used 
		// by the element above because of complete melting
		dT -= ql_Rest / (Edata.c[TEMPERATURE] * Edata.Rho * Edata.L);
		if (dT > 0.) {
			return;
		}
		// Determine the DECREASE in ice content and the INCREASE of water content
		A = (Edata.c[TEMPERATURE] * Edata.Rho) / ( Constants::density_ice * Constants::lh_fusion ); // Adapt A to calculate mass changes
		dth_i = A * dT;
		dth_w = - (Constants::density_ice / Constants::density_water) * dth_i;
		// It could happen that there is enough energy available to melt more ice than is present. You can only melt so much ice as is there ....
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
		// Calculate the new chemical concentrations
		for (i = 0; i < N_SOLUTES; i++) {
			if( dth_w > 0. ) {
				Edata.conc[WATER][i] = (Edata.theta[WATER] * Edata.conc[WATER][i] + dth_w * Edata.conc[ICE][i]) / (Edata.theta[WATER] + dth_w);
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
		// Calculate the bulk density
		Edata.Rho = Constants::density_ice * Edata.theta[ICE] + (Constants::density_water * Edata.theta[WATER] ) + (Edata.theta[SOIL] * Edata.soil[SOIL_RHO]);
		if ( (Edata.theta[SOIL] == 0.0) && !(Edata.Rho > 0 && Edata.Rho <= Constants::max_rho) ) {
			prn_msg(__FILE__, __LINE__, "err", -1., "Rho(snow)=%lf (SubSurfaceMelt)", Edata.Rho);
			throw IOException("Error in calcSubSurfaceMelt()", AT);
		}
	 	// Re-Calculate the specific heat capacity of our porous medium
		Edata.c[TEMPERATURE] = SnLaws::calcHeatCapacity(Edata);
		// Calculate MELTING ENERGY
		Edata.Qmf = (dth_i * Constants::density_ice * Constants::lh_fusion) / dt; // (W m-3)
		Edata.dth_w = dth_w;                                // (1)
		// Reset the element temperature
		Edata.Te += dT;
		if ( Edata.Te <= Constants::melting_tk ) {
			ql_Rest = 0.0;
		} else {
			ql_Rest = Edata.c[TEMPERATURE] * Edata.Rho * Edata.L * (Edata.Te - Constants::melting_tk);
			Edata.Te = Constants::melting_tk;
		}
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
 * -# The coefficient A is calculated.
 * -# The changes d th(i) and d th(w) in the volumetric ice and water content are calculated.
 * -# If there is  NOT  enough water to perform the calculated d th(w) then
 *    d_th(w) ( = th(w) ) as well as the changes d_th(i) and dT are corrected.
 * -# Check whether the element is frozen solid and edit a message.  Calculation of the
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
void PhaseChange::calcSubSurfaceFrze(SN_ELEM_DATA& Edata, const double& dt)
{
	double sum;   // sum of volumetric components
	double dT;    // Constants::freezing_tk - Te > 0
	double A;     // coefficient A ( see notes below )
	double dth_i; // change in volumetric ice content
	double dth_w; // change in volumetric water content
	double hsum;  // auxiliary variable
	int i;

	// First make sure the sum of all volumetric contents is near 1 (Can make a 1% error)
	for (sum = 0.0, i = 0; i < N_COMPONENTS; i++) {
		sum += Edata.theta[i];
	}
	if ( sum <= 0.99 || sum >= 1.01 ) {
		prn_msg(__FILE__, __LINE__, "wrn", -1., "SUM of volumetric contents = %1.4f", sum);
	}
	// The volumetric contents must remain positive ....
	if ( (Edata.theta[SOIL]  < -0.0001) || (Edata.theta[ICE]   < -0.01) || (Edata.theta[WATER] < -0.01) || (Edata.theta[AIR]   < -0.01) ) {
		prn_msg(__FILE__, __LINE__, "wrn", -1., "Negative volumetric content! [I/W/A/S] [%e/%e/%e/%e]", Edata.theta[ICE], Edata.theta[WATER], Edata.theta[AIR], Edata.theta[SOIL]);
	}
	/*
	 * Freezing within the snowpack can occur if (1) the temperature of the element is below 0
	 * and if water is present to be refrozen
	*/
	if ( (Edata.Te >= Constants::freezing_tk) || (Edata.theta[WATER] <= PhaseChange::theta_r) ) {
		return;
	} else {
		dT = Constants::freezing_tk - Edata.Te;
		// Adapt A to calculate mass changes
		A = (Edata.c[TEMPERATURE] * Edata.Rho) / ( Constants::density_ice * Constants::lh_fusion );
		// Calculate the change in volumetric ice and water contents
		dth_i = A * dT;
		dth_w = - (Constants::density_ice / Constants::density_water) * dth_i;
		// Make sure that there is enough water to refreeze
		if ( (Edata.theta[WATER] + dth_w) < PhaseChange::theta_r ) {
			dth_w = - fabs( Edata.theta[WATER] - PhaseChange::theta_r );
			dth_i = - (Constants::density_water / Constants::density_ice) * dth_w;
			dT = dth_i / A;
		}
		// See if the element is pure ICE
		hsum =  Edata.theta[ICE] + PhaseChange::theta_r + dth_i + Edata.theta[SOIL];
		if ( hsum >= 1.0 ) {
			dth_w = - fabs( Edata.theta[WATER] - PhaseChange::theta_r );
			dth_i = - (Constants::density_water / Constants::density_ice) * dth_w;
			Edata.theta[ICE] = 1.0 - Edata.theta[SOIL];
			Edata.theta[WATER] = PhaseChange::theta_r;
			Edata.theta[AIR] = 0.0;
		} else {
			// If not calculate the volumetric components
			for (i = 0; i < N_SOLUTES; i++) {
				if (dth_i > 0.) {
					Edata.conc[ICE][i] = ( Edata.theta[ICE] * Edata.conc[ICE][i] + dth_i * Edata.conc[WATER][i]) / ( Edata.theta[ICE] + dth_i);
				}
			}
			Edata.theta[ICE] += dth_i;
			Edata.theta[WATER] += dth_w;
			Edata.theta[AIR] = MAX(0., 1.0 - Edata.theta[ICE] - Edata.theta[WATER] - Edata.theta[SOIL]);
		}
		// State when the water content has disappeared (PERMAFROST)
		if ( Edata.theta[WATER] >= 1.0 ) {
			prn_msg(__FILE__, __LINE__, "msg+", -1., "Wet Element! (dth_w=%e)", dth_w);
			Edata.theta[WATER] = 1.0;
		}
		// Calculate the densities
		Edata.Rho = Constants::density_ice * Edata.theta[ICE] + (Constants::density_water * Edata.theta[WATER]) + (Edata.theta[SOIL] * Edata.soil[SOIL_RHO]);
		if( (Edata.theta[SOIL] == 0.0) && !(Edata.Rho > 0. && Edata.Rho <= Constants::max_rho) ) {
			prn_msg(__FILE__, __LINE__, "err", -1., "Rho(snow)=%lf (SubSurfaceFrze)", Edata.Rho);
			throw IOException("Error in calcSubSurfaceFrze()", AT);
		}
		// Re-Calculate the specific heat capacity of our porous medium
		Edata.c[TEMPERATURE] = SnLaws::calcHeatCapacity(Edata);
		// Calculate the REFREEZING Energy
		Edata.Qmf = dth_i * Constants::density_ice * Constants::lh_fusion / dt; // (W m-3)
		Edata.dth_w = dth_w;                              // (1)
		// Reset the element temperature
		Edata.Te += dT;
	}
}

/**
 * @brief Driving routine for  subsurface melting and  refreezing as well as surface melting.
 * The basic equation in both subsurface processes is: \n
 * d_th(i) = A * dt \n
 * with th(i) volumetric ice content (1), c_p(T) specific heat capacity of ice (J kg-1 K-1), \n
 * Q_f the freezing / melting energy (J kg-1), T the absolute temperature (K ), \n
 * and the coefficient: \n
 * A = c_p(T) * th(i) * Q_f
 * @param *Xdata
 */
void PhaseChange::runPhaseChange(const SN_SURFACE_DATA& Sdata, SN_STATION_DATA& Xdata)
{
	int e, i, nE;
	double ql_Rest, sum, Te_old; // Energy that is transferred from the upper element to the lower one in case of complete melting of the former
	double ColdContent0=Xdata.ColdContent, ColdContent1=0., S_Qmf=0., M0=0., M1=0.;
	SN_ELEM_DATA* EMS;

	nE = Xdata.getNumberOfElements(); EMS = &Xdata.Edata[0];  vector<SN_NODE_DATA>& NDS = Xdata.Ndata;
	// Initialize and Determine Energy Content
	for (e = 0; e < nE; e++) {
		EMS[e].dth_w = EMS[e].Qmf = 0.0;
		EMS[e].c[TEMPERATURE] = SnLaws::calcHeatCapacity(EMS[e]);
		M0 += EMS[e].theta[ICE] * EMS[e].L;
	}

	try {
		// In the first step check the density
		for (e = 0; e < nE; e++) {
			if ( EMS[e].theta[SOIL] == 0.0) {
				if ( !(EMS[e].Rho > 0. && EMS[e].Rho <= Constants::max_rho) ) {
					prn_msg(__FILE__, __LINE__, "wrn", -1., "Phase Change Begin: rho[%d]=%lf", e, EMS[e].Rho);
				}
			}
			// and make sure the sum of all volumetric contents is near 1 (Can make a 1% error)
			for (sum = 0.0, i = 0; i < N_COMPONENTS; i++) {
				sum += EMS[e].theta[i];
			}
			if ( sum <= 0.99 || sum >= 1.01 ) {
				prn_msg(__FILE__, __LINE__, "msg+", -1., "Phase Change Begin: sum of volumetric contents = %7.4f", sum);
				prn_msg(__FILE__, __LINE__, "msg", -1., "Element=%d, nE=%d  ICE %lf, Water %lf, Air %lf Soil %lf", 
					   e, nE, EMS[e].theta[ICE], EMS[e].theta[WATER], EMS[e].theta[AIR], EMS[e].theta[SOIL]);
			}
		}
		
		if ( Xdata.SubSurfaceMelt ) {
			double thresh_th_w;
			for (e = nE-1, ql_Rest = 0.; e >= 0; e--) {
				try {
					calcSubSurfaceMelt(EMS[e], sn_dt, ql_Rest);
				} catch(...) {
					prn_msg(__FILE__, __LINE__, "err", -1., "SubSurfaceMelt at element [%d], nE=%d", e, nE);
					throw;
				}

				// In soils, some water may still be liquid below freezing
				if ( EMS[e].theta[SOIL] > 0. ) {
					thresh_th_w = PhaseChange::theta_r;
				} else {
					thresh_th_w = 0.0;
				}
				// Now, if you have water in the element -- the nodal temperatures must be 0!!!
				if ( (EMS[e].theta[WATER] > thresh_th_w) && 
					((EMS[e].theta[ICE] > Constants::min_ice_content) || (EMS[e].Te <= Constants::melting_tk)) ) {
					NDS[e+1].T = Constants::melting_tk;
					// Bottom soil node temperature cannot be changed
					if ( e > 0 || !(EMS[e].theta[SOIL] > 0.) ) {
						NDS[e].T = Constants::melting_tk;
					}
				}
				// See whether hoar has melted
				if ( EMS[e].theta[WATER] > 0.0 ) {
					NDS[e].hoar = 0.0;
					NDS[e+1].hoar = 0.0;
				}
			}
		}
		// TODO If WATER_LAYER && ql_rest > 0, consider evaporating water left in the last element above soil!
		
		// Subsurface refreezing check
		if ( Xdata.SubSurfaceFrze ) {
			double thresh_th_w;
			for (e = 0; e < nE; e++) {
				Te_old = EMS[e].Te;

				try {
					calcSubSurfaceFrze(EMS[e], sn_dt);
				} catch(...) {
					prn_msg(__FILE__, __LINE__, "err", -1., "SubSurfaceFrze at element [%d], nE=%d", e, nE);
					throw;
				}

				// In soils, some water may still be liquid below freezing
				if ( EMS[e].theta[SOIL] > 0. ) {
					thresh_th_w = PhaseChange::theta_r;
				} else {
					thresh_th_w = 0.0;
				}
				// Now, if you have water in the element -- the nodal temperatures must be 0!!!
				if ( (EMS[e].theta[WATER] > thresh_th_w) && (EMS[e].theta[ICE] > Constants::min_ice_content) ) {
					NDS[e+1].T = Constants::melting_tk;
					// Bottom soil node temperature cannot be changed
					if ( e > 0 || !(EMS[e].theta[SOIL] > 0.) )
						NDS[e].T = Constants::melting_tk;
				}
				// Make sure all nodal temperatures are consistent
				if (EMS[e].Te < Te_old) {
					NDS[e].T += 0.5 * (EMS[e].Te - Te_old);
					NDS[e+1].T += 0.5 * (EMS[e].Te - Te_old);
				}
			}
		}

		// In the final step compute temperature and temperature gradient, check both density and mass balance
		for (e = 0; e < nE; e++) {
			EMS[e].gradT = (NDS[e+1].T - NDS[e].T) / EMS[e].L;
			EMS[e].Te = (NDS[e].T + NDS[e+1].T) / 2.0;
			if ( EMS[e].theta[SOIL] == 0.0 ) {
				if ( !(EMS[e].Rho > 0. && EMS[e].Rho <= Constants::max_rho) ) {
					prn_msg(__FILE__, __LINE__, "err", -1., "Phase Change End: rho[%d]=%lf", e, EMS[e].Rho);
					throw IOException("Runtime error in runPhaseChange()", AT);
				}
			}
			// Also make sure the sum of all volumetric contents is near 1 (Can make a 1% error)
			for (sum = 0.0, i = 0; i < N_COMPONENTS; i++) {
				sum += EMS[e].theta[i];
			}
			if ( sum <= 0.99 || sum >= 1.01 ) {
				prn_msg(__FILE__, __LINE__, "err", -1., "Phase Change End: sum of volumetric contents = %1.4f", sum);
				prn_msg(__FILE__, __LINE__, "msg", -1., "Element=%d, nE=%d  ICE %lf, Water %lf, Air %lf Soil %lf", 
					   e, nE, EMS[e].theta[ICE], EMS[e].theta[WATER], EMS[e].theta[AIR], EMS[e].theta[SOIL]);
				throw IOException("Runtime error in runPhaseChange()", AT);
		}
			M1 += EMS[e].theta[ICE] * EMS[e].L;
			ColdContent1 += EMS[e].c[TEMPERATURE] * EMS[e].Rho * (EMS[e].Te - Constants::melting_tk) * EMS[e].L;
			S_Qmf += EMS[e].Qmf * EMS[e].L;
		}
		if ( 0 && (S_Qmf > Constants::eps) ) {
			prn_msg(__FILE__, __LINE__, "msg+", -1., "Checking energy balance  (W/m2):");
			prn_msg(__FILE__, __LINE__, "msg+", -1., " E1: %lf   E0: %lf  E1-E0: %lf  S_Qmf: %lf  Surface EB : %lf", 
				   (ColdContent1) / sn_dt, (ColdContent0) / sn_dt, 
				   (ColdContent1 - ColdContent0) / sn_dt, S_Qmf, 
				   Sdata.qs + Sdata.ql + Sdata.lw_net + Sdata.qr + Sdata.qw);
		}
	} catch (exception& ex) {
		prn_msg(__FILE__, __LINE__, "err", -1, "Runtime error in runPhaseChange()");
		throw;
	}
}


/*
 * END of Phase Change
 */
