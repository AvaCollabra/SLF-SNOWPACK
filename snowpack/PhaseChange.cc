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

//Saturated Water Content, for now we say  1.0
const double PhaseChange::theta_s = 1.0;

/************************************************************
 * non-static section                                       *
 ************************************************************/

PhaseChange::PhaseChange(const mio::Config& cfg)
{
	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	sn_dt = M_TO_S(calculation_step_length);
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
 * @param *ql_Rest Latent heat flux balance (W m-2)
 */
void PhaseChange::compSubSurfaceMelt(ElementData& Edata, const unsigned int nSolutes, const double& dt, double& ql_Rest)
{
	double dT;    // Constants::freezing_tk - Te > 0
	double A;     // coefficient A ( see notes below )
	double dth_i; // change in volumetric ice content
	double dth_w; // change in volumetric water content

	Edata.checkVolContent();
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
		// Now we take into account that there might be some extra energy that could not
		// be used by the element above because of complete melting
		dT -= ql_Rest / (Edata.c[TEMPERATURE] * Edata.Rho * Edata.L);
		if (dT > 0.) {
			return;
		}
		// Determine the DECREASE in ice content and the INCREASE of water content
		A = (Edata.c[TEMPERATURE] * Edata.Rho) / ( Constants::density_ice * Constants::lh_fusion ); // Adapt A to compute mass changes
		dth_i = A * dT;
		dth_w = - (Constants::density_ice / Constants::density_water) * dth_i;
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
		if (Edata.Te <= Constants::melting_tk) {
			ql_Rest = 0.0;
			Edata.Te = Constants::melting_tk;
		} else {
			ql_Rest = Edata.c[TEMPERATURE] * Edata.Rho * Edata.L * (Edata.Te - Constants::melting_tk);
			Edata.Te = Constants::melting_tk;
		}
		Edata.Qmf = (dth_i * Constants::density_ice * Constants::lh_fusion) / dt; // (W m-3)
		Edata.dth_w = dth_w;                                // (1)
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
		Edata.Rho = Constants::density_ice * Edata.theta[ICE]
		                + (Constants::density_water * Edata.theta[WATER] )
		                    + (Edata.theta[SOIL] * Edata.soil[SOIL_RHO]);
		if ((Edata.theta[SOIL] == 0.0) && !(Edata.Rho > 0. && Edata.Rho <= Constants::max_rho)) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "Rho(snow)=%f (SubSurfaceMelt)", Edata.Rho);
			throw IOException("Error in compSubSurfaceMelt()", AT);
		}
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
void PhaseChange::compSubSurfaceFrze(ElementData& Edata, const unsigned int nSolutes, const double& dt)
{
	Edata.checkVolContent();
	/*
	 * Freezing within the snowpack can occur if (1) the temperature of the element is below 0
	 * and if water is present to be refrozen
	*/
	if ((Edata.Te >= Constants::freezing_tk) || (Edata.theta[WATER] <= PhaseChange::theta_r)) {
		return;
	} else {
		double dT = Constants::freezing_tk - Edata.Te;
		// Adapt A to compute mass changes
		double A = (Edata.c[TEMPERATURE] * Edata.Rho) / ( Constants::density_ice * Constants::lh_fusion );
		// Compute the change in volumetric ice and water contents
		double dth_i = A * dT;
		double dth_w = - (Constants::density_ice / Constants::density_water) * dth_i;
		// Make sure that there is enough water to refreeze
		if ((Edata.theta[WATER] + dth_w) < PhaseChange::theta_r) {
			dth_w = - fabs( Edata.theta[WATER] - PhaseChange::theta_r );
			dth_i = - (Constants::density_water / Constants::density_ice) * dth_w;
			dT = dth_i / A;
		}
		// See if the element is pure ICE
		if ((Edata.theta[ICE] + PhaseChange::theta_r + dth_i + Edata.theta[SOIL]) >= 1.0) {
			dth_w = - fabs( Edata.theta[WATER] - PhaseChange::theta_r );
			dth_i = - (Constants::density_water / Constants::density_ice) * dth_w;
			Edata.theta[ICE] = 1.0 - Edata.theta[SOIL];
			Edata.theta[WATER] = PhaseChange::theta_r;
			Edata.theta[AIR] = 0.0;
		} else {
			// If the new chemical concentrations
			for (unsigned int ii = 0; ii < nSolutes; ii++) {
				if (dth_i > 0.) {
					Edata.conc[ICE][ii] = (Edata.theta[ICE] * Edata.conc[ICE][ii] + dth_i * Edata.conc[WATER][ii]) / ( Edata.theta[ICE] + dth_i);
				}
			}
			Edata.theta[ICE] += dth_i;
			Edata.theta[WATER] += dth_w;
			Edata.theta[AIR] = MAX(0., 1.0 - Edata.theta[ICE] - Edata.theta[WATER] - Edata.theta[SOIL]);
		}
		// State when the water content has disappeared (PERMAFROST)
		if (Edata.theta[WATER] >= 1.0) {
			prn_msg(__FILE__, __LINE__, "msg+", Date(), "Wet Element! (dth_w=%e)", dth_w);
			Edata.theta[WATER] = 1.0;
		}
		Edata.Rho = Constants::density_ice * Edata.theta[ICE] + (Constants::density_water * Edata.theta[WATER]) + (Edata.theta[SOIL] * Edata.soil[SOIL_RHO]);
		if ((Edata.theta[SOIL] == 0.0) && !(Edata.Rho > 0. && Edata.Rho <= Constants::max_rho)) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "Rho(snow)=%f (SubSurfaceFrze)", Edata.Rho);
			throw IOException("Error in compSubSurfaceFrze()", AT);
		}
		Edata.heatCapacity();
		// Compute the volumetric refreezing power
		Edata.Qmf = dth_i * Constants::density_ice * Constants::lh_fusion / dt; // (W m-3)
		Edata.dth_w = dth_w;
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
 * A = c_p(T) * th(i) * Q_f \n
 * ql_Rest is the Energy that is transferred from the upper element to the lower one
 * in case of complete melting of the former
 * @param Sdata
 * @param Xdata
 */
void PhaseChange::compPhaseChange(const SurfaceFluxes& Sdata, SnowStation& Xdata)
{
	unsigned int e, nE;
	double ql_Rest;
	double cold_content_in=Xdata.ColdContent, cold_content_out=0., sum_Qmf=0.;
	ElementData* EMS;
	bool prn_CK = false;

	nE = Xdata.getNumberOfElements(); EMS = &Xdata.Edata[0];  vector<NodeData>& NDS = Xdata.Ndata;
	// Initialize and Determine Energy Content
	for (e = 0; e < nE; e++) {
		EMS[e].dth_w = EMS[e].Qmf = 0.0;
	}

	try {
		// In the first step check the density
		for (e = 0; e < nE; e++) {
			if (EMS[e].theta[SOIL] == 0.0) {
				if (!(EMS[e].Rho > 0. && EMS[e].Rho <= Constants::max_rho)) {
					prn_msg(__FILE__, __LINE__, "wrn", Date(), "Phase Change Begin: rho[%d]=%f", e, EMS[e].Rho);
				}
			}
			// and make sure the sum of all volumetric contents is near 1 (Can make a 1% error)
			if (!EMS[e].checkVolContent()) {
				prn_msg(__FILE__, __LINE__, "msg+", Date(),
				        "Phase Change Begin: Element=%d, nE=%d  ICE %f, Water %f, Air %f Soil %f",
				        e, nE, EMS[e].theta[ICE], EMS[e].theta[WATER], EMS[e].theta[AIR], EMS[e].theta[SOIL]);
			}
		}

		if (Xdata.SubSurfaceMelt) {
			double thresh_th_w;
			e = nE;
			ql_Rest = 0.;
			//for (e = nE-1, ql_Rest = 0.; e >= 0; e--) {
			while (e-- > 0) {
				try {
					compSubSurfaceMelt(EMS[e], Xdata.number_of_solutes, sn_dt, ql_Rest);
				} catch(...) {
					prn_msg(__FILE__, __LINE__, "err", Date(), "SubSurfaceMelt at element [%d], nE=%d", e, nE);
					throw;
				}

				// In soils, some water may still be liquid below freezing
				if (EMS[e].theta[SOIL] > 0.) {
					thresh_th_w = PhaseChange::theta_r;
				} else {
					thresh_th_w = 0.0;
				}
				// Now, if you have water in the element -- the nodal temperatures must be 0!!!
				if ((EMS[e].theta[WATER] > thresh_th_w) &&
					((EMS[e].theta[ICE] > Snowpack::min_ice_content) || (EMS[e].Te <= Constants::melting_tk))) {
					NDS[e+1].T = Constants::melting_tk;
					// Bottom soil node temperature cannot be changed
					if (e > 0 || !(EMS[e].theta[SOIL] > 0.)) {
						NDS[e].T = Constants::melting_tk;
					}
				}
// 				if (NDS[e].T > Constants::melting_tk || NDS[e+1].T > Constants::melting_tk)
// 					prn_msg(__FILE__, __LINE__, "wrn", Date(),
// 					        "Melt: Snow temperature above melting (%f) in element %d (nE=%d) %f %f",
// 					        EMS[e].Te, e, nE, NDS[e].T, NDS[e+1].T);
				// See whether hoar has melted
				if (EMS[e].theta[WATER] > 0.0) {
					NDS[e].hoar = 0.0;
					NDS[e+1].hoar = 0.0;
				}
			}
		}
		// TODO If WATER_LAYER && ql_rest > 0, consider evaporating water left in the last element above soil!

		// Subsurface refreezing check
		if (Xdata.SubSurfaceFrze) {
			double thresh_th_w;
			for (e = 0; e < nE; e++) {
				double Te_old = EMS[e].Te;

				try {
					compSubSurfaceFrze(EMS[e], Xdata.number_of_solutes, sn_dt);
				} catch(...) {
					prn_msg(__FILE__, __LINE__, "err", Date(), "SubSurfaceFrze at element [%d], nE=%d", e, nE);
					throw;
				}

				// In soils, some water may still be liquid below freezing
				if (EMS[e].theta[SOIL] > 0.) {
					thresh_th_w = PhaseChange::theta_r;
				} else {
					thresh_th_w = 0.0;
				}
				// Now, if you have water in the element -- the nodal temperatures must be 0!!!
				if ((EMS[e].theta[WATER] > thresh_th_w) && (EMS[e].theta[ICE] > Snowpack::min_ice_content)) {
					NDS[e+1].T = Constants::melting_tk;
					// Bottom soil node temperature cannot be changed
					if (e > 0 || !(EMS[e].theta[SOIL] > 0.))
						NDS[e].T = Constants::melting_tk;
				}
				// Make sure all nodal temperatures are consistent
				if (EMS[e].Te < Te_old) {
					NDS[e].T += 0.5 * (EMS[e].Te - Te_old);
					NDS[e+1].T += 0.5 * (EMS[e].Te - Te_old);
				}
// 				if (NDS[e].T > Constants::melting_tk || NDS[e+1].T > Constants::melting_tk)
// 					prn_msg(__FILE__, __LINE__, "wrn", Date(),
// 					        "Freeze: Snow temperature above freezing (%f) in element %d (nE=%d) %f %f",
// 					        EMS[e].Te, e, nE, NDS[e].T, NDS[e+1].T);
			}
		}

		// In the final step compute temperature and temperature gradient, check both density and mass balance
		for (e = 0; e < nE; e++) {
			EMS[e].gradT = (NDS[e+1].T - NDS[e].T) / EMS[e].L;
			EMS[e].Te = (NDS[e].T + NDS[e+1].T) / 2.0;
			if ((e > Xdata.SoilNode) && (e < nE-1) && ((EMS[e].Te - Constants::melting_tk) > 0.1))
				prn_msg(__FILE__, __LINE__, "wrn", Date(),
				        "Snow temperature above melting point (%f) in element %d (nE=%d) %f %f",
				        EMS[e].Te, e, nE, NDS[e].T, NDS[e+1].T);
			if (EMS[e].theta[SOIL] == 0.) {
				if (!(EMS[e].Rho > 0. && EMS[e].Rho <= Constants::max_rho)) {
					prn_msg(__FILE__, __LINE__, "err", Date(), "Phase Change End: rho[%d]=%f", e, EMS[e].Rho);
					throw IOException("Run-time error in compPhaseChange()", AT);
				}
			}
			// Also make sure the sum of all volumetric contents is near 1 (Can make a 1% error)
			if (!EMS[e].checkVolContent()) {
				prn_msg(__FILE__, __LINE__, "err", Date(),
				        "Phase Change End: Element=%d, nE=%d  ICE %f, Water %f, Air %f Soil %f",
				        e, nE, EMS[e].theta[ICE], EMS[e].theta[WATER], EMS[e].theta[AIR], EMS[e].theta[SOIL]);
				throw IOException("Run-time error in compPhaseChange()", AT);
			}
			cold_content_out += EMS[e].c[TEMPERATURE] * EMS[e].Rho * (EMS[e].Te - Constants::melting_tk) * EMS[e].L;
			sum_Qmf += EMS[e].Qmf * EMS[e].L;
		}
		if (prn_CK && (sum_Qmf > 0.)) { //HACK
			prn_msg(__FILE__, __LINE__, "msg+", Date(), "Checking energy balance  (W/m2):");
			prn_msg(__FILE__, __LINE__, "msg+", Date(), " E1: %f   E0: %f  E1-E0: %f  sum_Qmf: %f  Surface EB : %f",
			           (cold_content_out) / sn_dt, (cold_content_in) / sn_dt,
			               (cold_content_out - cold_content_in) / sn_dt, sum_Qmf,
			                   Sdata.qs + Sdata.ql + Sdata.lw_net + Sdata.qr + Sdata.qw);
		}
	} catch (const exception& ) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Run-time error in compPhaseChange()");
		throw;
	}
}
