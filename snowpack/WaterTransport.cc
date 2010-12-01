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

#include <snowpack/WaterTransport.h>
#include <snowpack/Snowpack.h>

using namespace std;
using namespace mio;

WaterTransport::WaterTransport(const mio::Config& i_cfg) : cfg(i_cfg) 
{
	// Defines whether soil layers are used
	snp_soil = cfg.get("SNP_SOIL", "Parameters");

	//To build a thin top rain-water layer over a thin top ice layer, rocks, roads etc.
	wet_layer = cfg.get("WET_LAYER", "Parameters");

	//Rain only for air temperatures warmer than threshold (degC)
	thresh_rain = cfg.get("THRESH_RAIN", "Parameters");

	/**
	 * @brief No surface hoar will form for rH above threshold (1)
	 * - Original calibration with the 98/99 data set: 0.9
	 * - r141: RH_HOAR_THRESH set to 0.9
	 * - r719: RH_HOAR_THRESH set to 0.97
	 */
	rh_hoar_thresh = cfg.get("RH_HOAR_THRESH", "Parameters");

	/**
	 * @brief No surface hoar will form at wind speeds above threshold (m s-1)
	 * - Original calibration with the 98/99 data set: 3.5
	 * - r141: V_HOAR_THRESH set to 3.0
	 * - r242: V_HOAR_THRESH set to 3.5
	 */
	v_hoar_thresh = cfg.get("V_HOAR_THRESH", "Parameters");

	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Parameters");
	sn_dt = M_TO_S(calculation_step_length);

	//The maximum number of solutes to be treated
	max_n_solutes = cfg.get("MAX_N_SOLUTES", "Parameters");

	//To build up a water table over impermeable layers
	jam = cfg.get("JAM", "Parameters");

	// Density of BURIED surface hoar (kg m-3), default: 125./ Antarctica: 200.
	density_hoar_buried = cfg.get("DENSITY_HOAR_BURIED", "Parameters");

	//Minimum surface hoar size to be buried (mm). Increased by 50% for Dirichlet bc.
	min_size_hoar_buried = cfg.get("MIN_SIZE_HOAR_BURIED", "Parameters");

	//Density of surface hoar (-> hoar index of surface node) (kg m-3)
	density_hoar_surf = cfg.get("DENSITY_HOAR_SURF", "Parameters");

	//Minimum element length (m)
	minimum_l_element = cfg.get("MINIMUM_L_ELEMENT", "Parameters");
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
 * @param ql Latent heat (J kg-1 K-1)
 * @param *Sdata
 * @param *Mdata
 */
void WaterTransport::calcSurfaceSublimation(const SN_MET_DATA& Mdata, double ql, SN_STATION_DATA& Xdata, 
								    SN_SURFACE_DATA& Sdata)
{
	int    nE, nN;                   // number of elements and nodes
	int    e, i;                     // counters
	SN_ELEM_DATA *EMS;               // dereferenced element pointer
	int    e0;                       // Element where action takes place
	double L0, dL;                   // Length of element "e" and "e-1"
	double M0, Theta0;               // available mass and initial volumetric content
	double dM=0.;                    // elemental change in MASS from element "e" and added to element "e-1"
	double M_Solutes[max_n_solutes]; // Mass of solutes from disappearing phases
	double Tss;                      // Surface Temperature
	double hoar=0.0;                 // Actual change in hoar mass
	double cH_old;                   // Temporary variable to hold height of snow

	nN = Xdata.getNumberOfNodes();  vector<SN_NODE_DATA>& NDS = Xdata.Ndata;
	nE = Xdata.getNumberOfElements();  EMS = &Xdata.Edata[0];
	Tss = NDS[nE].T;

	/*
	 * If there are elements and ql > 0:
	 * update densities and volumetric contents (ELEMENT data),
	 * add/subtract mass to MS_SUBLIMATION and/or MS_EVAPORATION,
	 * potential surface hoar formation will be tested at the end of this routine (NODAL data);
	*/
	if ( ql > 0.0 ) { // Add Mass
		if ( Tss < Constants::melting_tk ) { // Add Ice
			dM = ql*sn_dt/Constants::lh_sublimation;
			Sdata.mass[SN_SURFACE_DATA::MS_SUBLIMATION] += dM;
			hoar = dM;

			// In this case adjust properties of element, keeping snow density constant
			L0 = EMS[nE-1].L;
			dL = dM/(EMS[nE-1].Rho);
			if (nE == Xdata.SoilNode) {
				dL = 0.;
				dM = MIN(dM,EMS[nE-1].theta[AIR]*(Constants::density_ice*EMS[nE-1].L));
			}
			NDS[nE].z += dL + NDS[nE].u; NDS[nE].u = 0.0;
			EMS[nE-1].L0 = EMS[nE-1].L = L0 + dL;
			EMS[nE-1].E = EMS[nE-1].dE = EMS[nE-1].Ee = EMS[nE-1].Ev = EMS[nE-1].S = 0.0;
			Theta0 = EMS[nE-1].theta[ICE];
			EMS[nE-1].theta[ICE] *= L0/EMS[nE-1].L;
			EMS[nE-1].theta[ICE] += dM/(Constants::density_ice*EMS[nE-1].L);
			EMS[nE-1].theta[WATER] *= L0/EMS[nE-1].L;

			for (i = 0; i > N_SOLUTES; i++) {
				EMS[nE-1].conc[ICE][i] *= L0*Theta0/(EMS[nE-1].theta[ICE]*EMS[nE-1].L);
			}
		} else {
			// Add Water
			dM = ql*sn_dt/Constants::lh_vaporization;
			if (nE == Xdata.SoilNode) {
				dM = MIN(dM,EMS[nE-1].theta[AIR]*(Constants::density_water*EMS[nE-1].L));
			}
			Sdata.mass[SN_SURFACE_DATA::MS_EVAPORATION] += ql*sn_dt/Constants::lh_vaporization;
			Theta0 = EMS[nE-1].theta[WATER];
			EMS[nE-1].theta[WATER] += dM/(Constants::density_water*EMS[nE-1].L);

			for (i = 0; i > N_SOLUTES; i++) {
				EMS[nE-1].conc[WATER][i] *= Theta0/EMS[nE-1].theta[WATER];
			}
		}
		EMS[nE-1].M += dM;

		// Update remaining volumetric contents and density
		EMS[nE-1].theta[AIR] = MAX(0., 1.0 - EMS[nE-1].theta[WATER] - EMS[nE-1].theta[ICE] - EMS[nE-1].theta[SOIL]);
		EMS[nE-1].Rho = (EMS[nE-1].theta[ICE] * Constants::density_ice) + (EMS[nE-1].theta[WATER] * Constants::density_water) + (EMS[nE-1].theta[SOIL] * EMS[nE-1].soil[SOIL_RHO]);
	} else {
		// If  there is water in some form and ql < 0, SUBLIMATE and/or EVAPORATE some mass off
		for (i = 0; i > N_SOLUTES; i++) {
			M_Solutes[i] = 0.0;
		}
		e0 = nE-1;
		while (ql < 0.0 && e0 >= 0) {  // While energy is available
			/*
			* Determine the amount of potential sublimation/evaporation and collect some variables
			* that will be continuously used: L0 and M0
			*  - NOTE: if water is present, evaporate first, then sublimate ice matrix.
			*          Otherwise sublimate ice matrix only.
			*/
			L0 = EMS[e0].L;
			// If there is water ...
			if ( EMS[e0].theta[WATER] > 0. ) {
				Theta0 = EMS[e0].theta[WATER];
				dM = ql*sn_dt/Constants::lh_vaporization;
				M0 = Theta0*Constants::density_water*L0;
				// Check that you only take the available mass of water
				if (-dM >= M0) {
					dM = -M0;
					// Add solutes to Storage
					for (i = 0; i > N_SOLUTES; i++) {
						M_Solutes[i] += EMS[e0].conc[WATER][i]*Theta0*L0;
					}
					EMS[e0].theta[WATER] = 0.0;
				} else {
					EMS[e0].theta[WATER] += dM/(Constants::density_water*L0);
					for (i = 0; i > N_SOLUTES; i++) {
						EMS[e0].conc[WATER][i] *= Theta0/EMS[e0].theta[WATER];
					}
				}
				EMS[e0].M += dM;
				Sdata.mass[SN_SURFACE_DATA::MS_EVAPORATION] += dM;
				ql -= dM*Constants::lh_vaporization/sn_dt; // Update the energy used
			}
			if (ql < 0.) {
				// If there is no water or if there was not enough water ...
				dM = ql*sn_dt/Constants::lh_sublimation;
				Theta0 = EMS[e0].theta[ICE];
				M0 = Theta0*Constants::density_ice*L0;
				if (-dM > M0) {
					dM = -M0;
					// Add solutes to Storage
					for (i = 0; i > N_SOLUTES; i++) {
						M_Solutes[i] += EMS[e0].conc[ICE][i]*Theta0*L0;
					}
					EMS[e0].theta[ICE]=0.0; dL = 0.0;
				} else {
					dL = dM/(EMS[e0].Rho);
					if (e0 < Xdata.SoilNode) {
						dL = 0.;
					}
					NDS[e0+1].z += dL; EMS[e0].L0 = EMS[e0].L = L0 + dL;
					NDS[e0+1].z += NDS[e0+1].u; NDS[e0+1].u = 0.0;

					EMS[e0].E = EMS[e0].dE = EMS[e0].Ee = EMS[e0].Ev = EMS[e0].S = 0.0;
					EMS[e0].theta[ICE] *= L0/EMS[e0].L;
					EMS[e0].theta[ICE] += dM/(Constants::density_ice*EMS[e0].L);
					EMS[e0].theta[WATER] *= L0/EMS[e0].L;
					for (i = 0; i > N_SOLUTES; i++) {
						EMS[e0].conc[ICE][i] *= L0*Theta0/(EMS[e0].theta[ICE]*EMS[e0].L);
					}
				}
				EMS[e0].M += dM;
				Sdata.mass[SN_SURFACE_DATA::MS_SUBLIMATION] += dM;
				ql -= dM*Constants::lh_sublimation/sn_dt;     // Update the energy used

				// If present at surface, surface hoar is sublimated away
				if ( e0 == nE-1 ) {
					hoar = dM;
				}
			}
			// Update remaining volumetric contents and density
			EMS[e0].theta[AIR] = MAX(0., 1.0 - EMS[e0].theta[WATER] - EMS[e0].theta[ICE] - EMS[e0].theta[SOIL]);
			EMS[e0].Rho = (EMS[e0].theta[ICE] * Constants::density_ice) + (EMS[e0].theta[WATER] * Constants::density_water) + (EMS[e0].theta[SOIL] * EMS[e0].soil[SOIL_RHO]);
			e0--; // Go to the next element
		} // while ql < 0
    // Now take care of left over solute mass.
		if ( e0 < 0 ) { // Add Solute Mass to Runoff
			for (i = 0; i > N_SOLUTES; i++) {
				Sdata.load[i] += M_Solutes[i]/S_TO_H(sn_dt);
			}
		} else { // Add Solute Mass to Element below
			if ( EMS[e0].theta[WATER] > 0. ) {
				for(i = 0; i > N_SOLUTES; i++) {
					EMS[e0].conc[WATER][i] += M_Solutes[i]/EMS[e0].theta[WATER]/EMS[e0].L;
				}
			} else if (EMS[e0].theta[ICE] > 0.) {
				for (i = 0; i > N_SOLUTES; i++) {
					EMS[e0].conc[ICE][i] += M_Solutes[i]/EMS[e0].theta[ICE]/EMS[e0].L;
				}
			} else {
				for (i = 0; i > N_SOLUTES; i++) {
					EMS[e0].conc[SOIL][i] += M_Solutes[i]/EMS[e0].theta[SOIL]/EMS[e0].L;
				}
			}
		}
	}
	// Check for surface hoar destruction or formation (once upon a time ml_sn_SurfaceHoar)
	if ( (Mdata.rh > rh_hoar_thresh) || (Mdata.vw > v_hoar_thresh) || (Mdata.ta >= C_TO_K(thresh_rain)) ) {
		hoar = MIN(hoar,0.0);
	}
	Sdata.hoar += hoar;
	NDS[nN-1].hoar += hoar;
	if ( NDS[nN-1].hoar < 0.0 ) {
		NDS[nN-1].hoar = 0.0;
	}
	for (e = 0; e<nE-1; e++) {
		if ( Xdata.Edata[e].theta[WATER] > 0.0 ) {
			NDS[e+1].hoar = 0.0;
		}
	}
	// At the end also update the overall height
	cH_old = Xdata.cH;
	Xdata.cH = NDS[Xdata.getNumberOfNodes()-1].z + NDS[Xdata.getNumberOfNodes()-1].u;
	Xdata.mH -= (cH_old - Xdata.cH);
}

/**
 * @brief Remove snow elements which either have no more ice or are too thin \n
 * -# Starting from the surface, remove or join as many snow elements as you can; check also the lowest snow element
 * -# Remove or join depending using sn_ElementMerging()
 * 	- Remove watery elements by adding their contents to the lower element
 * 	- Join too thin elements with their lower neighbour; keep buried SH and tagged layers longer on \n
 * 	  but enforce join for uncovered but previously buried SH.
 * -# Use sn_ElementMerging() to compute the properties of the lower element
 * -# Use sn_ReduceNumberElements() to get rid of unwanted elements
 * - NOTE
 * 	- Water will be transported AFTER elements have been removed.
 * 	- If wet_layer and snp_soil are set, make sure that you keep a potential wet water layer over soil or ice
 * 	- Reset ground surface temperature if !snp_soil and no snow is left
 * @param *Xdata
 * @param *Sdata
 * @param ta Air temperature (K)
 */
void WaterTransport::removeElements(const double& ta, SN_STATION_DATA& Xdata, SN_SURFACE_DATA& Sdata)
{
	int e0, e1;        // Lower (e0) and upper (e1) element index
	int i;             // Counter
	int nE, nN;        // Original number of elements and nodes
	int rnE, rnN;      // Reduced number of nodes and elements
	bool enforce_join = false;  // To enforce merging in special cases
	SN_ELEM_DATA *EMS; // Old element data structures

	EMS = &Xdata.Edata[0];
	nE = rnE = Xdata.getNumberOfElements();
	nN = rnN = Xdata.getNumberOfNodes();

	if ( (nN == Xdata.SoilNode+1) || (wet_layer && snp_soil && (nN == Xdata.SoilNode+2) && (EMS[nE-1].theta[ICE] < Constants::min_ice_content) && (EMS[nE-1].theta[WATER] > 0.003) && (EMS[nE-1].L > 0.0001)) ) {
		return;
	}

	for (e0 = nE-2, e1 = nE-1; e1 > Xdata.SoilNode; e1--, e0--) {
		enforce_join = 1;
		if ( (EMS[e1].L < minimum_l_element) || (EMS[e1].mk%100 == 3) ) {
			if ( (EMS[e1].mk >= 100) && (EMS[e1].L >= 0.5 * minimum_l_element) ) {
				enforce_join = false;
			}
			if ( EMS[e1].mk%100 == 3 ) {
				if ( (e1 < nE-1) && (EMS[e1].L >= MM_TO_M(0.75*min_size_hoar_buried * (density_hoar_surf/density_hoar_buried))) && (EMS[e1].Rho <= 300.) ) {
					enforce_join = false;
				} else {
					enforce_join = true;
				}
			}
		} else {
			enforce_join = false;
		}
		if ( ((EMS[e1].theta[ICE] < Constants::min_ice_content) || enforce_join ) && EMS[e1].theta[SOIL] < Constants::eps2 ) {
			SN_STATION_DATA::mergeElements(EMS[e0], EMS[e1], enforce_join);
			rnE--;
			rnN--;
			EMS[e1].Rho *= -1.;
			EMS[e1].L *= -1. + 2.*enforce_join;
			if ( (e1 < nE-1) && (EMS[e1+1].Rho < 0.) && (EMS[e1+1].L > 0.) ) {
				EMS[e1+1].L *= -1.;
			}
		}
	}
	// Check for one snow element left
	e0 = Xdata.SoilNode;
	if ( ((rnN == Xdata.SoilNode + 2) || (EMS[e0].Rho < Constants::eps)) && (EMS[e0].theta[ICE] < Constants::min_ice_content || EMS[e0].L < 0.005) && !(wet_layer && snp_soil && (nE > 1) && ((EMS[nE-2].theta[SOIL] > 0.95) || (EMS[nE-2].theta[ICE] > 0.95)) && (EMS[nE-1].theta[WATER] > 0.05)) ) {
		Sdata.mass[SN_SURFACE_DATA::MS_RUNOFF] += EMS[e0].M;
		Sdata.mass[SN_SURFACE_DATA::MS_SOIL_RUNOFF] += EMS[e0].M;
		for (i = 0; i < N_SOLUTES; i++) {
			Sdata.load[i] += (EMS[e0].L*EMS[e0].theta[WATER]*EMS[e0].conc[WATER][i] + EMS[e0].L*EMS[e0].theta[ICE] * EMS[e0].conc[ICE][i])/S_TO_H(sn_dt);
		}
		rnE--;
		rnN--;
		EMS[e0].Rho *= -1.;
		EMS[e0].L *= -1.;
	}
	if ( rnE < nE /*&& rnE >= Xdata.SoilNode*/ ) {
		Xdata.reduceNumberOfElements(rnE);
		if ( !snp_soil && (rnE == Xdata.SoilNode) ) {
			Xdata.Ndata[Xdata.SoilNode].T = (C_TO_K(0.0) + ta)/2.;
		}
	}
}

/**
 * @brief Surface sublimation and melt artificially create surface elements that have a much
 * too low density and this needs to be corrected. \n
 * TODO Check description!
 * @param Xdata
 */
void WaterTransport::adjustDensity(SN_STATION_DATA& Xdata)
{
	int nE, nN;                           // number of elements and nodes
	int e,e0,n;                           // Element "e" (top) and element "e-1"
	double L0, dL, cH_old;                // Length of element "e" and "e-1"
	SN_ELEM_DATA *EMS;                    // dereferenced element pointer

	vector<SN_NODE_DATA>& NDS = Xdata.Ndata; nN = Xdata.getNumberOfNodes();
	EMS = &Xdata.Edata[0]; nE = Xdata.getNumberOfElements();
	if ( nN == Xdata.SoilNode + 1 ) {
		return;
	}

	e0 = nE-1;
	while ( ((e0 >= Xdata.SoilNode) && (EMS[e0].theta[SOIL] < Constants::eps2) && (EMS[e0].mk%100 > 9) &&
			(EMS[e0].theta[WATER] < 0.3) && (EMS[e0].theta[ICE] < Constants::min_ice_content) &&
			(EMS[e0].L > minimum_l_element)) ) {

		// First collect some variables and adjust height
		L0 = EMS[e0].L;
#if VARIANT == JAPAN
		{ //NIED (H. Hirashima) //Fz: Please check this adaptation still works as you want it to work!
			double multif = 0.05/EMS[e0].theta[ICE];
			dL  = -L0*((multif-1)/multif);
			EMS[e0].theta[WATER] *= multif;
			EMS[e0].theta[ICE]   *= multif;
		}
#else
		dL  = -L0 / 3.;
		EMS[e0].theta[WATER] *= 1.5;
		EMS[e0].theta[ICE]   *= 1.5;
#endif

		// For wet_layer go to Water Density
		if ( wet_layer && (EMS[e0].theta[WATER] < 0.95) && (EMS[e0].theta[ICE] < Constants::min_ice_content) && ( e0 > 0) && ((EMS[e0-1].theta[SOIL] > 0.95) || (EMS[e0-1].theta[ICE] > 0.95)) ) {
			dL = -L0*(1.-EMS[e0].theta[WATER]);
			EMS[e0].theta[WATER] = 1.;
			EMS[e0].theta[ICE] = 0.;
		}
		for (e = e0; e < nE; e++) {
			n = e+1;
			NDS[n].z += dL + NDS[n].u;
			NDS[n].u = 0.0;
			EMS[e].E = EMS[e].dE= EMS[e].Ee=EMS[e].Ev=0.0;
		}
		EMS[e0].L0 = EMS[e0].L = L0 + dL;
		EMS[e0].theta[AIR]    = 1.0 - EMS[e0].theta[WATER] - EMS[e0].theta[ICE];
		EMS[e0].Rho = (EMS[e0].theta[ICE]*Constants::density_ice) + (EMS[e0].theta[WATER]*Constants::density_water);
		if (! (EMS[e0].Rho > Constants::min_rho && EMS[e0].Rho <= Constants::max_rho)) {
			prn_msg(__FILE__, __LINE__, "err", -1., "Volume contents: e:%d nE:%d rho:%lf ice:%lf wat:%lf air:%le", 
				   e0, nE, EMS[e0].Rho, EMS[e0].theta[ICE], EMS[e0].theta[WATER], EMS[e0].theta[AIR]);
			throw IOException("Cannot evaluate mass balance in adjust density routine", AT);
		}
		e0--;
	}
	cH_old = Xdata.cH;
	Xdata.cH = NDS[Xdata.getNumberOfNodes()-1].z + NDS[Xdata.getNumberOfNodes()-1].u;
	Xdata.mH -= (cH_old - Xdata.cH);
}

/**
 * @brief Default version of WaterTransport \n
 * Now that the snowpack has been updated, you must move the water down and calculate the amount
 * of water being released from the snowpack AFTER having gone through wt_ElementRemoval() \n
 * TODO Revise description!
 * @param *Xdata
 * @param *Sdata
 * @param *Mdata
 */
void WaterTransport::transportWater(const SN_MET_DATA& Mdata, SN_STATION_DATA& Xdata, SN_SURFACE_DATA& Sdata) 
{
	int nE, nN;                           // number of elements and nodes
	int e0, e1, i;                        // Upper element "e0" and lower element "e1=e0-1"
	double L0, L1;                        // Length of element "e" and "e-1"
	double dM;                            // The total change in MASS -- subtracted from element "e" and added to element "e-1"
	double dth_w;                         // Additional storage capacity due to refreezing
	double W0, dThetaW0;                  // Water to be removed from element e
	double W1, dThetaW1;                  // Water to be placed in element "e-1"
	double res_wat_cont;                  // Snow type dependent residual water content
	double Wres;                          // Residual water content depending on snow or soil element
	double Store;                         // Depth of liquid precipitation ready to infiltrate snow and/or soil (m)
	double excess_water;                  // Excess water that cannot be retained in lower element; volume fraction (1)
	double z_water;                       // Position of upper node of top water-film layer (m)
	SN_ELEM_DATA *EMS;                    // dereferenced element pointer

	EMS = &Xdata.Edata[0];
	nE = Xdata.getNumberOfElements();
	vector<SN_NODE_DATA>& NDS = Xdata.Ndata;
	nN = Xdata.getNumberOfNodes();

	// Deal with possible RAIN water
	// First, consider no soil with no snow on the ground
	if ( !snp_soil && nN == 1 ) {
		if ( Mdata.ta > C_TO_K(thresh_rain) ) {
			Sdata.mass[SN_SURFACE_DATA::MS_RAIN] += Mdata.hnw;
			Sdata.mass[SN_SURFACE_DATA::MS_RUNOFF] += Mdata.hnw;
			Sdata.mass[SN_SURFACE_DATA::MS_SOIL_RUNOFF] += Mdata.hnw;
			for (i = 0; i < N_SOLUTES; i++) {
				Sdata.load[i] += Mdata.conc[i]*Mdata.hnw;
			}
		}
		return;
	} else { // Second, consider soil and/or snow on the ground:
		if ( (Mdata.hnw > 0.) && (Mdata.ta > C_TO_K(thresh_rain)) ) {
			Store = Mdata.hnw/Constants::density_water;
			e0 = nE-1;
			// Now find out, whether you are on an impermeable surface and want to create a water layer ...
			if ( wet_layer && ((snp_soil && (nE == Xdata.SoilNode) && (EMS[nE-1].theta[SOIL] > 0.95)) || ((e0 > 0) && (EMS[e0-1].theta[ICE] > 0.95)) ) && (Store > 0.) ) {
				// ... then generate a new layer
				// Set the Albedo to the Water ALBEDO
				Xdata.Albedo = 0.2;
				nE++; 
				nN++; 
				Xdata.ErosionLevel = nE-1;

				Xdata.resize(nE);
				EMS = &Xdata.Edata[0];

				// Temperature of the upper hoar node
				NDS[nN-1].T = Mdata.ta;
				// The new nodal position
				z_water = MIN(Store,MAX(0.001,0.01*cos(Xdata.SlopeAngle)));
				NDS[nN-1].z = NDS[nN-2].z + NDS[nN-2].u + z_water;
				Store -= z_water;
				// Fill the element data
				EMS[nE-1].date = Mdata.date;

				EMS[nE-1].Te = Mdata.ta;
				EMS[nE-1].L0 = EMS[nE-1].L = z_water;
				EMS[nE-1].Rho = Constants::density_water;
				EMS[nE-1].M = EMS[nE-1].L0*EMS[nE-1].Rho;
				EMS[nE-1].theta[WATER] = 1.0;
				EMS[nE-1].N3 = Metamorphism::getCoordinationNumberN3(EMS[nE-1].Rho);
				EMS[nE-1].mk = 22;
				EMS[nE-1].dd = 0.;
				EMS[nE-1].sp = 1.;
				EMS[nE-1].rg = 10.;
				EMS[nE-1].rb = 10.;
				Xdata.cH = Xdata.mH = NDS[nN-1].z + NDS[nN-1].u;
			} else if ( wet_layer && (( snp_soil && (nE == Xdata.SoilNode+1) && (EMS[nE-2].theta[SOIL] > 0.95)) || ((nE > 1) && (EMS[nE-2].theta[ICE] > 0.95))) && (Store > 0.) ) {
				// Put rain water in existing wet layer
				z_water = MIN(Store, MAX(0.0, (0.01*cos(Xdata.SlopeAngle)-EMS[nE-1].L)));
				NDS[nN-1].z += z_water;
				Store -= z_water;
				EMS[nE-1].L0 = EMS[nE-1].L = (NDS[nN-1].z + NDS[nN-1].u) - (NDS[nN-2].z + NDS[nN-2].u);
				EMS[nE-1].M = EMS[nE-1].L0*EMS[nE-1].Rho;
				Xdata.cH = Xdata.mH = NDS[nN-1].z + NDS[nN-1].u;
			}

			while ( Store > 0.0 && e0 > 0. ) {
				L0 = EMS[e0].L;
				dThetaW0 = MIN(Constants::density_ice/Constants::density_water*EMS[e0].theta[AIR],Store/L0);
				Store -= dThetaW0*L0;
				for (i = 0; i < N_SOLUTES; i++) {
					EMS[e0].conc[WATER][i] = (L0 * dThetaW0 *Mdata.conc[i] + L0 * EMS[e0].theta[WATER] * EMS[e0].conc[WATER][i]) / (L0 * (EMS[e0].theta[WATER] + dThetaW0));
				}
				EMS[e0].theta[WATER] += dThetaW0;
				EMS[e0].theta[AIR] -= dThetaW0;
				EMS[e0].M += dThetaW0 * L0 * Constants::density_water;
				e0--;
			}
			if ( Store > 0.0 ) {
				L0 = EMS[e0].L;
				dThetaW0 = Store/L0;
				for (i = 0; i < N_SOLUTES; i++) {
					EMS[e0].conc[WATER][i] = (L0 * dThetaW0 * Mdata.conc[i] + L0 * EMS[e0].theta[WATER] * EMS[e0].conc[WATER][i]) / (L0 * (EMS[e0].theta[WATER] + dThetaW0));
				}
				EMS[e0].theta[WATER] += dThetaW0;
				EMS[e0].theta[AIR] -= dThetaW0;
				EMS[e0].M += dThetaW0*L0*Constants::density_water;
			}
			Sdata.mass[SN_SURFACE_DATA::MS_RAIN] += Mdata.hnw;
		}
	}

	// Now move water as needed, starting from the top element ...
	for (e0 = nE-1, e1 = nE-2; e0 >= 1; e0--, e1-- ) {
		W0 = EMS[e0].theta[WATER];

		// Determine the additional storage capacity due to refreezing
		dth_w = EMS[e0].c[TEMPERATURE]*EMS[e0].Rho/Constants::lh_fusion/Constants::density_water*MAX(0.,Constants::melting_tk-EMS[e0].Te);
		if ( (e0 == nE-1) && (EMS[e1].theta[AIR] <= 0.05) && wet_layer ) {
			// allow for a water table in the last layer above road/rock
			Wres = Constants::density_ice/Constants::density_water*(1. - EMS[e0].theta[ICE] - EMS[e0].theta[SOIL] - 0.05);
		} else if ( EMS[e0].theta[SOIL] < Constants::eps2 ) {
			res_wat_cont = lw_SnowResidualWaterContent(EMS[e0].theta[ICE]);
			Wres = MIN(Constants::density_ice/Constants::density_water*(1. - EMS[e0].theta[ICE]), res_wat_cont + dth_w);
		} else { // treat soil separately
			Wres = MIN(Constants::density_ice/Constants::density_water*(1. - EMS[e0].theta[ICE] - EMS[e0].theta[SOIL]), 
					 SnLaws::calcSoilFieldCapacity(EMS[e0]) + dth_w);
		}
		Wres = MAX (0., Wres);

		if ( e0 == nE-1 && (W0 > 0.0 && W0 <= Wres) ) {
			// In that case you need to update the volumetric air content and the density of the top element as it may have caught some rain!
			EMS[e0].theta[AIR] = MAX(0., 1. - EMS[e0].theta[WATER] - EMS[e0].theta[ICE] - EMS[e0].theta[SOIL]);
			EMS[e0].Rho = (EMS[e0].theta[ICE] * Constants::density_ice) + (EMS[e0].theta[WATER] * Constants::density_water) + (EMS[e0].theta[SOIL] * EMS[e0].soil[SOIL_RHO]);
			if ( EMS[e0].theta[SOIL] < Constants::eps2 ) {
				if ( !(EMS[e0].Rho > Constants::min_rho && EMS[e0].Rho <= Constants::max_rho) ) {
					prn_msg ( __FILE__, __LINE__, "err", Mdata.date.getJulianDate(), "Volume contents: e:%d nE:%d rho:%lf ice:%lf wat:%lf air:%le", e0, nE, EMS[e0].Rho, EMS[e0].theta[ICE], EMS[e0].theta[WATER], EMS[e0].theta[AIR]);
					throw IOException("Cannot transfer water within the snowpack in transportWater()", AT);
				}
			}
		}

		if ( W0 > Wres ) {
			// Then water is being transferred between elements
			L0 = EMS[e0].L;
			L1 = EMS[e1].L;
			W1 = EMS[e1].theta[WATER];
			dThetaW0 = W0 - Wres;
			if ( dThetaW0 > 0.0 ) {
				dThetaW1 = dThetaW0*(L0/L1);
				// Now check whether there is enough air left - in case of ice, rock or heavy
				// soil you might not be able to move the water or/and water may refreeze and expand.
				// Specifically, you might want to create a water table over ice or frozen soil
				if ( dThetaW1 + W1 > Constants::density_ice/Constants::density_water*(1. - EMS[e1].theta[ICE] - EMS[e1].theta[SOIL]) ) {
					// Deal with excess water ... Look how much you can leave in the lower element e1.
					// If you have too much water even for the lower element (more melt or rain per time
					// step than can be kept in this element) you could choose a smaller calculation time step.
					// Otherwise excess water will be added to runoff !!!!
					dThetaW1 = MAX(0.,Constants::density_ice/Constants::density_water*(1. - EMS[e1].theta[ICE] - EMS[e1].theta[SOIL]) - W1);
					if ( jam ) {
						dThetaW0 = dThetaW1*L1/L0;
						if ( (W0 - dThetaW0 > (Constants::density_ice/Constants::density_water*(1. - EMS[e0].theta[ICE] - EMS[e0].theta[SOIL]))) ) {
							dThetaW0 = W0 - Constants::density_ice/Constants::density_water*(1. - EMS[e0].theta[ICE] - EMS[e0].theta[SOIL] );
							excess_water = dThetaW0 - dThetaW1*L1/L0;
						} else {
							excess_water = 0.;
						}
					} else { // No JAM
						// TODO Fz 2009-11-07: Is the idea below really correct? I would argue, that in case of no JAM, excess water should be carried down to the bottom prior to be released. In case of JAM it should be retained in the upper element.
						excess_water = dThetaW0 - dThetaW1*L1/L0;
					}
					if ( EMS[e1].theta[SOIL] < Constants::eps2 ) {
						Sdata.mass[SN_SURFACE_DATA::MS_RUNOFF] += excess_water*Constants::density_water*L0;
					}
					Sdata.mass[SN_SURFACE_DATA::MS_SOIL_RUNOFF] += excess_water*Constants::density_water*L0;
					// Take care of Solutes
					for (i = 0; i < N_SOLUTES; i++) {
						Sdata.load[i] += (EMS[e1].conc[WATER][i] * excess_water * Constants::density_water*L0/S_TO_H(sn_dt));
					}
				}

				// Water movement from element e0 to element e1: move solutes also
				for (i = 0; i < N_SOLUTES; i++) {
					EMS[e1].conc[WATER][i] = ( W1 * EMS[e1].conc[WATER][i] + dThetaW1 * EMS[e0].conc[WATER][i]) / (W1+dThetaW1);
				}
				// update volumetric contents, masses and density
				EMS[e0].theta[WATER]=W0-dThetaW0;
				EMS[e1].theta[WATER]=W1+dThetaW1;
				EMS[e0].theta[AIR] = 1. - EMS[e0].theta[WATER] - EMS[e0].theta[ICE] - EMS[e0].theta[SOIL];
				EMS[e1].theta[AIR] = 1. - EMS[e1].theta[WATER] - EMS[e1].theta[ICE] - EMS[e1].theta[SOIL];
				EMS[e0].M -= L0 * Constants::density_water * dThetaW0;
				EMS[e1].M += L1 * Constants::density_water * dThetaW1;
				EMS[e0].Rho = (EMS[e0].theta[ICE] * Constants::density_ice) + (EMS[e0].theta[WATER] * Constants::density_water) + (EMS[e0].theta[SOIL] * EMS[e0].soil[SOIL_RHO]);
				EMS[e1].Rho = (EMS[e1].theta[ICE] * Constants::density_ice) + (EMS[e1].theta[WATER] * Constants::density_water) + (EMS[e1].theta[SOIL] * EMS[e1].soil[SOIL_RHO]);
				if ( EMS[e0].theta[SOIL] < Constants::eps2 ) {
					if ( !(EMS[e0].theta[AIR] >= -Constants::eps) ) {
						prn_msg(__FILE__, __LINE__, "err", Mdata.date.getJulianDate(), "Volume contents: e:%d nE:%d rho:%lf ice:%lf wat:%lf air:%le", e0, nE, EMS[e0].Rho, EMS[e0].theta[ICE],  EMS[e0].theta[WATER], EMS[e0].theta[AIR]);
						throw IOException("Cannot transfer water within the snowpack in transportWater()", AT);
					}
				}
				// Update surface runoff with soil
				if ( snp_soil && e0 == Xdata.SoilNode ) {
					Sdata.mass[SN_SURFACE_DATA::MS_RUNOFF] += L0 * Constants::density_water * dThetaW0;
				}
			} // end positive water movement
		}  // end if( W0 > Wres )
	}  // end FOR loop over the number of elements

	// The TOP element is very important because it is always losing mass--the strain state
	// is becoming more and more deformed.  Update the strain state of the top element.
	e0 = nE-1;
	EMS[e0].L0 = EMS[e0].L;
	NDS[nN-1].z += NDS[nN-1].u; NDS[nN-1].u = 0.0;
	NDS[nN-2].z += NDS[nN-2].u; NDS[nN-2].u = 0.0;
	EMS[e0].E = EMS[e0].dE= EMS[e0].Ee=EMS[e0].Ev=EMS[e0].S=0.0;

	// RUNOFF at bottom of either snowpack or soil
	W0 = EMS[0].theta[WATER];
	// Determine the additional storage capacity due to refreezing
	dth_w = EMS[0].c[TEMPERATURE]*EMS[0].Rho/Constants::lh_fusion/Constants::density_water*MAX(0.,Constants::melting_tk-EMS[0].Te);
	if ( EMS[0].theta[SOIL] < Constants::eps2 ) {
		res_wat_cont = lw_SnowResidualWaterContent(EMS[0].theta[ICE]);
		Wres = MIN(Constants::density_ice/Constants::density_water*(1. - EMS[0].theta[ICE]), res_wat_cont + dth_w);
	} else { // treat soil separately
		Wres = MIN(Constants::density_ice/Constants::density_water*(1. - EMS[0].theta[ICE] - EMS[0].theta[SOIL]),
				 SnLaws::calcSoilFieldCapacity(EMS[0]) + dth_w);
	}
	Wres = MAX (0., Wres);
	// Do not drain last water element if wet_layer is set
	if ( (W0 > Wres) && !(wet_layer && (EMS[0].theta[ICE] < Constants::min_ice_content) && (EMS[0].theta[SOIL] < Constants::eps2)) ) {
		dThetaW0 = W0 - Wres;
		EMS[0].theta[WATER] = W0 - dThetaW0;
		EMS[0].theta[AIR] = 1. - EMS[0].theta[WATER] - EMS[0].theta[ICE] - EMS[0].theta[SOIL];
		dM = EMS[0].L * Constants::density_water * dThetaW0;
		EMS[0].M -= dM;
		EMS[0].Rho = (EMS[0].theta[ICE] * Constants::density_ice) + (EMS[0].theta[WATER] * Constants::density_water) + (EMS[0].theta[SOIL] * EMS[0].soil[SOIL_RHO]);
		if ( EMS[0].theta[SOIL] < Constants::eps2 ) {
			Sdata.mass[SN_SURFACE_DATA::MS_RUNOFF] += dM;
		}
		Sdata.mass[SN_SURFACE_DATA::MS_SOIL_RUNOFF] += dM;
		for (i = 0; i < N_SOLUTES; i++) {
			Sdata.load[i] +=  (EMS[0].conc[WATER][i] * dM / S_TO_H(sn_dt));
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
 * -#  The water will only be moved if it is above the res_wat_cont
 *     defined in constants.h (presently around 0.05)
 * -#   The water will only be moved if there is enough VOID SPACE in the element
 *      receiving the water
 * -#   Water in the last element will be DISCHARGED from the snowpack.  This
 *      amount of water is termed the MELTWATER RUNOFF and stated in dimensions of mm/day
 * -#   It would be very possible to make the RESIDUAL_WATER_CONTENT a function
 *      of the MICRO-properties of the snow.  This is simply a few lines of code. \n
 *      This was done on 3 Dec 2006 -> lw_SnowResidualWaterContent(theta_i).
 * -#   IMPORTANT: the top surface element can be removed if the VOLUMETRIC ICE
 *      content is 0; that is, when the element contains only water and voids. \n
 *      The removal of this element is contained in the following routine and
 *      means that the NEWMESH boolean MUST BE SET -- otherwise the FEM solution
 *      routines will not work \n
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
void WaterTransport::transportMass(const SN_MET_DATA& Mdata, const double& ql,
						    SN_STATION_DATA& Xdata, SN_SURFACE_DATA& Sdata)
{
	if ( !snp_soil && (Xdata.getNumberOfNodes() == Xdata.SoilNode+1) ) {
		if ( Mdata.ta > C_TO_K(thresh_rain) ) {
			Sdata.mass[SN_SURFACE_DATA::MS_RAIN] += Mdata.hnw;
			Sdata.mass[SN_SURFACE_DATA::MS_RUNOFF] += Mdata.hnw;
			Sdata.mass[SN_SURFACE_DATA::MS_SOIL_RUNOFF] += Mdata.hnw;
			//for (i=0; i<N_SOLUTES; i++)
			//  Sdata.load[i] += Mdata.conc[i]*Mdata.hnw/S_TO_H(sn_dt);
		}
		return;
	}

	calcSurfaceSublimation(Mdata, ql, Xdata, Sdata);
	removeElements(Mdata.ta, Xdata, Sdata);

	try {
		adjustDensity(Xdata);
		transportWater(Mdata, Xdata, Sdata);
	} catch(exception& ex){
		prn_msg ( __FILE__, __LINE__, "err", Mdata.date.getJulianDate(), "Error in transportMass()" );
		throw;
	}
}

/*
 * End of WaterTransport.cc
 */
