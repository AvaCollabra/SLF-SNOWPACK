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
/**
 * @file Snowpack.cc
 * @version 11.06
 * @bug     -
 * @brief This module contains the driving routines for the 1d snowpack model
 */

#include <snowpack/Snowpack.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/

//Uses an empirically determined size of deposited hydrometeors as new snow grain size (mm)
const bool Snowpack::hydrometeor = false;

//Warning is issued if depth of snowfall is larger than this amount (m)
const double Snowpack::snowfall_warning = 0.5;

const unsigned int Snowpack::new_snow_marker = 0;
const double Snowpack::new_snow_albedo = 0.9;

/// @brief Min volumetric ice content allowed
const double Snowpack::min_ice_content = SnLaws::min_hn_density / Constants::density_ice;

/************************************************************
 * non-static section                                       *
 ************************************************************/

Snowpack::Snowpack(const mio::Config& i_cfg) : cfg(i_cfg)
{
	cfg.getValue("VARIANT", "SnowpackAdvanced", variant);

	cfg.getValue("FIXED_ALBEDO", "SnowpackAdvanced", fixed_albedo);

	cfg.getValue("HN_DENSITY", "SnowpackAdvanced", hn_density);
	cfg.getValue("HN_DENSITY_MODEL", "SnowpackAdvanced", hn_density_model);

	/** Defines the management of the bottom boundary conditions with soil layers
	 * - 0 ==> Dirichlet, i.e fixed Temperature
	 * - 1 ==> Neumann, fixed geothermal heat flux GEO_HEAT */
	cfg.getValue("SOIL_FLUX", "Snowpack", soil_flux);

	cfg.getValue("GEO_HEAT", "Snowpack", geo_heat); //Constant geothermal heat flux at (great) depth (W m-2)

	/* Defines the management of the surface boundary conditions
	 * - 0: Neumann boundary conditions throughout
	 * - 1: Dirichlet if Tss < THRESH_CHANGE_BC, Neumann else */
	cfg.getValue("CHANGE_BC", "Snowpack", change_bc);
	cfg.getValue("THRESH_CHANGE_BC", "Snowpack", thresh_change_bc);

	//Should be 0 for data-sets which do not provide measured surface temperatures
	cfg.getValue("MEAS_TSS", "Snowpack", meas_tss);

	/**
	 * @brief Defines how the height of snow is going to be handled
	 * - 0: Depth of snowfall is determined from the water equivalent of snowfall (HNW)
	 * - 1: The measured height of snow is used to determine whether new snow has been deposited.
	 *      This setting MUST be chosen in operational mode. \n
	 *      This procedure has the disadvantage that if the snowpack settles too strongly
	 *      extra mass is added to the snowpack. \n
	 * New snow density is needed in both cases, either parameterized, measured, or fixed.
	 */
	cfg.getValue("ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack", enforce_measured_snow_heights);

	/**
	 * @brief Defines whether the canopy model is used \n
	 * NOTE: OUT_CANOPY must also be set to dump canopy parameters to file; see Constants_local.h
	 */
	cfg.getValue("CANOPY", "Snowpack", useCanopyModel);

	/**
	 * @brief Define the heights of the meteo measurements above ground (m) \n
	 * Required for surface energy exchange computation and for drifting and blowing snow.
	 */
	cfg.getValue("HEIGHT_OF_METEO_VALUES", "Snowpack", height_of_meteo_values);

	/**
	 * @brief Defines whether the measured shortwave radiation is incoming
	 * - 0 downward SW radiation is used
	 * - 1 reflected SW radiation is used
	 * - 2 both downward and reflected SW radiation is used \n
	 * @note { If SW_MODE == 2, the input must hold both fluxes! }
	 */
	cfg.getValue("SW_MODE", "Snowpack", sw_mode);
	sw_mode %= 10;

	/**
	 * @brief Height of new snow element (m) [NOT read from CONSTANTS_User.INI] \n
	 * Controls the addition of new snow layers. Set in qr_ReadParameters() \n
	 * The value depends on ENFORCE_MEASURED_SNOW_HEIGHTS:
	 * - 0: 2.0*MINIMUM_L_ELEMENT (value depends on VARIANT)
	 * - 1: 0.02
	 */
	cfg.getValue("HEIGHT_NEW_ELEM", "SnowpackAdvanced", height_new_elem);
	cfg.getValue("MINIMUM_L_ELEMENT", "SnowpackAdvanced", minimum_l_element);

	// Defines whether soil layers are used
	cfg.getValue("SNP_SOIL", "Snowpack", useSoilLayers);

	cfg.getValue("RESEARCH", "SnowpackAdvanced", research_mode);

	cfg.getValue("VISCOSITY_MODEL", "SnowpackAdvanced", viscosity_model);

	//Rain only for air temperatures warmer than threshold (degC)
	cfg.getValue("THRESH_RAIN", "SnowpackAdvanced", thresh_rain);

	/**
	 * @brief Precipitation only for rH above threshold (1)
	 * - default: 0.50
	 * 	- 2007-12-01: set THRESH_RH to 0.70 to be consistent with data range of ZWART new snow density model
	 * 	- 2008-01-21: set back THRESH_RH to 0.50 (IMIS sensor problem in operational mode)
	 * - Antarctica: 0.70
	 */
	cfg.getValue("THRESH_RH", "SnowpackAdvanced", thresh_rh);

	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	sn_dt = M_TO_S(calculation_step_length);
	meteo_step_length = cfg.get("METEO_STEP_LENGTH", "Snowpack");

	//Defines whether a multiband model is used for short wave radiation extinction
	cfg.getValue("MULTISTREAM", "SnowpackAdvanced", multistream);

	//Defines whether joining elements will be considered at all
	cfg.getValue("JOIN_ELEMENTS", "SnowpackAdvanced", join_elements);

	//Which boundary condition to use
	string boundary_condition; cfg.getValue("SURFACECODE", "SnowpackAdvanced", boundary_condition);
	if (boundary_condition == "NEUMANN_BC") {
		surfaceCode = Snowpack::NEUMANN_BC;
	} else if (boundary_condition == "DIRICHLET_BC") {
		surfaceCode = Snowpack::DIRICHLET_BC;
	} else {
		throw InvalidArgumentException("The SURFACECODE " + boundary_condition + " is unknown", AT);
	}

	//Warning is issued if snow tempeartures are out of bonds, that is, crazy
	cfg.getValue("T_CRAZY_MIN", "SnowpackAdvanced", t_crazy_min);
	cfg.getValue("T_CRAZY_MAX", "SnowpackAdvanced", t_crazy_max);

/** @brief Initial new snow parameters, see computeSnowFall()
 * - that rg and rb are equal to 0.5*gsz and 0.5*bsz, respectively. Both given in millimetres
 * - If VW_DENDRICITY is set, new snow dendricity is f(vw)
 * - BOND_FACTOR_RH new snow bonds get stronger for average winds >= SnLaws::event_wind_lowlim and
 *   mean relative humidity >= rh_lowlim
 */
	if (variant == "ANTARCTICA") {
		new_snow_dd = 0.5;
		new_snow_sp = 0.75;
		new_snow_dd_wind = 0.15;
		new_snow_sp_wind = 1.0;
		vw_dendricity = false;
		rh_lowlim = 0.7;
		bond_factor_rh = 3.0;
	} else {
		new_snow_dd = 1.0;
		new_snow_sp = 0.5;
		new_snow_dd_wind = 0.5;
		new_snow_sp_wind = 0.75;
		vw_dendricity = true;
		rh_lowlim = 1.0;
		bond_factor_rh = 1.0;
	}

	cfg.getValue("NEW_SNOW_GRAIN_RAD", "SnowpackAdvanced", new_snow_grain_rad);
	new_snow_bond_rad = 0.25 * new_snow_grain_rad;

	/**
	 * @name Thresholds for surface hoar formation and burial
	 * NOTE that the value of the parameter ROUGHNESS_LENGTH in CONSTANTS_User.INI is critical for surface hoar formation,
	 * particularly for Dirichlet boundary conditions. Value should be < 1 mm. Other considerations favor larger values.
	 * - 0.0007 m : original calibration with the 98/99 data set \n
	 * - 0.002  m : favored operational value with Dirichlet bc
	 */
	//Density of BURIED surface hoar (kg m-3), default: 125./ Antarctica: 200.
	cfg.getValue("HOAR_DENSITY_BURIED", "SnowpackAdvanced", hoar_density_buried);
	//Density of surface hoar (-> hoar index of surface node) (kg m-3)
	cfg.getValue("HOAR_DENSITY_SURF", "SnowpackAdvanced", hoar_density_surf);

	//Minimum surface hoar size to be buried (mm). Increased by 50% for Dirichlet bc.
	cfg.getValue("HOAR_MIN_SIZE_BURIED", "SnowpackAdvanced", hoar_min_size_buried);

}

/**
 * @brief Return rain/snow temperature threshold that Snowpack uses
 * @return rain/snow threshold temperature (K)
 */
double Snowpack::getThreshRain() const {
	return C_TO_K( thresh_rain );
}

/**
 * @brief This routine evaluates the element stiffness matrix and right hand side vector for the
 * creep solution process. That is, the routine evaluates [Ke], {Fc} (which are later placed
 * on the right hand side), {Fi}, the internal forces and finally the external forces {Fe}
 * meaning, the self-weight of the snowpack.
 * NOTE: For now we are assuming that the density remains CONSTANT over the time step. At the
 * end of the time step the density is updated. This improves both STABILITY and ACCURACY
 * @param *Edata
 * @param dt time step (s)
 * @param cos_sl Cosine of slope angle
 * @param Zn Height of element nodes (m)
 * @param Un Temperature of Nodes (K)
 * @param Se Element stiffness matrix
 * @param Fc Creep forces
 * @param Fi Internal forces
 * @param Fe Element right hand side vector
 * @return false on error, true if no error occurred
 */
bool Snowpack::compSnowForces(ElementData *Edata,  double dt, double cos_sl, double Zn[ N_OF_INCIDENCES ], double Un[ N_OF_INCIDENCES ], double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ], double Fc[ N_OF_INCIDENCES ], double Fi[ N_OF_INCIDENCES ], double Fe[ N_OF_INCIDENCES ])
{
	double L, L0; // Length, initial length
	double dVol;  // Change in volume
	double D;     // Modulus of elasticity
	double E;     // Strain
	double S;     // Stress
	double Sc;    // Psuedo Creep Stress
	double ddE;   // Change in axial strian, volumetric strain

	// First compute the new length and from the new length the total strain
	dVol = Edata->L;
	L0 = Edata->L0; L = Edata->L = Zn[1] + Un[1] - Zn[0] - Un[0];
	dVol /= L;
	if (L <= 0.) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Element length < 0.0!\n L0=%lf L=%lf Zn[0]=%lf Zn[1]=%lf Un[0]=%lf Un[1]=%lf,", L0, L, Zn[0], Zn[1], Un[0], Un[1]);
		return false;
	}
	// Compute the Natural Strain // Green Lagrange Strain
	E = Edata->E = log(L / L0);
	if (!(E >= -100. && E <= 1.e-5)) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "In strain E (Memory?)");
		return false;
	}
	/*
	 * The volume change conserves ice mass, i.e. the volumetric air content is
	 * changed and the volumetric ice and water contents remain the same.
	 * What is required is not the total volume change -- rather the increment
	 * in volume change. Otherwise cannot update the volumetric components (which
	 * might be changed elsewhere) correctly. Update density.
	*/
	ddE = (E - Edata->dE);
	Edata->dE = E;
	if (ddE <= -1.5) {
		ddE = -1.5;
		prn_msg(__FILE__, __LINE__, "wrn", Date(), "Time step is too Large!!!");
	}
	// Compute the volumetric components
	Edata->theta[ICE] = MIN (0.999999, Edata->theta[ICE] * dVol);
	Edata->theta[WATER] = MIN (0.999999, Edata->theta[WATER] * dVol);
	Edata->theta[AIR] = 1.0 - Edata->theta[ICE] - Edata->theta[WATER] - Edata->theta[SOIL];
	if (!(Edata->theta[AIR] <= 1.0 && Edata->theta[AIR] >= -0.05)) {
		prn_msg(__FILE__, __LINE__, "msg+", Date(), "ERROR AIR: %e (ddE=%e)", Edata->theta[AIR], ddE);
		prn_msg(__FILE__, __LINE__, "msg", Date(), "ELEMENT SIZE: L0=%e L=%e", Edata->L0, Edata->L);
		prn_msg(__FILE__, __LINE__, "msg", Date(), "DENSITY: rho=%e", Edata->Rho);
		prn_msg(__FILE__, __LINE__, "msg", Date(), "ThetaICE: ti=%e", Edata->theta[ICE]);
		prn_msg(__FILE__, __LINE__, "msg", Date(), "ThetaWATER: ti=%e", Edata->theta[WATER]);
		return false;
	}
	// Compute the self weight of the element (with the present mass)
	Fe[0] = Fe[1] = -(Edata->M * Constants::g * cos_sl) / 2.;
	// Compute the elastic strain and stress
	Edata->Ee = Edata->E - Edata->Ev;
	D = Edata->snowElasticity();
	S = Edata->S = D * Edata->Ee;
	// Compute the internal forces
	Fi[0] = -S;
	Fi[1] = S;
	// Compute the pseudo increment in creep stress
	Sc = D * Edata->EvDot * dt;
	// Compute the creep forces
	Fc[0] = -Sc;
	Fc[1] = Sc;
	/*
	 * Update the element force vector (NOTE: Assembled Fe[0..5]=0.0! if elastic
	 * istropic solution WITHOUT creep.)
	*/
	Fe[0] = Fe[0]+Fc[0];
	Fc[0] = Fe[0]-Fi[0];
	Fe[1] = Fe[1]+Fc[1];
	Fc[1] = Fe[1]-Fi[1];
	// Compute the stiffness matrix
	Se[0][0] = D/L ;  Se[0][1] = -D/L ;
	Se[1][0] =-D/L ;  Se[1][1] = D/L ;

	return true;
}

/**
 * @brief Snow creep
 * -# The Thing ain't settling any more in case of ice, soil or water only
 * -# Enhanced densification for wind slabs of Metamorphism::wind_slab_depth (m); see also mm_Metamorphism()
 *    dry snow AND strong wind AND near the surface => enhance densification \n
 * -# Normal densification
 * -# Empirism for surface hoar. Two different rates are used for either large or small SH.
 *    Implemented by Sascha Bellaire on 28.11.2006.
 * @todo name parameters for the computation of CDot
 * @param Xdata
 * @param Mdata
 */
void Snowpack::compSnowCreep(const CurrentMeteo& Mdata, SnowStation& Xdata)
{
	unsigned int e;         // Element counter
	double L0, dL, cH_old;  // Element length and change of length
	bool prn_WRN = false;

	unsigned int nN = Xdata.getNumberOfNodes();
	if (nN == (Xdata.SoilNode + 1))
		return;

	vector<NodeData>& NDS = Xdata.Ndata;
	unsigned int nE = Xdata.getNumberOfElements();
	vector<ElementData>& EMS = Xdata.Edata;
	e = nE;
	double SigC = 0.; // Cauchy stress
	double Sig0 = 0.; // "Sintering" stress
	const double SigC_fac = Constants::g * cos(DEG_TO_RAD(Xdata.meta.getSlopeAngle()));
	while (e-- > 0) {
		double oldStress = EMS[e].C;
		double age = MAX(0., Mdata.date.getJulianDate() - EMS[e].depositionDate.getJulianDate());
		if (e < nE-1)
			SigC -= (EMS[e+1].M / 2.) * SigC_fac;
		EMS[e].C = SigC -= (EMS[e].M / 2.) * SigC_fac;
		if (EMS[e].CDot / SigC > 0.05) {
			EMS[e].CDot *= exp(-0.037 * S_TO_D(sn_dt));
		} else {
			EMS[e].CDot = 0.;
		}
		if ((e < nE-1) && (age > Constants::eps)) {
			if ((SigC - oldStress) < 0.)
				EMS[e].CDot += (SigC - oldStress);
		}
	}

	for (e = Xdata.SoilNode; e < nE; e++) {
		double eta = SnLaws::smallest_viscosity; // snow viscosity
		if (EMS[e].Rho > 910. ||  EMS[e].theta[SOIL] > 0. || EMS[e].theta[ICE] < Constants::eps) {
			EMS[e].k[SETTLEMENT] = eta = 1.0e99;
		} else {
			EMS[e].k[SETTLEMENT] = eta = SnLaws::compSnowViscosity(variant, viscosity_model, EMS[e], Mdata.date);
			if (!(eta > 0.01 * SnLaws::smallest_viscosity && eta <= 1.e11 * SnLaws::smallest_viscosity)
			        && (EMS[e].theta[ICE] > 2. * Snowpack::min_ice_content) && (EMS[e].theta[ICE] < 0.6)) {
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date,
				          "Viscosity=%e out of range! e=%d nE=%d rg=%lf rb=%lf dd=%lf sp=%lf theta_i=%lf theta_w=%lf",
				            eta, e, nE, EMS[e].rg, EMS[e].rb, EMS[e].dd, EMS[e].sp,
				              EMS[e].theta[ICE], EMS[e].theta[WATER]);
			}
			if (eta < SnLaws::smallest_viscosity) {
				if (prn_WRN) //HACK
					prn_msg(__FILE__, __LINE__, "wrn", Mdata.date,
					        "Viscosity=%e reset to SMALLEST_VISCOSITY! e=%d nE=%d", eta, e, nE);
				EMS[e].k[SETTLEMENT] = eta = SnLaws::smallest_viscosity;
			}
		}
		Sig0 = SnLaws::compLoadingRateStress(viscosity_model, EMS[e], Mdata.date);
		L0 = EMS[e].L;

		if (!(EMS[e].mk%100 == 3)) { //ALL except SH
			double wind_slab=1.;
			double dz = NDS[nE].z - NDS[e].z;
			double dv = Mdata.vw - Metamorphism::wind_slab_vw;
			if ((EMS[e].theta[WATER] < 0.01)
			      && (Mdata.vw > Metamorphism::wind_slab_vw)
			        && ((dz < Metamorphism::wind_slab_depth) || (e == nE-1))) {
				if (variant == "ANTARCTICA") {
					// fits original parameterization at Metamorphism::wind_slab_vw + 0.6 m/s
					wind_slab += 2.7*Metamorphism::wind_slab_enhance
					                 * dv*dv*dv * (1. - dz / (1.25 * Metamorphism::wind_slab_depth));
				} else {
					// original parameterization by Lehning
					wind_slab += Metamorphism::wind_slab_enhance * dv;
				}
			}
			EMS[e].EvDot = wind_slab * (EMS[e].C + Sig0) / eta;
			dL = L0 * sn_dt * EMS[e].EvDot;
			if ((L0 + dL) < minimum_l_element)
				dL = MIN(0., minimum_l_element - L0);
		} else { //SH
			if (NDS[e+1].hoar > 0.006) { // TODO Large initial size, i.e., deposited hoar mass/HOAR_DENSITY_BURIED ??
				if ((Mdata.date.getJulianDate() - EMS[e].depositionDate.getJulianDate()) < 21.)
					dL = MM_TO_M(-0.391 * S_TO_D(sn_dt));
				else
					dL = MM_TO_M(-0.0807 * S_TO_D(sn_dt));
			} else {
				dL = MM_TO_M(-0.1107 * S_TO_D(sn_dt));
			}
			EMS[e].EvDot = dL / (L0 * sn_dt);
			if ((L0 + dL) < 0. )
				dL = 0.;
		}

		if (variant == "CALIBRATION") {
			NDS[e].f = (Sig0 - EMS[e].EDot) / eta; // sigMetamo
			EMS[e].EDot /= eta;                    // sigReac
			NDS[e].udot = EMS[e].C / eta;          // Deformation due to load alone
			EMS[e].S = EMS[e].CDot / EMS[e].C;     // ratio loadrate to load addLoad to load
		}

		EMS[e].theta[WATER] *= L0 / (L0 + dL);
		EMS[e].theta[ICE]   *= L0 / (L0 + dL);
		EMS[e].L0 = EMS[e].L = (L0 + dL);
		NDS[e+1].z = NDS[e].z + EMS[e].L;
		while (EMS[e].theta[ICE] + EMS[e].theta[WATER] + EMS[e].theta[SOIL] > 0.99) {
			EMS[e].theta[ICE] *= 0.99;
			EMS[e].theta[WATER] *= 0.99;
			EMS[e].M = EMS[e].L0 * ((EMS[e].theta[ICE] * Constants::density_ice)
			                           + (EMS[e].theta[WATER] * Constants::density_water));
		}
		EMS[e].theta[AIR] = 1.0 - EMS[e].theta[WATER] - EMS[e].theta[ICE] - EMS[e].theta[SOIL];
		EMS[e].Rho = (EMS[e].theta[ICE] * Constants::density_ice) + (EMS[e].theta[WATER]
		                *Constants::density_water) + (EMS[e].theta[SOIL]
		                  * EMS[e].soil[SOIL_RHO]);
		if (! (EMS[e].Rho > 0. && EMS[e].Rho <= Constants::max_rho)) {
			prn_msg(__FILE__, __LINE__, "err", Date(),
			          "Volume contents: e=%d nE=%d rho=%lf ice=%lf wat=%lf air=%le",
			            e, nE, EMS[e].Rho, EMS[e].theta[ICE], EMS[e].theta[WATER], EMS[e].theta[AIR]);
			throw IOException("Runtime Error in compSnowCreep()", AT);
		}
	} // end for settling-loop
	cH_old = Xdata.cH;
	Xdata.cH = NDS[nN-1].z + NDS[nN-1].u;
	Xdata.mH -= (cH_old - Xdata.cH);
}

/**
 * @brief Computes the element stiffness matrix and right hand side vector for a fully implicit time integration scheme \n
 * The matrices that must be evaluated are : \n
 * - [Se] = [Ce]/dt + [Ke] :  [Ce] is the element heat capacity matrix and [Ke] is the
 *                            element conductivity matrix.
 * - {Fe} = {Q} - [Ke]*{T0} : {Fe} is the element right-hand side vector containing the element heat flux. \n
 *                            {Q} arises from shortwave radiation -- the other heat fluxes are treated separately
 *                            when determining the meteo parameters.
 * @param *Edata
 * @param dt Calculation time step length (s)
 * @param dvdz Wind pumping velocity gradient (s-1)
 * @param T0 Initial nodal temperatures (K)
 * @param Tn Iterated nodal temperatures (K)
 * @param Se Element heat capacitity (stiffness) matrix
 * @param Fe Element right hand side vector
 * @param *SubSurfaceMelt Check for subsurface melt
 * @param *SubSurfaceFrze Check for subsurface freezing
 * @param VaporEnhance Vapor transport enhancement factor
 * @return false on error, true if no error occurred
 */
bool Snowpack::sn_ElementKtMatrix(ElementData *Edata, double dt, double dvdz, double T0[ N_OF_INCIDENCES ], double Tn[ N_OF_INCIDENCES ], double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ], double Fe[ N_OF_INCIDENCES ], char *SubSurfaceMelt, char *SubSurfaceFrze, double VaporEnhance)
{
	double k, c;    // conductivity and heat capacity
	double Keff;    // the effective thermal conductivity

	// Gather the data which is required to compute the element stiffness and force vector
	Fe[0] = Fe[1] = 0.0;
	if (Edata->L < 0.0) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Negative length L=%e", Edata->L);
		return false;
	}

	// Phase Change Check
	// NOTE Should it not preferably be done in PhaseChange.cc ??
	if ((Tn[0] > Constants::melting_tk || Tn[1] > Constants::melting_tk) && Edata->theta[ICE] > Constants::eps2)
		*SubSurfaceMelt = true;
	if ((Tn[0] < Constants::freezing_tk || Tn[1] < Constants::freezing_tk) && Edata->theta[WATER] > Constants::eps2)
		*SubSurfaceFrze = true;
	if (Edata->theta[WATER] > Constants::eps2 && Edata->theta[SOIL] < Constants::eps2)
		 T0[0] = T0[1] = Constants::melting_tk;

	// Find the conductivity of the element
	if (Edata->theta[SOIL] > 0.0) {
		Keff = SnLaws::compSoilThermalConductivity(*Edata, dvdz);
	} else if (Edata->theta[ICE] > 0.55 || Edata->theta[ICE] < min_ice_content) {
		Keff = Edata->theta[AIR] * Constants::conductivity_air + Edata->theta[ICE] * Constants::conductivity_ice +
		           Edata->theta[WATER] * Constants::conductivity_water + Edata->theta[SOIL] * Edata->soil[SOIL_K];
	} else {
		Keff = SnLaws::compSnowThermalConductivity(*Edata, dvdz);
	}
	// mimics effect of vapour transport if liquid water present in snowpack
	Keff *= VaporEnhance;
	Edata->k[TEMPERATURE] = Keff;
	k = Keff/Edata->L;   // Divide by the length to save from doing it during the matrix operations

	// Evaluate the stiffness matrix
	Se[0][0] = Se[1][1] = k; Se[0][1] = Se[1][0] = -k;
	// Heat the element via short-wave radiation
	if (Edata->sw_abs < 0.0) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "NEGATIVE Shortwave Radiation %e", Edata->sw_abs);
		return false;
	}
	Fe[1] += Edata->sw_abs;

	/*
	 * Phase change ennergy:
	 * double Qpc = Edata->Qmf * Edata->L / 4.0; // original: 2.0 Changed by ML: 15.02.2002 as a compromise between Fz and Pb
	 * Heat / Cool element via phasechanges (Edata->Qmf is computed in phase change module)
	 * Why should we do that ????
	 * Phase Change already USED this energy in the previous time step!!!
	 * ==> thus commented out 06.10.2001 Fz
	*/
	Fe[0] += 0.; // Qpc;
	Fe[1] += 0.; // Qpc;

	// Add the implicit time integration term to the right hand side
	Fe[0] -= (Se[0][0] * T0[0] + Se[0][1] * T0[1]);
	Fe[1] -= (Se[1][0] * T0[0] + Se[1][1] * T0[1]);

	// Now add the heat capacitity matrix
	Edata->heatCapacity();
	c = Edata->c[TEMPERATURE] * Edata->L * Edata->Rho / (6. * dt);
	Se[0][0] += 2. * c;
	Se[1][1] += 2. * c;
	Se[0][1] += c;
	Se[1][0] += c;

	return true;
}

/**
 * @brief Update all energy fluxes but solar irradiance \n
 * @par The fluxes are updated before entering the temperature routine sn_SnowTemperature.
 * These fluxes are then used to impose Neumann boundary conditions.
 * This splitting ensures that these initial fluxes are used for each iteration in the
 * semi-explicit solution while they are not altered by the implicit solution.
 *
 * - qs = alpha*(Ta - Tss) :     Sensible heat transfer \n
 *                               alpha is a coefficient based on the famous Monin - Obukhov theory.
 *                               It is strongly dependent on the wind speed. Ta and Tss are
 *                               air and surface temperatures (K), respectively
 * - ql = beta*(eA - eS) :       Latent heat transfer \n
 *                               eA and eS are the vapor pressures of air and snow, respectively. \n
 *                               Includes consideration of soil (one active element).
 * - qR = gamma*(T_air - Tss) :  Energy convected by rain
 * - q_lw = SE*SB*(e*Ta^4 - Ts^4) : Long wave radiation heat transfer \n
 *                               For the implicit solution,this will be treated as a convective
 *                               boundary condition, similar to the sensible heat exchange.
 *                               ES = emissivity_snow, SB = stefan-boltzmann constant.
 * @param *Mdata
 * @param *Xdata
 * @param *Bdata
 */
void Snowpack::updateMeteoHeatFluxes(const CurrentMeteo& Mdata, SnowStation& Xdata,
                                     BoundCond& Bdata)
{
	double alpha, gamma;		// The exchange coefficients
	double lw_in, lw_out;		// Incoming and outgoing long wave radiation
	const double& T_air = Mdata.ta;
	const double& Tss = Xdata.Ndata[Xdata.getNumberOfNodes()-1].T;

	alpha = SnLaws::compSensibleHeatCoefficient(Mdata, Xdata, height_of_meteo_values);
	Bdata.qs = MIN (350., MAX (-350., alpha * (T_air - Tss)));

	Bdata.ql = SnLaws::compLatentHeat_Rh(Mdata, Xdata, height_of_meteo_values);
	if ((Xdata.getNumberOfElements() > 0)
	        /*&& (Xdata.Edata[Xdata.getNumberOfElements()-1].theta[ICE] >= min_ice_content)*/) { //HACK: how should we handle large fluxes?
		Bdata.ql = MIN (250., MAX (-250., Bdata.ql));
	}

	if (T_air >= C_TO_K(thresh_rain)) {
		gamma = (Mdata.hnw / sn_dt) * Constants::specific_heat_water;
		Bdata.qr = gamma * (T_air - Tss);
	} else {
		Bdata.qr = 0.;
	}

	lw_in  = Constants::emissivity_snow * Constants::stefan_boltzmann * Mdata.ea * T_air*T_air*T_air*T_air;
	lw_out = Constants::emissivity_snow * Constants::stefan_boltzmann * Tss*Tss*Tss*Tss;
	Bdata.lw_out = lw_out;
	Bdata.lw_net = lw_in - lw_out;

}

/**
 * @brief Imposes Neumann boundary conditions at the surface. \n
 * The initial fluxes are computed in sn_MeteoHeatFluxes. \n
 * This splitting ensures that the initial fluxes are used for each iteration in the
 * semi-explicit solution while they are not altered by the implicit solution. \n
 * Note that long wave radiation heat transfer is treated as a convection boundary condition
 * for the implicit solution, similar to the sensible heat exchange: \n
 *            q_lw = delta*(e*Ta - T_iter) \n
 * where delta is a function of both Ta and T_iter and is computed by SnLaws::compLWRadCoefficient
 * @param Mdata
 * @param Bdata
 * @param Xdata
 * @param T_snow Initial (snow) surface temperature (K)
 * @param T_iter Iterated (snow) surface temperature (K)
 * @param Se Element stiffness matrix
 * @param Fe Element right hand side vector
 */
void Snowpack::neumannBoundaryConditions(const CurrentMeteo& Mdata, BoundCond& Bdata, const SnowStation& Xdata,
                                         const double& T_snow, const double& T_iter,
                                         double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ],
                                         double Fe[ N_OF_INCIDENCES ])
{
	double alpha, gamma, delta;  //  The heat exchange coefficients
	double S_eb;
	double T_air = Mdata.ta;

	// First zero out the interiour node contribution
	Se[0][0] = Se[0][1] = Se[1][0] = Se[1][1] = Fe[0] = Fe[1] = 0.0;

	/*
	 * Now branch between phase change cases (semi-explicit treatment) and dry snowpack dynamics
	 * (implicit treatment)
	*/
	if ((Xdata.Edata[Xdata.getNumberOfElements()-1].theta[WATER] > 0.) && (T_iter != T_snow)) {
		// Semi-explicit
		// Latent heat transfer and net longwave radiation
		Fe[1] += Bdata.ql + Bdata.lw_net;
		// Sensible heat transfer and advected rain energy
		S_eb = Bdata.qs + Bdata.qr;
		if ((T_air - T_iter) != 0.)
			alpha = MAX (0.0, S_eb / (T_air - T_iter));
		else
			alpha = 1.0;
		Se[1][1] += alpha;
		Fe[1] += alpha * T_air;
	} else { // Implicit
		// Sensible heat transfer: linear dependence on snow surface temperature
		alpha = SnLaws::compSensibleHeatCoefficient(Mdata, Xdata, height_of_meteo_values);
		Se[1][1] += alpha;
		Fe[1] +=alpha * T_air;
		// Latent heat transfer: NON-linear dependence on snow surface temperature
		//NOTE: should it not be linearized then?
		Fe[1] += Bdata.ql;
		// Advected rain energy: linear dependence on snow surface temperature
		if (T_air >= C_TO_K(thresh_rain)) {
			gamma = (Mdata.hnw / sn_dt) * Constants::specific_heat_water;
			Se[1][1] += gamma;
			Fe[1] += gamma * T_air;
		}
		// Net longwave radiation: NON-linear dependence on snow surface temperature
		delta = SnLaws::compLWRadCoefficient( T_iter, T_air, Mdata.ea);
		Se[1][1] += delta;
		Fe[1] += delta * pow( Mdata.ea, 0.25 ) * T_air;
	} // end else

	// Because of the implicit time integration, must subtract this term from the flux ....
	Fe[1] -= Se[1][1] * T_snow;
}

/**
 * @brief Imposes the Neumann boundary conditions at the bottom node. \n
 * Note that this can only be used if you have a deep enough soil,
 * because the heat flux is currently constant. For Ethan, we have had already a
 * version with upper and lower heat flux on a small snow sample and that worked perfectly.
 * @param flux Ground heat flux (W m-2)
 * @param T_snow Initial upper node temperature of bottom element (K)
 * @param Se Element stiffness matrix
 * @param Fe Element right hand side vector
 */
void Snowpack::neumannBoundaryConditionsSoil(const double& flux, const double& T_snow,
                                             double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ],
                                             double Fe[ N_OF_INCIDENCES ])
{
	double alpha;  //  The heat exchange coefficients
	double T_pseudo;

	// First zero out the interiour node contribution
	Se[0][0] = Se[0][1] = Se[1][0] = Se[1][1] = Fe[0] = Fe[1] = 0.0;

	// Use the numerical trick of an assumed temperature difference of 1 K over the boundary
	T_pseudo = T_snow - 1.;

	// For the implicit solution, we need to define a pseudo-exchange coefficient
	alpha = flux / (T_pseudo - T_snow);
	Se[1][1] += alpha;
	Fe[1] += alpha * T_pseudo;

	// Because of the implicit time integration, must subtract this term from the flux ....
	Fe[1] -= Se[1][1] * T_snow;
}

/**
 * @brief Computes the snow temperatures which are given by the following formula: \n
 * @par
 *             dT            d^2T
 *   rho*c(T)*----  -  k(T)*------   =  Q
 *             dt            dz^2
 * \n
 * with initial and Dirichlet and/or Neumann boundary conditions.
 *
 * T(z,t) = temperature (K);
 * rho = snow density (kg m-3); c = specific heat capacity (J K-1 kg-1); k = heat conductivity (W K-1 m-1); t = time (s);
 * \n
 * Note:  The equations are solved with a fully implicit time-integration scheme and the
 * system of finite element matrices are solved using a sparse matrix solver.
 * @param Xdata
 * @param Mdata
 * @param Bdata
 */
void Snowpack::compSnowTemperatures(SnowStation& Xdata, CurrentMeteo& Mdata, BoundCond& Bdata)
{
	size_t n, e;        // nodal and element counters
	int iteration;	// iteration counter (not really required)
	int NotConverged;	// = 1 if iteration not converged
	double MaxTDiff;	// maximum temperature difference for convergence
	double TDiff;		// temperature difference for convergence check
	double Big = Constants::big;	// big number (for DIRICHLET boundary conditions
	// Set the default solution routine convergence parameters
	int MaxItnTemp = 40;                         // maximum 40 iterations for temperature field
	double ConTolTemp = 0.01;	                   // solution convergence to within 0.01 degC

	int Ie[N_OF_INCIDENCES];                     // Element incidences
	double T0[N_OF_INCIDENCES];                  // Element nodal temperatures at t0
	double TN[N_OF_INCIDENCES];                  // Iterated element nodal temperatures
	double Se[N_OF_INCIDENCES][N_OF_INCIDENCES]; // Element stiffnes matrix
	double Fe[N_OF_INCIDENCES];                  // Element right hand side vector

	void *Kt;                                    // Avoids dereferencing the pointer
	double *U=NULL, *dU=NULL, *ddU=NULL;         // Solution vectors
	double I0;                                   // The net incoming shortwave irradiance
	double cAlb;                                 // Computed albedo
	double hs;                                   // Snow depth
	double d_pump, v_pump, dvdz;                 // Wind pumping parameters
	double VaporEnhance;                         // Enhancement factor for water vapor transport in wet snow

	// Dereference the pointers
	Kt = Xdata.Kt;
	vector<NodeData>& NDS = Xdata.Ndata;
	size_t nN = Xdata.getNumberOfNodes();
	vector<ElementData>& EMS = Xdata.Edata;
	size_t nE = Xdata.getNumberOfElements();

	// ABSORPTION OF SOLAR RADIATION WITHIN THE SNOWPACK
	// What snow depth should be used?
	if (enforce_measured_snow_heights && (Xdata.meta.getSlopeAngle() <= Constants::min_slope_angle))
		hs = Mdata.hs1;
	else
		hs = Xdata.cH;

	// Parameterized albedo (statistical model) including correct treatment of PLASTIC and WATER_LAYER
	if ((nE > Xdata.SoilNode)) { //Snow, glacier, ice, water, or plastic layer
		size_t eAlbedo = nE-1;
		size_t marker = EMS[eAlbedo].mk % 10;
		switch (marker) {
		case 9: // WATER_LAYER
			if (eAlbedo > Xdata.SoilNode)
				eAlbedo--;
		case 8: // Ice layer within the snowpack
			while ((eAlbedo > Xdata.SoilNode) && (marker == 8))
				eAlbedo--;
		default: // Snow, glacier ice, PLASTIC, or soil
			if (eAlbedo > Xdata.SoilNode && (EMS[eAlbedo].theta[SOIL] < Constants::eps2)) { // Snow, or glacier ice
				cAlb = SnLaws::parameterizedSnowAlbedo(fixed_albedo, EMS[eAlbedo], NDS[eAlbedo+1].T, Mdata);
			} else { // PLASTIC, or soil
				cAlb = Xdata.SoilAlb;
			}
			break;
		}
	} else { // Soil
		cAlb = Xdata.SoilAlb;
	}

	if (!(useCanopyModel && (Xdata.Cdata.height > 3.5))) {
		if (research_mode) { // Treatment of "No Snow" on the ground in research mode
			if ((hs < 0.02) || (NDS[nN-1].T > C_TO_K(3.5))
			                  || ((hs < 0.05) && (NDS[nN-1].T > C_TO_K(1.7)))) {
				cAlb = Xdata.SoilAlb;
			}
		}
		if (!ALPINE3D)
			cAlb = MAX(cAlb, Mdata.rswr / Constants::solcon);
		cAlb = MAX(Xdata.SoilAlb, MIN(0.95, cAlb));
	} else {
		cAlb = MAX(0.05, MIN(0.99, cAlb));
	}

	switch (sw_mode) {
	case 0: // need to assign rswr a correct value
		Mdata.rswr = Mdata.iswr * cAlb;
		break;
	case 1: // need to assign iswr a correct value
		Mdata.iswr = Mdata.rswr / cAlb;
		break;
		// ... while the ground is still snow covered according to HS measurements
	case 2: // use measured albedo ...
		if (Xdata.mAlbedo != Constants::undefined) {
			if ((!((Xdata.mAlbedo < 2.*Xdata.SoilAlb)
			        && ((Xdata.cH - Xdata.Ground) > 0.05))) && Xdata.mAlbedo <= 0.95)
				cAlb = Xdata.mAlbedo;
			else
				Mdata.rswr = Mdata.iswr * cAlb;
		} else {
			Mdata.rswr = Mdata.iswr = 0.;
		}
		break;
	default:
		prn_msg(__FILE__, __LINE__, "err", Mdata.date, "sw_mode = %d not implemented yet!", sw_mode);
		exit(EXIT_FAILURE);
		break;
	}
	Xdata.cAlbedo = cAlb; // Assign albedo to Xdata
	I0 = Mdata.iswr - Mdata.rswr; // Net irradiance perpendicular to slope
	if (I0 < 0.) {
		prn_msg(__FILE__, __LINE__, "err", Mdata.date, " iswr:%lf  rswr:%lf  cAlb:%lf", Mdata.iswr, Mdata.rswr, cAlb);
		exit(EXIT_FAILURE);
	}

	// Simple treatment of radiation absorption in snow: Beer-Lambert extinction (single or multiband).
	try {
		SnLaws::compShortWaveAbsorption(I0, useSoilLayers, multistream, Xdata);
	} catch(const exception&){
		prn_msg(__FILE__, __LINE__, "err", Mdata.date, "Runtime error in sn_SnowTemperature");
		throw;
	}

	// Take care of uppermost soil element
	if ((nE > Xdata.SoilNode+1) && (EMS[Xdata.SoilNode].sw_abs > EMS[Xdata.SoilNode+1].sw_abs)) {
		EMS[nE-1].sw_abs += (EMS[Xdata.SoilNode].sw_abs - EMS[Xdata.SoilNode+1].sw_abs);
		EMS[Xdata.SoilNode].sw_abs = EMS[Xdata.SoilNode+1].sw_abs;
	}

	// Set ground surface temperature
	if (nN == 1) {
		// No snow on the ground and !useSoilLayers
		if ((Mdata.ts0 > Constants::melting_tk) && ((Mdata.ts0 - Mdata.ta) > 10.))
			NDS[0].T = (Mdata.ts0 + Mdata.ta) / 2.;
		else
			NDS[0].T = Mdata.ts0;
		return;
	}

	if (Kt != 0)
		ds_Solve(ReleaseMatrixData, (MYTYPE*)Kt, 0);
	ds_Initialize(nN, 0, (MYTYPE**)&Kt);
	/*
	 * Define the structure of the matrix, i.e. its connectivity. For each element
	 * we compute the element incidences and pass the incidences to the solver.
	 * The solver assumes that the element incidences build a crique, i.e. the
	 * equations specified by the incidence set are all connected to each other.
	 * Initialize element data.
	*/
	for (e = 0; e < nE; e++) {
		int Nodes[2] = {e, e+1};
		ds_DefineConnectivity( (MYTYPE*)Kt, 2, Nodes , 1, 0 );
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
	ds_Solve(SymbolicFactorize, (MYTYPE*)Kt, 0);

	// Make sure that these vectors are always available for use ....
	errno=0;
	U=(double *) realloc(U, nN*sizeof(double));
	if (errno != 0) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "%s (allocating  solution vector U)", strerror(errno));
		throw IOException("Runtime error in sn_SnowTemperature", AT);
	}
	dU=(double *) realloc(dU, nN*sizeof(double));
	if (errno != 0) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "%s (allocating  solution vector dU)", strerror(errno));
		throw IOException("Runtime error in sn_SnowTemperature", AT);
	}
	ddU=(double *) realloc(ddU, nN*sizeof(double));
	if (errno != 0) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "%s (allocating  solution vector ddU)", strerror(errno));
		throw IOException("Runtime error in sn_SnowTemperature", AT);
	}

  // Make sure that the global data structures know where the pointers are for the next integration step after the reallocation ....
	Xdata.Kt = Kt;

	/*
	 * Set the temperature at the snowpack base to the prescribed value.
	 * NOTE if there is water in the base element then the base temperature MUST be 0 degC
	 * NOTE ts0 no longer set to 0 degC in Control.c
	*/
	if (!(useSoilLayers && soil_flux)) {
		if ((EMS[0].theta[ICE] >= min_ice_content)) {
			if ((EMS[0].theta[WATER] > 0.003)) {
				NDS[0].T = C_TO_K(0.0); /*Mdata.ts0 = C_TO_K(0.0);*/
			} else if (!useSoilLayers) {
				// To avoid temperatures above freezing while snow covered
				NDS[0].T = MIN(Mdata.ts0, C_TO_K(0.0));
			} else {
				NDS[0].T = Mdata.ts0;
			}
		} else {
			NDS[0].T = Mdata.ts0;
		}
	}
	if (surfaceCode == DIRICHLET_BC)
		NDS[nE].T = Mdata.tss;
	// Copy Temperature at time0 into First Iteration
	for (n = 0; n < nN; n++) {
		U[n] = NDS[n].T;
		dU[n] = 0.0;
		ddU[n] = 0.0;
		if (!(U[n] > 50. && U[n] < 500.)) {
			prn_msg(__FILE__, __LINE__, "err", Mdata.date, "Temperature out of bond at beginning of iteration!");
			prn_msg(__FILE__, __LINE__, "msg", Date(), "At snow node n=%d (nN=%d): T=%.2lf", n, nN, U[n]);
			free(U); free(dU); free(ddU);
			throw IOException("Runtime error in sn_SnowTemperature", AT);
		}
	}
	// Set the iteration counters, as well as the phase change boolean values
	iteration = 0;
	NotConverged = true;
	// Set the phase change booleans
	Xdata.SubSurfaceMelt = false;
	Xdata.SubSurfaceFrze = false;
	// Determine the displacement depth d_pump
	d_pump = SnLaws::compWindPumpingDisplacement(Xdata);

	// IMPLICIT INTEGRATION LOOP
	do {
		iteration++;
		// Reset the matrix data and zero out all the increment vectors
		ds_Solve(ResetMatrixData, (MYTYPE*)Kt, 0);
		for (n = 0; n < nN; n++) {
			ddU[n] = dU[n];
			dU[n] = 0.0;
		}
		// Reinialize the wind pumping speed at the surface
		if (nE > Xdata.SoilNode || SnLaws::wind_pump_soil)
			v_pump = SnLaws::compWindPumpingVelocity(Mdata, d_pump);
		else
			v_pump = 0.0;
		// Assemble matrix
		e = nE;
		while (e-- > 0) {
			EL_INCID( e, Ie );
			EL_TEMP( Ie, T0, TN, NDS, U );
			// Update the wind pumping velocity gradient
			dvdz = SnLaws::compWindGradientSnow(EMS[e], v_pump);
			// Update the water vapor transport enhancement factor
			VaporEnhance = MAX (1., SnLaws::compEnhanceWaterVaporTransportSnow(Xdata, e));
			// Now fill the matrix
			if (!sn_ElementKtMatrix(&EMS[e], sn_dt, dvdz, T0, TN, Se, Fe, &Xdata.SubSurfaceMelt,
			                        &Xdata.SubSurfaceFrze, VaporEnhance)) {
				prn_msg(__FILE__, __LINE__, "msg+", Mdata.date, "Error in sn_ElementKtMatrix @ element %d:", e);
				for (n = 0; n < nN; n++)
					fprintf(stdout, "U[%u]=%e K\n", (unsigned int)n, U[n]);
				free(U); free(dU); free(ddU);
				throw IOException("Runtime error in sn_SnowTemperature", AT);
			}
			ds_AssembleMatrix( (MYTYPE*)Kt, 2, Ie, 2,  (double*) Se );
			EL_RGT_ASSEM( dU, Ie, Fe );
		}

		/*
		 * The top element is special in that it handles the entire meteo conditions
		 * Several terms must be added to the global stiffness matrix Kt and flux
		 * right-hand side vector dU. Note:  Shortwave radiation --- since it is a body
		 * or volumetric force --- is computed in sn_ElementKtMatrix().
		*/
		if (surfaceCode == NEUMANN_BC) {
			EL_INCID(nE-1, Ie);
			EL_TEMP(Ie, T0, TN, NDS, U);
			neumannBoundaryConditions(Mdata, Bdata, Xdata, T0[1], TN[1], Se, Fe);
			ds_AssembleMatrix( (MYTYPE*)Kt, 2, Ie, 2,  (double*) Se );
			EL_RGT_ASSEM( dU, Ie, Fe );
		}
		if (surfaceCode == DIRICHLET_BC) {
			Ie[0] = nE;
			ds_AssembleMatrix( (MYTYPE*)Kt, 1, Ie, 1, &Big );
		}
		// Bottom node
		if (useSoilLayers && soil_flux) {
			// Neumann BC at bottom: The lower boundary is now a heat flux -- put the heat flux in dU[0]
			EL_INCID(0, Ie);
			EL_TEMP(Ie, T0, TN, NDS, U);
			Bdata.qg = geo_heat;
			neumannBoundaryConditionsSoil(Bdata.qg, T0[1], Se, Fe);
			ds_AssembleMatrix((MYTYPE*)Kt, 2, Ie, 2, (double*) Se);
			EL_RGT_ASSEM(dU, Ie, Fe);
		} else if ((Xdata.getNumberOfElements() < 3) && (Xdata.Edata[0].theta[WATER] >= 0.9
		               * Xdata.Edata[0].snowResidualWaterContent())) {
			dU[0] = 0.;
		} else {
			// Dirichlet BC at bottom: prescribed temperature value
			// NOTE Insert Big at this location to hold the temperature constant at the prescribed value.
			Ie[0] = 0;
			ds_AssembleMatrix( (MYTYPE*)Kt, 1, Ie, 1, &Big );
		}

		/*
		 * Solve the linear system of equation. The te_F vector is used first as right-
		 * hand-side vector for the linear system. The solver stores in this vector
		 * the solution of the system of equations, the new temperature.
		*/
		ds_Solve( ComputeSolution, (MYTYPE*)Kt, dU );
		// Update the solution vectors and check for convergence
		for (n = 0; n < nN; n++)
			ddU[n] = dU[n] - ddU[n];
		MaxTDiff = fabs(ddU[0]);
		for (n = 1; n < nN; n++) {
			TDiff = fabs(ddU[n]);
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
		if (U[nE] + ddU[nE] > Constants::melting_tk || EMS[nE-1].theta[WATER] > 0.) {
			ConTolTemp = 0.007;
			MaxItnTemp = 200; // NOTE originally 100;
		}
		NotConverged = ( MaxTDiff > ConTolTemp );
		if (iteration > MaxItnTemp) {
			prn_msg(__FILE__, __LINE__, "err", Mdata.date,
			          "Temperature did not converge (azi=%.0lf, slope=%.0lf)!",
			            Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
			prn_msg(__FILE__, __LINE__, "msg", Date(),
			          "%d iterations > MaxItnTemp=%d; ConTolTemp=%.4lf; nN=%d",
			            iteration, MaxItnTemp, ConTolTemp, nN);
			for (n = 0; n < nN; n++) {
				if (n > 0)
					prn_msg(__FILE__, __LINE__, "msg-", Date(),
					          "U[%03d]:%6.1lf K, ddU:%8.4lf K;  NDS.T(t-1)=%6.1lf K; EMS[n-1].th_w(t-1)=%.5lf",
					            n, U[n], ddU[n], NDS[n].T, EMS[n-1].theta[WATER]);
				else
					prn_msg(__FILE__, __LINE__, "msg-", Date(),
					          "U[%03d]:%6.1lf K, ddU:%8.4lf K;  NDS.T(t-1)=%6.1lf K;",
					            n, U[n], ddU[n], NDS[n].T);
			}
			prn_msg(__FILE__, __LINE__, "msg", Date(),
			          "Latent: %lf  Sensible: %lf  Rain: %lf  NetLong:%lf  NetShort: %lf",
			            Bdata.ql, Bdata.qs, Bdata.qr, Bdata.lw_net, I0);
			free(U); free(dU); free(ddU);
			throw IOException("Runtime error in sn_SnowTemperature", AT);
		}
		for (n = 0; n < nN; n++)
			U[n] += ddU[ n ];
	} while ( NotConverged ); // end Convergence Loop

	unsigned int crazy = 0;
	bool crazy_wrn = true;
	for (n = 0; n < nN; n++) {
		if ((U[n] > t_crazy_min) && (U[n] < t_crazy_max)) {
			NDS[n].T = U[n];
		} else { // Correct for - hopefully transient - crazy temperatures!
			NDS[n].T = 0.5*(U[n] + NDS[n].T);
			if (!ALPINE3D) {
				if (!crazy && (nN > Xdata.SoilNode + 3)) {
					prn_msg(__FILE__, __LINE__, "wrn", Mdata.date, "Crazy temperature(s)!");
					prn_msg(__FILE__, __LINE__, "msg-", Date(),
					        "T <= %5.1lf OR T >= %5.1lf; nN=%d; cH=%6.3lf m; azi=%.0lf, slope=%.0lf",
					        t_crazy_min, t_crazy_max, nN, Xdata.cH,
					        Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
					crazy_wrn = false;
				}
				if (useSoilLayers && (n <= Xdata.SoilNode)) {
					prn_msg(__FILE__, __LINE__, "msg-", Date(),
					        "SOIL node %3d: U(t)=%6.1lf  NDS.T(t-1)=%6.1lf K; EMS[n-1].th_w(t-1)=%.5lf",
					        n, U[n], NDS[n].T, EMS[n-1].theta[WATER]);
				} else if ((n > Xdata.SoilNode) && (nN > Xdata.SoilNode + 3)) {
					prn_msg(__FILE__, __LINE__, "msg-", Date(),
					        "SNOW node %3d: U(t)=%6.1lf  NDS.T(t-1)=%6.1lf EMS[n-1].th_w(t-1)=%.5lf",
					        n, U[n], NDS[n].T, EMS[n-1].theta[WATER]);
				}
			}
			crazy++;
		}
	}
	if (crazy) {
		if (crazy_wrn)
			prn_msg(__FILE__, __LINE__, "wrn", Mdata.date,
			        "%d crazy node(s) from total %d! azi=%.0lf, slope=%.0lf",
			        crazy, nN, Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
		else
			prn_msg(__FILE__, __LINE__, "msg-", Date(),
			        "%d crazy node(s) from total %d!",
			        crazy, nN);
		prn_msg(__FILE__, __LINE__, "msg-", Date(),
		        "Latent: %lf  Sensible: %lf  Rain: %lf  NetLong:%lf  NetShort: %lf",
		        Bdata.ql, Bdata.qs, Bdata.qr, Bdata.lw_net, I0);
	}

	for (e = 0; e < nE; e++) {
		EMS[e].Te = (NDS[e].T + NDS[e+1].T) / 2.0;
		EMS[e].heatCapacity();
		EMS[e].gradT = -(NDS[e].T - NDS[e+1].T) / EMS[e].L;
	}
	free(U); free(dU); free(ddU);
}

/**
 * @brief Determines whether new snow elements are added on top of the snowpack
 * - If enforce_measured_snow_heights=0 (research mode), new snow height corresponding to cumulated
 *   new snow water equivalent cumu_hnw must be greater than HEIGHT_NEW_ELEM to be added to Xdata->mH
 * - In case of virtual slopes, uses new snow depth and density from either flat field or luv slope
 * - The first thing is to compute the height of each element in the snow layer. For now,
 *   instead of trying to find an optimal number of elements, we will define the number of new
 *   elements as a constant nNnew. The constant is related to HEIGHT_NEW_ELEM, which is set in
 *   qr_ReadParameter(). The smaller HEIGHT_NEW_ELEM, the more elements we will use.
 * @param Mdata Meteorological data
 * @param Xdata Snow cover data
 * @param cumu_hnw cumulated amount of precipitation (kg m-2)
 */
void Snowpack::compSnowFall(const CurrentMeteo& Mdata, SnowStation& Xdata, double& cumu_hnw,
                            SurfaceFluxes& Sdata)
{
	bool force_layer = false, snow_fall = false, snowed_in = false;
	unsigned int nNewE, nHoarE;        // The number of elements to be added
	unsigned int nOldE, nE, nOldN, nN; // Old and new numbers of elements and nodes
	unsigned int e, n;                 // Element and node counters
	double z0;                         // Used to determine the z-location of new snowfall nodes
	double Ln;                         // Original new snow layer element length

	double rho_hn, t_surf, hn, hoar;   // New snow data
	double cos_sl, L0, dL, Theta0;     // Local values

	//Threshold for detection of the first snow fall on soil/canopy (grass/snow detection)
	const double TSS_threshold24=-1.5;			//deg Celcius of 24 hour average TSS
	const double TSS_threshold12_smallHSincrease=-0.5;	//deg Celcius of 12 hour average TSS in case of low rate of change of HS
	const double TSS_threshold12_largeHSincrease=1.0;	//deg Celcius of 12 hour average TSS in case of high rate of change of HS
	const double HS_threshold_smallincrease=0.01;		//low rate of change of HS (m/hour)
	const double HS_threshold_largeincrease=0.02;		//high rate of change of HS (m/hour)
	const double ThresholdSmallCanopy=1.;			//Set the threshold for the canopy height. Below this threshold, the canopy is considered to be small and snow will fall on top (like grass, or small bushes). Above this threshold, snow will fall through (like in forests). When the canopy is small, the measured snow height will be assigned to the canopy height in case of no snow pack on the ground.

	nOldN = Xdata.getNumberOfNodes();
	nOldE = Xdata.getNumberOfElements();
	cos_sl = cos(DEG_TO_RAD(Xdata.meta.getSlopeAngle()));
	if (surfaceCode == DIRICHLET_BC)
		t_surf = MIN(C_TO_K(-0.1), Mdata.tss);
	else
		t_surf = MIN(C_TO_K(-0.1), (Xdata.Ndata[nOldN-1].T + Mdata.ta)/2.);
	rho_hn = SnLaws::compNewSnowDensity(hn_density, hn_density_model,
	                                    Mdata, Xdata, t_surf, variant);
	if ((Sdata.cRho_hn < 0.) && (rho_hn != Constants::undefined))
		Sdata.cRho_hn = -rho_hn;

	if (!enforce_measured_snow_heights) {
		if (Mdata.ta < C_TO_K(thresh_rain)) {
			if ((cumu_hnw > 0.) && (rho_hn != Constants::undefined)) {
				if ((hn_density_model == "MEASURED")
				        || ((hn_density_model == "FIXED") && (rho_hn > SnLaws::max_hn_density))) {
					if ((meteo_step_length / sn_dt * Mdata.hnw) <= cumu_hnw) {
						Xdata.mH += (cumu_hnw/rho_hn);
						force_layer = 1;
					}
				} else if (cumu_hnw/rho_hn > height_new_elem*cos_sl) {
					Xdata.mH += (cumu_hnw/rho_hn);
					if (hn_density_model == "EVENT") {
						force_layer = 1;
					}
				}
			} else {
				return;
			}
		} else {
			// This is now very important to make sure that rain will not accumulate
			cumu_hnw -= Mdata.hnw;
		}
	}
	if (rho_hn == Constants::undefined)
		return;

	// To check thresholds for solid precipitation
	// NOTE No new snow during cloud free conditions
	snow_fall = (((Mdata.rh > thresh_rh) && (Mdata.ta < C_TO_K(thresh_rain)) && (Mdata.ta - Mdata.tss < 3.0))
                     || !enforce_measured_snow_heights || (Xdata.hn > 0.));
	// snowed_in is true if the ground is either already snowed in or snow will remain on it
	//  ... but check first for the availability of the following data. This is particularly important if enforce_measured_snow_heights == false.
	if (!enforce_measured_snow_heights || (Mdata.hs_change_rate == Constants::undefined)
		   || (Mdata.tss24 == Constants::undefined) || (Mdata.tss12 ==  Constants::undefined)) {
		snowed_in = true;
	} else {
		snowed_in = ( (Xdata.getNumberOfNodes() > Xdata.SoilNode+1)
		        || (Mdata.tss24!=IOUtils::nodata &&
		                Mdata.tss24 < C_TO_K(TSS_threshold24) && Mdata.hs_change_rate > HS_threshold_smallincrease)
		        || (Mdata.tss12!=IOUtils::nodata
		                && Mdata.tss12 < C_TO_K(TSS_threshold12_smallHSincrease)
		                && Mdata.hs_change_rate > HS_threshold_smallincrease)
		        || (Mdata.tss12!=IOUtils::nodata
		                && Mdata.tss12 < C_TO_K(TSS_threshold12_largeHSincrease)
		                && Mdata.hs_change_rate > HS_threshold_largeincrease)
		            );
	}

	//Now we check: we need snow fall AND ground which is snowed in or cold enough to maintain the snow pack. The latter condition is only relevant
	//if we are NOT in a slope! If we are in a slope, the slope should just get new snow when the flat field gets new snow:
	if (snow_fall && (snowed_in || (Xdata.meta.getSlopeAngle() > Constants::min_slope_angle))) {
		/* Now check if we have some canopy left. Then we reduce the canopy with the new snow height.
		 * This is to simulate the gradual sinking of the canopy under the weight of snow.
		 * We also adjust Xdata.mH, to have it reflect only measured snow, and not the canopy.
		 * (So if there is only a canopy, mH=0.)
		 * The third clause makes it a one way issue (only with snowfall, not with melt).
		 * The second clause limits this issue to small canopies only, to prevent problems
		 * with Alpine3D simulations in forests.
		 */
		if ((Xdata.Cdata.height > 0.) && (Xdata.Cdata.height < ThresholdSmallCanopy)
                && ((Xdata.cH + Xdata.Cdata.height) < Xdata.mH) && (Xdata.mH != Constants::undefined)) {
			/*First, reduce the Canopy height with the additional snow height. This makes the Canopy work
			 *like a spring. When increase in mh is 3 cm and the canopy height is 10 cm, the snow pack is
			 *assumed to be 6cm in thickness. When the total amount of snow measured (mh) equals the canopy
			 *height, the canopy is reduced to 0, and everything measured is assumed to be snow.
			 */
			Xdata.Cdata.height -= (Xdata.mH - (Xdata.cH + Xdata.Cdata.height));
			if (Xdata.Cdata.height < 0.) {
				Xdata.Cdata.height = 0.;
			} else {
				//Adjust measured snow height with canopy height ...
				Xdata.mH -= Xdata.Cdata.height;
				if (Xdata.mH < Xdata.Ground)
					Xdata.mH = Xdata.Ground;
			}
		}

		// In case of virtual slope use new snow depth and density from either flat field or luv slope
		if (((Xdata.cH + Xdata.Cdata.height) < (Xdata.mH - (height_new_elem * cos_sl)))
		        || (Xdata.hn > 0.) || force_layer) {
			if (Xdata.hn > 0. && (Xdata.meta.getSlopeAngle() > Constants::min_slope_angle)) {
				hn = Xdata.hn;
				rho_hn = Xdata.rho_hn;
			} else { // in case of flat field or PERP_TO_SLOPE
				hn = Xdata.mH - (Xdata.cH + Xdata.Cdata.height);
				// Store new snow depth and density
				if (!ALPINE3D) {
					Xdata.hn = hn;
					Xdata.rho_hn = rho_hn;
				}
			}
			if (hn > Snowpack::snowfall_warning)
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date,
				          "Large snowfall! hn=%.3lf cm (azi=%.0lf, slope=%.0lf)",
				            M_TO_CM(hn), Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
			nNewE = (int)(hn/(height_new_elem*cos_sl));
			if (nNewE < 1) {
				// Always add snow on virtual slope (as there is no storage variable available) and some other cases
				if (!ALPINE3D && ((Xdata.meta.getSlopeAngle() > Constants::min_slope_angle)
				                      || force_layer)) {
					nNewE = 1;
				} else {
					Xdata.hn = 0.;
					Xdata.rho_hn = Constants::undefined;
					return;
				}
			}
			cumu_hnw = 0.0;
			// Check whether surface hoar could be buried
			hoar = Xdata.Ndata[nOldN-1].hoar;
			if (nOldE > 0 && Xdata.Edata[nOldE-1].theta[SOIL] < Constants::eps2) {
				// W.E. of surface hoar must be larger than a threshold to be buried
				if (hoar > 1.5*MM_TO_M(hoar_min_size_buried)*hoar_density_surf) {
					nHoarE = 1;
				} else if (!(change_bc && meas_tss)
				               && (hoar > MM_TO_M(hoar_min_size_buried)*hoar_density_surf)) {
					nHoarE = 1;
				} else {
					nHoarE = 0;
				}
			} else { // Add surface hoar on ground to first new snow element(s)
				nHoarE = 0;
				hn += hoar/rho_hn;
				Xdata.Ndata[nOldN-1].hoar = 0.;
			}

			Xdata.cAlbedo = Snowpack::new_snow_albedo;

			nN = nOldN + nNewE + nHoarE;
			nE = nOldE + nNewE + nHoarE;
			Xdata.resize(nE);
			vector<NodeData>& NDS = Xdata.Ndata;
			vector<ElementData>& EMS = Xdata.Edata;

			// Create hoar layer
			if (nHoarE) {
				// Since mass of hoar was already added to element below, substract....
				// Make sure you don't try to extract more than is there
				hoar = MAX(0.,MIN(EMS[nOldE-1].M - 0.1,hoar));
				L0 = EMS[nOldE-1].L;
				dL = -hoar/(EMS[nOldE-1].Rho);
				EMS[nOldE-1].L0 = EMS[nOldE-1].L = L0 + dL;

				EMS[nOldE-1].E = EMS[nOldE-1].dE = EMS[nOldE-1].Ee = EMS[nOldE-1].Ev = EMS[nOldE-1].S = 0.0;
				Theta0 = EMS[nOldE-1].theta[ICE];
				EMS[nOldE-1].theta[ICE] *= L0/EMS[nOldE-1].L;
				EMS[nOldE-1].theta[ICE] += -hoar/(Constants::density_ice*EMS[nOldE-1].L);
				EMS[nOldE-1].theta[ICE] = MAX(EMS[nOldE-1].theta[ICE],0.);
				EMS[nOldE-1].theta[WATER] *= L0/EMS[nOldE-1].L;
				for (unsigned int ii = 0; ii < Xdata.number_of_solutes; ii++)
					EMS[nOldE-1].conc[ICE][ii] *= L0*Theta0/(EMS[nOldE-1].theta[ICE]*EMS[nOldE-1].L);
				EMS[nOldE-1].M -= hoar;
				EMS[nOldE-1].theta[AIR] = MAX(0., 1.0 - EMS[nOldE-1].theta[WATER]
				                                - EMS[nOldE-1].theta[ICE] - EMS[nOldE-1].theta[SOIL]);
				EMS[nOldE-1].Rho = (EMS[nOldE-1].theta[ICE] * Constants::density_ice)
				                      + (EMS[nOldE-1].theta[WATER] * Constants::density_water)
				                        + (EMS[nOldE-1].theta[SOIL]  * EMS[nOldE-1].soil[SOIL_RHO]);
				// Take care of old surface node
				NDS[nOldN-1].z += dL + NDS[nOldN-1].u;
				NDS[nOldN-1].u = 0.0;
				NDS[nOldN-1].hoar = 0.0;
				// Now fill nodal data for upper hoar node
				NDS[nOldN].T = t_surf;              // The temperature of the new node
				// The new nodal position;
				NDS[nOldN].z = NDS[nOldN-1].z + NDS[nOldN-1].u + hoar/hoar_density_buried;
				NDS[nOldN].u = 0.0;                 // Initial displacement is 0
				NDS[nOldN].hoar = hoar / hoar_density_buried;         // Surface hoar initial size
				NDS[nOldN].udot = 0.0;               // Settlement rate is also 0
				NDS[nOldN].f = 0.0;                 // Unbalanced forces is 0
				NDS[nOldN].S_n = INIT_STABILITY;
				NDS[nOldN].S_s = INIT_STABILITY;
			} else { // Make sure top node surface hoar mass is removed
				NDS[nOldN-1].hoar = 0.0;
			}
			// Fill the nodal data
			if (!useSoilLayers && (nOldN-1 == Xdata.SoilNode)) // New snow on bare ground w/o soil
				NDS[nOldN-1].T = (t_surf + Mdata.ta)/2.;
			Ln = (hn / nNewE);
			z0 = NDS[nOldN-1+nHoarE].z + NDS[nOldN-1+nHoarE].u + Ln;
			for (n = nOldN+nHoarE; n < nN; n++) {
				NDS[n].T = t_surf;                  // The temperature of the new node
				NDS[n].z = z0;                      // The new nodal position;
				NDS[n].u = 0.0;                     // Initial displacement is 0
				NDS[n].hoar = 0.0;                  // The new snow surface hoar is set to zero
				NDS[n].udot = 0.0;                  // Settlement rate is also 0
				NDS[n].f = 0.0;                     // Unbalanced forces is 0
				NDS[n].S_n = INIT_STABILITY;
				NDS[n].S_s = INIT_STABILITY;
				z0 += Ln;
			}

			// Fill the element data
			for (e = nOldE; e < nE; e++) {
				// Birthdate
				EMS[e].depositionDate = Mdata.date;
				// Temperature
				EMS[e].Te = t_surf;
				// Lengths
				EMS[e].L0 = EMS[e].L = (NDS[e+1].z + NDS[e+1].u) - (NDS[e].z + NDS[e].u);
				// Density
				EMS[e].Rho = rho_hn;
				if (nHoarE && e == nOldE)
					EMS[e].Rho = hoar_density_buried;
				// Mass
				EMS[e].M = EMS[e].L0*EMS[e].Rho;
				// Volumetric components
				EMS[e].theta[SOIL]  = 0.0;
				EMS[e].theta[ICE]   = EMS[e].Rho/Constants::density_ice;
				EMS[e].theta[WATER] = 0.0;
				EMS[e].theta[AIR]   = 1. - EMS[e].theta[ICE];
				for (unsigned int ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[e].conc[ICE][ii]   = Mdata.conc[ii]*Constants::density_ice/Constants::density_water;
					EMS[e].conc[WATER][ii] = Mdata.conc[ii];
					EMS[e].conc[AIR][ii]   = 0.0;
					EMS[e].conc[SOIL][ii]  = 0.0;
				}
				// Coordination number based on Bob's empirical function
				EMS[e].N3 = Metamorphism::getCoordinationNumberN3(EMS[e].Rho);
				// Constitutive Parameters
				EMS[e].heatCapacity();
				EMS[e].c[SEEPAGE] = EMS[e].c[SETTLEMENT]= 0.0;
				EMS[e].k[TEMPERATURE] = EMS[e].k[SEEPAGE] = EMS[e].k[SETTLEMENT]= 0.0;
				EMS[e].soil[SOIL_RHO] = EMS[e].soil[SOIL_K] = EMS[e].soil[SOIL_C] = 0.0;
				// Set the initial short wave radiation to zero
				EMS[e].sw_abs = 0.0;
				// Phase change variables:
				EMS[e].dth_w=0.0; // change of water content
				EMS[e].Qmf=0.0;   // change of energy due to phase changes
				Xdata.ColdContent += EMS[e].coldContent();
				// Total element strain (GREEN'S strains -- TOTAL LAGRANGIAN FORMULATION.
				EMS[e].dE = EMS[e].E = EMS[e].Ee = EMS[e].Ev = 0.0;
				// Total Strain Rate	(Simply, E/sn_dt)
				EMS[e].EDot = EMS[e].EvDot=0.0;
				// Total Element Stress
				EMS[e].S=0.0;
				EMS[e].CDot = 0.; // loadRate
				// NEW SNOW MICRO-STRUCTURE	 fitted from different data sources
				{
					double logit, value;
					double alpha, beta, gamma, delta, eta, phi;
					const double TA = K_TO_C(Mdata.ta);
					const double RH = Mdata.rh*100.;
					// Distinguish between Graupel and New Snow
					alpha = 49.6; beta = 0.857; gamma = -0.547;
					logit = alpha + beta*Mdata.vw + gamma*RH;
					value = exp(logit)/(1.+exp(logit));
					if (value > 1.0) { // Graupel
						EMS[e].mk = 4;
						EMS[e].dd = 0.;
						EMS[e].sp = 1.;
						EMS[e].rg = 0.6;
						EMS[e].rb = 0.2;
						// Because density and volumetric contents are already defined, redo it here
						EMS[e].Rho = 110.;
						EMS[e].theta[ICE] = EMS[e].Rho/Constants::density_ice;  // ice content
						EMS[e].theta[AIR] = 1. - EMS[e].theta[ICE];  // void content
					} else { // no Graupel
						EMS[e].mk = Snowpack::new_snow_marker;
						if (SnLaws::jordy_new_snow && (Mdata.vw > 2.9)
						        && ((hn_density_model == "LEHNING_NEW") || (hn_density_model == "LEHNING_OLD"))) {
							alpha = 1.87; beta = -0.04;
							EMS[e].dd = MAX(0.5,MIN(1.0,(alpha + beta*Mdata.vw)*(alpha + beta*Mdata.vw)));
							EMS[e].sp = new_snow_sp;
							alpha = 0.9; beta = 0.015; gamma = -0.0062; delta = -0.117; eta=0.0011; phi=-0.0034;
							EMS[e].rg = MIN(new_snow_grain_rad, MAX(0.3*new_snow_grain_rad,
							                  alpha + beta*TA + gamma*RH + delta*Mdata.vw
							                    + eta*RH*Mdata.vw + phi*TA*Mdata.vw));
							EMS[e].rb = 0.4*EMS[e].rg;
						} else {
							EMS[e].dd = new_snow_dd;
							EMS[e].sp = new_snow_sp;
							// Adapt dd and sp for blowing snow
							if ((Mdata.vw > 5.) && ((variant == "ANTARCTICA")
									|| (!SnLaws::jordy_new_snow && ((hn_density_model == "BELLAIRE")
										|| (hn_density_model == "LEHNING_NEW"))))) {
								EMS[e].dd = new_snow_dd_wind;
								EMS[e].sp = new_snow_sp_wind;
							} else if (vw_dendricity && ((hn_density_model == "BELLAIRE")
								           || (hn_density_model == "ZWART"))) {
								double vw;
								vw = MAX(0.05, Mdata.vw);
								// dd = f(vw)
								EMS[e].dd = (1 - pow(vw/10., 1.57));
								EMS[e].dd = MAX(0.2, EMS[e].dd);
							}
							if (Snowpack::hydrometeor) { // empirical
								alpha=1.4; beta=-0.08; gamma=-0.15; delta=-0.02;
								EMS[e].rg = 0.5*(alpha + beta*TA + gamma*Mdata.vw + delta*TA*Mdata.vw);
								EMS[e].rb = 0.25*EMS[e].rg;
							} else {
								EMS[e].rg = new_snow_grain_rad;
								EMS[e].rb = new_snow_bond_rad;
								if (((Mdata.vw_avg >= SnLaws::event_wind_lowlim) && (Mdata.rh_avg >= rh_lowlim))) {
									EMS[e].rb = MIN(bond_factor_rh*EMS[e].rb,
									                    Metamorphism::max_grain_bond_ratio*EMS[e].rg);
								}
							}
						}
					} // end no Graupel
					if (nHoarE &&	 e == nOldE) {
						EMS[e].mk = 3;
						EMS[e].dd = 0.;
						EMS[e].sp = 0.;
						EMS[e].rg = MAX(new_snow_grain_rad, 0.5*M_TO_MM(EMS[e].L0));
						// Note: L0 > hoar_min_size_buried/hoar_density_buried
						EMS[e].rb = EMS[e].rg/3.;
					}
					EMS[e].metamo = 0.;
				} // Treat all the initial snow types

				// Snow classification
				EMS[e].snowType();

				// Initialise the Stability Index for ml_st_CheckStability routine.)
				EMS[e].S_dr = INIT_STABILITY;
				EMS[e].hard = 0.;
			}   // End elements

			// Finally, update snowpack height
			Xdata.cH = Xdata.mH = NDS[nN-1].z + NDS[nN-1].u;
			Xdata.ErosionLevel = nE-1;
		} //  End of NEW_SNOWFALL Condition 2
	} else {
		// If there is no snowfall and no snowpack yet, we can assign the measured snow height to the canopy,
		// but only for small canopies, to prevent problems with Alpine3D simulations in forests.
		if (Xdata.Cdata.height < ThresholdSmallCanopy) {
			if ((Xdata.getNumberOfNodes() == Xdata.SoilNode+1) && (Xdata.mH != Constants::undefined)) {
				Xdata.Cdata.height = Xdata.mH - Xdata.Ground;
			}
		}
	}
}

/**
 * @brief The near future (s. below) has arrived on Wednesday Feb. 6, when it was finally snowing
 * in Davos and Sergey, Michael and Perry were working furiously on SNOWPACK again. Michael
 * prepared the coupling of the model to the energy balance model of Olivia and his own snow
 * drift model. He found it therefore necessary to move the time loop into the Main.c frame.
 * He also cleaned the data structure handling, i.e. made sure that all variables are handed
 * to the subroutines and that there is no longer a global referencing. In the meantime, Perry
 * and Sergey prepared the introduction of a vapor transport equation. While Sergey was sad
 * that his family had to go back, Perry was exchanging a lot of wet kisses with his new women
 * probably causing some of his confusion about vapor and heat transport.
 * Driving routine for the 1d-snowpack model.  The main program will probably be replaced
 * in future -- or, at least modified.  For now, it is responsible for several tasks:
 * -# In the next step the program enters the time-integration, or, time-stepping loop.
 *    Moreover, we find the solution at TimeN = Time0 + sn_dt.  For now, we will assume
 *    that sn_dt = 15min and Time0 = 0.0.  These restrictions can be later relaxed.
 *    The program must then "check" three modules before finding the temperature distribution
 *    and creep deformations, these are:
 *    -#   Determine the METEO parameters at time TimeN. (The METEO parameters include
 *             air temperature, relative humidity, short wave radiation, incoming longwave
 *             radiation, etc. as well as the ground temperature)
 *             >>>>>> METEO MODULE <<<<<<<
 *             Change: ML, 06.02.02: Mdata needs to be handed to the function
 *    -#   Determine if a SNOWFALL event as occcured. (This implies that the dynamically
 *             allocated data structures must be reallocated and the new material defined.
 *             The mesh connectivies (sn_NewMesh=TRUE) must be rebuilt for the numerical
 *             solution and the energy exchanges at the top surface initialized.)
 *             >>>>>> SNOWFALL MODULE <<<<<<<
 *    -#   Determine if SNOWDRIFT is occuring. Note that SNOWFALL distributes mass
 *             equally at ALL EXPOSITIONS.  SNOWDRIFT can occur in conjunction WITH SNOWFALL
 *             or WITHOUT SNOWFALL.  Mass is redistributed between the expositions. (Again,
 *             the internal data structures must be reinitialized.)
 *             >>>>>> SNOWDRIFT MODULE <<<<<<<
 * -# The phase change module can generate water; this water is moved downwards by the
 *    WATER TRANSPORT module.  At present, a simple, mass conserving scheme is employed to move
 *    the water; when the water content is greater than the residual water content then the excess
 *    water is moved to the next element.  If the element reaches the bottom of the snowpack
 *    water is removed -- this is called MELTWATER RUNOFF.  Also, elements which do not contain
 *    any ice are removed from the finite element mesh.
 * -# In the next step the TEMPERATURE distribution within the snowpack at each exposition
 *    is found.  If SNOWFALL and SNOWDRIFT has occured then the finite element matrices must
 *    be rebuilt.  >>>>> TEMPERATURE MODULE  <<<<<<<<
 * -# If the computed temperatures are above zero degrees then PHASECHANGES are occuring.
 *    This means that the volumetric contents of the elements must be updated. >>>>>> PHASE CHANGE MODULE <<<<<
 * -# In the next step the CREEPING deformations and 1d state of stress within the snowpack
 * at each exposition are found.  If SNOWFALL and SNOWDRIFT has occured then the finite
 * element matrices (connectivities) must be rebuilt.  >>>>CREEP MODULE<<<<<
 * @param Mdata
 * @param Xdata
 * @param cumu_hnw Variable to store amount of precipitation (kg m-2)
 * @param Bdata
 * @param Sdata
 */
void Snowpack::runSnowpackModel(CurrentMeteo& Mdata, SnowStation& Xdata, double& cumu_hnw,
                                BoundCond& Bdata, SurfaceFluxes& Sdata)
{
	WaterTransport watertransport(cfg);
	Metamorphism metamorphism(cfg);
	SnowDrift snowdrift(cfg);
	PhaseChange phasechange(cfg);

	try {
		// Adjust boundary condition
		if (change_bc && meas_tss) {
			if (Mdata.tss < C_TO_K(thresh_change_bc))
				surfaceCode = DIRICHLET_BC;
			else
				surfaceCode = NEUMANN_BC;
		}

		// If it is SNOWING, find out how much, prepare for new FEM data
		compSnowFall(Mdata, Xdata, cumu_hnw, Sdata);

		// Check to see if snow is DRIFTING, compute a simple snowdrift index and erode layers if
		// neccessary. Note that also the very important friction velocity is computed in this
		// routine and later used to compute the Meteo Heat Fluxes
		snowdrift.compSnowDrift(Mdata, Xdata, Sdata, cumu_hnw);

		// Reinitialize and compute the initial meteo heat fluxes
		memset((&Bdata), 0, sizeof(BoundCond));
		updateMeteoHeatFluxes(Mdata, Xdata, Bdata);

		// Find the temperature in the snowpack
		compSnowTemperatures(Xdata, Mdata, Bdata);

		//Good HACK (according to Charles, qui persiste et signe;-)... like a good hunter and a bad one...
		// If you switched from DIRICHLET to NEUMANN boundary conditions, correct
		//   for a possibly erroneous surface energy balance. The latter can be due e.g. to a lack
		//   of data on nebulosity leading to a clear sky assumption for incoming long wave.
		if ((change_bc && meas_tss) && (surfaceCode == NEUMANN_BC)
				&& (Xdata.Ndata[Xdata.getNumberOfNodes()-1].T < C_TO_K(thresh_change_bc))) {
			surfaceCode = DIRICHLET_BC;
			Xdata.Ndata[Xdata.getNumberOfNodes()-1].T = C_TO_K(thresh_change_bc/2.);
			compSnowTemperatures(Xdata, Mdata, Bdata);
		}

		// See if any SUBSURFACE phase changes are occuring
		phasechange.compPhaseChange(Sdata, Xdata);

		// Compute change of internal energy during last time step (J m-2)
		Xdata.compSnowpackInternalEnergyChange(sn_dt);
		// Compute the final meteo heat fluxes
		Sdata.ql += Bdata.ql; // Bad;-) HACK, needed because latent heat ql is not (yet)
		                      // linearized w/ respect to Tss and thus reamins unchanged
		                      // throughout the temperature iterations!!!
		updateMeteoHeatFluxes(Mdata, Xdata, Bdata);

		// The water transport routines must be placed here, otherwise the temperature
		// and creep solution routines will not pick up the new mesh boolean.
		watertransport.compTransportMass(Mdata, Bdata.ql, Xdata, Sdata);


		// Find the settlement of the snowpack.
		// HACK This routine was formerly placed here because the settlement solution MUST ALWAYS follow
		// computeSnowTemperatures where the vectors U, dU and dUU are allocated.
		compSnowCreep(Mdata, Xdata);

	} catch(const exception&) {
		prn_msg(__FILE__, __LINE__, "err", Mdata.date, "Snowpack computation not completed");
		throw;
	}

	metamorphism.runMetamorphismModel(Mdata, Xdata);

	if (join_elements)
		Xdata.joinElements(SnowStation::number_top_elements);
}
