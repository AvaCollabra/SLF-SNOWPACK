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
/**
 * @file Snowpack.c
 * @version 10.02
 * @bug     -
 * @brief This module contains the driving routines for the 1d snowpack model
 */

#include <snowpack/Snowpack.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/

const bool Snowpack::jordy_new_snow = false;

//Uses an empirically determined size of deposited hydrometeors as new snow grain size (mm)
const bool Snowpack::hydrometeor = false;

//Warning is issued if depth of snowfall is larger than this amount (m)
const double Snowpack::snowfall_warning = 0.5;

const int Snowpack::new_snow_marker = 0;
const double Snowpack::new_snow_albedo = 0.9;

double Snowpack::fixed_hn_density        = 100.0;
double Snowpack::max_hn_density          = 250.0;
double Snowpack::event_wind_lowlim       = 0.0;
double Snowpack::event_wind_highlim      = 0.0;
Snowpack::EventType	Snowpack::event_type = Snowpack::EVENT_DEFAULT;
Snowpack::NewSnowDensityModel Snowpack::hn_density_model = Snowpack::ZWART;

void Snowpack::initStaticData(const std::string& variant)
{
	/* Wind effect on new snow density, dendricity and grain size
	 * EVENT_... : thresholds for event driven new snow density and bond size */
	if (variant == "ANTARCTICA"){
		Snowpack::fixed_hn_density   = 300.0;
		Snowpack::max_hn_density     = 340.0;
		Snowpack::event_type         = Snowpack::EVENT_WIND;
		Snowpack::event_wind_lowlim  = 4.0;
		Snowpack::event_wind_highlim = 7.0;
		Snowpack::hn_density_model   = Snowpack::EVENT;
	}
}

/************************************************************
 * non-static section                                       *
 ************************************************************/

Snowpack::Snowpack(const mio::Config& i_cfg) : cfg(i_cfg)
{
	string tmp_variant = cfg.get("VARIANT", "Parameters");
	variant = tmp_variant;
	initStaticData(variant);

	/** Defines the management of the bottom boundary conditions with soil layers
	 * - 0 ==> Dirichlet, i.e fixed Temperature
	 * - 1 ==> Neumann, fixed geothermal heat flux GEO_HEAT */
	SOIL_FLUX = cfg.get("SOIL_FLUX", "Parameters");

	GEO_HEAT = cfg.get("GEO_HEAT", "Parameters"); //Constant geothermal heat flux at (great) depth (W m-2)

	/* Defines the management of the surface boundary conditions
	 * - 0: Neumann boundary conditions throughout
	 * - 1: Dirichlet if Tss < THRESH_CHANGE_BC, Neumann else */
	change_bc = cfg.get("CHANGE_BC", "Parameters");
	THRESH_CHANGE_BC = cfg.get("THRESH_CHANGE_BC", "Parameters");

	//Should be 0 for data-sets which do not provide measured surface temperatures
	meas_tss = cfg.get("MEAS_TSS", "Parameters");;

	/**
	 * @brief Defines how the height of snow is going to be handled
	 * - 0: Depth of snowfall is determined from the water equivalent of snowfall (HNW)
	 * - 1: The measured height of snow is used to determine whether new snow has been deposited.
	 *      This setting MUST be chosen in operational mode. \n
	 *      This procedure has the disadvantage that if the snowpack settles too strongly
	 *      extra mass is added to the snowpack. \n
	 * New snow density is needed in both cases, either parameterized, measured, or fixed.
	 */
	ENFORCE_MEASURED_SNOW_HEIGHTS = cfg.get("ENFORCE_MEASURED_SNOW_HEIGHTS", "Parameters");

	/**
	 * @brief Defines whether the canopy model is used \n
	 * NOTE: OUT_CANOPY must also be set to dump canopy parameters to file; see Constants_local.h
	 */
	CANOPY = cfg.get("CANOPY", "Parameters");

	/**
	 * @brief Define the heights of the meteo measurements above ground (m) \n
	 * Required for surface energy exchange calculation and for drifting and blowing snow.
	 */
	HEIGHT_OF_METEO_VALUES = cfg.get("HEIGHT_OF_METEO_VALUES", "Parameters");

	/// @brief Defines whether soil layers are used
	//SNP_SOIL = cfg.get("SNP_SOIL", "Parameters");

	/// @brief Time resolution of input data (min; oper: 30.)
	METEO_STEP_LENGTH = cfg.get("METEO_STEP_LENGTH", "Parameters");

	/**
	 * @brief Defines whether the measured shortwave radiation is incoming
	 * - (SW_REF 0 || 10) only incoming SW used
	 * - (SW_REF 1 || 11) only reflected SW used
	 * - (SW_REF = 2 || 12) both incoming and reflected SW used \n
	 * NOTE: If SW_REF == 2 || SW_REF >= 10, the input file must hold both fluxes!
	 */
	SW_REF  = cfg.get("SW_REF", "Parameters");

	/**
	 * @brief Height of new snow element (m) [NOT read from CONSTANTS_User.INI] \n
	 * Controls the addition of new snow layers. Set in qr_ReadParameters() \n
	 * The value depends on ENFORCE_MEASURED_SNOW_HEIGHTS:
	 * - 0: 2.0*MINIMUM_L_ELEMENT (value depends on VARIANT)
	 * - 1: 0.02
	 */
	hns_ne_height = cfg.get("HNS_NE_HEIGHT", "Parameters");

	// Defines whether soil layers are used
	SNP_SOIL = cfg.get("SNP_SOIL", "Parameters");


	string tmp_viscosity_model = cfg.get("VISCOSITY_MODEL", "Parameters"); 
	viscosity_model = tmp_viscosity_model;

	//Rain only for air temperatures warmer than threshold (degC)
	thresh_rain = cfg.get("THRESH_RAIN", "Parameters");

	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Parameters");
	sn_dt = M_TO_S(calculation_step_length);

	//Defines whether a multiband model is used for short wave radiation extinction
	multistream = cfg.get("MULTISTREAM", "Parameters");

	//Defines whether joining elements will be considered at all
	join_elements = cfg.get("JOIN_ELEMENTS", "Parameters");

	//Which boundary condition to use
	string boundarycondition = cfg.get("SURFACECODE", "Parameters");
	if (boundarycondition == "NEUMANN_BC"){
		surfaceCode = Snowpack::NEUMANN_BC;
	} else if (boundarycondition == "DIRICHLET_BC"){
		surfaceCode = Snowpack::DIRICHLET_BC;
	} else {
		throw InvalidArgumentException("The SURFACECODE " + boundarycondition + " is unknown", AT);
	}

	//Warning is issued if snow tempeartures are out of bonds, that is, crazy
	t_crazy_min = cfg.get("T_CRAZY_MIN", "Parameters");
	t_crazy_max = cfg.get("t_crazy_max", "Parameters");
	
	/* Initial new snow parameters, see sn_SnowFall()
	 * - that rg and rb are equal to 0.5*gsz and 0.5*bsz, respectively. Both given in millimetres
	 * - If VW_DENDRICITY is set, new snow dendricity is f(vw)
	 * - BOND_FACTOR_RH new snow bonds get stronger for average winds >= EVENT_WIND_LOWLIM and
	 *   mean relative humidity >= RH_LOWLIM */
	if (variant == "ANTARCTICA"){
		new_snow_dd = 0.5;
		new_snow_sp = 0.75;
		new_snow_dd_wind = 0.15;
		new_snow_sp_wind = 1.0;
		vw_dendricity    = false;
		rh_lowlim = 0.7;
		bond_factor_rh = 3.0;
	} else {
		new_snow_dd = 1.0;
		new_snow_sp = 0.5;
		new_snow_dd_wind = 0.5;
		new_snow_sp_wind = 0.75;
		vw_dendricity    = true;
		rh_lowlim = 1.0;
		bond_factor_rh = 1.0;
	}

	new_snow_grain_rad = cfg.get("NEW_SNOW_GRAIN_RAD", "Parameters");
	new_snow_bond_rad = 0.25 * new_snow_grain_rad;
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
 * @return ERROR if any, NO_ERROR if not
 */
int Snowpack::sn_SnowForces(SN_ELEM_DATA *Edata,  double dt, double cos_sl, double Zn[ N_OF_INCIDENCES ], double Un[ N_OF_INCIDENCES ], double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ], double Fc[ N_OF_INCIDENCES ], double Fi[ N_OF_INCIDENCES ], double Fe[ N_OF_INCIDENCES ])
{
	double L, L0; // Length, initial length
	double dVol;  // Change in volume
	double D;     // Modulus of elasticity
	double E;     // Strain
	double S;     // Stress
	double Sc;    // Psuedo Creep Stress
	double ddE;   // Change in axial strian, volumetric strain

	// First calculate the new length and from the new length the total strain
	dVol = Edata->L;
	L0 = Edata->L0; L = Edata->L = Zn[1] + Un[1] - Zn[0] - Un[0];
	dVol /= L;
	if ( L <= 0. ) {
		prn_msg(__FILE__, __LINE__, "err", -1., "Element length < 0.0!\n L0=%lf L=%lf Zn[0]=%lf Zn[1]=%lf Un[0]=%lf Un[1]=%lf,", L0, L, Zn[0], Zn[1], Un[0], Un[1]);
		return(ERROR);
	}
	// Calculate the Natural Strain // Green Lagrange Strain
	E = Edata->E = log(L / L0);
	if ( !(E >= -100. && E <= 1.e-5) ) {
		prn_msg(__FILE__, __LINE__, "err", -1., "In strain E (Memory?)");
		return(ERROR);
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
	if ( ddE <= -1.5 ) {
		ddE = -1.5;
		prn_msg(__FILE__, __LINE__, "wrn", -1., "Time step is too Large!!!");
	}
	// Calculate the volumetric components
	Edata->theta[ICE] = MIN (0.999999, Edata->theta[ICE] * dVol);
	Edata->theta[WATER] = MIN (0.999999, Edata->theta[WATER] * dVol);
	Edata->theta[AIR] = 1.0 - Edata->theta[ICE] - Edata->theta[WATER] - Edata->theta[SOIL];
	if ( !(Edata->theta[AIR] <= 1.0 && Edata->theta[AIR] >= -0.05) ) {
		prn_msg(__FILE__, __LINE__, "msg+", -1., "ERROR AIR: %e (ddE=%e)", Edata->theta[AIR], ddE);
		prn_msg(__FILE__, __LINE__, "msg", -1., "ELEMENT SIZE: L0=%e L=%e", Edata->L0, Edata->L);
		prn_msg(__FILE__, __LINE__, "msg", -1., "DENSITY: rho=%e", Edata->Rho);
		prn_msg(__FILE__, __LINE__, "msg", -1., "ThetaICE: ti=%e", Edata->theta[ICE]);
		prn_msg(__FILE__, __LINE__, "msg", -1., "ThetaWATER: ti=%e", Edata->theta[WATER]);
		return(ERROR);
	}
	// Calculate the self weight of the element (with the present mass)
	Fe[0] = Fe[1] = -(Edata->M * Constants::g * cos_sl) / 2.;
	// Calculate the elastic strain and stress
	Edata->Ee = Edata->E - Edata->Ev;
	D = lwsn_SnowElasticity(Edata->Rho);
	S = Edata->S = D * Edata->Ee;
	// Calculate the internal forces
	Fi[0] = -S;
	Fi[1] = S;
	// Calculate the pseudo increment in creep stress
	Sc = D * Edata->EvDot * dt;
	// Calculate the creep forces
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
	// Calculate the stiffness matrix
	Se[0][0] = D/L ;  Se[0][1] = -D/L ;
	Se[1][0] =-D/L ;  Se[1][1] = D/L ;
	return(NO_ERROR);
}

/**
 * @brief Snow creep
 * -# The Thing ain't settling any more in case of ice, soil or water only
 * -# Enhanced densification for wind slabs of Metamorphism::wind_slab_depth (m); see also mm_Metamorphism()
 *    dry snow AND strong wind AND near the surface => enhance densification \n
 * -# Normal densification
 * -# Empirism for surface hoar. Two different rates are used for either large or small SH.
 *    Implemented by Sascha Bellaire on 28.11.2006.
 * @param Xdata
 * @param Mdata
 */
void Snowpack::calcSnowCreep(const SN_MET_DATA& Mdata, SN_STATION_DATA& Xdata)
{
	int nE, nN;             // Number of elements and nodes
	int e;                  // Element counter
	double L0, dL, cH_old;  // Element length and change of length
	double Sig0, SigC, eta; // Initial and Cauchy stress, resp.; snow viscosity
	SN_ELEM_DATA *EMS;      // Dereferenced element pointer

	vector<SN_NODE_DATA>& NDS = Xdata.Ndata; nN = Xdata.getNumberOfNodes();
	EMS = &Xdata.Edata[0]; nE = Xdata.getNumberOfElements();
	for (e = nE-1, SigC = 0.0; e >= 0; e--) {
		SigC += -EMS[e].M * Constants::g * cos(Xdata.SlopeAngle);
		EMS[e].C = SigC;
	}
	if ( nN == Xdata.SoilNode + 1 ) {
		return;
	}

	for (e = Xdata.SoilNode; e < nE; e++) {
		if (EMS[e].Rho > 910. ||  EMS[e].theta[SOIL] > 0. || EMS[e].theta[ICE] < EPS) {
			EMS[e].k[SETTLEMENT] = eta = 1.0e99;
		} else {
			EMS[e].k[SETTLEMENT] = eta = lwsn_SnowViscosity(viscosity_model, EMS[e], Mdata.date);
			if ( !(eta > 0.01 * SMALLEST_VISCOSITY && eta <= 1.e11 * SMALLEST_VISCOSITY) && (EMS[e].theta[ICE] > 2. * MIN_ICE_CONTENT) && (EMS[e].theta[ICE] < 0.6) ) {
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date.getJulianDate(), "Viscosity=%e out of range! e=%d nE=%d rg=%lf rb=%lf dd=%lf sp=%lf theta_i=%lf theta_w=%lf", eta, e, nE, EMS[e].rg, EMS[e].rb, EMS[e].dd, EMS[e].sp, EMS[e].theta[ICE], EMS[e].theta[WATER]);
			}
			if ( eta < SMALLEST_VISCOSITY ) {
				if ( 0 ) {
					prn_msg(__FILE__, __LINE__, "wrn", Mdata.date.getJulianDate(), "Viscosity=%e reset to SMALLEST_VISCOSITY! e=%d nE=%d", eta, e, nE);
				}
				EMS[e].k[SETTLEMENT] = eta = SMALLEST_VISCOSITY;
			}
		}
		Sig0 = lwsn_InitialStress(viscosity_model, EMS[e]);
		L0 = EMS[e].L;

		if ( !(EMS[e].mk%100 == 3) ) {
			double wind_slab=1.;
			if ( (EMS[e].theta[WATER] < 0.01) && (Mdata.vw > Metamorphism::wind_slab_vw) && ((NDS[nE].z - NDS[e].z < Metamorphism::wind_slab_depth) || e==nE-1) ) {
				wind_slab += Metamorphism::wind_slab_enhance * (Mdata.vw - Metamorphism::wind_slab_vw);
			}
			EMS[e].EvDot = wind_slab * (EMS[e].C - Sig0 ) / eta;
			dL = L0 * sn_dt * EMS[e].EvDot;
			if ( (L0 + dL) < MINIMUM_L_ELEMENT ) {
				dL = MIN(0., MINIMUM_L_ELEMENT - L0);
			}
		} else { //SH
			if ( NDS[e+1].hoar > 0.006 ) { // TODO Large initial size, i.e., deposited hoar mass/DENSITY_HOAR_BURIED ??
				if ( (Mdata.date.getJulianDate() - EMS[e].date.getJulianDate()) < 21. ) {
					dL = MM_TO_M(-0.391 * S_TO_D(sn_dt));
				} else {
					dL = MM_TO_M(-0.0807 * S_TO_D(sn_dt));
				}
			} else {
				dL = MM_TO_M(-0.1107 * S_TO_D(sn_dt));
			}
			EMS[e].EvDot = dL / (L0 * sn_dt);
			if ( (L0 + dL) < 0. ) {
				dL = 0.;
			}
		}

		if (variant == "CALIBRATION"){
			NDS[e].udot = EMS[e].C/eta;
			EMS[e].EDot = -Sig0/EMS[e].C;
			NDS[e].f    = -Sig0/eta;
		}

		EMS[e].theta[WATER] *= L0 / (L0 + dL);
		EMS[e].theta[ICE]   *= L0 / (L0 + dL);
		while (EMS[e].theta[ICE] + EMS[e].theta[WATER] + EMS[e].theta[SOIL] > 0.99) {
			EMS[e].theta[ICE] *= 0.99;
			EMS[e].theta[WATER] *= 0.99;
			EMS[e].M = L0 * ((EMS[e].theta[ICE] * Constants::DENSITY_ICE) 
						  + (EMS[e].theta[WATER] * Constants::DENSITY_WATER));
		}
		EMS[e].theta[AIR] = 1.0 - EMS[e].theta[WATER] - EMS[e].theta[ICE] - EMS[e].theta[SOIL];
		EMS[e].Rho = (EMS[e].theta[ICE] * Constants::DENSITY_ICE) + (EMS[e].theta[WATER] * Constants::DENSITY_WATER) 
			+ (EMS[e].theta[SOIL] * EMS[e].soil[SOIL_RHO]);
		EMS[e].L0 = EMS[e].L = (L0 + dL);
		NDS[e+1].z = NDS[e].z + EMS[e].L;
		if (! (EMS[e].Rho > 0. && EMS[e].Rho <= MAX_RHO)) {
			prn_msg(__FILE__, __LINE__, "err", -1., "Volume contents: e=%d nE=%d rho=%lf ice=%lf wat=%lf air=%le", 
				   e, nE, EMS[e].Rho, EMS[e].theta[ICE], EMS[e].theta[WATER], EMS[e].theta[AIR]);

			throw IOException("Runtime Error in calcSnowCreep()", AT);
		}
	} // end for settling-loop
	cH_old = Xdata.cH;
	Xdata.cH = NDS[nN-1].z + NDS[nN-1].u;
	Xdata.mH -= (cH_old - Xdata.cH);
}

/**
 * @brief Calculates the element stiffness matrix and right hand side vector for a fully implicit time integration scheme \n
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
 * @return ERROR if any, NO_ERROR if not
 */
int Snowpack::sn_ElementKtMatrix(SN_ELEM_DATA *Edata, double dt, double dvdz, double T0[ N_OF_INCIDENCES ], double Tn[ N_OF_INCIDENCES ], double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ], double Fe[ N_OF_INCIDENCES ], char *SubSurfaceMelt, char *SubSurfaceFrze, double VaporEnhance)
{
	double k, c;    // conductivity and heat capacity
	double Keff;    // the effective thermal conductivity

	// Gather the data which is required to calculate the element stiffness and force vector
	Fe[0] = Fe[1] = 0.0;
	if ( Edata->L < 0.0 ) {
		prn_msg(__FILE__, __LINE__, "err", -1., "Negative length L=%e", Edata->L);
		return(ERROR);
	}

	// Phase Change Check
	if ( (Tn[0] > MELTING_TK || Tn[1] > MELTING_TK) && Edata->theta[ICE] > 0. ) {
		*SubSurfaceMelt = TRUE;
	}
	if ( (Tn[0] < FREEZING_TK || Tn[1] < FREEZING_TK) && Edata->theta[WATER] > 0. ) {
		*SubSurfaceFrze = TRUE;
	}
	// Calculate the element mean temperature and temperature gradient
	if ( Edata->theta[WATER] > 0. && Edata->theta[SOIL] == 0. ) {
		 T0[0] = T0[1] = MELTING_TK;
	}

	// Find the conductivity of the element
	if ( Edata->theta[SOIL] > 0.0 ) {
		Keff = lwsn_SoilThermalConductivity(Edata, dvdz);
	} else if ( Edata->theta[ICE] > 0.55 || Edata->theta[ICE] < MIN_ICE_CONTENT ) {
		Keff = Edata->theta[AIR] * CONDUCTIVITY_AIR + Edata->theta[ICE] * CONDUCTIVITY_ICE + Edata->theta[WATER] * CONDUCTIVITY_WATER + Edata->theta[SOIL] * Edata->soil[SOIL_K];
	} else {
		Keff = lwsn_SnowThermalConductivity(*Edata, dvdz);
	}
	// mimics effect of vapour transport if liquid water present in snowpack
	Keff *= VaporEnhance;
	Edata->k[TEMPERATURE] = Keff;
	k = Keff/Edata->L;   // Divide by the length to save from doing it during the matrix operations

	// Compute heat capacity of element
	Edata->c[TEMPERATURE] = lwsn_HeatCapacity(*Edata);
	c = Edata->c[TEMPERATURE] * Edata->L * Edata->Rho / (6. * dt);

	// Evaluate the stiffness matrix
	Se[0][0] = Se[1][1] = k; Se[0][1] = Se[1][0] = -k;
	// Heat the element via short-wave radiation
	if ( Edata->sw_abs < 0.0 ) {
		prn_msg(__FILE__, __LINE__, "err", -1., "NEGATIVE Shortwave Radiation %e", Edata->sw_abs);
		return(ERROR);
	}
	Fe[1] += Edata->sw_abs;

	/*
	 * Phase change ennergy:
	 * double Qpc = Edata->Qmf * Edata->L / 4.0; // original: 2.0 Changed by ML: 15.02.2002 as a compromise between Fz and Pb
	 * Heat / Cool element via phasechanges (Edata->Qmf is calculated in phase change module)
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
	Se[0][0] += 2. * c;
	Se[1][1] += 2. * c;
	Se[0][1] += c;
	Se[1][0] += c;
	return(NO_ERROR);
}

/**
 * @brief This routine updates all energy fluxes but solar irradiance \n
 * @par The fluxes are updated before entering the temperature routine sn_SnowTemperature.
 * These fluxes are then used to impose the boundary conditions in sn_Neuman.
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
 * - q_lw = SB*(e*Ta^4 - Ts^4) : Long wave radiation heat transfer \n
 *                               For the implicit solution,this will be treated as a convective
 *                               boundary condition, similar to the sensible heat exchange.
 *                               SB = STEFAN-BOLTZMANN constant.
 * @param *Mdata
 * @param *Bdata
 * @param *Sdata
 * @param *Xdata
 */
void Snowpack::updateMeteoHeatFluxes(const SN_MET_DATA& Mdata, const SN_STATION_DATA& Xdata, 
						    SN_BOUNDARY_DATA& Bdata, SN_SURFACE_DATA& Sdata)
{
	double alpha, gamma;		// The exchange coefficients
	double lw_in, lw_out;		// Incoming and outgoing long wave radiation
	const double& T_air = Mdata.ta;
	const double& Tss = Xdata.Ndata[Xdata.getNumberOfNodes()-1].T;

	alpha = lwsn_SensibleHeat(Mdata, Xdata, HEIGHT_OF_METEO_VALUES);
	Bdata.qs = MIN (350., MAX (-350., alpha * (T_air - Tss)));
	Sdata.qs += Bdata.qs;

	Bdata.ql = lwsn_LatentHeat_Rh(Mdata, Xdata, HEIGHT_OF_METEO_VALUES);
	if ( (Xdata.getNumberOfElements() > 0) && (Xdata.Edata[Xdata.getNumberOfElements()-1].theta[ICE] >= MIN_ICE_CONTENT) ) {
		Bdata.ql = MIN (250., MAX (-250., Bdata.ql));
	}
	Sdata.ql += Bdata.ql;

	if ( T_air > C_TO_K(thresh_rain) ) {
		gamma = Mdata.hnw / M_TO_S(METEO_STEP_LENGTH) * SPECIFIC_HEAT_WATER;
		Bdata.qr = gamma * (T_air - Tss);
	} else {
		Bdata.qr = 0.;
	}
	Sdata.qr += Bdata.qr;

	lw_in  = Constants::stefan_boltzmann * Mdata.ea * T_air*T_air*T_air*T_air;
	lw_out = Constants::stefan_boltzmann * Tss*Tss*Tss*Tss;
	Bdata.lw_out =  lw_out;
	Bdata.lw_net =  lw_in - lw_out;

	Sdata.lw_in  += lw_in;
	Sdata.lw_out += lw_out;
	Sdata.lw_net += Bdata.lw_net;
}

/**
 * @brief Imposes Neumann boundary conditions at the surface. \n
 * The initial fluxes are calculated in sn_MeteoHeatFluxes. \n
 * This splitting ensures that the initial fluxes are used for each iteration in the
 * semi-explicit solution while they are not altered by the implicit solution. \n
 * Note that long wave radiation heat transfer is treated as a convection boundary condition
 * for the implicit solution, similar to the sensible heat exchange: \n
 *            q_lw = delta*(e*Ta - Ts) \n
 * where delta is a function of both Ta and Ts and is calculated by lwsn_LongWave
 * @param Mdata
 * @param Bdata
 * @param Xdata
 * @param T_snow Initial (snow) surface temperature (K)
 * @param T_iter Iterated (snow) surface temperature (K)
 * @param Se Element stiffness matrix
 * @param Fe Element right hand side vector
 */
void Snowpack::sn_Neumann(const SN_MET_DATA& Mdata, SN_BOUNDARY_DATA& Bdata, const SN_STATION_DATA& Xdata, 
					 const double& T_snow, const double& T_iter, 
					 double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ], double Fe[ N_OF_INCIDENCES ])
{
	double alpha, gamma, delta;  //  The heat exchange coefficients
	double S_eb;
	double T_air = Mdata.ta;

	// First zero out the interiour node contribution
	Se[0][0] = Se[0][1] = Se[1][0] = Se[1][1] = Fe[0] = Fe[1] = 0.0;

	/*
	 * Now branch between Phase Change Cases (Semi Explicit Treatment) and Dry Snowpack Dynamics
	 * (Implicit Treatment)
	*/
	if ( (Xdata.Edata[Xdata.getNumberOfElements()-1].theta[WATER] > 0.) && (T_iter != T_snow) ) {
		// Semi-explicit
		// Latent heat transfer and Longwave Radiation have a non-linear dependence on T_snow
		Fe[1] += Bdata.ql + Bdata.lw_net;
		// Sensible and Rain Energy depend linearly on Surface Temperature
		S_eb = Bdata.qs + Bdata.qr;
		if ( (T_air - T_iter) != 0. ) {
			alpha = MAX (0.0, S_eb / (T_air - T_iter));
		} else {
			alpha = 1.0;
		}
		Se[1][1] += alpha;
		Fe[1] += alpha * T_air;
	} else { // Implicit
		// Sensible heat transfer
		alpha = lwsn_SensibleHeat(Mdata, Xdata, HEIGHT_OF_METEO_VALUES);
		Se[1][1] += alpha;
		Fe[1] +=alpha * T_air;
		// Latent heat transfer
		Fe[1] += Bdata.ql;
		// Rain Energy
		if ( T_air > C_TO_K(thresh_rain) ) {
			gamma = Mdata.hnw / M_TO_S(METEO_STEP_LENGTH) * SPECIFIC_HEAT_WATER;
			Se[1][1] += gamma;
			Fe[1] += gamma * T_air;
		}
		// Longwave Radiation
		delta = lwsn_LongWave( T_iter, T_air, Mdata.ea);
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
void Snowpack::sn_NeumannSoil(const double& flux, const double& T_snow, 
						double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ], double Fe[ N_OF_INCIDENCES ])
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
 * @brief Calculates the snow temperatures which are given by the following formula: \n
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
 * @param *Xdata
 * @param *Mdata
 * @param *Bdata
 * @param *mAlb Measured albedo (1)
 */
void Snowpack::sn_SnowTemperature(SN_STATION_DATA& Xdata, SN_MET_DATA& Mdata, SN_BOUNDARY_DATA& Bdata, double& mAlb)
{
	int n, e;        // nodal and element counters
	int iteration;	// iteration counter (not really required)
	int NotConverged;	// = 1 if iteration not converged
	double MaxTDiff;	// maximum temperature difference for convergence
	double TDiff;		// temperature difference for convergence check
	double Big = BIG;	// big number (for DIRICHLET boundary conditions
	// Set the default solution routine convergence parameters
	int MaxItnTemp = 40;                         // maximum 40 iterations for temperature field
	double ConTolTemp = 0.01;	                   // solution convergence to within 0.01 degC

	int Ie[N_OF_INCIDENCES];                     // Element incidences
	double T0[N_OF_INCIDENCES];                  // Element nodal temperatures at t0
	double TN[N_OF_INCIDENCES];                  // Iterated element nodal temperatures
	double Se[N_OF_INCIDENCES][N_OF_INCIDENCES]; // Element stiffnes matrix
	double Fe[N_OF_INCIDENCES];                  // Element right hand side vector

	void *Kt;                                    // Avoids dereferencing the pointer
	SN_ELEM_DATA *EMS;                           // Avoids dereferencing the pointer
	int nN, nE;                                  // Nodes and elements for exposition
	double *U=NULL, *dU=NULL, *ddU=NULL;         // Solution vectors
	double I0;                                   // The net incoming shortwave irradiance
	double Alb;                                  // Modeled albedo
	double hs;                                   // Snow depth
	double res_wat_cont;                         // Snow type dependent residual water content
	double d_pump, v_pump, dvdz;                 // Wind pumping parameters
	double VaporEnhance;                         // Enhancement factor for water vapor transport in wet snow


	// Dereference the pointers
	Kt = Xdata.Kt;
	vector<SN_NODE_DATA>& NDS = Xdata.Ndata;
	nN = Xdata.getNumberOfNodes();
	EMS = &Xdata.Edata[0];
	nE = Xdata.getNumberOfElements();

	// ABSORPTION OF SOLAR RADIATION WITHIN THE SNOWPACK
	// What snow depth should be used?
	if ( ENFORCE_MEASURED_SNOW_HEIGHTS && (Xdata.SlopeAngle < 0.017 * Constants::pi) ) {
		hs = Mdata.hs1;
	} else {
		hs = Xdata.cH;
	}

	// Measured albedo
	if ( (SW_REF > 1) && (Mdata.iswr > 5.) && (Mdata.rswr > 3.) ) {
		mAlb = Mdata.rswr / Mdata.iswr;
	} else {
		mAlb = 0.;
	}

	// Estimated albedo (statistical model) allowing correct treatment of PLASTIC or WET_LAYER
	if ( (nE > Xdata.SoilNode) && (EMS[nE-1].theta[SOIL] < 0.000001) ) {
		if ( (EMS[nE-1].theta[ICE] >= (MIN_ICE_CONTENT)) ) { // Snow on top
			Alb = lwsn_Albedo(EMS[nE-1], NDS[nN-1].T, Mdata, Mdata.date.getJulianDate() - EMS[nE-1].date.getJulianDate());
		} else if ( (nE > Xdata.SoilNode+1) && (EMS[nE-2].theta[ICE] >= (MIN_ICE_CONTENT)) ) {
			// Water over snow or ice
			Alb = lwsn_Albedo(EMS[nE-2], NDS[nN-2].T, Mdata, Mdata.date.getJulianDate() - EMS[nE-2].date.getJulianDate());
		} else {
			Alb = Xdata.SoilAlb;
		}
	} else {
		Alb = Xdata.SoilAlb;
	}

	if ( !(CANOPY && (Xdata.Cdata.height > 3.5)) ) {
		if ( RESEARCH ) { // Treatment of "No Snow" on the ground in research mode
			if ( (hs < 0.02) || (NDS[nN-1].T > C_TO_K(3.5)) ||
				((hs < 0.05) && (NDS[nN-1].T > C_TO_K(1.7))) ) {
				Alb = Xdata.SoilAlb;
			}
		}

		if (!ALPINE3D) {
			Alb = MAX (Alb, Mdata.rswr / Constants::solcon);
		}

		Alb = MAX (Xdata.SoilAlb, MIN(0.95, Alb));
	} else {
		Alb = MAX (0.05, MIN (0.99, Alb));
	}

	switch ( SW_REF ) {
		case 0: case 10: // need to assign rswr a correct value
			Mdata.rswr = Mdata.iswr * Alb;
			break;
		case 1: case 11: // need to assign iswr a correct value
			Mdata.iswr = Mdata.rswr / Alb;
			break;
			// ... while the ground is still snow covered according to HS measurements
		case 2: case 12: // use measured albedo ...
			if ( mAlb > 0. ) {
				if ( (!((mAlb < 2.*Xdata.SoilAlb) && ((Xdata.cH - Xdata.Ground) > 0.05)))
						&& mAlb <= 0.95 ) {
					Alb = mAlb;
				} else {
					Mdata.rswr = Mdata.iswr * Alb;
				}
			} else {
				Mdata.rswr = Mdata.iswr = 0.;
			}
			break;
		default:
			prn_msg(__FILE__, __LINE__, "err", Mdata.date.getJulianDate(), "SW_REF = %d not implemented yet!", SW_REF);
			exit(EXIT_FAILURE);
			break;
	}
	Xdata.Albedo = Alb; // Assign albedo to Xdata
	I0 = Mdata.iswr - Mdata.rswr; // Net irradiance perpendicular to slope
	if ( I0 < 0. ) {
		prn_msg(__FILE__, __LINE__, "err", Mdata.date.getJulianDate(), " iswr:%lf  rswr:%lf  Alb:%lf", Mdata.iswr, Mdata.rswr, Alb);
		exit(EXIT_FAILURE);
	}

	// Simple treatment of radiation absorption in snow: Beer-Lambert extinction (single or multiband).
	try {
		lwsn_ShortWaveAbsorption(Xdata, I0, SNP_SOIL, multistream);
	} catch(exception& ex){
		prn_msg(__FILE__, __LINE__, "err", Mdata.date.getJulianDate(), "Runtime error in sn_SnowTemperature");
		throw;
	}

	// Take care of uppermost soil element
	if ( (nE > Xdata.SoilNode+1) && (EMS[Xdata.SoilNode].sw_abs > EMS[Xdata.SoilNode+1].sw_abs) ) {
		EMS[nE-1].sw_abs += (EMS[Xdata.SoilNode].sw_abs - EMS[Xdata.SoilNode+1].sw_abs);
		EMS[Xdata.SoilNode].sw_abs = EMS[Xdata.SoilNode+1].sw_abs;
	}

	// Set ground surface temperature
	if ( nN == 1 ) {
		// No snow on the ground and !SNP_SOIL
		if ( (Mdata.ts0 - Mdata.ta) > 10. ) {
			NDS[0].T = (Mdata.ts0 + Mdata.ta) / 2.;
		} else {
			NDS[0].T = Mdata.ts0;
		}
		return;
	}

	if ( Kt != 0 ) {
		ds_Solve(ReleaseMatrixData, (MYTYPE*)Kt, 0);
	}
	ds_Initialize( nN, 0, (MYTYPE**)&Kt );
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
	U=(double *) realloc(U, nN*sizeof( double ));
	if ( errno != 0 ) {
		prn_msg(__FILE__, __LINE__, "err", -1., "%s (allocating  solution vector U)", strerror(errno));
		throw IOException("Runtime error in sn_SnowTemperature", AT);
	}
	dU=(double *) realloc(dU, nN*sizeof( double ));
	if ( errno != 0 ) {
		prn_msg(__FILE__, __LINE__, "err", -1., "%s (allocating  solution vector dU)", strerror(errno));
		throw IOException("Runtime error in sn_SnowTemperature", AT);
	}
	ddU=(double *) realloc(ddU, nN*sizeof( double ));
	if ( errno != 0 ) {
		prn_msg(__FILE__, __LINE__, "err", -1., "%s (allocating  solution vector ddU)", strerror(errno));
		throw IOException("Runtime error in sn_SnowTemperature", AT);
	}

  // Make sure that the global data structures know where the pointers are for the next integration step after the reallocation ....
	Xdata.Kt = Kt;

	/*
	 * Set the temperature at the snowpack base to the prescribed value.
	 * NOTE if there is water in the base element then the base temperature MUST be 0 degC
	 * NOTE ts0 no longer set to 0 degC in Control.c
	*/
	if ( !(SNP_SOIL && SOIL_FLUX) ) {
		if ( (EMS[0].theta[ICE] >= MIN_ICE_CONTENT) ) {
			if ( (EMS[0].theta[WATER] > 0.003) ) {
				NDS[0].T = C_TO_K(0.0); /*Mdata.ts0 = C_TO_K(0.0);*/
			} else if ( !SNP_SOIL ) {
				// To avoid temperatures above freezing while snow covered
				NDS[0].T = MIN(Mdata.ts0, C_TO_K(0.0));
			} else {
				NDS[0].T = Mdata.ts0;
			}
		} else {
			NDS[0].T = Mdata.ts0;
		}
	}
	if ( surfaceCode == DIRICHLET_BC ) {
		NDS[nE].T = Mdata.tss;
	}
	// Copy Temperature at time0 into First Iteration
	for (n = 0; n < nN; n++) {
		U[n] = NDS[n].T;
		dU[n] = 0.0;
		ddU[n] = 0.0;
		if( !(U[n] > 50. && U[n] < 500.) ) {
			prn_msg(__FILE__, __LINE__, "err", Mdata.date.getJulianDate(), "Temperature out of bond at beginning of iteration!");
			prn_msg(__FILE__, __LINE__, "msg", -1., "At snow node n=%d (nN=%d): T=%.2lf", n, nN, U[n]);
			free(U); free(dU); free(ddU);
			throw IOException("Runtime error in sn_SnowTemperature", AT);
		}
	}
	// Set the iteration counters, as well as the phase change boolean values
	iteration = 0;
	NotConverged = TRUE;
	// Set the phase change booleans
	Xdata.SubSurfaceMelt = FALSE;
	Xdata.SubSurfaceFrze = FALSE;
	// Determine the displacement depth d_pump
	d_pump = lwsn_WindPumpingDisplacement(Xdata);

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
		if ( nE > Xdata.SoilNode || WIND_PUMP_SOIL ) {
			v_pump = lwsn_WindPumpingVelocity(Mdata, d_pump);
		} else {
			v_pump = 0.0;
		}
		// Assemble matrix
		for (e = nE-1; e >= 0; e--) {
			EL_INCID( e, Ie );
			EL_TEMP( Ie, T0, TN, NDS, U );
			// Update the wind pumping velocity gradient
			dvdz = lwsn_WindGradientSnow(&v_pump, &EMS[e]);
			// Update the water vapor transport enhancement factor
			VaporEnhance = MAX (1., lwsn_EnhanceWaterVaporTransportSnow(Xdata, e));
			// Now fill the matrix
			if ( !sn_ElementKtMatrix(&EMS[e], sn_dt, dvdz, T0, TN, Se, Fe, &Xdata.SubSurfaceMelt, &Xdata.SubSurfaceFrze, VaporEnhance) ) {
				prn_msg(__FILE__, __LINE__, "msg+", Mdata.date.getJulianDate(), "Error in sn_ElementKtMatrix @ element %d:", e);
				for (n = 0; n < nN; n++) {
					fprintf(stdout, "U[%d]=%e K\n", n, U[n]);
				}
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
		 * or volumetric force --- is calculated in sn_ElementKtMatrix().
		*/
		if ( surfaceCode == NEUMANN_BC ) {
			EL_INCID(nE-1, Ie);
			EL_TEMP(Ie, T0, TN, NDS, U);
			sn_Neumann(Mdata, Bdata, Xdata, T0[1], TN[1], Se, Fe);

			ds_AssembleMatrix( (MYTYPE*)Kt, 2, Ie, 2,  (double*) Se );
			EL_RGT_ASSEM( dU, Ie, Fe );
		}
		if ( surfaceCode == DIRICHLET_BC ) {
			Ie[0] = nE;
			ds_AssembleMatrix( (MYTYPE*)Kt, 1, Ie, 1, &Big );
		}
		// Bottom node
		res_wat_cont = lw_SnowResidualWaterContent(Xdata.Edata[0].theta[ICE]);
		// Neumann boundary conditions at bottom: The lower boundary is now a heat flux -- put the heat flux in dU[0]
		if ( SNP_SOIL && SOIL_FLUX ) {
			EL_INCID(0, Ie);
			EL_TEMP(Ie, T0, TN, NDS, U);
			sn_NeumannSoil(GEO_HEAT, T0[1], Se, Fe);

			ds_AssembleMatrix( (MYTYPE*)Kt, 2, Ie, 2, (double*) Se );
			EL_RGT_ASSEM( dU, Ie, Fe );
		}
		/*
		 * Dirichlet boundary condition at bottom: prescribed temperature value
		 * (i.e. known) temperature.
		 * NOTE Insert Big at this location to hold the temperature constant at the prescribed value.
		*/
		else if ( (Xdata.getNumberOfElements() < 3) && (Xdata.Edata[0].theta[WATER] >= 0.9 * res_wat_cont) ) {
			dU[0] = 0.;
		} else {
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
		for (n = 0; n < nN; n++) {
			ddU[n] = dU[n] - ddU[n];
		}
		MaxTDiff = fabs(ddU[0]);
		for (n = 1; n < nN; n++) {
			TDiff = fabs(ddU[n]);
			if ( TDiff > MaxTDiff ) {
				MaxTDiff = TDiff;
			}
		}
		/*
		 * This is to increase the number of iterations for a phase changing uppermost element.
		 * Otherwise we would violate energy conservation because the implicit scheme
		 * leads to an adaptation of surface fluxes to the iterated surface temperature.
		 * In reality, the surface temperature will not change and therefore the fluxes
		 * must be constant. This means that the fluxes must be treated explicitely
		 * (see sn_Neumann)
		 */
		if ( U[nE] + ddU[nE] > MELTING_TK || EMS[nE-1].theta[WATER] > 0. ) {
			ConTolTemp = 0.007;
			MaxItnTemp = 200; // NOTE originally 100;
		}
		NotConverged = ( MaxTDiff > ConTolTemp );
		if ( iteration > MaxItnTemp ) {
			prn_msg(__FILE__, __LINE__, "err", Mdata.date.getJulianDate(), "Temperature did not converge (azi=%.0lf, slope=%.0lf)!", RAD_TO_DEG(Xdata.SlopeAzi), RAD_TO_DEG(Xdata.SlopeAngle));
			prn_msg(__FILE__, __LINE__, "msg", -1., "%d iterations > MaxItnTemp=%d; ConTolTemp=%.4lf; nN=%d", iteration, MaxItnTemp, ConTolTemp, nN);
			for (n = 0; n < nN; n++) {
				if ( n > 0 ) {
					prn_msg(__FILE__, __LINE__, "msg-", -1., "U[%03d]:%6.1lf K, ddU:%8.4lf K;  NDS.T(t-1)=%6.1lf K; EMS[n-1].th_w(t-1)=%.5lf", n, U[n], ddU[n], NDS[n].T, EMS[n-1].theta[WATER]);
				} else {
					prn_msg(__FILE__, __LINE__, "msg-", -1., "U[%03d]:%6.1lf K, ddU:%8.4lf K;  NDS.T(t-1)=%6.1lf K;", n, U[n], ddU[n], NDS[n].T);
				}
			}
			prn_msg(__FILE__, __LINE__, "msg", -1., "Latent: %lf  Sensible: %lf  Rain: %lf  NetLong:%lf  NetShort: %lf", Bdata.ql, Bdata.qs, Bdata.qr, Bdata.lw_net, I0);
			free(U); free(dU); free(ddU);
			throw IOException("Runtime error in sn_SnowTemperature", AT);
		}
		for (n = 0; n < nN; n++) {
			U[n] += ddU[ n ];
		}
	} while ( NotConverged ); // end Convergence Loop

	int crazy = 0;
	for (n = 0; n < nN; n++) {
		if ( (U[n] > t_crazy_min) && (U[n] < t_crazy_max) ) {
			NDS[n].T = U[n];
		} else { // Correct for - hopefully transient - crazy temperatures!
			if ( !crazy && !ALPINE3D && (nN > Xdata.SoilNode + 3) ) {
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date.getJulianDate(), "Crazy temperature(s)! T <= %5.1lf OR T >= %5.1lf; nN=%d; cH=%6.3lf m; azi=%.0lf, slope=%.0lf", t_crazy_min, t_crazy_max, nN, Xdata.cH, RAD_TO_DEG(Xdata.SlopeAzi), RAD_TO_DEG(Xdata.SlopeAngle));
			}
			if ( SNP_SOIL && (n <= Xdata.SoilNode) ) {
				if ( !ALPINE3D ) {
					prn_msg(__FILE__, __LINE__, "msg-", -1., "SOIL node n=%3d: U(t)=%6.1lf  NDS.T(t-1)=%6.1lf K; EMS[n-1].th_w(t-1)=%.5lf", n, U[n], NDS[n].T, EMS[n-1].theta[WATER]);
				}
      } else if ( (n > Xdata.SoilNode) && (nN > Xdata.SoilNode + 3) ) {
				if ( !ALPINE3D ) {
					prn_msg(__FILE__, __LINE__, "msg-", -1., "SNOW node n=%3d: U(t)=%6.1lf  NDS.T(t-1)=%6.1lf EMS[n-1].th_w(t-1)=%.5lf", n, U[n], NDS[n].T, EMS[n-1].theta[WATER]);
				}
			}
      NDS[n].T = 0.5*(U[n] + NDS[n].T);
      crazy++;
		}
	}
	if ( crazy ) {
		prn_msg(__FILE__, __LINE__, "msg-", -1., "%d crazy node(s) from total %d! Latent: %lf  Sensible: %lf  Rain: %lf  NetLong:%lf  NetShort: %lf", crazy, nN, Bdata.ql, Bdata.qs, Bdata.qr, Bdata.lw_net, I0);
	}

	for (e = 0; e < nE; e++) {
		EMS[e].Te = (NDS[e].T + NDS[e+1].T) / 2.0;
		EMS[e].gradT = -(NDS[e].T - NDS[e+1].T) / EMS[e].L;
	}
	free(U); free(dU); free(ddU);
}

/**
 * @brief Assign some results from sn_SnowTemperature() to Sdata.
 * Formerly assigned within sn_SnowTemperature(), but as the latter function may be
 * called twice in a row during one time step ... (Fierz, 15 Feb 2008)
 * @param Xdata
 * @param Mdata
 * @param Sdata
 * @param mAlb Measured albedo (1)
 */
void Snowpack::assignSomeFluxes(const SN_STATION_DATA& Xdata, const SN_MET_DATA& Mdata, const double& mAlb, 
						  SN_SURFACE_DATA& Sdata)
{
	double res_wat_cont; // Snow type dependent residual water content

	// 1) Short wave fluxes and Albedo.
	//     Depending on settings (SW_REF) and conditions,
	//     sw_in and sw_out may differ slightly from the original input
	Sdata.sw_in  += Mdata.iswr;
	Sdata.sw_out += Mdata.rswr;
	Sdata.qw     += Mdata.iswr - Mdata.rswr;

	Sdata.cA += Xdata.Albedo;
	Sdata.mA += mAlb;

	// 2) Ground heat fluxes
	if ( Xdata.getNumberOfElements() > 0 ) {
		res_wat_cont = lw_SnowResidualWaterContent(Xdata.Edata[0].theta[ICE]);

		// 2a) qg0: heat flux at soil/snow boundary
		if ( Xdata.getNumberOfElements() > Xdata.SoilNode ) {
			if ( (Xdata.getNumberOfElements() < 3) && (Xdata.Edata[0].theta[WATER] >= 0.9 * res_wat_cont) ) {
				Sdata.qg += 0.;
			} else {
				Sdata.qg0 += -Xdata.Edata[Xdata.SoilNode].k[TEMPERATURE] *
						Xdata.Edata[Xdata.SoilNode].gradT;
			}
		} else if ( Xdata.SoilNode > 0 ) {
				Sdata.qg0 += -Xdata.Edata[Xdata.SoilNode-1].k[TEMPERATURE] *
						Xdata.Edata[Xdata.SoilNode-1].gradT;
		} else {
			Sdata.qg0 += NODATA;
		}
		// 2b) qg : geothermal heat flux or heat flux at bottom of snow in case of no soil
		if ( SNP_SOIL && SOIL_FLUX ) {
			Sdata.qg += GEO_HEAT;
		} else if ( (Xdata.getNumberOfElements() < 3) && (Xdata.Edata[0].theta[WATER] >= 0.9 * res_wat_cont) ) {
			Sdata.qg += 0.;
		} else {
			Sdata.qg += -Xdata.Edata[0].k[TEMPERATURE] * Xdata.Edata[0].gradT;
		}
	} else {
		Sdata.qg0 += NODATA;
		Sdata.qg  += NODATA;
	}
}

/**
* @brief Various parameterizations of new-snow density
 * @param TA  Air temperature (K)
 * @param TSS Snow surface temperature (K)
 * @param RH  Relative air humidity (1)
 * @param VW  Mean wind velocity (m s-1)
 * @param HH  Altitude a.s.l. (m)
 * @param model Parameterization to be used
 * @return New snow density (kg m-3)
 */
double Snowpack::sn_NewSnowDensityPara(double TA, double TSS, double RH, double VW, double HH)
{
	double rho;

	TA  = K_TO_C(TA);
	TSS = K_TO_C(TSS);
	RH  *= 100.;
	HH  = floor(HH);

	switch(Snowpack::hn_density_model) {
		case LEHNING_OLD:
		{
			double alpha=70., beta=30., gamma=10., delta=0.4;
			double eta=30., phi=6.0, mu=-3.0, nu=-0.5;
			rho = alpha + beta*TA + gamma*TSS +  delta*RH + eta*VW + phi*TA*TSS + mu*TA*VW + nu*RH*VW;
			if ( Snowpack::jordy_new_snow && (VW > 2.9) ) {
				rho = sn_NewSnowDensityHendrikx(TA, TSS, RH, VW);
			}
			break;
		}
		case LEHNING_NEW:
		{
			double alpha=90., beta=6.5, gamma=7.5, delta=0.26;
			double eta=13., phi=-4.5, mu=-0.65, nu=-0.17, om=0.06;
			rho = alpha + beta*TA + gamma*TSS +  delta*RH + eta*VW + phi*TA*TSS + mu*TA*VW + nu*RH*VW + om*TA*TSS*RH;
			// Ad hoc temperature correction
			if (TA < -10.) {
				rho = MIN (rho, alpha*(1. + 0.05*(TA + 10.)));
			}
			// Increase snow density under snow transport conditions
			if ( (!Snowpack::jordy_new_snow) && (VW > 5.) ) {
				rho = 90. + (rho - 30.)*0.9;
			} else if ( Snowpack::jordy_new_snow && (VW > 2.9) ) {
				rho = sn_NewSnowDensityHendrikx(TA, TSS, RH, VW);
			}
			break;
		}
		case BELLAIRE:
		{
			double arg;
			double alpha=3.946, beta=0.07703, zeta=0.0001701, eta=0.02222, mu=-0.05371;
			// Transformations based on natural logarithm!!!
			arg = alpha + beta*TA + zeta*HH + eta*log(VW) + mu*TA*log(VW);
			rho = exp(arg);
			break;
		}
		case ZWART:
		{
			double arg;
			double beta01=3.28, beta1=0.03, beta02=-0.36, beta2=-0.75, beta3=0.3;
			VW = MAX (2., VW);
			RH = 0.8; // ori: MIN(1., RH/100.); see asin(sqrt()) below
			if (TA < -14.) {
				arg = beta01 + beta1*TA + beta2*asin(sqrt(RH)) + beta3*log10(VW);
			} else {
				arg = beta01 + beta1*TA + beta02 + beta2*asin(sqrt(RH)) + beta3*log10(VW); // += beta2*TA;
			}
			rho = pow(10., arg);
			break;
		}
		case PAHAUT:
		{
			rho = 109. + 6.*(C_TO_K(TA) - MELTING_TK) + 26.*sqrt(VW);
			break;
		}
		default:
		{
			prn_msg(__FILE__, __LINE__, "err", -1., "New snow density parameterization (%d) not available", 
				   Snowpack::hn_density_model);
			exit(EXIT_FAILURE);
			break;
		}
	}
	return(MIN(max_hn_density, MAX(MIN_HN_DENSITY, rho)));
}

/**
 * @brief Event driven new-snow density
 * - Antarctica (EVENT_WIND): rho = 250.3 kg m-3 @ 4 m s-1; rho = 338 kg m-3 @ 7 m s-1
 * @param Mdata  Meteorological input
 * @param event_type Wind, relative humidity, etc.)
 * @return New snow density
 */
double Snowpack::sn_NewSnowDensityEvent(const SN_MET_DATA& Mdata, const EventType& i_event_type)
{
	switch( i_event_type ) {
	case Snowpack::EVENT_WIND:
		{
			const double rho_0=361., rho_1=33.;
			if( (Mdata.vw_ave >= event_wind_lowlim) && (Mdata.vw_ave <= event_wind_highlim) ) {
				return (rho_0*log10(Mdata.vw_ave) + rho_1);
			} else {
				return NODATA;
			}
			break;
		}
	default:
		prn_msg(__FILE__, __LINE__,"err", -1.,"No parameterization for event driven nsd for event type %d", i_event_type);
		exit(EXIT_FAILURE);
		break;
	}
}

/**
 * @brief Jordy Hendrikx' new snow density parameterization for strong winds
 * - To be used with Lehning's models only!
 * @param ta  Air temperature (degC)
 * @param tss Snow surface temperature (degC)
 * @param rh  Relative air humidity (%)
 * @param vw  Mean wind velocity (m s-1)
 * @return New snow density
 */
double Snowpack::sn_NewSnowDensityHendrikx(const double ta, const double tss, const double rh, const double vw)
{
	const double alpha=91., beta=-35., gamma=-1.1, delta=49., eta=32.,  phi=4.6;
	return (alpha + beta*ta + gamma*rh +  delta*vw + eta*tss + phi*ta*vw);
}

/**
 * @brief Determines density of new snow
 * @param *Mdata Meteorological input
 * @param *Xdata Snow cover data
 * @param tss    Snow surface temperature (K)
 * @param hnw    Available amount of precipitation (kg m-2)
 * @param model  Model to be used
 * @return double
 */
double Snowpack::calculateNewSnowDensity(const SN_MET_DATA& Mdata, const SN_STATION_DATA& Xdata, const double& tss, const double& hnw, const NewSnowDensityModel& model)
{
	double rho;

	switch(model) {
	case FIXED:
		{
			return fixed_hn_density;
		}
	case MEASURED:
		{
			if ( Mdata.rho_hn > 1. ) { // HACK 1. instead of 0. used with current input for Antarctica
				return Mdata.rho_hn; // New snow density as read from input file
			} else if ( hnw > 0. ) {
				return Xdata.Edata[Xdata.getNumberOfElements()-1].Rho;
			} else {
				return NODATA;
			}
		}
	case EVENT:
		{
			rho = sn_NewSnowDensityEvent(Mdata, event_type);
			return rho;
		}
	default:
		{
			rho = sn_NewSnowDensityPara(Mdata.ta, tss, Mdata.rh, Mdata.vw, Xdata.Alt);
			return rho;
		}
	}
}

/**
 * @brief Determines whether new snow elements are added on top of the snowpack
 * - If ENFORCE_MEASURED_SNOW_HEIGHTS=0 (research mode), new snow height corresponding to cumulated
 *   new snow water equivalent cumu_hnw must be greater than HNS_NE_HEIGHT to be added to Xdata->mH
 * - In case of virtual slopes, uses new snow depth and density from either flat field or luv slope
 * - The first thing is to calculate the height of each element in the snow layer. For now,
 *   instead of trying to find an optimal number of elements, we will define the number of new
 *   elements as a constant nNnew. The constant is related to HNS_NE_HEIGHT, which is set in
 *   qr_ReadParameter(). The smaller HNS_NE_HEIGHT, the more elements we will use.
 * @param *Xdata Snow cover data
 * @param *cumu_hnw cumulated amount of precipitation (kg m-2)
 * @param *Mdata Meteorological data
 * @return ERROR if any, NO_ERROR if not
*/
void Snowpack::determineSnowFall(const SN_MET_DATA& Mdata, SN_STATION_DATA& Xdata, double& cumu_hnw)
{
	int nNewE, nHoarE;          // The number of new elements to be added to SN_STATION_DATA
	int i, e, nE, n, nN, nOldE, nOldN;    // Temporary values; node and element counters
	double z0;                  // Used to determine the z-location of new snowfall nodes
	SN_ELEM_DATA  *EMS;         // Avoids dereferencing the pointer
	double Ln;                  // Original new snow layer element length
	double SigC;                // Used to determine the Cauchy stresses

	double rho_hn=NODATA, t_surf, hn, hoar;  // New snow data
	double cos_sl, L0, dL, Theta0;  // Local values

	nOldN = Xdata.getNumberOfNodes();
	nOldE = Xdata.getNumberOfElements();
	cos_sl = cos(Xdata.SlopeAngle);
	if ( surfaceCode == DIRICHLET_BC ) {
		t_surf = MIN(C_TO_K(-0.1), Mdata.tss);
	} else {
		t_surf = MIN(C_TO_K(-0.1), (Xdata.Ndata[nOldN-1].T + Mdata.ta)/2.);
	}
	rho_hn = calculateNewSnowDensity(Mdata, Xdata, t_surf, cumu_hnw, Snowpack::hn_density_model);

	if ( !ENFORCE_MEASURED_SNOW_HEIGHTS ) {
		if ( Mdata.ta < C_TO_K(thresh_rain) ) {
			if ( (rho_hn != NODATA) && (cumu_hnw/rho_hn > hns_ne_height*cos_sl) ) {
				Xdata.mH += (cumu_hnw/rho_hn);
			}
		} else {
			// This is now very important to make sure that rain will not accumulate
			cumu_hnw -= Mdata.hnw;
		}
	}
	if ( rho_hn == NODATA ) {
		return;
	}
	// Check thresholds to see if a new snow layer has appeared
	// NOTE No new snow during cloud free conditions
	if ( ((Mdata.rh > THRESH_RH) && (Mdata.ta < C_TO_K(thresh_rain)) && (Mdata.ta - Mdata.tss < 3.0)) ||  !ENFORCE_MEASURED_SNOW_HEIGHTS || (Xdata.hn_slope > 0.) ) {
		if ( Xdata.cH < (Xdata.mH - (hns_ne_height*cos_sl)) || (Xdata.hn_slope > 0.) ) {
			// In case of virtual slope use new snow depth and density from either flat field or luv slope
		  if ( Xdata.rho_slope > 0. && (Xdata.SlopeAngle >= 0.017*Constants::pi) ) {
				hn = Xdata.hn_slope;
				rho_hn = Xdata.rho_slope;
			} else { // in case of flat field
				hn = Xdata.mH - Xdata.cH;
				// Store flat field new snow depth and density for virtual slope simulations
				if ( !ALPINE3D && (Xdata.SlopeAngle < 0.017*Constants::pi) ) {
					Xdata.hn_slope = hn;
					Xdata.rho_slope = rho_hn;
				}
			}
			if ( hn > Snowpack::snowfall_warning ) {
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date.getJulianDate(), "Large snowfall! hn=%.3lf cm (azi=%.0lf, slope=%.0lf)", M_TO_CM(hn), RAD_TO_DEG(Xdata.SlopeAzi), RAD_TO_DEG(Xdata.SlopeAngle));
			}
			nNewE = (int)(hn/(hns_ne_height*cos_sl));
			if ( nNewE < 1 ) {
				// Always add snow on virtual slope (as there is no storage variable available) and some other cases
				if ( !ALPINE3D && ((Xdata.SlopeAngle > 0.017*Constants::pi) || (Snowpack::hn_density_model == MEASURED) || (Snowpack::hn_density_model == EVENT)) ) {
					nNewE = 1;
				} else {
					Xdata.hn_slope = 0.;
					return;
				}
			}
			cumu_hnw = 0.0;
			// Check whether surface hoar could be buried
			hoar = Xdata.Ndata[nOldN-1].hoar;
			if ( nOldE > 0 && Xdata.Edata[nOldE-1].theta[SOIL] < 0.00001 ) {
				// W.E. of surface hoar must be larger than a threshold to be buried
				if ( hoar > 1.5*MM_TO_M(MIN_SIZE_HOAR_BURIED)*DENSITY_HOAR_SURF ) {
					nHoarE = 1;
				} else if ( !(change_bc && meas_tss) && (hoar > MM_TO_M(MIN_SIZE_HOAR_BURIED)*DENSITY_HOAR_SURF) ) {
					nHoarE = 1;
				} else {
					nHoarE = 0;
				}
			} else { // Add surface hoar on ground to first new snow element(s)
				nHoarE = 0;
				hn += hoar/rho_hn;
				Xdata.Ndata[nOldN-1].hoar = 0.;
			}

			nN = nOldN + nNewE + nHoarE;
			nE = nOldE + nNewE + nHoarE;

			Xdata.resize(nE);

			vector<SN_NODE_DATA>& NDS = Xdata.Ndata; EMS = &Xdata.Edata[0];

			Xdata.Albedo = Snowpack::new_snow_albedo;
			// Create hoar layer
			if ( nHoarE ) {
				// Since mass of hoar was already added to element below, substract....
				// Make sure you don't try to extract more than is there
				hoar = MAX(0.,MIN(EMS[nOldE-1].M - 0.1,hoar));
				L0 = EMS[nOldE-1].L;
				dL = -hoar/(EMS[nOldE-1].Rho);
				EMS[nOldE-1].L0 = EMS[nOldE-1].L = L0 + dL;

				EMS[nOldE-1].E = EMS[nOldE-1].dE = EMS[nOldE-1].Ee = EMS[nOldE-1].Ev = EMS[nOldE-1].S = 0.0;
				Theta0 = EMS[nOldE-1].theta[ICE];
				EMS[nOldE-1].theta[ICE] *= L0/EMS[nOldE-1].L;
				EMS[nOldE-1].theta[ICE] += -hoar/(Constants::DENSITY_ICE*EMS[nOldE-1].L);
				EMS[nOldE-1].theta[ICE] = MAX(EMS[nOldE-1].theta[ICE],0.);
				EMS[nOldE-1].theta[WATER] *= L0/EMS[nOldE-1].L;
				for (i = 0; i < N_SOLUTES; i++) {
					EMS[nOldE-1].conc[ICE][i] *= L0*Theta0/(EMS[nOldE-1].theta[ICE]*EMS[nOldE-1].L);
				}
				EMS[nOldE-1].M -= hoar;
				EMS[nOldE-1].theta[AIR] = MAX(0., 1.0 - EMS[nOldE-1].theta[WATER] - EMS[nOldE-1].theta[ICE] - EMS[nOldE-1].theta[SOIL]);
				EMS[nOldE-1].Rho = (EMS[nOldE-1].theta[ICE] * Constants::DENSITY_ICE) + (EMS[nOldE-1].theta[WATER] * Constants::DENSITY_WATER) + (EMS[nOldE-1].theta[SOIL]  * EMS[nOldE-1].soil[SOIL_RHO]);
				// Take care of old surface node
				NDS[nOldN-1].z += dL + NDS[nOldN-1].u;
				NDS[nOldN-1].u = 0.0;
				NDS[nOldN-1].hoar = 0.0;
				// Now fill nodal data for upper hoar node
				NDS[nOldN].T = t_surf;              // The temperature of the new node
				// The new nodal position;
				NDS[nOldN].z = NDS[nOldN-1].z + NDS[nOldN-1].u + hoar/DENSITY_HOAR_BURIED;
				NDS[nOldN].u = 0.0;                 // Initial displacement is 0
				NDS[nOldN].hoar = hoar / DENSITY_HOAR_BURIED;         // Surface hoar initial size
				NDS[nOldN].udot = 0.0;               // Settlement rate is also 0
				NDS[nOldN].f = 0.0;                 // Unbalanced forces is 0
				NDS[nOldN].S_n = INIT_STABILITY;
				NDS[nOldN].S_s = INIT_STABILITY;
			} else { // Make sure top node surface hoar mass is removed
				NDS[nOldN-1].hoar = 0.0;
			}
			// Fill the nodal data
			if ( !SNP_SOIL && (nOldN-1 == Xdata.SoilNode) ) { // New snow on bare ground w/o soil
				NDS[nOldN-1].T = (t_surf + Mdata.ta)/2.;
			}
			Ln = (hn / nNewE);
			z0 = NDS[nOldN-1+nHoarE].z + NDS[nOldN-1+nHoarE].u + Ln;
			for (n = nOldN+nHoarE; n < nN; n++) {
				NDS[n].T = t_surf;                  // The temperature of the new node
				NDS[n].z = z0;                      // The new nodal position;
				NDS[n].u = 0.0;                     //   initial displacement is 0
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
				EMS[e].date = Mdata.date;
				// Temperature
				EMS[e].Te = t_surf;
				// Lengths
				EMS[e].L0 = EMS[e].L = (NDS[e+1].z + NDS[e+1].u) - (NDS[e].z + NDS[e].u);
				// Density
				EMS[e].Rho = rho_hn;
				if ( nHoarE && e == nOldE ) {
					EMS[e].Rho = DENSITY_HOAR_BURIED;
				}
				// Mass
				EMS[e].M = EMS[e].L0*EMS[e].Rho;
				// Volumetric components
				EMS[e].theta[ICE]   = EMS[e].Rho/Constants::DENSITY_ICE; // ice content
				EMS[e].theta[WATER] = 0.0;                    // water content
				EMS[e].theta[AIR]   = 1. - EMS[e].theta[ICE]; // void content
				EMS[e].theta[SOIL]  = 0.0;                    // soil content
				for (i = 0; i < N_SOLUTES; i++) {
					EMS[e].conc[ICE][i]   = Mdata.conc[i]*Constants::DENSITY_ICE/Constants::DENSITY_WATER;
					EMS[e].conc[WATER][i] = Mdata.conc[i];
					EMS[e].conc[AIR][i]   = 0.0;
					EMS[e].conc[SOIL][i]  = 0.0;
				}
				// Coordination number based on Bob's empirical function
				EMS[e].N3 = Metamorphism::getCoordinationNumberN3(EMS[e].Rho);
				// Constitutive Parameters
				EMS[e].c[TEMPERATURE] = EMS[e].c[SEEPAGE] = EMS[e].c[SETTLEMENT]= 0.0;
				EMS[e].k[TEMPERATURE] = EMS[e].k[SEEPAGE] = EMS[e].k[SETTLEMENT]= 0.0;
				EMS[e].soil[SOIL_RHO] = EMS[e].soil[SOIL_K] = EMS[e].soil[SOIL_C] = 0.0;
				// Set the initial short wave radiation to zero
				EMS[e].sw_abs = 0.0;
				// Phase change variables:
				EMS[e].dth_w=0.0; // change of water content
				EMS[e].Qmf=0.0;   // change of energy due to phase changes
				// Total element strain (GREEN'S strains -- TOTAL LAGRANGIAN FORMULATION.
				EMS[e].dE = EMS[e].E = EMS[e].Ee = EMS[e].Ev = 0.0;
				// Total Strain Rate	(Simply, E/sn_dt)
				EMS[e].EDot = EMS[e].EvDot=0.0;
				// Total Element Stress
				EMS[e].S=0.0;
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
					if ( value > 1.0 ) { // Graupel
						EMS[e].mk = 4;
						EMS[e].dd = 0.;
						EMS[e].sp = 1.;
						EMS[e].rg = 0.6;
						EMS[e].rb = 0.2;
						// Because density and volumetric contents are already defined, redo it here
						EMS[e].Rho = 110.;
						EMS[e].theta[ICE] = EMS[e].Rho/Constants::DENSITY_ICE;  // ice content
						EMS[e].theta[AIR] = 1. - EMS[e].theta[ICE];  // void content
					} else { // no Graupel
						EMS[e].mk = Snowpack::new_snow_marker;
						if ( Snowpack::jordy_new_snow && (Mdata.vw > 2.9) && (Snowpack::hn_density_model < BELLAIRE) ) {
							alpha = 1.87; beta = -0.04;
							EMS[e].dd = MAX(0.5,MIN(1.0,(alpha + beta*Mdata.vw)*(alpha + beta*Mdata.vw)));
							EMS[e].sp = new_snow_sp;
							alpha = 0.9; beta = 0.015; gamma = -0.0062; delta = -0.117; eta=0.0011; phi=-0.0034;
							EMS[e].rg = MIN(new_snow_grain_rad, MAX(0.3*new_snow_grain_rad, alpha + beta*TA + gamma*RH + delta*Mdata.vw + eta*RH*Mdata.vw + phi*TA*Mdata.vw));
							EMS[e].rb = 0.4*EMS[e].rg;
						} else {
							EMS[e].dd = new_snow_dd;
							EMS[e].sp = new_snow_sp;
							// Adapt dd and sp for blowing snow
							if ( (Mdata.vw > 5.) && ((VARIANT == ANTARCTICA) || (!Snowpack::jordy_new_snow && ((Snowpack::hn_density_model == BELLAIRE) || (Snowpack::hn_density_model == LEHNING_NEW)))) ) {
								EMS[e].dd = new_snow_dd_wind;
								EMS[e].sp = new_snow_sp_wind;
							} else if ( vw_dendricity && ((Snowpack::hn_density_model == BELLAIRE) || (Snowpack::hn_density_model == ZWART)) ) {
								double vw;
								vw = MAX(0.05, Mdata.vw);
								// dd = f(vw)
								EMS[e].dd = (1 - pow(vw/10., 1.57));
								EMS[e].dd = MAX(0.2, EMS[e].dd);
							}
							if ( Snowpack::hydrometeor ) { // empirical
								alpha=1.4; beta=-0.08; gamma=-0.15; delta=-0.02;
								EMS[e].rg = 0.5*(alpha + beta*TA + gamma*Mdata.vw + delta*TA*Mdata.vw);
								EMS[e].rb = 0.25*EMS[e].rg;
							} else {
								EMS[e].rg = new_snow_grain_rad;
								EMS[e].rb = new_snow_bond_rad;
								if ( ((Mdata.vw_ave >= event_wind_lowlim) && (Mdata.rh_ave >= rh_lowlim)) ) {
									EMS[e].rb = MIN(bond_factor_rh*EMS[e].rb, Metamorphism::max_grain_bond_ratio*EMS[e].rg);
								}
							}
						}
					} // end no Graupel
					if ( nHoarE &&	 e == nOldE ) {
						EMS[e].mk = 3;
						EMS[e].dd = 0.;
						EMS[e].sp = 0.;
						EMS[e].rg = MAX(new_snow_grain_rad, 0.5*M_TO_MM(EMS[e].L0));
						// Note: L0 > MIN_SIZE_HOAR_BURIED/DENSITY_HOAR_BURIED
						EMS[e].rb = EMS[e].rg/3.;
					}
				} // Treat all the initial snow types

				// Classify the grain type
				EMS[e].type = ml_ag_Classify(EMS[e].dd, EMS[e].sp, 2.*EMS[e].rg,
						EMS[e].mk%100, EMS[e].theta[WATER], EMS[e].theta[ICE]);

				// Initialise the Stability Index for ml_st_CheckStability routine.)
				EMS[e].S_dr = INIT_STABILITY;
				EMS[e].hard = 0.;
			}   // End elements

			// Find the Cauchy stress
			for (e = nE-1, SigC = 0.0; e >= 0; e--) {
				SigC += -EMS[e].M*Constants::g*cos_sl;
				EMS[e].C = SigC;
			}
			// Finally, update snowpack height
			Xdata.cH = Xdata.mH = NDS[nN-1].z + NDS[nN-1].u;
			Xdata.ErosionLevel = nE-1;
		} //  End of NEW_SNOWFALL Condition 2
	} //  End of NEW_SNOWFALL Condition 1
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
 * -# If the calculated temperatures are above zero degrees then PHASECHANGES are occuring.
 *    This means that the volumetric contents of the elements must be updated. >>>>>> PHASE CHANGE MODULE <<<<<
 * -# In the next step the CREEPING deformations and 1d state of stress within the snowpack
 * at each exposition are found.  If SNOWFALL and SNOWDRIFT has occured then the finite
 * element matrices (connectivities) must be rebuilt.  >>>>CREEP MODULE<<<<<
 * @param *Mdata
 * @param *Xdata
 * @param *cumu_hnw Variable to store amount of precipitation (kg m-2)
 * @param *Bdata
 * @param *Sdata
 */
void Snowpack::runSnowpackModel(SN_MET_DATA& Mdata, SN_STATION_DATA& Xdata, double& cumu_hnw, 
						  SN_BOUNDARY_DATA& Bdata, SN_SURFACE_DATA& Sdata)
{
	double mAlb = 0.;    // To transfer measured albedo from sn_SnowTemperature() to Sdata
	WaterTransport watertransport(cfg);
	Metamorphism metamorphism(cfg);
	SnowDrift snowdrift(cfg);
	PhaseChange pc(cfg);

	try {

		// Adjust Boundary Condition
		if ( change_bc && meas_tss ) {
			if ( Mdata.tss < C_TO_K(THRESH_CHANGE_BC) ) {
				surfaceCode = DIRICHLET_BC;
			} else {
				surfaceCode = NEUMANN_BC;
			}
		}

		// If it is SNOWING, find out how much, prepare for new FEM data
		determineSnowFall(Mdata, Xdata, cumu_hnw);

		/*
		 * Check to see if snow is DRIFTING, calculate a simple snowdrift index and erode layers if
		 * neccessary. Note that also the very important friction velocity is calculated in this
		 * routine and later used to calculate the Meteo Heat Fluxes
		 */
		snowdrift.calcSnowDrift(Mdata, Xdata, Sdata, cumu_hnw);

		// Reinitialize and Calculate the Meteo Heat Fluxes
		memset((&Bdata), 0, sizeof(SN_BOUNDARY_DATA));
		updateMeteoHeatFluxes(Mdata, Xdata, Bdata, Sdata);
		
		// Find the temperature in the snowpack
		sn_SnowTemperature(Xdata, Mdata, Bdata, mAlb);

		//Good HACK (according to Charles)... like a good hunter and a bad one...
		// If you switched from DIRICHLET to NEUMANN boundary conditions, correct
		//   for a possibly erroneous surface energy balance. The latter can be due e.g. to a lack
		//   of data on nebulosity leading to a clear sky assumption for incoming long wave.
		if ( (change_bc && meas_tss) && (surfaceCode == NEUMANN_BC)
			&& (Xdata.Ndata[Xdata.getNumberOfNodes()-1].T < C_TO_K(THRESH_CHANGE_BC)) ) {

			surfaceCode = DIRICHLET_BC;
			Xdata.Ndata[Xdata.getNumberOfNodes()-1].T = C_TO_K(THRESH_CHANGE_BC/2.);

			sn_SnowTemperature(Xdata, Mdata, Bdata, mAlb);
		}

		// Assign some results from sn_SnowTemperature() to Sdata
		assignSomeFluxes(Xdata, Mdata, mAlb, Sdata);
		
		// See if any SUBSURFACE phase changes are occuring
		pc.runPhaseChange(Sdata, Xdata);

		/*
		 * The water transport routines must be place here, otherwise the temperature
		 * and creep solution routines will not pick up the NEW mesh boolean.
		 */
		watertransport.transportMass(Mdata, Bdata.ql, Xdata, Sdata);


		/*
		 * Find the settlement of the snowpack.
		 * HACK This routine was formerly placed here because the settlement solution MUST ALWAYS follow
		 * sn_SnowTemperatures where the vectors U, dU and dUU are allocated.
		 */
		calcSnowCreep(Mdata, Xdata);

	} catch(exception& ex) {
		prn_msg(__FILE__, __LINE__, "err", Mdata.date.getJulianDate(), "Snowpack Calculation Not Completed");
		throw;
	}

	metamorphism.runMetamorphismModel(Mdata, Xdata);
		
	if (join_elements)
		Xdata.joinElements(NUMBER_TOP_ELEMENTS);
}

/*
 * End of Snowpack.c
*/
