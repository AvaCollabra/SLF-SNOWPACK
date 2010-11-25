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
 * @file Metamorphism.c
 * @brief This module contains the snow metamorphism routines of the SLF one-dimensional OPERATIONAL and RESEARCH snowpack model \n
 * It represents a truly international research effort: Dr. Bob "Borolo" Brown of MONTANA STATE UNIVERSITY (USA)
 * provided a lot of the motivation and low temperature gradient micro-structure
 * physics; Dr. Michael "Give me a Girl" Lehning (GERMANY), when he was was not
 * busy on the telephone talking to his dentist in San Francisco, provided project
 * leadership, a sense of the practical and more importantly, the link to the
 * avalanche warning group, i.e. he told everybody what to do; Dr. Pramod Sataywali
 * "Porky" (INDIA) came up with the high temperature gradient micro-structure
 * routines, while, at the same time, missing his wife and child terribly, skiing
 * and dreaming of Indian spin bowlers. Dr. Perry Bartelt (PLANET PLUTO) was FORCED
 * to write the code and INTEGRATE it into a very sensible continuum mechanics
 * model, that works pretty well.  He does NOT accept any RESPONSIBLITY for the
 * VERY STRANGE physical contants that the Borolo Bob and Sataywali use: SNOW MICRO-
 * STRUCTURE is BULLSHIT, BULLSHIT, BULLSHIT, and more BULLSHIT.
 * Michael had written all the Metamorphism stuff up to document it and he started
 * to understand what Perry meant with his lines above and below. He started
 * cleaning the code the day he knew that Betty had received the price for the best
 * presentation and the most scientific content at the Innsbruck conference. He was
 * proud. The main task will be to get rid of all the trash in this routine such as
 * cgs thanks to Porky and use the functions that are already programmed such as
 * saturation vapor pressure. Of course, university professors can not know that
 * we need saturation vapor pressure not only for Metamorphism.......
 * Definition of some of the essential variables:
 *
 * PRIMARY micro-structure parameters calculated by Metamorphism routine:
 *
 * - rb : bond radius in [mm]
 * - rg : grain radius in [mm]
 * - dd : snow dendricity (=0 old snow; =1 new snow)
 * - sp : snow sphericity (=0 faceted; =1 rounded)
 * - mk : microstructure marker:
 *                            - 0 dendritic snow or neither faceted nor rounded
 *                            - 1 reached a sp of 0
 *                            - 2 reached a sp of 1
 *                            - 3 Surface Hoar SH
 *                            - 4 Graupel PPgp
 *                            - 5 not implemented yet --> thin crusts
 *                            - 8 Ice layer IFil
 *                            - mk < 10, mk=mk+10 : first complete wetting
 *                            - mk < 20, mk=mk+10 : first melt-freeze cycle completed
 *                            - mk / 100 >= 1     : tagged snow layer
 *
 * SECONDARY micro-structure parameters calculated by Metamorphism routine:
 * - N3   : coordination number (1)
 *
 * SECONDARY micro-structure parameters calculated by routines found in Laws_sn.c:
 * - keff : effective conductivity (W m-1 K-1)
 * - E    : modulus of elasticity  (Pa)
 * - eta  : viscosity (Pa s)
 *
 * These are the variables needed to calculate these values, contained in the
 * SN_ELEM_DATA of the sn_Xdata structure:
 * - theta[ICE]   : volumetric ice content (1)
 * - theta[WATER] : volumetric water  (1)
 * - Te           : temperature (K)
 * - dTdZ         : temperature gradient (K m-1)
 * - dPdZ         : vapor pressure gradient (bar m-1)
 * - Rho          : bulk density (kg m-3)
 * - S            : overburden stress (Pa)
 */
/* The code is dedicated to Perry's heros:  General Vo Nguyen Giap, author of
* "People's War, People's Army: The Viet Cong Insurrection Manual for
* Underdeveloped Countries", which Perry is now using as a self-help guide to
* to survive the SLF; the poets Allen Ginsberg and W.H. Auden and, of course,
* Brooks Robinson, who Perry saw on a beautiful September day in 1971 hit two
* journeyman singles to centerfield against the Cleveland Indians: I fell in love.
* The works of these great men will be cited throughout the code. We begin with
* the General: "Such was the essence of the strategic direction of the Dien Bien
* Phu campaign and of the WINTER-SPRING campaign as a whole.  This direction drew
* its inspiration from the principles of DYNAMISM, INITIATIVE, MOBILITY and
* RAPIDITY of decision in face of NEW SITUATIONS.  Its main objective was the
* DESTRUCTION of enemy MANPOWER.  It took full advantage of the CONTRADICTIONS in
* which the enemy was involved and developed to the utmost the spirit of active
* offensive of the revolutionary army. This CORRECT, CLEAR-SIGHTED and BOLD
* strategy enabled us to deprive the enemy of all possiblility of retrieving the
* initiative, and to create favourable conditions for us to fight a decisive battle
* on a battlefield CHOSEN and PREPARED for by US.  This strategic direction ensured
* the success of the whole WINTER-SPRING campaign which was CROWNED by the great
* victory of DIEN BIEN PHU."
                                                                                              |
* This code is OUR Dien Bien Phu:    a victory against FRENCH imperialism,
* JAPANESE fascism and AMERICAN interventionism, i.e bad French metamorphism
* ideas, funny Japanese viscosity functions and, of course, bad American university
* programmaying!!  (This guy has a serious problem, folks, criticizing the home of
* the (DE)brave (d) and land of the free (BIE), or whatever that studpid phrase is
* The french metamorphism routines were written in November 1995 by Perry Bartelt
* and Martin Schneebeli.  They were first used in the 2d snowpack code haefeli.
*/

#include <snowpack/Metamorphism.h>
#include <snowpack/Snowpack.h>

using namespace std;
using namespace mio;

/************************************************************
 * static section                                           *
 ************************************************************/

///Defines the vapor pressure gradient at which TG metamorphism begins (hPa m-1)
const double Metamorphism::mm_tg_dpdz = 5.;
const double Metamorphism::ba_g_fudge = 3.; ///< brief Defines Thorsten's Geometry FUDGE
const double Metamorphism::sa_g_fudge = 0.35; ///< Defines Satyawali's Geometry FUDGE 

/// @name Grain and bond growth
//@
const double Metamorphism::max_grain_growth = 5.; ///< A grain cannot grow more than 5.0 mm d-1

///Bond radius usually will not grow larger than Metamorphism::bond_size_stop * grain radius
const double Metamorphism::bond_size_stop = 0.75;

///Absolute limit of grain to bond ratio (SH, initialisation)
const double Metamorphism::max_grain_bond_ratio = 0.95;
//@}

///@name Thresholds for wind slab formation
//@
///For no action at all, set strength factor for wind slab formation to 0.0
const double Metamorphism::wind_slab_enhance = 5.;

///Wind slab forms for winds stronger than Metamorphism::wind_slab_vw (m s-1)
const double Metamorphism::wind_slab_vw = 3.;

/// Wind slab formation down to Metamorphism::wind_slab_depth (m)
const double Metamorphism::wind_slab_depth = 0.07;
//@}


map<string, MetaModelFn> Metamorphism::mapMetamorphismModel;
map<string, MetaSpRateFn> Metamorphism::mapSpRate;
const bool Metamorphism::__init = Metamorphism::initStaticData();

bool Metamorphism::initStaticData()
{
	mapMetamorphismModel["DEFAULT"] = &Metamorphism::metamorphismDEFAULT;
	mapMetamorphismModel["NIED"]    = &Metamorphism::metamorphismNIED;

	mapSpRate["DEFAULT"] = &Metamorphism::spRateDEFAULT;
	mapSpRate["NIED"]    = &Metamorphism::spRateNIED;

	return true;
}

/**
 * @brief This routine estimates the cross sectional pore area
 * @author Charles Fierz
 * @date 2009-12-23
 * @param Edata
 * @return area (mm2)
 */
double Metamorphism::csPoreArea(const SN_ELEM_DATA& Edata)
{
	const double rg = MM_TO_M(Edata.rg); // Grain and bond radius (m)

	return((Constants::pi * rg*rg / Edata.theta[ICE]) * (1 - Edata.theta[ICE] - Edata.theta[WATER]));
}

/**
 * @brief This routine estimates the coordination number as a function of the snow density,
 * used for both the FRENCH and MONTANA metamorphism models.  The MONTANA model
 * uses the coordination number to determine the thermal conductivity of snow. The
 * constants are based on EXPERIMENTAL results. This is  a TYPICAL MONTANA State
 * University routine -- a bunch of defines with a lot of unusual constants, a bunch
 * of multiplications and then a return.  And they expect us to believe this stuff.
 * @param Rho Snow density (kg m-3)
 * @return Coordination number (1)
 */
 /* This reminds me of a poem by Auden:  "The History of Tru.h>
 *
 * In that ago when being was believing,
 * TRUTH was the most of many CREDIBLES,  (UNLIKE this routine which is UNCREDIBLE)
 * More first, more always, than a bat-winged lion,
 * A fish-tailed dog or eagle-headed fish, (i.e, fantasy=this routine)
 * The least like mortals, doubted by their deaths. (the programmer should die)
 *
 * Truth was their MODEL as they strove to BUILD (we are doing this too)
 * A world of LASTING objects to believe in,(I would like to believe in this stuff)
 * Without believing earthernware and legend,
 * Archway and song, were truthful or untruthful,
 * The TRUTH was there already to be TRUE.
 *
 * This while when, PRACTICAL like paper-dishes, (This routine is, however, useful)
 * TRUTH is CONVERTIBLE to KILOWATTS,
 * Our last to do by is an ANTI-MODEL,
 * SOME UNTRUTH anyone can give the LIE to,
 * A NOTHING no one need believe is there.  (I wish this routine wasn't here)
 *
 * This poem basically sums up Perry's feeling about modelling micro-structure.
 */
double Metamorphism::getCoordinationNumberN3(const double& Rho)
{
	double N_0 = 1.4153;
	double N_1 = 7.5580e-5;
	double N_2 = 5.1495e-5;
	double N_3 = 1.7345e-7;
	double N_4 = 1.8082e-10;

	double R_2, R_3, R_4, N3;
	R_2 = Rho*Rho;
	R_3 = R_2*Rho;
	R_4 = R_3*Rho;

	// For Rho between 100 kg/m3 and 670 kg/m3, use the following.
	N3 = N_0 - N_1*Rho + N_2*R_2 - N_3*R_3 + N_4*R_4;

	// Outside these limits use the following.
	if ( Rho >= 670. ) {
		N3 = 10.5;             // upper limit on N3 near close-packing value for spheres.
	}
	if ( Rho <= 100.0 ) {
		N3 = 1.75*(Rho/100.);  // Decreases N3 to zero as density goes to zero.
	}

	return(N3);
}

/**
 * @brief The following routine calculates the rate of change of dendricity according to
 * the Montana State Model. This routine was programmed on a foggy day in middle
 * March 1998 by Perry and Bob in Michael's office.
 * @param Edata const SN_ELEM_DATA
 * @return Rate of change (d-1)
 */
double Metamorphism::ddRate(const SN_ELEM_DATA& Edata)
{
	double ddDot;
	double Tk, dTdZ;
	double c, f;

	// Get the data at time t from the element data structures
	dTdZ = fabs(Edata.gradT);
	Tk = Edata.Te;
	c = exp (-6000. / Tk); // original CROCUS: -6000.; set to -5800. by Bellaire 2004; back to ori 2007
	f = c * dTdZ / 7.; // Introduced by Ml on 11 Oct 2004; original CROCUS: c*pow(dTdZ, 0.4)

	// Calculate the change in dendricity
	if ( dTdZ < 5.0 ) {
		ddDot = -2.5e8*c;   // Ml: ori -2.0; set to -5.0e8 by Bellaire 2004, then to -2.5e8 2007;
	} else if( dTdZ < 15.0 ) {
		ddDot = -3.5e8*f;   // Ml: ori -4.0; set to -1.5e8 by Bellaire 2004, then to -3.5e8 2007;
	} else {
		ddDot = -3.0e8*f;
	}

	return(ddDot);
}

/************************************************************
 * non-static section                                       *
 ************************************************************/

Metamorphism::Metamorphism(const mio::Config& i_cfg) : cfg(i_cfg) 
{
	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Parameters");
	sn_dt = M_TO_S(calculation_step_length);

	new_snow_grain_rad = cfg.get("new_snow_grain_rad", "Parameters");

	string tmp_metamorphism_model = cfg.get("METAMORPHISM_MODEL", "Parameters");
	metamorphism_model = tmp_metamorphism_model;

	IOUtils::toUpper(metamorphism_model);

	map<string, MetaModelFn>::const_iterator it1 = mapMetamorphismModel.find(metamorphism_model);
	if (it1 == mapMetamorphismModel.end()) throw InvalidArgumentException("Unknown metamorphism model: "
															+metamorphism_model, AT);

	map<string, MetaSpRateFn>::const_iterator it2 = mapSpRate.find(metamorphism_model);
	if (it2 == mapSpRate.end()) throw InvalidArgumentException("Unknown metamorphism model: "+metamorphism_model, AT);
}

/**
 * @brief Rate of change for sphericity sp
 * Presently, it is based on the french parameterization (CROCUS)
 * @param Edata
 * @return Rate of change (d-1)
 */
double Metamorphism::spRateDEFAULT(const SN_ELEM_DATA& Edata)
{
	double spDot;
	double dTdZ;
	double c, f;

	dTdZ = fabs(Edata.gradT);
	c = exp(-6000./Edata.Te); // Original 6000.
	f = c*dTdZ/17.; // Introduced by ml on 2004-10-11; original CROCUS: c*pow(dTdZ, 0.4)

	if ( dTdZ < 5.0 ) {
		if ( Edata.sp < 0.49 ) {
			spDot = (5.0 - dTdZ)*0.2e8*c*Edata.sp;
		} else {
			spDot = (5.0 - dTdZ)*0.5e8*c;
		}
	} else {
		spDot = -((dTdZ - 5.)/dTdZ)*0.5e8*f; // ML: Ori 2.0
	}
	// Consider graupel separately
	if ( (Edata.mk%100) == 4 ) {
		spDot = 0.;
	}

	return(spDot);
}

/**
 * @brief Rate of change for sphericity sp
 * Based on the french parameterization (CROCUS) but adapted to Japanese conditions
 * @param *Edata
 * @return Rate of change (d-1)
 */
double Metamorphism::spRateNIED(const SN_ELEM_DATA& Edata)
{
	double spDot;
	double dTdZ;
	double c, f;

	dTdZ = fabs(Edata.gradT);
	c = exp(-6000./Edata.Te); // Original 6000.
	f = c*dTdZ/17.; // Introduced by ml on 2004-10-11; original CROCUS: c*pow(dTdZ, 0.4)

  if ( dTdZ < 15.0 ) { //NIED (H. Hirashima)
    spDot = (15.0 - dTdZ)*0.5e9*c;
	} else if ( (dTdZ > 15.0) && (dTdZ < 25.0) ) {
    spDot = -1.e8*f;
	} else {
    if ( Edata.sp > 0. ) {
      spDot =  -1.e8*f;
    } else {
      spDot = 0.0;
    }
  }
	// Now consider Graupel
	if ( (Edata.mk%100) == 4 ) {
		spDot = 0.;
	}

	return(spDot);
}

/**
 * @brief This is Satyawali's TG bond growth rate routine
 * @param *Edata
 * @return Bond radius growth rate (mm d-1)
 */
double Metamorphism::TGBondRate(const SN_ELEM_DATA& Edata)
{
	double A;          // average cross sectional area (m2)
	double TGradBond;  // micro temp gradient across bonds (K m-1)
	double flux;       // mass flux of vapor in the pore space - in cgs units
	double rbDot;      // Bond radius growth rate (mm d-1)

	const double rb = MM_TO_M(Edata.rb);    // initial bond radius (m)
	const double rg = MM_TO_M(Edata.rg);    // initial grain radius (m)
	const double TGrad = fabs(Edata.gradT); // absolute value of temp gradient within element (K m-1)

	A = 1./3. * (Constants::pi*(rb*rb + rg*rg) + csPoreArea(Edata));
	TGradBond = Edata.k[TEMPERATURE] / Constants::conductivity_ice * A / (Constants::pi * rb*rb) * (-TGrad);       // (K m-1) NOTE Why take TGrad neg.?
	flux = -Constants::diffusion_coefficient_in_air / (Constants::gas_constant * Edata.Te*Edata.Te) * (Constants::lh_sublimation / (Constants::gas_constant * Edata.Te) - 1) * TGradBond;
	flux *= lw_SaturationPressure(Edata.Te); // (kg s-1 m-2)
	// Bond radius growth rate (m s-1)
	rbDot = flux / Constants::density_ice * Metamorphism::sa_g_fudge;
	// Convert to mm d-1
	return(M_TO_MM(D_TO_S(rbDot)));
}

/**
 * @brief This function calculates the lattice constant (mm)
 * Note there is one source and one sink grain (1+1)
 * @param th_ice Volumetric ice content (1)
 * @return Lattice constant (mm)
*/
double Metamorphism::LatticeConstant0(const double& th_ice)
{
	double gsz0 = 2*new_snow_grain_rad;

	return( pow((1. + 1.)*Metamorphism::ba_g_fudge*gsz0*gsz0*gsz0/th_ice, 1./3.) );
}

/*
 * Kinetic grain growth routine
*/
/**
 * @brief Actual routine: mm_TGGrainRateImplicit from r7.7
 * Fierz' implicit version of Baunach's model
 * (see SnpR1.fz/Originaux/Metamorphism.fzV5.c) :
 * Formulation of graingrowth including intra- and layer-to-layer transport
 * --> see paper in Ann. Glaciol., 32 (IGS Innsbruck 2000).
 * lattice constant given as a function of both gsz-gsz0 and density (th_i),
 * assuming the overall snow density and thus theta_ice are not much
 * affected by L2L, which is quite true! (Ml2l ~ 1mg/element/time step.)
 * Note : gsz0 = 2*new_snow_grain_rad
 * @param *Edata
 * @param Tbot Lower node temperature (K)
 * @param Ttop Upper node temperature (K)
 * @param gradTSub Temperature gradient at lower node (K m-1)
 * @param gradTSup Temperature gradient at upper node (K m-1)
 * @return Grain radius growth rate (mm d-1)
 */
double Metamorphism::TGGrainRate(const SN_ELEM_DATA& Edata, const double& Tbot, const double& Ttop,
                                 const double& gradTSub, const double& gradTSup)
{
	double intraFlux;                      // Intra-layer vapor flux
	double botFlux;                        // Vapor flux into the element
	double topFlux;                        // Vapor flux out of the element
	double th_i;                           // Ice content
	double dFluxL2L;                       // Flux divergence due to L2L transport
	double Te;                             // Temperature (K)
	double gradTbot, gradT, gradTtop;      // Temperature gradients (K m-1)
	double gsz;                            // grain size (mm)
	double a;                              // lattice constant (mm)
	double hElem;                          // Element height (mm)
	double a0, a1;                         // Empirical regression coefficients to estimate the lattice constant, f(gsz)
	double reg0 = 0.15, reg1 = -0.00048;   // Empirical regression coefficients to estimate the lattice constant, f(th_i)
	double rgDot=0.;                       // Grain radius growth rate (mm d-1)

	// Collect the continuum values from the element data structures
	th_i       = Edata.theta[ICE];
	Te         = Edata.Te;
	gradT      = Edata.gradT;
	gradTbot   = (gradTSub + Edata.gradT)/2.;
	gradTtop   = (Edata.gradT + gradTSup)/2.;
	gsz        = 2 * (Edata.rg);
	hElem      = M_TO_MM(Edata.L);

	// Calculate the lattice constant a at time t but a <= hElem;  Units : mm
	a0 = LatticeConstant0( th_i );
	if ( gsz > 2*new_snow_grain_rad ) {
		// Use an empirical estimation of the lattice constant
		a1 = reg0 + reg1*(th_i * Constants::density_ice);
		a  = a0 + a1*(gsz - 2*new_snow_grain_rad);
	} else {
		a = a0;
	}
	a  = MIN (a, hElem);

	// Intra layer flux, where the direction of flow does not matter! Units: kg/(sm2)
	intraFlux =  fabs(Constants::diffusion_coefficient_in_snow / (Constants::gas_constant * Te*Te) * (Constants::lh_sublimation / (Constants::gas_constant * Te) - 1.) * gradT);
	intraFlux *= lw_SaturationPressure(Te);

	// Layer to layer flux, where the direction of flow DOES matter! Units: kg/(sm2)
	botFlux = - Constants::diffusion_coefficient_in_snow / (Constants::gas_constant * Tbot*Tbot) * (Constants::lh_sublimation / (Constants::gas_constant * Tbot) - 1.) * gradTbot;
	botFlux *= lw_SaturationPressure(Tbot);
	topFlux = - Constants::diffusion_coefficient_in_snow / (Constants::gas_constant * Ttop*Ttop) * (Constants::lh_sublimation / (Constants::gas_constant * Ttop) - 1.) * gradTtop;
	topFlux *= lw_SaturationPressure(Ttop);
	dFluxL2L = -(topFlux - botFlux);
	// Calculate the rate in m s-1
	rgDot = 0.5 * ( (intraFlux + dFluxL2L * (a / hElem) ) * a*a) / (2.0 * Metamorphism::ba_g_fudge * Constants::density_ice * (2 * new_snow_grain_rad) * gsz);
	// Conversion to mm d-1
	rgDot = M_TO_MM(D_TO_S(rgDot));

	return(rgDot);
}


/**
 * @return Below  is Borolo Bob's ET bond growth rate routine.  Determines the bond or neck
 * growth for low temperature gradients. Called from the routine
 * mm_Metamorphism. Note that the growth rate is converted from mm s-1 to
 * mm d-1 before returning rbDot.
 * @param *Edata
 * @return Bond radius growth rate (mm d-1)
 */
double Metamorphism::ETBondRate(const SN_ELEM_DATA& Edata)
{
	/*
	* B_1...B_3 are  constants calculated with Brown's advanced and sophisticated
	* mixture theory. Bartelt is so jealous of that fine piece of work.   Please note
	* hist sarcastic tirade later in this  unreadable program.
	*/
	const double B_1 = 0.1436e-3;         //  in mm/sec
	const double B_2 = -1.8850e-6;        //  in mm
	const double B_3 = 4.6690e+3;         //  deg K
	const double B_R = 273.;
	const double rc = lwsn_ConcaveNeckRadius(Edata.rg, Edata.rb);
	double rbDot; // Bond radius growth rate (mm s-1)

	if ( fabs(rc - Edata.rb) < Constants::eps ) {	//special case: thermodynamic neck radius rn is infinite
		rbDot  = B_1 * (1. - exp(B_2/Edata.rg) ) * exp( (B_3/B_R) - (B_3/Edata.Te)  );
	} else {
		const double rn = ( 2.*Edata.rb*rc) / (rc - Edata.rb);
		rbDot  = B_1 * (exp(B_2/ rn ) - exp(B_2/Edata.rg) ) * exp( (B_3/B_R) - (B_3/Edata.Te)  );
	}
	// Convert to mm d-1
	rbDot = D_TO_S(rbDot);

	return(rbDot);
}


/**
 * @brief Below is Borolo Bob's ET grain growth rate routine : Determines the grain growth
 * for low temperature gradients.  Called from the routine mm_Metamorphism.
 * Note that the growth ate is converted from mm/sec to mm/day before returning rgDot
 * @param *Edata
 * @return Grain radius growth rate (mm s-1)
 */
double Metamorphism::ETGrainRate(const SN_ELEM_DATA& Edata)
{
	double rgDot; // Grain radius growth rate (mm s-1)

	// These are the routine's FUDGE FACTORs
	const double C_1 = 9.403e-11;
	const double C_2 = 5.860e-9;
	const double C_3 = 2.900e3;
	const double C_R = 273.;

	rgDot = ((C_1 / Edata.rg) + C_2) * exp((C_3 / C_R) - (C_3 / Edata.Te));
	// Convert to mm d-1
	rgDot = D_TO_S(rgDot);

	return( rgDot );
}


/**
 * @brief Below  is Borolo Bob's UNIVERSAL LAW  (equivalent to Einstein's UNIVERSAL
 * FIELD LAW, Bob made  me  type this in to impress people.)
 * @param *Edata
 * @return Additional bond radius growth rate (mm s-1)
 */
double Metamorphism::PressureSintering(const SN_ELEM_DATA& Edata)
{
	double rbdot; // Bond radius growth rate (mm s-1)

	/*
	 * If there is little or no ice, then RETURN
	 *    There are several reasons for setting the limit at 0.05. First, if we have
	 *    dry snow with an ice fraction less than 5%, this is NEW snow which the
	 *    current program does not model and for which the dendricity is high and the
	 *    sphericity is low. If it is not new snow, the only other case would be
	 *    when it is WET, and by the time the ice fraction is this low, the bonding
	 *    between grains should be negligible, since grain-to-grain contacts of
	 *    spherical ice grains (by the time they have been wet long enough to reduce
	 *    the ice mass to this low value, the grains are spherical) has disappeared.
	*/
	if ( Edata.theta[ICE] <= Constants::min_ice_content ) {
		return (0.0);
	}

	/*
	 * If the temperatures are over zero, then also RETURN (this should not occur, but
	 * to be SAFE .... in fact the only time it can occur is when theta_i=00. (Actually
	 * it can occur that theta_i > 0 and T == 0. The next piece of code picks this
	 * case up....
	*/
	if ( Edata.Te > Constants::melting_tk ) {
		return (0.0);
	}

	/*
	 * Introduced by Ml 11 Oct 2004:
	 * Relate bond growth to total deformation calculated with actual (partly empirical) viscosity;
	 * Previously this was done with microstructure viscosity, which is only used for main part of
	 * settling, however.
	*/
	rbdot = -0.1 * Edata.rb * Edata.EvDot / lwsn_Neck2VolumetricStrain(Edata);

	// Convert from (mm s-1) to (mm d-1)
	rbdot = D_TO_S(rbdot);

	return(rbdot);
}


/**
 * @brief Main routine for Metamorphism model
 * Actual routine: bb_mm_Metamorphism from r7.7
 * @param Mdata
 * @param Xdata
 */
void Metamorphism::metamorphismDEFAULT(const SN_MET_DATA& Mdata, SN_STATION_DATA& Xdata)
{
	int e, nE;
	SN_ELEM_DATA *EMS;   // Avoids dereferencing the element data pointer

	double rgDot;        // Grain growth rate (mm d-1)
	double rbDot;        // Bond growth rate (mm d-1)
	double rgDotMax, rbDotMax;  // Maximum grain and bond growth rates
	double ddDot;        // Rate of dendricity change (d-1)
	double spDot;        // Rate of sphericity change (d-1)
	double splim1, splim2, splim3;  // Constants used to limit changes in sphericity after faceting
	double dDay;         // Time increment in days, specifically written to avoid confusion
	double P1, P2;       // Nodal pressures of element e
	double T1, T2;       // Nodal temperatures of element
	double dPdZ;         // Vapor pressure gradient within element
	const double a1 = 1.11e-3, a2 = 3.65e-5;  // mm3 day-1 Volumetric growth coefficients for wet snow
	int    marker;       // local variable

	// Dereference the element pointer containing micro-structure data
	EMS = &Xdata.Edata[0];
	//vector<SN_ELEM_DATA>& EMS = Xdata.Edata;
	vector<SN_NODE_DATA>& NDS = Xdata.Ndata;
	nE = Xdata.getNumberOfElements();

	const double cw = 1.e8 * exp(-6000. / 273.15);

	// Calculate for each snow element within the mesh the micro-structure changes
	for (e = Xdata.SoilNode; e < nE; e++) {
		// Set all rates of change to zero for element e
		ddDot = spDot = rbDot = rgDot = 0.0;

		if ( EMS[e].theta[ICE] < 0.00001 || EMS[e].theta[SOIL] > 0.00001 ) {
			continue;
		}

		// Determine the coordination number which is purely a function of the density
		EMS[e].N3 = getCoordinationNumberN3(EMS[e].Rho);

		// Compute local values
		double thetam_w = 1.e2 * (Constants::density_water * (EMS[e].theta[WATER]) / (EMS[e].Rho));

		splim1 = 20. * (new_snow_grain_rad - EMS[e].rg);
		if ( splim1 > 0.0 ) {
			splim1=0.0;
		}
		splim2 = EMS[e].sp / 0.5;
		if ( splim2 > 1.0 ) {
			splim2 = 1.0;
		}
		splim3 = -0.7;
		marker = EMS[e].mk%100;  // untag EMS[e].mk

		// Calculate the pressure gradient (kinetic or equilibrium growth metamorphism??)
		T1 = NDS[e].T;
		T2 = NDS[e+1].T;
		P1 = lw_SaturationPressure(T1);
		P2 = lw_SaturationPressure(T2);
		dPdZ = fabs((P2 - P1) / EMS[e].L) * 0.01;  //  Result is in hPa m-1

		// Equilibrium growth rates for old dry snow
		rgDot = ETGrainRate(EMS[e]);
		rbDot = ETBondRate(EMS[e]);

		// Kinetic growth rates
		// Since we need temperature gradients above and below the element we have to consider various cases for the kinetic grain growth
		if ( e > 0 && e < nE-1 ) { // inner element
			rgDotMax = TGGrainRate(EMS[e], T1, T2, EMS[e-1].gradT, EMS[e+1].gradT);
		} else if ( e == 0 ) {// bottom element: use twice EMS[e].gradT to avoid troubles if nE=1
			rgDotMax = TGGrainRate(EMS[e], T1, T2, EMS[e].gradT, EMS[e].gradT);
		} else {// top element
			rgDotMax = TGGrainRate(EMS[e], T1, T2, EMS[e-1].gradT, EMS[e].gradT);
		}
		rgDotMax = MAX(0.0, rgDotMax);
		rbDotMax = TGBondRate(EMS[e]);

		if ( (EMS[e].theta[WATER] < 0.01) && (Mdata.vw > Metamorphism::wind_slab_vw) && ((NDS[nE].z - NDS[e].z < Metamorphism::wind_slab_depth) || e == nE-1) ) {
			//if snow is dry AND wind strong AND we are near the surface => wind densification of snow
			// Introduce heuristic metamorphism for wind slabs of Metamorphism::wind_slab_depth (m)
			double wind_slab = 1.;

			// Enhancement factor; see also sn_SnowCreep()
			wind_slab += Metamorphism::wind_slab_enhance * (Mdata.vw - Metamorphism::wind_slab_vw);

			ddDot = wind_slab * ddRate(EMS[e]);
			spDot = wind_slab * (CALL_MEMBER_FN(*this, mapSpRate[metamorphism_model])(EMS[e]));
			rgDot = 0.;
			rbDot = 0.5 * rgDotMax; //HACK Why should it be half the kinetic rate?
		} else {
			//normal processes for snow
			// NEW SNOW
			if ( EMS[e].dd > 0.0 ) {
				// WET new snow
				if ( EMS[e].theta[WATER] > 0.01 ) { //NIED if(EMS[e].theta[WATER] > 0.1) CORRECTED SINCE version 7.4
					ddDot = -pow(thetam_w, 3.) / 16.;
					if ( (-ddDot) < cw ) {
						ddDot = -cw;
					}
					spDot = -0.5 * ddDot;

					rgDot = rbDot = 0.0; // no grain growth until dd <= 0.0
				} else {
					// DRY new snow
					ddDot = ddRate(EMS[e]);
					spDot = CALL_MEMBER_FN(*this, mapSpRate[metamorphism_model])(EMS[e]);
					if ( (EMS[e].dd < 0.8) && (dPdZ > 2.*Metamorphism::mm_tg_dpdz) ) {
						// TG metamorphism (new snow); mimicks the Japanese change according to Fierz
						rgDot += rgDotMax;
					} else {
						// no grain growth until dd <= 0.0 if gradient too small
						rgDot += 0.; // (dPdZ/2./Metamorphism::mm_tg_dpdz*rgDotMax);
					}
				}
			} else { // (OLD) SNOW
				// WET snow
				if ( EMS[e].theta[WATER] > 0.005 ) {
					ddDot = 0.0;
					spDot = pow(thetam_w, 3.) / 16.;
					if ( spDot < 2.*cw ) {
						spDot = 2.*cw;
					}
					// Faceted grains, dry and wet, need first to be rounded (sp > 0.5) before they grow due to the presence of liquid water.
					if ( (marker%10 == 2) || EMS[e].sp > 0.5 ) {
						rgDot = 1. / (4. * Constants::pi * EMS[e].rg * EMS[e].rg) * (a1 + a2 * pow(thetam_w, 3.));
						rbDot = 0.6 * rgDot;
					} else {
						rgDot = rbDot = 0.;
					}
				} else {
					// DRY snow
					ddDot = 0.0;
					spDot = CALL_MEMBER_FN(*this, mapSpRate[metamorphism_model])(EMS[e]);
					if ( fabs(EMS[e].gradT) < 5.0 ) {
						if ( (marker == 1 || marker == 3) && EMS[e].rg > 0.3 ) {
							spDot = spDot * exp(splim1);
						} else if ( (marker == 21 || marker == 23) && EMS[e].sp < 0.5 ) {
							spDot = spDot*(splim2);
						}
					} else if ( marker == 22 ) {
						spDot = spDot * exp(splim3);
					}
					if ( dPdZ > Metamorphism::mm_tg_dpdz ) { // HACK Why did we sum up both growth rates?
						rbDot = rbDotMax;
						rgDot = rgDotMax;
					} else { // HACK Why are grains growing more than predicted by Borolo Bob?
						rbDot += ( (dPdZ / Metamorphism::mm_tg_dpdz + 0.1) * rbDotMax );
						rgDot += ( (dPdZ / Metamorphism::mm_tg_dpdz + 0.1) * rgDotMax );
					}
				}
			}
		}

		// Pressure sintering: update bond growth rate, both wet and dry snow
		EMS[e].ps2rb = PressureSintering(EMS[e]);
		rbDot += EMS[e].ps2rb;

		// UPDATE THE MICROSTRUCTURE PARAMETERS
		// Time increment in days, specifically written to avoid confusion
		dDay = S_TO_D(sn_dt);
		// Update dendricity
		EMS[e].dd += ddDot * dDay;
		EMS[e].dd = MAX(0.0, MIN (1.0, EMS[e].dd));
		// Update sphericity
		EMS[e].sp += spDot * dDay;
		if ( (marker == 1) && (EMS[e].rg >= 0.4) ) {
			EMS[e].sp = MAX(0.0, MIN(0.5, EMS[e].sp)); // Limit effect of rounding on dry faceted grains
		} else {
			EMS[e].sp = MAX(0.0, MIN(1.0, EMS[e].sp));
		}
		// Update grain size ...
		rgDot = MIN(rgDot, Metamorphism::max_grain_growth);
		if ( marker != 3 ) {
			EMS[e].rg += rgDot*dDay;
		} else {
			//HACK ... but do not allow surface hoar to grow and limit its size to layer thickness.
			EMS[e].rg = MIN(EMS[e].rg, 0.5 * M_TO_MM(EMS[e].L));
		}
		// Update bond size and limit its growth to Metamorphism::bond_size_stop * EMS[e].rg
		rbDotMax = (Metamorphism::bond_size_stop * EMS[e].rg - EMS[e].rb) / dDay;
		rbDot = MAX(0., MIN(rbDot, rbDotMax));
		EMS[e].rb += rbDot * dDay;
		if ( marker == 3 ) { //HACK SH is only grain allowed to decrease its grain size!
			EMS[e].rb = MIN(EMS[e].rb, Metamorphism::max_grain_bond_ratio * EMS[e].rg);
		}

		// Compute proportion of grain bond growth due to pressure sintering
		if ( (EMS[e].dd < 0.005) && (rbDot > 0.) ) {
			EMS[e].ps2rb /= rbDot;
		} else {
			EMS[e].ps2rb = 0.0;
		}
		// Update the microstructure marker
		if ( EMS[e].dd < 0.001 ) { //NIED EMS[e].dd < 0.001
			if ( (EMS[e].sp < 0.1) && (marker % 10 == 0) ) {
				EMS[e].mk += 1;  // grains become fully faceted
			}
			if ( (EMS[e].sp > 0.999) && (marker % 10 == 0) ) {
				EMS[e].mk += 2;  // grains become fully rounded
			}
			if ( (EMS[e].theta[ICE] > 0.763) && (marker % 10 != 8) ) {
				EMS[e].mk = (EMS[e].mk / 10) * 10 + 8;
				// There is an ice layer forming for dry densities above 700 kg m-3!
			}
		}

		// First wetting
		if ( (EMS[e].theta[WATER] > 0.015) && (marker < 10) ) {
			// Non-dendritic snow: thrsh ori 0.3 changed by S.Bellaire to get thinner crusts (13.03.2006)
			// Dendritic snow: very rapid change to melt forms
			if ( (EMS[e].theta[WATER] > 0.35*lw_SnowResidualWaterContent(EMS[e].theta[ICE])) || (marker < 1) ) {
				EMS[e].mk += 10;
			}
		}
    // First melt-freeze cycle completed
		else if ( (marker < 20) && (marker >= 10) && (EMS[e].Te < Constants::melting_tk - 0.3) ) {
			EMS[e].mk += 10;
		}

		// Update the calculation of grain class.
		EMS[e].type = ml_ag_Classify(EMS[e].dd, EMS[e].sp, 2. * EMS[e].rg, EMS[e].mk % 100, EMS[e].theta[WATER], EMS[e].theta[ICE]);
	} // end of loop over snow elements
}

/**
 * @brief Main routine for Metamorphism model adapted according to NIED (H. Hirashima)
 * @param Mdata
 * @param Xdata
 */
void Metamorphism::metamorphismNIED(const SN_MET_DATA& Mdata, SN_STATION_DATA& Xdata)
{
	int e, nE;
	SN_ELEM_DATA *EMS;   // Avoids dereferencing the element data pointer

	double rgDot;        // Grain growth rate (mm d-1)
	double rbDot;        // Bond growth rate (mm d-1)
	double rgDotMax, rbDotMax;  // Maximum grain and bond growth rates
	double ddDot;        // Rate of dendricity change (d-1)
	double spDot;        // Rate of sphericity change (d-1)
	double splim1, splim2, splim3;  // Constants used to limit changes in sphericity after faceting
	double dDay;         // Time increment in days, specifically written to avoid confusion
	double P1, P2;       // Nodal pressures of element e
	double T1, T2;       // Nodal temperatures of element
	double dPdZ;         // Vapor pressure gradient within element
	const double a1 = 1.11e-3, a2 = 3.65e-5;  // mm3 day-1 Volumetric growth coefficients for wet snow
	double cw, thetam_w; // local variables
	// double res_wat_cont; //NIED local variable not used
	int    marker;       // local variable
	double dhfDot;       //NIED (H. Hirashima) Depth hoar factor ... //Fz HACK needs to be initialized! see line 894!
	double DenFact, Diffus, gradV; //NIED (H. Hirashima) //Fz HACK please describe variables

	// Dereference the element pointer containing micro-structure data
	EMS = &Xdata.Edata[0];
	vector<SN_NODE_DATA>& NDS = Xdata.Ndata;
	nE = Xdata.getNumberOfElements();

	// Calculate for each snow element within the mesh the micro-structure changes
	for (e = Xdata.SoilNode; e < nE; e++) {
		// Set all rates of change to zero for element e
		ddDot = spDot = rbDot = rgDot = 0.0;

		if ( EMS[e].theta[ICE] < 0.00001 || EMS[e].theta[SOIL] > 0.00001 ) {
			continue;
		}

		// Determine the coordination number which is purely a function of the density
		EMS[e].N3 = getCoordinationNumberN3(EMS[e].Rho);

		// Compute local values
		cw = 1.e8 * exp(-6000. / 273.15);
		thetam_w = 1.e2 * (Constants::density_water * EMS[e].theta[WATER] / EMS[e].Rho);
		splim1 = 20. * (new_snow_grain_rad - EMS[e].rg);
		if ( splim1 > 0.0 ) {
			splim1=0.0;
		}
		splim2 = EMS[e].sp / 0.5;
		if ( splim2 > 1.0 ) {
			splim2 = 1.0;
		}
		splim3 = -0.7;
		marker = EMS[e].mk % 100;  // untag EMS[e].mk

		// Calculate the pressure gradient (kinetic or equilibrium growth metamorphism??)
		T1 = NDS[e].T;
		T2 = NDS[e+1].T;
		P1 = lw_SaturationPressure(T1);
		P2 = lw_SaturationPressure(T2);
		dPdZ = fabs((P2 - P1) / EMS[e].L) * 0.01;  //  Result is in mbar m-1

		// Equilibrium growth rates for old dry snow
		rgDot = ETGrainRate(EMS[e]);
		rbDot = ETBondRate(EMS[e]);

		// Kinetic growth rates
		// Since we need temperature gradients above and below the element we have to consider various cases for the kinetic grain growth
		if ( e > 0 && e < nE-1 ) { // inner element
			rgDotMax = TGGrainRate(EMS[e], T1, T2, EMS[e-1].gradT, EMS[e+1].gradT);
		} else if ( e == 0 ) { // bottom element: use twice EMS[e].gradT to avoid troubles if nE=1
			rgDotMax = TGGrainRate(EMS[e], T1, T2, EMS[e].gradT, EMS[e].gradT);
		} else { // top element
			rgDotMax = TGGrainRate(EMS[e], T1, T2, EMS[e-1].gradT, EMS[e].gradT);
		}
		rgDotMax = MAX (0.0, rgDotMax);
		rbDotMax = TGBondRate(EMS[e]);

		if ( (EMS[e].theta[WATER] < 0.01) && (Mdata.vw > Metamorphism::wind_slab_vw) && ((NDS[nE].z - NDS[e].z < Metamorphism::wind_slab_depth) || e == nE-1) ) {
			//if snow is dry AND wind strong AND we are near the surface => wind densification of snow
			// Introduce heuristic metamorphism for wind slabs of Metamorphism::wind_slab_depth (m)
			double wind_slab = 1.;

			// Enhancement factor; see also sn_SnowCreep()
			wind_slab += Metamorphism::wind_slab_enhance * (Mdata.vw - Metamorphism::wind_slab_vw);

			ddDot = wind_slab * ddRate(EMS[e]);
			spDot = wind_slab * (CALL_MEMBER_FN(*this, mapSpRate[metamorphism_model])(EMS[e]));
			rgDot = 0.;
			rbDot = 0.5 * rgDotMax;
			//dhfdot = 0.; //Fz HACK You'd need to define dhfdot in this case also
		} else {
			//normal processes for snow
			// NEW SNOW
			if ( EMS[e].dd > 0.0 ) {
				// WET new snow
				if ( EMS[e].theta[WATER] > 0.01 ) { //NIED if(EMS[e].theta[WATER] > 0.1) CORRECTED SINCE version 7.4
					ddDot = -pow(thetam_w, 3.) / 16.;
					if ( (-ddDot) < cw ) {
						ddDot = -cw;
					}
					spDot = -0.5 * ddDot;

					rgDot = rbDot = 0.0; // no grain growth until dd <= 0.0
					dhfDot = -pow(thetam_w,3.)/16./86400.; //NIED (H. Hirashima)
				} else {
					// DRY new snow //NIED (H. Hirashima)
					ddDot = ddRate(EMS[e]);
					spDot = -(fabs(EMS[e].gradT) - 15.)/17.; //NIED (H. Hirashima)
					rgDot = rbDot = 0.0; // no grain growth until dd <= 0.05

					if (fabs(EMS[e].gradT)>150.) { //NIED (H. Hirashima)
						ddDot=-0.5;
						spDot=-0.25;
					}
					gradV=dPdZ*7.93E-4;  //NIED (H. Hirashima) hPa/m��kg/m2�ɕϊ�
					DenFact = -0.136*EMS[e].Rho+4.56;
					Diffus = MAX((2.23E-5*(1013.25/1013.25)*pow((EMS[e].Te)/273.15,1.78)),((0.78*(EMS[e].Te-273.15))+10.84)*1.0E-5); //NIED (H. Hirashima)
					dhfDot = fabs(-DenFact*Diffus*gradV*(1.0-EMS[e].dhf));
					if (fabs(EMS[e].gradT)<5.0) {
						dhfDot=-60000000.*exp(-6000./EMS[e].Te)/86400.;  //NIED (H. Hirashima)
					}
				}
			} else { // (OLD) SNOW
				// WET snow
				if ( EMS[e].theta[WATER] > 0.005 ) {
					ddDot = 0.0;
					spDot = pow(thetam_w, 3.) / 16.;
					if ( spDot < 2.*cw ) {
						spDot = 2.*cw;
						if (spDot<-(fabs(EMS[e].gradT) - 15.)/17.) { //NIED (H. Hirashima)
							spDot = -(fabs(EMS[e].gradT) - 15.)/17.;
						}
					}
					// Faceted grains, dry and wet, need first to be rounded (sp > 0.5) before they grow due to the presence of liquid water.
					if ( (marker%10 == 2) || EMS[e].sp > 0.5 ) {
						rgDot = 1. / (4. * Constants::pi * EMS[e].rg * EMS[e].rg) * (a1 + a2 * pow(thetam_w, 3.));
						rbDot = 0.6 * rgDot;
						dhfDot = -(pow(thetam_w,3.)/16./86400.);
						if ( dhfDot>-2.*cw/86400. ) {  //NIED (H. Hirashima)
							dhfDot=-2.*cw/86400.;
						}
					} else {
						rgDot = rbDot = 0.;
					}
				} else {
					// DRY snow  //NIED (H. Hirashima)
					ddDot = 0.0;
					spDot = CALL_MEMBER_FN(*this, mapSpRate[metamorphism_model])(EMS[e]);
					gradV=dPdZ*7.93E-4; //NIED (H. Hirashima) //hPa/m��kg/m2�ɕϊ�
					DenFact = -0.136*EMS[e].Rho+4.56;  //NIED (H. Hirashima)
					Diffus = MAX((2.23E-5*(1013.25/1013.25)*pow((EMS[e].Te)/273.15,1.78)),((0.78*(EMS[e].Te-273.15))+10.84)*1.0E-5); //NIED (H. Hirashima)
					dhfDot = fabs(-DenFact*Diffus*gradV*(1.0-EMS[e].dhf));
					if ( fabs(EMS[e].gradT)<5.0 ) {
						dhfDot=-500000000.0*exp(-6000.0/EMS[e].Te)*(5-fabs(EMS[e].gradT))/86400; //NIED (H. Hirashima)
					}
					if ( dPdZ > Metamorphism::mm_tg_dpdz ) {
						rbDot = TGBondRate( EMS[e] );
						// Since we need temperature gradients above and below the element we have to be careful for the grain growth
						if (e > 0 && e < nE-1) {
							rgDot = TGGrainRate( EMS[e], T1, T2, EMS[e-1].gradT, EMS[e+1].gradT );
						} else if ( e == 0 ) {
							rgDot = TGGrainRate( EMS[e], T1, T2, EMS[e].gradT, EMS[e+1].gradT );
						}	else {
							rgDot = TGGrainRate( EMS[e], T1, T2, EMS[e-1].gradT, EMS[e].gradT );
							// rgDot = mm_TGGrainRate( &EMS[e] );  Thorstens Formulation
						}
						if ( rgDot < 0.0 ) {
							rgDot = 0.0;
						}
						if ( fabs(EMS[e].gradT ) > 150. ) { //NIED (H. Hirashima)
							ddDot=-0.5;
							spDot=-0.25;
						}
						rgDot=(6.25E-12 * EMS[e].gradT + 6.48E-10)*86400.*1000.; //NIED (H. Hirashima)
					}
				}
			}
		}

		// Pressure sintering: update bond growth rate, both wet and dry  snow
		EMS[e].ps2rb = PressureSintering(EMS[e]);
		rbDot += EMS[e].ps2rb;

		// Time increment in days, specifically written to avoid confusion (see line 1012
		dDay = S_TO_D(sn_dt);

		// UPDATE THE MICROSTRUCTURE PARAMETERS
		if(EMS[e].theta[WATER] > 0.01 ) { //NIED (H. Hirashima)
			dhfDot = -(pow(thetam_w,3.)/16./86400.);
			if(dhfDot>-2.*cw/86400.) {
				dhfDot=-2.*cw/86400.;
			}
			if (EMS[e].dd == 0.) {
				dhfDot=dhfDot/2.; // HACK //Fz Hazardous comparison!
			}
		}
		EMS[e].dhf += dhfDot * sn_dt; //NIED (H. Hirashima) HACK //Fz use consistent units dDay instead of sn_dt
		EMS[e].dhf = MAX(0.0, MIN(1.0, EMS[e].dhf)); //NIED (H. Hirashima)
		// Update dendricity
		EMS[e].dd += ddDot * dDay;
		EMS[e].dd = MAX (0.0, MIN (1.0, EMS[e].dd));
		// Update sphericity
		EMS[e].sp += spDot * dDay;
		if ( (marker == 1) && (EMS[e].rg >= 2.) ) { //NIED (H. Hirashima)
			EMS[e].sp = MAX(0.0, MIN(0.5, EMS[e].sp)); // Limit effect of rounding on dry faceted grains
		} else {
			EMS[e].sp = MAX(0.0, MIN(1.0, EMS[e].sp));
		}
		// Update grain size ...
		//rgDotMax = Metamorphism::max_grain_growth;
		rgDot = MIN(rgDot, Metamorphism::max_grain_growth);
		if ( marker != 3 ) {
			EMS[e].rg += rgDot*dDay;
		} else {
			// ... but do not allow surface hoar to grow and limit its size to layer thickness.
			EMS[e].rg = MIN (EMS[e].rg, 0.5 * M_TO_MM(EMS[e].L));
		}
		// Update bond size
		rbDotMax = (Metamorphism::bond_size_stop * EMS[e].rg - EMS[e].rb) / dDay;
		rbDot = MAX (0., MIN (rbDot, rbDotMax));
		EMS[e].rb += rbDot * dDay;
		// Compute proportion of grain bond growth due to pressure sintering
		if ( (EMS[e].dd < 0.005) && (rbDot > 0.) ) {
			EMS[e].ps2rb /= rbDot;
		} else {
			EMS[e].ps2rb = 0.0;
		}
		// Update the microstructure marker
		if ( EMS[e].dd < 0.001 ) { //NIED EMS[e].dd < 0.001
			if ( (EMS[e].sp < 0.1) && (marker % 10 == 0) ) {
				EMS[e].mk += 1;  // grains become fully faceted
			}
			if ( (EMS[e].sp > 0.999) && (marker % 10 == 0) ) {
				EMS[e].mk += 2;  // grains become fully rounded
			}
			if ( (EMS[e].theta[ICE] > 0.763) && (marker % 10 != 8) ) {
				EMS[e].mk = (EMS[e].mk / 10) * 10 + 8;
				// There is an ice layer forming for dry densities above 700 kg m-3!
			}
		}

		// First wetting //NIED (H. Hirashima)
		if ( (marker < 10) && (EMS[e].theta[WATER] > 0.99*lw_SnowResidualWaterContent(EMS[e].theta[ICE])) ) {
			EMS[e].mk += 10;
		}
    // First melt-freeze cycle completed
		else if ( (marker < 20) && (marker >= 10) && (EMS[e].Te < Constants::melting_tk - 0.3) ) {
			EMS[e].mk += 10;
		}

		// Update the calculation of grain class.
		EMS[e].type = ml_ag_Classify(EMS[e].dd, EMS[e].sp, 2. * EMS[e].rg, EMS[e].mk % 100, EMS[e].theta[WATER], EMS[e].theta[ICE]);
//Fz Structure has no longer members named ng, nb; also note that mg and mb are both (hopefully) 0., so what!
//	  EMS[e].ng = (EMS[e].M - Constants::density_water * EMS[e].theta[WATER] * EMS[e].L) / (mg + 0.5 * EMS[e].N3 * mb);
//	  EMS[e].nb = 0.5 * N3 * EMS[e].ng;
	} // end of loop over snow elements
}

void Metamorphism::runMetamorphismModel(const SN_MET_DATA& Mdata, SN_STATION_DATA& Xdata) throw()
{
	CALL_MEMBER_FN(*this, mapMetamorphismModel[metamorphism_model])(Mdata, Xdata);
}

/*
 * End of Metamorphism.c
 */
