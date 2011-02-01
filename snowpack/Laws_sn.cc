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
 * @file Laws_sn.c
 * @version 10.03
 * @brief This module contains (ALL) the constitutive laws for the 1d snowpack model.
 *
 * This includes Michael Lehning's ALBEDO law as well as Bob Brown's VISCOSITY and HEAT CONDUCTION models.
 * It also contains the SENSIBLE and LATENT heat CONVECTION coefficients.  Moreover, this code
 * is used by ALL snowpack modules -- creep, temperature and water transport.  Many of the
 * laws, for example the ELASTICITY modulus were programmed up by Marc Christen and taken
 * directly from program HAEFELI, the 2d snowpack model.
 */
/* Warning: The following paragraph is a good example of how jealousy is able to impair a
 * otherwise fairly well working mind. The motivation for the paragraph stems from the fact
 * that NOBODY walked up to Perry after HIS presentation to tell HIM that HIS talk was the
 * best of the conference and that impressive work has been done (even though this
 * compliment was just as well directed at HIS and the whole party's address); not to mention
 * that no funky girl would even consider chatting HIM up whether Italian or Spanish.
 * (By the way, several people came up to Perry and told him, D. McClung, K. Lied, P. Foehn and
 * and and ..Michael himself.. )
 *
 * Perry cleaned the code up after returning from Norway where Michael presented the first
 * snowpack paper (Let's go to the first sheet!).  His presentation was good but could have
 * been better:  he spent the entire conference chatting up a good looking Italian girl -- NO,
 * sorry, a good looking Italian girl chatted him up -- which of course ruined
 * his concentration.  Michael, in a word, BROKE one of the serious plegdes of snowpack
 * programmers (see Main.c, PLEDGE no. 7). This piece of code, and also Michael's
 * behaviour, demonstrate that there is a serious difference between "the possible", (die
 * Moeglichkeit), e.g. constitutive laws that work and funky Italian girls, and "reality" (die
 * Wirklichkeit), e.g. a semi-stable code and steep granite walls without any borrebolte:
 *
 * MOEGLICHKEIT
 * "Hoeher als die Wirklichkeit steht die Moeglichkeit.  Das Moeglichsein, das je das Dasein
 * existenzial ist, unterscheidet sich ebensosehr von der LEEREN, logischen Moeglichkeit wie
 * von der Kontingenz eines Vorhandenen.  Die Moeglichkeit als Existenzial ist die
 * urspruenglischste und LETZTE POSITIVE ONTOLOGISCHE BESTIMMTHEIT des DASEINS."
 *
 * Basically, what Heidegger is saying is that possibilities are a lot of fun, and why playing
 * with different constitutive laws (and women) is central to human existence and also
 * snow research. Possibilites are important paraphenalia along the WAY (der Weg), for
 * example, to Norway:
 *
 * WEG (unterwegs)
 * "Dasein ist immer unterwegs, STEHEN und BLEIBEN sind nur Grenzfaelle dieses ausgerichteten
 * Unterwegs.  Es gilt, einen WEG zur Aufhellung der ontologischen Fundamentalfrage zu suchen
 * und zu GEHEN.  OB er der EINZIGE oder UEBERHAUPT der rechte ist, das kann erst NACH dem
 * GANG entschieden werden."
 *
 * So there is a lot of TRAIL and ERROR in this module, but I think in the end Heidegger would
 * approve of our method (WEG) and our many possibilities (FREEDOM).
 */

/*
 * INCLUDES
*/
#include <snowpack/Laws_sn.h>

using namespace std;
using namespace mio;

/************************************************************
 * static section                                           *
 ************************************************************/

std::string SnLaws::current_variant = "DEFAULT"; //the default variant for SnLaws

/// @brief Switch on or off wind pumping in snow
const bool SnLaws::wind_pump = true;
/// @brief Switch on or off wind pumping in soil
const bool SnLaws::wind_pump_soil = true;

/// @brief Defines whether the extinction coefficient for SW radiation in snow is based on 20 or less bands
const bool SnLaws::band20 = false;
unsigned int SnLaws::swa_nBands = 5;
std::vector<double> SnLaws::swa_k;
std::vector<double> SnLaws::swa_pc;
std::vector<double> SnLaws::swa_fb;


/// @name SOIL PARAMETERS
//@{
/**
 * @brief Define Method and Coefficents for the computation of the influence of soil water
 * content on Evaporation from Bare Soil Layers:
 * - 1 ==> Resistance Approach, see Laws_sn.c
 * - 0 ==> Relative Humidity Approach, see Snowpack.cc
 * - -1 ==> none, assume saturation pressure and no extra resistance
 */
const bool SnLaws::soil_evaporation = true;

/// @brief Minimum soil surface resistance, 50 sm-1 (van den Hurk et al, 2000)
const double SnLaws::rsoilmin = 50.0;

/// @brief Minimum relative saturation of top soil layer
const double SnLaws::relsatmin = 0.05;

/// @brief Ratio of Porosity to Tortuosity for Soil
const double SnLaws::alpha_por_tor_soil = 0.05;

/// @brief Pore length for surface soil for ventilation (m)
const double SnLaws::pore_length_soil = 0.01;

/// @brief Field capacity of the soil (1). Above this levels, water begins to drain
const double SnLaws::field_capacity_soil = 0.15;
//@}


/**
 * @name Thermal conductivity
 * @brief Defines the constants and parameters for computing snow and soil thermal conductivity
 */
//@{
/// @brief Factor controlling ice to ice conduction
const double SnLaws::montana_c_fudge = 0.13;

/// @brief Factor controlling increase in water vapor transport and thus energy transport in wet snow
const double SnLaws::montana_vapor_fudge = 2.5;
//@}

/// @name Wind pumping (for soil see soil parameters)
//@{
/// @brief Used to decribe advective flux attenuation
const double SnLaws::wind_ext_coef = 0.1;

/// @brief Integral of snow density corresponding to the maximal displacement depth d_pump (m)
const double SnLaws::displacement_coef = 0.7;

/// @brief Ratio of Porosity to Tortuosity
const double SnLaws::alpha_por_tor = 0.07;
//@}

/// @brief Defines the smallest allowable viscosity (Pa s) that a viscosity law will return.
/// Value is DAMM SMALL -- smaller values than this are pretty unrealistic.
const double SnLaws::smallest_viscosity = 1.0e6;

/**
 * @name Snow viscosity
 * @brief Defines the constants and parameters for computing the snow viscosity \n
 * NOTE The Japanese variant uses a viscosity parametrization by Kojima
 * - t_term Temperature dependence
 * - visc_time_fudge: Empirical constant related to age of snow
 * - visc_ice_fudge : Empirical constant related to volumetric ice content; ori: 0.4
 * - visc_sp_fudge  : Empirical constant related to sphericity of snow grains
 */	
double SnLaws::visc_time_fudge = 11.;
double SnLaws::visc_ice_fudge = 0.5;
double SnLaws::visc_sp_fudge = 0.3;
double SnLaws::visc_water_fudge = 31.;
bool SnLaws::setfix = false;

double SnLaws::albedoCage = -0.000575;

SnLaws::AlbedoModel SnLaws::currentAlbedoModel = SnLaws::alb_lehning_2;
SnLaws::ViscosityCalVersion SnLaws::visc_cal = SnLaws::visc_cal_default;
SnLaws::TempDependence SnLaws::t_term = SnLaws::default_term;

/************************************************************
 * non-static section                                       *
 ************************************************************/

const bool SnLaws::__init = SnLaws::setStaticData("DEFAULT");
bool SnLaws::setStaticData(const std::string& variant)
{
	current_variant = variant;

	if (current_variant == "ANTARCTICA") {
		t_term = SnLaws::arrhenius_critical;
		visc_time_fudge = 0.;
		visc_ice_fudge = 0.51;
		visc_sp_fudge = 1.19;
		visc_water_fudge = 0.;
		visc_cal = SnLaws::visc_cal_ant;
		setfix = false;
	} else if (current_variant == "CALIBRATION") {
		t_term = SnLaws::arrhenius_critical;
		visc_time_fudge = 0.;
		visc_ice_fudge = 0.51;
		visc_sp_fudge = 1.19; //see factors in Laws_sn.cc
		visc_water_fudge = 31.;
		visc_cal = SnLaws::visc_cal_2010;
		setfix = false;
	} else {
		t_term = SnLaws::default_term;
		visc_time_fudge = 11.;
		visc_ice_fudge = 0.5;
		visc_sp_fudge = 0.3;
		visc_water_fudge = 31.;
		visc_cal = SnLaws::visc_cal_default;
		setfix = false;
	}

	if (current_variant == "JAPAN") {
		currentAlbedoModel = alb_nied;
	} else {
		currentAlbedoModel = alb_lehning_2;
	}

	if (current_variant == "ANTARCTICA") {
		SnLaws::albedoCage = 0.;
	} else {
		SnLaws::albedoCage = -0.000575;
	}

	swa_nBands = 5; // use 5 bands per default
	if (SnLaws::band20){
		swa_nBands = 20;
		double k_init[20]  = {0.0304, 0.0438, 0.0752, 0.1799, 0.3817, 0.5398, 0.7229, 1.618, 3.0372, 6.1273,
							 8.0446, 9.9419, 17.0339, 23.3685, 16.5658, 21.7658, 29.504, 66.922, 130.125,
							 120.405};
		double pc_init[20] = {1.98, 5.46, 8.9, 16.8, 6.54, 4.48, 4.4, 8.6, 10.81, 3.48, 1.2, 1.09, 6.89,
							 6.31, 1.38, 0.43, 1.27, 4.95, 1.93, 4.65};
		double fb_init[20];
		for (unsigned int ii=0; ii<swa_nBands; ii++) fb_init[ii] = 10.;//fudge_bohren

		swa_k.resize(swa_nBands);
		swa_pc.resize(swa_nBands);
		swa_fb.resize(swa_nBands);
		std::copy(k_init, k_init+swa_nBands, swa_k.begin());
		std::copy(pc_init, pc_init+swa_nBands, swa_pc.begin());
		std::copy(fb_init, fb_init+swa_nBands, swa_fb.begin());
	} else {
		double k_init[5]  = {0.059, 0.180, 0.525, 4.75, 85.23};
		double pc_init[5] = {16.3, 16.8, 15.4, 31.0, 20.5};
		double fb_init[5] = {29., 15., 5., 9., 35.}; // fudge_bohren
		swa_k.resize(swa_nBands);
		swa_pc.resize(swa_nBands);
		swa_fb.resize(swa_nBands);
		std::copy(k_init, k_init+swa_nBands, swa_k.begin());
		std::copy(pc_init, pc_init+swa_nBands, swa_pc.begin());
		std::copy(fb_init, fb_init+swa_nBands, swa_fb.begin());
	}

	return true;
}

/**
 * @name SNOW ALBEDO (defined as an absolute value--and not as a rate of change)
 * @brief Statistical models of surface snow albedo based on measurements (SWin and SWout, K&Z CM21)
 * from the Weissfluhjoch study plot.
 * @author Michael Lehning
 * @version 9Y.mm
 * @param variant The variant to use (DEFAULT, JAPAN, ...)
 * @param Edata
 * @param Tss Snow surface temperature (K)
 * @param Mdata
 * @param tday Age of surface snow (d)
 * @return Snow surface albedo (1)
 */
double SnLaws::compAlbedo(const std::string& variant, const ElementData& Edata, const double& Tss,
                          const SN_MET_DATA& Mdata, const double& age)
{
	if (variant != SnLaws::current_variant)
		setStaticData(variant);

	double Alb, Alb1;
	const double Ta = Mdata.ta;

	switch (currentAlbedoModel) {
		case alb_lehning_0: {
			double sqrt_age, lwi;
			const double weight=0.1;
			const double a = -0.2, b = 1.3, c = -0.012, d = -0.011, e = 0.0024, f = 0.018;
			const double g = -7.8e-6, h = -3.1e-3, i = 3.5e-4, j = 1.6e-7, k = -2.4e-7;
			const double l = -9.2e-6, m = 6.7e-5, n = -7.2e-5, o = 2.0e-4;

			if ( age < 0.001 ) {
				sqrt_age=0.000001;
			} else {
				sqrt_age=sqrt(age);
			}
			lwi = Constants::stefan_boltzmann*Mdata.ea*Ta*Ta*Ta*Ta;
			Alb1 = exp(a + b*sqrt_age + c*Tss + d*Mdata.vw + e*Mdata.rswr + f*Mdata.rh
					+ g*sqrt_age*Ta*Ta + h*sqrt_age*Tss + i*sqrt_age*lwi
					+ j*Ta*Ta*Tss + k*Ta*Ta*lwi + l*Tss*Mdata.rswr
					+ m*Tss*lwi + n*Tss*Mdata.rh + o*Mdata.vw*Mdata.rh);
			Alb = weight * Edata.dd * Snowpack::new_snow_albedo + (1 - weight * Edata.dd) * Alb1;
			break;
		}
		case alb_lehning_1: {
			double mf = 0.;
			const double av = 0.77;
			const double Cta = -0.0052, Cv = 0.0056, Clwc = -3.0, Crho = -0.0003, Cmf = -0.032;
			const double Crb = 0.06, Cdd = 0.017, Csp = 0.021, Ctss = 0.0084, Cswout = -6.8e-5;
			const double Cta_tss = -1.1e-5;

			if ( Edata.mk%100 > 19 ) {
				mf = 1.;
				// av *= exp(-age/1700.);
			}
			Alb1 = Crho*Edata.Rho + Clwc*Edata.theta[WATER] + Cdd*Edata.dd + Csp*Edata.sp + Cmf*mf 
				  + Crb*Edata.rb +  Cta*Ta + Ctss*Tss + Cv*Mdata.vw + Cswout*Mdata.rswr + Cta_tss*Ta*Tss;
			Alb = av + log(1.0 + Alb1);

			break;
		}
		case alb_lehning_2: {
			//TODO: this perfoms very badly (if not completly wrong) for wet snowpack
			//for example, February 2007 in Davos that saw very warm weather leads to albedos of 0.3 ...
			const double av = 0.8042; // Value of original regression: av=0.8042
			const double inter = 1.442, Cta = -0.006, Cv = 0.00762, Clwc = -0.2735;
			const double Crho = -0.000056, Crh = 0.0333, Crb = -0.301, Crg = 0.175;
			const double Cdd = 0.064, Csp = -0.0736, Ctss = 0.00459, Cswout = -0.000101;

			Alb1 = inter + Crho*Edata.Rho + Clwc*Edata.theta[WATER] + Cdd*Edata.dd + Csp*Edata.sp 
				  + Crg*Edata.rg + Crb*Edata.rb + Cta*Ta + Ctss*Tss + Cv*Mdata.vw + Cswout*Mdata.rswr 
				  + Crh*Mdata.rh + SnLaws::albedoCage*age;
			if ( Alb1 > 0. ) {
				Alb = MAX(Constants::min_albedo, MIN(Constants::max_albedo, av + log(Alb1)));
			} else {
				Alb = Constants::min_albedo;
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date.getJulianDate(), "Alb1=%lf set Alb to %lf", Alb1, Alb);
			}
			break;
		}
		case alb_nied: { //NIED (H. Hirashima)
			const double av=0.75;
			const double inter=1.005, Cta=-0.000249*2.0, Cv=0.00578, Clwc=-2.15;
			const double Crho=-0.000047, Cage=-0.00016*10.0,Crh=0.129, Crb=-0.306, Crg=0.107;
			const double Cdd=0.076, Csp=0.00964, Ctss=-0.000166, Cswout=-1.8e-5;

			Alb1 = inter + Crho*Edata.Rho + Clwc*Edata.theta[WATER] + Cdd*Edata.dd + Csp*Edata.sp 
				+ Crg*Edata.rg + Crb*Edata.rb + Cta*Ta + Ctss*Tss + Cv*Mdata.vw 
				+ Cswout*Mdata.rswr + Crh*Mdata.rh + Cage*age;

			if ( Alb1 > 0 ) {
				Alb = MAX(Constants::min_albedo, MIN(Constants::max_albedo, av + log(Alb1)));
			} else {
				Alb = Constants::min_albedo;
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date.getJulianDate(), "Alb1=%lf set Alb to %lf", Alb1, Alb);
			}
			break;
		}
		default: {
			prn_msg ( __FILE__, __LINE__, "err", -1., "ALB=%d not implemented yet!", currentAlbedoModel);
			throw IOException("The albedo model configured not implemented yet!", AT);
			break;
		}
	}

	return(Alb);
}

/**
 * @brief Helens Solution to Radiation Transfer \n
 * NOTE on fudge_bohren (fb): Larger values increase extinction --> Energy stays on top;
 * originally not band dependent, set to 10.0 for Neumann and to 5.0 for Dirichlet BC
 * @author Helene le Vesconte
 * @version Y.mm
 * @param I0 net shortwave radiation (W m-2)
 * @param useSnowLayers
 * @param multistream
 * @param Xdata
 */
void SnLaws::compShortWaveAbsorption(const double& I0, const bool& useSnowLayers, const bool& multistream, SnowStation& Xdata)
{
	int nE, e, bottom_element;
	double I0_band, dI, Ks;

	ElementData *EMS; // Avoids dereferencing the pointer

	EMS = &Xdata.Edata[0];
	nE = Xdata.getNumberOfElements();
	if ( nE == 0 ) {
		return;
	}

	if ( useSnowLayers ){
		bottom_element = Xdata.SoilNode - 1;
	} else {
		bottom_element = 0;
	}

	for (e = 0; e < nE; e++) {
		EMS[e].sw_abs = 0.;
	}

	// Compute absorbed radiation
	if ( multistream ) {
		for (unsigned int i = 0; i < swa_nBands; i++) {
			I0_band = I0 * swa_pc[i] / 100.;
			for (e = nE-1; e > bottom_element; e--) {
				Ks = swa_fb[i] * 0.84 * sqrt( 1000. * swa_k[i] / (2. * EMS[e].rg) ) * EMS[e].Rho / Constants::density_ice;
				dI = I0_band * (1. - exp(-(Ks * (EMS[e].L))));  // Radiation absorbed by element e in wavelength band i
				EMS[e].sw_abs += dI;
				I0_band -= dI;
			}
			EMS[bottom_element].sw_abs += I0_band;
		}
	} else { // ad hoc "1-Band" model
		I0_band = I0;
		for (e = nE-1; e > bottom_element; e--) {
			dI = I0 * (1. - exp(-(EMS[e].extinction() * (EMS[e].L))));  // Radiation absorbed by element e
			EMS[e].sw_abs += dI;
			I0_band -= dI;
		}
		EMS[bottom_element].sw_abs += I0_band;
	}

	for (e = nE-1; e >= bottom_element; e--) {
		if ( EMS[e].sw_abs < 0. ) {
			prn_msg(__FILE__, __LINE__, "err", -1., "NEGATIVE Shortwave Radiation %lf absorbed by element %d (nE=%d)", 
				   EMS[e].sw_abs, e, nE);
			throw IOException("SnLaws::compShortWaveAbsorption did not complete successfully", AT);
		}
	}
}

/**
 * @brief Computes the displacement depth in case of ventilation
 * @author Michael Lehning \n Charles Fierz
 * @version 10.01
 * @param Xdata
 * @return Displacement depth (m)
 */
double SnLaws::compWindPumpingDisplacement(const SnowStation& Xdata)
{
	int e = Xdata.getNumberOfElements();
	double d_pump=0.;      // Displacement depth (m)
	double hw = 0.0;       // Cumulated water equivalent (kg m-2)

	while ( e > 0 && hw < SnLaws::displacement_coef) {
		e--;
		if ( (hw + Xdata.Edata[e].L * Xdata.Edata[e].Rho) < SnLaws::displacement_coef ) {
			d_pump += Xdata.Edata[e].L;
		} else {
			d_pump += (SnLaws::displacement_coef - hw) / Xdata.Edata[e].Rho;
		}
		hw += Xdata.Edata[e].L * Xdata.Edata[e].Rho;
	}
	return (d_pump);
}

/**
 * @brief Computes the wind pumping velocity at the surface
 * @param *Mdata
 * @author Michael Lehning \n Charles Fierz
 * @version 10.01
 * @param d_pump Displacement depth (m)
 * @return Wind pumping velocity (m s-1)
 */
double SnLaws::compWindPumpingVelocity(const SN_MET_DATA& Mdata, const double& d_pump)
{
	//TODO Limit v_pump?
	return (Mdata.ustar / 0.4 * log((Mdata.z0 + d_pump) / Mdata.z0));
}

/**
 * @brief Computes the wind gradient near the snow-cover surface
 * @author Michael Lehning \n Charles Fierz
 * @param Edata
 * @param v_pump Wind velocity at element depth (m s-1)
 * @return Wind pumping velocity gradient (s-1)
 */
double SnLaws::compWindGradientSnow(const ElementData& Edata, double& v_pump)
{
	double dv, v_EXt;

	v_EXt = SnLaws::wind_ext_coef * (Edata.Rho + 2.e4 * Edata.theta[WATER]);
	dv = v_pump * (1. - exp(-v_EXt * (Edata.L)));
	v_pump -= dv;

	return (dv / (Edata.L));
}

/**
 * @brief Heat conduction in soil
 * @author Michael Lehning
 * @version 9Y.mm
 * @param *Edata const ElementData
 * @param dvdz Wind velocity gradient (s-1)
 * @return Soil thermal conductivity (W K-1 m-1)
 */
double SnLaws::compSoilThermalConductivity(const ElementData& Edata, const double& dvdz)
{
	double C_eff_soil, C_eff_soil_max, weight;
	const double c_clay = 1.3, c_sand = 0.27;
	const double alpha1 = 0.389, alpha2 = 0.3567, alpha3 = 61.61;
	const double beta1 = 6., beta2 = 4.978, c_mineral = 2.9;

	/*
	 * This ugly line does completely ignore the pledge for
	 * esthetics, which modelers should always adhere to if
	 * they want to be happy. Michi is unhappy and Ingo claims
	 * that it is not his fault..... but Martina's and David's
	 * 0: means no soil.
	 * 10000: means rock, which is also no soil but Ingo seems not to understand this.
	*/
	if ( (Edata.rg > 0.) && (Edata.rg < 10000.) ) {
		C_eff_soil_max = Edata.theta[SOIL] * c_mineral + Edata.theta[WATER]
                           * Constants::conductivity_water + Edata.theta[ICE] * Constants::conductivity_ice;

		/*
		 * This nice formulation is based on some tedious curve fitting by
		 * Martina, who became more happy when she saw how nicely it could
		 * now be implemented. The data for the frozen part stems from Kersten
		 * in "Geotechnical Engeneering for Cold Regions" article by Harlan
		 * and Nixon, while the water influence was deduced from deVries and
		 * Afgan in "Heat and Mass Transfer in the Biosphere".
		*/
		weight = (c_clay - Edata.soil[SOIL_K]) / (c_clay - c_sand);
		C_eff_soil = (beta1 + weight * beta2) * Edata.theta[ICE];
		if (Edata.theta[WATER] > 0.0001) {
			C_eff_soil += MAX (0.27,(alpha1 + alpha2 * weight) * log(alpha3 * Edata.theta[WATER]));
		} else {
			C_eff_soil += 0.27;
		}
		C_eff_soil = MIN (C_eff_soil_max, C_eff_soil);
	} else {
		C_eff_soil = Edata.soil[SOIL_K] + Edata.theta[WATER] * Constants::conductivity_water 
                       + Edata.theta[ICE] * Constants::conductivity_ice;
	}

	// Now check for possible ERRORS
	if (!(C_eff_soil >= 0.1*Edata.soil[SOIL_K] && C_eff_soil < 10.))
		prn_msg ( __FILE__, __LINE__, "wrn", -1., "Thermal Conductivity of Soil: %lf", C_eff_soil );

	// In case of dry soil, simply use the given conductivity with a possible ventilation part
	if (SnLaws::wind_pump_soil){
		C_eff_soil += SnLaws::alpha_por_tor_soil * Constants::specific_heat_air 
                        * Edata.soil[SOIL_RHO] * SnLaws::pore_length_soil * SnLaws::pore_length_soil * dvdz;
	}

	return(C_eff_soil);
}

/**
 * @brief Computes the enhancement factor for water vapor transport in wet snow
 * @author Michael Lehning
 * @version 9Y.mm
 * @param *Xdata
 * @param e Element number
 * @return Enhancement factor
 */
double SnLaws::compEnhanceWaterVaporTransportSnow(const SnowStation& Xdata, const int& e)
{
	int e1 = e;
	double vapor_enhance = 1.;

	while ((e1 >= MAX(Xdata.SoilNode, e-7)) && (vapor_enhance < SnLaws::montana_vapor_fudge) ) {
		//TODO check limit on theta_water
		if ( (Xdata.Edata[e1].theta[WATER] > 0.001) && (Xdata.Edata[e].theta[ICE] < 0.7) ) {
			vapor_enhance = SnLaws::montana_vapor_fudge;
		}
		e1--;
	}

	if (vapor_enhance > 1.0) 
		vapor_enhance *= ((8. - (e - e1)) / 7.);
	
	return (vapor_enhance);
}

/**
 * @brief Heat conduction in snow \n
 * Actual version: k = C1*[C2 + C3 + C4 + C5] \n
 * Adams/Sato model. The following piece of code was programmed by Perry on a warm
 * June afternoon in Birmensdorf.  H.U. Gubler was viciously attacking Ammann, saying nasty
 * things about our beloved chef, (even if his hair is a bit funny), in the DAVOSER ZEITUNG.
 * Perry was dreaming about San Diego and .... and Michael was dreamy about rock climbing in
 * Italy with funky, spunky Italian girls, Borolo Bob was busy with his wife on his way to
 * Mallorca -- to watch sexy TOPLESS bikini-slip clad GERMAN secretaries on the BEACH.
 * As usual, if this code does not work it is his FAULT ....  (Hey, Bob, when are we going to
 * meet ED ADAMS?  And remember PLEDGE 7 of Metamorphism.c.  The theory behind this piece of
 * code can be found in Bob's ROUGH DRAFT on MICROSTRUCTURE. \n
 * Much much later, Michael re-programmed this code and introduced the effect of liquid water.
 * He extended the Adams/Sato model and hoped that he would not loose contact with Keegan, even
 * he was not going to fulfill the explicit wishes of his Grandma. \n
 * The model was now being used from Finland to the US. We will also conquer the southern hemis-
 * sphere.
 * @author Michael Lehning
 * @version 8.10 (??)
 * @param *Edata
 * @param dvdz Wind velocity gradient (s-1)
 * @return Thermal conductivity of snow (W K-1 m-1)
 */
double SnLaws::compSnowThermalConductivity(const ElementData& Edata, const double& dvdz)
{
	double C_eff, C1, C2, C3, C4, C5;
	double kap;                 // The conductivity including latent heat transfer (W K-1 m-1)
	double Ap, Aip, Aiw;        // Cross sectional areas of conduction paths (m2)
	double rg, rb;              // Grain and bond radius (m)
	double Te;                  // Element temperature (K)
	double Lh = Constants::lh_sublimation;
	double P = 950.;

	rg = MM_TO_M(Edata.rg);
	rb = MM_TO_M(Edata.rb);
	Te = MIN(Edata.Te, Constants::melting_tk);

	// Check for elements with no ice and assume they contain only water
	if ( Edata.theta[ICE] < Constants::min_ice_content ) {
		return(Constants::conductivity_water);
	}

	// Important are the conductivities and cross sectional areas.
	kap = Constants::conductivity_air + ((Lh*Lh* Constants::diffusion_coefficient_in_air * P 
		 * lw_SaturationPressure(Te)) 
		 / (Constants::gas_constant*Constants::gas_constant * Te*Te*Te * (P - lw_SaturationPressure(Te))));

	/*
	* Determine C1 = nca/ncl (nca:=number of grains per unit area; ncl:=number
	* of grains per unit length) = ncl since nca=ncl^2
	* Now we can also determine the series cross sectional areas
	*/
	C1 =  (3.0 * Edata.Rho / Constants::density_ice) / ( 4.0 * Constants::pi * rg*rg*rg);
	C1 =  pow(C1, 0.3333333);

	/*
	* Determine C2 = a funny fraction containing PI and the coordination number
	* (I believe that this is the solid ice conduction from one grain through a
	* bond/neck into another ice grain.)
	*/
	C2 = Constants::pi * Constants::pi * rb * Constants::conductivity_ice * Edata.N3 / 32.0;

	// Compute cross-sectional areas
	Ap = Metamorphism::csPoreArea(Edata); // (mm2)
	Aiw = MAX (0., Edata.theta[WATER] * (1. / C1 - rg) / C1 * (Ap + Constants::pi * rg*rg));
	Aip = MAX (0., Constants::pi * (rg*rg - rb*rb) - Aiw);

	/*
	 * Determine C3 = a funny fraction still containing PI but a lot of very funny
	 * k coefficients. (I believe that this is the series conduction in which
	 * heat goes through a grain, then across a pore and into another grain.)
	*/
	C3 = Constants::conductivity_ice * kap * Aip;
	C3 = C3 / (rg * kap + (1. / C1 - rg) * Constants::conductivity_ice);

	/*
	 * Determine C4 = a funny fraction containing one PI and  a few (one) funny
	 * k coefficients. (I believe that this is the series conduction purely
	 * across the pore, the most believable part.)  Note that the fraction is
	 * divided by 1/C1, so I multiply the whole thing ...
	*/
	C4 = C1 * kap * Ap;

	/*
	 * Determine C5 = a funny fraction still containing PI but a lot of very funny
	 * k coefficients. (I believe that this is the series conduction in which
	 * heat goes through a grain, then across water and into another grain.)
	*/
	C5 = Constants::conductivity_ice * Constants::conductivity_water * Aiw;
	C5 = C5 / (rg * Constants::conductivity_water  + (1./C1 - rg) * Constants::conductivity_ice);

	C_eff  = SnLaws::montana_c_fudge * C1 * (C2 + C3 + C4 + C5) * (2.0 - Edata.dd) * (1.0 + pow(Edata.theta[ICE], 1.7)) * (0.5 + Te*Te / (Constants::melting_tk*Constants::melting_tk));

	if ( !((C_eff < 5.*Constants::conductivity_ice) && (C_eff > 0.2*Constants::conductivity_air)) && !ALPINE3D ) {
		prn_msg(__FILE__, __LINE__, "wrn", -1., "Conductivity out of range (0.2*Constants::conductivity_air=%.3lf, 5.*Constants::conductivity_ice=%.3lf):", 0.2 * Constants::conductivity_air, 5. * Constants::conductivity_ice);
		prn_msg(__FILE__, __LINE__, "msg-", -1., "C_eff: %lf  C_1: %lf  C_2: %lf  C_3: %lf  C_4: %lf  C_5: %lf", C_eff, C1, C2, C3, C4, C5);
		prn_msg(__FILE__, __LINE__, "msg-", -1., "C_Air: %lf  C_Water: %lf  C_Ice: %lf", Constants::conductivity_air, Constants::conductivity_water, Constants::conductivity_ice);
		prn_msg(__FILE__, __LINE__, "msg-", -1., "Ap: %lf  kap: %lf  Thet_i: %lf  Thet_w: %lf  T: %lf", Ap, kap, Edata.theta[ICE], Edata.theta[WATER], Te );
		prn_msg(__FILE__, __LINE__, "msg-", -1., "type: %3d  rg: %.3lf mm  rb: %.3lf mm N3: %lf  Lel: %.3lf mm", Edata.type, Edata.rg, Edata.rb, Edata.N3, M_TO_MM(Edata.L));
	}

	if ( !(C_eff < Constants::conductivity_ice) ) {
		C_eff = Constants::conductivity_ice;
	}
	if ( !(C_eff > Constants::conductivity_air) ) {
		C_eff = Constants::conductivity_air;
	}
	// Now introduce the effect of windpumping
	if ( SnLaws::wind_pump /* && Edata.theta[WATER] < 0.001 */ ) {
		C_eff += SnLaws::alpha_por_tor * Constants::specific_heat_air * Edata.Rho  * dvdz / (C1*C1);
	}
	return (C_eff);
}

/**
 * @brief SENSIBLE HEAT EXCHANGE COEFFICIENT (Surface Energy Exchange)
 * @author Michael Lehning
 * @version 9Y.mm
 * @param Mdata
 * @param Xdata
 * @return Exchange coefficient for sensible heat (1)
*/
double SnLaws::compSensibleHeatCoefficient(const SN_MET_DATA& Mdata, const SnowStation& Xdata, const double& height_of_meteo_values)
{
	double c;
	double z, karman = 0.4, lrat;

	z = MAX (0.5, height_of_meteo_values - Xdata.cH + Xdata.Ground);
	if ( (Xdata.cH - Xdata.Ground) > 0.03 ) {
		//assert(Mdata.z0>0.);
		lrat = log(z / Mdata.z0);
	} else {
		//assert(Xdata.BareSoil_z0>0.);
		lrat = log(z / Xdata.BareSoil_z0);
	}

	c = karman * Mdata.ustar / (0.74 * MAX (0.7, lrat-Mdata.psi_s));

	return(c * Constants::density_air * Constants::specific_heat_air);
	//return(2.77e-3*Mdata.vw*DENSITY_AIR*Constants::specific_heat_air);
}

/**
 * @brief Latent heat flux including consideration of soil (one active element)
 * This method uses the Relative Humidity approach:
 * ql = beta*(eA - eS) Latent heat transfer. eA and eS are the vapor
 * pressures of air and snow, respectively.
 * @author Michael Lehning
 * @version 9Y.mm
 * @param Mdata
 * @param Xdata
 * @return Latent heat flux (W m-2)
 */
double SnLaws::compLatentHeat_Rh(const SN_MET_DATA& Mdata, SnowStation& Xdata, const double& height_of_meteo_values)
{
	double beta, eA, eS;
	double Vp1, Vp2;
	double th_w_ss;

	const double T_air = Mdata.ta;
	const double Tss = Xdata.Ndata[Xdata.getNumberOfElements()].T;

	// Vapor Pressures
	if ( Xdata.getNumberOfElements() == 0 ) {
		th_w_ss = 0.0;
	} else {
		th_w_ss = Xdata.Edata[Xdata.getNumberOfElements()-1].theta[WATER];
	}

	// TODO The part below needs to be rewritten in a more consistent way !!!
	//      In particular, look closely at the condition within compLatentHeat()
	eA = Mdata.rh * lw_SaturationPressure(T_air);
	Vp1 = lw_SaturationPressure(Tss);
	Vp2 = lw_SaturationPressure(Tss);

	// First, the case of no snow
	if ( Xdata.getNumberOfNodes() == Xdata.SoilNode + 1 && Xdata.getNumberOfElements() > 0 ) {
		if( Tss < Constants::melting_tk ) {
			eS = Vp1 ;
		} else {
			/*
			 * Soil evaporation can now be computed using the Relative Humidity approach below,
			 * or a Resistance approach modifying the ql value instead of the eS. The latter
			 * function is defined in compLatentHeat, and the Switch SnLaws::soil_evaporation is found
			 * in Laws_sn.h
			*/
			if (SnLaws::soil_evaporation && th_w_ss < Xdata.Edata[Xdata.SoilNode-1].soilFieldCapacity()){
				eS = Vp2 * 0.5 * ( 1 - cos (MIN (Constants::pi, th_w_ss * Constants::pi /
						(Xdata.Edata[Xdata.SoilNode-1].soilFieldCapacity() * 1.6))));
			} else {
				eS = Vp2;
			}
		}
	} else {
		// for snow assume saturation
		if ( Tss < Constants::melting_tk ) {
			eS = Vp1;
		} else {
			eS = Vp2;
		}
	}
	// Now the latent heat
	beta = SnLaws::compLatentHeat(Mdata, Xdata, height_of_meteo_values);

	return (beta * (eA - eS));
}

/**
 * @brief LATENT HEAT EXCHANGE (Surface Energy Exchange)
 * @author Michael Lehning
 * @version 9Y.mm
 * @param Mdata
 * @param Xdata
 * @return Latent heat flux (W m-2)
 */
double SnLaws::compLatentHeat(const SN_MET_DATA& Mdata, SnowStation& Xdata, const double& height_of_meteo_values)
{
	double c, eS, eA;
	double z, karman = 0.4, lrat;

	z = MAX (0.5, height_of_meteo_values - Xdata.cH + Xdata.Ground);
	if ( (Xdata.cH - Xdata.Ground) > 0.03 ) {
		//assert(Mdata.z0>0);
		lrat = log(z/Mdata.z0);
	} else {
		//assert(Xdata.BareSoil_z0>0);
		lrat = log(z / Xdata.BareSoil_z0);
	}
	c = karman * Mdata.ustar / (0.74 * MAX (0.7, lrat-Mdata.psi_s));

	/*
	 * Below, David Gustafsson (davidg@kth.se) has introduced a resistance approach for latent
	 * heat exchange from bare soil as alternative to the relative hyumidity (RH) approach
	 * implemented in Snowpack.c (line 473): \n
	 * An additional resistance, dependent on the relative saturation of the top soil layer,
	 * is used to reduce the heat exchange coefficient in the case of evaporation:
	 * c = 1/(Ra + Rsoil), where Ra = 1/c as computed above, and
	 * Rsoil = 50 [s/m] * field_capacity_soil / theta_soil. \n
	 * A new switch SnLaws::soil_evaporation is defined in Constants.h to select method.
	 * The resistance formulation originates from van den Hurk et al.(2000) "Offline validation
	 * of the ERA40 surface scheme": ECMWF Tech.Memo 295. \n
	 * A difference from the RH method is that the surface vapour pressure is always assumed
	 * to be at saturation at the surface. Thus, some unrealistic effects of the RH method in
	 * present form are avoided -> the RH approach tend to estimate negative vapour gradient
	 * between the surface and the atmosphere, causing large condensation, even if the top soil
	 * layer is saturated, and even if the soil surface is warmer than the atmosphere! If a RH
	 * method should work in a discretized model, it is important to consider the difference
	 * between vapour pressure at the surface and the average of the top soil layer. \n
	 * The soil resistance is only used for bare soil layers, when TSS >= 0C and eSurf >= eAtm
	*/
	if (Xdata.getNumberOfNodes() == Xdata.SoilNode + 1 && Xdata.getNumberOfElements() > 0 &&
		Xdata.Ndata[Xdata.getNumberOfElements()].T >= Constants::melting_tk && SnLaws::soil_evaporation == 1) {
		eA = Mdata.rh * lw_SaturationPressure(Mdata.ta);
		eS = lw_SaturationPressure(Xdata.Ndata[Xdata.getNumberOfElements()].T);
		if (eS >= eA) {
			c = 1. / c + SnLaws::rsoilmin / (MAX (SnLaws::relsatmin, MIN(1. ,Xdata.Edata[Xdata.getNumberOfElements()-1].theta[WATER]
                     / Xdata.Edata[Xdata.SoilNode-1].soilFieldCapacity())));
			c = 1. / c;
		}
	}
	return (c * 0.622 * Constants::lh_sublimation / Constants::gas_constant_air / Mdata.ta);
}

/**
 * @brief LONGWAVE RADIATION COEFFICIENT  (This routine might look a bit unusual: \n
 * Radiation is treated as a CONVECTIVE boundary condition, similar to the sensible and latent heat
 * exchanges.  The exchange coefficient, however, is not a constant, dependent on say the wind
 * velocity, rather it is dependent on the temperature.  This routine computes the "pseudo"
 * convective heat exchange coefficient for radiation.)
 * @author Perry Bartelt
 * @version 9Y.mm
 * @param T_snow Snow surface temperature (K)
 * @param T_atm Temperature of the atmosphere, i.e., air (K)
 * @param e_atm Emissivity of the atmosphere (1)
 * @return LW radiation coefficient (?)
 */
double SnLaws::compLWRadCoefficient(const double& t_snow, const double& t_atm, const double& e_atm)
{
	return (Constants::stefan_boltzmann 
		   * ((t_snow * t_snow) + (sqrt(e_atm) * t_atm * t_atm)) 
		   * ((t_snow) + (pow(e_atm, 0.25) * t_atm)));
}

/**
 * @brief NEW SNOW VISCOSITY (dendritic snow, i.e., dd > 0.) \n
 * Actual version : ml_lw_VS_Lehning from r7.7 \n
 * This is Michael's viscosity routine, which is not a function of micro-structure, but which
 * is nonetheless pretty important because it is numerically STABLE and does predict decent
 * settling rates, sometimes a bit too high. This routine was (is) used for NEW or WET snow.
 * @author Michael Lehning
 * @version 9Y.mm
 * @return Viscosity of new snow (Pa s)
 */
double SnLaws::NewSnowViscosityLehning(const ElementData& Edata)
{
	double rho_hn = MAX(Constants::min_hn_density, Edata.Rho);
	
	if (rho_hn > 913.) //upper boundary
		return (1.e9 * SnLaws::smallest_viscosity);

	if (Edata.theta[WATER] > 0.001) {
		//return (0.01*pow(rho_hn, 4.7));
		return(0.0001 * pow(rho_hn, 5.5));
	} else {
		return (0.007 * pow(rho_hn, (4.75 - K_TO_C(Edata.Te) / 40.)));
	}
}

/**
 * @brief Computes the additional stress during initial metamorphic settling process. Note the negative sign!
 * - TODO rewrite description!
 * - Once upon a time: \n
 *   meaning that the strain rate is purely a function of the overburden stress. Experiments, for
 *   example De Quervain's creep experiments, CLEARLY show that snow will settle in the ABSENCE
 *   of a pressure.  This effect is included in ALL the viscosity computations using an INITIAL
 *   STRESS.  This routine finds this stress which is superimposed on the overburden pressure.
 *   Helps predict settlement of NEW snow better.
 * - Now: \n
 *   Added stress sig0 due to surface tension is significant at LOW densities only and when
 *   applied stresses are very small. This formulation differs markedly from the original formulation
 *   which was proposed by Bob Brown and used up to research version r8.0 (see above).
 * @author Michael Lehning
 * @version 9Y.mm
 * @param i_viscosity_model Set the viscosity model to use ("VS_CALIBRATION", "DEFAULT" are relevant)
 * @param Edata
 * @param date current
 * @return Initial stress (Pa)
 */
double SnLaws::compInitialStress(const std::string& variant, ElementData& Edata, const mio::Date& date)
{
	if (variant == "CALIBRATION") {
		return initialStressCALIBRATION(Edata, date);
	} else {
		return initialStressDEFAULT(Edata, date);
	}
}

/// default r897: TODO name parameters?
double SnLaws::initialStressDEFAULT(ElementData& Edata, const mio::Date& date)
{
	//TODO CDot: uncomment
	//double sigReac=0., sigMetamo=0.;
	//const double age = MAX(0., date.getJulianDate() - Edata.date.getJulianDate());

	//sigReac = 15.9 * Edata.CDot * exp(-age / 101.);
	//if (Edata.dd > Constants::eps) {
	//	sigMetamo = 37.0e3 * Metamorphism::ddRate(Edata);
	//}
	//return sigReac + sigMetamo;

	// default r837
	//(void) date;
	const double sigTension = 0.11;  // Ice surface tension (N m-2)
	if ((Edata.dd < 0.9) && (Edata.dd > 0.3)) { // default r837, factor -3. from r712
		return (-3. * Metamorphism::ddRate(Edata) * sigTension / MM_TO_M(Edata.rg));
	} else {
		return 0.;
	}
}

double SnLaws::initialStressCALIBRATION(ElementData& Edata, const mio::Date& date)
{
	const unsigned int version = 0;
	const double sigTension = 0.11;  // Ice surface tension (N m-2)

	switch (version) {
	case 837: // as of revision 837
		if ((Edata.dd < 0.9) && (Edata.dd > 0.3)) {
			double facIS = -3.; // default r712
			if (SnLaws::visc_cal == SnLaws::visc_cal_steinkogler) facIS = -1.5; //-1.1; //-0.5; //
			return (facIS * Metamorphism::ddRate(Edata) * sigTension / MM_TO_M(Edata.rg)); //default
		} else {
			return 0.;
		}
		break;
	default: // new calibration, currently r987
		double sigReac=0., sigMetamo=0.;
		const double age = MAX(0., date.getJulianDate() - Edata.date.getJulianDate());
		//09-24: sig0 = 33.7 * Edata.CDot * exp(-age/553.); //3.5 //5. //2.7 //1.e3 * exp(-age/77.) * Edata.CDot / Edata.C;
		sigReac = 15.9 * Edata.CDot * exp(-age/101.); //tst2: 553. //tst1: 735. //
		//if ( 1 && (Edata.dd > 0.2) /*((Edata.dd < 0.9) && (Edata.dd > 0.3))*/ ) {
		if (true && (Edata.dd > Constants::eps) /*((Edata.dd < 0.9) && (Edata.dd > 0.3))*/) {
			//sig0 *= -10. * mm_ddRate(Edata);
			//sig0 += 73.5 * mm_ddRate(Edata) * sigTension / MM_TO_M(Edata.rg); // ori: -3.0; -23.
			//sigMetamo = 73.5 * mm_ddRate(Edata) * sigTension / MM_TO_M(Edata.rg); // ori: -3.0; -23.
			//sigMetamo = 139.7 * mm_ddRate(Edata) * sigTension / MM_TO_M(Edata.rg); // 2010-10-22
			//sigMetamo = 29.0 * mm_ddRate(Edata) / MM_TO_M(Edata.rg); // 2010-10-23
			sigMetamo = 37.0e3 * Metamorphism::ddRate(Edata); // 2010-10-23
			if (false && Edata.dd < 0.9 ) {
				//sigMetamo /= Edata.dd; // pow(Edata.dd, 1./2.);
				sigMetamo *= pow(Edata.dd, 1./2.);
			}
			Edata.EDot = sigReac;
			return (sigReac + sigMetamo);
		} else { //if ( age > 3000. ) {
			Edata.EDot = sigReac;
			return (sigReac);
		}
		break;
	}
}

/**
 * @brief Determines the fudge factor for viscosity \n
 * This fudge factor takes into account bond-ice imperfections and the effect of liquid water
 * @version 9.10
 * @param Edata
 * @param date current date
 * @return Fudge factor for snow viscosity
 */
double SnLaws::snowViscosityFudgeDEFAULT(const ElementData& Edata, const mio::Date& date)
{
	double visc_fudge=0., sp_fudge;
	double age = MAX(0., date.getJulianDate() - Edata.date.getJulianDate());
	const double thresh_rho1=0.5, thresh_rho2=0.7; // Thresholds for enhanced viscosity

	if ( Edata.mk%100 >= 20 && Edata.theta[WATER] < 0.005 ) {
		sp_fudge = 0.;
	} else {
		sp_fudge = SnLaws::visc_sp_fudge;
	}

	visc_fudge += SnLaws::visc_time_fudge * (1. - sqrt(MIN(1., age / 77.))) * (1. + MIN(0.3, (263.15 - Edata.Te) / 17.));
	visc_fudge += (SnLaws::visc_ice_fudge / Edata.theta[ICE]) + (sp_fudge * sqrt(Edata.sp))
		+ (3. * Edata.theta[WATER] / Edata.theta[ICE] * 0.5 * (1. - Edata.theta[ICE]));

	if ( Edata.theta[ICE] > thresh_rho1 ) {
		visc_fudge *= (1. - Edata.theta[ICE]);
	}
	if ( Edata.theta[ICE] > thresh_rho2 ) {
		visc_fudge *= (1. - Edata.theta[ICE]);
	}
	return visc_fudge;
}

/**
 * @brief Calibrate the fudge factor for viscosity \n
 * NOTE This is the fudge playground for calibrations used in SnowViscosityCALIBRATION()
 * @version 10.11
 * @param Edata
 * @param date current date
 * @return Fudge factor for snow viscosity
 */
double SnLaws::snowViscosityFudgeCALIBRATION(const ElementData& Edata, const mio::Date& date)
{
	double visc_fudge, sp_fudge;
	double age = MAX(0., date.getJulianDate() - Edata.date.getJulianDate());
	double thresh_rho1, thresh_rho2=0.7; // Thresholds for enhanced viscosity
	bool use_thresh=false;

	if ( Edata.mk%100 >= 20 && Edata.theta[WATER] < 0.005 ) {
		sp_fudge = 0.;
	} else {
		sp_fudge = SnLaws::visc_sp_fudge;
	}

	switch (SnLaws::visc_cal) {
		case SnLaws::visc_cal_new: {} // Calibration currently under test
		case SnLaws::visc_cal_2010: case SnLaws::visc_cal_ant: { // Calibration fall 2010 & adaptation to Antarctica TODO by Fierz
			double ice_fudge = SnLaws::visc_ice_fudge / Edata.theta[ICE];
			if ( SnLaws::visc_cal != SnLaws::visc_cal_ant ) {
				ice_fudge *= 1. - logisticFunction(Edata.theta[ICE], 0.431, 0.067) * Edata.theta[ICE];
			} else {  // TODO doit!
				ice_fudge *= 1. - logisticFunction(Edata.theta[ICE], 0.351, 0.097) * Edata.theta[ICE] * (1. - Edata.theta[ICE]);
			}
			sp_fudge *= (1. - logisticFunction(Edata.sp, 1.0, 0.077)) * sqrt(Edata.sp);
			visc_fudge = ice_fudge * (1. + sp_fudge) + 0.1;
			visc_fudge *= (1. + SnLaws::visc_water_fudge * Edata.theta[WATER] * (1. - Edata.theta[ICE]*Edata.theta[ICE]));
			break;
		}
		case SnLaws::visc_cal_steinkogler: { // Walter Steinkogler's playground; master thesis, September 2009
			visc_fudge += SnLaws::visc_time_fudge / exp(age / 35.);
			visc_fudge += (SnLaws::visc_ice_fudge / Edata.theta[ICE]) + (sp_fudge * sqrt(Edata.sp))
				+ (3. * Edata.theta[WATER] / Edata.theta[ICE] * 0.5 * (1. - Edata.theta[ICE]));
			thresh_rho1 = 0.5; // rho_dry > 458.5 kg m-3
			use_thresh=true;
			break;
		}
		default: { // calibrated version r837
			visc_fudge += SnLaws::visc_time_fudge * (1. - sqrt(MIN(1., age / 77.))) * (1. + MIN(0.3, (263.15 - Edata.Te) / 17.));
			visc_fudge += (SnLaws::visc_ice_fudge / Edata.theta[ICE]) + (sp_fudge * sqrt(Edata.sp))
				+ (3. * Edata.theta[WATER] / Edata.theta[ICE] * 0.5 * (1. - Edata.theta[ICE]));
			use_thresh=true;
			thresh_rho1 = 0.5; // rho_dry > 458.5 kg m-3
		}
	}

	if (use_thresh) {
		if (Edata.theta[ICE] > thresh_rho1) {
			visc_fudge *= (1. - Edata.theta[ICE]);
		}
		if (Edata.theta[ICE] > thresh_rho2) {
			visc_fudge *= (1. - Edata.theta[ICE]);
		}
	}
	return visc_fudge;
}

/**
 * @brief Computes the temperature term of viscosity
 * @author Charles Fierz
 * @version 9.10
 * @param Te Element temperature
 * @return Temperature term of snow viscosity
 */
double SnLaws::snowViscosityTemperatureTerm(const double& Te)
{
	const double Q = 67000.; // Activation energy for defects in ice J mol-1

	switch (SnLaws::t_term) {
		case SnLaws::arrhenius:
			return (1. / lw_ArrheniusLaw(Q, Te, 263.));
			break;
			case SnLaws::arrhenius_critical: {
				const double Q_fac=0.39; // Adjust Q to snow. from Sz: 0.24; Ml: 0.39; Fz: 0.31
				const double T_r=265.15; // Reference temperature (K); ori: 265.15 + 2.15
				return ((1. / lw_ArrheniusLaw(Q_fac * Q, Te, T_r)) * (0.3 * sqrt(Constants::melting_tk - Te) + 0.4));
				break;
			}
			case SnLaws::steinkogler: { // Master thesis, September 2009
				return (0.35 * pow((274.15 - Te), 0.5));
				break;
			}
			default: { // from r243
				return (9. - 8.7 * exp(0.015 * (Te - Constants::melting_tk)));
				break;
			}
	}
}

/**
 * @brief SNOW VISCOSITY (all types of snow)
 * TODO revise description
 * @param variant current
 * @param i_viscosity_model to use
 * @param Edata
 * @param date current
 */
double SnLaws::compSnowViscosity(const std::string& variant, const std::string& i_viscosity_model,
                                 ElementData& Edata, const mio::Date& date)
{
	if (variant != SnLaws::current_variant)
		setStaticData(variant);

	string viscosity_model(i_viscosity_model);
	IOUtils::toUpper(viscosity_model);

	if (viscosity_model == "VS_KOJIMA"){
		return snowViscosityKOJIMA(Edata, date);
	} else if (viscosity_model == "VS_CALIBRATION"){
		return snowViscosityCALIBRATION(Edata, date);
	} else if (viscosity_model == "DEFAULT"){
		return snowViscosityDEFAULT(Edata, date);
	} else {
		throw InvalidArgumentException("Unknown viscosity model: "+viscosity_model, AT);
	}
}

/**
 * @brief SNOW VISCOSITY (all types of snow)
 * TODO revise description
 * - bb_lw_VS_Montana from r7.7
 * - Fudge : according to Ml (8.1)
 * - Bob Brown's MICRO-STRUCTURE law Clearly the BEST law we have right now
 *   but also the most UNSTABLE:  note that the viscosity is not only a function of the grain
 *   dimensions, but also a function of the overburden stress.
 *   A series of equations that collectively give the viscosity Vis = S/eDot.
 *   The microstructure parameters rb and rg are obtained through Edata pointer and are in mm.
 *   The secondary microstructure parameters  L and rc are also computed in mm. This means that
 *   the dimensions of rg, rb, L, &  rc are in mm since they show up in the following  equations
 *   as ratios to give dimensionless numbers.
 * 	- LINEAR:
 *    This viscosity is not a function of stress and is therefore a linear viscosity.  Its value
 *    depends on rb, rg, N3, Edata.theta[ICE] and T. The expression  ((N3*Edata.theta[ICE])/(4.))*((rb*rb)/(rg*rg))
 *    determines the neck stress relative to the snow stress.
 * 	- NON-LINEAR:
 *    The units here are in m, s and Pa. Vis has dimensions of Pa s.
 *    Note that the viscosity is a function of the square of the applied stress.
 *    Density does not show up explicitly, but is buried implicitly in the equation
 *    because N3, L, rg and rb all change with time as the material deforms.
 * @author Michael Lehning \n Charles Fierz
 * @version 9.10
 * @param *Edata
 * @return Snow viscosity (Pa s)
 */
double SnLaws::snowViscosityDEFAULT(ElementData& Edata, const mio::Date& date)
{
	const double eps1Dot = 1.76e-7;    // Unit strain rate (at stress = 1 MPa) (s-1)
	const double sig1 = 0.5e6;         // Unit stress from Sinha's formulation (Pa)
	const double sigNeckYield = 0.4e6; // Yield stress for ice in neck (Pa)
	const double sig = -Edata.C;      // Overburden stress, that is, absolute value of Cauchy stress (Pa)
	double Te;                         // Element temperature (K)
	double visc_macro, visc_micro;           // Structure related multiplying factors
	double visc_fudge, visc_factor;          // Fit and reference parameters
	double eta;                        // Viscosity (Pa s)

	Te = MIN(Edata.Te, Constants::melting_tk);

	// Check needed while JAM set!
	if ( Edata.theta[WATER] > 0.3 ) {
		return (1.e9 * SnLaws::smallest_viscosity);
	}
	// Check that you are not in a ice or/and water layer
	if ( Edata.theta[ICE] + Edata.theta[WATER] > 0.99 ) {
		return (1.e9 * SnLaws::smallest_viscosity);
	}

	visc_fudge = SnLaws::snowViscosityFudgeDEFAULT(Edata, date);
	visc_factor = (sig1*sig1*sig1 / (eps1Dot * visc_fudge*visc_fudge*visc_fudge));
	visc_macro = Edata.neck2VolumetricStrain();
	visc_micro = Edata.neckStressEnhancement();
	eta = (1. / visc_macro) * SnLaws::snowViscosityTemperatureTerm(Te) * visc_factor;
	// NOT YIELDING, LINEAR; sigNeckYield = 0.4 MPa
	if ( (visc_micro * sig) <= 100. * sigNeckYield ) {
		eta /= visc_micro * sigNeckYield*sigNeckYield;
	// YIELDING, NON-LINEAR
	} else {
		eta /= visc_micro*visc_micro*visc_micro * sig*sig;
	}
	return eta;
}

/**
 * @brief SNOW VISCOSITY according to formulation by Kojima
 * @version 9.10
 * @param Edata
 * @param date current date
 * @return Snow viscosity (Pa s)
 */
double SnLaws::snowViscosityKOJIMA(const ElementData& Edata, const mio::Date& date)
{
	(void)(date);
	return (8.64e6 * exp(0.021*Edata.Rho));
}

/**
 * @brief Calibrate snow viscosity
 * NOTE This is the test or playground version for calibrating settling
 * @version 10.11
 * @param Edata
 * @param date current date
 * @return Snow viscosity (Pa s)
 */
double SnLaws::snowViscosityCALIBRATION(ElementData& Edata, const mio::Date& date)
{
	const double eps1Dot = 1.76e-7;    // Unit strain rate (at stress = 1 MPa) (s-1)
	const double sig1 = 0.5e6;         // Unit stress from Sinha's formulation (Pa)
	const double sigNeckYield = 0.4e6; // Yield stress for ice in neck (Pa)
	const double sig = -Edata.C;      // Overburden stress, that is, absolute value of Cauchy stress (Pa)
	double Te;                         // Element temperature (K)
	double visc_macro, visc_micro;           // Structure related multiplying factors
	double visc_fudge, visc_factor;          // Fit and reference parameters
	double eta;                        // Viscosity (Pa s)

	Te = MIN(Edata.Te, Constants::melting_tk);

	// TODO Check whether the two commented checks below are needed!
	// If the element length is SMALLER than the grain size then the thing aint settling ....
	//  if( Edata.L <= 2.*MM_TO_M(rg) )
	//    return(Constants::big);
	// Perry introduced this little check when the ice matrix is completely melted away -- in this case
	// set the viscosity to a high number to give the water in the element time to percolate away
	// NOTE MassTransport() is called before SnowCreep, thus this condition should never be met?
	//if ( Edata.theta[ICE] < Constants::eps ) {
	//	return(Constants::big);
	//}
	// Check needed while JAM set!
	if (true && Edata.theta[WATER] > 0.3) {
		return (1.e9 * SnLaws::smallest_viscosity);
	} else if (false && Edata.theta[WATER] >= 0.005) {
		return SnLaws::snowViscosityKOJIMA(Edata, date);
	}
	// Check that you are not in a ice or/and water layer
	if (Edata.theta[ICE] + Edata.theta[WATER] > 0.99) {
		return (1.e9 * SnLaws::smallest_viscosity);
	}

	visc_fudge = SnLaws::snowViscosityFudgeCALIBRATION(Edata, date);
	visc_factor = (sig1*sig1*sig1 / (eps1Dot * visc_fudge*visc_fudge*visc_fudge));
	visc_macro = Edata.neck2VolumetricStrain();
	visc_micro = Edata.neckStressEnhancement();
	eta = (1. / visc_macro) * SnLaws::snowViscosityTemperatureTerm(Te) * visc_factor;
	// NOT YIELDING, LINEAR; sigNeckYield = 0.4 MPa
	if ((visc_micro * sig) <= 100. * sigNeckYield) {
		eta /= visc_micro * sigNeckYield*sigNeckYield;
	// YIELDING, NON-LINEAR
	} else {
		eta /= visc_micro*visc_micro*visc_micro * sig*sig;
	}
  //ANT Quickfix for Antarctica only
	if (SnLaws::setfix && ((date.getJulianDate() - Edata.date.getJulianDate()) > 60.)) {
		eta /= 0.06;
	}
	return eta;
}
