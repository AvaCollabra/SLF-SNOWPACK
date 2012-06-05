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
 * @file Laws_sn.cc
 * @version 11.06
 * @brief This module contains (ALL) the constitutive laws for the 1d snowpack model \n
 * @par SnLaws is currently a merely static class (all functions and all variables are static members).
 *   In order to initialize the static member variables they are defined below and if dependent
 *   on the VARIANT in use initialized in setStaticData(const& std::string variant). This
 *   function also initializes static containers (swa_k, swa_fb, swa_pc). When calling
 *   compAlbedo(...) or compSnowViscosity(...) setStaticData might be called in case that
 *   the VARIANT has changed, thus some static variables may be adapted at runtime to the
 *   VARIANT specified.
 *
 * This includes Michael Lehning's ALBEDO law as well as Bob Brown's VISCOSITY and HEAT CONDUCTION models.
 * It also contains the SENSIBLE and LATENT heat CONVECTION coefficients.  Moreover, this code
 * is used by ALL snowpack modules -- creep, temperature and water transport.  Many of the
 * laws, for example the ELASTICITY modulus were programmed up by Marc Christen and taken
 * directly from program HAEFELI, the 2d snowpack model.
 */

#include <snowpack/Laws_sn.h>

using namespace std;
using namespace mio;

/************************************************************
 * static section                                           *
 ************************************************************/

std::string SnLaws::current_variant = "DEFAULT"; // The default variant for SnLaws

/// @brief Switch on or off wind pumping in snow
const bool SnLaws::wind_pump = true;
/// @brief Switch on or off wind pumping in soil
const bool SnLaws::wind_pump_soil = true;

/// @brief Snow extinction coefficients for absorbtion of SW radiation (swa)
//@{
size_t SnLaws::swa_nBands = 5;      ///< Number of bands used
std::vector<double> SnLaws::swa_k;  ///< mean extinction coefficient for pure ice per band
std::vector<double> SnLaws::swa_pc; ///< fraction of sun power spectrum per band
std::vector<double> SnLaws::swa_fb; ///< fudge_bohren
//@}

/**
 * @name SOIL PARAMETERS
 *
 * @brief Define Method and Coefficents for the computation of the influence of soil water
 * content on Evaporation from Bare Soil Layers:
 * - 1 ==> Resistance Approach, see Laws_sn.c
 * - 0 ==> Relative Humidity Approach, see Snowpack.cc
 * - -1 ==> none, assume saturation pressure and no extra resistance
 */
//@{
const bool SnLaws::soil_evaporation = true;

/// @brief Minimum soil surface resistance, 50 sm-1 (van den Hurk et al, 2000)
const double SnLaws::rsoilmin = 50.0;

/// @brief Minimum relative saturation of top soil layer
const double SnLaws::relsatmin = 0.05;

/// @brief Ratio of porosity to tortuosity for Soil
const double SnLaws::alpha_por_tor_soil = 0.05;

/// @brief Pore length for surface soil for ventilation (m)
const double SnLaws::pore_length_soil = 0.01;

/// @brief Field capacity of the soil (1). Above this levels, water begins to drain
const double SnLaws::field_capacity_soil = 0.15;
//@}

/**
 * @name THERMAL CONDUCTIVITY
 * @brief Defines the constants and parameters for computing snow and soil thermal conductivity
 */
//@{
/// @brief Factor controlling ice to ice conduction
const double SnLaws::montana_c_fudge = 0.13;

/// @brief Factor controlling increase in water vapor transport and thus energy transport in wet snow
const double SnLaws::montana_vapor_fudge = 2.5;
//@}

/**
 * @name WIND PUMPING
 * @note soil parameters
 */
//@{
/// @brief Used to decribe advective flux attenuation
const double SnLaws::wind_ext_coef = 0.1;

/// @brief Integral of snow density corresponding to the maximal displacement depth d_pump (m)
const double SnLaws::displacement_coef = 0.7;

/// @brief Ratio of porosity to tortuosity
const double SnLaws::alpha_por_tor = 0.07;
//@}

/// @brief To use J. Hendrikx's parameterization for wind speeds > 2.9 m s-1
const bool SnLaws::jordy_new_snow = false;

/// @brief Defines the smallest allowable viscosity (Pa s) that a viscosity law will return \n
/// Value is DAMM SMALL -- smaller values than this are pretty unrealistic.
const double SnLaws::smallest_viscosity = 1.0e6;

/**
 * @name SNOW VISCOSITY
 * @brief Defines the constants and parameters for computing the snow viscosity.
 * @note {
 * @par
 * The Japanese variant uses a viscosity parametrization by Kojima
 *
 * @par
 * These static variables are only defined below, if you want to change them
 * for your purposes you need to do so in the function SnLaws::setStaticData
 * where these parameters are set according to the VARIANT used
 * }
 * - t_term           : Temperature dependence.
 * 	- _arrhenius_critical, _arrhenius, t_term_837, _stk
 * - visc_ice_fudge   : Empirical constant related to volumetric ice content.
 * - visc_sp_fudge    : Empirical constant related to sphericity of snow grains.
 * - visc_water_fudge : Empirical constant related to volumetric water content.
 * - visc_time_fudge  : Empirical constant related to age of snow (deprecated).
 * - visc_*           : viscosity version:
 * 	- _dflt, _cal, _ant, _897, _837, _stk
 * - setfix           : Quickfix to "harden" antarctic snow (use with CALIBRATION variant only)
 */
//@{
SnLaws::TempDependence SnLaws::t_term = SnLaws::t_term_arrhenius_critical;
SnLaws::ViscosityVersion SnLaws::visc = SnLaws::visc_cal;
double SnLaws::visc_ice_fudge = 1.;
double SnLaws::visc_sp_fudge = 1.;
double SnLaws::visc_water_fudge = 0.;
bool   SnLaws::setfix = false;
//@}

/**
 * @name SNOW ALBEDO
 * @par {
 * These static variables are only defined below, if you want to change them
 * for your purposes you need to do so in the function SnLaws::setStaticData
 * where these parameters are set according to the VARIANT used
 * }
 * - currentAlbedoModel : Albedo model to be used
 * - albedoCage         : Empirical constant related to age of snow, set to zero in Antarctic variant
 */
//@{
SnLaws::AlbedoModel SnLaws::currentAlbedoModel = SnLaws::alb_lehning_2;
bool SnLaws::ageAlbedo = true;
//@}

/**
 * @name Event driven density of new snow
 * @par {
 * These static variables are only defined below, if you want to change them
 * for your purposes you need to do so in the function SnLaws::setStaticData
 * where these parameters are set according to the VARIANT used
 * }
 * - event : type of event, currently only event_wind implemented
 * - 	event_wind_lowlim : lower limit for event_wind
 * - 	event_wind_highlim : upper limit for event_wind
 */
//@{
SnLaws::EventType SnLaws::event = SnLaws::event_none;
double SnLaws::event_wind_lowlim = 0.0;
double SnLaws::event_wind_highlim = 0.0;
//@}
double SnLaws::min_hn_density = 30.;
double SnLaws::max_hn_density = 250.0;

const bool SnLaws::__init = SnLaws::setStaticData("DEFAULT");

/**
 * @brief This function is used to give default values to a bunch of static members of SnLaws
 *        it is called when the helper variable __init is initialized (compile-time)
 *        and everytime the user changes the variant and invokes either compAlbedo(..) or
 *        compSnowViscosity(...)
 * @return always true
 */
bool SnLaws::setStaticData(const std::string& variant)
{
	current_variant = variant;

	if (current_variant == "ANTARCTICA") {
		t_term = t_term_arrhenius_critical;
		visc = visc_dflt;
		visc_ice_fudge = 9.45;
		visc_sp_fudge = 16.5;
		//visc_water_fudge is set to zero by default
		setfix = false;

		event = event_wind;
		event_wind_lowlim = 4.0;
		event_wind_highlim = 7.0;

	} else if (current_variant == "CALIBRATION") {
		// actual calibration; see factors in Laws_sn.cc
		t_term = t_term_arrhenius_critical;
		visc = visc_cal; //visc_897; //visc_ant; //
		// Revision 837 & Steinkogler
// 		t_term = t_term_837; //_stk; // settings from r837 or WSteinkogler (stk)
// 		visc = visc_837; //visc_stk; //
		visc_ice_fudge = 1.;
		visc_sp_fudge = 1.; //see factors in Laws_sn.cc
		visc_water_fudge = 1.;
// 		visc_time_fudge = 11.;
// 		visc_ice_fudge = 0.5; //0.5 // 837=0.5, stk=0.51
// 		visc_sp_fudge = 0.3; //0.3 // 837=0.3, stk=1.19
// 		visc_water_fudge = 31.;
		setfix = false;

	} else {
		t_term = t_term_arrhenius_critical;
		visc = visc_dflt;
		visc_ice_fudge = 9.45;
		visc_sp_fudge = 16.5;
		visc_water_fudge = 33.;
		setfix = false;
	}

	if (current_variant == "JAPAN")
		currentAlbedoModel = alb_nied;
	else
		currentAlbedoModel = alb_lehning_2;

	if (current_variant == "ANTARCTICA")
		SnLaws::ageAlbedo = false;
	else
		SnLaws::ageAlbedo = true;

	// snow extinction coefficients
	double k_init[5]  = {0.059, 0.180, 0.525, 4.75, 85.23}; // values in use since r140
	double fb_init[5] = {29., 15., 5., 9., 35.}; // values in use since r140
	double pc_init[5] = {16.3, 16.8, 15.4, 31.0, 20.5}; // values in use since r140
// 	double k_init[5]  = {0.00125, 0.0198, 0.179, 3.74, 178.}; // 31 Jul 2008
// 	double fb_init[5] = {0.5, 0.05, 0.02, 3., 5.}; // fbn1 01 Aug 2008
	//double fb_init[5] = {0.01, 0.05, 0.05, 3., 3.}; // fbn0 31 Jul 2008
	swa_k.resize(swa_nBands);
	swa_pc.resize(swa_nBands);
	swa_fb.resize(swa_nBands);
	std::copy(k_init, k_init+swa_nBands, swa_k.begin());
	std::copy(pc_init, pc_init+swa_nBands, swa_pc.begin());
	std::copy(fb_init, fb_init+swa_nBands, swa_fb.begin());

	return true;
}

/**
 * @name THERMAL CONDUCTIVITY OF ICE
 * @brief Based on master thesis of Tobias Hipp, who used relationships by Ling & Yhang (2005).
 * @version 11.03
 * @param Temperature Temperature (K)
 * @return Thermal conductivity of ice
 */
double SnLaws::conductivity_ice(const double& Temperature)
{
	double ki = 0.4685 + (488.19)/(Temperature);
	return ki;
}

/**
 * @name THERMAL CONDUCTIVITY OF WATER
 * @brief Based on master thesis of Tobias Hipp, who used relationships by Ling & Yhang (2005).
 * @version 11.03
 * @param Temperature Temperature (K)
 * @return Thermal conductivity of water
 */
double SnLaws::conductivity_water(const double& Temperature)
{
	double kw = 0.11455 + 1.6318E-3 * (Temperature);
	return kw;
}

/**
 * @name SNOW ALBEDO (defined as an absolute value--and not as a rate of change)
 * @brief Various parameterizations for snow albedo
 * @version 11.05
 * currentAlbedoModel:
 * 	- alb_lehning_[012] : Statistical models of surface snow albedo based on measurements
 *      from the Weissfluhjoch study plot (SWin and SWout, K&Z CM21)
 * 	- alb_nied : The Japanese version of alb_lehning_2
 * @param i_fixed_albedo value
 * @param Edata
 * @param Tss Snow surface temperature (K)
 * @param Mdata
 */
double SnLaws::parameterizedSnowAlbedo(const double& i_fixed_albedo, const ElementData& Edata,
                                       const double& Tss, const CurrentMeteo& Mdata)
{
	double Alb = Constants::min_albedo, Alb1;
	const double Ta = Mdata.ta;
	double age = Mdata.date.getJulianDate() - Edata.depositionDate.getJulianDate();

	if (i_fixed_albedo != Constants::undefined) {
		Alb = i_fixed_albedo;
	} else if ((SnLaws::ageAlbedo && (age > 365.)) || (Edata.mk % 10 == 7)) {
		Alb = Constants::glacier_albedo;
	} else {
		switch (currentAlbedoModel) {
		case alb_lehning_0: {
			double sqrt_age, lwi;
			const double weight=0.1;
			const double a = -0.2, b = 1.3, c = -0.012, d = -0.011, e = 0.0024, f = 0.018;
			const double g = -7.8e-6, h = -3.1e-3, i = 3.5e-4, j = 1.6e-7, k = -2.4e-7;
			const double l = -9.2e-6, m = 6.7e-5, n = -7.2e-5, o = 2.0e-4;

			if (age < 0.001)
				sqrt_age=Constants::eps;
			else
				sqrt_age=sqrt(age);
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

			if (Edata.mk%100 > 19) {
				mf = 1.;
				// av *= exp(-age/1700.);
			}
			Alb1 = Crho*Edata.Rho + Clwc*Edata.theta[WATER] + Cdd*Edata.dd + Csp*Edata.sp
			                        + Cmf*mf + Crb*Edata.rb +  Cta*Ta + Ctss*Tss
			                          + Cv*Mdata.vw+ Cswout*Mdata.rswr + Cta_tss*Ta*Tss;
			Alb = av + log(1.0 + Alb1);
			break;
		}
		case alb_lehning_2: {
			//TODO: this perfoms very badly (if not completly wrong) for (very?) wet snowpack
			//for example, February 2007 in Davos with very warm weather resulting in (measured?) albedos of 0.3 ...
			double av = 0.8042; // Value of original regression
			if (!SnLaws::ageAlbedo) { // NOTE clean antarctic snow
				age = 0.;
				av = 0.7542; // estimated from comparison with measurements at Dome C
			} else {
				age = MIN(30., age);
			}
			const double inter = 1.442;
			const double Cage = -0.000575, Cta = -0.006, Cv = 0.00762, Clwc = -0.2735;
			const double Crho = -0.000056, Crh = 0.0333, Crb = -0.301, Crg = 0.175;
			const double Cdd = 0.064, Csp = -0.0736, Ctss = 0.00459, Cswout = -0.000101;
			Alb1 = inter + Cage*age + Crho*Edata.Rho + Clwc*Edata.theta[WATER]
			           + Cdd*Edata.dd + Csp*Edata.sp + Crg*Edata.rg + Crb*Edata.rb
			               + Cta*Ta + Ctss*Tss + Cv*Mdata.vw + Cswout*Mdata.rswr
			                   + Crh*Mdata.rh;
			if (Alb1 > 0.) {
				Alb = MAX(Constants::min_albedo, MIN(Constants::max_albedo, av + log(Alb1)));
			} else {
				Alb = Constants::min_albedo;
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date, "Alb1=%lf set Alb to %lf", Alb1, Alb);
			}
			if (age > 30.) {
				prn_msg(__FILE__, __LINE__, "msg+", Mdata.date,
						"albedo=%f : age=%f rho=%f th_w=%f ta=%f rswr=%f rb=%f sp=%f type=%d",
						Alb, age, Edata.Rho, Edata.theta[WATER], Ta, Mdata.rswr, Edata.rb, Edata.sp, Edata.type);
			}
			break;
		}
		case alb_nied: { // by H. Hirashima (NIED, Nagaoka, Japan)
			const double av = 0.75;
			const double inter = 1.005;
			const double Cage = -0.00016*10.0, Cta = -0.000249*2.0, Cv = 0.00578, Clwc = -2.15;
			const double Crho = -0.000047, Crh = 0.129, Crb = -0.306, Crg = 0.107;
			const double Cdd = 0.076, Csp = 0.00964, Ctss = -0.000166, Cswout = -1.8e-5;

			Alb1 = inter + Crho*Edata.Rho + Clwc*Edata.theta[WATER] + Cdd*Edata.dd + Csp*Edata.sp
			           + Crg*Edata.rg + Crb*Edata.rb + Cta*Ta + Ctss*Tss + Cv*Mdata.vw
			               + Cswout*Mdata.rswr + Crh*Mdata.rh + Cage*age;

			if (Alb1 > 0) {
				Alb = MAX(Constants::min_albedo, MIN(Constants::max_albedo, av + log(Alb1)));
			} else {
				Alb = Constants::min_albedo;
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date, "Alb1=%lf set Alb to %lf", Alb1, Alb);
			}
			break;
		}
		default:
			prn_msg(__FILE__, __LINE__, "err", Date(), "Albedo model %d not implemented yet!", currentAlbedoModel);
			throw IOException("The required snow albedo model is not implemented yet!", AT);
			break;
		}
	}
	return(Alb);
}

/**
 * @brief Helens Solution to Radiation Transfer \n
 * NOTE on fudge_bohren (fb): Larger values increase extinction --> Energy stays on top;
 * originally not band dependent, set to 10.0 for Neumann and to 5.0 for Dirichlet BC
 * @version 12.04
 * @param Xdata
 * @param I0 net shortwave radiation (W m-2)
 * @param multistream
 */
void SnLaws::compShortWaveAbsorption(SnowStation& Xdata, const double& I0, const bool& multistream)
{
	size_t nE, e, bottom_element;
	double I0_band, dI, Ks;

	ElementData *EMS; // Avoids dereferencing the pointer

	EMS = &Xdata.Edata[0];
	nE = Xdata.getNumberOfElements();
	if (nE == 0)
		return;
	else if (Xdata.SoilNode > 0)
		bottom_element = Xdata.SoilNode - 1;
	else
		bottom_element = Xdata.SoilNode;

	for (e = bottom_element; e < nE; e++)
		EMS[e].sw_abs = 0.;

	// Compute absorbed radiation
	if (multistream) {
		for (size_t ii = 0; ii < swa_nBands; ii++) {
			I0_band = I0 * swa_pc[ii] / 100.;
			for (e = nE-1; e > bottom_element; e--) {
				Ks = swa_fb[ii] * 0.84 * sqrt( 1000. * swa_k[ii] / (2. * EMS[e].rg) ) * EMS[e].Rho
				       / Constants::density_ice;
				if (EMS[e].mk%10 != 9)
					dI = I0_band * (1. - exp(-(Ks * (EMS[e].L))));  // Radiation absorbed by element e in band i
				else
					dI = 0.; // Water is transparent
				EMS[e].sw_abs += dI;
				I0_band -= dI;
			}
			EMS[bottom_element].sw_abs += I0_band;
		}
	} else { // ad hoc "1-Band" model
		I0_band = I0;
		for (e = nE-1; e > bottom_element; e--) {
			if (EMS[e].mk%10 != 9)
				dI = I0 * (1. - exp(-(EMS[e].extinction() * (EMS[e].L))));  // Radiation absorbed by element e
			else
				dI = 0.; // Water is transparent
			EMS[e].sw_abs += dI;
			I0_band -= dI;
		}
		EMS[bottom_element].sw_abs += I0_band;
	}

	for (e = bottom_element; e < nE; e++) {
		if (EMS[e].sw_abs < 0.) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "NEGATIVE Shortwave Radiation %lf absorbed by element %d (nE=%d)",
				   EMS[e].sw_abs, e, nE);
			throw IOException("SnLaws::compShortWaveAbsorption did not complete successfully", AT);
		}
	}
}

/**
 * @brief Computes the displacement depth in case of ventilation
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
		if ((hw + Xdata.Edata[e].L * Xdata.Edata[e].Rho) < SnLaws::displacement_coef)
			d_pump += Xdata.Edata[e].L;
		else
			d_pump += (SnLaws::displacement_coef - hw) / Xdata.Edata[e].Rho;
		hw += Xdata.Edata[e].L * Xdata.Edata[e].Rho;
	}
	return (d_pump);
}

/**
 * @brief Computes the wind pumping velocity at the surface
 * @version 10.01
 * @param Mdata
 * @param d_pump Displacement depth (m)
 * @return Wind pumping velocity (m s-1)
 */
double SnLaws::compWindPumpingVelocity(const CurrentMeteo& Mdata, const double& d_pump)
{
	//TODO Limit v_pump?
	return (Mdata.ustar / 0.4 * log((Mdata.z0 + d_pump) / Mdata.z0));
}

/**
 * @brief Computes the wind gradient near the snow-cover surface
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
 * @version 11.03: thermal conductivity made temperature dependent.
 * @param Edata
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
	if ((Edata.rg > 0.) && (Edata.rg < 10000.)) {
		C_eff_soil_max = Edata.theta[SOIL] * c_mineral + Edata.theta[WATER]
				* SnLaws::conductivity_water(Edata.Te) + Edata.theta[ICE]
					* SnLaws::conductivity_ice(Edata.Te);

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
		if (Edata.theta[WATER] > 0.0001)
			C_eff_soil += MAX(0.27,(alpha1 + alpha2 * weight) * log(alpha3 * Edata.theta[WATER]));
		else
			C_eff_soil += 0.27;
		C_eff_soil = MIN(C_eff_soil_max, C_eff_soil);
	} else {
		C_eff_soil = Edata.soil[SOIL_K] + Edata.theta[WATER] * SnLaws::conductivity_water(Edata.Te)
                       + Edata.theta[ICE] * SnLaws::conductivity_ice(Edata.Te);
	}

	// Now check for possible ERRORS
	if (!(C_eff_soil >= 0.1*Edata.soil[SOIL_K] && C_eff_soil < 10.)) //HACK
		prn_msg(__FILE__, __LINE__, "wrn", Date(), "Thermal Conductivity of Soil: %lf", C_eff_soil);

	// In case of dry soil, simply use the given conductivity with a possible ventilation part
	if (SnLaws::wind_pump_soil)
		C_eff_soil += SnLaws::alpha_por_tor_soil * Constants::specific_heat_air
		                * Edata.soil[SOIL_RHO] * SnLaws::pore_length_soil * SnLaws::pore_length_soil * dvdz;
	return(C_eff_soil);
}

/**
 * @brief Computes the enhancement factor for water vapor transport in wet snow
 * @version 9Y.mm
 * @param Xdata
 * @param e Element number
 * @return Enhancement factor
 */
double SnLaws::compEnhanceWaterVaporTransportSnow(const SnowStation& Xdata, const int& e)
{
	int e1 = e;
	double vapor_enhance = 1.;

	while ((e1 >= MAX(signed(Xdata.SoilNode), e-7)) && (vapor_enhance < SnLaws::montana_vapor_fudge) ) {
		//TODO check limit on theta_water
		if ((Xdata.Edata[e1].theta[WATER] > 0.001) && (Xdata.Edata[e].theta[ICE] < 0.7))
			vapor_enhance = SnLaws::montana_vapor_fudge;
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
 * @version 8.10 (??)
 * @param Edata
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
	Te = MIN(Edata.Te, Edata.melting_tk);

	// Check for elements with no ice and assume they contain only water
	if (Edata.theta[ICE] < Snowpack::min_ice_content)
		return(Constants::conductivity_water);

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

	C_eff  = SnLaws::montana_c_fudge * C1 * (C2 + C3 + C4 + C5) * (2.0 - Edata.dd) * (1.0 + pow(Edata.theta[ICE], 1.7)) * (0.5 + Te*Te / (Edata.melting_tk*Edata.melting_tk));

	if (!((C_eff < 5.*Constants::conductivity_ice) && (C_eff > 0.2*Constants::conductivity_air)) && !ALPINE3D) {
		prn_msg(__FILE__, __LINE__, "wrn", Date(), "Conductivity out of range (0.2*Constants::conductivity_air=%.3lf, 5.*Constants::conductivity_ice=%.3lf):", 0.2 * Constants::conductivity_air, 5. * Constants::conductivity_ice);
		prn_msg(__FILE__, __LINE__, "msg-", Date(), "C_eff: %lf  C_1: %lf  C_2: %lf  C_3: %lf  C_4: %lf  C_5: %lf", C_eff, C1, C2, C3, C4, C5);
		prn_msg(__FILE__, __LINE__, "msg-", Date(), "C_Air: %lf  C_Water: %lf  C_Ice: %lf", Constants::conductivity_air, Constants::conductivity_water, Constants::conductivity_ice);
		prn_msg(__FILE__, __LINE__, "msg-", Date(), "Ap: %lf  kap: %lf  Thet_i: %lf  Thet_w: %lf  T: %lf", Ap, kap, Edata.theta[ICE], Edata.theta[WATER], Te );
		prn_msg(__FILE__, __LINE__, "msg-", Date(), "type: %3d  rg: %.3lf mm  rb: %.3lf mm N3: %lf  Lel: %.3lf mm", Edata.type, Edata.rg, Edata.rb, Edata.N3, M_TO_MM(Edata.L));
	}

	if (!(C_eff < Constants::conductivity_ice))
		C_eff = Constants::conductivity_ice;
	if (!(C_eff > Constants::conductivity_air))
		C_eff = Constants::conductivity_air;
	// Now introduce the effect of windpumping
	if (SnLaws::wind_pump /* && Edata.theta[WATER] < 0.001 */)
		C_eff += SnLaws::alpha_por_tor * Constants::specific_heat_air * Edata.Rho * dvdz / (C1*C1);
	return (C_eff);
}

/**
 * @brief SENSIBLE HEAT EXCHANGE COEFFICIENT (Surface Energy Exchange)
 * @version 9Y.mm
 * @param Mdata
 * @param Xdata
 * @param height_of_meteo_values Height at which meteo parameters are measured
 * @return Exchange coefficient for sensible heat (1)
*/
double SnLaws::compSensibleHeatCoefficient(const CurrentMeteo& Mdata, const SnowStation& Xdata, const double& height_of_meteo_values)
{
	double c;
	double z, karman = 0.4, lrat;

	z = MAX (0.5, height_of_meteo_values - Xdata.cH + Xdata.Ground);
	if ((Xdata.cH - Xdata.Ground) > 0.03) {
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
 * @version 9Y.mm
 * @param Mdata
 * @param Xdata
 * @param height_of_meteo_values Height at which meteo parameters are measured
 * @return Latent heat flux (W m-2)
 */
double SnLaws::compLatentHeat_Rh(const CurrentMeteo& Mdata, SnowStation& Xdata, const double& height_of_meteo_values)
{
	double beta, eA, eS;
	double Vp1, Vp2;
	double th_w_ss;

	const double T_air = Mdata.ta;
	const double Tss = Xdata.Ndata[Xdata.getNumberOfElements()].T;

	// Vapor Pressures
	if (Xdata.getNumberOfElements() == 0)
		th_w_ss = 0.0;
	else
		th_w_ss = Xdata.Edata[Xdata.getNumberOfElements()-1].theta[WATER];

	// TODO The part below needs to be rewritten in a more consistent way !!!
	//      In particular, look closely at the condition within compLatentHeat()
	eA = Mdata.rh * lw_SaturationPressure(T_air);
	Vp1 = lw_SaturationPressure(Tss);
	Vp2 = lw_SaturationPressure(Tss);

	// First, the case of no snow
	if (Xdata.getNumberOfNodes() == Xdata.SoilNode + 1 && Xdata.getNumberOfElements() > 0) {
		if ( Tss < Xdata.Edata[Xdata.getNumberOfElements()-1].melting_tk) {
			eS = Vp1 ;
		} else {
			/*
			 * Soil evaporation can now be computed using the Relative Humidity approach below,
			 * or a Resistance approach modifying the ql value instead of the eS. The latter
			 * function is defined in compLatentHeat, and the Switch SnLaws::soil_evaporation is found
			 * in Laws_sn.h
			*/
			if (SnLaws::soil_evaporation && th_w_ss < Xdata.Edata[Xdata.SoilNode-1].soilFieldCapacity()) {
				eS = Vp2 * 0.5 * ( 1 - cos (MIN (Constants::pi, th_w_ss * Constants::pi
				         / (Xdata.Edata[Xdata.SoilNode-1].soilFieldCapacity() * 1.6))));
			} else {
				eS = Vp2;
			}
		}
	} else {
		// for snow assume saturation
		const double melting_tk = (Xdata.getNumberOfElements()>0)? Xdata.Edata[Xdata.getNumberOfElements()-1].melting_tk : Constants::melting_tk;
		if (Tss < melting_tk)
			eS = Vp1;
		else
			eS = Vp2;
	}
	// Now the latent heat
	beta = SnLaws::compLatentHeat(Mdata, Xdata, height_of_meteo_values);

	return (beta * (eA - eS));
}

/**
 * @brief LATENT HEAT EXCHANGE (Surface Energy Exchange)
 * @version 9Y.mm
 * @param Mdata
 * @param Xdata
 * @param height_of_meteo_values Height at which meteo parameters are measured
 * @return Latent heat flux (W m-2)
 */
double SnLaws::compLatentHeat(const CurrentMeteo& Mdata, SnowStation& Xdata, const double& height_of_meteo_values)
{
	double c, eS, eA;
	double z, karman = 0.4, lrat;

	z = MAX (0.5, height_of_meteo_values - Xdata.cH + Xdata.Ground);
	if ((Xdata.cH - Xdata.Ground) > 0.03) {
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
	if ((Xdata.getNumberOfNodes() == Xdata.SoilNode + 1) && (Xdata.getNumberOfElements() > 0)
		    && (Xdata.Ndata[Xdata.getNumberOfElements()].T >= Xdata.Edata[Xdata.getNumberOfElements()-1].melting_tk)
		      && (SnLaws::soil_evaporation == 1)) {
		eA = Mdata.rh * lw_SaturationPressure(Mdata.ta);
		eS = lw_SaturationPressure(Xdata.Ndata[Xdata.getNumberOfElements()].T);
		if (eS >= eA) {
			c = 1. / c + SnLaws::rsoilmin / MAX(SnLaws::relsatmin, MIN(1.,
			                                     Xdata.Edata[Xdata.getNumberOfElements()-1].theta[WATER]
			                                         / Xdata.Edata[Xdata.SoilNode-1].soilFieldCapacity()));
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
 * @version 9Y.mm
 * @param t_snow Snow surface temperature (K)
 * @param t_atm Temperature of the atmosphere, i.e., air (K)
 * @param e_atm Emissivity of the atmosphere (1)
 * @return LW radiation coefficient (?)
 */
double SnLaws::compLWRadCoefficient(const double& t_snow, const double& t_atm, const double& e_atm)
{
	return (Constants::emissivity_snow * Constants::stefan_boltzmann
	           * ((t_snow*t_snow) + (sqrt(e_atm) * t_atm*t_atm))
	               * (t_snow + (pow(e_atm, 0.25) * t_atm)));

}

/**
 * @brief Event driven new-snow density
 * @param i_event:
 * - event_wind: rho = 250.3 kg m-3 @ 4 m s-1; rho = 338 kg m-3 @ 7 m s-1 Antarctica
 * @param Mdata  Meteorological input
 */
double SnLaws::newSnowDensityEvent(const std::string& variant, const SnLaws::EventType& i_event,
                                   const CurrentMeteo& Mdata)
{
	if (variant != SnLaws::current_variant)
		setStaticData(variant);

	switch (i_event) {
	case event_wind: {
		const double rho_0=361., rho_1=33.;
		if ((Mdata.vw_avg >= event_wind_lowlim) && (Mdata.vw_avg <= event_wind_highlim))
			return (rho_0*log10(Mdata.vw_avg) + rho_1);
		else
			return Constants::undefined;
		break;
	}
	default:
		prn_msg(__FILE__, __LINE__,"err", Date(),
		        "No new snow density parameterization for event type %d", i_event);
		throw IOException("Event type not implemented yet!", AT);
		break;
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
double SnLaws::newSnowDensityPara(const std::string& i_hn_model,
                                  double TA, double TSS, double RH, double VW, double HH)
{
	double rho_hn;

	TA  = K_TO_C(TA);
	TSS = K_TO_C(TSS);
	RH  *= 100.;
	HH  = floor(HH);

	if (i_hn_model == "LEHNING_OLD") {
		double alpha=70., beta=30., gamma=10., delta=0.4;
		double eta=30., phi=6.0, mu=-3.0, nu=-0.5;
		rho_hn = alpha + beta*TA + gamma*TSS +  delta*RH + eta*VW + phi*TA*TSS + mu*TA*VW + nu*RH*VW;
		if (jordy_new_snow && (VW > 2.9))
			rho_hn = newSnowDensityHendrikx(TA, TSS, RH, VW);

	} else if (i_hn_model == "LEHNING_NEW") {
		double alpha=90., beta=6.5, gamma=7.5, delta=0.26;
		double eta=13., phi=-4.5, mu=-0.65, nu=-0.17, om=0.06;
		rho_hn = alpha + beta*TA + gamma*TSS +  delta*RH + eta*VW + phi*TA*TSS + mu*TA*VW + nu*RH*VW + om*TA*TSS*RH;
		// Ad hoc temperature correction
		if (TA < -10.)
			rho_hn = MIN(rho_hn, alpha*(1. + 0.05*(TA + 10.)));
		// Increase snow density under snow transport conditions
		if ((!jordy_new_snow) && (VW > 5.)) {
			rho_hn = 90. + (rho_hn - 30.)*0.9;
		} else if (jordy_new_snow && (VW > 2.9)) {
			rho_hn = newSnowDensityHendrikx(TA, TSS, RH, VW);
		}

	} else if (i_hn_model == "BELLAIRE") {
		double arg;
		double alpha=3.946, beta=0.07703, zeta=0.0001701, eta=0.02222, mu=-0.05371;
		// Transformations based on natural logarithm!!!
		VW = MAX(1., VW);
		arg = alpha + beta*TA + zeta*HH + eta*log(VW) + mu*TA*log(VW);
		rho_hn = exp(arg);

	} else if (i_hn_model == "ZWART") {
		double arg;
		double beta01=3.28, beta1=0.03, beta02=-0.36, beta2=-0.75, beta3=0.3;
		VW = MAX(2., VW);
		RH = 0.8; // ori: MIN(1., RH/100.); see asin(sqrt()) below
		if (TA < -14.)
			arg = beta01 + beta1*TA + beta2*asin(sqrt(RH)) + beta3*log10(VW);
		else
			arg = beta01 + beta1*TA + beta02 + beta2*asin(sqrt(RH)) + beta3*log10(VW); // += beta2*TA;
		rho_hn = pow(10., arg);

	} else if (i_hn_model == "PAHAUT") {
		rho_hn = 109. + 6.*(C_TO_K(TA) - Constants::melting_tk) + 26.*sqrt(VW);

	} else {
		prn_msg(__FILE__, __LINE__, "err", Date(),
		        "New snow density parameterization '%s' not available",
		        i_hn_model.c_str());
		exit(EXIT_FAILURE);
	}
	return(MIN(max_hn_density, MAX(min_hn_density, rho_hn)));
}

/**
 * @brief Jordy Hendrikx' new snow density parameterization for strong winds (> 2.9 m s-1)
 * @note To be used with Lehning's models only!
 * @param ta  Air temperature (degC)
 * @param tss Snow surface temperature (degC)
 * @param rh  Relative air humidity (%)
 * @param vw  Mean wind velocity (m s-1)
 * @return New snow density
 */
double SnLaws::newSnowDensityHendrikx(const double ta, const double tss, const double rh, const double vw)
{
	const double alpha=91., beta=-35., gamma=-1.1, delta=49., eta=32.,  phi=4.6;
	return (alpha + beta*ta + gamma*rh +  delta*vw + eta*tss + phi*ta*vw);
}

/**
 * @brief Computes the density of new snow. The options are:
 * - PARAMETERIZED:
 * 	- ZWART: Costijn Zwart's model (elaborated 2006; in use since 4 Dec 2007
 * 	- LEHNING_NEW: Improved model by M. Lehning, incl. ad-hoc wind & temperature effects (used until 06/07)
 * 	- LEHNING_OLD: First model by M. Lehning
 *       @note {models by M. Lehning can be augmented with a parameterization for winds > 2.9 m s-1
 *              worked out by J. Hendrikx => set jordy_new_snow in Laws_sn.cc}
 * 	- BELLAIRE: Sascha Bellaire's model (elaborated 2007; used summer/fall 2007)
 * 	- PAHAUT: Edmond Pahaut's model, introduced Sep 1995 in CROCUS by G. Giraud
 * - EVENT: Driven by event type, that is,
 * 	- event_wind: Implemented 2009 by Christine Groot Zwaaftink for Antarctic variant
 * - MEASURED: Use measured new snow density read from meteo input (RHO_HN must be set)
 * - fixed: Fixed new snow density by assigning HN_DENSITY a number > 0.
 * @param i_hn_density type of density computation
 * @param i_hn_density_model to use
 * @param Mdata Meteorological input
 * @param Xdata Snow cover data
 * @param tss Snow surface temperature (K)
 * @param variant which is currently used
 */
double SnLaws::compNewSnowDensity(const std::string& i_hn_density, const std::string& i_hn_density_model,
                                  const CurrentMeteo& Mdata, const SnowStation& Xdata, const double& tss,
                                  const std::string& variant)
{
	double rho;

	if (i_hn_density == "PARAMETERIZED") {
		rho = newSnowDensityPara(i_hn_density_model,
		                         Mdata.ta, tss, Mdata.rh, Mdata.vw,
		                         Xdata.meta.position.getAltitude());
	} else if (i_hn_density == "EVENT") {
		rho = newSnowDensityEvent(variant, event, Mdata);
	} else if (i_hn_density == "MEASURED") {
		if (Mdata.rho_hn != Constants::undefined) {
			rho = Mdata.rho_hn; // New snow density as read from input file
		} else if (Mdata.hnw > 0. ) {
			if (i_hn_density_model == "SURFACE_SNOW") {
				rho = Xdata.Edata[Xdata.getNumberOfElements()-1].Rho;
			} else {
				rho = newSnowDensityPara(i_hn_density_model,
				                         Mdata.ta, tss, Mdata.rh, Mdata.vw,
				                         Xdata.meta.position.getAltitude());
			}
		} else {
			rho = Constants::undefined;
		}
	} else { // "FIXED"
		if (!IOUtils::convertString(rho, i_hn_density, std::dec))
			throw ConversionFailedException("Cannot convert  '"+i_hn_density+"' to double", AT);
	}

	return rho;
}

/**
 * @brief NEW SNOW VISCOSITY (dendritic snow, i.e., dd > 0.) \n
 * Actual version : ml_lw_VS_Lehning from r7.7 \n
 * This is Michael's viscosity routine, which is not a function of micro-structure, but which
 * is nonetheless pretty important because it is numerically STABLE and does predict decent
 * settling rates, sometimes a bit too high. This routine was (is) used for NEW or WET snow.
 * @version 9Y.mm
 * @return Viscosity of new snow (Pa s)
 */
double SnLaws::NewSnowViscosityLehning(const ElementData& Edata)
{
	double rho_hn = MAX(min_hn_density, Edata.Rho);

	if (rho_hn > 913.) //upper boundary
		return (1.e9 * smallest_viscosity);

	if (Edata.theta[WATER] > 0.001)
		//return (0.01*pow(rho_hn, 4.7));
		return(0.0001 * pow(rho_hn, 5.5));
	else
		return (0.007 * pow(rho_hn, (4.75 - K_TO_C(Edata.Te) / 40.)));
}

/**
 * @brief Computes the temperature term of viscosity
 * @version 11.06
 * @param Te Element temperature (K)
 * @return Temperature term of snow viscosity
 */
double SnLaws::snowViscosityTemperatureTerm(const double& Te)
{
	const double Q = 67000.; // Activation energy for defects in ice (J mol-1)

	switch (SnLaws::t_term) {
	case t_term_arrhenius_critical: {
		const double Q_fac = 0.39; // Adjust Q to snow; from Schweizer et al. (2004): 0.24
		const double criticalExp = 0.7; //0.5; //0.3; //
		const double T_r = 265.15; // Reference temperature (K), from Schweizer et al. (2004)
		return ((1. / lw_ArrheniusLaw(Q_fac * Q, Te, T_r))
		             * (0.3 * pow((Constants::melting_tk - Te), criticalExp) + 0.4));
		break;
	}
	case t_term_arrhenius:
		return (1. / lw_ArrheniusLaw(Q, Te, 263.));
		break;
	case t_term_stk: // Master thesis, September 2009
		return (0.35 * sqrt(274.15 - Te));
		break;
	default: // as of revision 243, used up to revision 837 (deprecated)
		return (9. - 8.7 * exp(0.015 * (Te - Constants::melting_tk)));
		break;
	}
}

/**
 * @brief Computes the additional stress due to loading rate
 * @version 11.06
 * @param i_viscosity_model see compSnowViscosity()
 * @param Edata
 * @param date current
 * @return Initial stress (Pa); note it is a negative value!
 */
double SnLaws::compLoadingRateStress(const std::string& i_viscosity_model, ElementData& Edata, const mio::Date& date)
{
	if (i_viscosity_model == "CALIBRATION")
		return loadingRateStressCALIBRATION(Edata, date);
	else
		return loadingRateStressDEFAULT(Edata, date);
}

double SnLaws::loadingRateStressDEFAULT(ElementData& Edata, const mio::Date& date)
{
	const double age = MAX(0., date.getJulianDate() - Edata.depositionDate.getJulianDate());

	double sigReac = 15.5 * Edata.CDot * exp(-age/101.);
	if (Edata.theta[WATER] > Constants::eps)
		sigReac *= 0.37 * (1. + Edata.theta[WATER]);
	return sigReac;
}

double SnLaws::loadingRateStressCALIBRATION(ElementData& Edata, const mio::Date& date)
{
	const double sigTension = 0.11;  // Ice surface tension (N m-2)

	Edata.EDot = 0.;
	switch (visc) {
	case visc_dflt: case visc_cal: case visc_ant:  { // new calibration
		double sigMetamo = 0.;
		const double age = MAX(0., date.getJulianDate() - Edata.depositionDate.getJulianDate());
		double sigReac = 15.5 * Edata.CDot * exp(-age/101.);
		if (Edata.theta[WATER] > Constants::eps)
			sigReac *= 0.37 * (1. + Edata.theta[WATER]); // 0.2 ; 0.37
// 		sigReac = 15.5 * Edata->CDot * exp(-age/101.)
// 		              * MAX(0.1, 1. - 9.*Edata.theta[WATER]);
// 		sigReac = 15.5 * Edata->CDot * exp(-age/101.)
// 		              * MAX(0., 1. - pow(Edata.theta[WATER], 1./3.));
// 		sigReac = 15.9 * Edata->CDot * exp(-age/101.);
// 		sigReac = 17.3 * Edata->CDot * exp(-age/101.); //tst2: 553. //tst1: 735. //
		Edata.EDot = sigReac;
		if (false && (Edata.dd > Constants::eps)) {
			sigMetamo = 3.7e2 * Metamorphism::ddRate(Edata);
			if (false && Edata.dd < 0.9 ) {
				sigMetamo *= pow(Edata.dd, 1./2.);
			}
		}
		return sigReac;
// 		return (sigReac + sigMetamo);
// 		return 0.;
		break;
	}
	case visc_897: { // r897
		double sigReac = 0., sigMetamo = 0.;
		const double age = MAX(0., date.getJulianDate() - Edata.depositionDate.getJulianDate());
		sigReac = 15.9 * Edata.CDot * exp(-age/101.); //tst2: 553. //tst1: 735. //
		Edata.EDot = sigReac;
		//if ( 1 && (Edata->dd > 0.2) /*((Edata.dd < 0.9) && (Edata.dd > 0.3))*/ ) {
		if (1 && (Edata.dd > Constants::eps) /*((Edata->dd < 0.9) && (Edata->dd > 0.3))*/) {
			//sigMetamo = 73.5 * mm_ddRate(Edata) * sigTension / MM_TO_M(Edata.rg); // ori: -3.0; -23.
			//sigMetamo = 139.7 * mm_ddRate(Edata) * sigTension / MM_TO_M(Edata.rg); // 2010-10-22
			//sigMetamo = 29.0 * mm_ddRate(Edata) / MM_TO_M(Edata.rg); // 2010-10-23
			sigMetamo = 37.0e3 * Metamorphism::ddRate(Edata); // 2010-10-23
			if (0 && Edata.dd < 0.9) {
				//sigMetamo /= Edata.dd; // pow(Edata.dd, 1./2.);
				sigMetamo *= pow(Edata.dd, 1./2.);
			}
		}
		return (sigReac + sigMetamo);
		break;
	}
	case visc_837: case visc_stk: { // as of revision 837
		double sig0 = 0.;
		if ((Edata.dd < 0.9) && (Edata.dd > 0.3)) {
			double facIS = 3.; // default r712
			if (SnLaws::visc == SnLaws::visc_stk)
				facIS = -1.5; //-1.1; //-0.5; //
			sig0 = facIS * Metamorphism::ddRate(Edata) * sigTension / MM_TO_M(Edata.rg);
		}
		return sig0;
		break;
	}
	default:
		prn_msg(__FILE__, __LINE__, "err", Date(),
				"visc=%d not a valid choice for loadingRateStress!", visc);
		throw IOException("Choice not implemented yet!", AT);
		break;
	}
}

/**
 * @brief Determines the fudge factor for viscosity \n
 * This fudge factor takes into account bond-ice imperfections and the effect of liquid water
 * @version 11.06
 * @param Edata
 * @return Fudge factor for snow viscosity
 */
double SnLaws::snowViscosityFudgeDEFAULT(const ElementData& Edata)
{
	double visc_fudge, sp_fudge;

	double ice_fudge = SnLaws::visc_ice_fudge / Edata.theta[ICE];
	ice_fudge *= (1. - logisticFunction(Edata.theta[ICE], 0.019, 0.15))
	                 * pow(Edata.theta[ICE], 0.77);

	if (Edata.mk%100 >= 20 && Edata.theta[WATER] < 0.005)
		sp_fudge = 0.;
	else
		sp_fudge = SnLaws::visc_sp_fudge;
	sp_fudge *= (1. - logisticFunction(Edata.sp, 0.37, 0.21))
	                 * Edata.sp*Edata.sp*Edata.sp;

	visc_fudge = ice_fudge * (1. + sp_fudge) + 0.1;

	visc_fudge *= (1. + SnLaws::visc_water_fudge * Edata.theta[WATER]
	                   * (1. - Edata.theta[ICE]));
	return visc_fudge;
}

/**
 * @brief To calibrate the fudge factor for viscosity
 * @note This is the fudge playground for calibrations used in SnowViscosityCALIBRATION()
 * @version 10.11
 * @param Edata
 * @param date current date
 * @return Fudge factor for snow viscosity
 */
double SnLaws::snowViscosityFudgeCALIBRATION(const ElementData& Edata, const mio::Date& date)
{
	double visc_fudge, sp_fudge;
	double const age = MAX(0., date.getJulianDate() - Edata.depositionDate.getJulianDate());
	double thresh_rho1 = 1., thresh_rho2 = 1.; // Thresholds for enhanced viscosity
	bool use_thresh = false;

	if (Edata.mk%100 >= 20 && Edata.theta[WATER] < 0.005)
		sp_fudge = 0.;
	else
		sp_fudge = SnLaws::visc_sp_fudge;

	switch (visc) {
	case visc_cal: case visc_dflt: { // Calibration currently under test
		double ice_fudge = (0.67*14.1*SnLaws::visc_ice_fudge / Edata.theta[ICE]); // 0.67; 0.59; r897: 0.51
		ice_fudge *= (1. - logisticFunction(Edata.theta[ICE], 0.019, 0.15))
		                  * pow(Edata.theta[ICE], 0.77);
// 		ice_fudge *= 0.67*13.7 * (1. - lwsn_logisticFunction(Edata.theta[ICE], 0.019, 0.15))
// 		ice_fudge *= 0.67*23.7*(1. - lwsn_logisticFunction(Edata.theta[ICE], 0.009, 0.11))
// 		                 * Edata.theta[ICE];

		sp_fudge *= 0.67*24.7;  // 0.67 ; 0.77 ; r897: 1.19 //
		sp_fudge *= (1. - logisticFunction(Edata.sp, 0.37, 0.21))
		                 * Edata.sp*Edata.sp*Edata.sp;
// 		sp_fudge *= 0.67*24.7*(1. - lwsn_logisticFunction(Edata.sp, 0.35, 0.17))
// 		                * Edata.sp*Edata.sp*Edata.sp;
// 		sp_fudge *= 0.67*(1. - lwsn_logisticFunction(Edata.sp, 1.0, 0.077)) * sqrt(Edata.sp);

		visc_fudge = ice_fudge * (1. + sp_fudge) + 0.1;

		visc_fudge *= (1. + 33.*visc_water_fudge * Edata.theta[WATER]
		                   * (1. - Edata.theta[ICE])); // 31. ; 33. ; 35.
// 		visc_fudge *= (1. + 31.*visc_water_fudge * Edata.theta[WATER]
// 		                   * (1. - Edata.theta[ICE]*Edata.theta[ICE])); // ori r897
// 		                   * (0.67 - 0.29*Edata.theta[ICE]));
// 		                   * (1.37 - Edata.theta[ICE]));
// 		                   * (0.8 - 0.6*pow(Edata.theta[ICE], 0.5)));
// 		                   * MAX(0.2, MIN(0.8, 1.16 - 1.2*Edata.theta[ICE])));
// 		                   * (MAX(0., 0.8 - Edata.theta[ICE]*Edata.theta[ICE])));
		break;
	}
	case visc_897: case visc_ant: { // Calibration fall 2010 & Antarctica
		double ice_fudge = 0.51*SnLaws::visc_ice_fudge / Edata.theta[ICE]; //
		if (SnLaws::visc != SnLaws::visc_ant) {
			ice_fudge *= 1. - logisticFunction(Edata.theta[ICE], 0.431, 0.067) * Edata.theta[ICE];
		} else {  // antarctic
			ice_fudge *= 1. - logisticFunction(Edata.theta[ICE], 0.351, 0.097) * Edata.theta[ICE]
			               * (1. - Edata.theta[ICE]);
		}

		sp_fudge *= 1.19; //
		sp_fudge *= (1. - logisticFunction(Edata.sp, 1.0, 0.077)) * sqrt(Edata.sp);

		visc_fudge = ice_fudge * (1. + sp_fudge) + 0.1;

		visc_fudge *= (1. + 31.*SnLaws::visc_water_fudge * Edata.theta[WATER]
		                   * (1. - Edata.theta[ICE]*Edata.theta[ICE]));
		break;
	}
	case visc_stk: { // Walter Steinkogler's playground; master thesis, September 2009
		const double visc_time_fudge = 8.;
		visc_fudge = visc_time_fudge / exp(age / 35.);
		visc_fudge += (0.5*SnLaws::visc_ice_fudge / Edata.theta[ICE])
		                 + (0.3*sp_fudge * sqrt(Edata.sp))
		                   + (3. * Edata.theta[WATER] / Edata.theta[ICE] * 0.5
		                          * (1. - Edata.theta[ICE]));
		use_thresh=true;
		thresh_rho1 = 0.5; // rho_dry > 458.5 kg m-3
		thresh_rho2 = 0.7;
		break;
	}
	case visc_837: { // as of revision 712, used up to r837 (deprecated)
		const double visc_time_fudge = 11.;
		visc_fudge = visc_time_fudge * (1. - sqrt(MIN(1., age / 77.)))
		                  * (1. + MIN(0.3, (263.15 - Edata.Te) / 17.));
		visc_fudge += (0.5*SnLaws::visc_ice_fudge / Edata.theta[ICE])
		                 + (0.3*sp_fudge * sqrt(Edata.sp))
		                   + (3. * Edata.theta[WATER] / Edata.theta[ICE]
		                          * 0.5 * (1. - Edata.theta[ICE]));
		use_thresh=true;
		thresh_rho1 = 0.5; // rho_dry > 458.5 kg m-3
		thresh_rho2 = 0.7;
		break;
	}
	default:
		prn_msg(__FILE__, __LINE__, "err", Date(),
		        "visc=%d not a valid choice for SnowViscosityFudge!", visc);
		throw IOException("Choice not implemented yet!", AT);
		break;
	}

	if (use_thresh) { // Default false
		if (Edata.theta[ICE] > thresh_rho1)
			visc_fudge *= (1. - Edata.theta[ICE]);
		if (Edata.theta[ICE] > thresh_rho2)
			visc_fudge *= (1. - Edata.theta[ICE]);
	}
	return visc_fudge;
}

/**
 * @name SNOW VISCOSITY
 * @brief Computes snow viscosity according to the following options:
 * - DEFAULT      : default model (corresponding to last calibration)
 * - KOJIMA       : according to parameterization by Kojima
 * - CALIBRATION  : used for calibration purposes
 * @param variant see VARIANT
 * @param i_viscosity_model to use
 * @param Edata
 * @param date current
 */
double SnLaws::compSnowViscosity(const std::string& variant, const std::string& i_viscosity_model,
                                 ElementData& Edata, const mio::Date& date)
{
	if (variant != SnLaws::current_variant)
		setStaticData(variant);

	if (i_viscosity_model == "KOJIMA") {
		return snowViscosityKOJIMA(Edata);
	} else if (i_viscosity_model == "CALIBRATION") {
		return snowViscosityCALIBRATION(Edata, date);
	} else if (i_viscosity_model == "DEFAULT") {
		return snowViscosityDEFAULT(Edata);
	} else {
		throw InvalidArgumentException("Unknown viscosity model: "+i_viscosity_model, AT);
	}
}

/**
 * @brief SNOW VISCOSITY (all types of snow)
 * @todo revise description
 * @version 11.02
 * @param Edata
 * @return Snow viscosity (Pa s)
 */
double SnLaws::snowViscosityDEFAULT(ElementData& Edata)
{
	const double eps1Dot = 1.76e-7;    // Unit strain rate (at stress = 1 MPa) (s-1)
	const double sig1 = 0.5e6;         // Unit stress from Sinha's formulation (Pa)
	const double sigNeckYield = 0.4e6; // Yield stress for ice in neck (Pa)
	const double sig = -Edata.C;       // Overburden stress, that is, absolute value of Cauchy stress (Pa)
	double visc_macro, visc_micro;     // Structure related multiplying factors
	double visc_fudge, visc_factor;    // Fit and reference parameters
	double eta;                        // Viscosity (Pa s)

	// Check needed while JAM set!
	if (Edata.theta[WATER] > 0.3)
		return (1.e9 * SnLaws::smallest_viscosity);
	// Check that you are not in a ice or/and water layer
	if (Edata.theta[ICE] + Edata.theta[WATER] > 0.99)
		return (1.e9 * SnLaws::smallest_viscosity);

	visc_fudge = SnLaws::snowViscosityFudgeDEFAULT(Edata);
	visc_factor = (sig1*sig1*sig1 / (eps1Dot * visc_fudge*visc_fudge*visc_fudge));
	visc_macro = Edata.neck2VolumetricStrain();
	visc_micro = Edata.neckStressEnhancement();
	double Te = MIN(Edata.Te, Edata.melting_tk);
	eta = (1. / visc_macro) * SnLaws::snowViscosityTemperatureTerm(Te) * visc_factor;
	// HACK multiply sigNeckYield by 100. to avoid yielding on purpose
	if ((visc_micro * sig) <= 100. * sigNeckYield) // NOT YIELDING, LINEAR
		eta /= visc_micro * sigNeckYield*sigNeckYield;
	else // YIELDING, NON-LINEAR
		eta /= visc_micro*visc_micro*visc_micro * sig*sig;
	return eta;
}

/**
 * @brief SNOW VISCOSITY according to formulation by Kojima
 * @version 9.10
 * @param Edata
 * @return Snow viscosity (Pa s)
 */
double SnLaws::snowViscosityKOJIMA(const ElementData& Edata)
{
	return (8.64e6 * exp(0.021*Edata.Rho));
}

/**
 * @brief Calibrate snow viscosity
 * NOTE This is the test or playground version for calibrating settling
 * @version 10.11
 * @param Edata
 * @param date current
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

	Te = MIN(Edata.Te, Edata.melting_tk);

	// TODO Check whether the two commented checks below are needed!
	// If the element length is SMALLER than the grain size then the thing aint settling ....
	//  if( Edata.L <= 2.*MM_TO_M(rg) )
	//    return(Constants::big);
	// Perry introduced this little check when the ice matrix is completely melted away -- in this case
	// set the viscosity to a high number to give the water in the element time to percolate away
	// NOTE MassTransport() is called before SnowCreep, thus this condition should never be met?
	//if (Edata.theta[ICE] < Constants::eps) {
	//	return(Constants::big);
	//}
	// Check needed while JAM set!
	if (true && (Edata.theta[WATER] > 0.3)) {
		return (1.e9 * SnLaws::smallest_viscosity);
	} else if (false && (Edata.theta[WATER] >= 0.005)) {
		return SnLaws::snowViscosityKOJIMA(Edata);
	}
	// Check that you are not in a ice or/and water layer
	if (Edata.theta[ICE] + Edata.theta[WATER] > 0.99)
		return (1.e9 * SnLaws::smallest_viscosity);

	visc_fudge = SnLaws::snowViscosityFudgeCALIBRATION(Edata, date);
	visc_factor = (sig1*sig1*sig1 / (eps1Dot * visc_fudge*visc_fudge*visc_fudge));
	visc_macro = Edata.neck2VolumetricStrain();
	visc_micro = Edata.neckStressEnhancement();
	eta = (1. / visc_macro) * SnLaws::snowViscosityTemperatureTerm(Te) * visc_factor;
	// HACK multiply sigNeckYield by 100. to avoid yielding on purpose
	if ((visc_micro * sig) <= 100. * sigNeckYield) // NOT YIELDING, LINEAR
		eta /= visc_micro * sigNeckYield*sigNeckYield;
	else // YIELDING, NON-LINEAR
		eta /= visc_micro*visc_micro*visc_micro * sig*sig;
	//ANT Quickfix for Antarctica only
	if (SnLaws::setfix && ((date.getJulianDate() - Edata.depositionDate.getJulianDate()) > 60.))
		eta /= 0.06;

	return eta;
}
