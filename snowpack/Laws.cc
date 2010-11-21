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
#include <snowpack/Laws.h>
#include <snowpack/Constants.h>

/**
 * @file Laws.c
 * @version 10.02
 * @bug
 * @brief This contains some physical laws or parametrizations that are used by Snowpack but are independent of its data structures.
 * They can therefore be used outside of Snowpack.
 */

/**
 * @brief Calculates the emissivity of a body
 * @author Mathias Bavay
 * @version 10.02
 * @param lwr longwave radiation emitted by the body (W m-2)
 * @param T   surface temperature of the body (K)
 */
double lw_emissivity(const double lwr, const double T)
{
	return ( lwr / (Constants::stefan_boltzmann * (T*T*T*T)) );
}

/**
 * @brief Calculates the long wave radiation coming from a body knowing its emissivity
 * @author Mathias Bavay
 * @version Y.mm
 * @param ea emissivity of the body (0-1)
 * @param T   surface temperature of the body (K)
 */
double lw_lwr(const double ea, const double T)
{
	return ( ea * (Constants::stefan_boltzmann * (T*T*T*T)) );
}

/**
 * @brief Extrapolate air temperature
 * @author Mathias Bavay
 * @version Y.mm
 * @param ta_ref const double
 * @param alti_ref const double
 * @param altitude const double
 * @return double
 */
double lw_TairLapseRate(const double ta_ref, const double alti_ref, const double altitude)
{
	double ta;
	const double lapse_rate = 0.0065; // K m-1

	ta = ta_ref - (altitude - alti_ref) * lapse_rate;

	return(ta);
}

/**
 * @brief Equation for standard atmosphere; Stull "Meteorology for Scientists and Engineers", p.12-13
 * standard pressure and standard temperature together with height as geopotential height
 * with standard sea-level temperature 288.15 K and molecular scale temperature gradient
 * -6.5 K/km
 * @author Mathias Bavay
 * @version Y.mm
 * @param altitude const double
 * @return double
 */
double lw_AirPressure(const double altitude)
{
	double p;
	const double p0 = 101325.; 		// Air and standard pressure in Pa
	const double lapse_rate = 0.0065;	// K m-1
	const double sea_level_temp = 288.15;	// K
	const double expo = Constants::g / (lapse_rate * Constants::gas_constant_air);
	const double R0 = 6356766.0;		// Earth's radius in m

	p = p0 * pow( 1. - ( (lapse_rate * R0 * altitude) / (sea_level_temp * (R0 + altitude)) ), expo );

	return(p);
}

double lw_WetBulbTemperature(const double L, const double T, const double RH, const double altitude)
{
	// L:	latent heat of vaporisation
	const double mixing_ratio = 0.622;	//gas_constant(dry_air) / gas_constant(water_vapor)
	const double p = lw_AirPressure(altitude);
	const double Vp = lw_SaturationPressure(T);

	return ( T - (RH*Vp - Vp) * mixing_ratio * L / p / Constants::specific_heat_air );

}

/**
 * @brief Saturation pressure
 * This function computes the saturation vapor pressure (in Pa). The saturation pressures are
 * determined ( approximated ) using the empirical formulas of "Beziehung von Murray"
 * for flat ice or water surfaces (others: Magnus formula, the precise formula after Golff-Gratch):
 * e  = p0 * exp ( Qlh * ( T - T0 ) / ( R * T * T0 ) )
 * where
 * - p0, T0:	triple point pressure and temperature of water
 * - Qlh:	latent heat of evaporation, (water) or sublimation (ice)
 * Note: the zero point of the Celsius scale of Murray is at 273.16 K!!!
 * from F.W. Murray (1967), J. of Applied Meteorology (transformed equation of Teten (1930) who gave the Magnus
 * formula)
 * @author Mathias Bavay \n Michael Lehning
 * @version Y.mm
 * @param T air temperature in K (const double)
 * @return double
 */
double lw_SaturationPressure(const double T)
{
	double exp_p_sat; // exponent
	const double p0 = 610.78; // triple point pressure of water
	double c2, c3; // varying constants

	if ( T < 273.16 ) { // for a flat ice surface
		c2 = 21.88;
		c3 = 7.66;
	} else { // for a flat water surface
		c2 = 17.27;
		c3 = 35.86;
	}
	exp_p_sat = c2 *  (T - 273.16) / (T - c3);

	return( p0 * exp( exp_p_sat ) );
}

/**
 * @brief Conversion from relative humidity to dew point temperature
 * This function convert a relative humidity to a dew point temperature
 * @author Mathias Bavay \n Michael Lehning
 * @version Y.mm
 * @param RH relative humidity, between 0 and 1 (double)
 * @param TA air temperature in K (double)
 * @param force_water switch to force computation above water (if set to 1)
 * This is needed in order to work with measurements, since WMO defines RH measurements
 * as always above water (for any temperature) (const short int)
 * @return double
 */
double RhtoDewPoint(double RH, double TA, const bool& force_water)
{
	//Convert a Relative Humidity into a dew point temperature
	//TA is in Kelvins, RH between 0 and 1, returns Td in Kelvins
	TA = K_TO_C(TA);
	double Es, E, Tdw, Tdi; //saturation and current water vapor pressure
	const double Aw = 611.21, Bw = 17.502, Cw = 240.97;	//parameters for water
	const double Ai = 611.15, Bi = 22.452, Ci = 272.55;	//parameters for ice
	const double Tfreeze = 0.;			//freezing temperature
	const double Tnucl = -16.0;			//nucleation temperature
	const double di = 1. / ((TA - Tnucl) * (TA - Tnucl) + 1e-6);		//distance to pure ice
	const double dw = 1. / ((Tfreeze - TA) * (Tfreeze - TA) + 1e-6);	//distance to pure water

	//in order to avoid getting NaN if RH=0
	RH += 0.0001;
	if ((TA >= Tfreeze) || force_water) {//above freezing point, water
		Es = Aw * exp( (Bw * TA) / (Cw + TA) );
		E = RH * Es;
		Tdw = ( Cw * log(E / Aw) ) / ( Bw - log(E / Aw) );
		return C_TO_K(Tdw);
	}
	if (TA < Tnucl) { //below nucleation, ice
		Es = Ai * exp( (Bi * TA) / (Ci + TA) );
		E = RH * Es;
		Tdi = ( Ci * log(E / Ai) ) / ( Bi - log(E / Ai) );
		return C_TO_K(Tdi);
	}

	//no clear state, we do a smooth interpolation between water and ice
	Es = Ai * exp( (Bi*TA) / (Ci + TA) );
	E = RH * Es;
	Tdi = ( Ci * log(E / Ai) ) / ( Bi - log(E / Ai) );

	Es = Aw * exp( (Bw * TA) / (Cw + TA) );
	E = RH * Es;
	Tdw = ( Cw * log(E / Aw) ) / ( Bw - log(E / Aw) );

	return C_TO_K( (di / (di + dw) * Tdi + dw / (di + dw) * Tdw) );
}

/**
 * @brief Conversion from dew point temperature to relative humidity
 * This function converts a dew point temperature to a relative humidity between 0 and 1
 * @author Mathias Bavay \n Michael Lehning
 * @version Y.mm
 * @param TD dew point temperature (double)
 * @param TA air temperature in K (double)
 * @param force_water switch to force computation above water (if set to 1)
 * This is needed in order to work with measurements, since WMO defines RH measurements
 * as always above water (for any temperature) (const short int)
 * @return double
 */
double DewPointtoRh(double TD, double TA, const bool& force_water)
{
	//Convert a dew point temperature into a Relative Humidity
	//TA, TD are in Kelvins, RH is returned between 0 and 1
	TA = K_TO_C(TA);
	TD = K_TO_C(TD);
	double Es, E, Rhi, Rhw, Rh;			//saturation and current water vapro pressure
	const double Aw = 611.21, Bw = 17.502, Cw = 240.97;	//parameters for water
	const double Ai = 611.15, Bi = 22.452, Ci = 272.55;	//parameters for ice
	const double Tfreeze = 0.;			//freezing temperature
	const double Tnucl = -16.0;			//nucleation temperature
	const double di = 1. / ((TA - Tnucl) * (TA - Tnucl) + 1e-6);		//distance to pure ice
	const double dw = 1. / ((Tfreeze - TA) * (Tfreeze - TA) + 1e-6);	//distance to pure water

	if ((TA >= Tfreeze) || force_water) {
		//above freezing point, water
		Es = Aw * exp( (Bw * TA) / (Cw + TA) );
		E  = Aw * exp( (Bw * TD) / (Cw + TD) );
		Rhw = (E / Es);
		if (Rhw > 1.) {
			return 1.;
		} else {
			return Rhw;
		}
	}
	if (TA < Tnucl) {
		//below nucleation, ice
		Es = Ai * exp( (Bi * TA) / (Ci + TA) );
		E  = Ai * exp( (Bi * TD) / (Ci + TD) );
		Rhi = (E / Es);
		if (Rhi > 1.) {
			return 1.;
		} else {
			return Rhi;
		}
	}

	//no clear state, we do a smooth interpolation between water and ice
	Es = Ai * exp( (Bi * TA) / (Ci + TA) );
	E  = Ai * exp( (Bi * TD) / (Ci + TD) );
	Rhi = E / Es;

	Es = Aw * exp( (Bw * TA) / (Cw + TA) );
	E  = Aw * exp( (Bw * TD) / (Cw + TD) );
	Rhw = E / Es;

	Rh = (di / (di + dw) * Rhi + dw / (di + dw) * Rhw);
	if(Rh > 1.) {
		return 1.;
	} else {
		return Rh;
	}
}

/**
 * @brief LONGWAVE ATMOSPHERIC EMISSIVITY ESTIMATION
 * Estimate the longwave radiation using the formula from Brutsaert -- "On a Derivable
 * Formula for Long-Wave Radiation From Clear Skies", Journal of Water Resources
 * Research, Vol. 11, No. 5, October 1975, pp 742-744.
 * Alternative: Satterlund (1979): Water Resources Research, 15, 1649-1650.
 * @author Perry Bartelt ??
 * @version Y.mm
 * @param e0 Water vapor saturation pressure (Pa)
 * @param ta Air temperature (K)
 * @return double
 */
double lw_LW_Brutsaert(const double e0, const double ta)
{
	return (1.24 * pow((0.01 * e0 / ta), (1. / 7.)));
}


/**
 * @brief Estimate the longwave radiation using the formula from Omstedt, 1990.
 * The air pressure e0 is in Pa and the cloud fraction is in (1)
 * @author Michael Lehning
 * @version Y.mm
 * @param e0 Air pressure in Pa (const double)
 * @param cloud_frac cloud fraction, between 0 and 1 (const double)
 * @return double
 */
double lw_Omstedt(const double e0, const double cloud_frac)
{
	return (0.97 * (0.68 + 0.0036 * sqrt(e0)) * (1. + 0.18 * cloud_frac * cloud_frac));
}


/**
 * @brief Estimate the residual water content RWC by Vol from work by Coleou and Lesaffre, 1998.
 * Ann. Glaciol., 26, 64-68.
 * Experimental range:
 * - density unsoaked: 235 to 580
 * - density soaked: 328 to 589 kg m-3
 * - RWC by Vol 0.049 to 0.029
 * Note that function will limit range to 0.0264 to 0.08 RWC by Vol
 * @author Charles Fierz
 * @version 8.11
 * @param theta_i Volumetric fraction of ice in the snow
 * @return double
 */
double lw_SnowResidualWaterContent(const double theta_i)
{
	const double max_content = 0.08; //TODO: is it really needed?
	double res_wat_cont;

	if ( theta_i > 0.23 ) {
		res_wat_cont = 0.0264 + 0.0099 * (1. - theta_i) / theta_i;
	} else {
		res_wat_cont = 0.08 - 0.1023 * (theta_i - 0.03);
	}

	//maximum returned value: max_content
	if(res_wat_cont > max_content) {
		return max_content;
	} else {
		return res_wat_cont;
	}
}

/**
 * @brief Computes air emissivity (1) \n
 * Uses either incoming long wave radiation or either Brutsaert (clear sky) or Omstedt (cloudiness) parametrization \n
 * NOTE ta and rh must be checked values, best containing no NODATA!
 * @author Mathias Bavay et al.
 * @version 10.04
 * @param input
 * - input between 0 and 1 : fractional cloud cover
 * - input greater than 1 : incoming longwave radiation (W m-2)
 * - input less than or equal to 0 : use clear sky parameterization of Brutsaert
 * @param ta Checked air temperature (K)
 * @param rh Checked relative humidity (1)
 * @return air emissivity in range [MIN_AIR_EMISSIVITY,1.] (1)
 */
double lw_AirEmissivity(const double input, const double ta, const double rh, const double min_air_emissivity)
{
	double ea;

	if(input > 1.) {
		ea = input/(Constants::stefan_boltzmann*ta*ta*ta*ta);
	} else {
		if((ta == NODATA) || (rh == NODATA)) {
			return NODATA;
		}
		const double e0 = rh * lw_SaturationPressure(ta); // (Pa)
		if(input <= 0.) {
			ea = lw_LW_Brutsaert(e0, ta);
		} else {
			ea = lw_Omstedt(e0, input);
		}
	}
	if(ea < min_air_emissivity) {
		return min_air_emissivity;
	} else if(ea < 1.) {
		return ea;
	} else {
		return 1.;
	}
}

/**
 * @brief Computes an Arrhenius-type temperature dependency
 * @author Charles Fierz
 * @version 9.11
 */
double lw_ArrheniusLaw(const double ActEnergy, const double T, const double T_ref)
{
	const double R = 8.31;		// gas constant J mol-1 K-1

	return exp((ActEnergy / R) * (1. / T_ref - 1. / T));
}

/*
 * End of Laws.c
*/
