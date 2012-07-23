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
 * @file Meteo.cc
 * @author Michael Lehning and others
 * @version 9.x
 * @date -
 * @bug -
 * @brief Computes missing meteorological information such as friction velocity and roughness length
 * - 29.10.2002: Michael Lehning implements Micromet()
 * - 15.03.2005: Andy and Michi implement stability correction for turbulent fluxes in the hope
 *               that this will also improve the little bit too strong melting of the version 8.1
 */
#include <meteoio/MeteoIO.h>
using namespace mio;

#include <snowpack/Meteo.h>

/************************************************************
* non-static section                                       *
************************************************************/

Meteo::Meteo(const mio::Config& cfg) : canopy(cfg), research_mode(false), useCanopyModel(false), alpine3d(false)
{
	cfg.getValue("ALPINE3D", "SnowpackAdvanced", alpine3d);

	/**
	 * @brief Defines the way to deal with atmospheric stability:
	 * -    0: Standard MO iteration with Paulson and Stearns & Weidner (can be used with BC_CHANGE=0)
	 * -    1: Assume neutral stratification. Should be used with BC_CHANGE=1, i.e., Dirichlet bc but also
	 *         recommended with Neumann b.c., i.e., BC_CHANGE=0
	 * - (-1): Simplified Richardson number stability correction
	 */
	cfg.getValue("NEUTRAL", "Snowpack", neutral);

	/**
	 * @brief Initial estimate of the roughness length for the site; will be adjusted iteratively. \n
	 * Default value and operational mode: 0.002 m
	 */
	cfg.getValue("ROUGHNESS_LENGTH", "Snowpack", roughness_length);

	/**
	 * @brief Defines whether the canopy model is used \n
	 * NOTE: OUT_CANOPY must also be set to dump canopy parameters to file; see Constants_local.h
	 */
	cfg.getValue("CANOPY", "Snowpack", useCanopyModel);

	/**
	 * @brief Define the heights of the meteo measurements above ground (m) \n
	 * Required for surface energy exchange computation and for drifting and blowing snow.
	 */
	cfg.getValue("HEIGHT_OF_WIND_VALUE", "Snowpack", height_of_wind_value);

	cfg.getValue("RESEARCH", "SnowpackAdvanced", research_mode);
}

/**
 * @brief set the atmosphere stability to a given value
 * @param i_neutral stability model (see possible values in constructor)
 */
void Meteo::setStability(const int& i_neutral)
{
	neutral = i_neutral;
}

/**
 * @brief get the atmosphere stability
 * @return stability model (see possible values in constructor)
 */
int Meteo::getStability() const
{
	return neutral;
}

/**
 * @brief Projects precipitations and snow height perpendicular to slope
 * @param hs Height of snow (m)
 * @param precips precipitations (kg m-2)
 * @param slope_angle (deg)
 */
void Meteo::projectPrecipitations(const double& slope_angle, double& precips, double& hs)
{
	precips *= cos(DEG_TO_RAD(slope_angle));
	hs *= cos(DEG_TO_RAD(slope_angle));
}

/**
 * @brief Atmospheric stability correction for wind values.
 * This makes an iteration to find z0 and ustar at the same time
 * @param Mdata
 * @param Xdata
 * @param adjust_VW_height if set to false, assumes a constant measurement height for wind values (default: true, ie.
 * take into account the snow height decreasing the sensor height above the surface)
 */
void Meteo::MicroMet(const SnowStation& Xdata, CurrentMeteo& Mdata, const bool& adjust_VW_height)
{
	const int max_iter = 100;

	// Ideal approximation of pressure and vapor pressure
	const double p0 = Atmosphere::stdAirPressure(Xdata.meta.position.getAltitude());
	const double sat_vap = Atmosphere::waterSaturationPressure(Mdata.ta);
	const double vw = MAX(0.3, Mdata.vw);

	// Initialize snow surface temperature as well as virtual temperatures for stability
	const double t_surf = Xdata.Ndata[Xdata.getNumberOfElements()].T;
	const double ta_v = Mdata.ta * (1. + 0.377 * sat_vap / p0);
	const double t_surf_v = t_surf * (1. + 0.377 * sat_vap / p0);

	// Adjust for snow height if fixed_height_of_wind=false
	const double zref = (adjust_VW_height)? MAX(0.5, height_of_wind_value - (Xdata.cH - Xdata.Ground)) : height_of_wind_value ;
	// In case of ventilation ... Wind pumping displacement depth (m)
	const double d_pump = (SnLaws::wind_pump)? SnLaws::compWindPumpingDisplacement(Xdata) : 0.;

	// Iterate to find atmospheric stability, possibly adjusting z0 to drifting snow and ventilation
	// initial guess (neutral)
	const double eps1 = 1.e-3, eps2 = 1.e-5, a2 = 0.16;
	double z0 = roughness_length;
	double ustar_old, z0_old;
	double psi_m = 0., psi_s = 0., ustar = 0.4 * vw / log((zref - d_pump) / z0);
	double stab_ratio = 0.;
	int iter = 1;
	do {
		iter++;
		if (iter > max_iter) {
			Mdata.z0 = roughness_length;
			Mdata.ustar = 0.4 * vw / log((zref - d_pump) / z0);
			Mdata.psi_s = 0.;
			prn_msg(__FILE__, __LINE__, "wrn", Mdata.date,
			          "Stability correction did not converge (azi=%.0lf, slope=%.0lf) --> assume neutral",
			            Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
			return;
		}
		ustar_old = ustar;
		z0_old = z0;
		z0 = 0.9 * z0_old + 0.1 * (a2 * ustar*ustar / 2. / Constants::g); //update z0
		const double z_ratio = log((zref - d_pump) / z0);

		// Stability corrections
		if (neutral < 0) { // Switch for Richardson
			const double Ri = Constants::g / t_surf_v * (ta_v - t_surf_v) * zref / vw / vw;
			if (Ri < 0.2) { // neutral and unstable
				stab_ratio = Ri;
			} else {
				stab_ratio = Ri/(1.-5.*Ri);
			}
			if (Ri < 0.) { // unstable
				stab_ratio = Ri;
				const double dummy = pow((1. - 15. * stab_ratio), 0.25);
				psi_m = log((0.5 * (1 + dummy*dummy)) * (0.5 * (1 + dummy)) * (0.5 * (1 + dummy)))
				            - 2. * atan(dummy) + 0.5 * Constants::pi;
				psi_s = 2. * log(0.5 * (1 + dummy*dummy));
			} else if (Ri < 0.1999) { // stable
				stab_ratio = Ri / (1. - 5. * Ri);
				psi_m = psi_s = -5. * stab_ratio;
			} else {
				stab_ratio = Ri / (1. - 5. * 0.1999);
				psi_m = psi_s = -5. * stab_ratio;
			}
			ustar = 0.4 * vw / (z_ratio - psi_m);

		} else if (neutral == 0 || (!research_mode && (Mdata.tss > 273.) && (Mdata.ta > 277.))) { // MO Iteration
			ustar = 0.4 * vw / (z_ratio - psi_m);
			const double Tstar = 0.4 * (t_surf_v - ta_v) / (z_ratio - psi_s);
			stab_ratio = -0.4 * zref * Tstar * Constants::g / (t_surf * ustar*ustar);

			if (stab_ratio > 0.) { // stable
				// Stearns & Weidner, 1993
				const double dummy1 = pow((1. + 5. * stab_ratio), 0.25);
				psi_m = log(1. + dummy1) * log(1. + dummy1) + log(1. + dummy1*dummy1)
				            - 1. * atan(dummy1) - 0.5 * dummy1*dummy1*dummy1 + 0.8247; // Original 2.*atan(dummy1) - 1.3333
				// Launiainen and Vihma, 1990
				//psi_m = -17. * (1. - exp(-0.29 * stab_ratio));

				// Holtslag and DeBruin (1988) prepared from Ed Andreas
				//psi_m = psi_s = -(0.7 * stab_ratio + 0.75 * (stab_ratio - 14.28)
				//                    * exp(-0.35 * stab_ratio) + 10.71);

				// Stearns & Weidner, 1993, for scalars
				const double dummy2 = sqrt(1. + 5. * stab_ratio);
				psi_s = log(1. + dummy2) * log(1. + dummy2)
				            - 1. * dummy2 - 0.3 * dummy2*dummy2*dummy2 + 1.2804; // Ori: 2. * dummy2 - 0.66667 * ...
			} else {
				// Stearns & Weidner, 1993 - Must be an ERROR somewhere NOTE maybe - -1. below ;-)
				//const double dummy0 = pow((1.-15. * stab_ratio),0.25);
				//psi_m = log(1. - dummy0) * log(1. - dummy0) + log(1. + dummy0*dummy0)
				//            - 2.*atan(dummy0) - -1. + dummy0 - 0.5086;

				// Paulson - the original
				const double dummy1 = pow((1. - 15. * stab_ratio), 0.25);
				psi_m = 2. * log(0.5 * (1. + dummy1)) + log(0.5 * (1. + dummy1*dummy1))
				            - 2. * atan(dummy1) + 0.5 * Constants::pi;

				// Stearns & Weidner, 1993, for scalars
				const double dummy2 = pow((1. - 22.5 * stab_ratio), 0.33333);
				psi_s = pow(log(1. + dummy2 + dummy2*dummy2), 1.5) - 1.732 * atan(0.577 * (1. + 2. * dummy2)) + 0.1659;
			}
		} else { // NEUTRAL
			psi_m = 0.;
			psi_s = 0.;
		}
	} while ( (fabs(ustar_old - ustar) > eps1) && (fabs(z0_old - z0) > eps2) );

	// Save the values in the global Mdata data structure to use it later
	Mdata.ustar = ustar;
	Mdata.z0 = z0;
	Mdata.psi_s = psi_s;
}

/**
 * @brief Compute measured snow depth change rate to detect growing grass (canopy) vs. snowfall on bare ground
 * @param Mdata
 * @param Xdata
 * @param hs_a3hl6 snow depth average from t_now - 6 h to t_now - 3 h
 * @return whether grass should be detected
 */
bool Meteo::compHSrate(CurrentMeteo& Mdata, const SnowStation& Xdata, const double hs_a3hl6)
{
	if (Xdata.getNumberOfNodes() == Xdata.SoilNode+1) { //Detect only when there is no snow pack yet.
		if ((hs_a3hl6 != Constants::undefined) && (Mdata.hs_a3h != Constants::undefined)) {
			// NOTE we compare two consecutive time spans of 3 hours and take the rate from
			//      the "middle" of the two time spans. hs_rate is in m h-1.
			Mdata.hs_rate = (Mdata.hs_a3h - hs_a3hl6) / 3.;
			return true;
		} else {
			Mdata.hs_rate = Constants::undefined;
			return false;
		}
	} else {
		Mdata.tss_a12h = Constants::undefined;
		Mdata.tss_a24h = Constants::undefined;
		Mdata.hs_rate = Constants::undefined;
		return false;
	}
}

/**
 * @brief
 * \li with CANOPY set:
 * 		In case of an existing canopy, call canopy routine, which computes precipitation, radiation,
 * 		friction velocity and reference temperature for the surface below the canopy.
 * 		Note that solar radiation may change also in dg_cn_Canopy(). \n
 * 		- Mdata->iswr  incoming global solar radiation (direct + diffuse), adapted to canopy
 * 		- Mdata->rswr  reflected global solar radiation (diffuse), adapted to canopy
 * 		- Mdata->ustar friction velocity, adapted to canopy
 * 		- Mdata->z0    roughness length, adapted to canopy
 * 		- Mdata->ea    atmospheric emissivity below canopy, i.e., to give correct
 * 		               longwave radiation as function of air temperature, however
 * 		               modified to include effect of canopy
 * \li without canopy (CANOPY is not set):
 * 		For bare soil as well as snowed-in canopy or some other problems, compute the roughness
 * 		length z0, the friction velocity ustar as well as the atmospheric stability correction
 * 		psi_s for scalar heat fluxes
 * 		- Mdata->ustar friction velocity
 * 		- Mdata->z0    roughness length
 * 		- psi_s        stability correction for scalar heat fluxes
 * @param *Mdata
 * @param *Xdata
 */
void Meteo::compMeteo(CurrentMeteo *Mdata, SnowStation *Xdata)
{
	if (useCanopyModel)
		canopy.runCanopyModel(Mdata, Xdata, roughness_length, height_of_wind_value, alpine3d);

	if (!(useCanopyModel) || Xdata->Cdata.zdispl < 0.) {
		if(alpine3d) MicroMet(*Xdata, *Mdata, false); // for Alpine3D: do not adjust sensor height for snow height
		else MicroMet(*Xdata, *Mdata, true);
	}
}

