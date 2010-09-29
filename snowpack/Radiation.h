/* **********************************************************************************************/
/*                                        VERSION 9.x                                          */
/*                               Derived from RESEARCH VERSION 9.0                             */
/* **********************************************************************************************/
/* **********************************************************************************/
/*  Copyright WSL Institute for Snow and Avalanche Research    SLF-DAVOS           */
/* **********************************************************************************/
/* This file is part of Snowpack.
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

#ifndef __RADIATION_H__
#define __RADIATION_H__

#include <snowpack/Constants.h>
#include <snowpack/Snowpack.h>
#include <snowpack/Laws.h>

#include <meteoio/MeteoIO.h>

/*
 * These classes contain the solar parameters needed to split and project the global
 * incoming irradiance.
*/
class PositionSun {
	public:
		PositionSun() : ecc_corr(0.), decl(0.), eq_time(0.), solar_time(0.), hr_angle(0.), azi_Sacw(0.),
		                azi_Ncw(0.), elev(0.), zen(0.), sunx(0.), suny(0.), sunz(0.), ang_inc(0.) {}

		double ecc_corr;         ///< Correction due to the eccentricity of the earth's orbit
		double decl;             ///< Solar declination: angle between the vector earth/sun and the equatorial plane
		double eq_time;          ///< Equation of time
		double solar_time;       ///< True solar time
		double hr_angle;         ///< Angular displacement of the sun E or W of the local(!) meridian
		
		double azi_Sacw;         ///< Solar azimuth, anticlockwise from South
		double azi_Ncw;          ///< Solar azimuth, clockwise from North
		double elev;             ///< Solar elevation
		double zen;              ///< Solar zenith angle
		
		double sunx, suny, sunz; ///< Components of sun(-earth) vector
		double ang_inc;          ///< Solar incidence (rad)
		
};

class RadiationData {
	public:
		RadiationData() : global_hor(0.), Md(0.), dir_hor(0.), dir_slope(0.), diffsky(0.),
		                  toa_h(0.), pot_dir(0.), pot_diffsky(0.) {}

		double global_hor;       ///< Mean global incoming short wave radiation at the meteo station
		double Md;               ///< Radiation splitting parameter
		double dir_hor;          ///< Incoming direct irradiance at the meteo station
		double dir_slope;        ///< Incoming direct irradiance perpendicular to slope
		double diffsky;          ///< Incoming diffuse irradiance from the sky
		
		double toa_h;            ///< Irradiance on surface normal to the sun vector on top of the atmosphere
		double pot_dir;          ///< Potential incoming direct irradiance perpendicular to the ground
		double pot_diffsky;      ///< Potential incoming diffuse irradiance from the sky
};

/**
 * @class Radiation
 * @author Malcolm Connolly and Michael Lehning (2002/03) and Nora Helbig
 * @author Charles Fierz 2004: provide a common basis for ALPINE3D and SNOWPACK (radiation)
 * @author Nora Helbig, Charles Fierz, and Mathias Bavay 2008: adapt radiation part to A3D-developments by Nora
 * @version 9.x
 * @date -
 * @bug -
 * @brief PURPOSE:
 * - split radiation into direct and diffuse fractions
 * - project radiation onto slope
 * - project precipitation onto slope \n
 * NOTE: Up to Aug 2008, this part of the code was located in Meteo.c
 * .
 * - References \n
 *   TODO complete reference list!
 *   - Bird and Hulstrom (1980, 1981)
 *   - Bourges ???
 *   - Erbs, D.G., S.A. Klein and J.A. Duffie. 1982. \n
 *           Estimation of the diffuse radiation fraction for hourly, daily and monthly-average global radiation.
 *           Sol. Energy, 28(4), 293-304.
 *   - Duffie, J.A. ???. 2006. ???
 *   - Iqbal, M. 1983. \n
 *           An Introduction to Solar Radiation.
 *           Academic Press, Toronto, 390pp.
 *   - Kasten and Young (1989)
 *   - Reindl et al. 1990.
 *   - Spencer, J.W. 1971. \n
 *           Fourier series representation of the position of the sun.
 *           Search, 2(5), p.172.
 *   - Oke, T.R. 1987. \n
 *           Boundary Layer Climates.
 *           2nd ed, Routledge, London, 435pp.
 */
class Radiation {

	public:

		Radiation(const mio::Config& i_cfg);

		void radiationOnSlope(const SN_STATION_DATA& Xdata, SN_MET_DATA& Mdata, SN_SURFACE_DATA& Sdata, 
		                      PositionSun& Psolar, RadiationData& Rdata);

		void flatFieldRadiation(const SN_STATION_DATA& Xdata, SN_MET_DATA& Mdata, 
		                        PositionSun& Psolar, RadiationData& Rdata);

	private:

		double ProjectToHorizontal(const double& slope_component, const double& ang_inc, const double& sunx, 
							  const double& suny, const double& sunz);

		//if slope_angle>0, uses slope_angle, slope_azi for computation, otherwise using sx,sy
		void angleOfIncidence(const double& sx, const double& sy, const double& slope_angle,
		                      const double& slope_azi, PositionSun& Psolar); 

		//day_number=day of year day_number_equi=days since March equinox
		void computeDayNumbers(const mio::Date& date, double& day_number, double& day_number_equi);

		//filling PositionSun's daily parameters
		void computeSolarDailyParameters(const mio::Date& date, PositionSun& Psolar);

		//filling PositionSun's hourly parameters
		void computePositionSun(const double& local_time, const double& Lat, const double& Lon, PositionSun& Psolar);

		//diff/direct balance for RadiationData
		void computeSplittingCoefficient(const PositionSun& solar, const double& thresh_elev, RadiationData& Rdata); 

		void computePotentialRadiation(const PositionSun& Psolar, const double& mean_alb, const double& altitude,
		                               const double& pressure, const double& rh, const double& ta, RadiationData& Rdata);

		//projection of incoming radiation for a given slope
		void projectRadiationOnSlope(const PositionSun& Psolar, const double& Alb, SN_MET_DATA& Mdata, RadiationData& Rdata); 
	private:

		mio::Config cfg;
		int sw_ref;

		static const double thresh_sun_elevation;
};

#endif //END of Radiation.h
