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

#ifndef __SNOWDRIFT_H__
#define __SNOWDRIFT_H__

#include <meteoio/MeteoIO.h>
#include <snowpack/Constants.h>
#include <cmath>
#include <snowpack/Saltation.h>
#include <snowpack/Snowpack.h>
#include <vector>

class Saltation;
struct ElementData;
struct SurfaceFluxes;
struct CurrentMeteo;
struct SnowStation;

/**
 * @class SnowDrift
 * @brief This class contains the computation of local snow drift and the associated erosion
 * @version 10.02
 */
class SnowDrift {

	public:
		SnowDrift(const mio::Config& i_cfg);

		void compSnowDrift(const CurrentMeteo& Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata, double& cumu_hnw);

 	private:
		double compMassFlux(const ElementData& Edata, const double& ustar, const double& slope_angle);

		bool enforce_measured_snow_heights, snow_redistribution; // Will be read from cfg object
		Saltation saltation; // The saltation model used
		double sn_dt;        //Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
		int nSlopes;

		static const double schmidt_drift_fudge;
		static const bool msg_erosion;
}; //End class SnowDrift

#endif //#ifndef __SNOWDRIFT_H__

