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
 * @file WaterTransport.h
 * @version 10.02
 */

#ifndef __WATERTRANSPORT_H__
#define __WATERTRANSPORT_H__

#include <snowpack/Constants.h>
#include <snowpack/DataClasses.h>
#include <snowpack/Laws_sn.h>

#include <meteoio/MeteoIO.h>

/**
 * @class WaterTransport
 * @version 10.02
 * @bug Prone to bugs at any changes! Be aware!
 * @brief This module contains water transport routines for the 1d snowpack model
 */
class WaterTransport {

	public:
		WaterTransport(const mio::Config& i_cfg);

		void compTransportMass(const CurrentMeteo& Mdata, const double& ql, SnowStation& Xdata, SurfaceFluxes& Sdata);

	private:

		void compSurfaceSublimation(const CurrentMeteo& Mdata, double ql, SnowStation& Xdata, SurfaceFluxes& Sdata);

		void removeElements(SnowStation& Xdata, SurfaceFluxes& Sdata);

		void adjustDensity(SnowStation& Xdata);

		void transportWater(const CurrentMeteo& Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata);

		std::string variant;
		double thresh_rain;
		double sn_dt;
		double hoar_thresh_rh;
		double hoar_thresh_vw;
		double hoar_density_buried, hoar_density_surf, hoar_min_size_buried;
		double minimum_l_element;
		bool useSoilLayers, water_layer, jam;
};
#endif //End of WaterTransport.h
