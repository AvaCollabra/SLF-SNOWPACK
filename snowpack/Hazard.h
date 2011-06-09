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

/**
 * @file Hazard.h
 * @version 10.02
 * This module contains the hazard computation routines and structures
*/

#ifndef __HAZARD_H__
#define __HAZARD_H__

#include <snowpack/Constants.h>
#include <snowpack/Laws.h>
#include <snowpack/DataClasses.h>
#include <snowpack/Stability.h>
#include <meteoio/MeteoIO.h>

#include <vector>

class Hazard {
	public:
		Hazard(const mio::Config& i_cfg, const double duration);

		void initializeHazard(double *vecDrift, double slope_angle,
		                      std::vector<ProcessDat>& Hdata, std::vector<ProcessInd>& Hdata_ind);

		static double driftIndex(double *old_drift, double drift, const double rho, const int nHours,
		                         double slope_angle, const int shift);

		void getDriftIndex(ProcessDat& Hdata, ProcessInd& Hdata_ind,
                       double *old_drift, double& drift, double slope_angle);

		void getHazardData(ProcessDat& Hdata, ProcessInd& Hdata_ind,
		                   const double& delta_hs6, const double& delta_hs24,
		                   CurrentMeteo& Mdata, SurfaceFluxes& Sdata, SN_ZWISCHEN_DATA& Zdata,
		                   SnowStation& Xdata_station, SnowStation& Xdata_south,
		                   const unsigned int& nSlopes, const bool& virtual_slope);

		static const double typical_slope_length, wind_slab_density;

	private:
		double compDewPointDeficit(double TA, double TSS, double RH);

		double compHoarIndex(double *OldHoar, double new_hoar, int nhour, int new_step);

		static void compMeltFreezeCrust(const SnowStation& Xdata, ProcessDat& Hdata, ProcessInd& Hdata_ind);

		void compHazard(ProcessDat& Hdata, ProcessInd& Hdata_ind, const double& d_hs6, const double& d_hs24,
		                const CurrentMeteo& Mdata, SurfaceFluxes& Sdata, SN_ZWISCHEN_DATA& Zdata,
		                const SnowStation& Xdata);

		const mio::Config& cfg;
		bool research_mode;
		bool enforce_measured_snow_heights;
		bool force_rh_water;
		double sn_dt;
		double i_time_zone;
		int hazard_steps_between;
		int nHz;
		double hoar_density_surf, hoar_min_size_surf;
		static const double minimum_drift, maximum_drift;
};

#endif
