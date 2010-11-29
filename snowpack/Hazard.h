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
 * This module contains the hazard calculation routines and structures
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
		Hazard(const mio::Config& i_cfg);

		void initializeHazard(const double TimeEnd, double *OldDrift, double SlopeAngle, 
						  std::vector<Q_PROCESS_DAT>& Hdata, std::vector<Q_PROCESS_IND>& Hdata_ind);
	
		static double driftIndex(double *OldDrift, double Drift, double rho, int nhour, int shift);
 
		double hoarIndex(double *OldHoar, double new_hoar, int nhour, int new_step);

		static void calculateMeltFreezeCrust(const SN_STATION_DATA& Xdata, 
									  Q_PROCESS_DAT& Hdata, Q_PROCESS_IND& Hdata_ind);

		void calculateHazard(const double& d_hs6, const double& d_hs24, const SN_STATION_DATA& Xdata,
						 const SN_MET_DATA& Mdata, const int& nAvg, SN_ZWISCHEN_DATA& Zdata, 
						 Q_PROCESS_DAT& Hdata, Q_PROCESS_IND& Hdata_ind, SN_SURFACE_DATA&  Sdata);

		static const double typical_slope_length, wind_slab_density, minimum_drift, maximum_drift;
		static const bool r_in_n;

	private:		
		double calcDewPointDeficit(double TA, double TSS, double RH);

		mio::Config cfg;
		bool research_mode;
		bool enforce_measured_snow_heights;
		bool force_rh_water;
		double sn_dt, calculation_step_length;
		int hazard_steps_between;
		double min_size_hoar_surf, density_hoar_surf;
};

#endif
