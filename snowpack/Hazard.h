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

#ifndef HAZARD_H
#define HAZARD_H

#include <snowpack/Constants.h>
#include <snowpack/Laws.h>
#include <snowpack/DataClasses.h>
#include <snowpack/Stability.h>
#include <meteoio/MeteoIO.h>

#include <vector>

struct SN_ZWISCHEN_DATA;
struct SN_STATION_DATA;
struct SN_MET_DATA;
struct SN_SURFACE_DATA;

/// Structure of double values to supplement snowpack model
typedef struct Q_PROCESS_DAT {
	mio::Date date;
	char stat_abbrev[16];
	int  loc_for_snow;
	int  loc_for_wind;
	// Version used, date, user, ...
	char sn_version[MAX_STRING_LENGTH];          ///< SNOWPACK version
	char sn_computation_date[MAX_STRING_LENGTH];
	double sn_jul_computation_date;
	mio::Date sn_compile_date;
	char sn_user[MAX_STRING_LENGTH];             ///< SNOWPACK user
	int nHz;               ///< Number of hazard steps

	// Snow depth
	double ch;             ///< cm
	// SWE, total liquid water content and runoff
	double swe;
	double tot_lwc;
	double runoff;
	// Surface hoar index
	double dewpt_def;      ///< degC
	double hoar_size;      ///< mm
	double hoar_ind6;      ///< kg m-2
	double hoar_ind24;     ///< kg m-2
	// Drift index
	double wind_trans;     ///< cm
	double wind_trans24;   ///< cm
	// New snow depths
	double hns_half_hour;  ///< cm
	double hns3;           ///< cm
	double hns6;           ///< cm
	double hns12;          ///< cm
	double hns24;          ///< cm
	double hns72;          ///< cm
	double hns72_24;       ///< cm;
	// New snow water equivalents
	double wc_half_hour;   ///< kg m-2
	double wc3;            ///< kg m-2
	double wc6;            ///< kg m-2
	double wc12;           ///< kg m-2
	double wc24;           ///< kg m-2
	double wc72;           ///< kg m-2
	// Stability indices
	int stab_class1;       ///< Stability class 1,3,5
	int stab_class2;       ///< Profile type 0..10
	double stab_index1;
	double stab_index2;
	double stab_index3;
	double stab_index4;
	double stab_index5;
	double stab_height1;   ///< cm
	double stab_height2;   ///< cm
	double stab_height3;   ///< cm
	double stab_height4;   ///< cm
	double stab_height5;   ///< cm
	// Special parameters
	double crust;          ///< cm
	double en_bal;         ///< kJ m-2
	double sw_net;         ///< kJ m-2
	double t_top1;         ///< degC
	double t_top2;         ///< degC
} Q_PROCESS_DAT;

typedef struct Q_PROCESS_IND {
	short stat_abbrev;
	short loc_for_snow;
	short loc_for_wind;
	// Data
	short ch;

	short swe;
	short tot_lwc;
	short runoff;

	short dewpt_def;
	short hoar_size;
	short hoar_ind6;
	short hoar_ind24;

	short wind_trans;
	short wind_trans24;

	short hns3;
	short hns6;
	short hns12;
	short hns24;
	short hns72;
	short hns72_24;

	short wc3;
	short wc6;
	short wc12;
	short wc24;
	short wc72;

	short stab_class1;
	short stab_class2;
	short stab_index1;
	short stab_index2;
	short stab_index3;
	short stab_index4;
	short stab_index5;
	short stab_height1;
	short stab_height2;
	short stab_height3;
	short stab_height4;
	short stab_height5;

	short crust;
	short en_bal;
	short sw_net;
	short t_top1;
	short t_top2;
} Q_PROCESS_IND;

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
		mio::Config cfg;
		bool enforce_measured_snow_heights;
		bool force_rh_water;
		double sn_dt, calculation_step_length;
		int hazard_steps_between;
};
#endif //End of Hazard.h
