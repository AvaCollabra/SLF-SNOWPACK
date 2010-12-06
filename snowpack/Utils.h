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
 * @file Queries.h
 * @version 10.02
 * @date -
 * @brief This header file contains the definition of the structures required to handle I/O
 */

#ifndef __UTILS_H__
#define __UTILS_H__

#include <snowpack/DataClasses.h>
#include <snowpack/Constants.h> 
#include <snowpack/Laws.h>
#include <meteoio/MeteoIO.h>

#include <cstdarg> // needed for va_list
#include <string>
#include <cstring>
#include <vector>

/*
 * FUNCTION PROTOTYPES
 */
#ifdef GNU	//in this case, GCC can check the format arguments for types, number, ...
void prn_msg(const char *theFile, int theLine, const char *msg_type, double JulianDate, const char *format, ...)
__attribute__ ((format (printf, 5, 6)));
#else
void prn_msg(const char *theFile, int theLine, const char *msg_type, double JulianDate, const char *format, ...);
#endif

int qr_BooleanTime(const double& JulianDate, double days_between, 
			    const double& start, const double& calculation_step_length);

void deleteOldOutputFiles(const std::string& outdir, const std::string& experiment, 
					 const std::string& station, const int& number_expo);

void qr_VersionUserRuntime(char *version, char *computation_date, double *jul_computation_date, 
					  char *user, mio::Date& date);

void qr_AverageFluxTimeSeries(const int& n_steps, const bool& useCanopyModel, 
						SN_SURFACE_DATA& Sdata, SN_STATION_DATA& Xdata);

void qr_TypeToCode(int *F1, int *F2, int *F3, int type);

bool qr_MassBalanceCheck(const SN_STATION_DATA& Xdata, const SN_SURFACE_DATA& Sdata, double& tot_mass_in);

double getModelledTemperature(const double& z, const SN_STATION_DATA& Xdata);

int findUpperNode(const double& z, const std::vector<SN_NODE_DATA>& Ndata, const int& nN);

double qro_ForcedErosion(const double hs1, SN_STATION_DATA *Xdata);

void qro_DeflateInflate(const SN_MET_DATA *Mdata, SN_STATION_DATA *Xdata, double *dhs_corr, double *mass_corr);

int ml_ag_Classify(const double& dendricity, const double& sphericity, 
			    const double& grain_dia, const int& marker, const double& theta_w, const double& theta_i);

#endif //End of Queries.h
