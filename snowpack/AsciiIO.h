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

#ifndef __ASCIIIO_H__
#define __ASCIIIO_H__

#include <snowpack/Constants.h> 
#include <snowpack/SnowpackIOInterface.h> 
#include <snowpack/Hazard.h>
#include <meteoio/MeteoIO.h>
#include <snowpack/Canopy.h>

#include <sys/time.h> // time manipulation functions
#include <unistd.h>   // used for getlogin
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cstdarg> // needed for va_list
#include <time.h>
#include <string>

/// @brief Defines whether hardness (R) is output either in N (Swiss scale) or steps
#define R_IN_N 1

/// Defines whether surface temperature is included in comparison; If no values are provided in Master-file, they will be extrapolated
#define T_SRF 0
/// Defines whether snow/ground temperature is included in comparison; If no values are provided in Master-file, they will be extrapolated
#define T_GND 0

class AsciiIO : public SnowpackIOInterface {

	public:
		AsciiIO(const mio::Config& i_cfg);

		virtual void readSnowCover(const std::string& station, SN_SNOWSOIL_DATA& SSdata, SN_ZWISCHEN_DATA& Zdata);

		virtual void writeSnowCover(const mio::Date& date, const std::string& station, const SN_STATION_DATA& Xdata, 
							   const SN_ZWISCHEN_DATA& Zdata, const bool& forbackup=false);
		
		virtual void writeTimeSeries(const std::string& station, const SN_STATION_DATA& Xdata, 
							    const SN_SURFACE_DATA& Sdata, const SN_MET_DATA& Mdata, 
							    const Q_PROCESS_DAT& Hdata);
		
		virtual void writeProfile(const mio::Date& date, const std::string& station, const unsigned int& expo,
							 const SN_STATION_DATA& Xdata, const Q_PROCESS_DAT& Hdata);

		virtual void writeHazardData(const std::string& station, 
							    const std::vector<Q_PROCESS_DAT>& Hdata, const int& num);

	private:
		std::string getFilenamePrefix(const std::string& stationname);

		bool checkHeader(const char *fnam, const char *first_string, const Q_PROCESS_DAT* Hdata, 
					  const std::string& station, const char *ext, ...);

		void writeFreeProfilesDEFAULT(const SN_STATION_DATA& Xdata, FILE *fout);
		void writeFreeProfilesCALIBRATION(const SN_STATION_DATA& Xdata, FILE *fout);

		int writeTemperatures(FILE *fout, const double& z_vert, const double& T, 
						  const int& i, const SN_STATION_DATA& Xdata);

		double calcPerpPosition(const double& z_vert, const double& hs_ref, 
						    const double& ground, const double& angleSlope);
		double checkMeasuredTemperature(const double& T, const double& z, const double& mH);
		
		int findTaggedElement(const int& tag, const SN_STATION_DATA& Xdata);
		int writeHeightTemperatureTag(FILE *fout,const int& tag,const SN_MET_DATA& Mdata,const SN_STATION_DATA& Xdata);
		
		void writeFreeSeriesDEFAULT(const SN_STATION_DATA *Xdata, const SN_SURFACE_DATA *Sdata, 
							   const SN_MET_DATA *Mdata, const Q_PROCESS_DAT *Hdata, FILE *fout);
		void writeFreeSeriesANTARCTICA(const SN_STATION_DATA *Xdata, const SN_SURFACE_DATA *Sdata, 
								 const SN_MET_DATA *Mdata, const Q_PROCESS_DAT *Hdata, FILE *fout);
		void writeFreeSeriesCALIBRATION(const SN_STATION_DATA *Xdata, const SN_SURFACE_DATA *Sdata, 
								  const SN_MET_DATA *Mdata, const Q_PROCESS_DAT *Hdata, FILE *fout);

		void writeProfileDB(const mio::Date& date, const std::string& station, 
						const SN_STATION_DATA& Xdata, const Q_PROCESS_DAT& Hdata);

		mio::Config cfg;

		double calculation_step_length, hazard_steps_between;
		bool useCanopyModel, useSnowLayers, research_mode;
		int t_internal, sw_ref;
		int NUMBER_EXPO, NUMBER_SENSORS;//Actual number of "sensors" that are monitored, including tags in advanced mode
		bool OUT_HEAT, OUT_LW, OUT_SW, OUT_METEO, OUT_HAZ, OUT_MASS, OUT_T, OUT_LOAD, OUT_STAB, OUT_CANOPY;
		bool PERP_TO_SLOPE;

		std::string variant, experiment, outpath, inpath;

		//Defines heights of fixed sensors or/and initial depths of sensors with fixed settling rates
		std::vector<double> depth_of_sensors;

		int CHANGE_BC, MEAS_TSS;
		static const std::string profile_filename;
};

#endif //End of AsciiIO.h
