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
#include <ctime>
#include <string>
#include <set>

class AsciiIO : public SnowpackIOInterface {

	public:
		AsciiIO(const mio::Config& i_cfg);

		virtual void readSnowCover(const std::string& i_snowfile, const std::string& stationID,
		                           SN_SNOWSOIL_DATA& SSdata, SN_ZWISCHEN_DATA& Zdata);

		virtual void writeSnowCover(const mio::Date& date, const SnowStation& Xdata, const SN_SNOWSOIL_DATA& SSdata,
                                const SN_ZWISCHEN_DATA& Zdata, const bool& forbackup=false);
		
		virtual void writeTimeSeries(const SnowStation& Xdata, const SurfaceFluxes& Sdata, const CurrentMeteo& Mdata,
		                             const ProcessDat& Hdata, const double wind_trans24);
		
		virtual void writeProfile(const mio::Date& date, SnowStation& Xdata, const ProcessDat& Hdata);

		virtual bool writeHazardData(const std::string& station, const std::vector<ProcessDat>& Hdata,
                                 const std::vector<ProcessInd>& Hdata_ind, const int& num);

	private:
		bool appendFile(const std::string& filename, const mio::Date& startdate, const std::string& ftype);
		bool parseMetFile(const char& eoln, const mio::Date& start_date, std::istream& fin, std::ostream& ftmp);
		bool parseProFile(const char& eoln, const mio::Date& start_date, std::istream& fin, std::ostream& ftmp);

		std::string getFilenamePrefix(const std::string& fnam, const std::string& path, const bool addexp=true);

		bool checkHeader(const char *fnam, const char *first_string, const ProcessDat& Hdata, const char *ext, ...);

		void writeFreeProfileDEFAULT(SnowStation& Xdata, FILE *fout);
		void writeFreeProfileCALIBRATION(SnowStation& Xdata, FILE *fout);

		int writeTemperatures(FILE *fout, const double& z_vert, const double& T,
                          const int& i, const SnowStation& Xdata);

		double compPerpPosition(const double& z_vert, const double& hs_ref,
                            const double& ground, const double& slope_angle);
		double checkMeasuredTemperature(const double& T, const double& z, const double& mH);
		
		int findTaggedElement(const int& tag, const SnowStation& Xdata);
		int writeHeightTemperatureTag(FILE *fout,const int& tag,const CurrentMeteo& Mdata,const SnowStation& Xdata);
		
		void writeFreeSeriesDEFAULT(const SnowStation& Xdata, const SurfaceFluxes& Sdata,
                                const CurrentMeteo& Mdata, const double crust, const unsigned int nCalcSteps, FILE *fout);
		void writeFreeSeriesANTARCTICA(const SnowStation& Xdata, const SurfaceFluxes& Sdata,
                                   const CurrentMeteo& Mdata, const double crust, const unsigned int nCalcSteps,
                                   FILE *fout);
		void writeFreeSeriesCALIBRATION(const SnowStation& Xdata, const SurfaceFluxes& Sdata,
                                    const CurrentMeteo& Mdata, const double crust, const unsigned int nCalcSteps,
                                    FILE *fout);

		static const bool r_in_n, t_srf, t_gnd;
		
		const mio::Config& cfg;
		double time_zone; ///< input data time zone

		double calculation_step_length, hazard_steps_between, ts_days_between;
		double avgsum_time_series;
		bool useCanopyModel, useSoilLayers, research_mode;
		int sw_mode;
		bool out_heat, out_lw, out_sw, out_meteo, out_haz, out_mass, out_t, out_load, out_stab, out_canopy;
		bool perp_to_slope;
		double min_depth_subsurf;
		int number_meas_temperatures, number_fixed_heights, number_fixed_rates, max_number_sensors;
		int number_sensors; //Actual number of "sensors" that are monitored, including tags in advanced mode
		double hoar_density_surf, hoar_min_size_surf;
		//int number_of_solutes;

		std::string variant, experiment, station_name;
		std::string inpath, snowfile, i_snopath, outpath, o_snopath;

		//Defines heights of fixed sensors or/and initial depths of sensors with fixed settling rates
		std::vector<double> fixed_sensor_depths;

		int CHANGE_BC, MEAS_TSS;

		std::set<std::string> setAppendableFiles;

		//std::ofstream fout;//Output file streams
};

#endif //End of AsciiIO.h
