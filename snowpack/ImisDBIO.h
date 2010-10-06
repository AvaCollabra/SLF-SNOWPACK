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

#ifndef __IMISDBIO_H__
#define __IMISDBIO_H__

#include <snowpack/SnowpackIOInterface.h>
#include <snowpack/DataClasses.h>
#include <meteoio/MeteoIO.h>
#include <occi.h>

class ImisDBIO : public SnowpackIOInterface{

	public:
		ImisDBIO(const mio::Config& i_cfg);

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
		void parseStationName(const std::string& stationName, std::string& stName, std::string& stNumber);

		void deleteHdata(const std::string& stationName, const std::string& stationNumber,
					  const mio::Date& dateStart, const mio::Date& dateEnd, 
					  oracle::occi::Environment*& env, oracle::occi::Connection*& conn);

		void insertHdata(const std::string& stationName, const std::string& stationNumber,
					  const std::vector<Q_PROCESS_DAT>& Hdata, const int& num, 
					  oracle::occi::Environment*& env, oracle::occi::Connection*& conn);

		mio::Config cfg;
		std::string oracleDB, oraclePassword, oracleUser;

		static const double in_tz; //All IMIS data is in gmt+1
		static const std::string sqlDeleteHdata; //Delete statement for Hdata from snowpack.ams_pmod
		static const std::string sqlInsertHdata; //Insert statement for Hdata to snowpack.ams_pmod
};

#endif //End of ImisDBIO.h
