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

#ifndef __SMET_IO_H__
#define __SMET_IO_H__

#include <meteoio/MeteoIO.h>
#include <snowpack/Constants.h>
#include <snowpack/SnowpackIOInterface.h>
#include <snowpack/Hazard.h>
#include <snowpack/Canopy.h>

#include <string>

class SmetIO : public SnowpackIOInterface {

	public:
		SmetIO(const mio::Config& i_cfg);

		virtual void readSnowCover(const std::string& i_snowfile, const std::string& stationID,
		                           SN_SNOWSOIL_DATA& SSdata, SN_ZWISCHEN_DATA& Zdata);

		virtual void writeSnowCover(const mio::Date& date, const SnowStation& Xdata, const SN_SNOWSOIL_DATA& SSdata,
		                            const SN_ZWISCHEN_DATA& Zdata, const bool& forbackup=false);

		virtual void writeTimeSeries(const SnowStation& Xdata, const SurfaceFluxes& Sdata, const CurrentMeteo& Mdata,
		                             const ProcessDat& Hdata, const double wind_trans24);

		virtual void writeProfile(const mio::Date& date, SnowStation& Xdata, const ProcessDat& Hdata);

		virtual bool writeHazardData(const std::string& stationID, const std::vector<ProcessDat>& Hdata,
		                             const std::vector<ProcessInd>& Hdata_ind, const int& num);

	private:
		std::string getFilenamePrefix(const std::string& fnam, const std::string& path, const bool addexp=true);
		void setBasicHeader(const SnowStation& Xdata, const std::string& fields, smet::SMETWriter& smet_writer) const;
		void setSnoSmetHeader(const SnowStation& Xdata, const SN_SNOWSOIL_DATA& SSdata, const mio::Date& date,
		                      smet::SMETWriter& smet_writer) const;
		void setFormatting(const size_t& nr_solutes,
		                   std::vector<int>& vec_width, std::vector<int>&  vec_precision) const;
		void writeSnoFile(const std::string& filename, const mio::Date& date, const SnowStation& Xdata,
		                  const SN_SNOWSOIL_DATA& SSdata, const SN_ZWISCHEN_DATA& Zdata) const;
		void writeHazFile(const std::string& filename, const mio::Date& date,
		                  const SnowStation& Xdata, const SN_ZWISCHEN_DATA& Zdata) const;
		double get_doubleval(const smet::SMETReader& reader, const std::string& keyname) const;
		int get_intval(const smet::SMETReader& reader, const std::string& keyname) const;
		mio::Date read_snowsmet(const std::string& snofilename, const std::string& stationID, SN_SNOWSOIL_DATA& SSdata);
		mio::Date read_snosmet_header(const smet::SMETReader& sno_reader, const std::string& stationID,
		                              SN_SNOWSOIL_DATA& SSdata);
		mio::Date read_hazsmet(const std::string& hazfilename, SN_ZWISCHEN_DATA& Zdata);

	private:
		std::string outpath, o_snopath, snowpath, experiment, inpath, i_snopath;
		double time_zone;
		int sw_mode;
		bool perp_to_slope, useSoilLayers;
};

#endif //End of SmetIO.h
