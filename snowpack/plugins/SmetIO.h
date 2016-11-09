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

#ifndef SMET_IO_H
#define SMET_IO_H

#include <meteoio/MeteoIO.h>
#include <snowpack/Constants.h>
#include <snowpack/plugins/SnowpackIOInterface.h>
#include <snowpack/Hazard.h>

#include <string>

class SmetIO : public SnowpackIOInterface {

	public:
		SmetIO(const SnowpackConfig& i_cfg, const RunInfo& run_info);
		SmetIO(const SmetIO&);
		
		~SmetIO();
		
		SmetIO& operator=(const SmetIO&); ///<Assignement operator, required because of pointer member

		virtual bool snowCoverExists(const std::string& i_snowfile, const std::string& stationID) const;

		virtual void readSnowCover(const std::string& i_snowfile, const std::string& stationID,
		                           SN_SNOWSOIL_DATA& SSdata, ZwischenData& Zdata);

		virtual void writeSnowCover(const mio::Date& date, const SnowStation& Xdata,
		                            const ZwischenData& Zdata, const bool& forbackup=false);

		virtual void writeTimeSeries(const SnowStation& Xdata, const SurfaceFluxes& Sdata, const CurrentMeteo& Mdata,
		                             const ProcessDat& Hdata, const double wind_trans24);

		virtual void writeProfile(const mio::Date& date, const SnowStation& Xdata);

		virtual bool writeHazardData(const std::string& stationID, const std::vector<ProcessDat>& Hdata,
		                             const std::vector<ProcessInd>& Hdata_ind, const size_t& num);

		mio::Date read_hazsmet(const std::string& hazfilename, ZwischenData& Zdata);
		void writeHazFile(const std::string& hazfilename, const mio::Date& date,
		                         const SnowStation& Xdata, const ZwischenData& Zdata) const;

	private:
		std::string getFilenamePrefix(const std::string& fnam, const std::string& path, const bool addexp=true) const;
		void setBasicHeader(const SnowStation& Xdata, const std::string& fields, smet::SMETWriter& smet_writer) const;
		static void setSnoSmetHeader(const SnowStation& Xdata, const mio::Date& date, smet::SMETWriter& smet_writer);
		static void setFormatting(const size_t& nr_solutes,
		                   std::vector<int>& vec_width, std::vector<int>&  vec_precision);
		void writeSnoFile(const std::string& snofilename, const mio::Date& date, const SnowStation& Xdata,
		                         const ZwischenData& Zdata) const;
		bool keyExists(const smet::SMETReader& reader, const std::string& key) const;
		double get_doubleval(const smet::SMETReader& reader, const std::string& keyname) const;
		int get_intval(const smet::SMETReader& reader, const std::string& keyname) const;
		mio::Date read_snosmet(const std::string& snofilename, const std::string& stationID, SN_SNOWSOIL_DATA& SSdata);
		mio::Date read_snosmet_header(const smet::SMETReader& sno_reader, const std::string& stationID,
		                              SN_SNOWSOIL_DATA& SSdata);
		double compPerpPosition(const double& z_vert, const double& hs_ref,
		                        const double& ground, const double& cos_sl);
		std::string getFieldsHeader();
		void writeTimeSeriesHeader(const SnowStation& Xdata);
		void writeTimeSeriesData(const SnowStation& Xdata, const SurfaceFluxes& Sdata, const CurrentMeteo& Mdata, const ProcessDat& Hdata, const double &wind_trans24);

	private:
		std::vector<double> fixedPositions;
		std::string outpath, o_snowpath, snowpath, experiment, inpath, i_snowpath, sw_mode;
		const RunInfo info;
		smet::SMETWriter *ts_smet_writer;
		double in_dflt_TZ;
		double calculation_step_length, hazard_steps_between, ts_days_between;
		double min_depth_subsurf;
		bool avgsum_time_series, useCanopyModel, useSoilLayers, research_mode, perp_to_slope;
		bool out_heat, out_lw, out_sw, out_meteo, out_haz, out_mass, out_t, out_load, out_stab, out_canopy, out_soileb;
};

#endif
