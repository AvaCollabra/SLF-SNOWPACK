/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of Snowpack.
MeteoIO is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MeteoIO is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef __CAAMLIO_H__
#define __CAAMLIO_H__

#include <meteoio/MeteoIO.h>
#include <snowpack/Constants.h>
#include <snowpack/plugins/SnowpackIOInterface.h>
#include <snowpack/Hazard.h>
#include <snowpack/Canopy.h>
#include <snowpack/SmetIO.h>

#include <string>

#include <libxml/parser.h>
#include <libxml/xpath.h>

/**
 * @class CaaMLIO
 * @brief Reading snow profile data in CAAML format.
 * Reads in CAAML snow profile data, the CAA-IACS international standard
 *
 * @author Charles Fierz (Mathias Bavay)
 * @date   2014
 */
class CaaMLIO : public SnowpackIOInterface {

	public:
		CaaMLIO(const SnowpackConfig& i_cfg, const RunInfo& run_info);
		CaaMLIO(const CaaMLIO&);
		~CaaMLIO() throw();

		CaaMLIO& operator=(const CaaMLIO&); ///<Assignement operator, required because of pointer member

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

	private:
		void cleanup() throw();
		void init(const SnowpackConfig& cfg);
		void openIn_CAAML(const std::string& in_snowfile);
		void closeIn_CAAML() throw();
// 		bool parseStationData(const std::string& station_id, const xmlXPathContextPtr& xpathCtx, mio::StationData &sd);

// 		bool parseCaamlData(const mio::Date& dateStart, const mio::Date& dateEnd, const std::string& station_id,
// 		                    const mio::StationData& sd, const xmlXPathContextPtr& xpathCtx, std::vector<CaamlData> &vecCaaml) const;

		void setBasicHeader(const SnowStation& Xdata, const std::string& fields, smet::SMETWriter& smet_writer) const;
		void setSnoSmetHeader(const SnowStation& Xdata, const SN_SNOWSOIL_DATA& SSdata, const mio::Date& date, smet::SMETWriter& smet_writer) const;
		void setFormatting(const size_t& nr_solutes, std::vector<int>& vec_width, std::vector<int>&  vec_precision) const;

		std::string getFilenamePrefix(const std::string& fnam, const std::string& path, const bool addexp=true) const;
		bool read_snocaaml(const std::string& snofilename, const std::string& stationID, SN_SNOWSOIL_DATA& SSdata);
		void writeSnowFile(const std::string& snofilename, const mio::Date& date, const SnowStation& Xdata,
		                   const ZwischenData& Zdata) const;

		const RunInfo info;
		std::string i_snowpath, sw_mode, o_snowpath, experiment;
		bool useSoilLayers, perp_to_slope;
		/*static const*/ double in_tz; //plugin specific time zones
		std::string snow_prefix, snow_ext; //for the file naming scheme
		double caaml_nodata; //plugin specific no data value

		xmlDocPtr in_doc;
		xmlXPathContextPtr in_xpathCtx;
		xmlCharEncoding in_encoding;
// 		static const xmlChar* xml_attribute;
		static const xmlChar* xml_namespace;
		static const xmlChar* xml_namespace_abrev;
		static const std::string StationMetaData_xpath, SnowData_xpath;
};

#endif //End of CAAMLIO.h
