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

#ifndef __SNOWPACKIO_H__
#define __SNOWPACKIO_H__

#include <meteoio/MeteoIO.h>

#include <snowpack/DataClasses.h>
#include <snowpack/plugins/SnowpackIOInterface.h>
#include <snowpack/plugins/SmetIO.h>
#include <snowpack/plugins/AsciiIO.h>
#ifdef CAAMLIO
	#include <snowpack/plugins/CaaMLIO.h>
#endif
#ifdef IMISDBIO
	#include <snowpack/plugins/ImisDBIO.h>
#endif

/**
 * @page snowpackio Snowpack formats
 * Snowpack has the ability to read various format for its meteorological input data through the 
 * <a href="https://models.slf.ch/p/meteoio">MeteoIO</a> pre-processing library, so please check into
 * MeteoIO's documentation, in the <i>"Available plugins and usage"</i> section for the applicable formats.
 *
 * @section available_single_profile_plugins Available plugins for single snow profiles
 * The Snowpack specific data are supported directly in Snowpack and the formats listed in the table below 
 * are available, both for input and output of snow profiles with the <b>"SNOW"</b> keyword. 
 * Please read the documentation for each plugin in order to know the plugin-specific keywords!
 * <center><table border="1">
 * <tr><th>Key</th><th>Description</th><th>Extra requirements</th></tr>
 * <tr><td>\subpage ascii "SNOOLD"</td><td>legacy Snowpack profile</td><td></td></tr>
 * <tr><td>\subpage smet "SMET"</td><td>SMET based profile, recommended</td><td></td></tr>
 * <tr><td>\subpage caaml "CAAML"</td><td>CAAML profile</td><td><A HREF="http://xmlsoft.org/">libxml</A></td></tr>
 * </table></center>
 * 
 * @section available_profile_ts_plugins Available plugins for snow profiles time series
 * The Snowpack specific data are supported directly in Snowpack and the formats listed in the table below 
 * are available for output of snow profiles time series with the <b>"PROFILE_FORMAT"</b> keyword. 
 * Please read the documentation for each plugin in order to know the plugin-specific keywords!
 * <center><table border="1">
 * <tr><th>Key</th><th>Description</th><th>Extra requirements</th></tr>
 * <tr><td>\subpage ascii "PRO"</td><td>legacy Snowpack profile time series</td><td></td></tr>
 * <tr><td>\subpage ascii "PRF"</td><td>easier to parse profile time series</td><td></td></tr>
 * <tr><td>\subpage profile_imis "IMIS"</td><td>write profile time series to the IMIS database</td><td></td></tr>
 * </table></center>
 *
 */
class SnowpackIO : public SnowpackIOInterface {

	public:
		SnowpackIO(const SnowpackConfig& cfg);
		SnowpackIO(const SnowpackIO& source);
		~SnowpackIO();

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

		SnowpackIO& operator=(const SnowpackIO& source);

	private:
#ifdef IMISDBIO
		ImisDBIO *imisdbio;
#endif
#ifdef CAAMLIO
		CaaMLIO *caamlio;
#endif
		SmetIO *smetio;
		AsciiIO *asciiio;
		bool input_snow_as_smet, output_snow_as_smet;
		bool input_snow_as_caaml, output_snow_as_caaml;
		bool input_snow_as_ascii, output_snow_as_ascii;
		bool output_prf_as_ascii, output_prf_as_caaml, output_prf_as_imis;
		bool output_ts_as_ascii, output_haz_as_imis;
};

#endif //End of SnowpackIO.h
