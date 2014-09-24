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

#include <snowpack/DataClasses.h>
#include <meteoio/MeteoIO.h>

#include <snowpack/plugins/SnowpackIOInterface.h>

#include <snowpack/plugins/SmetIO.h>
#include <snowpack/plugins/AsciiIO.h>
#ifdef CAAMLIO
#include <snowpack/plugins/CaaMLIO.h>
#endif
#ifdef IMISDBIO
#include <snowpack/plugins/ImisDBIO.h>
#endif

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
