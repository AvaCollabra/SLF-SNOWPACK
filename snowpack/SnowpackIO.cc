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

#include <snowpack/SnowpackIO.h>

using namespace std;
using namespace mio;

#ifdef IMISDBIO
SnowpackIO::SnowpackIO(const SnowpackConfig& cfg) : imisdbio(NULL),
#else
SnowpackIO::SnowpackIO(const SnowpackConfig& cfg) :
#endif
              asciiio(NULL), smetio(NULL),
              input_snow_as_smet(false), output_snow_as_smet(false),
              output_prf_as_ascii(false), output_ts_as_ascii(false),
              output_prf_as_imis(false), output_haz_as_imis(false)
              
{
	//Format of initial snow profile:
	const string in_snow = cfg.get("SNOW", "Input", IOUtils::nothrow);
	if (in_snow == "SMET") { //TODO: document Input::SNOW = SMET or SNOOLD
		input_snow_as_smet = true;
	}
	//Format of transitional and final snow profile(s):
	const string out_snow = cfg.get("SNOW", "Output", IOUtils::nothrow);
	if (out_snow == "SMET") { //TODO: document ouput::SNOW = SMET or SNOOLD
		output_snow_as_smet = true;
	}
	/* Profiles may also be dumped in up to 4 formats specified by the key PROF_FMT in section [Output]:
	 * VISU     : ASCII-format for visulization with SN_GUI
	 * FULL_PRF : full profiles in tabular ASCII-format
	 * AGGR_PRF : aggregated profiles in tabular ASCII-format
	 * IMIS     : aggregated profiles for upload to the SLF database sdbo
	 */
	std::vector<string> vecProfileFmt = cfg.get("PROF_FMT", "Output", IOUtils::nothrow);
	if (vecProfileFmt.size() > 4) {
		throw InvalidArgumentException("The key PROF_FMT in section [Output] can have four values at most", AT);
	} else {
		for (size_t ii=0; ii<vecProfileFmt.size(); ii++) {
			if (vecProfileFmt[ii] == "VISU") {
				output_prf_as_ascii = true;
			} else if (vecProfileFmt[ii] == "FULL_PRF") {
				output_prf_as_ascii  = true;
			} else if (vecProfileFmt[ii] == "AGGR_PRF") {
				output_prf_as_ascii  = true;
			} else if (vecProfileFmt[ii] == "IMIS") {
				output_prf_as_imis  = true;
			} else {
				throw InvalidArgumentException("Key PROF_FMT in section [Output] takes only VISU, FULL_PRF, AGGR_PRF or IMIS values", AT);
			}
		}
	}
	//Format of meteo time series:
    output_ts_as_ascii = cfg.get("TS_WRITE", "Output", IOUtils::nothrow);

//set the "plugins" pointers
	RunInfo run_info;
	if(input_snow_as_smet || output_snow_as_smet) smetio = new SmetIO(cfg, run_info);
	if(output_prf_as_ascii || output_ts_as_ascii) asciiio = new AsciiIO(cfg, run_info);
#ifdef IMISDBIO
	output_haz_as_imis = output_prf_as_imis;
	if(output_prf_as_imis || output_haz_as_imis) imisdbio = new ImisDBIO(cfg, run_info);
#endif
}

SnowpackIO::SnowpackIO(const SnowpackIO& source) :
#ifdef IMISDBIO
          imisdbio(source.imisdbio),
#endif
          asciiio(source.asciiio), smetio(source.smetio),
          input_snow_as_smet(source.input_snow_as_smet), output_snow_as_smet(source.output_snow_as_smet),
          output_prf_as_ascii(source.output_prf_as_ascii), output_ts_as_ascii(source.output_ts_as_ascii),
          output_prf_as_imis(source.output_prf_as_imis), output_haz_as_imis(source.output_haz_as_imis)
{}

SnowpackIO::~SnowpackIO()
{
	if(smetio != NULL) delete smetio;
	if(asciiio != NULL) delete asciiio;
#ifdef IMISDBIO
	if(imisdbio != NULL) delete imisdbio;
#endif
}

bool SnowpackIO::snowCoverExists(const std::string& i_snowfile, const std::string& stationID) const
{
	if (input_snow_as_smet){
		return smetio->snowCoverExists(i_snowfile, stationID);
	} else {
		return asciiio->snowCoverExists(i_snowfile, stationID);
	}
}

void SnowpackIO::readSnowCover(const std::string& i_snowfile, const std::string& stationID,
                               SN_SNOWSOIL_DATA& SSdata, ZwischenData& Zdata)
{
	if (input_snow_as_smet){
		smetio->readSnowCover(i_snowfile, stationID, SSdata, Zdata);
	} else {
		asciiio->readSnowCover(i_snowfile, stationID, SSdata, Zdata);
	}
}

void SnowpackIO::writeSnowCover(const mio::Date& date, const SnowStation& Xdata, const SN_SNOWSOIL_DATA& SSdata,
                                const ZwischenData& Zdata, const bool& forbackup)
{
	if (output_snow_as_smet){
		smetio->writeSnowCover(date, Xdata, SSdata, Zdata, forbackup);
	} else {
		asciiio->writeSnowCover(date, Xdata, SSdata, Zdata, forbackup);
	}
}

void SnowpackIO::writeTimeSeries(const SnowStation& Xdata, const SurfaceFluxes& Sdata, const CurrentMeteo& Mdata,
                                 const ProcessDat& Hdata, const double wind_trans24)
{
	if(output_ts_as_ascii)
		asciiio->writeTimeSeries(Xdata, Sdata, Mdata, Hdata, wind_trans24);
}

void SnowpackIO::writeProfile(const mio::Date& date, const SnowStation& Xdata)
{
	if (output_prf_as_ascii)
		asciiio->writeProfile(date, Xdata);

#ifdef IMISDBIO
	if (output_prf_as_imis)
		imisdbio->writeProfile(date, Xdata);
#endif
}

#ifdef IMISDBIO
bool SnowpackIO::writeHazardData(const std::string& stationID, const std::vector<ProcessDat>& Hdata,
                                 const std::vector<ProcessInd>& Hdata_ind, const int& num)
{
	if(output_haz_as_imis)
		return imisdbio->writeHazardData(stationID, Hdata, Hdata_ind, num);
	return false;
}
#else
bool SnowpackIO::writeHazardData(const std::string& /*stationID*/, const std::vector<ProcessDat>& /*Hdata*/,
                                 const std::vector<ProcessInd>& /*Hdata_ind*/, const int& /*num*/)
{
	return false;
}
#endif

SnowpackIO& SnowpackIO::operator=(const SnowpackIO& source)
{
	if(this != &source) {
#ifdef IMISDBIO
		imisdbio = source.imisdbio;
#endif
		asciiio = source.asciiio;
		smetio = source.smetio;
		output_prf_as_ascii = source.output_prf_as_ascii;
		output_prf_as_imis = source.output_prf_as_imis;
		output_snow_as_smet = source.output_snow_as_smet;
		input_snow_as_smet = source.input_snow_as_smet;
		output_ts_as_ascii = source.output_ts_as_ascii;
		output_haz_as_imis = source.output_haz_as_imis;
	}
	return *this;
}
