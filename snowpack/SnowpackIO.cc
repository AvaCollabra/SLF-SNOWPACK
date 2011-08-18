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

#include <snowpack/SnowpackIO.h>

using namespace std;
using namespace mio;

#ifdef IMISDBIO
SnowpackIO::SnowpackIO(const mio::Config& i_cfg) : cfg(i_cfg), asciiio(cfg), smetio(cfg), imisdbio(cfg)
{
#else
SnowpackIO::SnowpackIO(const mio::Config& i_cfg) : cfg(i_cfg), asciiio(cfg), smetio(cfg)
{
#endif
	//Actual constructor code following
	outputprofile_as_ascii = false;
	outputprofile_as_imis  = false;

	//The profiles may be dumped either in ASCII format or in another ASCII format for upload to the DB
	//The user can switch the desired mode on by specifying "ASCII" or "IMIS" or both in the io.ini
	vector<string> vecProfileOutput = cfg.get("PROFILE", "Output", Config::nothrow);
	if (vecProfileOutput.size() == 0) {
		outputprofile_as_ascii = true;
		outputprofile_as_imis  = false;
	} else if (vecProfileOutput.size() > 2) {
		throw InvalidArgumentException("The key PROFILE in section OUTPUT can have two values at most", AT);
	} else {
		for (unsigned int ii=0; ii<vecProfileOutput.size(); ii++){
			if (vecProfileOutput[ii] == "ASCII"){
				outputprofile_as_ascii = true;
			} else if (vecProfileOutput[ii] == "IMIS") {
				outputprofile_as_imis  = true;
			} else {
				throw InvalidArgumentException("Key PROFILE / section OUTPUT: only values ASCII or IMIS expected", AT);
			}
		}
	}

	//Initial snow profile, where to read it from:
	string in_snow;
	cfg.getValue("SNOW", "Input", in_snow, Config::nothrow);

	if (in_snow == "SMET"){
		input_snow_as_smet = true;
	} else {
		input_snow_as_smet = false;
	}
}

void SnowpackIO::readSnowCover(const std::string& i_snowfile, const std::string& stationID,
                               SN_SNOWSOIL_DATA& SSdata, SN_ZWISCHEN_DATA& Zdata)
{
	if (input_snow_as_smet){
		smetio.readSnowCover(i_snowfile, stationID, SSdata, Zdata);
	} else {
		asciiio.readSnowCover(i_snowfile, stationID, SSdata, Zdata);
	}
}

void SnowpackIO::writeSnowCover(const mio::Date& date, const SnowStation& Xdata, const SN_SNOWSOIL_DATA& SSdata,
                                const SN_ZWISCHEN_DATA& Zdata, const bool& forbackup)
{
	asciiio.writeSnowCover(date, Xdata, SSdata, Zdata, forbackup);
	smetio.writeSnowCover(date, Xdata, SSdata, Zdata, forbackup);
}

void SnowpackIO::writeTimeSeries(const SnowStation& Xdata, const SurfaceFluxes& Sdata, const CurrentMeteo& Mdata,
                                 const ProcessDat& Hdata, const double wind_trans24)
{
	asciiio.writeTimeSeries(Xdata, Sdata, Mdata, Hdata, wind_trans24);
}

void SnowpackIO::writeProfile(const mio::Date& date, SnowStation& Xdata, const ProcessDat& Hdata)
{
	if (outputprofile_as_ascii)
		asciiio.writeProfile(date, Xdata, Hdata);

	if (outputprofile_as_imis){
#ifdef IMISDBIO
		imisdbio.writeProfile(date, Xdata, Hdata);
#endif
	}
}

#ifdef IMISDBIO
bool SnowpackIO::writeHazardData(const std::string& stationID, const std::vector<ProcessDat>& Hdata,
                                 const std::vector<ProcessInd>& Hdata_ind, const int& num)
{
	if (imisdbio.writeHazardData(stationID, Hdata, Hdata_ind, num)){
		return true;
	} else {
		return false;
	}
}
#else
bool SnowpackIO::writeHazardData(const std::string& /*stationID*/, const std::vector<ProcessDat>& /*Hdata*/,
                                 const std::vector<ProcessInd>& /*Hdata_ind*/, const int& /*num*/)
{
	return false;
}
#endif
