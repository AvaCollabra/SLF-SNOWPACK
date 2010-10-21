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
SnowpackIO::SnowpackIO(const mio::Config& i_cfg) : cfg(i_cfg), asciiio(cfg), imisdbio(cfg)
{
#else
SnowpackIO::SnowpackIO(const mio::Config& i_cfg) : cfg(i_cfg), asciiio(cfg)
{
#endif
	//Actual constructor code following
	outputprofile_as_ascii = false;
	outputprofile_as_imis  = false;

	//The profiles may be dumped either in ASCII format or in another ASCII format for upload to the DB
	//The user can switch the desired mode on by specifying "ASCII" or "IMIS" or both in the io.ini
	vector<string> vecProfileOutput = cfg.get("PROFILE", "Output", Config::nothrow);
	if (vecProfileOutput.size() == 0){
		outputprofile_as_ascii = true;
		outputprofile_as_imis  = false;
	} else if (vecProfileOutput.size() > 2){
		throw InvalidArgumentException("The key PROFILE in section OUTPUT can have two values at most", AT);
	} else {
		for (unsigned int ii=0; ii<vecProfileOutput.size(); ii++){
			if (vecProfileOutput[ii] == "ASCII"){
				outputprofile_as_ascii = true;
			} else if (vecProfileOutput[ii] == "IMIS"){
				outputprofile_as_imis  = true;
			} else {
				throw InvalidArgumentException("Key PROFILE / section OUTPUT: only values ASCII or IMIS expected", AT);
			}
		}
	}
}

void SnowpackIO::readSnowCover(const std::string& station, SN_SNOWSOIL_DATA& SSdata, SN_ZWISCHEN_DATA& Zdata)
{
	asciiio.readSnowCover(station, SSdata, Zdata);
}

void SnowpackIO::writeSnowCover(const mio::Date& date, const std::string& station, const SN_STATION_DATA& Xdata, 
				const SN_ZWISCHEN_DATA& Zdata, const bool& forbackup)
{
	asciiio.writeSnowCover(date, station, Xdata, Zdata, forbackup);
}
	
void SnowpackIO::writeTimeSeries(const std::string& station, const SN_STATION_DATA& Xdata, 
				 const SN_SURFACE_DATA& Sdata, const SN_MET_DATA& Mdata, const Q_PROCESS_DAT& Hdata)
{
	asciiio.writeTimeSeries(station, Xdata, Sdata, Mdata, Hdata);
}
	
void SnowpackIO::writeProfile(const mio::Date& date, const std::string& station, const unsigned int& expo,
						const SN_STATION_DATA& Xdata, const Q_PROCESS_DAT& Hdata)
{
	if (outputprofile_as_ascii)
		asciiio.writeProfile(date, station, expo, Xdata, Hdata);

	if (outputprofile_as_imis){
#ifdef IMISDBIO
		imisdbio.writeProfile(date, station, expo, Xdata, Hdata);
#endif
	}
}

void SnowpackIO::writeHazardData(const std::string& station, const std::vector<Q_PROCESS_DAT>& Hdata, 
						   const std::vector<Q_PROCESS_IND>& Hdata_ind, const int& num)
{
#ifdef IMISDBIO
	imisdbio.writeHazardData(station, Hdata, Hdata_ind, num);
#endif
}
