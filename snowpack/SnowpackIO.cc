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

SnowpackIO::SnowpackIO(const mio::Config& i_cfg) : cfg(i_cfg), asciiio(cfg) {}

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
	
void SnowpackIO::writeProfile(const mio::Date& date, const std::string& station,
						const SN_STATION_DATA& Xdata, const Q_PROCESS_DAT& Hdata)
{
	asciiio.writeProfile(date, station, Xdata, Hdata);
}
