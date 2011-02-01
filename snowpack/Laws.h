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
/**
 * @file Laws.h
 * @version 10.02
 */

#ifndef __LAWS_H__
#define __LAWS_H__

#include <math.h>

double lw_emissivity(const double lwr, const double T);

double lw_lwr(const double ea, const double T);

double lw_TairLapseRate(const double ta, const double ref_alti, const double altitude);

double lw_AirPressure(const double altitude);

double lw_WetBulbTemperature(const double L, const double T, const double RH, const double altitude);

double lw_SaturationPressure(const double T);

double RhtoDewPoint(double RH, double TA, const bool& force_water);

double DewPointtoRh(double TD, double TA, const bool& force_water);

double lw_LW_Brutsaert(const double e0, const double ta);

double lw_Omstedt(const double e0, const double cloud_frac);

/* minimum observed air emissivity
 * - default: 0.55 (from 1993 data at Weissfluhjoch)
 * - Antarctica: 0.31 (from 2006/2007 data of Dome C) */
double lw_AirEmissivity(const double input, const double ta, const double rh, const double min_air_emissivity=0.55);

double lw_ArrheniusLaw(const double ActEnergy, const double T, const double T_ref);

#endif //End of Laws.h
