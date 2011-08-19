/***********************************************************************************/
/*  Copyright 2009-2010 WSL Institute for Snow and Avalanche Research    SLF-DAVOS */
/***********************************************************************************/
/* This file is part of MeteoIO.
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
#ifndef __MAINPAGE_H__
#define __MAINPAGE_H__

 /**
 * @mainpage Welcome to SNOWPACK
 * @section intro_sec Introduction
 * SNOWPACK is the operational model of the Swiss avalanche warning service and is available as an integrated software package, made of a library that can be reused
 * in third party applications and a standalone model. It simulates the evolution of the snow cover based on meteorological input data, based on the physical modeling
 * of the various processes taking place. It requires the following meteorological parameters:
 * - air temperature
 * - relative humidity
 * - wind speed
 * - incoming short wave radiation
 * - icoming long wave radiation OR surface temperature
 * - snow height OR precipitation
 * - ground temperature (if available)
 * - snow temperatures at various depth (if available and only for comparisons)
 *
 * These parameters MUST be available at least at a hourly time step. Please keep in mind that any inaccuracy on these parameters will have an impact on
 * the quality of the simulation. Since the modelling is physically-based, manually re-constructing missing data can lead to completely wrong results if not
 * done with great care (the model performing various checks on the physical consistency of the input data, it WILL exclude data points that are not consistent
 * with the other parameters. For example, precipitation occuring simultaneously with quite dry air will be refused).

 * International intercomparison studies show that SNOWPACK is successfully applied to alpine, arctic, maritime and continental snow covers.
 *
 * This library is available under LPGL version 3 or above, see <a href="http://www.gnu.org/licenses/lgpl.txt">www.gnu.org</a>.
 *
 */

#endif
