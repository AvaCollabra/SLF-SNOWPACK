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
 * @file Constants.h
 * @version 10.02
 * This module defines all the snow constants needed for the 1d snowpack model \n
 * Most external declarations are repeated and documented in Constants.c
 */

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

/// @brief Defines whether a distributed or 1D use of SNOWPACK is foreseen
#define ALPINE3D 0

/// @brief Default research variant
#define DEFLT 0
#define JAPAN 1
#define ANTARCTICA 2
#define CALIBRATION 3

#define VARIANT DEFLT //CALIBRATION //ANTARCTICA //JAPAN //

/**
 * @brief SVNREV is given as a compilation flag to tell us what is the version number \n
 * Please only use SN_VERSION in the code
 */
#ifndef SVNREV
#define SVNREV "20101006.001" //"20100616.837-845Mpc" //
#endif
#ifndef SN_VERSION
	#define SN_VERSION SVNREV
#endif

/*
 * INCLUDES
 */
#include <string.h>
#include <errno.h>

/*
 * CONSTANTS & ENUMERATIONS
 */
#define MAX_STRING_LENGTH 256
#define MAX_LINE_LENGTH 6000


/*****************************************************************************/
/*****************************************************************************/

/// @name Flags and numbers
//@{
#define ON 1
#define OFF 0
#define ERROR 0
#define NO_ERROR 1

#define NODATA -999.
#define INODATA 999
#define SNOWPACK_UNDEFINED -77777777
//@}

/**
 * @name INTERNAL MODEL PARAMETERS AND CONSTANTS
 * @brief These parameters should only be changed by advanced users.
 */
/// @brief Initial value for stability parameter
#define INIT_STABILITY 999.
/// @brief Switch on or off wind pumping in snow
#define WIND_PUMP ON
/// @brief Switch on or off wind pumping in soil
#define WIND_PUMP_SOIL ON

/// @brief Time step control parameters
typedef struct {
  double TimeN;      ///< Time of present calculation (s)
  int    nStep;      ///< Time step number
  double TimeEnd;    ///< Calculate up to TimeEnd that can be less than last time of last record (s)
  int    TsDump;     ///< Flag for time series dump
  int    nAvg;       ///< Number of calculation time steps to average fluxes etc.
  int    HzStep;     ///< Hazard step number (should be half of nStep in operational mode)
  int    HzDump;     ///< Calculation of hazard information will be performed
  int    PrDump;     ///< Flag for profile dump
  int    XdataDump;  ///< Backup of Xdata will be performed
  int    TaglayDump; ///< Flag for tagged layer series dump
  int    PrTabDump;  ///< Flag for tabular profile dump
} MN_TIME_DATA;

namespace Constants {

	const double pi = 3.14159265358979323846;
	const double g  = 9.80665;
	const double stefan_boltzmann = 5.67051e-8;

	const double solcon = 1366.1; ///< Total Solar Irradiance (W m-2) (Froehlich, 2006)

	const double gas_constant = 461.9; ///< (J mol-1 K-1)
	const double gas_constant_air = 287.0; ///< for air (J kg-1 K-1)
	const double gas_constant_mol = 8.31;  ///< (J mol-1 K-1)

	/// @name Density (kg m-3)
	//@{
	const double density_air = 1.1; ///< Approximation: use ideal gas law
	const double density_ice = 917.0; ///< At T = 0 degC
	const double density_water = 1000.0; ///<  At T = 0 degC
	const double density_vapor = 1.0; ///< Approximation: use ideal gas law
	/// @name Max and min densities (kg m-3)
	//@{
	const double min_rho = 5.0;
	const double max_rho = 1000.0;
	/// @brief Min volumetric ice content allowed
	const double min_hn_density = 30.0;
	const double min_ice_content = min_hn_density / density_ice;
	//@}
	//@}

	///@name Specific heat (J kg-1 K-1)
	//@{
	const double specific_heat_ice = 2100.0; ///< at T = 0 degC
	const double specific_heat_water = 4190.0; ///< at T = 0 degC
	const double specific_heat_air = 1004.67; ///< see Stull "Meteorology for scientists and engineers" p44
	//@}

	///@name Thermal conductivity of snow components ice, water and air (W m-1 K-1)
	//@{
	const double conductivity_ice = 2.2; ///< (W m-1 K-1)
	const double conductivity_water = 0.598; ///< (W m-1 K-1)
	const double conductivity_air = 0.026; ///< (W m-1 K-1)
	//@}

	///@name Vapor Diffusion in Air and Snow (m2 s-1)
	//@{
	///@brief This value was taken from: Colbeck, S.C., 1993. The vapor diffusion coefficient for snow, Water Resources Research, 29(1)
	const double diffusion_coefficient_in_air = 22.0e-6;  
	const double diffusion_coefficient_in_snow = 85.0e-6; ///< It is larger (see Colbeck) according to work prior to 2008!
	//@}

	///@name Phase change constants
	//@{
	const double melting_tk = 273.15; ///< (K)
	const double freezing_tk = 273.15; ///< (K)
	const double triple_point_t = 273.16; ///< (K)
	const double triple_point_p = 611.657; ///< (Pa)
	const double lh_sublimation = 2.838e6; ///< (J kg-1) (solid to vapor)
	const double lh_vaporization = 2.504e6; ///< (J kg-1) (liquid to vapor)
	const double lh_fusion = 3.34e5; ///< (J kg-1) (solid to liquid)
	//@}

	///@name Numerical Constants
	//@{
	/// @brief Small numbers for comparisons and to avoid divisions by zero
	const double eps = 1.e-6;
	const double eps2 = eps*eps;
	const double big = 1.0e200; ///< used for Dirichlet boundary conditions
	//@}
}


/**
 * @name MACRO definitions (to be replaced by proper functions some day)
 * TODO Move this definitions to Util.h, as functions to Util.c
 */
//@{
/// @name Elementary math functions
//@{
#define MAX(x,y)    (((x) < (y)) ? (y) : (x))
#define MIN(x,y)    (((x) < (y)) ? (x) : (y))
//@}
/// @name Conversion macros
//@{
/// @brief For temperatures
#define C_TO_K( T ) ( T + 273.15 )   // Celsius to Kelvin
#define K_TO_C( T ) ( T - 273.15 )   // Kelvin to Celsius
/// @brief For lengths
#define CM_TO_M( l ) ( (l) / 100. )  // centimeters to meters
#define M_TO_CM( l ) ( (l) * 100. )  // meters to centimeters
#define M_TO_MM( l ) ( (l) * 1000. ) // meters to millimeters
#define MM_TO_M( l ) ( (l) / 1000. ) // millimeters to meters
#define MM_TO_CM( l ) ( (l) / 10. )  // millimeters to meters
/// @brief For time
#define D_TO_H( t ) ( (t) * 24. )    // days to hours
#define D_TO_M( t ) ( (t) * 1440. )  // days to minutes
#define D_TO_S( t ) ( (t) * 86400. ) // days to seconds
#define H_TO_D( t ) ( (t) / 24.0 )   // hours to days
#define H_TO_M( t ) ( (t) * 60. )    // hours to minutes
#define H_TO_S( t ) ( (t) * 3600. )  // hours to seconds
#define M_TO_D( t ) ( (t) / 1440. )  // minutes to days
#define M_TO_H( t ) ( (t) / 60. )    // minutes to hours
#define M_TO_S( t ) ( (t) * 60. )    // minutes to seconds
#define S_TO_D( t ) ( (t) / 86400. ) // seconds to days
#define S_TO_H( t ) ( (t) / 3600. )  // seconds to hours
#define S_TO_M( t ) ( (t) / 60. )    // seconds to minutes
/// @brief For angles
#define DEG_TO_RAD(deg) ( deg *Constants::pi/180. )  // degree to radian
#define RAD_TO_DEG(rad) ( rad *180./Constants::pi )  // radian to degree
//@}
//@}

#endif
/*
 * End of Constants.h
*/
