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
/**
 * @name USER SETTABLE CONSTANTS \n
 */
//@{
/// @brief String lengths
#define MAX_STRING_LENGTH 256
#define MAX_LINE_LENGTH 6000
//@}

/// @brief The maximum number of solutes to be treated
#define MAX_N_SOLUTES 10


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
//@{
/**
 * @brief Minimum element length (m)
 * - default: 0.0025
 * - Antarctica: 0.0001
 */
#if VARIANT == ANTARCTICA
	#define MINIMUM_L_ELEMENT 0.0001
#else
	#define MINIMUM_L_ELEMENT 0.0025
#endif

/// @brief Initial value for stability parameter
#define INIT_STABILITY 999.
/// @brief Switch on or off wind pumping in snow
#define WIND_PUMP ON
/// @brief Switch on or off wind pumping in soil
#define WIND_PUMP_SOIL ON

/**
 * @name Thresholds for surface hoar formation and burial
 * NOTE that the value of the parameter ROUGHNESS_LENGTH in CONSTANTS_User.INI is critical for surface hoar formation,
 * particularly for Dirichlet boundary conditions. Value should be < 1 mm. However, other considerations favor larger values.
 * - 0.0007 m : original calibration with the 98/99 data set \n
 * - 0.002  m : favored operational value with Dirichlet bc
 */
//@{
/// @brief Minimum size to show surface hoar on surface (mm)
#define MIN_SIZE_HOAR_SURF    0.5
/// @brief Density of surface hoar (-> hoar index of surface node) (kg m-3)
#define DENSITY_HOAR_SURF   100.
/// @brief Minimum surface hoar size to be buried (mm). Increased by 50% for Dirichlet bc.
#define MIN_SIZE_HOAR_BURIED  2.0
/**
 * @brief Density of BURIED surface hoar (kg m-3)
 * - default: 125.
 * - Antarctica: 200.
 */
#if VARIANT == ANTARCTICA
	#define DENSITY_HOAR_BURIED 200.
#else
	#define DENSITY_HOAR_BURIED 125.
#endif
//@}

/**
 * @name User parameters
 * @brief External declarations for user parameters set in CONSTANTS_User.INI (documented in Constants.c)
 */
//@{
//extern double INITIAL_HS;
//@}

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

  //Density (kg m-3)
  const double density_air = 1.1; ///< Approximation: use ideal gas law
  const double density_ice = 917.0; ///< At T = 0 degC
  const double density_water = 1000.0; ///<  At T = 0 degC
  const double density_vapor = 1.0; ///< Approximation: use ideal gas law

  const double solcon = 1366.1; ///< Total Solar Irradiance (W m-2) (Froehlich, 2006)

  const double gas_constant = 461.9; ///< (J mol-1 K-1)
  const double gas_constant_air = 287.0; ///< for air (J kg-1 K-1)
  const double gas_constant_mol = 8.31;  ///< (J mol-1 K-1)

  const double specific_heat_ice = 2100.0; ///< at T = 0 degC
  const double specific_heat_water = 4190.0; ///< at T = 0 degC
  const double specific_heat_air = 1004.67; ///< see Stull "Meteorology for scientists and engineers" p44
}

/// @name PHYSICAL CONSTANTS
//@{
/// @name Constants
//@{
/// @brief 
//#define SOLCON  1366.1
/// @brief (m s-2)
//#define GRAVITY	9.80665
/// @brief (W m-2 K-4)
//#define STEFAN_BOLTZMANN 5.67051e-8
/// @brief For air (J kg-1 K-1)
//#define GAS_CONSTANT_AIR 287.
/// @brief For water vapor (J kg-1 K-1)
//#define GAS_CONSTANT 461.9
/// @brief (J mol-1 K-1)
//#define GAS_CONSTANT_MOL 8.31
//@}
/// @name Density (kg m-3)
//@{
/// @brief  At T = 0 degC
//#define DENSITY_ICE 917.0
/// @brief At T = 0 degC
//#define DENSITY_WATER 1000.0
/// @brief Approximation: use ideal gas law
//#define DENSITY_VAPOR 1.0
/// @brief Approximation: use ideal gas law
//#define DENSITY_AIR 1.1
/// @name Max and min densities (kg m-3)
//@{
#define MIN_RHO 5.0
#define MAX_RHO 1000.
/// @brief Min volumetric ice content allowed
#define MIN_HN_DENSITY     30.
#define MIN_ICE_CONTENT MIN_HN_DENSITY/Constants::density_ice
//@}
//@}
/// @name Specific heat (J kg-1 K-1)
//@{
/// @brief At T = 0 degC
//#define SPECIFIC_HEAT_ICE 2100.0
/// @brief At T = 0 degC
//#define SPECIFIC_HEAT_WATER 4190.0
/// @brief See Stull "Meteorology for scientists and engineers" p44
//#define SPECIFIC_HEAT_AIR 1004.67
//@}
/// @name Thermal conductivity of snow components ice, water and air (W m-1 K-1)
//@{
#define CONDUCTIVITY_ICE 2.2
#define CONDUCTIVITY_WATER 0.598
#define CONDUCTIVITY_AIR 0.026
//@}
/// @brief Vapor Diffusion in Air and Snow (m2 s-1)
//@{
/// @brief This value was taken from: Colbeck, S.C., 1993. The vapor diffusion coefficient for snow, Water Resources Research, 29(1)
#define	DIFFUSION_COEFFICIENT_IN_AIR  22.0e-6
/// @brief It is larger (see Colbeck) according to work prior to 2008!!!
#define	DIFFUSION_COEFFICIENT_IN_SNOW 85.0e-6
//@}
/// @name Phasechange constants
//@{
/// @brief (degC)
#define MELTING_TC 0.0
/// @brief (K)
#define MELTING_TK 273.15
/// @brief (K)
#define FREEZING_TK 273.15
/// @brief (K)
#define TRIPLE_POINT_T 273.16
/// @brief (Pa)
#define TRIPLE_POINT_P 610.5
/// @brief (J kg-1)	( Solid to Vapor )
#define LH_SUBLIMATION 2.838e6
/// @brief (J kg-1)	( Liquid to Vapor )
#define LH_VAPORIZATION 2.504e6
/// @brief (J kg-1)	( Solid to Liquid )
#define LH_FUSION 3.34e5
//@}
/// @name Numerical Constants
//@{
/// @brief Small numbers for comparisons and to avoid divisions by zero
#define EPS 1.e-6
#define EPS2 EPS*EPS
/// @brief Used for Dirichlet boundary conditions
#define BIG 1.0e200
//#define PI  3.141593
//@}
//@}

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
