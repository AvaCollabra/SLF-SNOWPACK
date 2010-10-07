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
 * @file Laws_sn.h
 * @version 10.02
 * @date 2009-10-21
 */

#ifndef __LAWS_SN_H__
#define __LAWS_SN_H__

#include <string>
#include <meteoio/MeteoIO.h>
#include <snowpack/DataClasses.h>
#include <snowpack/Snowpack.h> //some constants are necessary
/*
 * CONSTANTS & ENUMERATIONS
 */
/**
 * @name Albedo models
 * @brief Statistical models based on measured data
 * - ALB_LEHNING: data from the SLF study plot at Weissfluhjoch, 2540 m a.s.l.;
 *                *_2 corresponds to the latest regression. Variant ANTARCTICA sets coefficient Cage to 0.0
 * - ALB_NIED: Japanese adaptation of model ALB_LEHNING_2
 */
//@{
typedef enum {
 ALB_LEHNING_0,
 ALB_LEHNING_1,
 ALB_LEHNING_2,
 ALB_NIED,
 N_ALBM
} ALBM;
/// @brief Albedo model to be used
#if VARIANT == JAPAN
	#define ALBEDO_MODEL ALB_NIED
#else
	#define ALBEDO_MODEL ALB_LEHNING_2
#endif
//@}

/// @brief Defines whether the extinction coefficient for SW radiation in snow is based on 20 or less bands
#define BAND20 0

/**
 * @name Thermal conductivity
 * @brief Defines the constants and parameters for computing snow and soil thermal conductivity
 */
//@{
/// @brief Factor controlling ice to ice conduction
#define MONTANA_C_FUDGE 0.13
/// @brief Factor controlling increase in water vapor transport and thus energy transport in wet snow
#define MONTANA_VAPOR_FUDGE 2.5
//@}
/// @name Wind pumping (for soil see soil parameters)
//@{
/// @brief Used to decribe advective flux attenuation
#define WIND_EXT_COEF 0.1
/// @brief Integral of snow density corresponding to the maximal displacement depth d_pump (m)
#define DISPLACEMENT_COEF 0.7
/// @brief Ratio of Porosity to Tortuosity
#define ALPHA_POR_TOR 0.07
//@}

/// @name SOIL PARAMETERS
//@{
/**
 * @brief Define Method and Coefficents for the calculation of the influence of soil water
 * content on Evaporation from Bare Soil Layers:
 * - 1 ==> Resistance Approach, see Laws_sn.c
 * - 0 ==> Relative Humidity Approach, see Snowpack.c
 * - -1 ==> none, assume saturation pressure and no extra resistance
 */
#define SOIL_EVAPORATION 1
/// Minimum soil surface resistance, 50 sm-1 (van den Hurk et al, 2000)
#define RSOILMIN 50.0
/// Minimum relative saturation of top soil layer
#define RELSATMIN 0.05
/// @brief Ratio of Porosity to Tortuosity for Soil
#define ALPHA_POR_TOR_SOIL 0.05
/// @brief Pore length for surface soil for ventilation (m)
#define PORE_LENGTH_SOIL 0.01
/// Field capacity of the soil (1). Above this levels, water begins to drain
#define FIELD_CAPACITY_SOIL    0.15
//@}

/**
 * @name Snow viscosity
 * @brief Defines the constants and parameters for computing the snow viscosity \n
 * NOTE The Japanese variant uses a viscosity parametrization by Kojima
 * - T_TERM Temperature dependence
 * - V_TIME_FUDGE: Empirical constant related to age of snow
 * - V_ICE_FUDGE : Empirical constant related to volumetric ice content; ori: 0.4
 * - V_SP_FUDGE  : Empirical constant related to sphericity of snow grains
 * @brief Switches to test settling routine with lwsn_SnowViscosityCALIBRATION()
 * - VISC_CAL defines the calibration version to be used
 * - SETFIX is a quickfix for Antarctica, reducing settling to 6% of calculated value for snow older than 2 months
 */
//@{
typedef enum {
 VS_ANT=777,
 VS_STEINKOGLER=999,
 N_VISC_CAL
} VISC_CAL;
typedef enum {
 ARRHENIUS=1,
 ARRHENIUS_CRITICAL,
 STEINKOGLER,
 N_T_TERM
} T_TERM;
/// @brief T_term to be used
#if VARIANT == ANTARCTICA
	#define T_TERM ARRHENIUS_CRITICAL
	#define V_TIME_FUDGE 11.
	#define V_ICE_FUDGE 0.5
	#define V_SP_FUDGE 0.3
	#define VISC_CAL VS_ANT
	#define SETFIX 0
#elif VARIANT == CALIBRATION
	#define T_TERM ARRHENIUS_CRITICAL //STEINKOGLER //DEFAULT //
	#define V_TIME_FUDGE 11. //8. //
	#define V_ICE_FUDGE 0.5
	#define V_SP_FUDGE 0.3
	#define VISC_CAL VS_ANT //VS_STEINKOGLER //DEFAULT //
	#define SETFIX 0
#else
	#define T_TERM DEFLT
	#define V_TIME_FUDGE 11.
	#define V_ICE_FUDGE 0.5
	#define V_SP_FUDGE 0.3
	#define VISC_CAL DEFLT
	#define SETFIX 0
#endif
/// @brief Defines the smallest allowable viscosity (Pa s) that a viscosity law will return.
/// Value is DAMM SMALL -- smaller values than this are pretty unrealistic.
#define SMALLEST_VISCOSITY 1.0e6
//@}

enum {
	VS_KOJIMA=1,
	VS_CALIBRATION,
	N_VS_MODEL
};


/*
class LawsSnowpack {

	public:
		LawsSnowpack(const mio::Config& i_cfg);
*/
		double lwsn_InitialStress(const std::string& i_viscosity_model, const SN_ELEM_DATA& Edata);

		double lwsn_Albedo(const SN_ELEM_DATA& Edata, const double& Tss, const SN_MET_DATA& Mdata, const double& age);

		double lwsn_Extinction(const SN_ELEM_DATA& Edata);

		void lwsn_ShortWaveAbsorption(SN_STATION_DATA& Xdata, const double& I0, const int& SNP_SOIL, 
								const bool& multistream);

		double lwsn_SnowpackInternalEnergy(SN_STATION_DATA& Xdata);

		double lwsn_SoilFieldCapacity(const SN_ELEM_DATA *Edata);

		double lwsn_WindPumpingDisplacement(const SN_STATION_DATA& Xdata);

		double lwsn_WindPumpingVelocity(const SN_MET_DATA& Mdata, const double& d_pump);

		double lwsn_WindGradientSnow(double *v_pump, const SN_ELEM_DATA *Edata);

		double lwsn_SoilThermalConductivity(const SN_ELEM_DATA *Edata, const double dvdz);

		double lwsn_EnhanceWaterVaporTransportSnow(const SN_STATION_DATA& Xdata, const int& e);
		
		double lwsn_SnowThermalConductivity(const SN_ELEM_DATA& Edata, const double& dvdz);
		
		double lwsn_HeatCapacity(const SN_ELEM_DATA& Edata);
 
		double lwsn_SensibleHeat(const SN_MET_DATA& Mdata, const SN_STATION_DATA& Xdata, 
							const double& height_of_meteo_values);
 
		double lwsn_LatentHeat_Rh(const SN_MET_DATA& Mdata, const SN_STATION_DATA& Xdata, 
							 const double& height_of_meteo_values);

		double lwsn_LatentHeat(const SN_MET_DATA& Mdata, const SN_STATION_DATA& Xdata, 
						   const double& height_of_meteo_values);

		double lwsn_LongWave(const double T_snow, const double T_atm, const double e_atm);

		double lwsn_SnowElasticity(const double rho);

		double lwsn_NewSnowViscosityLehning(const SN_ELEM_DATA Edata);

		double lwsn_ConcaveNeckRadius(const double rg, const double rb);

		double lwsn_NeckLength(const double rg, const double rc);

		double lwsn_Neck2VolumetricStrain(const SN_ELEM_DATA& Edata);

		double lwsn_NeckStressEnhancement(const SN_ELEM_DATA& Edata);


		double lwsn_SnowViscosityTemperatureTerm(const double Te);

		double lwsn_SnowViscosity(const std::string& i_viscosity_model, 
							 const SN_ELEM_DATA& Edata, const mio::Date& date);
		double lwsn_SnowViscosityFudgeDEFAULT(const SN_ELEM_DATA& Edata, const mio::Date& date);
		double lwsn_SnowViscosityFudgeCALIBRATION(const SN_ELEM_DATA& Edata, const mio::Date& date);
		double lwsn_SnowViscosityDEFAULT(const SN_ELEM_DATA& Edata, const mio::Date& date);
		double lwsn_SnowViscosityKOJIMA(const SN_ELEM_DATA& Edata, const mio::Date& date);
		double lwsn_SnowViscosityCALIBRATION(const SN_ELEM_DATA& Edata, const mio::Date& date);
		
		double lwsn_InitialStressDEFAULT(const SN_ELEM_DATA& Edata);
		double lwsn_InitialStressCALIBRATION(const SN_ELEM_DATA& Edata);

/*
	private:
		mio::Config cfg;
		};
*/
#endif //End of Laws_sn.h
