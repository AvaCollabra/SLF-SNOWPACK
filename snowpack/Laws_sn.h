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
#include <snowpack/Constants.h>
#include <snowpack/Laws.h>
#include <snowpack/Snowpack.h> //some constants are necessary
#include <snowpack/Metamorphism.h>

class SnLaws {

	public:
		/**
		 * @brief Switches to test settling routine with lwsn_SnowViscosityCALIBRATION()
		 * - VISC_CAL defines the calibration version to be used
		 * - SETFIX is a quickfix for Antarctica, reducing settling to 6% of calculated value for snow older than 2 months
		 */
		enum ViscosityCalVersion {
			vs_default=0,
			vs_ant=777,
			vs_steinkogler=999
		};

		enum TempDependence {
			default_term,
			arrhenius,
			arrhenius_critical,
			steinkogler,
			n_t_term
		};


		/**
		 * @name Albedo models
		 * @brief Statistical models based on measured data
		 * - ALB_LEHNING: data from the SLF study plot at Weissfluhjoch, 2540 m a.s.l.;
		 *                *_2 corresponds to the latest regression. Variant ANTARCTICA sets coefficient Cage to 0.0
		 * - ALB_NIED: Japanese adaptation of model ALB_LEHNING_2
		 */
		enum AlbedoModel {
			alb_lehning_0,
			alb_lehning_1,
			alb_lehning_2,
			alb_nied,
			n_albm
		};

	public:
		static double calcHeatCapacity(const SN_ELEM_DATA& Edata);
		static double calcNeckStressEnhancement(const SN_ELEM_DATA& Edata);
		static double calcConcaveNeckRadius(const double& rg, const double& rb);
		static double calcNeckLength(const double& rg, const double& rc);
		static double calcNeck2VolumetricStrain(const SN_ELEM_DATA& Edata);
		static double calcWindPumpingDisplacement(const SN_STATION_DATA& Xdata);

		static double calcSoilFieldCapacity(const SN_ELEM_DATA& Edata);
		static double calcSnowViscosityTemperatureTerm(const double& Te);
		static double calcSnowElasticity(const double& rho);

		static double calcLWRadCoefficient(const double& t_snow, const double& t_atm, const double& e_atm);
				
		static double calcSensibleHeat(const SN_MET_DATA& Mdata, const SN_STATION_DATA& Xdata, 
								 const double& height_of_meteo_values); 
		static double calcLatentHeat_Rh(const SN_MET_DATA& Mdata, const SN_STATION_DATA& Xdata, 
								  const double& height_of_meteo_values);
		static double calcLatentHeat(const SN_MET_DATA& Mdata, const SN_STATION_DATA& Xdata, 
							    const double& height_of_meteo_values);

		static double calcSoilThermalConductivity(const SN_ELEM_DATA& Edata, const double& dvdz);
		static double calcSnowThermalConductivity(const SN_ELEM_DATA& Edata, const double& dvdz);

		static double calcEnhanceWaterVaporTransportSnow(const SN_STATION_DATA& Xdata, const int& e);

		static double calcNewSnowViscosityLehning(const SN_ELEM_DATA& Edata);

		static double calcWindPumpingVelocity(const SN_MET_DATA& Mdata, const double& d_pump);
		static double calcWindGradientSnow(const SN_ELEM_DATA& Edata, double& v_pump);

		static double calcAlbedo(const std::string& variant, const SN_ELEM_DATA& Edata, const double& Tss, 
                                   const SN_MET_DATA& Mdata, const double& age);

		static double calcExtinction(const SN_ELEM_DATA& Edata);

		static void calcShortWaveAbsorption(const double& I0, const bool& useSnowLayers,
									 const bool& multistream,SN_STATION_DATA& Xdata);

		static double calcSnowpackInternalEnergy(SN_STATION_DATA& Xdata);

		static double calcInitialStress(const std::string& variant, const std::string& i_viscosity_model, 
								  const SN_ELEM_DATA& Edata);


		static double calcSnowViscosity(const std::string& variant, const std::string& i_viscosity_model, 
								  const SN_ELEM_DATA& Edata, const mio::Date& date);
		static double calcSnowViscosityFudgeDEFAULT(const SN_ELEM_DATA& Edata, const mio::Date& date);
		static double calcSnowViscosityFudgeCALIBRATION(const SN_ELEM_DATA& Edata, const mio::Date& date);
		static double calcSnowViscosityDEFAULT(const SN_ELEM_DATA& Edata, const mio::Date& date);
		static double calcSnowViscosityKOJIMA(const SN_ELEM_DATA& Edata, const mio::Date& date);
		static double calcSnowViscosityCALIBRATION(const SN_ELEM_DATA& Edata, const mio::Date& date);
		
		static const double smallest_viscosity;
		static const bool wind_pump, wind_pump_soil;

	private:
		static bool setStaticData(const std::string& variant);
		static const bool __init;
		static std::string current_variant;

		static AlbedoModel currentAlbedoModel;
		static double albedoCage;

		static ViscosityCalVersion visc_cal;
		static TempDependence t_term;
		static double v_time_fudge, v_ice_fudge, v_sp_fudge;
		static bool setfix;

		static const bool band20;
		static unsigned int swa_nb;
		static std::vector<double> swa_k, swa_pc, swa_fb;

		static const int soil_evaporation;
		static const double rsoilmin, relsatmin, alpha_por_tor_soil, pore_length_soil, field_capacity_soil;
		static const double montana_c_fudge, montana_vapor_fudge;
		static const double wind_ext_coef, displacement_coef, alpha_por_tor;
};

#endif //End of Laws_sn.h
