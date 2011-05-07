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

		//@{
		/// Types of events for computing new snow density
		enum EventType {
			event_default,
			event_wind,   ///< Wind driven deposition of snow (Antarctica)
			nEventType	
		};
		/// Defines which version to use while calibrating snow viscosity
		enum ViscosityCalVersion {
			visc_cal_new=111,         ///< version currently under test
			visc_cal_897=333,         ///< calibration fall 2010 by Fierz
			visc_cal_ant=911,         ///< calibration fall 2010 adapted to Antarctica
			visc_cal_steinkogler=933, ///< calibration 2009 by Walter Steinkogler (MSc thesis)
			visc_cal_837=999,         ///< as of revision 837 (deprecated)
			nViscosityCalVersion
		};
		/// Defines temperature dependence of snow viscosity
		enum TempDependence {
			t_term_arrhenius_critical, ///< Arrhenius type multiplied by critical function near melting point
			t_term_arrhenius,          ///< pure Arrhenius type w/ excitation energy of ice
			t_term_steinkogler=933,    ///< calibration 2009 by Walter Steinkogler (MSc thesis)
			t_term_837=999,            ///< as of revision 837 (deprecated)
			nTempDependence
		};
		/// Defines which snow albedo parameterization to use
		enum AlbedoModel {
			alb_lehning_0,
			alb_lehning_1,
			alb_lehning_2,
			alb_nied,      ///< Japanese version of alb_lehning_2
			nAlbedoModel
		};
		//@}
		
		static double conductivity_ice(const double& Temperature);
		static double conductivity_water(const double& Temperature);
		static double conductivity_air(void);


		static double compWindPumpingDisplacement(const SnowStation& Xdata);
		static double compWindPumpingVelocity(const CurrentMeteo& Mdata, const double& d_pump);
		static double compWindGradientSnow(const ElementData& Edata, double& v_pump);

		static double compSensibleHeatCoefficient(const CurrentMeteo& Mdata, const SnowStation& Xdata,
		                                          const double& height_of_meteo_values);
		static double compLatentHeat_Rh(const CurrentMeteo& Mdata, SnowStation& Xdata,
		                                const double& height_of_meteo_values);
		static double compLatentHeat(const CurrentMeteo& Mdata, SnowStation& Xdata,
		                             const double& height_of_meteo_values);

		static double compSoilThermalConductivity(const ElementData& Edata, const double& dvdz);
		static double compSnowThermalConductivity(const ElementData& Edata, const double& dvdz);

		static double compEnhanceWaterVaporTransportSnow(const SnowStation& Xdata, const int& e);

		static double compLWRadCoefficient(const double& t_snow, const double& t_atm, const double& e_atm);

		static double compSnowAlbedo(const std::string& variant, const double& i_fixed_albedo,
		                             const ElementData& Edata, const double& Tss, const CurrentMeteo& Mdata,
		                             const double& age);
		static void compShortWaveAbsorption(const double& I0, const bool& useSoilLayers,
		                                    const bool& multistream,SnowStation& Xdata);

		static double compNewSnowDensity(const std::string& i_hn_density_model, const double& hn_fixed_density,
		                                 const CurrentMeteo& Mdata, const SnowStation& Xdata,
		                                 const double& tss, const double& hnw);
		static double NewSnowViscosityLehning(const ElementData& Edata);

		static double compInitialStress(const std::string& variant, ElementData& Edata, const mio::Date& date);
		static double initialStressDEFAULT(ElementData& Edata, const mio::Date& date);
		static double initialStressCALIBRATION(ElementData& Edata, const mio::Date& date);
		static double snowViscosityFudgeDEFAULT(const ElementData& Edata);
		static double snowViscosityFudgeCALIBRATION(const ElementData& Edata, const mio::Date& date);
		static double snowViscosityTemperatureTerm(const double& Te);
		static double compSnowViscosity(const std::string& variant, const std::string& i_viscosity_model,
		                                ElementData& Edata, const mio::Date& date);
		static double snowViscosityDEFAULT(ElementData& Edata);
		static double snowViscosityKOJIMA(const ElementData& Edata);
		static double snowViscosityCALIBRATION(ElementData& Edata, const mio::Date& date);

		static double min_hn_density, event_wind_lowlim;
		static const double smallest_viscosity, field_capacity_soil;
		static const bool jordy_new_snow, wind_pump, wind_pump_soil;

	private:
		static const mio::Config& cfg;
		static bool setStaticData(const std::string& variant);

		static double newSnowDensityPara(const std::string& i_hn_model,
		                                 double TA, double TSS, double RH, double VW, double HH);
		static double newSnowDensityEvent(const SnLaws::EventType& i_event_type, const CurrentMeteo& Mdata);
		static double newSnowDensityHendrikx(const double ta, const double tss, const double rh, const double vw);
		static double max_hn_density, event_wind_highlim;
		static EventType event_type;

		static double sn_dt; //Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
		static const bool __init;
		static std::string current_variant;
		static AlbedoModel currentAlbedoModel;
		static double albedoCage;
		static ViscosityCalVersion visc_cal;
		static TempDependence t_term;
		static double visc_time_fudge, visc_ice_fudge, visc_sp_fudge, visc_water_fudge;
		static bool setfix;
		static const bool band20;
		static unsigned int swa_nBands;
		static std::vector<double> swa_k, swa_pc, swa_fb;
		static const bool soil_evaporation;
		static const double rsoilmin, relsatmin, alpha_por_tor_soil, pore_length_soil;
		static const double montana_c_fudge, montana_vapor_fudge;
		static const double wind_ext_coef, displacement_coef, alpha_por_tor;
};

#endif
