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
 * @file Canopy.h
 * @version 10.02
 */

#ifndef __CANOPY_H__
#define __CANOPY_H__

#include <snowpack/Constants.h>
#include <snowpack/DataClasses.h>
#include <snowpack/Hazard.h>
#include <snowpack/Utils.h>
#include <snowpack/Laws.h>
#include <snowpack/Laws_sn.h>

/*
 * CONSTANTS & ENUMERATIONS
 */
/**
 * @name CANOPY PARAMETERS
 * @brief History of changed values:
 * - 2007-12-20: update based on data from all SnowMIP2 sites, and calibration using Alptal data
 */
//@{
/// @brief INTERCEPTION
/**
 * @brief Specific interception capacity for snow (i_LAI) (mm/LAI) \n
 * Please note that this parameter is further multiplied with (0.27+46/new_snow_density[Ta]) following (Pomeroy et al, Hydr. Proc. 1998)
 * - 5.9 Spruce and 6.6 Pine (Schmidt&Glums,CanJForRes,1991)
 */
#define INT_CAP_SNOW 5.9
/// Specific interception capacity for rain (I_LAI) (mm/LAI)
#define INT_CAP_RAIN 0.3
/// Coef in interception function, see (Pomeroy et al,1998) where a value of 0.7 was found to be appropriate for hourly time-step, but smaller time steps require smaller values, 0.5 was found reasoanble by using the SnowMIP2 data (2007-12-09)
#define INTERCEPTION_TIMECOEF 0.5
/// @brief RADIATION BALANCE
#define CAN_ALB_DRY 0.09            // Albedo of dry canopy (calibr: 0.09, Alptal)
#define CAN_ALB_WET 0.09            // Albedo of wet canopy (calibr: 0.09, Alptal)
#define CAN_ALB_SNOW 0.35           // Albedo of snow covered albedo (calibr: 0.35, Alptal)
#define KRNT_LAI 0.7               // Radiation transmissivity parameter, in the range 0.4-0.8
                                   //     (calibr: 0.7, Alptal)
#define CANOPYTRANSMISSION 1       // (1=on, 0=off) optional radiation transfer model taking solar elevation into account
#define CAN_DIAMETER 1.0           // average canopy (tree) diameter [m], parameter in the new radiation transfer model

/// @brief TURBULENT HEAT EXCHANGE
/// @brief Stab. corr. aerodyn. resist. above and below canopy: 0=off and 1=on (Monin-Obukhov formulation)
#define CANOPY_STABILITYCORRECTION 1
/// @brief Ratio between canopy height and roughness length
#define ROUGHMOM_TO_CANOPYHEIGHT_RATIO 0.10
/// @brief As above for displacement height
#define DISPL_TO_CANOPYHEIGHT_RATIO 0.6667
/**
 * @brief Fractional increase of aerodynamic resistance for evaporation of intercepted snow.
 * - 10.0 from Koivusalo and Kokkonen (2002)
 * - 8.0 calibration with Alptal data
 */
#define RAINCREASE_SNOW 8.0
/// @brief Maximum allowed canopy temperature change (K hr-1)
#define CANOPYTEMP_MAXCHANGE_PERHOUR 10.0
/// @brief (~=1, but Not allowed to be exactly 1)
#define ROUGHHEAT_TO_ROUGHMOM_RATIO 0.9999
/// @brief minimum heat exchange (Wm-2K-1) at zero wind
#define CAN_CH0 3.0
/// @brief 1+CAN_RS_MULT = maximum factor to increase Cdata->rs below canopy
#define CAN_RS_MULT 3.0
/// @brief TRANSPIRATION
/// @brief Minimum canopy surface resistance, 500 (sm-1) is for needle leaf treas van den Hurk et al (2000) *75% Gustafsson et al (2003)
#define RSMIN 375.0
/**
 * @brief gd (Pa-1) parameter for canopy surface resistance response to vapour pressure:
 * - 0.0003 = trees (needle or broadleafs)
 * - 0=crops, grass, tundra etc
 */
#define F3_GD 0.0003
/// @brief Root depth, determining the soil layers influenced by root water uptake
#define ROOTDEPTH 1.0
/// @brief Wilting point, defined as a fraction of water content at field capacity (-)
#define WP_FRACTION 0.17
//@}

class Canopy {

 	public:
		Canopy(const mio::Config& i_cfg);

		static void cn_DumpCanopyData(FILE *OutFile, const SN_CANOPY_DATA *Cdata, const SN_SURFACE_DATA *Sdata, const double cos_sl);
		
		void runCanopyModel(SN_MET_DATA *Mdata, SN_STATION_DATA *Xdata, double roughness_length, 
						double height_of_wind_val);

 	private:

		double cn_f1(const double& ris);

		double cn_RootFraction(const double& zupper, const double& zlower);

		void cn_SoilWaterUptake(const int& SoilNode, const double& transpiration, SN_ELEM_DATA* EMS);

		double cn_f4(const double& tempC);

		double cn_f2f4(const int& SoilNode, SN_ELEM_DATA* EMS);

		double cn_f3(const double& vpd);

		double cn_IntCapacity(const double& tair, const double& density_of_new_snow, const double& lai);

		double cn_IntUnload(const double& capacity, const double& storage);

		double cn_IntRate(const double& capacity, const double& storage, const double& prec,
					   const double& direct, const double& tair);

		double cn_CanopyAlbedo(const double& tair, const double& wetfrac);

		double cn_TotalAlbedo(double CanAlb, double sigf, double SurfAlb, double SkyViewFraction, 
						  double CanopyClosureDirect, double RadFracDirect, double sigfdirect);

		double cn_CanopyShadeSoilCover(const double& HEIGHT, const double& COVER, const double& ELEV);

		double cn_CanopyWetFraction(const double& capacity, const double& storage);

		double cn_CanopyTransmissivity(const double& lai, const double& elev);

		void cn_LineariseNetRadiation(const SN_MET_DATA& Mdata,const SN_CANOPY_DATA& Cdata, const SN_STATION_DATA& Xdata,
								double& iswrac, double& rsnet, double& ilwrac, double& r0,double& r1, 
								const double& canopyalb, double& CanopyClosureDirect, double& RadFracDirect, 
								const double& sigfdirect, double& r1p);

		void cn_LineariseSensibleHeatFlux(const double& ch_canopy, const double& tair, double& h0, double& h1);

		double cn_DSaturationPressureDT(const double& L, const double& T);

		void cn_LineariseLatentHeatFlux(const double& ce_canopy, const double& tc_old, const double& vpair,
								  double& le0, double& le1);

		void cn_CanopyEnergyBalance(const double& h0, const double& h1, const double& le0, 
							   const double& le1, const double& vpair, const double& ce_canopy, 
							   const double& ce_condensation, const double& r1p, const double& CanopyClosure, 
							   double& r0, double& r1, double& TCANOPY, double& RNCANOPY, 
							   double& HCANOPY, double& LECANOPY);

		void cn_CanopyEvaporationComponents(double ce_canopy, //double ce_interception,
									 double ce_transpiration, double *LECANOPY,
									 double ta,double vpair,double I, double DT,
									 double *CanopyEvaporation,
									 double *INTEVAP, double *TRANSPIRATION,
									 double *RNCANOPY, double *HCANOPY,double *TCANOPY,
									 double *r0, double *r1, double h0, double h1, //double le0,double le1,
									 double *LECANOPYCORR, double r1p, double CanopyClosure, double wetfraction);

		double cn_psim(const double& xi);

		double cn_psih(const double& xi);

		double cn_RichardsonToAeta(double za, double TempAir, double DiffTemp, double Windspeed, double zom, double zoh, int maxitt);

		void cn_CanopyTurbulentExchange(const SN_MET_DATA& Mdata, const double& refheight, const double& zomg, 
								  const double& wetfraction, SN_STATION_DATA& Xdata, double& ch_canopy, 
								  double& ce_canopy, double& ce_transpiration, 
								  double& ce_interception, double& ce_condensation);
		
		void cn_CanopyRadiationOutput(SN_STATION_DATA *Xdata, SN_MET_DATA *Mdata, double ac, 
								double *iswrac, double *rswrac,
								double *iswrbc, double *rswrbc, double *ilwrac,
								double *rlwrac, double *ilwrbc, double *rlwrbc,
								double CanopyClosureDirect, double RadFracDirect, double sigfdirect);
		
		mio::Config cfg;
		int snp_soil;
		double thresh_rain, calculation_step_length;
};

#endif //END of Canopy.h
