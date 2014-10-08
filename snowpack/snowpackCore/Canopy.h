/*
 *  SNOWPACK stand-alone
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

#ifndef __CANOPY_H__
#define __CANOPY_H__

#include <snowpack/Constants.h>
#include <snowpack/DataClasses.h>
#include <snowpack/Hazard.h>
#include <snowpack/Utils.h>
#include <snowpack/Laws_sn.h>

class Canopy {

 	public:
		Canopy(const SnowpackConfig& i_cfg);

		static void cn_DumpCanopyData(FILE *OutFile, const CanopyData *Cdata, const SurfaceFluxes *Sdata, const double cos_sl);
		void runCanopyModel(CurrentMeteo &Mdata, SnowStation &Xdata, double roughness_length,
		                    double height_of_wind_val, const bool& alpine3d=false);
		static void cn_writeTimeSeriesAdd2LCanopy(FILE *OutFile, const CanopyData *Cdata);
		static const double can_alb_dry, can_alb_wet, can_alb_snow, krnt_lai; //public constants

 	private:
		double cn_f1(const double& ris);
		double cn_RootFraction(const double& zupper, const double& zlower);
		void cn_SoilWaterUptake(const size_t& SoilNode, const double& transpiration, ElementData* EMS);
		double cn_f4(const double& tempC);
		double cn_f2f4(const size_t& SoilNode, ElementData* EMS);
		double cn_f3(const double& vpd);
		double cn_IntCapacity(const double& tair, const double& density_of_new_snow, const double& lai);
		double cn_IntCapacitySnowMIP2(const double& tair, const double& density_of_mixed, const double& lai, double& hnws);
		double cn_IntUnload(const double& capacity, const double& storage);
		double cn_IntRate(const double& capacity, const double& storage, const double& prec,
		                  const double& direct);

		double cn_CanopyAlbedo(const double& tair, const double& wetfrac);
		double cn_TotalAlbedo(double CanAlb, double sigf, double SurfAlb, double DirectThroughfall,
		                      double CanopyClosureDirect, double RadFracDirect, double sigfdirect);

		double cn_CanopyShadeSoilCover(const double& HEIGHT, const double& COVER, const double& ELEV);
		double cn_CanopyWetFraction(const double& capacity, const double& storage);
		double cn_CanopyTransmissivity(const double& lai, const double& elev);

		void cn_LineariseNetRadiation(const CurrentMeteo& Mdata,const CanopyData& Cdata, const SnowStation& Xdata,
		                              double& iswrac, double& rsnet, double& ilwrac, double& r0,double& r1,
		                              const double& canopyalb, double& CanopyClosureDirect, double& RadFracDirect,
		                              const double& sigfdirect, double& r1p);
		void cn_LineariseNetRadiation2L(const CurrentMeteo& Mdata, const CanopyData& Cdata, const SnowStation& Xdata,
                                      double& iswrac, double& rsnet, double& ilwrac, double& r0,double& r1, double& r2,
                                      double& rt0, double& rt1, double& rt2, const double& canopyalb, double& CanopyClosureDirect, double& RadFracDirect,
                                      const double& sigfdirect, const double& sigftrunkdirect, double& r1p, double& r2p);
		void cn_TrunkEnergyBalance(double r2, double rt0, double rt1, double rt2, double ht0, double ht1, double let0, double let1,
                                          double HMt0, double HMt1, double &TT0, double &TT1, double TCANOPY, double &Ttrunk);
		void cn_LineariseSensibleHeatFlux(const double& ch_canopy, const double& tair, double& h0, double& h1, double scalingfactor);

		double cn_DSaturationPressureDT(const double& L, const double& T);
		void cn_LineariseLatentHeatFlux(const double& ce_canopy, const double& tc_old, const double& vpair,
		                                double& le0, double& le1, double scalingfactor);
		void cn_CalculateHeatMass(const double& height, const double& BasalArea, double& lai ,double& HMLeaves,  double& HMTrunks);

		void cn_LineariseConductiveHeatFlux(const double& tc_old, const double& HM, double& HM0, double& HM1,  double DT, double scalingfactor);

                void cn_CanopyEnergyBalance(const double& h0, const double& h1, const double& le0,
                                                         const double& le1, const double& HM0,  const double& HM1,
                                                         const double& ce_canopy,
                                                         const double& ce_condensation,
                                                         double& r0, double& r1, double& TCANOPY, double& RNCANOPY,
                                                         double& HCANOPY, double& LECANOPY);

		void cn_CanopyEnergyBalance2L(double& h0, double& h1, double& le0,
                                                         double& le1, double& HM0, double& HM1, double& TT0, double& TT1,
					                 const double& ce_canopy,
                                                         const double& ce_condensation,
                                                         double& r0, double& r1, double& r2, double& TCANOPY, double& Ttrunk, double& RNCANOPY,
                                                         double& HCANOPY, double& LECANOPY);

		void cn_CanopyEvaporationComponents(double& ce_canopy,
                                      double& ce_transpiration, double& LECANOPY,
                                      double& ta,double& I, double DT,
                                      double& CanopyEvaporation,
                                      double& INTEVAP, double& TRANSPIRATION,
                                      double& RNCANOPY, double& HCANOPY,double& TCANOPY,
                                      double& r0, double& r1, double& h0, double& h1,
                                      double& LECANOPYCORR,
                                      double& wetfraction, double& HM0, double& HM1);

		void cn_CanopyEvaporationComponents2L(double& ce_canopy,
							double& ce_transpiration, double& LECANOPY,
							double& ta, double& I, double DT,
							double& CanopyEvaporation,
							double& INTEVAP, double& TRANSPIRATION,
							double& RNCANOPY, double& HCANOPY,double& TCANOPY, double& Ttrunk,
							double& TT0, double& TT1,
							double& r0, double& r1, double& r2, double& h0, double& h1,
							double& LECANOPYCORR,
							double& wetfraction,
							double& HM0, double& HM1);
		double cn_psim(const double& xi);
		double cn_psih(const double& xi);
		double cn_RichardsonToAeta(double za, double TempAir, double DiffTemp, double Windspeed, double zom, double zoh, int maxitt);

		void cn_CanopyTurbulentExchange(const CurrentMeteo& Mdata, const double& refheight, const double& zomg,
								  const double& wetfraction, SnowStation& Xdata, double& ch_canopy,
								  double& ce_canopy, double& ce_transpiration,
								  double& ce_interception, double& ce_condensation);

		void cn_CanopyRadiationOutput(SnowStation& Xdata, CurrentMeteo& Mdata, double ac,
								double *iswrac, double *rswrac,
								double *iswrbc, double *rswrbc, double *ilwrac,
								double *rlwrac, double *ilwrbc, double *rlwrbc,
								double CanopyClosureDirect, double RadFracDirect, double sigfdirect, double sigftrunkdirect);

		static const double int_cap_snow, int_cap_rain, interception_timecoef;
		static const bool canopy_stabilitycorrection;
		static const double can_diameter, roughmom_to_canopyheight_ratio, displ_to_canopyheight_ratio, raincrease_snow;
		static const double canopytemp_maxchange_perhour, roughheat_to_roughmom_ratio, can_ch0, can_rs_mult, rsmin;
		static const double f3_gd, rootdepth, wp_fraction;

		std::string hn_density, hn_density_parameterization, variant;
		double hn_density_fixedValue, thresh_rain, calculation_step_length;
		bool useSoilLayers;
		// variables for canopy heat mass and 2-layer canopy
                bool CanopyHeatMass;
                bool Twolayercanopy;
                bool canopytransmission;
                bool forestfloor_alb;
                static const double biomass_heat_capacity, biomass_density, lai_frac_top_default, trunk_frac_height, trunkalb, et;
};

#endif //END of Canopy.h
