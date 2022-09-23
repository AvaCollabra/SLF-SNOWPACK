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

#ifndef CANOPY_H
#define CANOPY_H

#include <snowpack/DataClasses.h>

#include <string>
#include <fstream>

/**
 * @brief Computes interception of precipitation and radiation, and reduction of windspeed
 * in a canopy layer above thesnow or soil surface.
 * This has been published in Gouttevin, i., M. Lehning, T. Jonas, D. Gustafsson, and Meelis Mölder,
 * <i>"A two-layer canopy model with thermal inertia for an improved snowpack energy balance below needleleaf forest
 * (model SNOWPACK, version 3.2. 1, revision 741)."</i>, Geoscientific Model Development <b>8.8</b>, pp 2379-2398, 2015.
 *
 * @section canopy_modeling Canopy modeling
 * -# 2layer canopy model:
 *	- key: TWO_LAYER_CANOPY = true [Snowpack]
 *	- logical in the code: Twolayercanopy
 *	- content: canopy is divided into a trunk layer (interception radiations with factor sigftrunk) and
 *		a leaf-layer (intercepting radiations with factor sigf). SW radiations reaching the ground are modified
 *		accordingly. An energy balance is computed for each layer, producing leaf-layer and trunk layer
 *		temperatures (TC and Ttrunk) that affect LW radiations to the ground.
 *		Optionally, trunks can get direct solar insolation (important for sparse canopies), look for
 *		CanClosDirTrunks in the code
 *	- further details: Gouttevin et al. (2014): A two-layer canopy with thermal inertia for an improved modelling
 *		of the sub-canopy snowpack energy-balance (in prep).
 *
 * -# canopy heat mass:
 *	- key: CANOPY_HEAT_MASS = true [Snowpack]
 *	- logical in the code: CanopyHeatMass
 *	- content: the canopy gets an heat mass (whole, or separate between trunks and leaves if Twolayercanopy) that
 *		adds a biomass heat flux term in the canopy energy balance.
 *	- an additionnal parameter is now required in the input/station.snoold file : CanopyBasalArea (m2/m2),
 *		to be placed after CanopyLeafAreaIndex. Example value for closed canopies like Alptal : 0.004.
 *	- further details: Gouttevin et al. (2014): A two-layer canopy with thermal inertia for an improved modelling
 *		of the sub-canopy snowpack energy-balance (in prep).
 *
 * -# forest-floor albedo:
 *	- key: FORESTFLOOR_ALB = true [Snowpack]
 *	- logical in the code: forestfloor_alb
 *	- content: Litter falling on the forest floor can reduce albedo. This effect is currently parameterized
 *		through an exponential decay of the albedo to a value of 0.3 with a time-constant of 7 days,
 *		based on parameterizations commonly used in Land-Surface models.
 *		There is room for improvement !
 *
 * @section canopy_comments Important comments:
 *	- Snowpack can take precipitation phase (relying on the psum_ph variable) for applications such as the
 *	  SnowMIP experiments (Rutter et al., 2009).
 *	- an additionnal parameter is now required in the input/station.snoold file : CanopyBasalArea (m2/m2),
 *	  to be placed after CanopyLeafAreaIndex.
 *	- Some cleaning was done to suppressed outputs that can be easily derived from other outputs.
 *	  There is now space for outputs specific to the 2layer model, which are written if variant = 2L_CANOPY in
 *	  [SnowpackAdvanced] (Canopy::writeTimeSeriesAdd2LCanopy).
 *
 */

class Canopy {

 	public:
		Canopy(const SnowpackConfig& i_cfg);

		static void DumpCanopyHeader(std::ofstream &fout);
		static void DumpCanopyUnits(std::ofstream &fout);
		static void DumpCanopyData(std::ofstream &fout, const CanopyData *Cdata,
                          const SurfaceFluxes *Sdata, const double cos_sl);
		bool runCanopyModel(CurrentMeteo &Mdata, SnowStation &Xdata,
                          const double roughness_length, const double height_of_wind_val,
                          const bool adjust_VW_height=true);
		static void writeTimeSeriesAdd2LCanopy(std::ofstream &fout, const CanopyData *Cdata);

 	private:

		void initStochasticUnload(const SnowpackConfig& cfg);

		double get_f1(const double ris) const;

		double RootFraction(const double zupper, const double zlower, const double rootdepth) const;

		void SoilWaterUptake(const size_t SoilNode, const double transpiration,
                         ElementData* ems, const double wp_fraction,
                         const double rootdepth, const double h_wilt);

		double get_f4(const double tempC) const;

		double get_f2f4(const size_t SoilNode, const ElementData* const ems,
                          const double wp_fraction,const double rootdepth) const;

		double get_f3(const double vpd, const double f3_gd) const;

		double IntCapacity(const CurrentMeteo& Mdata, const SnowStation& Xdata,
                        const bool force_rain=false) const;

		double IntUnload(const double capacity, const double storage) const;

		double StochasticUnload(const CurrentMeteo& Mdata, const double storage, const double solidfraction) const;

		void updateStorageAndUnloadElements(const double unload, ElementData& unloadedSnow, ElementData& snowStored);

		double IntRate(const double capacity, const double storage, const double prec,
		               const double direct, const double interception_timecoef) const;

		void updateInterceptionLayer(double interception, ElementData& snowStored, double density_new_snow,
		                             CurrentMeteo& Mdata);

		double CanopyAlbedo(const double tair, const double wetfrac, const SnowStation& Xdata) const;

		double TotalAlbedo(double CanAlb, double sigf, double SurfAlb, double DirectThroughfall,
		                   double CanopyClosureDirect, double RadFracDirect, double sigfdirect) const;

		void compactStoredSnow(ElementData& snowStored, double age, const CurrentMeteo& Mdata);

		double CanopyShadeSoilCover(const double height, const double cover, const double elev, const double can_diameter) const;

		double CanopyWetFraction(const double capacity, const double storage) const;

		double CanopyTransmissivity(const double lai, const double elev, const double krnt_lai) const;

		void LineariseNetRadiation(const CurrentMeteo& Mdata,const CanopyData& Cdata, const SnowStation& Xdata,
		                              double& iswrac, double& rsnet, double& ilwrac, double& r0,double& r1,
		                              const double canopyalb, double& CanopyClosureDirect, double& RadFracDirect,
		                              const double sigfdirect, double& r1p);

		void LineariseNetRadiation2L(const CurrentMeteo& Mdata, const CanopyData& Cdata, const SnowStation& Xdata,
                                      double& iswrac, double& rsnet, double& ilwrac, double& r0,double& r1, double& r2,
                                      double& rt0, double& rt1, double& rt2, const double canopyalb, double& CanopyClosureDirect, double& RadFracDirect,
                                      const double sigfdirect, const double sigftrunkdirect, double& r1p, double& r2p);

		void LineariseSensibleHeatFlux(const double chCanopy, const double tair, double& h0, double& h1, double scalingfactor);

		double DSaturationPressureDT(const double lh, const double temperature);

		void LineariseLatentHeatFlux(const double ceCanopy, const double tc_old, const double vpair,
		                                double& le0, double& le1, double scalingfactor);

		void CalculateHeatMass(const double height, const double BasalArea, double& lai ,double& HMLeaves,
						double& HMTrunks, const double biomass_density, const double biomass_heat_capacity);

		void LineariseConductiveHeatFlux(const double tc_old, const double hm, double& hm0, double& hm1,  const double dt, const double scalingfactor);

		void CanopyEnergyBalance(const double h0, const double h1, const double le0,
                                                         const double le1, const double hm0,  const double hm1,
                                                         const double ceCanopy,
                                                         const double ceCondensation,
                                                         double& r0, double& r1, double& TCanopy, double& RnCanopy,
                                                         double& HCanopy, double& LeCanopy);

		void CanopyEnergyBalance2L(double& h0, double& h1, double& le0,
                                                         double& le1, double& hm0, double& hm1, double& tt0, double& tt1,
                                                         const double ceCanopy,
                                                         const double ceCondensation,
                                                         double& r0, double& r1, double& r2, double& TCanopy, double& Ttrunk, double& RnCanopy,
                                                         double& HCanopy, double& LeCanopy);

		void CanopyEvaporationComponents(const double ceCanopy,
                                      const double ceTranspiration, double& LeCanopy,
                                      const double ta, const double i, const double dt,
                                      double& CanopyEvaporation,
                                      double& intevap, double& transpiration,
                                      double& RnCanopy, double& HCanopy,double& TCanopy,
                                      const double r0, const double r1, const double h0, const double h1,
                                      double& LECanopyCorr,
                                      const double wetfraction, const double hm0, const double hm1);

		void CanopyEvaporationComponents2L(const double ceCanopy,
                                      const double ceTranspiration, double& LeCanopy,
                                      const double ta, const double i, const double dt,
                                      double& CanopyEvaporation,
                                      double& intevap, double& transpiration,
                                      double& RnCanopy, double& HCanopy,double& TCanopy, double& Ttrunk,
                                      const double tt0, const double tt1,
                                      const double r0, const double r1, const double r2, const double h0, const double h1,
                                      double& LECanopyCorr,
                                      const double wetfraction,
                                      const double hm0, const double hm1);

		double get_psim(const double xi) const;

		double get_psih(const double xi) const;

		double RichardsonToAeta(double za, double TempAir, double DiffTemp, double Windspeed, double zom, double zoh, int maxitt) const;

		void CanopyTurbulentExchange(const CurrentMeteo& Mdata, const double refheight, const double zomg,
								  const double wetfraction, SnowStation& Xdata, double& chCanopy,
								  double& ceCanopy, double& ceTranspiration,
								  double& ce_interception, double& ceCondensation);

		void CanopyRadiationOutput(SnowStation& Xdata, const CurrentMeteo& Mdata, double ac,
								double &iswrac, double &rswrac,
								double &iswrbc, double &rswrbc, double &ilwrac,
								double &rlwrac, double &ilwrbc, double &rlwrbc,
								double CanopyClosureDirect, double RadFracDirect, double sigfdirect, double sigftrunkdirect);

		double GetStochasticUnloadProb1(double vw, double ta) const;
		double GetStochasticUnloadProb2(double vw, double ta) const;
		double GetStochasticUnloadProb3(double vw, double ta) const;
		double GetStochasticUnloadProb4(double vw, double ta) const;

		std::string hn_density, hn_density_parameterization, variant, watertransportmodel_soil;
		double hn_density_fixedValue, calculation_step_length;
		bool useSoilLayers;
		// variables for canopy heat mass and 2-layer canopy
		bool CanopyHeatMass;
		bool Twolayercanopy, Twolayercanopy_user;
		bool canopytransmission;
		bool forestfloor_alb;
		bool useUnload;
		double min_unload;
		bool stochasticUnload;
		double StochasticUnloadFrac;
		std::vector<double> stochasticParams;
		size_t stochasticDegree;
		std::function<double(double, double)> GetStochasticUnloadProb;


};

#endif
