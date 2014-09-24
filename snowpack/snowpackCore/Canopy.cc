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
/**
 * @file Canopy.cc
 * @version 10.02
 * @author David Gustafsson (davidg@kth.se), Michael Lehning
 * @bug     -
 * @brief Computes interception of precipitation and radiation, and reduction of windspeed
 * in a canopy layer above thesnow or soil surface.

 * - 2013-10-23 bis (I. Gouttevin, M. Bavay): simplification of the canopy energy Balance (EB) calculation, based on
		(i) suppression of the limitation of TC change by temp_maxchange_per_hour;
		    if temp_maxchange_per_hour is exceeded, turbulent coefficients
		    are re-calculated
		(ii) increased EB iterations (from 3 to 7)
		(iii) no radiation and latent fluxes updates within each EB iteration (pb of EB closure with the output variables if convergence is not reached)
		(iv) printing out the real (instead of the potential) canopy evaporation.

 * - 2013-10-23 (I. Gouttevin, M. Lehning, M. Bavay, D. Gustafsson): bug correction in r0change, psi, RAV and RAG calculations; re-calculation of the turbulent exchange coefficient in the EB-loop; a scaling factor is suggested for turbulent fluxes within the canopy EB but kept to 1. for the moment.

 * - 2007-12-20 The changes 071122 and 071022 included in rev 260 and commited to SVN.

 * - 2007-11-22 Changes in the Canopy following on the movement of SOLdata
                   members diff and elev to Mdata. Final commit of Note 2007-10-22
                   to SVN.

 * - 2007-11-22 Reversed order of notes so the latest is in the top

 * - 2007-10-22 Updated SVN version of Canopy model with the updates at
                   2007-05-14 that were never committed to svn.
                   Also note that the updates to canopy turbulence earlier this
                   summer were refined when David visited Michi in July. Important
                   now is that the Canopy model returns Cdata.zdispl = -0.7 if
                   canopy height is smaller than snow depth,which means that no
                   canopy computations were made.

 * - 2007-06-28 Made some additional checks in the canopy turbulence scheme to avoid
                   negative ustar, which happens if the roughness of the canopy is
                   smaller than the roughness of the snow surface (may happen when canopy
                   with small height is about to be covered by snow)

 * - 2007-05-14 Developments of the Canopy submodel for the Snowmip2 project have now been
                   included in the official svn version (also in ALpine3D). However,
                   the hack to consider liquid and frozen precipitation is not included,
                   thus, it is mainly the new radiation transfer scheme that differ from earlier
                   SVN revisions. The new canopy radiation scheme includes seperate transmission
                   coefficients for diffuse and direct radiation, where diffuse shortwave and
                   longwave has the same value depending on LAI and extinction coefficent only,
                   and direct shortwave is also dependent on solar elevation. The new scheme can
                   be switched on/off with the parameter CANOPYTRANSMISSION in Constants.c

 * - 2006-11-20 Introducing a variable canopy transmissivity for shortwave radiation
                   as a function of solar elevation (Gryning et al, 1999), which becomes
                   important at higher latitudes with low solar elevation angle. \n
                   Also making a hack for Snowmip2, to consider liquid and frozen
                   precipitation, as well as throughfall, separately.

 * - 2006-10-11 Changed default value of INTERCEPTION_TIMECOEF from 1 to 0.7 to obtain full
                   similarity with (Pomeroy et al., 1998: Hydr Proc)
                   Changed default values of most parameters in accordance with latest
                   findings from Alptal data (Gustafsson et al, 2007: Eur J For Res)

 * - 2006-10-09 Implemented (again) simplified linear scaling of canopy height
                   as model for roughness length and displacement height (taking the old
                   already existing parameters in use again, and DELETING parameter CAN_Z0M_COEF).

 * - 2006-06-02 Additional parameter to increase of canopy roughness length (we should
                   replace the Shaw Perreria function by a simple fraction of canopy height again,
                   thus only one parameter instead of ....) \n
                   Further more, some modifications to allow Canopy model to run without soil

 * - 2006-03-01  Introduced new definition of reference height above canopy
                   for Alpine3D applications. The reference height is now equal
                   to the Canopy height + HEIGHT_OF_WIND_VALUE if the model is
                   runned in ALPINE3D mode, and only HEIGHT_OF_WIND_VALUE in
                   SNOWPACK stand-alone mode. In both cases, the reference height
                   is minimized with Canopy height + 2.0m. \n
                   Furthermore, the final computation of the friction velocity
                   below the canopy was previously computed using the reference
                   height determined for the canopy - which was correct whenever the
                   HEIGHT_OF_WIND_VALUE > Canopy height + 2.0. However, in other cases
                   and especially with the new definition for ALPINE3D-mode,
                   this computation has to be done with the reference height
                   that actually is used in Meteo.c. This is now corrected.

 * - 2006-02-14  Actually, the rs_mult is still needed, however, it has a
                    new meaning... \n
                    The definition is now Cdata->rs = Cdata->rs * (1 + f(LAI))
                    , where f(LAI) = rs_mult * (1-exp(-LAI)).
                    In other words, the aerodynamic resistance below the canopy is
                    increased due to the LAI of the canopy, with a maximum factor
                    of 1+rs_mult. This is much better than before, when the
                    rs_mult was a constant (and hidden in the source code).
                    Now however, the increase of the resistance is also scaled by
                    the input LAI, so hopefully, there is less need for the users
                    to chose their own rs_mult. Anyway, rs_mult is added to the
                    CANOPY block in Constants.c \n
                    Furthermore, all those little special parameters related to
                    turbulent exchange from and below the canopy have been
                    calibrated using the ALPTAL data (see Gustafsson et al....),
                    and the new values are added to Constants.c

 * - 2006-01-09  Introduced a stability correction below the canopy (after
                    having realized that the global switch NEUTRAL doesnt have
                    there is effect when the CANOPY is switched on).
                    Now the rs_mult can be removed again...

 * - 2005-09-29: Introduce a multiplicative factor to increase the aerodynamic
                    resistance from the snow surface to the canopy level.
                    This factor seems to be around 2 for the Alptal data.
                    It can be a sign for the need of a stabililty correction
                    below the canopy (cf Koivusalo or Storck).
                    However, new stabililty functions for the snow surface
                    (switch NEUTRAL), didnt have there is effect on the exaggerated
                    contribution of turbulanet heat fluxes to snowmelt below
                    the canopy. Exaggerated longwave radiation below the canopy
                    can also be an explanaition.

 * - 2005-05-25: Cleaning up the main module cn_Canopy to enable average
                    (cumulated) flux output

 * - 2005-03-31: New simplified canopy turbulence scheme following
                    dual-source model presented by (Blyth et al)

 * - 2004-08-20: Field capacity is defined as the plant available water
                    in Alpine3D, therefore the wilting point is set to 0,
                    in order to get a reasonable soil water control on
                    transpiration.

                    Root distribution changed so that a root depth can be
                    specified: see functions SoilWaterUptake

 */

#include <snowpack/snowpackCore/Canopy.h>
#include <snowpack/snowpackCore/Snowpack.h>
#include <meteoio/MeteoIO.h>

#include <assert.h>

using namespace std;
using namespace mio;

/************************************************************
 * static section                                           *
 ************************************************************/

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
const double Canopy::int_cap_snow = 5.9;
/// Specific interception capacity for rain (I_LAI) (mm/LAI)
const double Canopy::int_cap_rain = 0.3;
/** Coef in interception function, see (Pomeroy et al,1998) where a value of 0.7 was
 * found to be appropriate for hourly time-step, but smaller time steps require smaller
 * values, 0.5 was found reasoanble by using the SnowMIP2 data (2007-12-09)
 */
const double Canopy::interception_timecoef = 0.5;

/// @brief RADIATION BALANCE
const double Canopy::can_alb_dry = 0.09;  // Albedo of dry canopy (calibr: 0.09, Alptal)
const double Canopy::can_alb_wet = 0.09;  // Albedo of wet canopy (calibr: 0.09, Alptal)
const double Canopy::can_alb_snow = 0.35; // Albedo of snow covered albedo (calibr: 0.35, Alptal)
const double Canopy::krnt_lai = 0.7;      // Radiation transmissivity parameter, in the range 0.4-0.8
                                          // (calibr: 0.7, Alptal)

const bool Canopy::canopytransmission = true; //optional radiation transfer model taking solar elevation into account
const double Canopy::can_diameter = 1.0; // average canopy (tree) diameter [m], parameter in the new radiation transfer model

/// @brief TURBULENT HEAT EXCHANGE
/// @brief Stab. corr. aerodyn. resist. above and below canopy: 0=off and 1=on (Monin-Obukhov formulation)
const bool Canopy::canopy_stabilitycorrection = true;
/// @brief Ratio between canopy height and roughness length
const double Canopy::roughmom_to_canopyheight_ratio = 0.10;
/// @brief As above for displacement height
const double Canopy::displ_to_canopyheight_ratio = 0.6667;

/**
 * @brief Fractional increase of aerodynamic resistance for evaporation of intercepted snow.
 * - 10.0 from Koivusalo and Kokkonen (2002)
 * - 8.0 calibration with Alptal data
 */
const double Canopy::raincrease_snow = 8.0;
/// @brief Maximum allowed canopy temperature change (K hr-1)
const double Canopy::canopytemp_maxchange_perhour = 7.0;
/// @brief (~=1, but Not allowed to be exactly 1)
const double Canopy::roughheat_to_roughmom_ratio = 0.9999;
/// @brief minimum heat exchange (Wm-2K-1) at zero wind
const double Canopy::can_ch0 = 3.0;
/// @brief 1+CAN_RS_MULT = maximum factor to increase Cdata->rs below canopy
const double Canopy::can_rs_mult = 3.0;
/// @brief TRANSPIRATION
/// @brief Minimum canopy surface resistance, 500 (sm-1) is for needle leaf treas van den Hurk et al (2000) *75% Gustafsson et al (2003)
const double Canopy::rsmin = 375.0;

/**
 * @brief gd (Pa-1) parameter for canopy surface resistance response to vapour pressure:
 * - 0.0003 = trees (needle or broadleafs)
 * - 0=crops, grass, tundra etc
 */
const double Canopy::f3_gd = 0.0003;
/// @brief Root depth, determining the soil layers influenced by root water uptake
const double Canopy::rootdepth = 1.0;
/// @brief Wilting point, defined as a fraction of water content at field capacity (-)
const double Canopy::wp_fraction = 0.17;
//@}


/**
 * @brief Dump 28 canopy parameters to Outfile
 * @param *OutFile Dump file
 * @param *Cdata
 * @param *Sdata
 * @param cos_sl Cosine of slope angle
 */
void Canopy::cn_DumpCanopyData(FILE *OutFile, const CanopyData *Cdata, const SurfaceFluxes *Sdata, const double cos_sl)
{
	fprintf(OutFile, ",%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
	// PRIMARY "STATE" VARIABLES
	Cdata->storage/cos_sl,      // intercepted water (mm or kg m-2)
	K_TO_C(Cdata->temp),        // temperature (degC)

	// SECONDARY "STATE" VARIABLES
	Cdata->canopyalb,           // albedo (1)
	Cdata->wetfraction,         // wet fraction
	Cdata->intcapacity/cos_sl,  // interception capacity (kg m-2)

	// RADIATIVE FLUXES (W m-2)
	Cdata->rsnet,               // net shortwave radiation, canopy
	Cdata->rlnet,               // net longwave radiation, canopy
	Cdata->rsnet+Cdata->rlnet,  // net radiation, canopy

	// HEAT FLUXES CANOPY (W m-2)
	-Cdata->sensible,           // sensible heat flux, canopy
	-Cdata->latentcorr,             // latent heat flux, canopy

	// WATER FLUXES CANOPY (kg m-2)
	Cdata->transp/cos_sl,       // transpiration
	Cdata->intevap/cos_sl,      // interception evaporation
	Cdata->interception/cos_sl, // interception
	Cdata->throughfall/cos_sl,  // throughfall
	Cdata->snowunload/cos_sl,   // unload of snow

	// TOTAL SURFACE FLUXES,EVAPORATION, ETC
	-Cdata->sensible+Sdata->qs, // sensible heat exchange of total surface
	-Cdata->latent+Sdata->ql,   // latent heat exchange of total surface
	Cdata->rlwrac,              // upward longwave radiation ABOVE canopy
	Cdata->ilwrac,              // downward longwave radiation ABOVE canopy
	Cdata->rlnet+Sdata->lw_net, // net longwave radiation
	Cdata->rswrac,              // upward shortwave above canopy
	Cdata->iswrac,              // downward shortwave radiation above canopy
	Cdata->rsnet+Sdata->qw,     // net shortwave radiation
	Cdata->totalalb,            // total albedo [-]
	Cdata->rlnet+Sdata->lw_net+Cdata->rsnet+Sdata->qw, // net radiation of total surface
	K_TO_C(pow(Cdata->rlwrac/Constants::stefan_boltzmann, 0.25)), // surface (radiative) temperature
	// precipitation (mm h-1) (surface variable is actually throughfall in canopy version)
	(Cdata->interception+Cdata->throughfall-Cdata->snowunload)/cos_sl,
	// evapotranspiration of total surface (mm h-1)
	(Cdata->transp+Cdata->intevap-(Sdata->mass[SurfaceFluxes::MS_SUBLIMATION]+Sdata->mass[SurfaceFluxes::MS_EVAPORATION]))/cos_sl);
}

/************************************************************
 * non-static section                                       *
 ************************************************************/
Canopy::Canopy(const SnowpackConfig& cfg)
        : hn_density(), hn_density_parameterization(), variant(),
          hn_density_fixedValue(Constants::undefined), thresh_rain(0.), calculation_step_length(0.), useSoilLayers(false)
{
	cfg.getValue("VARIANT", "SnowpackAdvanced", variant);

	cfg.getValue("SNP_SOIL", "Snowpack", useSoilLayers);

	cfg.getValue("THRESH_RAIN", "SnowpackAdvanced", thresh_rain);

	cfg.getValue("CALCULATION_STEP_LENGTH", "Snowpack", calculation_step_length);

	cfg.getValue("HN_DENSITY", "SnowpackAdvanced", hn_density);
	cfg.getValue("HN_DENSITY_PARAMETERIZATION", "SnowpackAdvanced", hn_density_parameterization);
	cfg.getValue("HN_DENSITY_FIXEDVALUE", "SnowpackAdvanced", hn_density_fixedValue);
}

/**
 * @brief multiplicative increase of canopy surface resistance as
 * a function of downward solar radiation (van den Burk et al (2000):
 * Offline validation of the ERA40 surface scheme, ECMWF Tech.Mem.295)
 * @param ris
 * @return double
 */
double Canopy::cn_f1(const double& ris)
{
	const double a = 0.81;
	const double b = 0.004;
	const double c = 0.05;
	const double f1 = ( a * ( 1.0 + b * ris ) ) / ( b * ris + c );
	if (f1 < 1.0) {
		return 1.0;
	}
	return (f1);

} // end cn_f1

/**
 * @brief Computes the FRACTION OF ROOTS in a soil layer between zupper and zlower meters
 * below the soil surface (van den Burk et al (2000): Offline validation of the ERA40
 * surface scheme, ECMWF Tech.Mem.295)
 * @param zupper
 * @param zlower
 * @return double
 */
double Canopy::cn_RootFraction(const double& zupper, const double& zlower)
{
	double rf = 0.0;

	// Constants.h: Canopy::rootdepth, default 0.5
	if ( zupper < Canopy::rootdepth ) {
		const double ar = 6.706; // evergreen needleleaf trees
		const double br = 2.175; // evergreen needleleaf trees
		// fraction of roots below root depth (according to exponential distribution)
		const double tail = 0.5 * (exp(-ar * Canopy::rootdepth)+ exp(-br * Canopy::rootdepth));
		// multiplicative factor to distribute tail on layers above root depth
		rf = ( ( 1. + tail / ( 1. - tail ) ) * 0.5 *
			(exp(-ar * zupper) + exp(-br * zupper)
			-exp(-ar * MIN(zlower, Canopy::rootdepth))
			-exp(-br * MIN(zlower, Canopy::rootdepth))));
	}

	return (rf);
}


/**
 * @brief Computes ROOT WATER UPTAKE in soil layers above rootdepth.
 * Theory:
	- 1) Transpiration is partitioned on soil layers according to fraction of roots.
	- 2) Root water uptake is limited to the plant available water in each layer,
	- defined by the wilting point and the field capacity.
	- 3) Wilting point and field capacity is dependent on soil type:
	- Fieldcapacity = ElementData::soilFieldCapacity
	- Wilting point = WP_FRACTION * Fieldcapacity
 * Last update: David Gustafsson, 2005-03-16.
 * @param SoilNode
 * @param *EMS
 * @param transpiration
 */
void Canopy::cn_SoilWaterUptake(const size_t& SoilNode, const double& transpiration, ElementData* EMS)
{
	// transpiration [mm]
	if ( transpiration == 0. )
		return;

	// Mass of water [kg m-2] that is to be extracted from the soil
	double waterresidual = transpiration;
	double waterresidual_real = transpiration;

	// Loop over soil layers above rootdepth
	double zupper = 0.;
	size_t RootLayer = SoilNode;
	for( size_t e = SoilNode-1; e --> 0; ) {
		// fraction of roots in layer
		const double rootfr = cn_RootFraction(zupper, zupper + EMS[e].L);
		const double water = transpiration;
		if( rootfr > 0.0 ){
			// Index of last layer with roots
			RootLayer = e;

			// Change in volumetric water content in layer
			const double d_theta_l = MIN ( MAX (0., ( EMS[e].theta[WATER] -
			                         Canopy::wp_fraction * EMS[e].soilFieldCapacity() )),
			                         rootfr*water / ( Constants::density_water * EMS[e].L ) );

			// residual water to be extracted in layers below
			waterresidual -= rootfr * water;
			waterresidual_real -= d_theta_l * Constants::density_water * EMS[e].L;

			// Update volumetric water content in layer
			EMS[e].theta[WATER] -= d_theta_l;
			assert(EMS[e].theta[WATER] >= -Constants::eps);
			EMS[e].theta[AIR] += d_theta_l;
			assert(EMS[e].theta[AIR] >= -Constants::eps);
		}
		// Depth of the upper edge of layer below
		zupper += EMS[e].L;
	}// End of loop

	// Extract the residual water uptake from first layer below rootzone
	if ( RootLayer > 0 ) {
		// modify by Moustapha if there is a problem .
		RootLayer -= 1;
	}
	const double d_theta = MIN ( MAX (0., ( EMS[RootLayer].theta[WATER] -
	                       Canopy::wp_fraction * EMS[RootLayer].soilFieldCapacity() ) ),
	                       waterresidual / ( Constants::density_water * EMS[RootLayer].L ) );

	EMS[RootLayer].theta[WATER] -= d_theta;
	assert(EMS[RootLayer].theta[WATER] >= -Constants::eps);
	EMS[RootLayer].theta[AIR] += d_theta;
	assert(EMS[RootLayer].theta[AIR] >= -Constants::eps);
	waterresidual_real -= d_theta * Constants::density_water * EMS[RootLayer].L;

	// Check if water content is below wilting point in last layer
	if ( waterresidual_real > 0.5 ) {
		// modify by Moustapha if there is problem .
		prn_msg(__FILE__, __LINE__, "wrn", mio::Date(), "Transpiration Error [mm]: %lf", waterresidual_real);
	}
}


/**
 * @brief multiplicative increase of canopy surface resistance as
 * a function of soil temperature, based on �gren et al (1976)
 * (A=0.8 and B=0.8 implies 1/f4=1 at 10oC)
 * Last update: David Gustafsson, 2005-03-16
 * Last Update: David Gustafsson, 2006-11-22>> A=1.75,B=0.5, implies that
		- 1/f4 increases much faster towards 1. 1/f4 = 0.9 already at 2oC, but 1/f4=1 still
		- at about 10oC.
 * @param tempC
 * @return double
 */
double Canopy::cn_f4(const double& tempC)
{
	const double F4_A = 1.75;
	const double F4_B = 0.5;

	// OBS tempC in C
	const double f4 = 1.0 / ( 1.0 - exp( -F4_A * pow( MAX( 0.00001, tempC ), F4_B ) ) );
	return (f4);
}


/**
 * @brief multiplicative increase of canopy surface resistance as
 * a function of liquid water content[1] and soil temperature [2]
 * in the root zone: [1] van den Hurk et al (2000): Offline validation
 * of the ERA40 surface scheme, ECMWF Tech.Mem.295, [2] �gren (1976)/Mellander (2005)
 * @param SoilNode
 * @param *EMS
 * @return double
 */
double Canopy::cn_f2f4(const size_t& SoilNode, ElementData* EMS)
{
	double f2_wpwp; double f2_wcap;
	double thet_act;
	double rootresidual = 1.;
	double f2 = 0.0; double f4 = 0.0;
	size_t RootLayer = SoilNode;

	// loop over layers:
	double zupper = 0.;
	for( size_t e = SoilNode-1; e --> 0; ) {
		// 1) root fraction in layer
		const double rootfr = cn_RootFraction(zupper, zupper + EMS[e].L);
		if( rootfr > 0.0 ){
			RootLayer = e;
			// 2) Field Capacity in layer
			f2_wcap = EMS[e].soilFieldCapacity();
			// 3) Wilting point in layer
			f2_wpwp = f2_wcap * Canopy::wp_fraction;
			// 4) Soil water content in layer (from a root's point of view)
			thet_act = MAX(f2_wpwp, EMS[e].theta[WATER]);
			// 4) Inversed soilwater stress weighted by root fractin in layer
			f2 += rootfr * (thet_act-f2_wpwp) / (f2_wcap - f2_wpwp);
			// 5) Soil temperature stress weighted by root fraction in layer
			f4 += cn_f4(K_TO_C(EMS[e].Te)) * rootfr;
			// 6) Update rootresidual and depth of upper edge of next layer
			rootresidual -= rootfr;
		}
		zupper += EMS[e].L;
	}
	// End of loop and now do the bottom layer
	if ( RootLayer > 0 ){
		RootLayer -= 1;
	}
	f2_wcap = EMS[RootLayer].soilFieldCapacity();
	f2_wpwp = f2_wcap * Canopy::wp_fraction;
	thet_act = MAX(f2_wpwp, EMS[RootLayer].theta[WATER]);
	f2 += rootresidual * (thet_act - f2_wpwp) / (f2_wcap - f2_wpwp);
	f4 += cn_f4(K_TO_C(EMS[RootLayer].Te)) * rootresidual;

	// Limit inverse of the f2 function between 0 and 1
	f2 = MAX( 0.00001, MIN( 1., f2 ) );
	// invert f2
	f2 = 1.0 / f2;
	// return maximum of f2 and f4
	return ( MAX(f2, f4) );
}


/**
 * @brief multiplicative increase of canopy surface resistance as
 * a function of atmospheric vapor pressure deficit (van den Burk et al (2000):
 * Offline validation of the ERA40 surface scheme, ECMWF Tech.Mem.295)
 * @param vpd
 * @return double
 */
double Canopy::cn_f3(const double& vpd)
{
	/*
	 * double F3_GD=0.0003; => now defined in Constants.h
	 * gd [Pa-1] value is for trees (needle or bradleafs), other veg. as crops,
	 * grass tundra etc should have gd=0;
	 */
	const double f3 = 1.0 / exp( -Canopy::f3_gd * vpd );
	return (f3);
}


double Canopy::cn_IntCapacity(const double& tair, const double& density_of_new_snow, const double& lai)
{
	// in the future, temperature threshold might be abandoned
	if ( K_TO_C(tair) < thresh_rain && density_of_new_snow > 0 ){
		return (Canopy::int_cap_snow * lai * ( 0.27+46.0 / density_of_new_snow ));
	} else {
		return (Canopy::int_cap_rain * lai);
	}
}

/**
 * @brief Intercepted snow or rain above the capacity is unloaded imediately
 * @param capacity
 * @param storage
 * @param *unload
 */
double Canopy::cn_IntUnload(const double& capacity, const double& storage)
{
	if ( storage > capacity ){
		return (storage - capacity);
	} else{
		return 0.;
	}
}

/**
 * @brief interception rate according to exponential function from Ashton ()
 * as formulated by i.e. Pomeroy and Hedstrom (1998)
 * @param capacity
 * @param storage
 * @param prec
 * @param *interception
 * @param direct
 */
double Canopy::cn_IntRate(const double& capacity, const double& storage, const double& prec, const double& direct)
{
	const double interception = MIN ( ( 1.0 - direct ) * prec,
                                Canopy::interception_timecoef * ( capacity - storage)*
                                ( 1.0 - exp( -(1.0 - direct) * prec / capacity ) ) );
	if ( interception < 0.0)
		return 0.;

	return interception;
}


double Canopy::cn_CanopyAlbedo(const double& tair, const double& wetfrac)
{
	// Albedo of partly "wet" canopy = weighted average of dry and wet parts
	if( tair > Constants::melting_tk ) {
		return (wetfrac * Canopy::can_alb_wet + (1. - wetfrac) * Canopy::can_alb_dry);
	} else {
		return (wetfrac * Canopy::can_alb_snow + (1. - wetfrac) * Canopy::can_alb_dry);
	}
}

/**
 * @brief Function returning the total surface albedo of a canopy covered snow or soil surface
 * @param CanAlb
 * @param sigf
 * @param SurfAlb
 * @param SkyViewFraction
 * @param CanopyClosureDirect
 * @param RadFracDirect
 * @param sigfdirect
 * @return double
 */
double Canopy::cn_TotalAlbedo(double CanAlb, double sigf, double SurfAlb, double SkyViewFraction,
			  double CanopyClosureDirect, double RadFracDirect, double sigfdirect)
{
	// Total surface albedo (diffuse fraction)
	const double albedo_diff = ( 1.0 - RadFracDirect ) * ( (sigf * CanAlb + SurfAlb * Optim::pow2(1.0 - sigf) /
			(1.0 - sigf * CanAlb * SurfAlb) ) * (1. - SkyViewFraction) + SurfAlb * SkyViewFraction);
	// Total surface albedo (direct fraction)
	const double albedo_dir = RadFracDirect * ( (sigfdirect * CanAlb + SurfAlb * Optim::pow2(1.0 - sigfdirect) /
			(1.0 - sigfdirect * CanAlb * SurfAlb) ) * CanopyClosureDirect + SurfAlb *
			(1.0 - CanopyClosureDirect) );
	return albedo_diff+albedo_dir;
}

/**
 * @brief Function returning the soil cover fraction of the canopy shade as a
 * function of canopy height, canopy diameter, and solar elevation
 * Computes the canopy shade soil cover fraction, as function of canopy height, crown diameter,
 * vertical canopy soil cover, and solar elevation angle. Derivation can be found in Gryning et al. (2001).
 * Boundary-Layer Meteorology, 99 (3), 465-488.
 * @param height
 * @param cover
 * @param elev in radiants
 * @return double
 */
double Canopy::cn_CanopyShadeSoilCover(const double& height, const double& cover, const double& elev)
{
	if ( elev > 0.0 ) {
		return MIN(1.0, cover * (1.0 + 4.0 * height / (Constants::pi * Canopy::can_diameter * tan(elev))));
	} else {
		return 1.0;
	}
}

/**
 * @brief fraction of the canopy surface covered by intercepted rain or snow
 * @param capacity
 * @param storage
 */
double Canopy::cn_CanopyWetFraction(const double& capacity, const double& storage)
{
	if ( storage > 0. ) {
		// limit the wet fraction to minimum 0.01 otherwise it will never completely dry
		return MAX(0.01, MIN(1.0, pow(storage / capacity, 2./3.)) );
	} else {
		return 0.;
	}
}

/**
 * @brief Transmissivity is now (optionally) a function of solar elevation (Chen et al, 1997)
 * @param *sigf
 * @param lai
 * @param elev
 * @return sigf
 */
double Canopy::cn_CanopyTransmissivity(const double& lai, const double& elev)
{
	const double pai = 0.; // pai [additional plant area index] could be a function of interception storage
	return (1. - exp(-Canopy::krnt_lai * (lai + pai) / MAX(sin(elev), 0.0001))); // Beer-Lambert type of law
}

/**
 * @brief This routine estimates the radiation balance of a vegetation covered grid
 * inputs:
 * incoming shortwave (RG) and longwave (RA) radiation
 * vegetation temperature (TV) and surface temperature of the ground below (TG)
 * vegetation albedo (AV) and ground albedo (AG)
 * vegetation shielding coefficient (SIGF) [shortwave and longwave]
 * emissivity of vegetation (EV) and ground (EG)
 * outputs:
 * Net longwave and shortwave radiation of vegetation (RGV,RAV) and ground (RGV,RGG)
 * Total grid albedo (AGRID) and grid surface radiation temperature (TGRID)
 * @param Mdata
 * @param Cdata
 * @param Xdata
 * @param iswrac
 * @param rsnet
 * @param ilwrac
 * @param r0
 * @param r1
 * @param canopyalb
 * @param CanopyClosureDirect
 * @param RadFracDirect
 * @param sigfdirect
 * @param r1p
 */
void Canopy::cn_LineariseNetRadiation(const CurrentMeteo& Mdata, const CanopyData& Cdata, const SnowStation& Xdata,
                                      double& iswrac, double& rsnet, double& ilwrac, double& r0,double& r1,
                                      const double& canopyalb, double& CanopyClosureDirect, double& RadFracDirect,
                                      const double& sigfdirect, double& r1p)
{
	// Variables used a lot
	const bool snow = (Xdata.getNumberOfElements()>Xdata.SoilNode);
	const double Tsfc = (snow)? Xdata.Ndata[Xdata.getNumberOfElements()].T : Mdata.ta;
	const double ag = (snow)? Xdata.Albedo : Xdata.SoilAlb;

	/*
	 * Canopy Closure = Canopy Soil Cover Fraction, is made a function of solar elevation for direct shortwave
	 * First, check whether the solar elevation and splitted radiation data makes there is sense
	 */

	const double elev = Mdata.elev;
	const double diffuse = Mdata.diff;
	const double direct = Mdata.iswr - diffuse;
	double RadFracDiffuse;
	if ( direct > 0.0 ) {
		RadFracDirect = direct / (diffuse + direct);
		RadFracDiffuse = 1.0 - RadFracDirect;
	} else {
		RadFracDirect = 0.0;
		RadFracDiffuse = 1.0;
	}
	const double sigf = Cdata.sigf;
	// Canopy Closure for diffuse shortwave and longwave
	const double CanopyClosure = 1. - Xdata.Cdata.direct_throughfall; //HACK: we already pass Cdata

	// Canopy Closure for direct shortwave
	if ( Canopy::canopytransmission ){
		CanopyClosureDirect = cn_CanopyShadeSoilCover(Cdata.height, CanopyClosure, elev);
	} else{
		CanopyClosureDirect = CanopyClosure;
	}

	// Shortwave radiation fluxes above and absorbed by canopy above
	iswrac = Mdata.iswr;

	// net absorbed by canopy

	// first only diffuse fraction
	rsnet = RadFracDiffuse * iswrac * (1. - canopyalb) * sigf *
		(1. + ag * (1. - sigf) / (1. - sigf * ag * canopyalb));
	// Longwave radiation above canopy:
	ilwrac = Mdata.ea * Constants::stefan_boltzmann * (Mdata.ta * Mdata.ta * Mdata.ta * Mdata.ta);

	// Longwave absorbed by canopy: auxiliary variables
	const double eg = 1.0; // emissivity of ground assumed to be 1
	const double star = 1. - sigf * (1. - Cdata.ec) * (1. - eg);
	const double psi = (1. - sigf) * (1. - eg) * Cdata.ec;


	// RNC = RNSC + RNLC: r0p and r1p correpsonds to RNC(t) = r0p + r1p * TC(t)^4
	const double r0p = rsnet + sigf * ((Cdata.ec + psi / star) *
	ilwrac + Cdata.ec * eg * Constants::stefan_boltzmann * Optim::pow4(Tsfc) / star);
	r1p = -sigf * (Cdata.ec * Constants::stefan_boltzmann + Cdata.ec * eg * Constants::stefan_boltzmann /
		star + psi * Constants::stefan_boltzmann / star);

	/*
	 * Linearise RNC arond TC(t) by using TC(t)=TC(t-1)^4+4*TC(t-1)^3*(TC(t)-TC(t-1)),
	 * which gives us r0 and r1, correpsonding to RNC(t)=r0+r1*TC(t)
	 */
	const double TC_old = Cdata.temp;

	r0 = r0p - 3. * r1p * Optim::pow4(TC_old);
	r1 = 4.* r1p * Optim::pow3(TC_old);

	// Scaling by CanopyClosure (= 1-SkyViewFraction)
	rsnet *= CanopyClosure;
	r0 *= CanopyClosure;
	r1 *= CanopyClosure;

	// Now, add the direct component with different CanopyClosure
	const double rsnetdir = CanopyClosureDirect * RadFracDirect * iswrac *
		(1. - canopyalb) * sigfdirect * (1. + ag * (1. - sigfdirect) / (1. - sigfdirect * ag * canopyalb));

	rsnet += rsnetdir;
	r0 += rsnetdir;

}

/**
 * @brief the sensible heat flux is already a linear function of TC(t),
 * whould be different with stability corrections.
 * @param ch_canopy
 * @param tair
 * @param *h0
 * @param *h1
 */
void Canopy::cn_LineariseSensibleHeatFlux(const double& ch_canopy, const double& tair, double& h0, double& h1, double scalingfactor)
{
	h1 = scalingfactor*ch_canopy;
	h0 = -scalingfactor*ch_canopy * tair;
}

/**
 * @brief Temperature derivative of the saturation pressure function
 * @param L
 * @param T
 * @return double
 */
double Canopy::cn_DSaturationPressureDT(const double& L, const double& T)
{
	double c2, c3;

	if ( L != Constants::lh_sublimation ) {
		//c1 = 610.780;
		c2 = 17.08085;
		c3 = 234.175 ;
	} else {
		//c1 = 610.714;
		c2 = 22.44294;
		c3 = 272.440 ;
	}

	//const double dpdt = lw_SaturationPressure(T) * L / (Constants::gas_constant * T * T);
	const double dpdt =  Atmosphere::waterSaturationPressure(T) * c2 * c3 / ((c3 + K_TO_C(T)) * (c3 + K_TO_C(T)));

	return(dpdt);
}

/**
 * @brief the latent heat flux LE = ce_canopy * (esat(TC(t)) - eair) is linearised around TC(t) by:
 * applying esat(TC(t)) = esat(TC(t - 1)) + Desat(TC(t - 1)) / DT * (TC(t) - TC(t - 1))
 * @param ce_canopy
 * @param tc_old
 * @param vpair
 * @param le0
 * @param le1
 */
void Canopy::cn_LineariseLatentHeatFlux(const double& ce_canopy, const double& tc_old, const double& vpair,
                                        double& le0, double& le1, double scalingfactor)
{
	if(tc_old > 273.15) {
		le1 = scalingfactor*ce_canopy * cn_DSaturationPressureDT(Constants::lh_vaporization, tc_old);
		le0 = scalingfactor*ce_canopy * (Atmosphere::waterSaturationPressure(tc_old) - vpair) - (le1) * tc_old;
	} else {
		le1 = scalingfactor*ce_canopy * cn_DSaturationPressureDT(Constants::lh_sublimation, tc_old);
		le0 = scalingfactor*ce_canopy * (Atmosphere::waterSaturationPressure(tc_old) - vpair) - (le1) * tc_old;
	}
}

/**
 * @brief Last update: 2007-05-10, David Gustafsson
 * @param r0
 * @param r1
 * @param h0
 * @param h1
 * @param le0
 * @param le1
 * @param TCANOPY
 * @param RNCANOPY
 * @param HCANOPY
 * @param LECANOPY
 * @param ce_canopy
 * @param ce_condensation
 */
void Canopy::cn_CanopyEnergyBalance(const double& h0, const double& h1, const double& le0,
							 const double& le1, const double& ce_canopy,
							 const double& ce_condensation,
							 double& r0, double& r1, double& TCANOPY, double& RNCANOPY,
							 double& HCANOPY, double& LECANOPY)

{
	/*
 	 * Introduced filter to reduce sensitivity of canopy energy balance solution:
	 * Maximum allowed change of canopy temperature per hour
	 * 1. infer TCANOPY from (R0 + R1 * TCANOPY) = (H0 + H1 * TCANOPY) + (LE0 + LE1 * TCANOPY)
	 */
	double TC_CHANGE = (h0 + le0 - r0) / (r1 - h1 - le1) - TCANOPY;

	TCANOPY += TC_CHANGE;

	// 3. and re-compute Rn, H, and LE
//      Previously, r0 and r1 were updated after each TC change computed in the EB loop.
//      With only 3 iterations of the EB, this prevented the closure of the canopy EB when looking at the output variables,
//      because TC had not completely converged.
//	The suggestion is to increase the iterations of the EB to 7 (instead of 3) and get rid off these artefacts.
// 	Similarly, LECANOPY is put to its computed value.
	RNCANOPY = r0 + r1  * TCANOPY ;
	HCANOPY = h0 + h1 * TCANOPY;
	LECANOPY = le0 + le1 * TCANOPY;

	// 3b. re-compute in case of condensation/sublimation on canopy
	if( LECANOPY < 0.0 ) {
		TCANOPY -= TC_CHANGE;
		TC_CHANGE = (h0 + le0 * ce_condensation / ce_canopy - r0) /
				(r1 - h1 - le1 * ce_condensation / ce_canopy) - TCANOPY;
		TCANOPY += TC_CHANGE;
	RNCANOPY = r0 + r1  * TCANOPY ;
	HCANOPY = h0 + h1 * TCANOPY;
        LECANOPY = le0* ce_condensation / ce_canopy + le1* ce_condensation / ce_canopy * TCANOPY;
	}
}


/**
 * @brief Partition latent heat flux on interception and transpiration
 * @param ce_canopy
 * @param ce_transpiration
 * @param *LECANOPY
 * @param ta
 * @param I
 * @param DT
 * @param *CanopyEvaporation
 * @param *INTEVAP
 * @param *TRANSPIRATION
 * @param *RNCANOPY
 * @param *HCANOPY
 * @param *TCANOPY
 * @param *r0
 * @param *r1
 * @param h0
 * @param h1
 * @param *LECANOPYCORR
 * @param wetfraction
 */
void Canopy::cn_CanopyEvaporationComponents(double ce_canopy, //double ce_interception,
				      double ce_transpiration, double& LECANOPY,
				      double ta,double I, double DT,
				      double& CanopyEvaporation,
				      double& INTEVAP, double& TRANSPIRATION,
				      double& RNCANOPY, double& HCANOPY,double& TCANOPY,
				      double& r0, double& r1, double h0, double h1, //double le0,double le1,
				      double& LECANOPYCORR,
                                      double wetfraction)
{
	if ( ta > Constants::freezing_tk ) {
		CanopyEvaporation = DT * 3600.0 * LECANOPY / Constants::lh_vaporization; // [mm]
	} else {
		CanopyEvaporation = DT * 3600.0 * LECANOPY / Constants::lh_sublimation;  // [mm]
	}

	if ( CanopyEvaporation <= 0.0 ) {
		INTEVAP = CanopyEvaporation; // [mm]
		TRANSPIRATION = 0.0;            // [mm]
		LECANOPYCORR = LECANOPY;
	} else {
		TRANSPIRATION = CanopyEvaporation * ce_transpiration * (1.0 - wetfraction) / ce_canopy;
		INTEVAP = CanopyEvaporation - TRANSPIRATION;
                if ( INTEVAP > I ) {
			INTEVAP = I;
			CanopyEvaporation = INTEVAP + TRANSPIRATION;
			if ( ta > Constants::freezing_tk ) {
				LECANOPYCORR = CanopyEvaporation * Constants::lh_vaporization / (DT * 3600.0);
			} else {
				LECANOPYCORR = CanopyEvaporation * Constants::lh_sublimation / (DT * 3600.0);
			}
			// re-compute TCANOPY from (R0 + R1 * TCANOPY) = (H0 + H1 * TCANOPY) + LECANOPYCORR
			TCANOPY  = (LECANOPYCORR + h0 - r0) / (r1 - h1);
			// Re-compute RNCANOPY, HCANOPY, and LECANOPY with new temperature
			RNCANOPY = r0 + r1  * TCANOPY ;
			HCANOPY  = h0 + h1 * TCANOPY;
			LECANOPY = LECANOPYCORR;
		} else {
			LECANOPYCORR = LECANOPY;
		}
	}
}

/**
 * @brief STABILITY FUNCTIONS FOR MOMENTUM
 * @param xi
 * @return double
 */
double Canopy::cn_psim(const double& xi)
{
	if ( xi <= 0.0 ) {
		// unstable case from Paulsen et al 1970
		const double x = pow((1. - 19. * xi), 0.25); // 19 from H�gstr�m, 1996
		return log((1. + x) * (1. + x) * (1. + x * x) / 8.) - 2 * atan(x) + mio::Cst::PI / 2.;
	} else {
		// stable case from Holstlag and Bruin, following Beljaars & Holtslag 1991
		const double a = 1.;
		const double b = 2./3.;
		const double c = 5.;
		const double d = 0.35;
		return -(a * xi + b * (xi - c/d) * exp(-d * xi) + b * c/d);
	}
}

/**
 * @brief STABILITY FUNCTIONS FOR HEAT
 * @param xi
 * @return double
 */
double Canopy::cn_psih(const double& xi)
{
	if ( xi <= 0) {
		// Unstable case. Integrated by Paulsen, 1970 from phi-functions found by Businger et al, 1971
		const double x = pow((1. - 11.6 * xi), 0.25);   // 11.6 from H�gstr�m, 1996
		return (2. * log((1. + x*x) / 2.) );
	} else {
		// Stable case, func=1, equation from Holtslag and De Bruin following Beljaars & Holstlag, 1991
		const double a = 1.;
		const double b = 2. / 3.;
		const double c = 5.;
		const double d = 0.35;
		return -(pow((1. + 2. / 3. * a * xi), 3. / 2.) + b * (xi - c/d) * exp(-d * xi) + b * c/d - 1.);
	}
}

/**
 * @brief SOLVES EQUATION eta=Ri*Ri2eta(eta) ITERATIVELY
 * @param za
 * @param TempAir
 * @param DiffTemp
 * @param Windspeed
 * @param zom
 * @param zoh
 * @param maxitt
 * @return double
 */
double Canopy::cn_RichardsonToAeta(double za, double TempAir, double DiffTemp,
			      double Windspeed, double zom, double zoh, int maxitt)
{
	// CALCULATE RICHARDSON NUMBER
	const double Ri = 9.81 * DiffTemp * za / (TempAir * Windspeed * Windspeed);
	// STEP 1: Compute function Ri2Eta(Eta)
	double Eta = 0.0;
	double Ri2eta = (log(za / zom) - cn_psim(Eta) + cn_psim(Eta * zom / za)) *
		(log(za / zom) - cn_psim(Eta) + cn_psim(Eta * zom / za)) /
		(log(za / zoh) - cn_psih(Eta) + cn_psih(Eta*zoh/za));
	// STEP 2: Compute error in terms of Ri using etaOld and Ri2eta(etaOld)
	double Error = Eta / Ri2eta - Ri;
	// STEP 3: solve iteratively
	const double acc = 0.0001;
	int itt=1;
	while ( fabs(Error) > acc && itt <= maxitt ) {
		// 3.1 new Eta
		Eta = Ri2eta * Ri;
		const double divider = (log(za / zoh) - cn_psih(Eta) + cn_psih(Eta * zoh / za)); //HACK: check with Davide
		if(divider!=0.) {
			Ri2eta = (log(za / zom) - cn_psim(Eta) + cn_psim(Eta * zom / za)) *
			(log(za / zom) - cn_psim(Eta) + cn_psim(Eta * zom / za)) /
			divider;
		} else {
			Ri2eta = 1.e12;
		}
		// 3.2 error in terms of Richardson number
		Error = Eta / Ri2eta - Ri;
		// 3.3 update number of iterations
		itt = itt + 1;
	}
	// STEP 4: Return new eta when convergance criteria is fullfilled
	Eta = Ri2eta * Ri;
	// check minimum Monin-Obukhov length (there isthing between 2 to some 100 meters in literature)
	if ( Eta > 0.0 ) {
		const double L = MAX (za, za / Eta);
		Eta = za / L;
	} else if ( Eta < 0.0 ) {
		const double L = MAX (-1000.0, za / Eta);
		Eta = za / L;
	}
	return (Eta);
}

/**
 * @brief Turbulent exchange coefficients canopy model
 * The Shuttelworth and Wallace (1985) approach contains rather complicated
 * formulas compared to the general low level of understanding of within canopy
 * turbulent exchange (which was rightly commented by the referees to our
 * paper Towards an improved...). Also, the old stability correction very
 * easily caused oscillations in the canopy energy balance.
 * To meet the referees and to get rid of the oscillations, David introduced
 * a more simple method, based only on roughness lengths (from Blyth, 1999)
 * Roughness lengths are estimated from canopy height and leaf area index.
 * A new stability correction method based on Monin-Obukhov theory was also
 * instead of the Richardsson number formulation
 * Last update: David, 2005-03-31
 * Stability correction now also below the canopy (using M-O) / DG, 2006-01-09
 * @param *Mdata
 * @param *Xdata
 * @param refheight
 * @param zomg
 * @param wetfraction
 * @param *ch_canopy
 * @param *ce_canopy
 * @param *ce_transpiration
 * @param *ce_interception
 * @param *ce_condensation
 */
void Canopy::cn_CanopyTurbulentExchange(const CurrentMeteo& Mdata, const double& refheight, const double& zomg,
                                        const double& wetfraction, SnowStation& Xdata, double& ch_canopy,
                                        double& ce_canopy, double& ce_transpiration,
                                        double& ce_interception, double& ce_condensation)
{
	const double karman = 0.4;
	CanopyData *Cdata = &Xdata.Cdata;
	const size_t nE = Xdata.getNumberOfElements();

	// check wind speed to be at least 0.1 m/s
	const double vw_local = (Mdata.vw>0.3)? Mdata.vw : 0.3;

	// canopy height above snow surface
	const double zcan = Cdata->height - (Xdata.cH - Xdata.Ground);

	/*
	 * 1. displacement and roughness (mom) according to Shaw and Perreira (1981)
	 * zdisplcan = 0.803 + 0.108 * CanDensMax - (0.462 - 0.086 * CanDensMax) *->
	 *-> exp(-(0.163 + 0.283 * CanDensMax) * Cdata->lai);
	 * zdisplcan = MAX (0., MIN (refheight - 0.5, zdisplcan * zcan));

	 * 1.3 roughness length
	 * const double EQ1 = (0.175 - 0.098 * CanDensMax) + (-0.098 + 0.045 * CanDensMax) * log10(Cdata->lai);
	 * const double EQ2 = (0.150 - 0.025 * CanDensMax) + (0.122 - 0.0135 * CanDensMax) * log10(Cdata->lai);
	 * zomc = MIN(RoughLmax, MAX(zcan * MIN(EQ1, EQ2), RoughLmin)) * CAN_Z0M_COEF;

	 * 1. displacement and roughness as simple scaling of canopy height.
	 * Please note:
	 * 1) Canopy roughness is not allowed to be smaller than roughness of
	 * snow surface as defined by zomg, otherwise ustar may be negative!
	* This is now guaranteed by the computation of RoughLmin = MAX(0.01,zomg)
	 * 2) refheight is already given as height relative the snow surface
	 */

	// Shaw Perreira parameters
	// double CanDensMax = 0.7;
	const double RoughLmin = 0.01;
	const double RoughLmax = 100.;
	const double zdisplcan = MAX (0., MIN (refheight - 0.5, Canopy::displ_to_canopyheight_ratio * zcan));
	const double zomc = MAX (MAX (RoughLmin, zomg), MIN (RoughLmax, Canopy::roughmom_to_canopyheight_ratio * zcan));

	/*
	 * 2. aerodynamic resistances simple approach (Blyth, 1999)
	 * 2.1 roughness length for scalars (heat and vapour)
	 */
	const double zohc = Canopy::roughheat_to_roughmom_ratio * zomc;
	const double zohg = Canopy::roughheat_to_roughmom_ratio * zomg;

	// update Cdata variables
	Cdata->z0m = zomc;
	Cdata->z0h = zohc;
	Cdata->zdispl = zdisplcan;

	// 2.2 Stability correction (adopted from Beljaars and Holtslag, 1991)
	double psim = 0.0;
	double psih = 0.0;

	if ( Canopy::canopy_stabilitycorrection ) {
		/*
		 * 2.2.1 Get Aeta = Monin-Obukhov stabilityparameter from Richardson number
		 */
		const double aeta = cn_RichardsonToAeta(refheight - zdisplcan, Mdata.ta,
			Mdata.ta - Cdata->temp, vw_local, zomc, zohc, 5);
		psih = -cn_psih(aeta) + cn_psih(aeta * zohc / (refheight - zdisplcan));
		psim = -cn_psim(aeta) + cn_psim(aeta * zomc / (refheight - zdisplcan));
	}

	// 2.3 FRICTION VELOCITY ABOVE CANOPY
	const double ustar = vw_local * karman / (log((refheight - zdisplcan) / zomc) + psim);

	// 2.4 TRANSFER COEFFICIENT FOR SCALARS ABOVE CANOPY
	double ch = Canopy::can_ch0 / (Constants::density_air * Constants::specific_heat_air)
		+ ustar * karman / (log((refheight - zdisplcan) / zohc) + psih);
	const double ch_e = ustar * karman / (log((refheight - zdisplcan) / zohc) + psih);

	// 2.5 AERODYNAMIC RESISTANCE ABOVE CANOPY
	Cdata->ra = 1. / ch;
	const double ra_e = 1. / ch_e;

	// 2.6 CANOPY TO CANOPY LEVEL RESISTANCE
	if ( log(zomc / zohc) > 0.0 ) {
		Cdata->rc = (log(zomc/zohc))/(karman * ustar );
	} else {
		Cdata->rc = 0.0;
	}

	// 2.7 SURFACE TO CANOPY LEVEL RESISTANCE
	if (log(zomc / zohg) > 0.) {
		Cdata->rs = (log(zomc / zohg)) / (karman * ustar ) * (1. + Canopy::can_rs_mult * (1. - exp(-Cdata->lai)));
	} else {
		Cdata->rs = 0.;
	}

	// 2.8 a stability correction is needed for the surface to canopy level resistance
	if ( Canopy::canopy_stabilitycorrection && (Cdata->rs > 0.) ) {
		double aeta_g = 0.;
		int i = 0;
		double rs_change = 1.;
		while( (i < 100) && (fabs(rs_change) > 0.0001) ) {
			i++;
			// 1. estimate ustar and ua(zdisplcan) above surface from ras and zomg, zohg, and zref = zdisplcan
			const double ustar_below1 = (1. / Cdata->rs) / karman * (log(zdisplcan / zohg)
			               - cn_psih(aeta_g) + cn_psih(aeta_g * zohg / (zdisplcan)));
			const double vw_zdisplcan = ustar_below1 / karman * (log(zdisplcan / zomg) -
			               cn_psim(aeta_g) + cn_psim(aeta_g * zomg / (zdisplcan)));
			// 2. estimate aeta above surface
			const double Tsup = (nE>0)? Xdata.Ndata[nE].T : Mdata.ta;
			aeta_g = cn_RichardsonToAeta(zdisplcan, Cdata->temp, Cdata->temp -
			                             Tsup, vw_zdisplcan, zomg, zohg, 5);
			// 3. new guess of ustar based on uadisplcan and new aeta_g
			const double ustar_below2 = vw_zdisplcan * karman / (log((zdisplcan)/zomg) -
			              cn_psim(aeta_g) + cn_psim(aeta_g * zomg / (zdisplcan)));

			// 4. TRANSFER COEFFICIENT FOR SCALARS below CANOPY
			ch = ustar_below2 * karman / (log((zdisplcan) / zohg) -
			     cn_psih(aeta_g) + cn_psih(aeta_g * zohg / (zdisplcan)));

			// 5. new guess for AERODYNAMIC RESISTANCE below CANOPY
			rs_change = 1. / ch - Cdata->rs;
			Cdata->rs = 1. / ch;
		}
	}

	 /*
	  * Surface resistance for transpiration (van den Hurk et al, 2000)
	  * In case there is no soil data, use the air temperature as guess for the soil temperature,
	  * and skip soil moisture function
	  */
	if ( useSoilLayers ) {
		Cdata->rstransp = Canopy::rsmin * cn_f1(Cdata->iswrac)*cn_f2f4(Xdata.SoilNode, &Xdata.Edata[0]) *
		                  cn_f3((1. - Mdata.rh) * Atmosphere::waterSaturationPressure(Mdata.ta)) / Cdata->lai;
	} else {
		const double Temp = (nE>0)? 0. : K_TO_C(Mdata.ta);
		Cdata->rstransp = Canopy::rsmin * cn_f1(Cdata->iswrac) * cn_f4(Temp) * cn_f3((1. - Mdata.rh) *
		                  Atmosphere::waterSaturationPressure(Mdata.ta)) / Cdata->lai;
	}
	// Exchange coefficients sensible heat
	ch_canopy = Constants::density_air * Constants::specific_heat_air / (Cdata->ra + Cdata->rc);

	// latent heat interception
	if ( Mdata.ta < 273.15 ) {
		ce_condensation  = 0.622 * Constants::lh_sublimation / (Constants::gas_constant_air * Mdata.ta
		                   * Canopy::raincrease_snow * (ra_e + Cdata->rc));// * MAX(0.1,wetfraction);
		ce_interception  = 0.622 * Constants::lh_sublimation / (Constants::gas_constant_air * Mdata.ta
		                   * Canopy::raincrease_snow * (ra_e + Cdata->rc));// * wetfraction;
		ce_transpiration = 0.0;
	} else {
		ce_condensation  = 0.622 * Constants::lh_vaporization / (Constants::gas_constant_air * Mdata.ta
		                   * (ra_e + Cdata->rc));// * MAX(0.1,wetfraction);
		ce_interception  = 0.622 * Constants::lh_vaporization / (Constants::gas_constant_air * Mdata.ta
		                   * (ra_e + Cdata->rc));// * wetfraction;
		ce_transpiration = 0.622 * Constants::lh_vaporization / (Constants::gas_constant_air * Mdata.ta
		                   * (ra_e + Cdata->rstransp + Cdata->rc));// * (1.0-wetfraction);
	}

	ce_canopy = ce_interception * MAX(0.001, wetfraction) + ce_transpiration * (1.0 - wetfraction);
}

/**
 * @brief Computes upward and downward radiation below and above canopy
 * @param Xdata
 * @param Mdata
 * @param ac
 * @param *iswrac
 * @param *rswrac
 * @param *iswrbc
 * @param *rswrbc
 * @param *ilwrac
 * @param *rlwrac
 * @param *ilwrbc
 * @param *rlwrbc
 * @param CanopyClosureDirect
 * @param RadFracDirect
 * @param sigfdirect
 */
void Canopy::cn_CanopyRadiationOutput(SnowStation& Xdata, CurrentMeteo& Mdata, double ac, double *iswrac, double *rswrac, double *iswrbc, double *rswrbc, double *ilwrac, double *rlwrac, double *ilwrbc, double *rlwrbc, double CanopyClosureDirect, double RadFracDirect, double sigfdirect)
{
	const bool snow = (Xdata.getNumberOfElements() > Xdata.SoilNode);
	const double Tsfc4 = (snow)? Optim::pow4(Xdata.Ndata[Xdata.getNumberOfElements()].T) : Optim::pow4(Mdata.ta);
	const double ag = (snow)? Xdata.Albedo : Xdata.SoilAlb;

	const double  TC4 = Optim::pow4(Xdata.Cdata.temp);
	const double  sigf = Xdata.Cdata.sigf;
	const double  ec = Xdata.Cdata.ec;
	const double  eg = 1.0;

	// Diffuse Shortwave radiation fluxes above and below canopy
	const double  rswrac_loc = *iswrac * (sigf * ac + ag * (1.0 - sigf) * (1.0 - sigf) / (1.0 - sigf * ac * ag));
	const double  iswrbc_loc = *iswrac * (1. - sigf) / (1.0 - sigf * ac * ag);
	const double  rswrbc_loc = iswrbc_loc * ag;

	// Direct Shortwave radiation fluxes above and below canopy
	const double  rswrac_loc2 = *iswrac * (sigfdirect * ac + ag * (1.0 - sigfdirect) * (1.0 - sigfdirect) / (1.0 - sigfdirect * ac * ag));
	const double  iswrbc_loc2 = *iswrac * (1. - sigfdirect) / (1.0 - sigfdirect * ac * ag);
	const double  rswrbc_loc2 = iswrbc_loc2 * ag;

	// Longwave radiation fluxes above and below canopy:
	const double  RAG = eg *(-Constants::stefan_boltzmann*Tsfc4+((1.-sigf)*(*ilwrac)+ sigf*ec*Constants::stefan_boltzmann*TC4 + eg*sigf*(1.-ec)*Constants::stefan_boltzmann*Tsfc4)/(1. - sigf * (1. - ec) * (1. - eg)));
	const double  RAV = sigf * ec * ( (*ilwrac) - 2.0 * Constants::stefan_boltzmann * TC4 + (Constants::stefan_boltzmann * ( eg * Tsfc4 + ec * sigf * TC4 * (1.-eg) ) + (1.0 - sigf) * (1.0 - eg) * (*ilwrac)) / (1.0 - sigf * (1.0 - ec)* ( 1.0 - eg)));
	*ilwrbc = RAG / eg + Constants::stefan_boltzmann * Tsfc4;
	*rlwrbc = (1 - eg)* (*ilwrbc) + eg * Constants::stefan_boltzmann * Tsfc4;
	*rlwrac = *ilwrac - RAG - RAV;

	// Scaling of results with CanopyClosureDiffuse and CanopyClosureDirect
	const double  CanopyClosureDiffuse = 1. - Xdata.Cdata.direct_throughfall;

	// Shortwave fluxes (diffuse)
	*rswrac = (rswrac_loc * CanopyClosureDiffuse + (*iswrac) * ag * (1.0 - CanopyClosureDiffuse)) * (1.0 - RadFracDirect);
	*iswrbc = (iswrbc_loc * CanopyClosureDiffuse + (*iswrac) * (1.0 - CanopyClosureDiffuse)) * (1.0 - RadFracDirect);
	*rswrbc = (rswrbc_loc * CanopyClosureDiffuse + *iswrac * ag * (1.0 - CanopyClosureDiffuse)) * (1.0 - RadFracDirect);

	// Shortwave fluxes (direct)
	*rswrac += (rswrac_loc2 * CanopyClosureDirect + (*iswrac) * ag * (1.0 - CanopyClosureDirect)) * RadFracDirect;
	*iswrbc += (iswrbc_loc2 * CanopyClosureDirect + (*iswrac) * (1.0 - CanopyClosureDirect)) * RadFracDirect;
	*rswrbc += (rswrbc_loc2 * CanopyClosureDirect + (*iswrac) * ag * (1.0 - CanopyClosureDirect)) *RadFracDirect;

	// Longwave fluxes (treat as diffuse)
	*rlwrac = *rlwrac * CanopyClosureDiffuse + Constants::stefan_boltzmann * eg * Tsfc4 * (1.0-CanopyClosureDiffuse);
	*ilwrbc = *ilwrbc * CanopyClosureDiffuse + *ilwrac * (1.0 - CanopyClosureDiffuse);
	*rlwrbc = *rlwrbc * CanopyClosureDiffuse + Constants::stefan_boltzmann * eg * Tsfc4 * (1.0-CanopyClosureDiffuse);
}

/**
 * @brief MAIN CANOPY FUNCTION CALLED BY Meteo.c
 * This routine computes interception of precipitation and radiation,
 * and reduction of turbulent exchange in a canopy layer above the ground.
 * Computations are made in the following order:
 * 1. Preliminar mass balance (interception and throughfall)
 * 2. Canopy surface energy balance (net radiation, sensible and latent
 * heat fluxes)
 * 3. Final mass balance (evaporation of intercepted water, and
 * transpiration
 * @param Mdata CurrentMeteo
 * @param Xdata Profile
 * @param roughness_length
 * @param height_of_wind_val
 * @param alpine3d changes the computation of zref. No idea what is the reason //HACK
 */
 //HACK: put an explicit name instead of "alpine3d"
void Canopy::runCanopyModel(CurrentMeteo &Mdata, SnowStation &Xdata, double roughness_length, double height_of_wind_val, const bool& alpine3d)
{
	const double hs = Xdata.cH - Xdata.Ground;
	const size_t nE = Xdata.getNumberOfElements();
	// First check, whether there is Canopy above the snow, i.e. whether s.th. needs to be done here
	if ( (Xdata.Cdata.height - 0.01) < hs ) {
		Xdata.Cdata.zdispl = -0.7;
		return;
	}

	// Check that some important initial values are within reasonable bounds
	if ( Xdata.Cdata.temp < 203.15 ){
		Xdata.Cdata.temp = 273.15;
	}
	if ( Xdata.Cdata.storage < 0.0 ) {
		Xdata.Cdata.storage = 0.0;
	}
	if ( Xdata.Cdata.lai <= 0.0 ) {
		Xdata.Cdata.lai = 0.0;
		return; //abort function execution, there is no canopy at this point
	}
	if ( Xdata.Cdata.height <= 0.0 ) {
		Xdata.Cdata.height = 0.0;
	}

	/*
	 * 1.1 compute the interception capacity [mm m-2]
	 * 1.1a Always new snow density as estimate of density in intercepted storage
	 */
	const double density_of_new_snow = SnLaws::compNewSnowDensity(hn_density, hn_density_parameterization, hn_density_fixedValue,
	                                                              Mdata, Xdata, Xdata.Cdata.temp, variant);

	// 1.1b Determine interception capacity [mm] as function of density of intercepted snow/rain
	const double intcapacity = cn_IntCapacity(Mdata.ta, density_of_new_snow, Xdata.Cdata.lai);

	// 1.2 compute direct unload [mm timestep-1], update storage [mm]
	double unload = cn_IntUnload(intcapacity, Xdata.Cdata.storage);

	if ( unload < 0.0 ) {
		prn_msg(__FILE__, __LINE__, "wrn", Mdata.date, "Negative unloading!!!");
		unload = 0.0;
	}
	Xdata.Cdata.storage -= unload;

	// 1.3 compute the interception [mm timestep-1] and update storage [mm]
	const double precipitation = Mdata.hnw;
	const double interception = cn_IntRate(intcapacity, Xdata.Cdata.storage, precipitation, Xdata.Cdata.direct_throughfall);

	Xdata.Cdata.storage += interception;

	// 1.4 compute the throughfall [mm timestep-1]
	const double throughfall = precipitation - interception + unload;

	Mdata.hnw = throughfall; // Please give the total amount for the time step

	// 2.1 prepare for canopy energy balance

	// Wetfraction update is moved to canopy energy balance loop  - use old value first
	double wetfrac = Xdata.Cdata.wetfraction;

	/*
	 * Radiation Transmissivity
	 * (could possibly be a function of interception - but is constant for the moment)
	 * Firstly, transmissivity of diffuse (and longwave) radiation
	 */
	Xdata.Cdata.sigf = cn_CanopyTransmissivity(Xdata.Cdata.lai, Constants::pi / 2.0);

	// Secondly, transmissivity of direct solar radiation
	const double sigfdirect = (Canopy::canopytransmission)? cn_CanopyTransmissivity(Xdata.Cdata.lai, Mdata.elev) : Xdata.Cdata.sigf;

	/*
	 * Reference Height [m above snow surface] for meteo input, at least 2 m above canopy height above snow surface
	 * 2006-03-01: introduced new definition of reference height above canopy for Alpine3D applications
	 * 2006-06-02: this should work also without soil data, since Xdata->Ground is initialized to 0.0
	*/
	const double ref_height = (alpine3d)? Xdata.Cdata.height+height_of_wind_val : height_of_wind_val;
	double zref = MAX( 2.0 + (Xdata.Cdata.height - hs), ref_height - hs );
	const double z0m_ground = (hs>0.03)? roughness_length : Xdata.BareSoil_z0;
	Mdata.z0 = (hs>0.03)? roughness_length : Xdata.BareSoil_z0;

	/*
	 * Turbulent Transport Coefficients:
	 * ch_canopy        = sensible heat
	 * ce_transpiration = latent heat from dry fraction, incl stomatal control.
	 * ce_interception  = latent heat from wet fraction
	 * ce_canopy        = latent heat total
	 * cm_canopy        = momentum (through the canopy, i.e to estimate wind speed below)
	*/
	double ch_canopy, ce_transpiration, ce_interception, ce_canopy, ce_condensation;
	cn_CanopyTurbulentExchange(Mdata, zref, z0m_ground, wetfrac, Xdata, ch_canopy, ce_canopy,
	                           ce_transpiration, ce_interception, ce_condensation);

	/*
	 * 2.2 Energy balance of the canopy
	 * The main purpose is to estimate the evapotranspiration loss,
	 * and the radiation balance of the canopy, which influences the
	 * snowpack below and the reflection/emittance to the atmosphere.
	 * Method:
	 * The energy balance of the canopy is assumed to be equal to,
	 * (1)   RnCanopy=HCanopy+LECanopy
	 * where RnCanopy is a function of a) incoming shortwave and longwave
	 * radiation, b) albedo, transmissivity, emissivity and temperature of
	 * the canopy, and c) albedo, emissivity and temperature of the ground
	 * below, taking multiple reflection into account (see documentation).
	 * Sensible and latent heat fluxes are computed using the standard
	 * bulk formulations, assuming an logarithmic wind profile above the
	 * canopy. The stomatal control of transpiration is represented by an
	 * additional surface resistance, estimated as a function of incoming
	 * solar radiation, atmospheric vapour pressure deficit and soil water
	 * content
	 * Numerical solution  following SiSPAT manual (Isabelle Braud, 2000):
	 * Equation (1) is linearised around the canopy temperature at time t
	 * using the temperature from the previous timestep t-1, so that:
	 * (2) RnCanopy  = r0  + r1  * TempCanopy(t)
	 * (3) HCanopy   = h0  + h1  * TempCanopy(t)
	 * (4) LECanopy  = le0 + le1 * TempCanopy(t)
	 * TempCanopy(t) is given by inserting eq (2)-(4) into (1).
	 * See the functions for r0,r1, h0,h1, le0, and le1 for details
	 * Alternative (to be implemented shortly):
	 * (1) RnCanopy = HCanopy + LECanopy + DQCanopy
	 * where DQCanopy = change of heat content of Canopy
	 * (5) DQCanopy = dq0+dq1 * TempCanopy(t)
	 * and dq/dt = HeatCapacity * (TempCanopy(t)-TempCanopy(t-1)
	 * HeatCapacity = HeatCapacityDry + IntStorage(t-1) * L
	 * INPUT: emissivity eground, surface temperature TGROUND, and albedo ground ALBGROUND
	 * Start of Energy Balance Loop, 3 iterations should be enough in most cases
	*/

	// local energy flux variables
	double RNCANOPY=0., HCANOPY=0., LECANOPY=0., LECANOPYCORR=0.;
	double iswrac, rswrac, iswrbc, rswrbc, ilwrac, rlwrac, ilwrbc, rlwrbc, rsnet=IOUtils::nodata;

	// local auxiliary variables
	double canopyalb=IOUtils::nodata;
	double r0, r1, h0, h1, le0, le1;
	double canopyclosuredirect=IOUtils::nodata, radfracdirect=IOUtils::nodata, r1p;

	double CanopyEvaporation=0., INTEVAP=0., TRANSPIRATION=0.;

	for ( int ebalitt = 0; ebalitt < 7; ebalitt++ ) {
		const double TC_OLD = Xdata.Cdata.temp; // Cdata.temp is updated in the iteration...

		// update ce_canopy as function of wetfraction
		ce_canopy = MAX (0.001 * ce_interception, ce_interception * wetfrac + ce_transpiration * (1.0 - wetfrac));
		ce_condensation = ce_interception * MAX (0.1, wetfrac);

		// canopy albedo
		canopyalb = cn_CanopyAlbedo(Mdata.ta, wetfrac);

		// compute properties r0 and r1 in eq (2) (and downward lw and sw for snowpack model)
		cn_LineariseNetRadiation(Mdata, Xdata.Cdata, Xdata, iswrac, rsnet, ilwrac, r0, r1,
		                         canopyalb, canopyclosuredirect, radfracdirect, sigfdirect, r1p);

		// compute properties h0 and h1 in eq (3)
		// NOTE: for sparse canopies turbulent fluxes should be scaled in the
		// canopy EB calculation; for the moment scalingfactor is 1
		cn_LineariseSensibleHeatFlux(ch_canopy, Mdata.ta, h0, h1, 1.);

		// compute properties le0 and le1 in eq (4)
		cn_LineariseLatentHeatFlux(ce_canopy, Xdata.Cdata.temp, Mdata.rh*Atmosphere::waterSaturationPressure(Mdata.ta), le0, le1, 1.);
		/* final canopy energy balance */
		cn_CanopyEnergyBalance(h0, h1, le0, le1, ce_canopy, ce_condensation,
		                       r0, r1, Xdata.Cdata.temp, RNCANOPY, HCANOPY, LECANOPY);

		/*
		 * Partition latent heat flux on interception and transpiration
		 * and correct energy balance for overestimated interception evaporation
		*/
		cn_CanopyEvaporationComponents(ce_canopy, ce_transpiration, LECANOPY, Mdata.ta, Xdata.Cdata.storage,
		                               M_TO_H(calculation_step_length), CanopyEvaporation, INTEVAP, TRANSPIRATION,
		                               RNCANOPY, HCANOPY, Xdata.Cdata.temp, r0, r1, h0, h1, LECANOPYCORR, wetfrac);
		const double newstorage = Xdata.Cdata.storage - INTEVAP;

		// wet surface fraction
		wetfrac = cn_CanopyWetFraction(intcapacity, newstorage);
		// Changes of temperature induce changes in stability correction.
		// re-computation of turbulent exchange coefficient is needed in case of big changes in TC.
		if (fabs(Xdata.Cdata.temp - TC_OLD) > Canopy::canopytemp_maxchange_perhour * M_TO_H(calculation_step_length)){
		cn_CanopyTurbulentExchange(Mdata, zref, z0m_ground, wetfrac, Xdata, ch_canopy, ce_canopy,
                    ce_transpiration, ce_interception, ce_condensation);
		}
		Xdata.Cdata.temp = (Xdata.Cdata.temp+TC_OLD)*0.5;
		wetfrac = (Xdata.Cdata.wetfraction+wetfrac)*0.5;
	} // End of Energy Balance Loop

	// Now REDUCE WaterContent in the Soil Elements --- Could also be part of WaterTransport.c
	if (useSoilLayers)
		cn_SoilWaterUptake(Xdata.SoilNode, TRANSPIRATION, &Xdata.Edata[0]);

	// final adjustment of interception storage due to evaporation
	Xdata.Cdata.storage = Xdata.Cdata.storage - INTEVAP;


	/*
	 * Preparation of output variables using += sign to allow for cumulated or averaged output
	 * (remember to reset these variables to 0 in Main.c before next integration step)
	*/
	// radiation above and below canopy
	cn_CanopyRadiationOutput(Xdata, Mdata, canopyalb, &iswrac, &rswrac, &iswrbc,&rswrbc,&ilwrac,
	                         &rlwrac,&ilwrbc,&rlwrbc,canopyclosuredirect,radfracdirect,sigfdirect);

	// longwave and shortwave radiation components
	Xdata.Cdata.iswrac += iswrac;
	Xdata.Cdata.rswrac += rswrac;
	Xdata.Cdata.iswrbc += iswrbc;
	Xdata.Cdata.rswrbc += rswrbc;
	Xdata.Cdata.ilwrac += ilwrac;
	Xdata.Cdata.rlwrac += rlwrac;
	Xdata.Cdata.ilwrbc += ilwrbc;
	Xdata.Cdata.rlwrbc += rlwrbc ;

	// Net longwave and shortwave radiation of canopy [W m-2]
	Xdata.Cdata.rlnet += RNCANOPY-rsnet;
	Xdata.Cdata.rsnet += rsnet;

	// atmospheric emissivity as seen by surface below canopy with regard to air temperature
	Mdata.ea = ilwrbc / (Constants::stefan_boltzmann * Optim::pow4(Mdata.ta));

	// downward and reflected shortwave below canopy
	Mdata.rswr = rswrbc;
	Mdata.iswr = iswrbc;

	// Adjust friction velocity below canopy using the same reference height as in Meteo.c
	zref = MAX(0.5, height_of_wind_val - hs);
	Mdata.ustar = 0.74 * log(zref / z0m_ground) / 0.4 / (Xdata.Cdata.ra + Xdata.Cdata.rs);

	// Canopy turbulent heat fluxes [W m-2]
	Xdata.Cdata.sensible += HCANOPY;
	Xdata.Cdata.latent += LECANOPY;
	Xdata.Cdata.latentcorr += LECANOPYCORR;

	// Canopy evaporative fluxes [mm]
	Xdata.Cdata.transp += TRANSPIRATION;
	Xdata.Cdata.intevap += INTEVAP;

	// Canopy mass fluxes [mm]

	Xdata.Cdata.interception += interception;
	Xdata.Cdata.throughfall += throughfall;
	Xdata.Cdata.snowunload += unload;

	// Canopy auxiliaries
	Xdata.Cdata.wetfraction = wetfrac;
	Xdata.Cdata.intcapacity += intcapacity;
	Xdata.Cdata.canopyalb += canopyalb;
	const double albedo = (nE>Xdata.SoilNode)? Xdata.Albedo : Xdata.SoilAlb;
	Xdata.Cdata.totalalb +=  cn_TotalAlbedo(canopyalb, Xdata.Cdata.sigf, albedo,
	                          Xdata.Cdata.direct_throughfall, canopyclosuredirect, radfracdirect, sigfdirect);
}
