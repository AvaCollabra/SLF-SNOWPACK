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
 * @file Canopy.cc
 * @version 10.02
 * @author David Gustafsson (davidg@kth.se) \n Michael Lehning
 * @bug     -
 * @brief Computes interception of precipitation and radiation, and reduction of windspeed
 * in a canopy layer above thesnow or soil surface.
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

#include <snowpack/Canopy.h>
#include <snowpack/Snowpack.h>

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
const double Canopy::canopytemp_maxchange_perhour = 10.0;
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
	-Cdata->latent,             // latent heat flux, canopy

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

Canopy::Canopy(const mio::Config& cfg)
{
	// Defines whether soil layers are used
	cfg.getValue("SNP_SOIL", "Snowpack", snp_soil);

	//Rain only for air temperatures warmer than threshold (degC)
	cfg.getValue("THRESH_RAIN", "SnowpackAdvanced", thresh_rain);

	cfg.getValue("CALCULATION_STEP_LENGTH", "Snowpack", calculation_step_length);

	cfg.getValue("HN_DENSITY", "SnowpackAdvanced", hn_density);
	cfg.getValue("HN_DENSITY_MODEL", "SnowpackAdvanced", hn_density_model);
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
	double f1;
	f1 = ( a * ( 1.0 + b * ris ) ) / ( b * ris + c );
	if (f1 < 1.0) {
		f1 = 1.0;
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
	const double ar = 6.706; // evergreen needleleaf trees
	const double br = 2.175; // evergreen needleleaf trees
	double tail;   // fraction of roots below rootdepth (according to exponential distribution)

	// Constants.h: Canopy::rootdepth, default 0.5
	if ( zupper < Canopy::rootdepth ) {
		// fraction of roots below root depth (according to exponential distribution)
		tail = 0.5 * (exp(-ar * Canopy::rootdepth)+ exp(-br * Canopy::rootdepth));
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
void Canopy::cn_SoilWaterUptake(const int& SoilNode, const double& transpiration, ElementData* EMS)
{

	// Constants.h: WP_FRACTION, default 0.01
	double rootfr, d_theta;
	double zupper = 0.;
	double water, waterresidual, waterresidual_real;
	int e, RootLayer;

	// transpiration [mm]
	if ( transpiration == 0. )
		return;

	// Mass of water [kg m-2] that is to be extracted from the soil
	water = waterresidual = waterresidual_real = transpiration;

	RootLayer = SoilNode;
	// Loop over soil layers above rootdepth
	for ( e = SoilNode-1; e > 0; e-- ) {
		// fraction of roots in layer
		rootfr = cn_RootFraction(zupper, zupper + EMS[e].L);
		if( rootfr > 0.0 ){
			// Index of last layer with roots
			RootLayer = e;

			// Change in volumetric water content in layer
			d_theta = MIN ( MAX (0., ( EMS[e].theta[WATER] -
					Canopy::wp_fraction * EMS[e].soilFieldCapacity() )),
						 rootfr*water / ( Constants::density_water * EMS[e].L ) );

			// residual water to be extracted in layers below
			waterresidual -= rootfr * water;
			waterresidual_real -= d_theta * Constants::density_water * EMS[e].L;

			// Update volumetric water content in layer
			EMS[e].theta[WATER] -= d_theta;
			EMS[e].theta[AIR] += d_theta;
		}
		// Depth of the upper edge of layer below
		zupper += EMS[e].L;
	}// End of loop

	// Extract the residual water uptake from first layer below rootzone
	if ( RootLayer > 0 ) {
		// modify by Moustapha if there is a problem .
		RootLayer -= 1;
	}
	d_theta = MIN ( MAX (0., ( EMS[RootLayer].theta[WATER] -
			Canopy::wp_fraction * EMS[RootLayer].soilFieldCapacity() ) ),
			waterresidual / ( Constants::density_water * EMS[RootLayer].L ) );

	EMS[RootLayer].theta[WATER] -= d_theta;
	EMS[RootLayer].theta[AIR] += d_theta;
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
	double f4;
	const double F4_A = 1.75;
	const double F4_B = 0.5;

	// OBS tempC in C
	f4 = 1.0 / ( 1.0 - exp( -F4_A * pow( MAX( 0.00001, tempC ), F4_B ) ) );
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
double Canopy::cn_f2f4(const int& SoilNode, ElementData* EMS)
{
	double f2_wpwp; double f2_wcap;
	double rootfr, thet_act;
	double zupper = 0.;
	double rootresidual = 1.;
	double f2 = 0.0; double f4 = 0.0;
	int e, RootLayer;

	RootLayer = SoilNode;
	// loop over layers:
	for ( e = SoilNode-1; e > 0; e-- ) {
		// 1) root fraction in layer
		rootfr = cn_RootFraction(zupper, zupper + EMS[e].L);
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
	double f3;
	f3 = 1.0 / exp( -Canopy::f3_gd * vpd );
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
 * @param tair
 */
double Canopy::cn_IntRate(const double& capacity, const double& storage, const double& prec,
                          const double& direct, const double& tair)
{
	double interception;

	if ( K_TO_C(tair) < thresh_rain ){
		interception = MIN ( ( 1.0 - direct ) * prec,
						 Canopy::interception_timecoef * ( capacity - storage )*
						 ( 1.0 - exp( -(1.0 - direct) * prec / capacity ) ) );
	} else{
		interception = MIN ( ( 1.0 - direct ) * prec,
						 Canopy::interception_timecoef * ( capacity - storage)*
						 ( 1.0 - exp( -(1.0 - direct) * prec / capacity ) ) );
	}

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
	double albedo_net;
	// Total surface albedo (diffuse fraction)
	albedo_net = ( 1.0 - RadFracDirect ) * ( (sigf * CanAlb + SurfAlb * (1.0 - sigf) * (1.0 - sigf) /
			(1.0 - sigf * CanAlb * SurfAlb) ) * (1. - SkyViewFraction) + SurfAlb * SkyViewFraction);
	// Total surface albedo (direct fraction)
	albedo_net += RadFracDirect * ( (sigfdirect * CanAlb + SurfAlb * (1.0 - sigfdirect) *
			(1.0 - sigfdirect) / (1.0 - sigfdirect * CanAlb * SurfAlb) ) * CanopyClosureDirect + SurfAlb *
			(1.0 - CanopyClosureDirect) );
	return albedo_net;
}

/**
 * @brief Function returning the soil cover fraction of the canopy shade as a
 * function of canopy height, canopy diameter, and solar elevation
 * Computes the canopy shade soil cover fraction, as function of canopy height, crown diameter,
 * vertical canopy soil cover, and solar elevation angle. Derivation can be found in Gryning et al. (2001).
 * Boundary-Layer Meteorology, 99 (3), 465-488.
 * @param HEIGHT
 * @param COVER
 * @param ELEV
 * @return double
 */
double Canopy::cn_CanopyShadeSoilCover(const double& height, const double& cover, const double& elev)
{
	double CanopyShadeSoilCover;

	if ( elev > 0.0 ){
		CanopyShadeSoilCover = MIN (1.0, cover * (1.0 + 4.0 * height / (Constants::pi * Canopy::can_diameter * tan(elev))));
	} else{
		CanopyShadeSoilCover = 1.0;
	}

	return CanopyShadeSoilCover;
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
		return MAX (0.01, MIN (1.0, pow(storage / capacity, 2. / 3.)));
	} else{
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
 * @param *Mdata
 * @param *Cdata
 * @param *Xdata
 * @param *iswrac
 * @param *rsnet
 * @param *ilwrac
 * @param *r0
 * @param *r1
 * @param canopyalb
 * @param *CanopyClosureDirect
 * @param *RadFracDirect
 * @param sigfdirect
 * @param *r1p
 */
void Canopy::cn_LineariseNetRadiation(const CurrentMeteo& Mdata, const CanopyData& Cdata, const SnowStation& Xdata,
							   double& iswrac, double& rsnet, double& ilwrac, double& r0,double& r1,
							   const double& canopyalb, double& CanopyClosureDirect, double& RadFracDirect,
							   const double& sigfdirect, double& r1p)
{
	double TC_old, TG, eg, ag;
	double star; double psi; double r0p;
	double CanopyClosure, elev, diffuse, direct, RadFracDiffuse, rsnetdir, sigf;

	// Variables used a lot
	if ( Xdata.getNumberOfElements() > 0 ) {
		TG = Xdata.Ndata[Xdata.getNumberOfElements()].T;	// ground surface temperature
		ag = Xdata.Albedo;
	} else {
		TG = Mdata.ta;
		ag = Xdata.SoilAlb;
	}
	/*
	 * Canopy Closure = Canopy Soil Cover Fraction, is made a function of solar elevation for direct shortwave
	 * First, check whether the solar elevation and splitted radiation data makes there is sense
	 */

	elev = Mdata.elev;
	diffuse = Mdata.diff;
	direct = Mdata.iswr - diffuse;
	if ( direct > 0.0 ) {
		RadFracDirect = direct / (diffuse + direct);
		RadFracDiffuse = 1.0 - RadFracDirect;
	} else {
		RadFracDirect = 0.0;
		RadFracDiffuse = 1.0;
	}
	sigf = Cdata.sigf;
	// Canopy Closure for diffuse shortwave and longwave
	CanopyClosure = 1. - Xdata.Cdata.direct_throughfall; //HACK: we already pass Cdata

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
	rsnet = RadFracDiffuse * iswrac * (1. - canopyalb) * Cdata.sigf *
		(1. + ag * (1. - Cdata.sigf) / (1. - Cdata.sigf * ag * canopyalb));
	// Longwave radiation above canopy:
	ilwrac = Mdata.ea * Constants::stefan_boltzmann * (Mdata.ta * Mdata.ta * Mdata.ta * Mdata.ta);

	// Longwave absorbed by canopy: auxiliary variables
	eg = 1.0; // emissivity of ground assumed to be 1
	star = 1. - Cdata.sigf * (1. - Cdata.ec) * (1. - eg);
	psi = (1. - Cdata.sigf) * (1. - Cdata.ec) * Cdata.ec;


	// RNC = RNSC + RNLC: r0p and r1p correpsonds to RNC(t) = r0p + r1p * TC(t)^4
	r0p = rsnet + Cdata.sigf * ((Cdata.ec + psi / star) *
		ilwrac + Cdata.ec * eg * Constants::stefan_boltzmann * (TG * TG * TG * TG) / star);
	r1p = -Cdata.sigf * (Cdata.ec * Constants::stefan_boltzmann + Cdata.ec * eg * Constants::stefan_boltzmann /
		star + psi * Constants::stefan_boltzmann / star);

	/*
	 * Linearise RNC arond TC(t) by using TC(t)=TC(t-1)^4+4*TC(t-1)^3*(TC(t)-TC(t-1)),
	 * which gives us r0 and r1, correpsonding to RNC(t)=r0+r1*TC(t)
	 */

	TC_old = Cdata.temp;

	r0 = r0p - 3. * r1p * (TC_old*TC_old*TC_old*TC_old);
	r1 = 4.* r1p * (TC_old * TC_old * TC_old);

	// Scaling by CanopyClosure (= 1-SkyViewFraction)
	rsnet *= CanopyClosure;
	r0 *= CanopyClosure;
	r1 *= CanopyClosure;

	// Now, add the direct component with different CanopyClosure
	rsnetdir = CanopyClosureDirect * RadFracDirect * iswrac *
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
void Canopy::cn_LineariseSensibleHeatFlux(const double& ch_canopy, const double& tair, double& h0, double& h1)
{
	h1 = ch_canopy;
	h0 = -ch_canopy * tair;
}

/**
 * @brief Temperature derivative of the saturation pressure function
 * @param L
 * @param T
 * @return double
 */
double Canopy::cn_DSaturationPressureDT(const double& L, const double& T)
{
	double dpdt;
	double c1, c2, c3;

	// dpdt = lw_SaturationPressure(T) * L / (Constants::gas_constant * T * T);
	if ( L != Constants::lh_sublimation ) {
		c1 = 610.780;
		c2 = 17.08085;
		c3 = 234.175 ;
	} else {
		c1 = 610.714;
		c2 = 22.44294;
		c3 = 272.440 ;
	}

	dpdt =  lw_SaturationPressure(T) * c2 * c3 / ((c3 + K_TO_C(T)) * (c3 + K_TO_C(T)));

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
								double& le0, double& le1)
{
	if(tc_old > 273.15) {
		le1 = ce_canopy * cn_DSaturationPressureDT(Constants::lh_vaporization, tc_old);
		le0 = ce_canopy * (lw_SaturationPressure(tc_old) - vpair) - (le1) * tc_old;
	} else {
		le1 = ce_canopy * cn_DSaturationPressureDT(Constants::lh_sublimation, tc_old);
		le0 = ce_canopy * (lw_SaturationPressure(tc_old) - vpair) - (le1) * tc_old;
	}
}

/**
 * @brief Last update: 2007-05-10, David Gustafsson
 * @param *r0
 * @param *r1
 * @param h0
 * @param h1
 * @param le0
 * @param le1
 * @param vpair
 * @param *TCANOPY
 * @param *RNCANOPY
 * @param *HCANOPY
 * @param *LECANOPY
 * @param ce_canopy
 * @param ce_condensation
 * @param r1p
 * @param CanopyClosure
 */
void Canopy::cn_CanopyEnergyBalance(const double& h0, const double& h1, const double& le0,
							 const double& le1, const double& vpair, const double& ce_canopy,
							 const double& ce_condensation, const double& r1p, const double& CanopyClosure,
							 double& r0, double& r1, double& TCANOPY, double& RNCANOPY,
							 double& HCANOPY, double& LECANOPY)

{
	// local variables
	double TC_CHANGE, TC_OLD, r0change, r1change;
	double MaxChange = Canopy::canopytemp_maxchange_perhour * M_TO_H(calculation_step_length);
	TC_OLD = TCANOPY;

	const double TC_OLD_POW3 = pow(TC_OLD, 3);
	const double TC_OLD_POW4 = TC_OLD_POW3 * TC_OLD;

	/*
 	 * Introduced filter to reduce sensitivity of canopy energy balance solution:
	 * Maximum allowed change of canopy temperature per hour
	 * 1. infer TCANOPY from (R0 + R1 * TCANOPY) = (H0 + H1 * TCANOPY) + (LE0 + LE1 * TCANOPY)
	 */
	TC_CHANGE = (h0 + le0 - r0) / (r1 - h1 - le1) - TCANOPY;

	// 2. minimize the rate of change to CANOPYTEMP_MAXCHANGE_PERHOUR [K hr - 1]
	TC_CHANGE = MAX (MIN (TC_CHANGE, MaxChange), -1.0 * MaxChange);
	TCANOPY += TC_CHANGE;

	// 3. and re-compute Rn, H, and LE
	r0change = 3.0 * r1p * CanopyClosure * (pow(TCANOPY, 4) - TC_OLD_POW4);
	r1change = pow(TCANOPY, 3) / TC_OLD_POW3;

	RNCANOPY = r0 + r0change + r1  * TCANOPY * r1change;
	HCANOPY = h0 + h1 * TCANOPY;
	LECANOPY = ce_canopy * (lw_SaturationPressure(TCANOPY) - vpair);

	// 3b. re-compute in case of condensation/sublimation on canopy
	if( LECANOPY < 0.0 ) {
		TCANOPY -= TC_CHANGE;
		TC_CHANGE = (h0 + le0 * ce_condensation / ce_canopy - r0) /
				(r1 - h1 - le1 * ce_condensation / ce_canopy) - TCANOPY;
		TC_CHANGE = MAX (MIN (TC_CHANGE, MaxChange), -1.0 * MaxChange);
		TCANOPY += TC_CHANGE;
		r0change = 3.0 * r1p * CanopyClosure * (pow(TCANOPY, 4) - TC_OLD_POW4);
		r1change = pow(TCANOPY, 3) / TC_OLD_POW3;
		RNCANOPY = r0 + r0change + r1  * TCANOPY * r1change;
		HCANOPY = h0 + h1 * TCANOPY;
		LECANOPY = ce_condensation * (lw_SaturationPressure(TCANOPY) - vpair);
	}
	r1 *= r1change;
	r0 += r0change;
}

/**
 * @brief Partition latent heat flux on interception and transpiration
 * @param ce_canopy
 * @param ce_transpiration
 * @param *LECANOPY
 * @param ta
 * @param vpair
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
 * @param r1p
 * @param CanopyClosure
 * @param wetfraction
 */
void Canopy::cn_CanopyEvaporationComponents(double ce_canopy, //double ce_interception,
				      double ce_transpiration, double *LECANOPY,
				      double ta,double vpair,double I, double DT,
				      double *CanopyEvaporation,
				      double *INTEVAP, double *TRANSPIRATION,
				      double *RNCANOPY, double *HCANOPY,double *TCANOPY,
				      double *r0, double *r1, double h0, double h1, //double le0,double le1,
				      double *LECANOPYCORR, double r1p, double CanopyClosure,
                                      double wetfraction)
{
	double TC_OLD, r0change, r1change;

	r1change = 1.0;
	r0change = 0.0;

	if ( ta > 273.15 ) {
		*CanopyEvaporation = DT * 3600.0 * (*LECANOPY) / Constants::lh_vaporization; // [mm]
	} else {
		*CanopyEvaporation = DT * 3600.0 * (*LECANOPY) / Constants::lh_sublimation;  // [mm]
	}

	if ( *CanopyEvaporation <= 0.0 ) {
		*INTEVAP = *CanopyEvaporation; // [mm]
		*TRANSPIRATION = 0.0;            // [mm]
		*LECANOPYCORR = *LECANOPY;
	} else {
		// *INTEVAP = *CanopyEvaporation * ce_interception * MAX(0.1,wetfraction) / ce_canopy;
		*TRANSPIRATION = *CanopyEvaporation * ce_transpiration * (1.0 - wetfraction) / ce_canopy;
		*INTEVAP = *CanopyEvaporation - (*TRANSPIRATION);
		if ( *INTEVAP > I ) {
			*INTEVAP = I;
			*CanopyEvaporation = *INTEVAP + (*TRANSPIRATION);
			if ( ta > 273.15 ) {
				*LECANOPYCORR = *CanopyEvaporation * Constants::lh_vaporization / (DT * 3600.0);
			} else {
				*LECANOPYCORR = *CanopyEvaporation * Constants::lh_sublimation / (DT * 3600.0);
			}
			TC_OLD = *TCANOPY;

			// re-compute TCANOPY from (R0 + R1 * TCANOPY) = (H0 + H1 * TCANOPY) + LECANOPYCORR
			*TCANOPY  = (*LECANOPYCORR + h0 - (*r0)) / (*r1 - h1);

			// Re-compute RNCANOPY, HCANOPY, and LECANOPY with new temperature
			r0change  = 3* r1p * CanopyClosure * (*TCANOPY * (*TCANOPY) *
					(*TCANOPY) * (*TCANOPY) - (TC_OLD * TC_OLD * TC_OLD * TC_OLD));
			r1change  = *TCANOPY * (*TCANOPY) * (*TCANOPY) / (TC_OLD * TC_OLD * TC_OLD);
			*RNCANOPY = *r0 + r0change + (*r1)  * (*TCANOPY) * r1change;
			*HCANOPY  = h0 + h1 * (*TCANOPY);
			*LECANOPY = ce_canopy * (lw_SaturationPressure(*TCANOPY) - vpair);
		} else {
			*LECANOPYCORR = *LECANOPY;
		}
	}
	*r1 *= r1change;
	*r0 += r0change;
}

/**
 * @brief STABILITY FUNCTIONS FOR MOMENTUM
 * @param xi
 * @return double
 */
double Canopy::cn_psim(const double& xi)
{
	double psim; double x;
	// UNSTABLE CASE FROM PAULSEN ET AL 1970
	if ( xi <= 0.0 ) {
		x = pow((1. - 19. * xi), 0.25); // 19 from H�gstr�m, 1996
		psim = log((1. + x) * (1. + x) * (1. + x * x) / 8.) - 2 * atan(x) + 3.141593 / 2.;
	} else {
		// stable case from Holstlag and Bruin, following Beljaars & Holtslag 1991
		double a = 1.;
		double b = 2. / 3.;
		double c = 5.;
		double d = 0.35;
		psim = -(a * xi + b * (xi - c/d) * exp(-d * xi) + b * c/d);
	}
	return(psim);
}

/**
 * @brief STABILITY FUNCTIONS FOR HEAT
 * @param xi
 * @return double
 */
double Canopy::cn_psih(const double& xi)
{
	double psih, x;
	// Unstable case. Integrated by Paulsen, 1970 from phi-functions found by Businger et al, 1971
	if ( xi <= 0) {
		x = pow((1. - 11.6 * xi), 0.25);   // 11.6 from H�gstr�m, 1996
		psih = 2. * log((1. + x * x) / 2.);
	} else {
		// Stable case, func=1, equation from Holtslag and De Bruin following Beljaars & Holstlag, 1991
		double a = 1.;
		double b = 2. / 3.;
		double c = 5.;
		double d = 0.35;
		psih = -(pow((1. + 2. / 3. * a * xi), 3. / 2.) + b * (xi - c/d) * exp(-d * xi) + b * c/d - 1.);
	}
	return(psih);
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
	double Ri2eta, Ri;
	double Eta = 0.0;
	double Error = 0.0;
	double acc = 0.0001;
	int itt = 1;
	double L;
	double divider;

	// CALCULATE RICHARDSON NUMBER
	Ri = 9.81 * DiffTemp * za / (TempAir * Windspeed * Windspeed);
	// STEP 1: Compute function Ri2Eta(Eta)
	Ri2eta = (log(za / zom) - cn_psim(Eta) + cn_psim(Eta * zom / za)) *
		(log(za / zom) - cn_psim(Eta) + cn_psim(Eta * zom / za)) /
		(log(za / zoh) - cn_psih(Eta) + cn_psih(Eta*zoh/za));
	// STEP 2: Compute error in terms of Ri using etaOld and Ri2eta(etaOld)
	Error = Eta / Ri2eta - Ri;
	// STEP 3: solve iteratively
	while ( fabs(Error) > acc && itt <= maxitt ) {
		// 3.1 new Eta
		Eta = Ri2eta * Ri;
		divider = (log(za / zoh) - cn_psih(Eta) + cn_psih(Eta * zoh / za)); //HACK: check with Davide
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
		L = MAX (za, za / Eta);
		Eta = za / L;
	} else if ( Eta < 0.0 ) {
		L = MAX (-1000.0, za / Eta);
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
	double vw_local;
	double karman = 0.4;
	double zcan;
	double zdisplcan;
	//  double EQ1;
	//  double EQ2;
	double zomc;
	double zohc;
	double zohg;
	double psim;
	double psih;
	double aeta;
	double ustar;
	double ch;
	double ch_e,ra_e;
	// Shaw Perreira parameters
	// double CanDensMax = 0.7;
	double RoughLmin = 0.01;
	double RoughLmax = 100.;

	// some new variables for the stability correction below canopy
	double ustar_below;
	double aeta_g;
	double rs_change;
	double vw_zdisplcan;
	int i;

	CanopyData *Cdata;

	// Dereference Cdata from Xdata
	Cdata = &Xdata.Cdata;

	// check wind speed to be at least 0.1 m/s
	vw_local = Mdata.vw;
	if ( vw_local < 0.3 ) {
		vw_local=0.3;
	}

	// canopy height above snow surface
	zcan = Cdata->height - (Xdata.cH - Xdata.Ground);

	/*
	 * 1. displacement and roughness (mom) according to Shaw and Perreira (1981)
	 * zdisplcan = 0.803 + 0.108 * CanDensMax - (0.462 - 0.086 * CanDensMax) *->
	 *-> exp(-(0.163 + 0.283 * CanDensMax) * Cdata->lai);
	 * zdisplcan = MAX (0., MIN (refheight - 0.5, zdisplcan * zcan));

	 * 1.3 roughness length
	 * EQ1 = (0.175 - 0.098 * CanDensMax) + (-0.098 + 0.045 * CanDensMax) * log10(Cdata->lai);
	 * EQ2 = (0.150 - 0.025 * CanDensMax) + (0.122 - 0.0135 * CanDensMax) * log10(Cdata->lai);
	 * zomc = MIN(RoughLmax, MAX(zcan * MIN(EQ1, EQ2), RoughLmin)) * CAN_Z0M_COEF;

	 * 1. displacement and roughness as simple scaling of canopy height.
	 * Please note:
	 * 1) Canopy roughness is not allowed to be smaller than roughness of
	 * snow surface as defined by zomg, otherwise ustar may be negative!
	* This is now guaranteed by the computation of RoughLmin = MAX(0.01,zomg)
	 * 2) refheight is already given as height relative the snow surface
	 */

	zdisplcan = MAX (0., MIN (refheight - 0.5, Canopy::displ_to_canopyheight_ratio * zcan));
	zomc = MAX (MAX (RoughLmin, zomg), MIN (RoughLmax, Canopy::roughmom_to_canopyheight_ratio * zcan));

	/*
	 * 2. aerodynamic resistances simple approach (Blyth, 1999)
	 * 2.1 roughness length for scalars (heat and vapour)
	 */
	zohc = Canopy::roughheat_to_roughmom_ratio * zomc;
	zohg = Canopy::roughheat_to_roughmom_ratio * zomg;

	// update Cdata variables
	Cdata->z0m = zomc;
	Cdata->z0h = zohc;
	Cdata->zdispl = zdisplcan;

	// 2.2 Stability correction (adopted from Beljaars and Holtslag, 1991)
	psim = 0.0;
	psih = 0.0;
	aeta = 0.0;

	if ( Canopy::canopy_stabilitycorrection ) {
		/*
		 * 2.2.1 Get Aeta = Monin-Obukhov stabilityparameter from Richardson number
		 */
		aeta = cn_RichardsonToAeta(refheight - zdisplcan, Mdata.ta,
			Mdata.ta - Cdata->temp, vw_local, zomc, zohc, 5);
		psih = -cn_psih(aeta) + cn_psih(aeta * zohc / (refheight - zdisplcan));
		psim = -cn_psim(aeta) + cn_psim(aeta * zomc / (refheight - zdisplcan));
	}

	// 2.3 FRICTION VELOCITY ABOVE CANOPY
	ustar = vw_local * karman / (log((refheight - zdisplcan) / zomc) + psim);

	// 2.4 TRANSFER COEFFICIENT FOR SCALARS ABOVE CANOPY
	ch = Canopy::can_ch0 / (Constants::density_air * Constants::specific_heat_air)
		+ ustar * karman / (log((refheight - zdisplcan) / zohc) + psih);
	ch_e = ustar * karman / (log((refheight - zdisplcan) / zohc) + psih);

	// 2.5 AERODYNAMIC RESISTANCE ABOVE CANOPY
	Cdata->ra = 1. / ch;
	ra_e = 1. / ch_e;

	// 2.6 CANOPY TO CANOPY LEVEL RESISTANCE
	if ( log(zomc / zohc) > 0.0 ) {
		Cdata->rc = (log(zomc/zohc))/(karman * ustar );
	} else {
		Cdata->rc = 0.0;
	}

	// 2.7 SURFACE TO CANOPY LEVEL RESISTANCE
	if (log(zomc / zohg) > 0.0) {
		Cdata->rs = (log(zomc / zohg)) / (karman * ustar ) * (1.0 + Canopy::can_rs_mult * (1 - exp(-Cdata->lai)));
	} else {
		Cdata->rs = 0.0;
	}

	// 2.8 a stability correction is needed for the surface to canopy level resistance
	if ( Canopy::canopy_stabilitycorrection && (Cdata->rs > 0.0) ) {
		aeta_g = 0;
		i = 0;
		rs_change = 1;
		while( (i < 100) && (fabs(rs_change) > 0.0001) )
		{
			i = i+1;
			// 1. estimate ustar and ua(zdisplcan) above surface from ras and zomg, zohg, and zref = zdisplcan
			ustar_below = (1. / Cdata->rs) / karman * (log(zdisplcan / zohg)
					- cn_psih(aeta_g) + cn_psih(aeta_g * zohg / (zdisplcan)));
			vw_zdisplcan = ustar_below / karman * (log(zdisplcan / zomg) -
					cn_psim(aeta_g) + cn_psim(aeta_g * zomg / (zdisplcan)));
			// 2. estimate aeta above surface
			if ( Xdata.getNumberOfElements() > 0){
				aeta_g = cn_RichardsonToAeta(zdisplcan, Cdata->temp, Cdata->temp -
					Xdata.Ndata[Xdata.getNumberOfElements()].T, vw_zdisplcan, zomg, zohg, 5);
			} else {
				aeta_g = cn_RichardsonToAeta(zdisplcan, Cdata->temp, Cdata->temp -
						Mdata.ta, vw_zdisplcan, zomg, zohg, 5);
			}
			// 3. new guess of ustar based on uadisplcan and new aeta_g
			ustar_below = vw_zdisplcan * karman / (log((zdisplcan)/zomg) -
					cn_psim(aeta_g) + cn_psim(aeta_g * zomg / (zdisplcan)));

			// 4. TRANSFER COEFFICIENT FOR SCALARS below CANOPY
			ch = ustar_below * karman / (log((zdisplcan) / zohg) -
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
	if ( snp_soil ) {
		Cdata->rstransp = Canopy::rsmin * cn_f1(Cdata->iswrac)*cn_f2f4(Xdata.SoilNode,&Xdata.Edata[0]) *
				cn_f3((1. - Mdata.rh) * lw_SaturationPressure(Mdata.ta)) / Cdata->lai;
	} else {
		if ( Xdata.getNumberOfElements() > 0 ) {
			Cdata->rstransp = Canopy::rsmin * cn_f1(Cdata->iswrac) * cn_f4(0.0) * cn_f3((1. - Mdata.rh) *
						lw_SaturationPressure(Mdata.ta)) / Cdata->lai;
		} else {
			Cdata->rstransp = Canopy::rsmin * cn_f1(Cdata->iswrac) * cn_f4(K_TO_C(Mdata.ta)) *
			cn_f3((1. - Mdata.rh) * lw_SaturationPressure(Mdata.ta)) / Cdata->lai;
		}
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
	double TC, TG, ag, sigf, ec, eg, RAG, RAV, CanopyClosureDiffuse, rswrac_loc;
	double rswrbc_loc, rswrac_loc2, iswrbc_loc2, rswrbc_loc2, iswrbc_loc;

	// Variables used a lot
	if ( Xdata.getNumberOfElements() > 0 ) {
		TG = Xdata.Ndata[Xdata.getNumberOfElements()].T;	// ground surface temperature
		ag=Xdata.Albedo;
	} else {
		TG = Mdata.ta;
		ag = Xdata.SoilAlb;
	}

	TC = Xdata.Cdata.temp;
	sigf = Xdata.Cdata.sigf;
	ec = Xdata.Cdata.ec;
	eg = 1.0;

	// Diffuse Shortwave radiation fluxes above and below canopy
	rswrac_loc = *iswrac * (sigf * ac + ag * (1.0 - sigf) * (1.0 - sigf) / (1.0 - sigf * ac * ag));
	iswrbc_loc = *iswrac * (1. - sigf) / (1.0 - sigf * ac * ag);
	rswrbc_loc = iswrbc_loc * ag;

	// Direct Shortwave radiation fluxes above and below canopy
	rswrac_loc2 = *iswrac * (sigfdirect * ac + ag * (1.0 - sigfdirect) * (1.0 - sigfdirect) / (1.0 - sigfdirect * ac * ag));
	iswrbc_loc2 = *iswrac * (1. - sigfdirect) / (1.0 - sigfdirect * ac * ag);
	rswrbc_loc2 = iswrbc_loc2 * ag;

	// Longwave radiation fluxes above and below canopy:
	RAG = (1. - sigf) * eg * ( *ilwrac - Constants::stefan_boltzmann * TG * TG * TG * TG) - eg * ec * sigf * Constants::stefan_boltzmann * (TG*TG*TG*TG - TC*TC*TC*TC) / (1. - sigf * (1. - ec) * (1. - eg));
	RAV = sigf * (ec * ( *ilwrac - Constants::stefan_boltzmann * TC*TC*TC*TC) + (Constants::stefan_boltzmann * ec * eg * (TG*TG*TG*TG - TC*TC*TC*TC) + (1.0 - sigf) * (1.0 - eg) * ec * ( *ilwrac - Constants::stefan_boltzmann * TC*TC*TC*TC)) / (1.0 - sigf * (1.0 - ec)* ( 1.0 - eg)));
	*ilwrbc = RAG / eg + Constants::stefan_boltzmann * TG*TG*TG*TG;
	*rlwrbc = - (1 - eg)* (*ilwrbc) + eg * Constants::stefan_boltzmann * TG*TG*TG*TG;
	*rlwrac = *ilwrac - RAG - RAV;

	// Scaling of results with CanopyClosureDiffuse and CanopyClosureDirect
	CanopyClosureDiffuse = 1. - Xdata.Cdata.direct_throughfall;

	// Shortwave fluxes (diffuse)
	*rswrac = (rswrac_loc * CanopyClosureDiffuse + (*iswrac) * ag * (1.0 - CanopyClosureDiffuse)) * (1.0 - RadFracDirect);
	*iswrbc = (iswrbc_loc * CanopyClosureDiffuse + (*iswrac) * (1.0 - CanopyClosureDiffuse)) * (1.0 - RadFracDirect);
	*rswrbc = (rswrbc_loc * CanopyClosureDiffuse + *iswrac * ag * (1.0 - CanopyClosureDiffuse)) * (1.0 - RadFracDirect);

	// Shortwave fluxes (direct)
	*rswrac += (rswrac_loc2 * CanopyClosureDirect + (*iswrac) * ag * (1.0 - CanopyClosureDirect)) * RadFracDirect;
	*iswrbc += (iswrbc_loc2 * CanopyClosureDirect + (*iswrac) * (1.0 - CanopyClosureDirect)) * RadFracDirect;
	*rswrbc += (rswrbc_loc2 * CanopyClosureDirect + (*iswrac) * ag * (1.0 - CanopyClosureDirect)) *RadFracDirect;

	// Longwave fluxes (treat as diffuse)
	*rlwrac = *rlwrac * CanopyClosureDiffuse + Constants::stefan_boltzmann * eg * TG * TG * TG * TG * (1.0-CanopyClosureDiffuse);
	*ilwrbc = *ilwrbc * CanopyClosureDiffuse + *ilwrac * (1.0 - CanopyClosureDiffuse);
	*rlwrbc = *rlwrbc * CanopyClosureDiffuse + Constants::stefan_boltzmann * eg * TG * TG * TG * TG * (1.0-CanopyClosureDiffuse);
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
 */

void Canopy::runCanopyModel(CurrentMeteo *Mdata, SnowStation *Xdata, double roughness_length, double height_of_wind_val)
{

	// local mass flux variables
	double precipitation, throughfall, unload, interception;
	double CanopyEvaporation=0., INTEVAP=0., TRANSPIRATION=0.;

	// local energy flux variables
	double RNCANOPY=0., HCANOPY=0., LECANOPY=0., LECANOPYCORR=0.;
	double iswrac, rswrac, iswrbc, rswrbc, ilwrac, rlwrac, ilwrbc, rlwrbc, rsnet;

	// local auxiliary variables
	double intcapacity, wetfrac, canopyalb;
	double ch_canopy, ce_transpiration, ce_interception, ce_canopy, ce_condensation;
	double r0, r1, h0, h1, le0, le1;
	double zref, z0m_ground;
	double density_of_new_snow, TC_OLD, newstorage;
	double canopyclosuredirect, radfracdirect, sigfdirect, r1p;
	int ebalitt;

	// First check, whether there is Canopy above the snow, i.e. whether s.th. needs to be done here

	if ( (Xdata->Cdata.height - 0.01) < (Xdata->cH - Xdata->Ground) ) {
		Xdata->Cdata.zdispl = -0.7;
		return;
	}

	// Check that some important initial values are within reasonable bounds
	if ( Xdata->Cdata.temp < 203.15 ){
		Xdata->Cdata.temp = 273.15;
	}
	if ( Xdata->Cdata.storage < 0.0 ) {
		Xdata->Cdata.storage = 0.0;
	}
	if ( Xdata->Cdata.lai <= 0.0 ) {
		Xdata->Cdata.lai = 0.0;
		return; //abort function execution, there is no canopy at this point
	}
	if ( Xdata->Cdata.height <= 0.0 ) {
		Xdata->Cdata.height = 0.0;
	}
	/*
	 * 1.1 compute the interception capacity [mm m-2]
	 * 1.1a Always new snow density as estimate of density in intercepted storage
	 */
	density_of_new_snow = SnLaws::compNewSnowDensity(hn_density, hn_density_model,
	                                                 *Mdata, *Xdata, Xdata->Cdata.temp, -1.);

	// 1.1b Determine interception capacity [mm] as function of density of intercepted snow/rain
	intcapacity = cn_IntCapacity(Mdata->ta, density_of_new_snow, Xdata->Cdata.lai);

	// 1.2 compute direct unload [mm timestep-1], update storage [mm]
	unload = cn_IntUnload(intcapacity, Xdata->Cdata.storage);

	if ( unload < 0.0 ) {
		prn_msg(__FILE__, __LINE__, "wrn", Mdata->date, "Negative unloading!!!");
		unload = 0.0;
	}
	Xdata->Cdata.storage -= unload;

	// 1.3 compute the interception [mm timestep-1] and update storage [mm]
	precipitation = Mdata->hnw;
	interception = cn_IntRate(intcapacity, Xdata->Cdata.storage, precipitation, Xdata->Cdata.direct_throughfall, Mdata->ta);

	Xdata->Cdata.storage += interception;

	// 1.4 compute the throughfall [mm timestep-1]
	throughfall = precipitation - interception + unload;

	Mdata->hnw = throughfall; // Please give the total amount for the time step

	// 2.1 prepare for canopy energy balance

	// Wetfraction update is moved to canopy energy balance loop  - use old value first
	wetfrac = Xdata->Cdata.wetfraction;

	/*
	 * Radiation Transmissivity
	 * (could possibly be a function of interception - but is constant for the moment)
	 * Firstly, transmissivity of diffuse (and longwave) radiation
	 */
	Xdata->Cdata.sigf = cn_CanopyTransmissivity(Xdata->Cdata.lai, Constants::pi / 2.0);

	// Secondly, transmissivity of direct solar radiation
	if ( Canopy::canopytransmission ) {
		sigfdirect = cn_CanopyTransmissivity(Xdata->Cdata.lai, Mdata->elev);
	} else {
		sigfdirect = Xdata->Cdata.sigf;
	}

	/*
	 * Reference Height [m above snow surface] for meteo input, at least 2 m above canopy height above snow surface
	 * 2006-03-01: introduced new definition of reference height above canopy for Alpine3D applications
	 * 2006-06-02: this should work also without soil data, since Xdata->Ground is initialized to 0.0
	*/
	if ( ALPINE3D ) {
		zref = MAX (2.0 + (Xdata->Cdata.height - (Xdata->cH - Xdata->Ground)), Xdata->Cdata.height + height_of_wind_val - (Xdata->cH - Xdata->Ground));
	} else {
		zref = MAX(2.0 + (Xdata->Cdata.height - (Xdata->cH - Xdata->Ground)), height_of_wind_val - (Xdata->cH - Xdata->Ground));
	}

	if ( (Xdata->cH - Xdata->Ground) > 0.03 ) {
		Mdata->z0 = roughness_length;
		z0m_ground = roughness_length;
	} else {
		z0m_ground = Xdata->BareSoil_z0;
		Mdata->z0 = Xdata->BareSoil_z0;
	}
	/*
	 * Turbulent Transport Coefficients:
	 * ch_canopy        = sensible heat
	 * ce_transpiration = latent heat from dry fraction, incl stomatal control.
	 * ce_interception  = latent heat from wet fraction
	 * ce_canopy        = latent heat total
	 * cm_canopy        = momentum (through the canopy, i.e to estimate wind speed below)
	*/
	cn_CanopyTurbulentExchange(*Mdata, zref, z0m_ground, wetfrac,
						  *Xdata, ch_canopy, ce_canopy, ce_transpiration, ce_interception, ce_condensation);


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
	for ( ebalitt = 0; ebalitt < 3; ebalitt++ ) {
		TC_OLD = Xdata->Cdata.temp; // Cdata.temp is updated in the iteration...

		// update ce_canopy as function of wetfraction
		ce_canopy = MAX (0.001 * ce_interception, ce_interception * wetfrac + ce_transpiration * (1.0 - wetfrac));
		ce_condensation = ce_interception * MAX (0.1, wetfrac);

		// canopy albedo
		canopyalb = cn_CanopyAlbedo(Mdata->ta, wetfrac);

		// compute properties r0 and r1 in eq (2) (and downward lw and sw for snowpack model)
		cn_LineariseNetRadiation(*Mdata, Xdata->Cdata, *Xdata, iswrac, rsnet, ilwrac, r0, r1,
							canopyalb, canopyclosuredirect, radfracdirect, sigfdirect, r1p);

		// compute properties h0 and h1 in eq (3)
		cn_LineariseSensibleHeatFlux(ch_canopy, Mdata->ta, h0, h1);

		// compute properties le0 and le1 in eq (4)
		cn_LineariseLatentHeatFlux(ce_canopy, Xdata->Cdata.temp, Mdata->rh*lw_SaturationPressure(Mdata->ta), le0, le1);

		/* final canopy energy balance */
		cn_CanopyEnergyBalance(h0, h1, le0, le1, Mdata->rh * lw_SaturationPressure(Mdata->ta),
						   ce_canopy, ce_condensation, r1p, 1. - Xdata->Cdata.direct_throughfall,
						   r0, r1, Xdata->Cdata.temp, RNCANOPY, HCANOPY, LECANOPY);

		/*
		 * Partition latent heat flux on interception and transpiration
		 * and correct energy balance for overestimated interception evaporation
		*/
		cn_CanopyEvaporationComponents(ce_canopy, ce_transpiration, &LECANOPY, Mdata->ta,
								 Mdata->rh * lw_SaturationPressure(Mdata->ta), Xdata->Cdata.storage,
								 M_TO_H(calculation_step_length), &CanopyEvaporation, &INTEVAP, &TRANSPIRATION,
								 &RNCANOPY, &HCANOPY, &Xdata->Cdata.temp, &r0, &r1, h0, h1, &LECANOPYCORR,
								 r1p, 1. - Xdata->Cdata.direct_throughfall ,wetfrac);

		newstorage = Xdata->Cdata.storage - INTEVAP;

		// wet surface fraction
		wetfrac = cn_CanopyWetFraction(intcapacity, newstorage);

		Xdata->Cdata.temp = (Xdata->Cdata.temp+TC_OLD)*0.5;
		wetfrac=(Xdata->Cdata.wetfraction+wetfrac)*0.5;
	} // End of Energy Balance Loop

	// Now REDUCE WaterContent in the Soil Elements --- Could also be part of WaterTransport.c
	if (snp_soil)
		cn_SoilWaterUptake(Xdata->SoilNode, TRANSPIRATION, &Xdata->Edata[0]);

	// final adjustment of interception storage due to evaporation
	Xdata->Cdata.storage = Xdata->Cdata.storage - INTEVAP;


	/*
	 * Preparation of output variables using += sign to allow for cumulated or averaged output
	 * (remember to reset these variables to 0 in Main.c before next integration step)
	*/
	// radiation above and below canopy
	cn_CanopyRadiationOutput(*Xdata, *Mdata, canopyalb, &iswrac, &rswrac, &iswrbc,&rswrbc,&ilwrac,
	                         &rlwrac,&ilwrbc,&rlwrbc,canopyclosuredirect,radfracdirect,sigfdirect);

	// longwave and shortwave radiation components
	Xdata->Cdata.iswrac += iswrac;
	Xdata->Cdata.rswrac += rswrac;
	Xdata->Cdata.iswrbc += iswrbc;
	Xdata->Cdata.rswrbc += rswrbc;
	Xdata->Cdata.ilwrac += ilwrac;
	Xdata->Cdata.rlwrac += rlwrac;
	Xdata->Cdata.ilwrbc += ilwrbc;
	Xdata->Cdata.rlwrbc += rlwrbc ;

	// Net longwave and shortwave radiation of canopy [W m-2]
	Xdata->Cdata.rlnet += RNCANOPY-rsnet;
	Xdata->Cdata.rsnet += rsnet;

	// atmospheric emissivity as seen by surface below canopy with regard to air temperature
	Mdata->ea = ilwrbc / (Constants::stefan_boltzmann * (Mdata->ta * Mdata->ta * Mdata->ta * Mdata->ta));

	// downward and reflected shortwave below canopy
	Mdata->rswr = rswrbc;
	Mdata->iswr = iswrbc;

	// Adjust friction velocity below canopy using the same reference height as in Meteo.c
	zref = MAX(0.5, height_of_wind_val - (Xdata->cH - Xdata->Ground));
	Mdata->ustar = 0.74 * log(zref / z0m_ground) / 0.4 / (Xdata->Cdata.ra + Xdata->Cdata.rs);

	// Canopy turbulent heat fluxes [W m-2]
	Xdata->Cdata.sensible += HCANOPY;
	Xdata->Cdata.latent += LECANOPY;
	Xdata->Cdata.latentcorr += LECANOPYCORR;

	// Canopy evaporative fluxes [mm]
	Xdata->Cdata.transp += TRANSPIRATION;
	Xdata->Cdata.intevap += INTEVAP;

	// Canopy mass fluxes [mm]

	Xdata->Cdata.interception += interception;
	Xdata->Cdata.throughfall += throughfall;
	Xdata->Cdata.snowunload += unload;

	// Canopy auxiliaries
	Xdata->Cdata.wetfraction = wetfrac;
	Xdata->Cdata.intcapacity += intcapacity;
	Xdata->Cdata.canopyalb += canopyalb;
	if ( Xdata->getNumberOfElements() > 0 ) {
		Xdata->Cdata.totalalb +=  cn_TotalAlbedo(canopyalb, Xdata->Cdata.sigf, Xdata->Albedo,
					Xdata->Cdata.direct_throughfall, canopyclosuredirect, radfracdirect, sigfdirect);
	} else {
		Xdata->Cdata.totalalb +=  cn_TotalAlbedo(canopyalb, Xdata->Cdata.sigf, Xdata->SoilAlb,
					Xdata->Cdata.direct_throughfall, canopyclosuredirect, radfracdirect, sigfdirect);
	}
}
