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

#include <snowpack/Radiation.h>

using namespace std;
using namespace mio;


/************************************************************
 * static section                                           *
 ************************************************************/
/**
 * Sun elevation threshold (rad) below which radiation
 * splitting is no longer working; ori Nora: 1.e-2
 */
const double Radiation::thresh_sun_elevation = 1.e-2; 

/************************************************************
 * non-static section                                       *
 ************************************************************/

Radiation::Radiation(const mio::Config& cfg)
{
	cfg.getValue("SW_MODE", "Snowpack", sw_mode);
	sw_mode %= 10;
}

double Radiation::ProjectToHorizontal(const double& slope_component, const double& ang_inc, 
                                      const double& sunx, const double& suny, const double& sunz)
{ // Project a radiation component that was computed on slope back to horizontal
	//HACK: this method does not work for angles of incidence = 90 deg
	//HACK: do it using local slope/azimuth in order to be safer
	const double rad_to_sun_vector = slope_component / cos( ang_inc );
	const double sun_vector_to_hor = sunz / sqrt( sunx*sunx + suny*suny + sunz*sunz );

	return ( rad_to_sun_vector * sun_vector_to_hor );
}

/**
 * @brief Compute the solar incidence (rad), i.e. the angle between the incident sun beam
 *   and the normal to the slope
 * @param sx
 * @param sy
 * @param slope_angle (deg)
 * @param azi (deg)
 * @param Psolar PositionSun
 */
void Radiation::angleOfIncidence(const double& sx, const double& sy, const double& slope_angle,
                                 const double& azi, PositionSun& Psolar)
{
	double cos_ang_inc;

	// TODO Is there a better way ?
	if ( slope_angle > 0. ) {
		// From Oke pp.345-346, use solar azimuth clockwise from N
		if ( Psolar.elev > 0. ) {
			cos_ang_inc = (cos(DEG_TO_RAD(slope_angle))*sin( Psolar.elev )) +
					(sin(DEG_TO_RAD(slope_angle))*cos( Psolar.elev )*cos( Psolar.azi_Ncw - DEG_TO_RAD(azi) ));
		} else {
			cos_ang_inc = 0.;
		}
	} else {
			// A3D: simplified according to Funk (1984 pp.139-140)
		cos_ang_inc = (sx * Psolar.sunx + sy * Psolar.suny + Psolar.sunz) /
				sqrt( sx * sx + sy * sy + 1. );
	}

	cos_ang_inc = MAX(0., cos_ang_inc);
	Psolar.ang_inc = acos( cos_ang_inc );
} // End of AngleOfIncidence

/**
 * @brief This function evaluates day numbers of the year
 * day_number:      1. == 1 Jan
 * day_number_equi: 1. == Spring equinox, that is, 20 March
 * @param *day_number
 * @param *day_number_equi
 * @param date_in
 */
void Radiation::computeDayNumbers(const mio::Date& date_in, double& day_number, double& day_number_equi)
{
	int    YYYY,MM,DD,HH,MI;
	double spring_equi;

	// Compute the day-of-the-year (DOY) from the Julian Date
	date_in.getDate(YYYY, MM, DD, HH, MI);
	Date year_end(YYYY-1, 12, 31, 0, 0, date_in.getTimeZone());
	day_number = (int)(date_in.getJulianDate() - year_end.getJulianDate());

	// Compute the day-of-the-year (DOY) from the Julian Date starting at spring equinox
	// spring equinox time in days from the beginning of the year (Bourges (1985))
	// spring_equi = 78.801 + 0.2422 * (YYYY - 1969) - (int) ( 0.25 * (YYYY - 1969) );
	// starting from spring equinox at the year 2000 (20.3. 8:30(UTC+1-> 79.3542)
	if ( YYYY < 2000 ) {
		spring_equi = 78.8010 + 0.2422 * (YYYY - 1969) - (int) ( 0.25 * (YYYY - 1969) );
	} else {
		spring_equi = 79.3542 + 0.2422 * (YYYY - 2000) - (int) ( 0.25 * (YYYY - 2000) );
	}

	// actual day number transformed that the time in days starts from spring equinox (Bourges (1985))
	if ( (day_number_equi = (int)(day_number - 0.5 - spring_equi)) < 0 ) {
		day_number_equi = (int)(day_number_equi + 365.2425);
		// because day_number_equi is lower than zero
		// for days before spring equinox in the year using the Gregorian calendar with
		// a mean year length of 365.2425 days
	}
} // End of ComputeDayNumbers

/**
 * @brief This function evaluates daily values of:
 * the correction due to the eccentricity of the earth's orbit
 * the solar declination
 * the equation of time
 * @param *Psolar
 * @param date_in
 * @return int
 */
void Radiation::computeSolarDailyParameters(const mio::Date& date_in, PositionSun& Psolar)
{
	double day_angle, day_angle_decl;
	double day_number, day_number_equi;

	computeDayNumbers(date_in, day_number, day_number_equi);
	// day_angle_decl in radians analog to Iqbal (1983); with day_number_equi as day_number
	// starting from spring equinox
	// using the Gregorian calendar mean number of days in a year = 365.2425
	day_angle_decl = DEG_TO_RAD( ((day_number_equi - 1.) * 360. / 365.2425) );

	// day_angle in radians (Iqbal (1983) p. 3)
	day_angle = DEG_TO_RAD( (day_number - 1.) * 360. / 365.2425 );

	// correction due to the eccentricity of the earth's orbit (Spencer (1971), Iqbal (1983), p.3)
	Psolar.ecc_corr = (1.000110 + 0.034221 * cos( day_angle ) + 0.001280 * sin( day_angle )
				+ 0.000719 * cos( 2. * day_angle ) + 0.000077 * sin( 2. * day_angle ));

	// solar declination after Bourges (1985) with a Fourier series approximation
	// and day number starting from 20.3. = 0 ...
	// and errors smaller than 0.02 degree
	// solar declination in degree
	Psolar.decl = 0.3723 + 23.2567 * sin ( day_angle_decl ) - 0.758   * cos ( day_angle_decl )
			+ 0.1149  * sin ( 2. * day_angle_decl ) + 0.3656  * cos ( 2. * day_angle_decl )
			- 0.1712  * sin ( 3. * day_angle_decl ) + 0.0201  * cos ( 3. * day_angle_decl );

	// transformation to radians
	Psolar.decl = DEG_TO_RAD( Psolar.decl );

	// equation of time in hours (Spencer (1971), Iqbal (1983), p.11)
	// 229.18 = 24*60/2PI (conversion in minutes)
	Psolar.eq_time = 229.18 * (0.000075 + 0.001868 * cos( day_angle ) - 0.032077 * sin( day_angle )
			- 0.014615 * cos( 2. * day_angle ) - 0.040849 * sin( 2. * day_angle )) / 60.;
} // End of ComputeSolarDailyParameters

/**
 * @brief This function evaluates hourly values of:
 *  - the solar elevation and zenith (rad)
 *  - the solar azimuth (rad)
 *  Note: We use a right handed system (y : N positive, x : E positive);
 * @param *Psolar PositionSun
 * @param local_time double
 * @param Lat double
 * @param Lon double
 * @return int
 */
void Radiation::computePositionSun(const double& local_time, const double& Lat, const double& Lon, PositionSun& Psolar)
{
	double azi;   // Solar azimuth in radians, clockwise from North

	// true SOLAR TIME or local apparent time in hours (Oke (1987), p.340 or Iqbal (1983), p.12):
	// local mean time                                                           + equation of time
	// local standard time (not daylight saving time) + longitude correction     + equation of time
	// Duffie (2006): Lon in degrees west: 0 < Lon < 360 -> (360 - Lst - (360 - Llo))
	// local standard time + 4min * (Lon_local - Lon_standameridian)          + equation of time
	// 4min for every degree, accounts for the difference between the local and standard meridians (longitude)
	// standard meridian: longitude of the center of my local time zone, 15*time difference to UTC
	// local standard time:
	// meteo input always has to be in local standard time (winter in Davos: UTC+1) to avoid confusion
	// with daylight saving time and for the computation of the hour_angle

	if ( Lon >= 0. ) {
		Psolar.solar_time = local_time
				+ (4. / 60.) * (Lon - (int)( Lon / 15. + 0.5 ) * 15.)
				+ Psolar.eq_time;
	} else {
			// for the western hemisphere
		Psolar.solar_time = local_time + (4. / 60.) * (Lon - (int)( Lon / 15. - 0.5 ) * 15.)
					+ Psolar.eq_time;
	}

	// hour angle in radians (Oke (1987), p.340 or Iqbal (1983), p.15)
	// zero at noon, positive before solar noon and negative afterwards
	// angular displacement of the sun east or west of the local(!) meridian
	Psolar.hr_angle = DEG_TO_RAD( 15. * (12. - Psolar.solar_time) );

	// solar zenith angle in radians Oke (1987), p.339
	Psolar.zen = acos( sin( Psolar.decl ) * sin( DEG_TO_RAD( Lat ) )
			+ cos( Psolar.decl ) * cos( DEG_TO_RAD( Lat ) ) * cos( Psolar.hr_angle ) );
	// solar elevation in radians
	Psolar.elev = Constants::pi / 2. - Psolar.zen;

	// solar azimuth in radians, measured clockwise from North for the meteo station (Oke (1987), p.339)
	azi = (sin( Psolar.decl ) * cos( DEG_TO_RAD( Lat ) )
			- cos( Psolar.decl ) * sin( DEG_TO_RAD( Lat ) )
			* cos( Psolar.hr_angle )) / sin( Psolar.zen );

	Psolar.azi_Ncw = acos( azi );

	if ( Psolar.solar_time > 12. ) {
		Psolar.azi_Ncw = 2. * Constants::pi - Psolar.azi_Ncw;
	}

	// Convert to angle measured from South, counterclockwise (rad)
	if ( Psolar.azi_Ncw <= Constants::pi ) {
		Psolar.azi_Sacw = Constants::pi - Psolar.azi_Ncw;
	} else {
		Psolar.azi_Sacw = 3. * Constants::pi - Psolar.azi_Ncw;
	}

	// sun vector S
	// derived as shown in Funk (1984) p. 107, but for a y-coordinate increasing northwards
	Psolar.sunx =  sin( Psolar.azi_Sacw ) * cos( Psolar.elev ); 	// left handed coordinate system
	// inverse y-coordinate because of reference system chosen
	Psolar.suny = -cos( Psolar.azi_Sacw ) * cos( Psolar.elev );
	Psolar.sunz =  sin( Psolar.elev );
} // End of ComputePositionSun

/**
 * @brief Computes the hourly potential direct and diffuse irradiance on an arbitrary inclined
 * surface according to Iqbal (1983), p.188 (Bird and Hulstrom (1980, 1981) model)
 * The top of atmosphere radiation Rdata->toa_h is multiplied with transmittances according
 * to the atmospheric composition (e.g. water vapor, ozon); furthermore the position of the
 * sun, that is, the angle of incidence is taken into account.
 * @param *Rdata RadiationData
 * @param *Psolar const PositionSun
 * @param mean_alb const double
 * @param altitude const double
 * @param pressure const double
 * @param rh const double
 * @param ta const double
 * @return int
 */
void Radiation::computePotentialRadiation(const PositionSun& Psolar, const double& mean_alb, const double& altitude,
                                          const double& pressure, const double& rh, const double& ta, RadiationData& Rdata)
{
	double e_stern;               // saturation water vapour pressure
	double ma;                    // actual optical air mass
	double mr;                    // relative optical air mass
	double taur;                  // broadband transmittance by Rayleigh scattering
	double tauoz;                 // broadband transmittance by ozone
	const double olt = 0.32;      // ozone layer thickness (cm) U.S.standard = 0.34 cm
	double u3;                    // ozone relative optical path length
	double taug;                  // broadband transmittance by mixed gases
	double tauw;                  // broadband transmittance by water vapour
	double precw;                 // precipitable water-vapor thickness
	double u1;                    // pressure-corrected relative optical path length of precipitable water
	double taua;                  // broadband transmittance by aerosols due to total attenuation
                                //   (scattering+absorption)
	double tauaa;                 // broadband transmittance by aerosols due to absorption only
	const double w0 = 0.9;        // fraction of energy scattered to total attenuation by aerosols
				//   (Bird and Hulstrom(1981))
	double tauas;                 // broadband transmittance by aerosols due to scattering only
	double alb_sky;               // albedo of the cloudless sky
	const double fc = 0.84;       // fraction of forward scattering to total scattering (Bird and Hulstrom(1981))
	const double alpha = 1.3;     // wavelength exponent (Iqbal(1983) p.118)
	const double beta = 0.03;
	double ka, ka1, ka2;          // coefficients of aerosol transmittance (beta<0.4)
	double beta_z;                // correction term for increased transmittance with altitude

 	 // relative optical air mass Kasten and Young (1989)
	mr = 1. / (cos( Psolar.zen ) + 0.50572 * pow( 96.07995 - RAD_TO_DEG( Psolar.zen ),-1.6364 ));

	// actual air mass: because mr is applicable for standard pressure
	// it is modified for other pressures (in Iqbal (1983), p.100)
	// pressure in Pa
	ma = mr * (pressure / 101325.);

	// the equations for all the transmittances of the individual atmospheric constituents
	// are from Bird and Hulstrom (1980, 1981) and can be found summarized in Iqbal (1983) 
	// on the quoted pages

  	// broadband transmittance by Rayleigh scattering (Iqbal (1983), p.189)
	taur = exp( -0.0903 * pow( ma,0.84 ) * (1 + ma - pow( ma,1.01 )) );

  	// broadband transmittance by ozone (Iqbal (1983), p.189)
	u3 = olt * mr;
	tauoz = 1. - (0.1611 * u3 * pow( 1. + 139.48 * u3,-0.3035 ) -
			0.002715 * u3 * pow( 1. + 0.044  * u3 + 0.0003 * u3 * u3,-1 ));

  	// broadband transmittance by uniformly mixed gases (Iqbal (1983), p.189)
	taug = exp( -0.0127 * pow( ma,0.26 ) );

  	// broadband transmittance by water vapor (in Iqbal (1983), p.189):

  	// saturation vapor pressure in Pa 
	e_stern = lw_SaturationPressure(ta);

	// Leckner (1978); pressure and temperature correction not necessary since it is
	// included in its numerical constant (in Iqbal (1983), p.94)
	precw = 0.493 * rh * e_stern / ta;

  	// pressure corrected relative optical path length of precipitable water (Iqbal (1983), p.176)
	u1 = precw * mr;

  	// broadband transmittance by water vapor (in Iqbal (1983), p.189)
	tauw = 1 - 2.4959 * u1 * (1. / (pow( 1.0 + 79.034 * u1,0.6828 ) + 6.385 * u1));

	// broadband total transmittance by aerosols (in Iqbal (1983), pp.189-190)
	// using Angstroem's turbidity formula Angstroem (1929, 1930) for the aerosol thickness
	// in Iqbal (1983), pp.117-119
	// aerosol optical depth at wavelengths 0.38 and 0.5 micrometer
	ka1 = beta * pow( 0.38,-alpha );
	ka2 = beta * pow( 0.5 ,-alpha );

	  // broadband aerosol optical depth:
	ka  = 0.2758 * ka1 + 0.35 * ka2;

  	// total aerosol transmittance function for the two wavelengths 0.38 and 0.5 micrometer:
	taua = exp( -pow( ka,0.873 ) * (1. + ka - pow( ka,0.7088 )) * pow( ma,0.9108 ) );

  	// Iqbal (1983) p. 190
	tauaa = 1. - (1. - w0) * (1. - ma + pow( ma,1.06 )) * (1. - taua);

	// broadband transmittance function due to aerosols scattering only 
	// Iqbal (1983) p. 146 (Bird and Hulstrom (1981))
	tauas = taua / tauaa;

  	// cloudless sky albedo Bird and Hulstrom (1980, 1981) (in Iqbal (1983) p. 190)
	  // alb_sky = alb_rayleighscattering + alb_aerosolscattering 
	alb_sky = 0.0685 + (1. - fc) * (1. - tauas);

	// direct normal solar irradiance in range 0.3 to 3.0 micrometer (Iqbal (1983) ,p.189)
	// 0.9751 is for the wavelength range ??
	// Bintanja (1996) (see Corripio (2002)) introduced a correction beta_z for increased 
	// transmittance with altitude that is linear up to 3000 m and than fairly constant up to 5000 - 6000 m
	if( altitude < 3000. ) {
		beta_z = 2.2 * 1.e-5 * altitude;
	} else {
		beta_z = 2.2 * 1.e-5 * 3000.;
	}
	Rdata.pot_dir = 0.9751 * Constants::solcon * Psolar.ecc_corr * (taur * tauoz * taug * tauw * taua + beta_z);
	Rdata.pot_dir *= cos( Psolar.ang_inc );

	// Diffuse radiation from the sky
	// Rayleigh-scattered diffuse radiation after the first pass through atmosphere (Iqbal (1983), p.190)
	Rdata.pot_diffsky = 0.79 * Rdata.toa_h * tauoz * taug * tauw * tauaa
			* 0.5 * (1. - taur ) / (1. - ma + pow( ma,1.02 ));

  	// aerosol scattered diffuse radiation after the first pass through atmosphere (Iqbal (1983), p.190)
	Rdata.pot_diffsky += 0.79 * Rdata.toa_h * tauoz * taug * tauw * tauaa
			* fc  * (1. - tauas) / (1. - ma + pow( ma,1.02 ));

  	// multiple reflected diffuse radiation between surface and sky (Iqbal (1983), p.154)
	Rdata.pot_diffsky += (Rdata.pot_diffsky + Rdata.pot_dir * cos( Psolar.zen )) * mean_alb * alb_sky /
			(1. - mean_alb * alb_sky);
} // End of ComputePotentialRadiation

/**
 * @brief Evaluate the coefficient Psolar->Md to split incoming global solar radiation into direct
 * and diffuse components. Splitting is based on "clearness of the sky", i.e. the ratio of
 * measured incoming global radiation to top of the atmosphere radiation toa_h.g
 * Erbs et al., Iqbal p.269e
 * @param *Rdata RadiationData
 * @param *Psolar const PositionSun
 * @param thresh_elev const double
 * @return int
 */
void Radiation::computeSplittingCoefficient(const PositionSun& Psolar, const double& thresh_elev, RadiationData& Rdata)
{
	//The clearsky index Mt is not usable for very small positive solar elevations:
	//   Mt is getting too large (>> 1) pretending a clear sky day,
	//   but actually almost all radiation is diffuse!
	if ( Psolar.elev < thresh_elev ) {
		Rdata.Md = 1.0;
		return;
	}
	
	// clear sky index (ratio global measured to top of atmosphere radiation)
	const double Mt = Rdata.global_hor / Rdata.toa_h;
	
	// diffuse fraction: hourly ratio of diffuse to global radiation incident on a horizontal surface
	// splitting according to a combination of Reindl et al.(1990)'s models (Mt-model and Mt&Psolar.elev-model):
	if ( Mt >= 0.78 ) {				// Mt in [0.78;1] -> clear day
		Rdata.Md = 0.147;
	} else {
		if ( Mt <= 0.3 ) {			// Mt in [0;0.3] -> overcast
			Rdata.Md = 1.02 - 0.248 * Mt;
			Rdata.Md = MIN(1.0, Rdata.Md);
		} else {					// Mt in ]0.3;0.78[ -> cloudy
			Rdata.Md = 1.4 - 1.749 * Mt + 0.177 * sin( Psolar.elev );
			Rdata.Md = MIN(0.97, Rdata.Md);
			Rdata.Md = MAX(0.1 , Rdata.Md);
		}
	}
} // End of ComputeSplittingCoefficient

/**
 * @brief Used in Snowpack for projecting direct SW-radiation on (virtual) slopes (see snowpack_opera)
 * => Direct solar radiation received on horizontal surface is projected onto slope while
 * diffuse component is assumed to be isotropic.
 * Global solar radiation received on slope is then the sum of projected direct solar
 * radiation and diffuse radiation.
 * Use is made of equation A1.8 found in Oke, p. 346, but for sun elevation larger than
 * 9 deg only (see caveat in Iqbal on p. 77)
 * @param *Mdata CurrentMeteo
 * @param *Rdata RadiationData
 * @param *Psolar const PositionSun
 * @param Alb const double
 * @return int
 */
void Radiation::projectRadiationOnSlope(const PositionSun& Psolar, const double& Alb, 
                                        CurrentMeteo& Mdata, RadiationData& Rdata)
{
	if ( Psolar.elev > 0.157 ) { // 9 deg = 0.157 rad
		Rdata.dir_slope = Rdata.dir_hor*cos( Psolar.ang_inc )/sin( Psolar.elev );
	} else {
		Rdata.dir_slope = Rdata.dir_hor*cos( Psolar.ang_inc )/sin( 0.157 );
	}
	Rdata.dir_slope = MAX(0., Rdata.dir_slope);

  	// Global radiation received on slope
	if ( Rdata.global_hor > 0. ) {
		Mdata.iswr = MIN(Rdata.dir_slope + Rdata.diffsky, Constants::solcon);
		Mdata.rswr = Alb*Mdata.iswr;
	} else {
		Mdata.iswr = 0.;
		Mdata.rswr = 0.;
	}
} // End ProjectRadiationOnSlope

/**
 * @brief Compute global, direct and diffuse radiation on flat field, including potential radiation
 * ATTENTION: This function should NOT be called by Alpine3D !!!
 * @param *Mdata CurrentMeteo
 * @param *Xdata Profile
 * @param *Psolar PositionSun
 * @param *Rdata RadiationData
 * @return int
 */
void Radiation::flatFieldRadiation(const SnowStation& Xdata, CurrentMeteo& Mdata,
                                   PositionSun& Psolar, RadiationData& Rdata)
{
	// SOLAR PARAMETERS have to be computed
	double local_time = Mdata.date.getJulianDate() - floor(Mdata.date.getJulianDate()) - 0.5; // in julian days
	// Compute daily solar parameters once a day (shortly after midnight)
	// HACK what's the purpose of the condition on ecc_corr?
	if ((local_time > 0.005 && local_time < 0.016) || !(Psolar.ecc_corr > 0.)) {
		computeSolarDailyParameters(Mdata.date, Psolar);
	}
	// Hourly position of the sun
	computePositionSun(D_TO_H(local_time), Xdata.meta.position.getLat(), Xdata.meta.position.getLon(), Psolar);
	// Angle of incidence
	angleOfIncidence(0., 0., Xdata.meta.getSlopeAngle(), Xdata.meta.getAzimuth(), Psolar);

	// Incoming global solar radiation on horizontal surface
	if (sw_mode == 1) {
		Rdata.global_hor = Mdata.rswr/Xdata.cAlbedo;
	} else {
		Rdata.global_hor = Mdata.iswr;
	}

	// Top of atmosphere radiation
	if ((Rdata.toa_h = Constants::solcon * Psolar.ecc_corr * cos( Psolar.zen )) < 0.) {
		Rdata.toa_h = Constants::nodata;
	}
	// Compute potential radiation only for Psolar.elev > THRESH_SUN_ELEVATION ...
	if (Psolar.elev >= Radiation::thresh_sun_elevation) {
		computePotentialRadiation(Psolar, Xdata.cAlbedo, Xdata.meta.position.getAltitude(), lw_AirPressure(Xdata.meta.position.getAltitude()), Mdata.rh, Mdata.ta, Rdata);
	} else {// ... because radiation is only diffuse otherwise
		Rdata.pot_dir = 0.;
		Rdata.pot_diffsky = Rdata.global_hor;
	}

	// Split flat field radiation into direct and diffuse components
	computeSplittingCoefficient(Psolar, Radiation::thresh_sun_elevation, Rdata);
	if ((Rdata.global_hor > 0.) && (Rdata.pot_dir > 0.)) {
		Rdata.dir_hor = (1. - Rdata.Md)*Rdata.global_hor;
		Rdata.diffsky = Rdata.Md*Rdata.global_hor;
	} else {
		if (Rdata.global_hor > 0.) {
			Rdata.dir_hor = 0.;
			Rdata.diffsky = MAX(Rdata.Md*Rdata.global_hor, Rdata.pot_diffsky);
		} else {
			Rdata.dir_hor = 0.;
			Rdata.diffsky = 0.;
		}
	}

	// Reset global radiation values
	if (sw_mode == 1) {
		Mdata.rswr = (Rdata.dir_hor + Rdata.diffsky)*Xdata.cAlbedo;
	} else {
		Mdata.iswr = Rdata.dir_hor + Rdata.diffsky;
	}
} // End FlatFieldRadiation

/**
 * @brief Compute radiation and precipitation on slopes (not for A3D !!!)
 * This is meant to be used for Snowpack:
 *  radiation and precipitation are measured on a flat field, radiation split into direct and
 *  diffuse components, and then projected onto a slope
 *  Note: Assume the reference station is not shaded (horizont) at any time!
 *        Otherwise one should use a mask (not implemented yet TODO)
 * ATTENTION: This function should NOT be called by Alpine3D !!!
 * @param Mdata CurrentMeteo
 * @param Xdata Profile
 * @param Sdata SurfaceFluxes
 * @param Psolar PositionSun
 * @param Rdata RadiationData
 */
void Radiation::radiationOnSlope(const SnowStation& Xdata, CurrentMeteo& Mdata, SurfaceFluxes& Sdata,
                                 PositionSun& Psolar, RadiationData& Rdata)
{
	if ( Xdata.meta.getSlopeAngle() > Constants::min_slope_angle ) {
		// Angle of incidence
		angleOfIncidence(0., 0., Xdata.meta.getSlopeAngle(), Xdata.meta.getAzimuth(), Psolar);
		// Project radiation
		projectRadiationOnSlope(Psolar, Xdata.cAlbedo, Mdata, Rdata);
	} else {
		Rdata.dir_slope = Rdata.dir_hor;
	}
	
	// Assign radiation values to Sdata
	Sdata.sw_hor  += Rdata.global_hor;
	Sdata.sw_dir  += Rdata.dir_slope;
	Sdata.sw_diff += Rdata.diffsky;
} // End ProjectOnSlope

/*
 * End of Radiation.cc
 */
