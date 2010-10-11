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
 * @file Meteo.c
 * @author Michael Lehning and others
 * @version 9.x
 * @date -
 * @bug -
 * @brief Calculates missing meteorological information such as friction velocity and roughness length
 * - 29.10.2002: Michael Lehning implements Micromet()
 * - 15.03.2005: Andy and Michi implement stability correction for turbulent fluxes in the hope
 *               that this will also improve the little bit too strong melting of the version 8.1
 */

#include <snowpack/Meteo.h>

Meteo::Meteo(const mio::Config& i_cfg) : cfg(i_cfg), canopy(cfg) 
{
	/**
	 * @brief Defines the way to deal with atmospheric stability:
	 * -    0: Standard MO iteration with Paulson and Stearns & Weidner (can be used with BC_CHANGE=0)
	 * -    1: Assume neutral stratification. Should be used with BC_CHANGE=1, i.e., Dirichlet bc but also
	 *         recommended with Neumann b.c., i.e., BC_CHANGE=0
	 * - (-1): Simplified Richardson number stability correction
	 */
	neutral = cfg.get("NEUTRAL", "Parameters");

	/**
	 * @brief Initial estimate of the roughness length for the site; will be adjusted iteratively. \n
	 * Default value and operational mode: 0.002 m
	 */
	ROUGHNESS_LENGTH = cfg.get("ROUGHNESS_LENGTH", "Parameters");

	/**
	 * @brief Defines whether the canopy model is used \n
	 * NOTE: OUT_CANOPY must also be set to dump canopy parameters to file; see Constants_local.h
	 */
	useCanopyModel = cfg.get("CANOPY", "Parameters");


	/**
	 * @brief Define the heights of the meteo measurements above ground (m) \n
	 * Required for surface energy exchange calculation and for drifting and blowing snow.
	 */
	HEIGHT_OF_WIND_VALUE = cfg.get("HEIGHT_OF_WIND_VALUE", "Parameters");

	research_mode = cfg.get("RESEARCH", "Parameters");
}


/**
 * @brief Projects precipitations and snow height perpendicular to slope
 * @param precips 
 * @param hs Height of snow (m)
 * @param SlopeAngle Slope angle
 */
void Meteo::projectPrecipitations(const double& SlopeAngle, double& precips, double& hs)
{
	precips *= cos(SlopeAngle);
	hs  *= cos(SlopeAngle);
} // End of ProjectPrecipitations

/**
 * @brief Make an iteration to find z0 and ustar at the same time
 * @param Mdata
 * @param Xdata
 */
void Meteo::MicroMet(const SN_STATION_DATA& Xdata, SN_MET_DATA& Mdata)
{
	int e, iter = 1, max_iter = 100;
	const double eps1 = 1.e-3, eps2 = 1.e-5;
	double ustar, z0 = ROUGHNESS_LENGTH, zref, a2 = 0.16 , vw, z0_old, ustar_old;
	double d_pump; // Wind pumping displacement depth (m)

	// New variables for stability correction
	double psi_m = 0., psi_s = 0., Tstar = 0.;
	double tss, tss_v, ta_v, dummy, p0, sat_vap, LH;
	double z_ratio = 1., stab_ratio = 0.;
	double Ri;  // Richardson number for simple stability correction

	// Initialize Virtual Temperatures for Stability
	tss = Xdata.Ndata[Xdata.getNumberOfElements()].T;

	// Ideal approximation of pressure and vapor pressure
	p0 = lw_AirPressure(Xdata.Alt);
	if (Mdata.ta > MELTING_TK) {
		LH = LH_VAPORIZATION;
	} else {
		LH = LH_SUBLIMATION;
	}

	sat_vap = lw_SaturationPressure(Mdata.ta);
	
	ta_v = Mdata.ta * (1. + 0.377 * sat_vap / p0);
	tss_v = tss * (1. + 0.377 * sat_vap / p0);

	/*
	 * Now start the real thing - iteratively determining stability and possibly adjusting z0 to
	 * drifting snow and ventilation
	*/
	e = Xdata.getNumberOfElements();
	vw = MAX(0.3, Mdata.vw);
	// Adjust for snow height
	zref = MAX (0.5, HEIGHT_OF_WIND_VALUE - (Xdata.cH - Xdata.Ground));
	
	// In case of ventilation ...
	if ( WIND_PUMP ) {
		d_pump = lwsn_WindPumpingDisplacement(Xdata);
	} else {
		d_pump = 0.;
	}

	// Iterate to find atmospheric stability
	// initial guess (neutral)
	ustar = 0.4 * vw / log((zref - d_pump) / z0);
	do {
		iter++;
		if ( iter > max_iter ) {
			Mdata.z0 = z0 = ROUGHNESS_LENGTH;
			Mdata.ustar = 0.4 * vw / log((zref - d_pump) / z0);
			Mdata.psi_s = 0.;
			prn_msg ( __FILE__, __LINE__, "wrn", Mdata.date.getJulianDate(), "Stability correction did not converge (azi=%.0lf, slope=%.0lf) --> assume neutral", RAD_TO_DEG(Xdata.SlopeAzi), RAD_TO_DEG(Xdata.SlopeAngle));
			return;
		}
		ustar_old = ustar;
		z0_old = z0;
		// Update z0
		z0 = 0.9 * z0_old + 0.1 * (a2 * ustar*ustar / 2. / Constants::g);
		z_ratio = log((zref - d_pump) / z0);
		
		// Prepare Values for Richardson
		if ( neutral < 0 ) { // Switch for Richardson
			Ri = Constants::g / tss_v * (ta_v - tss_v) * zref / vw / vw;
			if ( Ri < 0.2 ) {// neutral and unstable
				stab_ratio = Ri;
			} else {
				stab_ratio = Ri/(1.-5.*Ri);
			}			
			if ( Ri < 0. ) { // unstable
				stab_ratio = Ri;
				dummy = pow((1. - 15. * stab_ratio), 0.25);
				psi_m = log((0.5 * (1 + dummy*dummy)) * (0.5 * (1 + dummy)) * (0.5 * (1 + dummy))) - 2. * atan(dummy) + 0.5 * Constants::pi;
				psi_s = 2. * log(0.5 * (1 + dummy*dummy));
			} else if ( Ri < 0.1999 ) { // stable
				stab_ratio = Ri / (1. - 5. * Ri);
				psi_m = psi_s = -5. * stab_ratio;
			} else {
				stab_ratio = Ri / (1. - 5. * 0.1999);
				psi_m = psi_s = -5. * stab_ratio;
			}

			// Calculate ustar
			ustar = 0.4 * vw / (z_ratio - psi_m);
		} else if ( neutral == 0 || (!research_mode && (Mdata.tss > 273.) && (Mdata.ta > 277.)) ) { // MO Iteration
			// Calculate ustar
			ustar = 0.4 * vw / (z_ratio - psi_m);
			// Calculate Tstar
			Tstar = 0.4 * (tss_v - ta_v) / (z_ratio - psi_s);
			// Calculate Stability Parameter
			stab_ratio = -0.4 * zref * Tstar * Constants::g / (tss * ustar*ustar);
		
			if ( stab_ratio > 0. ) { // stable
				/* Stearns & Weidner, 1993 */
				dummy = pow((1. + 5. * stab_ratio), 0.25);
				psi_m = (log(1. + dummy) * log(1. + dummy) + log(1. + dummy*dummy) - 1. * atan(dummy) - 0.5 * dummy*dummy*dummy + 0.8247); // Original 2.*atan(dummy) - 1.3333
				// Launiainen and Vihma, 1990
				// psi_m = -17*(1.-exp(-0.29*stab_ratio));
			
				// Stearns & Weidner, 1993 for scalars
				dummy = sqrt(1. + 5. * stab_ratio);
				psi_s = (log(1. + dummy) * log(1. + dummy) - 1. * dummy - 0.3 * dummy*dummy*dummy + 1.2804); // Original 2.*dummy - 0.66667*...

				// Holtslag and DeBruin (1988) prepared from Ed Andreas
				// psi_m = psi_s = -(0.7*stab_ratio + 0.75*(stab_ratio - 14.28)*
				// exp(-0.35*stab_ratio) + 10.71);
			} else {
				// Stearns & Weidner, 1993 - Must be an ERROR somewhere
				//	dummy = pow((1.-15.*stab_ratio),0.25);
				//	psi_m = log(1.-dummy)*log(1.-dummy) + log(1.+dummy*dummy) - 2.*atan(dummy) -
				// 	        -1. + dummy - 0.5086;
				// Paulson - the original
				dummy = pow((1. - 15. * stab_ratio), 0.25);
				psi_m = 2. * log(0.5 * (1. + dummy)) +  log(0.5 * (1. + dummy*dummy)) - 2. * atan(dummy) + 0.5 * Constants::pi;
			
				// Stearns & Weidner, 1993 for scalars
				dummy = pow((1. - 22.5 * stab_ratio), 0.33333);
				psi_s = pow(log(1. + dummy + dummy*dummy), 1.5) - 1.732 * atan(0.577 * (1. + 2. * dummy)) + 0.1659;
			}
		} else { // NEUTRAL
			psi_m = 0.;
			psi_s = 0.;
		} 
	} while ( (fabs(ustar_old - ustar) > eps1) && (fabs(z0_old - z0) > eps2) );
	// Save the values in the global sn_Mdata data structure to use it later
	Mdata.ustar = ustar;
	Mdata.z0 = z0;
	Mdata.psi_s = psi_s;
	if ( (log(zref / z0) - psi_s) < 0.01 ) {
		psi_s = log(zref / z0) - 0.01; // Prevent contragradient fluxes
	}
} // End MicroMet

/**
 * @brief
 * \li with CANOPY set:
 * 		In case of an existing canopy, call canopy routine, which calculates precipitation, radiation,
 * 		friction velocity and reference temperature for the surface below the canopy.
 * 		Note that solar radiation may change also in dg_cn_Canopy(). \n
 * 		- Mdata->iswr  incoming global solar radiation (direct + diffuse), adapted to canopy
 * 		- Mdata->rswr  reflected global solar radiation (diffuse), adapted to canopy
 * 		- Mdata->ustar friction velocity, adapted to canopy
 * 		- Mdata->z0    roughness length, adapted to canopy
 * 		- Mdata->ea    atmospheric emissivity below canopy, i.e., to give correct
 * 		               longwave radiation as function of air temperature, however
 * 		               modified to include effect of canopy
 * \li without canopy (CANOPY is not set):
 * 		For bare soil as well as snowed-in canopy or some other problems, compute the roughness
 * 		length z0, the friction velocity ustar as well as the atmospheric stability correction
 * 		psi_s for scalar heat fluxes
 * 		- Mdata->ustar friction velocity
 * 		- Mdata->z0    roughness length
 * 		- psi_s        stability correction for scalar heat fluxes
 * @param *Mdata
 * @param *Xdata
 */
void Meteo::calculateMeteo(SN_MET_DATA *Mdata, SN_STATION_DATA *Xdata)
{
	if (useCanopyModel)
		canopy.runCanopyModel(Mdata, Xdata, ROUGHNESS_LENGTH, HEIGHT_OF_WIND_VALUE);

	if (!useCanopyModel || Xdata->Cdata.zdispl < 0.)
		MicroMet(*Xdata, *Mdata);
}

