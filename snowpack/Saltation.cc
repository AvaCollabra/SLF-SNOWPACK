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

#include <snowpack/Saltation.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/
/**
 * Height relative to maximum jump height taken
 * as hs, i.e. height where the concentration
 * is determined as BC for suspension
 */
const double Saltation::hs_frac = 1.0;

const double Saltation::karman = 0.4; ///< Von Karman constant
const double Saltation::elas = 0.5; ///<Coefficient of Elasticity 0.5

const double Saltation::angle_ej = 25.; // (deg)
const double Saltation::ratio_ve_ustar = 3.1; ///< Original Value by Judith: 2.9

const int Saltation::strong = 1;
const int Saltation::weak = 0;

///Saltation Z0 Preliminary values will be changed later Judith original: 0.00098
const double Saltation::z0_salt = 0.0017;
const double Saltation::salt_height = 0.07;

/*
 * The follwing const variables are no longer needed, their values shall remain documented:
 * const double Saltation::c_red = 1.0; //Defines red. of saltation conc. in suspension sol.
 * const double Saltation::grain_size = 0.000680;
 * const double Saltation::tau_thresh = 0.094; //Original Value by Judith: 0.094
 * const double Saltation::z0 = 0.01; //Wind Field Z0 - includes larger surface features
 */

/************************************************************
 * non-static section                                       *
 ************************************************************/

Saltation::Saltation(const mio::Config& cfg) 
{
	//Use Doorschots Saltation model instead of Sorenson (much slower)
	cfg.getValue("DOORSCHOT", "SnowpackAdvanced", doorschot);
}

/**
 * @brief Returns the wind profile
 * @param z
 * @param tauA
 * @param tauS
 * @param z0
 * @param u_start
 * @param slope_angle (deg)
 * @return u
 */
double Saltation::sa_vw(const double& z, const double& tauA, const double& tauS, const double& z0, 
                        const double& u_start, const double& slope_angle)
{
	double u, B, hs;
	
	hs = z0 / 100 + 0.6 * pow(u_start, 2) / (2. * Constants::g * cos(DEG_TO_RAD(slope_angle)));
	//  hs = DSQR(0.6*u_start)/(2.*Constants::g * cos(DEG_TO_RAD(slope_angle)));
	B = (tauS - tauA) / log(hs / z0);
	if (B > 0.0) {
		u = (pow(tauA + B * log(z / z0), 1.5) - pow(tauA, 1.5))
			/ (1.5 * B * Saltation::karman * sqrt(Constants::density_air));
	} else {
		u = sqrt(tauS / Constants::density_air) / Saltation::karman * log(z / z0);
	}
	
	return (u);
}

/**
 * @brief Returns the wind profile (alternative)
 * @param z
 * @param tauA
 * @param tauS
 * @param z0
 * @param u_start
 * @param slope_angle (deg)
 * @return u
 */
double Saltation::sa_vw2(const double& z, const double& tauA, const double& tauS, const double& z0, 
												 const double& u_start, const double& slope_angle)
{
	double u, hs, r, ustar, ustarz, dudz, z_act;
	double dz = 0.00002;
	
	u = 0; z_act = z0; 
	
	//  hs = z0 + DSQR(0.6*u_start)/(2.*Constants::g * cos(DEG_TO_RAD(slope_angle)));
	hs = pow(0.6 * u_start, 2) / (2. * Constants::g * cos(DEG_TO_RAD(slope_angle)));
	
	r = tauA / tauS;
	ustar = sqrt(tauS / Constants::density_air);
	//  fprintf(stdout, "\n z_act=%lf, r=%lf, ustar=%l",z_act,r,ustar);
	
	while (z_act < z) {
		z_act += dz; 
		//    fprintf(stdout, "\n z_act=%lf",z_act);
		ustarz = ustar * (1 - (1 - sqrt(r)) *  exp(-z_act / hs));
		dudz = ustarz / Saltation::karman / z_act;
		u += dudz * dz; 
	}
	//    fprintf(stdout, "\n z_act=%lf, ustarz=%lf, dudz=%lf, u=%lf",z_act,ustarz,dudz,u); }
	
	return(u);
}

/**
 * @brief Computes a trajectory and returns the characteristics
 * @param u0
 * @param angle_e_rad (rad)
 * @param slope_angle (deg)
 * @param dg
 * @param tauA
 * @param tauS
 * @param z0
 * @param ubar
 * @param u_i
 * @param angle_i_rad (rad)
 * @param t_i
 * @param z_max
 * @return bool
 */
bool Saltation::sa_Traject(const double& u0, const double& angle_e_rad, const double& slope_angle, const double& dg,
                          const double& tauA, const double& tauS, const double& z0, 
                          double& ubar, double& u_i, double& angle_i_rad, double& t_i, double& z_max)
{
	int nt = 0;
	const double DT = 0.0005;
	const double vis = 1.74e-5;
	double x, z, xdot, zdot, xdotdot, zdotdot, Ur, u;
	double Cd, Re;
	
	// Initialize velocities of particle and position
	xdot = u0 * cos(angle_e_rad);
	zdot = u0 * sin(angle_e_rad);
	x = xdot * DT;
	z = z0 + zdot * DT; 
	
	ubar = xdot;
	z_max = z0;
	
	if (ubar == 0.0) return false;
	
	// Make a simple forward time integration
	do {
		// Wind velocity
		u = sa_vw(z, tauA, tauS, z0, u0*sin(angle_e_rad), slope_angle);
		//  fprintf(stdout, "u=%lf \n",u);
		// Relative velocity
		Ur = sqrt(pow((xdot - u), 2) + pow(zdot, 2));
		// Reynolds number and Drag coefficient
		Re = dg * Ur / vis;
		Cd = 24. / Re + 4. / pow(Re,0.33);
		Cd = 24. / Re + 6 / (1 + sqrt(Re)) + 0.4; 
		//    fprintf(stdout, "\n Cd=%lf",Cd);
		
		// Accelerations
		xdotdot = -0.75 * Constants::density_air / Constants::density_ice * Ur / dg * Cd * (xdot - u) - Constants::g * sin(DEG_TO_RAD(slope_angle));
		zdotdot = -0.75 * Constants::density_air / Constants::density_ice * Ur / dg * Cd * zdot - Constants::g * cos(DEG_TO_RAD(slope_angle));
		
		// Velocities
		xdot += xdotdot * DT;
		zdot += zdotdot * DT;
		
		// Positions
		x += xdot * DT;
		z += zdot * DT;
		
		ubar += xdot;
		if (z > z_max) {
			z_max = z;
		}
		nt++;
	} while (z > z0);
	
	if (!((nt > 1) && (xdot > 0.0)) ) {
		return false;
	} else {
		// Mean horizontal velocity
		ubar /= nt;
		// Impact velocity
		u_i = sqrt(pow(xdot, 2) + pow(zdot, 2));
		// Impact angle
		angle_i_rad = atan(-zdot / xdot);
		// Flight time
		t_i = nt * DT;
		
		return true;
	}
}

/**
 * @brief Computes the saltation flux for one bottom element
 * @param z0e
 * @param tauS
 * @param tauA
 * @param slope
 * @param dg
 * @param tau_th
 * @param z_max
 * @param ubar
 * @param cs
 * @return mass flux
 */
double Saltation::sa_MassFlux(const double& z0, const double& tauS, const double& tauA, const double& slope_angle,
                              const double& dg, const double& tau_th, double& z_max, double& ubar, double& cs)
{
  	//  double tau_th=0.04;
	// pa preliminary until coupled to SNOWPACK
	double u0, angle_e_rad;
	double u_i, angle_i_rad, t_i;
	double Nreb, tau1G_reb;
	double mass;
	double hsalt, udisturb, ulog;
	
	// Initialize mass: Model Mass Flus is extremely sensitive to grain mass
	mass = Constants::density_ice * 4. / 3. * Constants::pi * pow( (dg / 2.), 3. );
	// Compute trajectories until stationary
	//  angle_e_rad = Constants::pi/6.;
	angle_e_rad = DEG_TO_RAD(Saltation::angle_ej);
	//  u0 = 0.63/sin(angle_e_rad)*sqrt(tauS/DENSITY_AIR);
	u0 = Saltation::ratio_ve_ustar * sqrt((tauS - tau_th) / Constants::density_air);
	//  u0 = 3.7*sqrt(tauS-tau_th);
	
	u_i = u0;
	// Iterate until stationary trajectory
	while (u_i == u0) {
		if (!sa_Traject(u0, angle_e_rad, slope_angle, dg, tauA, tauS, z0, ubar, u_i, angle_i_rad, t_i, z_max)) {
			cs = 0.;
			return 0.;
		}
	}
	
	// Shear stress from one grain
	tau1G_reb = mass * ( u_i * cos(angle_i_rad) - u0 * cos(angle_e_rad) ) / t_i;
	//  tau1G_reb = mass*( u_i*cos(angle_i_rad) - u0*cos(angle_e_rad) );
	
	// Number of rebounding grains
	Nreb = (tauS - tauA) / tau1G_reb;
	
	// Now compute the concentration
	if ( (u0 > 0.0) && (ubar > 0.0) ) {
		hsalt = z0 + Saltation::hs_frac * 0.5 * pow((u0 * cos(angle_e_rad)), 2) / 2. / Constants::g;
		udisturb = sa_vw(hsalt, tauA, tauS, z0, u0 * sin(angle_e_rad), slope_angle);
		ulog = sqrt(tauS / Constants::density_air) / Saltation::karman * log(hsalt / z0);
		cs = MIN (0.02, Constants::density_air * pow((ulog - udisturb) / (ubar), 2));
	} else {
		cs = 0.0;
	}
	
	if ( !(mass*Nreb*ubar > 0. && mass*Nreb*ubar < 1e10) ) {
		prn_msg(__FILE__, __LINE__, "msg+", Date(),
			"Infinite Flux from Mass Flux  mass:%lf NrebG:%lf ubar:%lf", mass, Nreb, ubar);
	}
	
	return (mass*Nreb*ubar);
} // End of MassFlux

/**
 * @brief Computes the saltation flux for one bottom element
 * @param z0
 * @param tauS
 * @param slope_angle (deg)
 * @param dg
 * @param tau_th
 * @param flux
 * @param z_max
 * @param ubar
 * @param cs
 * @return surface shear stress
*/
double Saltation::sa_AeroEntrain(const double& z0, const double& tauS, const double& slope_angle, const double& dg, 
                                 const double& tau_th, double& flux, double& z_max, double& ubar, double& cs)
{
	int    iter = 0, maxit = 40;
	double tauA, u0, angle_e_rad;
	double u_i, angle_i_rad, t_i;
	double eps = 0.001, tauA_old;
	double Nae, n_ae;
	double mass;
	double hsalt, udisturb, ulog;
	
	// Initialize mass, entrainment number and surface shear stress
	mass = Constants::density_ice * 4. / 3. * Constants::pi * pow((dg / 2.), 3.);
	angle_e_rad = DEG_TO_RAD(Saltation::angle_ej);
	//  n_ae = 1./(1.09*mass*sqrt(tauS/DENSITY_AIR));
	n_ae  = 0.5 / 8. / Constants::pi / dg / dg; 
	//  n_ae = 200000;
	tauA = (tauS + tau_th) / 2;
	
	// Compute trajectories until stationary
	//  u0 = 0.63/sin(angle_e_rad)*sqrt(tauS/DENSITY_AIR);
	u0 = Saltation::ratio_ve_ustar * sqrt((tauS - tau_th) / Constants::density_air);
	//  u0 = 3.7*sqrt(tauS-tau_th);
	
	iter=0;
	do {
		if (!sa_Traject(u0, angle_e_rad, slope_angle, dg, tauA, tauS, z0, ubar, u_i, angle_i_rad, t_i, z_max)) {
			flux = 0.;
			cs = 0.;
			return tauS;
		}
		tauA_old = tauA;
		Nae = n_ae * (tauA - tau_th);
		tauA = tauS - Nae *mass * (u_i * cos(angle_i_rad) - u0 * cos(angle_e_rad)) / t_i;
		//    tauA = tauS - Nae*mass*(u_i*cos(angle_i_rad) - u0*cos(angle_e_rad));
		if(tauA < tau_th) tauA = tau_th; 
		tauA = (tauA+tauA_old)/2.; 
		//     fprintf(stdout, "n_ae=%lf, tauA=%lf \n",n_ae, tauA);
		iter++;
	} while ( (fabs((tauA - tauA_old)/tauA) > eps) && (iter < maxit) );
	
	if (iter == maxit) {
		prn_msg(__FILE__, __LINE__, "wrn", Date(), "Airborne shear stress Iteration did not converge");
		tauA = 0.5 * tauS;
		Nae = n_ae * (tauA - tau_th);
	}
	if (tauA < tau_th) {
		prn_msg(__FILE__, __LINE__, "wrn", Date(), "Airborne shear stress smaller than threshold");
		Nae = 0.; 
	}
	
	flux = mass * Nae * ubar;
	
	if ( !((flux > 0.) && (flux < 1e10)) ) {
		prn_msg(__FILE__, __LINE__, "msg+", Date(),
			"Infinite Flux from Aero Entrain  mass:%lf Nae:%lf ubar:%lf",mass,Nae,ubar);
	}
	
	// Now compute the concentration
	if ( (u0 > 0.0) && (ubar > 0.0) ){
		hsalt = z0 + Saltation::hs_frac * 0.5 * pow((u0 * cos(angle_e_rad)), 2) / 2. / Constants::g;
		udisturb = sa_vw(hsalt, tauA, tauS, z0, u0 * sin(angle_e_rad), slope_angle);
		ulog = sqrt(tauS / Constants::density_air) / Saltation::karman * log(hsalt / z0);
		cs = MIN (0.02, Constants::density_air * pow((ulog - udisturb) / ubar, 2));
	} else {
		cs = 0.0;
	}
	
	return(tauA);

} // End of sa_AeroEntrain

/**
 * @brief Computes the saltation flux for one bottom element
 * @param z0
 * @param tauS
 * @param tauA
 * @param slope_angle (deg)
 * @param dg
 * @param tau_th
 * @param z_max
 * @param ubar
 * @return saltation weak or strong
 */
int Saltation::sa_TestSaltation(const double& z0, const double& tauS, const double& tauA, const double& slope_angle,
                                const double& dg, const double& tau_th, double& z_max, double& ubar)
{
	int j;
	double u0, angle_e_rad;
	double u_i, angle_i_rad, t_i;
	double hsalt1, hsalt2, mass;
	
	// Initialize mass
	mass = Constants::density_ice * 4. / 3. * Constants::pi * pow((dg / 2.), 3);
	// Compute trajectories until stationary
	angle_e_rad = DEG_TO_RAD(Saltation::angle_ej);
	//  u0 = 0.63/sin(angle_e_rad)*sqrt(tauS/DENSITY_AIR);
	u0 = Saltation::ratio_ve_ustar * sqrt((tauS - tau_th) / Constants::density_air);
	
	// Compute the first trajectory
	sa_Traject(u0, angle_e_rad, slope_angle, dg, tauA, tauS, z0, ubar, u_i, angle_i_rad, t_i, z_max);
	hsalt1 = z_max;
	u0 = Saltation::elas * u_i;
	// u0 = sqrt(ELAS*DSQR(u_i)) ;
	
	// Compute three trajectories to see whether they are growing
	for (j = 0; j < 3; j++) {
		sa_Traject(u0, angle_e_rad, slope_angle, dg, tauA, tauS, z0, ubar, u_i, angle_i_rad, t_i, z_max);
		hsalt2 = z_max;
		u0 = Saltation::elas * u_i;
		//    u0 = sqrt(ELAS*DSQR(u_i));
	}
		
	if (hsalt2 > hsalt1) {
		return (Saltation::strong);
	} else {
		return (Saltation::weak);
	}

} // End of sa_TestSaltation

/**
 * @brief Computes the saltation flux for one bottom element
 * @param i_tauS
 * @param tau_th
 * @param slope_angle
 * @param dg
 * @param massflux
 * @param c_salt
 * @return bool
*/
bool Saltation::compSaltation(const double& i_tauS, const double& tau_th, const double& slope_angle, const double& dg,
                                  double& massflux, double& c_salt)
{
	int    j, k = 5;
	double tauS = i_tauS;
	double eps = 1e-5, ustar, ustar_thresh;
	double tauA, tauA_left, tauA_right, tauA_middle, tau_r, tau_e = 0.;
	double taumean, taumax, taustep, Cp, tau_j, Ptau_j, ubar = 0., z_lower = 0.;
	double flux_mean = 0., cs_mean = 0., flux, cs;  // What we finally want
	
	if (!doorschot) {
		// Sorensen
		ustar = sqrt(tauS / Constants::density_air);
		ustar_thresh = sqrt(tau_th / Constants::density_air);
		if (ustar > ustar_thresh) {
			massflux = 0.0014 * Constants::density_air * ustar * (ustar - ustar_thresh) * (ustar + 7.6*ustar_thresh + 205);
			c_salt = massflux / ustar*0.001;
			// Arbitrary Scaling to match Doorschot concentration
		} else {
			massflux = 0.;
			c_salt = 0.;
		}
	} else {
		// Judith
		// Initialize Shear Stress Distribution
		taumean = tauS;
		taumax = 15.* tauS;
		taustep = (taumax - tau_th) / k;
		Cp = 1. / taumean;
	
	
		for (j = 0; j < k; j++) {
			tau_j = tau_th + ((double)(j) + 0.5) * taustep;
			Ptau_j = exp(-Cp * (tau_j - 0.5 * taustep)) - exp(-Cp * (tau_j+ 0.5 * taustep));
	
			tauS = tau_j;
	
			if(tauS > tau_th) {
	
				// First test for large rebound thresholds
				if (sa_TestSaltation(Saltation::z0_salt, tauS, tauS, slope_angle, dg, tau_th, z_lower, ubar) == Saltation::weak) {
					tauA = sa_AeroEntrain(Saltation::z0_salt, tauS, slope_angle, dg, tau_th, flux, z_lower, ubar, cs);
				} else {
					// Use an iterative method to determine the rebound threshold at the ground
					tauA_right = tauS;
					tauA_left = 0.0;
					do {
						tauA_middle = (tauA_left + tauA_right) / 2.;
						// fprintf(stdout, "tauA=%lf \n",tauA_middle);
						if (sa_TestSaltation(Saltation::z0_salt, tauS, tauA_middle, slope_angle,
									      dg, tau_th, z_lower, ubar) == Saltation::strong) {
							tauA_right = tauA_middle;
						} else {
							tauA_left = tauA_middle;
						}
					} while (tauA_right - tauA_left > eps);
					tau_r = tauA_middle;
	
					/*
					* Distinguish the different possibilities after Judith and compute
					* the flux; Start with computation of tau_e: The surface shear stress
					* for the hypothetical case of aerodynamic entrainment only given a
					* certain overall shear stress, tauS
					*/
					tau_e = sa_AeroEntrain(Saltation::z0_salt, tauS, slope_angle, dg, tau_th, flux, z_lower, ubar, cs);
					if (tau_e < tau_th) {
						tauA = tau_r;
						flux = sa_MassFlux(Saltation::z0_salt, tauS, tauA, slope_angle, dg, tau_th, z_lower, ubar, cs);
					} else if (tau_e < tau_r) {
						tauA = tau_e;
					} else {
						tauA = tau_r; // Flux computation is wrong, must be redone
						flux = sa_MassFlux(Saltation::z0_salt, tauS, tauA, slope_angle, dg, tau_th, z_lower, ubar, cs);
					}
				} // else large rebound threshold
				cs_mean += cs * Ptau_j;
				flux_mean += flux * Ptau_j;
	
			} // if there is s.th. to do
	
		} // for all shear stress classes
	
		// Fill return Values
		massflux = flux_mean;
		c_salt = cs_mean;
	}
	
	//  printf("Saltation Conc. %lf", c_salt);
	
	return true;
}

