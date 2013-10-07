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

#include <snowpack/SnowDrift.h>
#include <assert.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/

///Deviation from geometrical factors defined by Schmidt
const double SnowDrift::schmidt_drift_fudge = 1.0;

///Enables erosion notification
const bool SnowDrift::msg_erosion = false;


/************************************************************
 * non-static section                                       *
 ************************************************************/

SnowDrift::SnowDrift(const SnowpackConfig& cfg) : saltation(cfg),
                     enforce_measured_snow_heights(false), snow_redistribution(false), alpine3d(false),
                     sn_dt(0.), nSlopes(0)
{
	cfg.getValue("ALPINE3D", "SnowpackAdvanced", alpine3d);

	/**
	 * @brief Defines how the height of snow is going to be handled
	 * - false: Depth of snowfall is determined from the water equivalent of snowfall (HNW)
	 * - true : The measured height of snow is used to determine whether new snow has been deposited.
	 *      This setting MUST be chosen in operational mode. \n
	 *      This procedure has the disadvantage that if the snowpack settles too strongly
	 *      extra mass is added to the snowpack. \n
	 * New snow density is needed in both cases, either parameterized, measured, or fixed.
	 */
	cfg.getValue("ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack", enforce_measured_snow_heights);

	/*
	 * Number of aspects incl. the real flat field: at least 1, either 5 or 9 for SNOW_REDISTRIBUTION
	 * - 1 real simulation on flat field (or one slope, see PERP_TO_SLOPE)
	 * - 5 real simulation on flat field plus 4 virtual slopes
	 * - 9 real simulation on flat field plus 8 virtual slopes
	 */
	cfg.getValue("NUMBER_SLOPES", "Snowpack", nSlopes);

	// Defines whether real snow erosion and redistribution should happen under
	// blowing snow conditions. Set in operational mode.
	if(nSlopes>1)
		cfg.getValue("SNOW_REDISTRIBUTION", "Snowpack", snow_redistribution);

	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	const double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	sn_dt = M_TO_S(calculation_step_length);

}

/**
 * @brief Computes the local mass flux of snow
 * @bug Contribution from suspension not considered yet!
 * @param *Edata
 * @param ustar Shear wind velocity (m s-1)
 * @param angle Slope angle (rad)
 * @return Saltation mass flux (kg m-1 s-1)
 */
double SnowDrift::compMassFlux(const ElementData& Edata, const double& ustar, const double& slope_angle)
{
	// Compute basic quantities that are needed: friction velocity, z0, threshold vw
	// For now assume logarithmic wind profile; TODO change this later
	const double weight = 0.02 * Constants::density_ice * (Edata.sp + 1.) * Constants::g * MM_TO_M(Edata.rg);
	// weight = Edata.Rho*(Edata.sp + 1.)*Constants::g*MM_TO_M(Edata.rg);
	const double sig = 300.;
	const double binding = 0.0015 * sig * Edata.N3 * Optim::pow2(Edata.rb/Edata.rg);
	const double tau_thresh = SnowDrift::schmidt_drift_fudge * (weight + binding);  // Original value for fudge: 1. (Schmidt)
	//const double ustar_thresh = sqrt(tau_thresh / Constants::density_air);
	// fprintf(stdout, "weight: %lf   binding:%lf  tau_th:%lf\n",weight, binding, tau_thresh);
	const double tau = Constants::density_air * Optim::pow2(ustar);

	// First, look whether there is any transport at all: use formulation of Schmidt
	if ( tau_thresh > tau ) {
		return (0.0);
	}

	// Compute the saltation mass flux (after Pomeroy and Gray)
	double Qsalt = 0., Qsusp = 0., c_salt; // The mass fluxes in saltation and suspension (kg m-1 s-1)
	if (!saltation.compSaltation(tau, tau_thresh, slope_angle, MM_TO_M(2.*Edata.rg), Qsalt, c_salt)) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Saltation computation failed");
		throw IOException("Saltation computation failed", AT);
	}

	if (Qsalt > 0.) {
		Qsusp = 0.; // TODO What about sd_IntegrateFlux(ustar, ustar_thresh, Qsalt, z0); ???
	} else {
		Qsalt = 0.;
		Qsusp = 0.;
	}

	return (Qsalt + Qsusp);
}

/**
 * @brief Erodes Elements from the top and computes the associated mass flux
 * @brief Even so the code is quite obscure, it should cover all of the following cases:
 * -# externally provided eroded mass given with forced_massErode
 * -# Flat field simulation, where two flux values are computed from vw and vw_drift,
 *    using either the true snow surface or the ErosionLevel marker (virtual erosion layer), respectively.
 * -# Windward virtual slope simulations with vw_drift and the true snow surface
 * -# Slope simulations in research mode, using vw and the true snow surface
 * @param Mdata
 * @param Xdata
 * @param Sdata
 * @param forced_massErode if greater than 0, force the eroded mass to the given value (instead of computing it)
*/
void SnowDrift::compSnowDrift(const CurrentMeteo& Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata, double& forced_massErode)
{
	int nErode = 0;                           // number of eroded elements and erosion level
	double real_flux = 0., virtual_flux = 0.; // mass flux, either real or virtual

	size_t nE = Xdata.getNumberOfElements();
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;

	const bool no_snow = ((nE < Xdata.SoilNode+1) || (EMS[nE-1].theta[SOIL] > 0.));
	const bool no_wind_data = (Mdata.vw_drift == mio::IOUtils::nodata);
	if (no_snow || no_wind_data) {
		Xdata.ErosionLevel = Xdata.SoilNode;
		Xdata.ErosionMass = 0.;
		if (no_snow) {
			Sdata.drift = 0.;
		} else {
			Sdata.drift = Constants::undefined;
		}
		return;
	}

	if (Xdata.ErosionLevel > nE-1) {
		prn_msg(__FILE__, __LINE__, "wrn", Mdata.date, "ErosionLevel=%d did get messed up (nE-1=%d)", Xdata.ErosionLevel, nE-1);
	}

	bool windward = false; // Set to true for windward slope
	if (!alpine3d && (Xdata.meta.getSlopeAngle() > Constants::min_slope_angle) && Xdata.windward) {
		windward = true;
	}
	// Scale ustar_max
	const double ustar_max = (Mdata.vw>0.1) ? Mdata.ustar * Mdata.vw_drift / Mdata.vw : 0.;

	// Evaluate possible real erosion
	double massErode = 0.; // eroded mass loss due to erosion
	if (alpine3d) {
		massErode = MAX(0., -forced_massErode); //negative mass is erosion
	} else {
		try {
			if (enforce_measured_snow_heights && !windward) {
				real_flux = compMassFlux(EMS[nE-1], Mdata.ustar, Xdata.meta.getSlopeAngle()); // Flat Field, vw && nE-1
			} else {
				real_flux = compMassFlux(EMS[nE-1], ustar_max, Xdata.meta.getSlopeAngle()); // Windward slope && vw_drift && nE-1
			}
		} catch(const exception&){
			prn_msg(__FILE__, __LINE__, "err", Mdata.date, "SnowDrift");
			throw;
		}
		massErode = (real_flux * sn_dt / Hazard::typical_slope_length);
	}

	// Real erosion either on flat field or on windward slope or in Alpine3D
	if ((snow_redistribution && (((Xdata.mH + 0.02) < (Xdata.cH - Xdata.Ground)) && (Xdata.meta.getSlopeAngle() <= Constants::min_slope_angle))) || alpine3d || windward) {
		Xdata.ErosionMass = 0.;
		// Erode at most one element with a maximal error of +- 5 % on mass ...
		if (massErode >= 0.95 * EMS[nE-1].M) {
			if (windward) {
				Xdata.rho_hn = EMS[nE-1].Rho;
			}
			nE--;
			Xdata.cH -= EMS[nE].L;
			NDS[nE].hoar = 0.;
			massErode -= EMS[nE].M;
			Xdata.ErosionMass = EMS[nE].M;
			Xdata.ErosionLevel = MIN(nE-1, Xdata.ErosionLevel);
			nErode++;
			forced_massErode = -massErode;
		} else if (massErode > 0.) { // ... or take away massErode from top element - partial real erosion
			if (fabs(EMS[nE-1].L * EMS[nE-1].Rho - EMS[nE-1].M) > 0.001) {
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date, "Inconsistent Mass:%lf   L*Rho:%lf", EMS[nE-1].M,EMS[nE-1].L*EMS[nE-1].Rho);
				EMS[nE-1].M = EMS[nE-1].L * EMS[nE-1].Rho;
				assert(EMS[nE-1].M>=0.); //mass must be positive
			}
			if (windward) {
				Xdata.rho_hn = EMS[nE-1].Rho;
			}
			const double dL = -massErode / (EMS[nE-1].Rho);
			NDS[nE].z += dL;
			EMS[nE-1].L0 = EMS[nE-1].L = EMS[nE-1].L + dL;
			Xdata.cH += dL;
			NDS[nE].z += NDS[nE].u;
			NDS[nE].u = 0.0;
			NDS[nE].hoar = 0.;
			EMS[nE-1].M -= massErode;
			assert(EMS[nE-1].M>=0.); //mass must be positive
			Xdata.ErosionMass = massErode;
			nErode = -1;
			forced_massErode = 0.;
		}
		/* ... or check whether you can do a virtual erosion for drift index generation
		 * in case of no (or small) real erosion for all cases except Alpine3D or virtual slopes
		 */
	} else if (((nSlopes == 1) || (!snow_redistribution && (Xdata.meta.getSlopeAngle() < Constants::min_slope_angle))) && (Xdata.ErosionLevel > Xdata.SoilNode)) {
		virtual_flux = compMassFlux(EMS[Xdata.ErosionLevel], ustar_max, Xdata.meta.getSlopeAngle()*mio::Cst::to_rad);
		// Convert to eroded snow mass
		if ((massErode = virtual_flux*sn_dt / Hazard::typical_slope_length) > 0.) {
			nErode = -1;
			Sdata.mass[SurfaceFluxes::MS_WIND] = massErode;
		}
		// Add (negative) value stored in Xdata.ErosionMass
		if ( Xdata.ErosionMass < 0. ) {
			massErode -= Xdata.ErosionMass;
		}
		// Now keep track of mass that either did or did not lead to erosion of full layer
		if (massErode > EMS[Xdata.ErosionLevel].M) {
			massErode -= EMS[Xdata.ErosionLevel].M;
			Xdata.ErosionLevel--;
		}
		Xdata.ErosionMass = (-massErode);
		Xdata.ErosionLevel = MAX(Xdata.SoilNode, MIN(Xdata.ErosionLevel, nE-1));
	}

	// If real or virtual erosion took place, take the corresponding flux the other being zero
	if (nErode != 0 && (massErode > 0.)) {
		Sdata.drift = MAX(real_flux, virtual_flux);
		nErode = MAX(0, nErode);
	}

	if ( nErode > 0 ) {
		if (SnowDrift::msg_erosion && (Xdata.ErosionMass > 0.) && !alpine3d) { //reduce number of warnings for Alpine3D
			if ( windward ) {
				prn_msg(__FILE__, __LINE__, "msg+", Mdata.date, "Eroding %d layer(s) w/ total mass %.3lf kg/m2 (windward: azi=%.1lf, slope=%.1lf)", nErode, Xdata.ErosionMass, Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
			} else {
				prn_msg(__FILE__, __LINE__, "msg+", Mdata.date, "Eroding %d layer(s) w/ total mass %.3lf kg/m2 (azi=%.1lf, slope=%.1lf)", nErode, Xdata.ErosionMass, Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
			}
		}
		Xdata.resize(nE);
	}
}

/*
 * End of SnowDrift.cc
 */
