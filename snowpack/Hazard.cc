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
 * @file Hazard.c
 * @version 10.02
 * This module contains the hazard calculation routines
*/

#include <snowpack/Hazard.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/

//Lee slope length (m) used to convert mass flux to drift index deposition depth (cm/24h)
const double Hazard::typical_slope_length = 70.0;

//Predefined snow density (kg m-3) used to convert mass flux to drift index deposition depth (cm/24h)
const double Hazard::wind_slab_density = 77.0;

//At least that mass flux (kg m-1 h-1) must have summed up for the drift index to be larger than 0.
const double Hazard::minimum_drift = 0.0;

//Largest additional accumulation allowed (cm h-1)
const double Hazard::maximum_drift = 5.0;

//Selects whether hand hardness index R is output either in N (Swiss scale) or steps
const bool Hazard::r_in_n = true;

/************************************************************
 * non-static section                                       *
 ************************************************************/

Hazard::Hazard(const mio::Config& i_cfg) : cfg(i_cfg) 
{
	/**
	 * @brief Defines how the height of snow is going to be handled
	 * - 0: Depth of snowfall is determined from the water equivalent of snowfall (HNW)
	 * - 1: The measured height of snow is used to determine whether new snow has been deposited.
	 *      This setting MUST be chosen in operational mode. \n
	 *      This procedure has the disadvantage that if the snowpack settles too strongly
	 *      extra mass is added to the snowpack. \n
	 * New snow density is needed in both cases, either parameterized, measured, or fixed.
	 */
	enforce_measured_snow_heights = cfg.get("ENFORCE_MEASURED_SNOW_HEIGHTS", "Parameters");

	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Parameters");
	sn_dt = M_TO_S(calculation_step_length);

	/*
	 * Hazard data interval in units of CALCULATION_STEP_LENGTH \n
	 * WARNING: In operational mode, this has to result in a 30 min interval!
	 * It is a matter of consitency. If you change this, a big mess will result!!!
	 */
	hazard_steps_between = cfg.get("HAZARD_STEPS_BETWEEN", "Parameters");

	/* Dew point relative to water or ice
	 * - default: 1
	 * - Antarctica: 0 */
	force_rh_water = cfg.get("FORCE_RH_WATER", "Parameters");

	research_mode = cfg.get("RESEARCH", "Parameters");

	//Density of surface hoar (-> hoar index of surface node) (kg m-3)
	density_hoar_surf = cfg.get("DENSITY_HOAR_SURF", "Parameters");

	//Minimum size to show surface hoar on surface (mm)
	min_size_hoar_surf = cfg.get("MIN_SIZE_HOAR_SURF", "Parameters");
}

/**
 * @brief Allocates and initializes Hazard data (flat field station only)
 * - Fills the snowpack version number, date of computation, user, ...
 * - Computes a zeroth order drift index for the first time step w/o shifting OldDrift!
 * @author Michael Lehning
 * @version 10.03
 * @param TimeEnd
 * @param *OldDrift
 * @param SlopeAngle
 * @return Number of hazard steps
 */
void Hazard::initializeHazard(const double TimeEnd, double *OldDrift, double SlopeAngle, 
						std::vector<Q_PROCESS_DAT>& Hdata, std::vector<Q_PROCESS_IND>& Hdata_ind)
{
	int nHz = (int)(TimeEnd / (hazard_steps_between * M_TO_S(calculation_step_length))) + 2;
	if (nHz <= 0) nHz = 1;

	Hdata = vector<Q_PROCESS_DAT>(nHz);
	Hdata_ind = vector<Q_PROCESS_IND>(nHz);

	memset(&Hdata[0], 0, sizeof(Q_PROCESS_DAT)*nHz);
	memset(&Hdata_ind[0], 0, sizeof(Q_PROCESS_IND)*nHz);

	Hdata[0].nHz = nHz;
	Hdata[nHz-1].nHz = nHz;

	qr_VersionUserRuntime(Hdata[0].sn_version, Hdata[0].sn_computation_date, 
					  &Hdata[0].sn_jul_computation_date, Hdata[0].sn_user, Hdata[0].sn_compile_date);

	Hdata[0].wind_trans = driftIndex(OldDrift, 0., Hazard::wind_slab_density, 6, -1) / cos(SlopeAngle);
	Hdata[0].wind_trans24 = driftIndex(OldDrift, 0., Hazard::wind_slab_density, 24, -1) / cos(SlopeAngle);
}

/**
 * @brief Determines a nhour drift index
 * - shift
 * 	-  1 : shift OldDrift, overwrite OldDrift[0]
 * 	-  0 : overwrite OldDrift[0] only
 * 	- -1 : no action at all
 * @author Michael Lehning
 * @version Y.mm
 * @param *OldDrift double
 * @param Drift double
 * @param rho double
 * @param nhour double
 * @param shift double
 */
double Hazard::driftIndex(double *OldDrift, double Drift, double rho, int nhour, int shift)
{
	int    i;
	double sumindex = 0., flux, ero_depo;

	switch ( shift ) {
		case 1: // Shift drift data
			for(i=47; i>0; i--) {
				OldDrift[i] = OldDrift[i-1];
			}
		case 0: // Overwrite OldDrift[0]
			OldDrift[0] = Drift;
		default:
			break;
	}

	for (i = 0; i < 2*nhour; i++ ) {
		sumindex += OldDrift[i];
	}

	flux = H_TO_S ( MAX (0,(sumindex - Hazard::minimum_drift)) / (2. * nhour) ); // kg m-1 h-1
	ero_depo = M_TO_CM ( flux * nhour / (Hazard::typical_slope_length * rho) );

	if ( ero_depo > (nhour * Hazard::maximum_drift) ) {
		ero_depo = nhour * Hazard::maximum_drift;
	}
	return(ero_depo);
}

/**
 * @brief Determines hoar mass index for last n hours in kg m-2
 * @author Michael Lehning
 * @version Y.mm
 * @param *OldHoar double
 * @param new_hoar double
 * @param nhour double
 * @param new_step double
 * @return double
 */
double Hazard::hoarIndex(double *OldHoar, double new_hoar, int nhour, int new_step)
{
	int i;
	double hoar_ind = 0.;

	// Shift hoar data
	if ( new_step )
	for (i = 47; i > 0; i-- ) {
		OldHoar[i] = OldHoar[i-1];
	}
	OldHoar[0] = new_hoar;

	// Determine hoar_ind
	for (i = 0; i < 2*nhour; i++ ) {
		hoar_ind += OldHoar[i];
	}
	return(hoar_ind);
}

void Hazard::calculateMeltFreezeCrust(const SN_STATION_DATA& Xdata, Q_PROCESS_DAT& Hdata, Q_PROCESS_IND& Hdata_ind)
{
	double crust_dep=0., crust_thick=0.;
	double cos_sl = cos(Xdata.SlopeAngle);
	int e = Xdata.getNumberOfElements()-1;

	while( (e > Xdata.SoilNode) && (crust_dep <= 0.03) ) {
		if ( (Xdata.Edata[e].type == 772) || (Xdata.Edata[e].type == 880) ) {
			crust_thick += Xdata.Edata[e].L/cos_sl;
			if ( (Xdata.Edata[e-1].type != 772) && (Xdata.Edata[e-1].type != 880) ) {
				break;
			}
		} else {
			crust_dep += Xdata.Edata[e].L/cos_sl;
		}
		e--;
	}
	if ( (crust_thick >= 0.) && (crust_thick <= Xdata.cH/cos_sl) ) {
		Hdata.crust = M_TO_CM(crust_thick);
	} else {
		Hdata_ind.crust = -1;
	} 
}

/***********************************************************************************************/
double Hazard::calcDewPointDeficit(double TA, double TSS, double RH)
/*---------------------------------------------------------------------------------------------+
 | Determines the dew point deficit in degC                                                    |
 +---------------------------------------------------------------------------------------------*/
{ 
  double b=9.5, c=265.5;
  
  TA = K_TO_C(TA);
  TSS = K_TO_C(TSS);

  return(TSS - (c*(log10(RH) + b*TA/(c+TA))/
                (b - log10(RH) - b*TA/(c+TA))));
  
} // End of ml_hz_DewPointDeficit


/**
 * @brief Given the "Zwischen" data containing the depth-hoar index, the wind-drift index and the
 * three and twenty-four newsnowfall rates, this routine calculates the Hdata.
 * @param d_hs6
 * @param d_hs24
 * @param Zdata
 * @param Hdata
 * @param Hdata_ind
 * @param Xdata
 * @param Mdata
 * @param Sdata
 * @param nAvg
*/
void Hazard::calculateHazard(const double& d_hs6, const double& d_hs24, const SN_STATION_DATA& Xdata,
					    const SN_MET_DATA& Mdata, const int& nAvg, SN_ZWISCHEN_DATA& Zdata, 
					    Q_PROCESS_DAT& Hdata, Q_PROCESS_IND& Hdata_ind, SN_SURFACE_DATA&  Sdata)
{
	int    e, k, l, nE;
	double t_hns[9], sumhns = 0., sum_we = 0., hns[9], swe[9];
	double cos_sl, hs, h_top1, h_top2;

	const SN_ELEM_DATA *EMS;  // Pointer to element data

	// Dereference Pointers
	nE  = Xdata.getNumberOfElements();
	EMS = &Xdata.Edata[0];

	// Initialization plus dummy values for "RESEARCH VERSION"
	// cos of slope angle to convert Hdata to vertical
	cos_sl = cos( Xdata.SlopeAngle );
	hs = Xdata.cH / cos_sl;

	Hdata.date = Mdata.date;

	Hdata.dewpt_def = 21.7;    Hdata_ind.dewpt_def  = 0;
	Hdata.hoar_ind6 = 21.7;    Hdata_ind.hoar_ind6 = 0;
	Hdata.hoar_ind24 = 21.7;   Hdata_ind.hoar_ind24 = 0;
	Hdata.hoar_size = 21.7;    Hdata_ind.hoar_ind24 = 0;

	Hdata.wind_trans = 21.7;   Hdata_ind.wind_trans = 0;
	Hdata.wind_trans24 = 21.7; Hdata_ind.wind_trans24 = 0;

	Hdata.hns3 = 21.7;         Hdata_ind.hns3 = 0;
	Hdata.hns6 = 21.7;         Hdata_ind.hns6 = 0;
	Hdata.hns12 = 21.7;        Hdata_ind.hns12 = 0;
	Hdata.hns24 = 21.7;        Hdata_ind.hns24 = 0;
	Hdata.hns72 = 21.7;        Hdata_ind.hns72 = 0;
	Hdata.hns72_24 = 21.7;     Hdata_ind.hns72_24 = 0;
	Hdata.wc3  = 21.7;         Hdata_ind.wc3  = 0;
	Hdata.wc6  = 21.7;         Hdata_ind.wc6  = 0;
	Hdata.wc12 = 21.7;         Hdata_ind.wc12 = 0;
	Hdata.wc24 = 21.7;         Hdata_ind.wc24 = 0;
	Hdata.wc72 = 21.7;         Hdata_ind.wc72 = 0;

	Hdata.stab_class1 = 0;     Hdata_ind.stab_class1 = 0;
	Hdata.stab_class2 = 5;     Hdata_ind.stab_class2 = 0;

	Hdata.stab_index1 = 6.;    Hdata_ind.stab_index1 = 0;
	Hdata.stab_index2 = 6.;    Hdata_ind.stab_index2 = 0;
	Hdata.stab_index3 = 6.;    Hdata_ind.stab_index3 = 0;
	Hdata.stab_index4 = 6.;    Hdata_ind.stab_index4 = 0;
	Hdata.stab_index5 = 6.;    Hdata_ind.stab_index5 = 0;
	Hdata.stab_height1 = hs;   Hdata_ind.stab_height1 = 0;
	Hdata.stab_height2 = hs;   Hdata_ind.stab_height2 = 0;
	Hdata.stab_height3 = hs;   Hdata_ind.stab_height3 = 0;
	Hdata.stab_height4 = hs;   Hdata_ind.stab_height4 = 0;
	Hdata.stab_height5 = hs;   Hdata_ind.stab_height5 = 0;

	Hdata.ch = M_TO_CM(hs);    Hdata_ind.ch = 0;

	Hdata.swe     = 0.;        Hdata_ind.swe     = 0;
	Hdata.tot_lwc = 0.;        Hdata_ind.tot_lwc = 0;
	Hdata.runoff  = 0.;        Hdata_ind.runoff  = 0;

	Hdata.crust  =  0.0;       Hdata_ind.crust  = 0;
	Hdata.en_bal = 21.7;       Hdata_ind.en_bal = 0;
	Hdata.sw_net = 21.7;       Hdata_ind.sw_net = 0;
	Hdata.t_top1 = 21.7;       Hdata_ind.t_top1 = 0;
	Hdata.t_top2 = 21.7;       Hdata_ind.t_top2 = 0;

	// Calculate new snow depths for given time intervals
	t_hns[0] = 0.5;
	t_hns[1] = 3.;     // Caution! Originally 3 was assigned to t_hns[0].
	t_hns[2] = 6.;
	t_hns[3] = 12.;
	t_hns[4] = 24.;
	t_hns[5] = 72.;
	for (k = 0, e = nE-1; k <= 5; k++) {
		if ( e < 0 ) {
			hns[k] = 0.0;
			swe[k] = 0.0;
		} else {
			while ( (e >= Xdata.SoilNode) && ((Mdata.date.getJulianDate() - EMS[e].date.getJulianDate()) < (H_TO_D(t_hns[k]))) ) {
				sumhns += EMS[e].L;
				sum_we += EMS[e].L * EMS[e].Rho;
				e--;
			}
			hns[k] = sumhns;
			swe[k] = sum_we;
		}
	}
	Hdata.hns_half_hour = M_TO_CM(hns[0] / cos_sl);
	Hdata.hns3 =  M_TO_CM(hns[1] / cos_sl);
	Hdata.hns6 =  M_TO_CM(hns[2] / cos_sl);
	Hdata.hns12 =  M_TO_CM(hns[3] / cos_sl);
	Hdata.hns24 =  M_TO_CM(hns[4] / cos_sl);
	Hdata.hns72 =  M_TO_CM(hns[5] / cos_sl);
	Hdata.wc_half_hour = swe[0] / cos_sl;
	Hdata.wc3 =  swe[1] / cos_sl;
	Hdata.wc6 =  swe[2] / cos_sl;
	Hdata.wc12 =  swe[3] / cos_sl;
	Hdata.wc24 =  swe[4] / cos_sl;
	Hdata.wc72 =  swe[5] / cos_sl;

	// Calculate 72h sum of 24h new snow depths for a total of 3 days
	for (l = 143; l > 0; l--) {
		Zdata.hns24[l] = Zdata.hns24[l-1];
	}
	Zdata.hns24[0] = hns[4];
	Hdata.hns72_24 =  M_TO_CM((Zdata.hns24[0] + Zdata.hns24[48] + Zdata.hns24[96]) / cos_sl);

	if ( enforce_measured_snow_heights ) {
		// Check for CONSISTENCY and set indicator to -1 if necessary.
		//   The following is a cheap fix for InfoManager consistency if the deviation
		//   is small, i.e. most probably caused by the smoothing
		//   (4 cm for hns6, 10 cm for hns24).
		if ( ((Hdata.hns6 > 1.0)  &&  (Hdata.hns6 < M_TO_CM(d_hs6))) || ((Hdata.hns6 <= 1.0) && (M_TO_CM(d_hs6) > 3.0)) ) {
			if ( (Hdata.hns6 - M_TO_CM(d_hs6)) < 4.0 ) {
				Hdata.hns6 = M_TO_CM(d_hs6) * 1.1;
				Hdata.wc6 = Hdata.hns6 * 0.8;
				Hdata.hns3 = Hdata.hns6 * 0.7;
				Hdata.wc3 = Hdata.hns3 * 0.8;
				Hdata.hns12 = MAX (Hdata.hns6, Hdata.hns12);
				Hdata.wc12 = Hdata.hns12 * 0.8;
			} else {
				Hdata_ind.hns3  = -1;
				Hdata_ind.hns6  = -1;
				Hdata_ind.hns12 = -1;
				Hdata_ind.wc3   = -1;
				Hdata_ind.wc6   = -1;
				Hdata_ind.wc12  = -1;
			}
		}
		if ( ((Hdata.hns24 > 2.0)  &&  (Hdata.hns24 < M_TO_CM(d_hs24))) || ((Hdata.hns24 <= 2.0) && (M_TO_CM(d_hs24) > 5.0)) ) {
			if ( (Hdata.hns24 - M_TO_CM(d_hs24)) < 10.0) {
				Zdata.hns24[0] = d_hs24 * 1.1;
				Hdata.hns24 = M_TO_CM(Zdata.hns24[0]);
				Hdata.wc24 = Hdata.hns24 * 0.8;
				Hdata.hns72 = MAX (Hdata.hns24 * 1.2, Hdata.hns72);
				Hdata.wc72 = Hdata.hns72 * 0.8;
				Hdata.hns72_24 = M_TO_CM((Zdata.hns24[0] + Zdata.hns24[48] + Zdata.hns24[96]) / cos_sl);
			} else {
				Hdata_ind.hns24 = -1;
				Hdata_ind.hns72 = -1;
				Hdata_ind.hns72_24 = -1;
				Hdata_ind.wc24 = -1;
				Hdata_ind.wc72 = -1;
			}
		}
	}

	// INSTANTANEOUS DEWPOINT DEFICIT between TSS and Td(air)
	if (research_mode){
		Hdata.dewpt_def = K_TO_C(Xdata.Ndata[Xdata.getNumberOfNodes()-1].T) 
			- RhtoDewPoint(Mdata.rh, K_TO_C(Mdata.ta), force_rh_water);
	} else {
		Hdata.dewpt_def = calcDewPointDeficit(Mdata.ta, Xdata.Ndata[Xdata.getNumberOfNodes()-1].T, Mdata.rh);
	}

	if ( !((Hdata.dewpt_def > -50.) && (Hdata.dewpt_def < 50.)) ) {
		Hdata_ind.dewpt_def = -1;
	}

	// HOAR SIZE, size in mm assuming DENSITY_HOAR_SURF at surface
	Hdata.hoar_size = M_TO_MM(Xdata.Ndata[nE].hoar / density_hoar_surf);
	// Check for lower size limit
	if ( Hdata.hoar_size <= min_size_hoar_surf ) {
		Hdata.hoar_size = 0.;
	}
	if ( !((Hdata.hoar_size >= 0.) && (Hdata.hoar_size < 100.)) ) {
		Hdata_ind.hoar_size = -1;
	}
	// HOAR INDEX (24h and 6h), mass in kg m-2
	Hdata.hoar_ind24 = hoarIndex(&(Zdata.hoar24[0]), Sdata.hoar, 24, 1);
	if ( !((Hdata.hoar_ind24 > -10.) && (Hdata.hoar_ind24 < 10.)) ) {
		Hdata_ind.hoar_ind24 = -1;
	}
	Hdata.hoar_ind6  = hoarIndex(&(Zdata.hoar24[0]), Sdata.hoar,  6, 0);
	if ( !((Hdata.hoar_ind6 > -10.) && (Hdata.hoar_ind6 < 10.)) ) {
		Hdata_ind.hoar_ind6 = -1;
	}

	if ( !((Sdata.drift >= 0.0) && (Sdata.drift < 5.0)) ) {
		Sdata.drift = 0.0;
		Hdata_ind.wind_trans = -1;
		Hdata_ind.wind_trans24 = -1;
	}
	Hdata.wind_trans = driftIndex(&(Zdata.drift24[0]), Sdata.drift, Hazard::wind_slab_density, 6, 1) / cos_sl;
	if ( !((Hdata.wind_trans >= 0.) && (Hdata.wind_trans <= 6.*Hazard::maximum_drift)) ) {
		Hdata_ind.wind_trans = -1;
	}
	Hdata.wind_trans24 = driftIndex(&(Zdata.drift24[0]), Sdata.drift, Hazard::wind_slab_density, 24, 0) / cos_sl;
	if ( !((Hdata.wind_trans24 >= 0.) && (Hdata.wind_trans24 <= 24.*Hazard::maximum_drift)) ) {
		Hdata_ind.wind_trans24 = -1;
	}

	// SWE, total liquid water content and runoff, computed outside of Hazard.c
	Hdata.swe = Sdata.mass[SN_SURFACE_DATA::MS_SWE];
	Hdata.tot_lwc = Sdata.mass[SN_SURFACE_DATA::MS_WATER];
	Hdata.runoff = Sdata.mass[SN_SURFACE_DATA::MS_RUNOFF];

	// Profile type
	if ( (Xdata.S_class1 <= 10) && (Xdata.S_class1 >= 0) ) {
		Hdata.stab_class1 = Xdata.S_class1;
	} else {
		Hdata_ind.stab_class1 = -1;
	}

	// Stability class
	if ( (Xdata.S_class2 <= 5) && (Xdata.S_class2 >= 1) ) {
		Hdata.stab_class2 = Xdata.S_class2;
	} else {
		Hdata_ind.stab_class2 = -1;
	}
	// Stability index: Deformation index
	if ( (Xdata.S_d < (Stability::max_stability + Constants::eps)) && (Xdata.S_d > 0.) ) {
		Hdata.stab_index1 = Xdata.S_d;
		Hdata.stab_height1 = M_TO_CM(Xdata.z_S_d / cos_sl);
	} else {
		Hdata_ind.stab_index1 = -1;
		Hdata_ind.stab_height1 = -1;
	}
	// Natural stability index Sn38
	if ( (Xdata.S_n < (Stability::max_stability + Constants::eps)) && (Xdata.S_n > 0.) ) {
		Hdata.stab_index2 = Xdata.S_n;
		Hdata.stab_height2 = M_TO_CM(Xdata.z_S_n / cos_sl);
	} else {
		Hdata_ind.stab_index2 = -1;
		Hdata_ind.stab_height2 = -1;
	}
	// Skier stability index Sk38
	if ( (Xdata.S_s < (Stability::max_stability + Constants::eps)) && (Xdata.S_s > 0.) ) {
		Hdata.stab_index3 = Xdata.S_s;
		Hdata.stab_height3 = M_TO_CM(Xdata.z_S_s / cos_sl);
	} else {
		Hdata_ind.stab_index3 = -1;
		Hdata_ind.stab_height3 = -1;
	}
	// Structural stability index SSI
	if ( (Xdata.S_4 < (Stability::max_stability + Constants::eps)) && (Xdata.S_4 > 0.) ) {
		Hdata.stab_index4 = Xdata.S_4;
		Hdata.stab_height4 = M_TO_CM(Xdata.z_S_4 / cos_sl);
	} else {
		Hdata_ind.stab_index4 = -1;
		Hdata_ind.stab_height4 = -1;
	}
	// ???Index???
	if ( (Xdata.S_5 < (Stability::max_stability + Constants::eps)) && (Xdata.S_5 > 0.) ) {
		Hdata.stab_index5 = Xdata.S_5;
		Hdata.stab_height5 = M_TO_CM(Xdata.z_S_5 / cos_sl);
	} else {
		Hdata_ind.stab_index5 = -1;
		Hdata_ind.stab_height5 = -1;
	}

	// Surface crust [type == 772] calculated for southerly aspect outside Hazard.c

	// Energy input ... (kJ m-2)
	if ( nE > Xdata.SoilNode ) {
		Hdata.en_bal = Sdata.dIntEnergy * hazard_steps_between / 1000.;
	} else {
		Hdata.en_bal = ( ((Sdata.qw + Sdata.lw_net + Sdata.qs + Sdata.ql + Sdata.qr) / nAvg) * sn_dt * hazard_steps_between ) / 1000.;
  }
	if ( !((Hdata.en_bal > -3000.) && (Hdata.en_bal < 3000.)) ) {
		Hdata_ind.en_bal = -1;
	}

	// Net SW energy at surface (kJ m-2)
	if ( Sdata.sw_in > 0. ) {
		Hdata.sw_net = (Sdata.qw * sn_dt * hazard_steps_between) / 1000.;
		if ( !((Hdata.sw_net > -3000.) && (Hdata.sw_net < 3000.)) ) {
			Hdata_ind.sw_net = -1;
		}
	} else {
		Hdata.sw_net = 0.;
	}

	// Snow temperatures t_top1 and t_top2 in degC at 5 cm and 10 cm below the surface, respectively
	h_top1 = hs - 0.05;
	Hdata.t_top1 = getModelledTemperature(h_top1, Xdata);
	if ( !((Hdata.t_top1 > -50.) && (Hdata.t_top1 <= 0.)) ) {
		Hdata_ind.t_top1 = -1;
	}
	h_top2 = hs - 0.10;
	Hdata.t_top2 = getModelledTemperature(h_top2, Xdata);
	if ( !((Hdata.t_top2 > -50.) && (Hdata.t_top2 <= 0.)) ) {
		Hdata_ind.t_top2 = -1;
	}
}

/*
 * End of Hazard.c
 */
