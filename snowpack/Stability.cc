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

#include <snowpack/Stability.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/

const double Stability::psi_ref = 38.0; ///< Reference slope angle
const double Stability::max_stability = 6.; ///< Upper stability limit

///Minimum slab thickness for natural and deformation stability index (m) 
const double Stability::minimum_slab = 0.1; 

///The first GROUND_ROUGH m of snow will not be unstable due to ground roughness
const double Stability::ground_rough = 0.2; 

///MIN_DEPTH_SSI m of snow must be left after discarding penetration depth
const double Stability::min_depth_ssi = 0.1;

///Skiers will not trigger failures SKIER_DEPTH m below penetration depth
const double Stability::skier_depth = 1.0;

///Minimum thickness for a supporting melt-freeze crust (perp to slope, in m)
const double Stability::min_thick_crust = 0.03;

///Defines regression model for surface hoar shear strength
const int Stability::sh_mod = 2;

///Maximum number of structural instabilities looked at ("lemons")
const int Stability::nmax_lemon = 2;

/**
 * @brief Defines classification scheme for snow profiles
 * - 0: Based on Skier Stability Index (SSI) thresholds (3 classes)
 * - 1: Based on re-analysis by Schweizer/Bellaire of SSI and SK38 (2006)
 * - 2: Nov 2007: re-analysis after recalibration of settling (see rev 250/251)
 * - 3: According to Schweizer and Wiesinger (5 classes)
 */
const int Stability::prof_classi = 2;

map<string, StabMemFn> Stability::mapHandHardness;
map<string, StabFnShearStrength> Stability::mapShearStrength;
const bool Stability::__init = Stability::initStaticData();

bool Stability::initStaticData()
{
	mapHandHardness["DEFAULT"]  = &Stability::st_HandHardnessDEFAULT;
	mapHandHardness["ASARC"]    = &Stability::st_HandHardnessASARC;
	
	mapShearStrength["DEFAULT"] = &Stability::st_ShearStrengthDEFAULT;
	mapShearStrength["NIED"]    = &Stability::st_ShearStrengthSTRENGTH_NIED;

	return true;
}

/************************************************************
 * non-static section                                       *
 ************************************************************/

Stability::Stability(const mio::Config& i_cfg) : cfg(i_cfg) 
{
	string tmp_strength_model = cfg.get("STRENGTH_MODEL", "Parameters");
	strength_model = tmp_strength_model;

	string tmp_hardness_model = cfg.get("HARDNESS_MODEL", "Parameters");
	hardness_model = tmp_hardness_model;

	IOUtils::toUpper(strength_model);
	IOUtils::toUpper(hardness_model);

	map<string, StabMemFn>::const_iterator it1 = mapHandHardness.find(hardness_model);
	if (it1 == mapHandHardness.end()) throw InvalidArgumentException("Unknown hardness model: "+hardness_model, AT);

	map<string, StabFnShearStrength>::const_iterator it2 = mapShearStrength.find(strength_model);
	if (it2 == mapShearStrength.end()) throw InvalidArgumentException("Unknown strength model: "+strength_model, AT);

	//To build a sandwich with a non-snow layer (plastic or wood chips) on top; 
	//originally introduced for snow farming applications
	plastic = cfg.get("PLASTIC", "Parameters");

	// Density of BURIED surface hoar (kg m-3), default: 125./ Antarctica: 200.
	density_hoar_buried = cfg.get("DENSITY_HOAR_BURIED", "Parameters");
}


/**
 * @brief Assign hardness to snow types according to density, Swiss version
 * @author Implemented by C. Fierz: Regression by Sascha Bellaire January 2005 (all types except MFcr).
 * The original Swiss regression has been modified for PP, RG and FC to get a better agreement
 * to observed hardness. If there is a new settlement formulation in future, and therefore a
 * better agreement to observed density,it must be checked whether the new or the original
 * regression is the better choice. (18 September 2005; Fierz / S. Bellaire)
 * Original regression values are added as comments where needed.
 * @param *Edata
 * @return hand hardness index (1)
 */
double Stability::st_HandHardnessDEFAULT(const SN_ELEM_DATA& Edata)
{
	int    F1, F2, F3; // grain shape
	double hardness;
	double gsz;

	// Dereference some values
	gsz = 2.*Edata.rg;

	// Decompose type in its constituents
	qr_TypeToCode(&F1, &F2, &F3, Edata.type);

	if ( (Edata.mk%100) < 20 ) { // all types except MFcr (hardness 5)
		double A, B;
		switch( F1 ) {
			case 0: { // Graupel PPgp; introduced by Yamaguchi & Fierz, Feb 2004
				A = 0.0078;
				B = 0.0105;
				break;
			}
			case 1: { // Precipitation Particles PP; ori: A = 0.7927; B = 0.0038;
				A = 0.7927;
				B = 0.0036;
				break;
			}
			case 2: { // Decomposing and Fragmented precipitation particles DF
				A = 0.4967;
				B = 0.0074;
				break;
			}
			case 3: { // Rounded Grains RG; ori: A = 0.2027; B = 0.0092;
				A = 0.2027;
				B = 0.0072;
				break;
			}
			case 4: { // Faceted Crystals FC; ori: A = 0.3867; B = 0.0071;
				A = 0.3867;
				B = 0.0083;
				break;
			}
			case 5: { // Depth Hoar DH
				A = -0.00249;
				B = 0.0072;
				break;
			}
			case 6: { // Surface hoar SH; empirical: index 1 to 2 from DENSITY_HOAR_BURIED to 250 kg m-3
				A = 1. - density_hoar_buried/(250. - density_hoar_buried);
				B = 1./(250. - density_hoar_buried);
				break;
			}
			case 7: { // Melt Forms MF
				A = 0.5852;
				B = 0.0056;
				break;
			}
			case 8: { // Ice layer IFil
				A = 6.;
				B = 0.;
				break;
			}
			case 9: { // Rounding faceted particles FCxr
				A = -0.5226;
				B = 0.0104;
				break;
			}
			default: {
				A = NODATA;
				B = 0.;
				break;
			}
		}
		hardness = A + B*Edata.Rho;
		// Large surface hoar stays longer unstable! 1 dec 2007 (sb)
		if ( (F1 == 6) && (gsz >= 5.) ) {
			hardness = 1;
		} else if ((F1 == 6 ) && (gsz < 5.)) {
			hardness = MIN(hardness, 2.);
		}
	} else if ( Edata.theta[ICE] <= 0.7 ) { // Melt-freeze crust MFcr
		double res_wat_cont;
		res_wat_cont = lw_SnowResidualWaterContent(Edata.theta[ICE]);
		if ( Edata.theta[WATER] < 0.3*res_wat_cont ) {
			hardness = 5.;
		} else if ( Edata.theta[WATER] < 0.6*res_wat_cont ) {
			hardness = 4.5;
		} else if ( Edata.theta[WATER] < 0.85*res_wat_cont ) {
			hardness = 4.;
		} else {
			hardness = 3.;
		}
	} else { // Ice Formations IF
		hardness = 6.;
	}
	// Limit to range {1, 6}
	hardness = MAX(1., MIN(6., hardness));
	return(hardness);
}

/**
 * @brief Assign hand hardness to snow types according to density and grain size, original Canadian version
 * @author Implemented by C. Fierz: Regression from ASARC database by Bruce Jamieson on 2002-08-14
 * @param *Edata
 * @return hand hardness index (1)
 */
double Stability::st_HandHardnessASARC(const SN_ELEM_DATA& Edata)
{
	int    F1, F2, F3;
	double A=0., B=0., C=0.;
	double hardness;
	double gsz;

	// Dereference some values
	gsz     = 2.*Edata.rg;

	// Decompose type in its constituents
	qr_TypeToCode(&F1, &F2, &F3, Edata.type);

	// all types except MFcr
	if( Edata.mk%100 < 20 ) {
		switch ( F1 ) {
			case 0: { // Graupel PPgp; empirical!
				A = 1.5;
				B = 0.;
				C = 0.;
				break;
			}
			case 1: { // Precipitation Particles PP
				A = 0.45;
				B = 0.0068;
				C =  0.;
				break;
			}
			case 2: { // Decomposing and Fragmented precipitation particles DF
				A =  0.;
				B = 0.0140;
				C =  0.;
				break;
			}
			case 3: { // Rounded Grains RG
				A =  1.94;
				B = 0.0073;
				C = -0.192;
				break;
			}
			case 4: {
				if ( F2 != 9 ) { // Faceted Crystals FC
					A =  0.;
					B = 0.0138;
					C = -0.284;
				} else { // Rounding faceted particles FCxr, because F1=9 does not occur in SNOWPACK
					A =  1.29;
					B = 0.0094;
					C = -0.350;
				};
				break;
			}
			case 5: {
				if ( gsz > 1.5 ) { // Depth hoar DH, small dataset (N=41) !!!
					A = -0.80;
					B = 0.0150;
					C = -0.140;
				} else { // use FC values for small depth hoar grains
					A =  0.00;
					B = 0.0138;
					C = -0.284;
				};
				break;
			}
			case 6: { // Surface hoar SH; empirical: index 1 to 2 from DENSITY_HOAR_BURIED to 250 kg m-3
				A = 1. - density_hoar_buried/(250. - density_hoar_buried);
				B = 1./(250. - density_hoar_buried);
				C = 0.;
				break;
			}
			case 7: { // Melt Forms MF
				if ( Edata.theta[WATER] < 0.001 ) { // MF, dry
					A = 2.14;
					B = 0.0048;
					C =  0.;
				} else { // MF, wet (LWC > 0.)
					A = 3.00;
					B = 0.0000;
					C =  0.;
				};
				break;
			}
			case 8: { // Ice layer IFil
				A =  6.;
				B = 0.;
				C =  0.;
				break;
			}
			case 9: { // Rounding faceted particles FCxr
				A =  1.29;
				B = 0.0094;
				C = -0.350;
				break;
			}
			default: {
				A = NODATA;
				B = 0.;
				C = 0.;
				break;
			}
		}
	} else if ( Edata.theta[ICE] <= 0.7 ) { // Melt-freeze crust MFcr
		double res_wat_cont;
		res_wat_cont = lw_SnowResidualWaterContent(Edata.theta[ICE]);
		if ( Edata.theta[WATER] < 0.3*res_wat_cont ) {
			A = 5.;
		} else if ( Edata.theta[WATER] < 0.6*res_wat_cont ) {
			A = 4.5;
		} else if ( Edata.theta[WATER] < 0.85*res_wat_cont ) {
			A = 4.;
		} else {
			A = 3.;
		}
	} else { // Ice Formations IF
		A = 6.;
	}
	hardness = A + B*Edata.Rho + C*gsz;
	if ( F1 == 6 ) {
		hardness = MIN(hardness, 2.);
	}
	// Limit to range {1, 6}
	hardness = MAX(1., MIN(6., hardness));
	return(hardness);
}

/*
 * START OF STABILITY SECTION
*/
/**
 * @brief Returns the critical stress state of a layer given the temperature and plastic strain rate.
 * @param epsNeckDot Neck strain rate (s-1)
 * @param Ts Temperature of layer (K)
 * @return Critical stress (Pa)
 */
double Stability::st_CriticalStress(const double& epsNeckDot, const double& Ts)
{
	double sigBrittle=1.e7;   // Brittle fracture stress of ice (Pa)
	double Pm;                // Hydrostatic pressure that induces melting (Pa)
	double phi;               // Function of strain rate dependent failure surface
	double epsa;              // Absolute value of plastic strain rate

	const double C1=-6.6249;     // Constant
	const double C2=6.0780e-2;   // Constant
	const double C3=-1.3380e-4;  // Constant
	const double P1=70.000;      // Constant (Pa)

	// Find the rate dependent friction angle phi
	epsa = fabs(epsNeckDot);
	phi = P1*pow(epsa, 0.23)*Constants::pi/180.;

	// Hydrostatic melting pressure
	Pm = (C1 + C2*(Ts) + C3*(Ts)*(Ts)) * 1.e9;
	// Find the critical stress
	return(Pm * tan(phi) * sqrt(1. - (Pm/(Pm + sigBrittle)))); // TODO check that argument of sqrt is correctly written
} // End of st_CriticalStress

/**
 * @brief Returns the layer stability index
 * NOTE:  The intra-layer stability criteria is given by the ratio S_f = S_c/S_n where
 * S_n is the neck stress and S_c is the critical stress.  The critical stress is determined
 * in the function st_CriticalStress. This function might get a little more involved as
 * time goes on.
 * @param *Edata
 * @return Deformation rate index
 */
double Stability::st_DeformationRateIndex(const SN_ELEM_DATA& Edata)
{
	const double eps1Dot = 1.76e-7; // Unit strain rate (at stress = 1 MPa) (s-1)
	const double sig1 = 0.5e6;      // Unit stress from Sinha's formulation (Pa)
	const double sig = -Edata.C;   // Overburden stress, that is, absolute value of Cauchy stress (Pa)
	double sigNeck;                 // Neck stress (Pa)
	double epsNeckDot;              // Total strain rate in the neck (s-1)
	double Te;                      // Element temperature (K)

	Te = MIN(Edata.Te, Constants::melting_tk);
	// If you have less than 5% ice then say you know you have something unstable
	if ( Edata.theta[ICE] < 0.05 ) {
		return(0.1);
	}
	// First find the absolute neck stress
	sigNeck = SnLaws::calcNeckStressEnhancement(Edata) * (sig);
	// Now find the strain rate in the neck
	epsNeckDot =  eps1Dot * lwsn_SnowViscosityTemperatureTerm(Te) * (sigNeck/sig1)*(sigNeck/sig1)*(sigNeck/sig1);
	// Return the stability index
	return(MAX(0.1, MIN(st_CriticalStress(epsNeckDot, Te) / sigNeck, 6.)));
} // End of st_DeformationRateIndex

/**
 * @brief Initializes strenght parameters
 * @param STpar
 * @param *Xdata
 * @param SIdata
 * @param psi_ref Reference slope angle (deg)
 */
void Stability::initStability(const double& psi_ref, StabilityData& STpar,
                              SN_STATION_DATA& Xdata, vector<InstabilityData>& SIdata)
{
	int e;

	STpar.cos_psi_ref = cos(psi_ref*Constants::pi/180.);
	STpar.sin_psi_ref = sin(psi_ref*Constants::pi/180.);
	STpar.sig_n = 999.;
	STpar.sig_s = 999.;
	STpar.alpha_max = 54.3*Constants::pi/180.; // alpha_max(38.) = 54.3 deg (J. Schweizer, IB 712, SLF)

	for (e=Xdata.getNumberOfNodes()-1; e>=Xdata.SoilNode; e--) {
		SIdata[e].ssi        = Stability::max_stability;
		Xdata.Ndata[e].S_n  = Stability::max_stability;
		Xdata.Ndata[e].S_s  = Stability::max_stability;
		if (e < Xdata.getNumberOfNodes()-1 ) {
			Xdata.Edata[e].S_dr = Stability::max_stability;
		}
	}
}

/**
 * @brief Returns the skier's penetration depth (Adapted from Canadian parameterization)
 * @param *Xdata
 */
double Stability::st_PenetrationDepth(const SN_STATION_DATA& Xdata)
{
	double rho_Pk=0., dz_Pk=0.;            // Penetration depth Pk, from mean slab density
	double cos_sl;                         // Cosine of slope angle
	double top_crust=0., thick_crust=0.;   // Crust properties
	int    crust=0;                        // Checks for crust
	int    e=Xdata.getNumberOfElements()-1, e_crust=-99; // Counters

	cos_sl = cos(Xdata.SlopeAngle);
	while ( (e >= Xdata.SoilNode) && (((Xdata.cH - (Xdata.Ndata[e].z + Xdata.Ndata[e].u))/cos_sl < 0.3)) ) {
		rho_Pk += Xdata.Edata[e].Rho*Xdata.Edata[e].L;
		dz_Pk  += Xdata.Edata[e].L;
		// Test for strong mf-crusts MFcr.
		// Look for the first (from top) with thickness perp to slope > 3cm
		if ( !crust ) {
			if ( (Xdata.Edata[e].mk%100 >= 20) && (Xdata.Edata[e].Rho > 500.) ) {
				if ( e_crust < 0 ) {
					e_crust = e;
					top_crust = (Xdata.Ndata[e+1].z + Xdata.Ndata[e+1].u)/cos_sl;
					thick_crust += Xdata.Edata[e].L;
				} else if ( ((e_crust - e) < 2) ) {
					thick_crust += Xdata.Edata[e].L;
					e_crust = e;
				}
			} else if ( e_crust > 0 ) {
				if ( thick_crust > Stability::min_thick_crust ) {
					crust = 1;
				} else {
					e_crust = INODATA;
					top_crust = 0.;
					thick_crust = 0.;
				}
			}
		}
		e--;
	}
	rho_Pk /= (dz_Pk + 1.e-12);

  // NOTE Pre-factor 0.8 introduced May 2006 by S. Bellaire
	return(MIN(0.8*43.3/(rho_Pk + 1.e-12), (Xdata.cH/cos_sl - top_crust)));
}

/**
 * @brief Computes normal and shear stresses (kPa) reduced to psi_ref
 * @param STpar
 * @param stress Overload perpendicular to slope (Pa)
 * @param cos_sl Cosine of slope angle (1)
 */
void Stability::calcReducedStresses(const double& stress, const double& cos_sl, StabilityData& STpar)
{
	STpar.sig_n = -stress*(STpar.cos_psi_ref*STpar.cos_psi_ref)/(cos_sl*cos_sl)/1000.;
	STpar.sig_s = (STpar.sig_n)*STpar.sin_psi_ref/STpar.cos_psi_ref;
}

/**
 * @brief DEFAULT: Estimates the critical shear stress based on appropriate parameterisations
 * @param *Edata Xdata->Edata[e+1]
 * @param *Ndata Xdata->Ndata[e+1]
 * @param STpar
 * @param cH Calculated height of snow (m)
 * @param cos_sl Cosine of slope angle (1)
 * @param date
 * @return return false on error, true otherwise
 */
bool Stability::st_ShearStrengthDEFAULT(const double& cH, const double& cos_sl, const mio::Date& date,
                                        SN_ELEM_DATA& Edata, SN_NODE_DATA& Ndata, StabilityData& STpar)
{
	int    F1, F2, F3;                    // Grain shape
	double Sig_cC, Sig_c2, Sig_c3;        // Critical shear stress (kPa)
	double phi;                           // Normal load correction
	double rho_ri;                        // Snow density relative to ice

	// Snow density relative to ice
	rho_ri = Edata.Rho/Constants::density_ice;
	// Determine majority grain shape
	qr_TypeToCode(&F1, &F2, &F3, Edata.type);

	// Determine critical shear stress of element (kPa)
	// 1. Conway
	Sig_cC = 19.5*rho_ri*rho_ri;
		// 2. Grain Type dependent mostly from Jamieson,
		//    Ann. Glaciol., 26, 296-302 (2001) and Ann. Glaciol., 32, 59-69 (1998)
	phi = 0.;
	Sig_c2 = -1.0;
	Sig_c3 = -1.0;
	switch( F1 ) {
		case 0: // Graupel, from O. Abe, Ann. Glaciol. 38 (2004), size-effect corrected
			Sig_c2 = 0.65*(82.*pow(rho_ri, 2.8));
			phi = 0.08*Sig_c2 + 0.224;
			break;
		case 1: // PP
			Sig_c2 = 2.85*exp(1.13*log(rho_ri));
			phi = 0.08*Sig_c2 + 0.056 + 0.022*STpar.sig_n;
			break;
		case 2: // DF
			Sig_c2 = 8.75*exp(1.54*log(rho_ri));
			phi = 0.08*Sig_c2 + 0.224;
			break;
		case 3: // RG
			Sig_c2 = 7.39*exp(1.20*log(rho_ri));
			phi = 0.08*Sig_c2 + 0.224;
			break;
		case 6: // SH
			switch(Stability::sh_mod) {
				case 0: // original T. Chalmers
					Sig_c2 = 0.336 + 0.0139*(date.getJulianDate() - Edata.date.getJulianDate()) +
							1.18*STpar.sig_n/(STpar.cos_psi_ref*STpar.cos_psi_ref) - 0.625*(cH -
							(Ndata.z + Ndata.u))/cos_sl + 0.0804 *
							cH/cos_sl - 28.7*Edata.L/cos_sl +
							0.0187*K_TO_C(Edata.Te) + 0.0204*Edata.rg;
					break;
				case 1: // original T. Chalmers & accounting for Emin as 2*rg (ml 13 Feb 2003)
					Sig_c2 = 0.336 + 0.0139*(date.getJulianDate() - Edata.date.getJulianDate()) +
							1.18*STpar.sig_n/(STpar.cos_psi_ref*STpar.cos_psi_ref) - 0.625*(cH -
							(Ndata.z + Ndata.u))/cos_sl + 0.0804 *
							cH/cos_sl - 28.7*Edata.L/cos_sl +
							0.0187*K_TO_C(Edata.Te) + 0.0204*2.*Edata.rg;
					break;
				case 2: // New regression by Bruce Jamieson w/o Emin (14 Feb 2003)
					Sig_c2 = 0.429 + 0.0138*(date.getJulianDate() - Edata.date.getJulianDate()) +
							1.12*STpar.sig_n/(STpar.cos_psi_ref*STpar.cos_psi_ref) - 0.596*(cH -
							(Ndata.z + Ndata.u))/cos_sl + 0.0785 *
							cH/cos_sl - 27.1*Edata.L/cos_sl +
							0.0202*K_TO_C(Edata.Te);
					break;
				default:
					Sig_c2 = 1.0;
					break;
			}
			Sig_c2 = MAX(0.1, Sig_c2);
			Sig_c3 = 84.*exp(2.55*log(rho_ri));
			break;
		case 7: // MF
			Sig_c2 = 21.*exp(1.24*log(rho_ri));
			phi = 0.08*Sig_c2 + 0.224;
			break;
		default: // FC, DH, FCmx
			Sig_c2 = 18.5*exp(2.11*log(rho_ri));
			Sig_c3 = 1.36*exp(0.55*log(STpar.sig_n/STpar.cos_psi_ref));
			// phi = 0.08*Sig_c2 + 0.224; // Correction not used by B. Jamieson, see IGS03
			break;
	}

		// Hack for MFCs
	if ( Edata.mk >= 20 ) {
		Sig_c2 = Sig_c3 = 4.;
	}

	// Final assignements
	Edata.s_strength = Sig_c2;
	Sig_c2 = MIN(Sig_c2, STpar.strength_up);
	STpar.Sig_c2 = Sig_c2;
	STpar.phi = phi;

	// Warning message may be enabled for large differences in snow shear stength models
	if ( 0 && (((fabs(Sig_c2-Sig_cC)/Sig_cC) > 10.) || ((Sig_c3 > 0.) && ((fabs(Sig_c3-Sig_cC)/Sig_cC > 10.)))) ) {
		prn_msg( __FILE__, __LINE__, "wrn", date.getJulianDate(),"Large difference in Snow Shear Stength (type=%d)", F1);
		prn_msg(__FILE__, __LINE__, "msg-", -1., "Conway: %lf Sig_c2: %lf Sig_c3: %lf\n", Sig_cC, Sig_c2, Sig_c3);
		return false;
	} else {
		return true;
	}
} // End st_ShearStrengthDEFAULT

/**
 * @brief STRENGTH_NIED: Estimates the critical shear stress based on appropriate parameterisations adapted for Japan
 * @param *Edata Xdata->Edata[e+1]
 * @param *Ndata Xdata->Ndata[e+1]
 * @param STpar
 * @param cH Calculated height of snow (m)
 * @param cos_sl Cosine of slope angle (1)
 * @param date
 * @return return false on error, true otherwise
 */
bool Stability::st_ShearStrengthSTRENGTH_NIED(const double& cH, const double& cos_sl, const mio::Date& date,
                                              SN_ELEM_DATA& Edata, SN_NODE_DATA& Ndata, StabilityData& STpar)
{
	int    F1, F2, F3;             // Grain shape
	double Sig_cC, Sig_c2, Sig_c3; // Critical shear stress (kPa)
	double Sig_ET, Sig_DH;         //NIED (H. Hirashima)
	double phi;                    // Normal load correction
	double rho_ri;                 // Snow density relative to ice

	// Snow density relative to ice
	rho_ri = Edata.Rho/Constants::density_ice;
	// Determine majority grain shape
	qr_TypeToCode(&F1, &F2, &F3, Edata.type);

	// Determine critical shear stress of element (kPa)
  // 1. Conway
	Sig_cC = 19.5*rho_ri*rho_ri;
		// 2. Grain Type dependent mostly from Jamieson,
		//    Ann. Glaciol., 26, 296-302 (2001) and Ann. Glaciol., 32, 59-69 (1998)
	phi = 0.;
	Sig_c2 = -1.0;
	Sig_c3 = -1.0;
	switch ( F1 ) {
		case 0: // Graupel, from O. Abe, Ann. Glaciol. 38 (2004), size-effect corrected
			Sig_c2 = 0.65*(82.*pow(rho_ri, 2.8));
			phi = 0.08*Sig_c2 + 0.224;
			break;
		case 1: // PP //NIED (H. Hirashima)
			Sig_c2= 9.4*0.0001*pow(Edata.theta[ICE]*Constants::density_ice,2.91)*exp(-0.235*Edata.theta[WATER]*100.)/1000.;
			phi = 0.08*Sig_c2 + 0.056 + 0.022*STpar.sig_n;
			break;
		case 2: case 3: // DF & RG //NIED (H. Hirashima)
			Sig_c2= 9.4*0.0001*pow(Edata.theta[ICE]*Constants::density_ice,2.91)*exp(-0.235*Edata.theta[WATER]*100.)/1000.;
			phi = 0.08*Sig_c2 + 0.224;
			break;
		case 6: // SH
			switch(Stability::sh_mod) {
				case 0: // original T. Chalmers
					Sig_c2 = 0.336 + 0.0139*(date.getJulianDate() - Edata.date.getJulianDate()) +
							1.18*STpar.sig_n/(STpar.cos_psi_ref*STpar.cos_psi_ref) - 0.625*(cH -
							(Ndata.z + Ndata.u))/cos_sl + 0.0804 *
							cH/cos_sl - 28.7*Edata.L/cos_sl +
							0.0187*K_TO_C(Edata.Te) + 0.0204*Edata.rg;
					break;
				case 1: // original T. Chalmers & accounting for Emin as 2*rg (ml 13 Feb 2003)
					Sig_c2 = 0.336 + 0.0139*(date.getJulianDate() - Edata.date.getJulianDate()) +
							1.18*STpar.sig_n/(STpar.cos_psi_ref*STpar.cos_psi_ref) - 0.625*(cH -
							(Ndata.z + Ndata.u))/cos_sl + 0.0804 *
							cH/cos_sl - 28.7*Edata.L/cos_sl +
							0.0187*K_TO_C(Edata.Te) + 0.0204*2.*Edata.rg;
					break;
				case 2: // New regression by Bruce Jamieson w/o Emin (14 Feb 2003)
					Sig_c2 = 0.429 + 0.0138*(date.getJulianDate() - Edata.date.getJulianDate()) +
							1.12*STpar.sig_n/(STpar.cos_psi_ref*STpar.cos_psi_ref) - 0.596*(cH -
							(Ndata.z + Ndata.u))/cos_sl + 0.0785 *
							cH/cos_sl - 27.1*Edata.L/cos_sl +
							0.0202*K_TO_C(Edata.Te);
					break;
				default:
					Sig_c2 = 1.0;
					break;
			}
			Sig_c2 = MAX (0.1, Sig_c2);
			Sig_c3 = 84.*exp(2.55*log(rho_ri));
			break;
		case 7: // MF //NIED (H. Hirashima)
			Sig_c2= 4.97*0.0001*pow(Edata.theta[ICE]*Constants::density_ice,2.91)*exp(-0.235*Edata.theta[WATER]*100.)/1000.;
			phi = 0.08*Sig_c2 + 0.224;
			break;
		default: // FC, DH, FCmx
			Sig_c2 = 18.5*exp(2.11*log(rho_ri));
			//Sig_c2 = 0.0391*exp(0.0141*Edata.theta[ICE]*Constants::density_ice); //NIED (H. Hirashima)
			Sig_c3 = 1.36*exp(0.55*log(STpar.sig_n/STpar.cos_psi_ref));
			// phi = 0.08*Sig_c2 + 0.224; // Correction not used by B. Jamieson, see IGS03
			break;
	}

	// Hack for MFCs; not used by //NIED (H. Hirashima)

	// Final assignements
	//NIED (H. Hirashima)
	Sig_ET = 9.4*0.0001*pow(Edata.theta[ICE]*Constants::density_ice,2.91)*exp(-0.235*Edata.theta[WATER]*100.)/1000.;
	Sig_DH = 2.3*0.0001*pow(Edata.theta[ICE]*Constants::density_ice,2.78)*exp(-0.235*Edata.theta[WATER]*100.)/1000.;
	Ndata.Sigdhf = Sig_ET - Edata.dhf*(Sig_ET - Sig_DH);
	Ndata.S_dhf = (Ndata.Sigdhf + phi*STpar.sig_n)/STpar.sig_s;
	// original SNOWPACK
	Edata.s_strength = Sig_c2;
	Sig_c2 = MIN(Sig_c2, STpar.strength_up);
	STpar.Sig_c2 = Sig_c2;
	STpar.phi = phi;

	// Warning message may be enabled for large differences in snow shear stength models
	if ( 0 && (((fabs(Sig_c2-Sig_cC)/Sig_cC) > 10.) || ((Sig_c3 > 0.) && ((fabs(Sig_c3-Sig_cC)/Sig_cC > 10.)))) ) {
		prn_msg( __FILE__, __LINE__, "wrn", date.getJulianDate(),"Large difference in Snow Shear Stength (type=%d)", F1);
		prn_msg(__FILE__, __LINE__, "msg-", -1., "Conway: %lf Sig_c2: %lf Sig_c3: %lf\n", Sig_cC, Sig_c2, Sig_c3);
		return false;
	} else {
		return true;
	}
} // End st_ShearStrengthSTRENGTH_NIED

/**
 * @brief Returns the natural stability index Sn
 * The classic natural stability index Sn, that is, the ratio of shear stress to shear strength (static)
 * @param STpar
 */
double Stability::st_NaturalStabilityIndex(const StabilityData& STpar)
{
	// Limit natural stability index to range {0.05, Stability::max_stability}
	return(MAX(0.05, MIN(((STpar.Sig_c2 + STpar.phi*STpar.sig_n)/STpar.sig_s), Stability::max_stability)));
}

/**
 * @brief Returns the skier stability index Sk reduced to psi_ref (usually 38 deg => Sk_38)
 * The classic skier stability index Sk(psi_ref), using P. Foehn's formula
 * (IAHS No162, 1987, p201) for the skier induced shear stress,
 * corrected for skier penetration depth Pk.
 * @param depth_lay Depth of layer to investigate (m)
 * @param STpar
 */
double Stability::st_SkierStabilityIndex(const double& depth_lay, const StabilityData& STpar)
{
	double delta_sig; // Skier contribution to shear stress (kPa) at psi_ref (usually 38 deg)

	if ( depth_lay > Constants::eps ) {
		delta_sig = (2.*0.5*cos(STpar.alpha_max)*sin(STpar.alpha_max)*sin(STpar.alpha_max)*sin(STpar.alpha_max + (STpar.psi_ref*Constants::pi/180.)));
		delta_sig /= Constants::pi*STpar.cos_psi_ref*depth_lay;
		// Limit skier stability index to range {0.05, Stability::max_stability}
		return(MAX(0.05, MIN(((STpar.Sig_c2 + STpar.phi*STpar.sig_n)/(STpar.sig_s + delta_sig)), Stability::max_stability)));
	} else {
		return(Stability::max_stability); // strictly speaking, Sk is not defined
	}
}

/**
 * @brief Returns the structural stability index SSI
 * Adds one lemon to Sk for each structural instability found, presently hardness and grain size differences
 * above a given threshold
 * @param Edata_low Xdata->Edata[e]
 * @param Edata_up Xdata->Edata[e+1]
 * @param Sk Skier stability index Sk (Xdata->Ndata[e+1].S_s)
 * @param SIdata [e+1]
 */
void Stability::setStructuralStabilityIndex(const SN_ELEM_DATA& Edata_low, const SN_ELEM_DATA& Edata_up,
                                            const double& Sk, InstabilityData& SIdata)
{
	const double thresh_dhard=1.5, thresh_dgsz=0.50; // Thresholds for structural instabilities

	SIdata.n_lemon = 0;
	SIdata.dhard = fabs(Edata_low.hard - Edata_up.hard);
	if ( SIdata.dhard > thresh_dhard ) {
		SIdata.n_lemon++;
	}
	SIdata.dgsz = 2.*fabs(Edata_low.rg - Edata_up.rg);
	if ( SIdata.dgsz > thresh_dgsz ) {
		SIdata.n_lemon++;
	}
	// Skier Stability Index (SSI)
	SIdata.ssi = (Stability::nmax_lemon - SIdata.n_lemon) + Sk;
	// Limit stability indices to range {0.05, Stability::max_stability}
	SIdata.ssi = MAX(0.05, MIN (SIdata.ssi, Stability::max_stability));
}

/**
 * @brief Returns the Profile Stability Classification (Schweizer-Wiesinger Method)
 * @param Xdata
 * @return false if error, true otherwise
 */
bool Stability::classifyProfileStability(SN_STATION_DATA& Xdata)
{
	int S = 5, S0, e, nE, i, i_weak, count=0;
	int F1, F2, F3;                            // Grain shape
	double mH = 0., thH, maxH = 0., minH = 7.; // Mean Hardness and Threshold
	double mH_u, delta_H;
	double h_Slab;
	double cos_sl;

	SN_ELEM_DATA *EMS;   // Avoids dereferencing the element data pointer

	// Dereference the element pointer containing micro-structure data
	EMS = &Xdata.Edata[0]; nE = Xdata.getNumberOfElements();
	vector<SN_NODE_DATA>& NDS = Xdata.Ndata;

	// Initialize
	S0 = S;
	cos_sl = cos(Xdata.SlopeAngle);
	h_Slab = EMS[nE-1].L/cos_sl;

	// Classify only for Snowpacks thicker than Stability::minimum_slab (vertically)
	if ( (NDS[nE].z+NDS[nE].u)/cos_sl < Stability::minimum_slab ) {
		Xdata.S_class2 = 5;
		return true;
	}

	// First, find mean, maximum and minimum hardness
	for (e = Xdata.SoilNode; e < nE; e++) {
		maxH = MAX (maxH, EMS[e].hard);
		minH = MIN (minH, EMS[e].hard);
		mH += EMS[e].hard;
	}

	// Find mean hardness of profile and determine threshold for critical transitions
	mH /= nE;
	thH = MIN (0.5*mH, 1.);

	// Now make the classification for all critical transitions and keep track ....
	for (e = nE-2; e > Xdata.SoilNode; e--) {
		h_Slab += EMS[e].L/cos_sl;
		delta_H = EMS[e+1].hard - EMS[e].hard;
		if ( fabs(delta_H) > thH ) {
			count++;
			if ( delta_H < 0. ) {
				i_weak = e+1;
			} else {
				i_weak = e;
			}

			// Decompose grain type to determine majority shape F1
			qr_TypeToCode(&F1, &F2, &F3, EMS[i_weak].type);

			// Remember that S is initialized to 5!!!
			// First consider wet weak layer
			if ( EMS[i_weak].theta[WATER] > 0.75*lw_SnowResidualWaterContent(EMS[i_weak].theta[ICE]) ) {
				if ( EMS[i_weak].mk%100 < 20 ) {
					S = 1;
				}
			} else { // Then do some stuff for dry snow
				mH_u = 0.;
				for (i = i_weak; i < nE; i++) {
					mH_u += EMS[e].hard;
				}
				mH_u /= (nE - i_weak);
				if ( mH > 2. ) {
					if ( delta_H < 0. ) {
						// Proposal Fz; (see original in version 7.4
						if ( (mH_u < 2.5) || (h_Slab < 0.7) ) {
							if ( (mH_u > 2.) && (h_Slab > 0.5) ) {
								S = MIN(S,4);
							} else {
								S = MIN(S,3);
							}
						} else {
							if ( minH > 2. ) {
								S = MIN(S,4);
							} else {
								S = MIN(S,3);
							}
						}
					}
				} else if ( mH < 1.5 ) {
					if ( (EMS[i_weak].rg > 0.5) && ((F1 > 3) && (F1 < 7)) ) {
						if ( (EMS[i_weak].rg > 0.75) && ((mH_u < 1.5) && (maxH < 2.5)) ) {
							S = 1;
						} else {
							S = MIN (S, 2);
						}
					} else {
						S = MIN (S, 3);
					}
				} else {
					S = MIN (S, 3);
				}
			} // end dry snow
			if ( S < S0 ) {
				S0 = S;
			}
		} // if weak layer found; also end of loop over elements
	} // end for

	if ( count == 0 ) {
		if ( mH > 2. ) {
			S = MIN (S, 4);
		} else if ( mH < 1.5) {
			if ( maxH > 2.3 ) {
				S = MIN (S, 2);
			} else {
				S = 1;
			}
		} else {
			S = MIN (S, 3);
		}
	}

	Xdata.S_class2 = S;

	if ( (S > 5) || (S < 1) ) {
		return false;
	} else {
		return true;
	}
}  // End classifyProfileStability

/**
 * @brief "Pattern recognition" of 10 profile types according to Schweizer, J. and M. Luetschg (2001).
 * "Characteristics of human-triggered avalanches." Cold Reg. Sci. Technol. 33(2-3): 147-162.
 *  Note that analysis is done on vertical snow height.
 * @param *Xdata
 * @param date
 * @return false on error, true otherwise
 */
bool Stability::recognizeProfileType(const mio::Date& date, SN_STATION_DATA& Xdata)
{
	SN_ELEM_DATA *EMS;                              // Element pointer

	int    prf_type=-1;                             // Profile type
	vector<double> z_el;                            // Vertical element heigth (m)
	vector<double> L_el;                            // Vertical element thickness (m)
	vector<double> hard;                            // Hardness in N
	vector<double> red_hard;                        // Reduced hardness in N
	vector<double> deltaN;                          // Difference in hardness between layers in N

	// Temporary variables
	int    e, e_min, e_el, e_max, nE_s;             // Counters
	int    mf_base, weak_base;                      // Booleans
	int    n_window=5;                              // Window half-width in number of elements
	double cos_sl;                                  // cos of slope angle
	double cH;                                      // Vertical snow depth
	double L_base_0=0.2, L_base, L_sum;             // Lengths (m)
	double min_hard=19.472, slope_hard=150.;        // Constants to calculate reduced hardness,
                                                     // (N) and (N m-1), respectively
	double thresh_hard;                             // Hardness threshold (N)
	double mean_hard, mean_red_hard, mean_gsz;      // Means
	double sum_red_hard;
	double red_hard_min, red_hard_max, deltaN_max;             // Extremes
	double z_red_hard_min, z_red_hard_max, z_deltaN_max;       // Extremes' abs. heights
	double pos_min, pos_max, pos_max_deltaN;                   // Extremes' rel. heights

	// cos of slope angle to convert height and thickness to vertical values
	cos_sl = cos(Xdata.SlopeAngle);
	// Vertical snow depth
	cH = (Xdata.cH - Xdata.Ground)/cos_sl;

	// Check for snow profile shallower than 1.5*L_base_0 m (not classifiable)
	if ( cH <= 1.5*L_base_0 ) {
		Xdata.S_class1 = -1;
		return true;
	} else {
		nE_s = Xdata.getNumberOfElements() - Xdata.SoilNode;
	}

	// Dereference element and node pointers
	EMS = &Xdata.Edata[0];
	vector<SN_NODE_DATA>& NDS = Xdata.Ndata;

	try { // Allocate space for temporary vectors
		z_el.resize(nE_s, 0.0);
		L_el.resize(nE_s, 0.0);
		hard.resize(nE_s, 0.0);
		red_hard.resize(nE_s, 0.0);
		deltaN.resize(nE_s, 0.0);
	} catch(exception& e){
		prn_msg(__FILE__, __LINE__, "err", date.getJulianDate(),
			   "Cannot allocate space for temporary objects in Stability::recognizeProfileType");
		throw IOException(e.what(), AT); //this will catch all allocation exceptions
	}

	// Absolute and reduced hardness profiles (N)
	for (e = nE_s-1; e >= 0; e--) {
		z_el[e] = (((NDS[e+Xdata.SoilNode].z + NDS[e+Xdata.SoilNode].u) +
				(NDS[e+Xdata.SoilNode+1].z + NDS[e+Xdata.SoilNode+1].u))/2.)/cos_sl;
		L_el[e] = EMS[e+Xdata.SoilNode].L/cos_sl;
		hard[e] = min_hard*pow(EMS[e+Xdata.SoilNode].hard, 2.3607);
		red_hard[e] = hard[e] - (min_hard + slope_hard*(cH - z_el[e]));
		if ( e == nE_s-1 ) {
			deltaN[e] = fabs(red_hard[e] - min_hard);
		} else {
			deltaN[e] = fabs(red_hard[e] - red_hard[e+1]);
		}
	}

	// Check for base strength (bottom L_base_0 m of snow cover)
	// not considering basal melt-freeze crusts
	L_base = L_base_0;
	L_sum = 0.;
	mean_hard = mean_red_hard = mean_gsz = 0.;
	thresh_hard = 19.472*pow(4., 2.3607);
	mf_base = (hard[0] > thresh_hard);
	e = 0;
	while ( L_sum <= L_base ) {
		if ( mf_base && (hard[e] < thresh_hard) ) {
			mf_base = 0;
			L_base -= L_sum;
			L_sum = 0.;
			mean_hard = mean_red_hard = mean_gsz = 0.;
		}
		L_sum += L_el[e];
		mean_hard += L_el[e]*hard[e];
		mean_red_hard += L_el[e]*red_hard[e];
		mean_gsz += L_el[e]*(2.*EMS[e+Xdata.SoilNode].rg);
		e++;
	}

	// Averages
	mean_hard /= L_sum;
	mean_red_hard /= L_sum;
	mean_gsz /= L_sum;

	// Weak or strong base?
	weak_base = ((mean_hard <= 275.) && (mean_gsz > 0.9));

	// Seek extremes over profile depth
	// Initialise
	e_min = e_el = MIN (e, nE_s - 1);
	e_max = MIN (e_el + n_window, nE_s - 1);
	L_sum = 0.;
	sum_red_hard = 0.;
	red_hard_max = -9999.;
	red_hard_min = 9999.;
	deltaN_max = -9999.;
	z_red_hard_min = z_red_hard_max = z_deltaN_max = 9999.;

	// First evaluation
	for (e = e_min; e <= e_max; e++) {
		L_sum += L_el[e];
		sum_red_hard += L_el[e]*red_hard[e];
	}

	// Use window width of 2*n_window+1 elements
	while ( e_el <= nE_s-1 ) {
		if ( (e_el - e_min) > n_window ) {
			L_sum -= L_el[e_min];
			sum_red_hard -= L_el[e_min]*red_hard[e_min];
			e_min++;
		}
		if ( (e_max < nE_s-1) && ((e_max - e_el) < n_window) ) {
			e_max++;
			L_sum += L_el[e_max];
			sum_red_hard += L_el[e_max]*red_hard[e_max];
		}
		// Find extremes ...
		if ( sum_red_hard/L_sum > red_hard_max ) {
			red_hard_max = sum_red_hard/L_sum;
			z_red_hard_max = z_el[e_el];
		}
		if ( sum_red_hard/L_sum < red_hard_min ) {
			red_hard_min = sum_red_hard/L_sum;
			z_red_hard_min = z_el[e_el];
		}
		e_el++;
	}

	// Find extremes for deltaN (no window required)
	e = 0;
	while ( (z_el[e] < L_base_0) && (e < nE_s-1) ) {
		e++;
	}
	for (; e < nE_s-1; e++) {
		if ( deltaN[e] > deltaN_max ) {
			deltaN_max = deltaN[e];
			z_deltaN_max = z_el[e];
		}
	}

	if ( !((red_hard_max > (-150.*cH)) && (red_hard_min < 1500.)) ) {
		Xdata.S_class1 = -1;
		return false;
	}

	// Classify
	double hard_max; // Max. Hardness at position of maximum

	hard_max = red_hard_max + (min_hard + slope_hard*(cH - z_red_hard_max));
	if ( weak_base ) {
		// Position of extremes
		pos_max = (z_red_hard_max - L_base)/(cH - L_base);
		pos_min = (z_red_hard_min - L_base)/(cH - L_base);
		// Assign weak profile type
		if ( red_hard_max < 50. ) {
			prf_type = 1;
		} else if ( pos_max <= 0.3 ) {
			prf_type = 4;
		} else if ( pos_max <= 0.7 ) {
			prf_type = 3;
		} else if ( pos_max <= 0.9 ) {
			prf_type = 2;
		} else if ( (pos_max) > 0.9 && (hard_max > thresh_hard) ) {
			prf_type = 5;
		} else {
			prf_type = 4;
		}
	} else {// strong base
		// Position of extremes
		pos_max = (z_red_hard_max - L_base)/(cH - L_base);
		pos_min = (z_red_hard_min - L_base)/(cH - L_base);
		pos_max_deltaN = (z_deltaN_max - L_base)/(cH - L_base);

		// Assign strong profile type
		if ( (pos_max_deltaN > 0.85) && (hard_max > thresh_hard) ) {
			prf_type = 9;
		} else if ( (deltaN_max > 150.) && (pos_max_deltaN > pos_min) ) {
			prf_type = 7;
		} else if ( pos_max < 0.3 ) {
			prf_type = 6;
		} else if ( hard_max > thresh_hard ) {
			if ( fabs(red_hard_max - red_hard_min) < 50. ) {
				prf_type = 10;
			} else if ( red_hard_max < 50. ) {
				prf_type = 8;
			} else {
				prf_type =6;
			}
		} else {
			prf_type = 0;
		}
	}
	// end of classify

	Xdata.S_class1 = prf_type;

	return true;
} // End of recognizeProfileType

/*
 *  END OF CLASSIFICATION SECTION
*/

/**
 * @brief On a beautiful morning in September, with Foehn winds outside and incredibly fresh
 * colors, Michael finally started to implement the Stability thing. He is not con-
 * vinced that it is going to work but still dreams of lying in the sun at the Davos
 * lake. The stability information will be based on a very empirical principle. First
 * a distinction is made between "direct action" and "slab" situations. The former
 * have to do with strain weakening during heavy snowfalls or during melt situations.
 * The original Bob intra-layer stability will be adapted for this situation and
 * complemented by the Conway approach. The latter will be handled by an adaptation
 * of the Schweizer - Wiesinger profile classification in combination with a more
 * conventional stability index based on critical shear strength values. Halleluja.
 * @param Mdata SN_MET_DATA
 * @param Xdata SN_STATION_DATA
 */
void Stability::checkStability(const SN_MET_DATA& Mdata, SN_STATION_DATA& Xdata)
{
	int    e, nE, nN;                // Nodal and element counters
	int    Swl_lemon;                // Temporary lemon counter
	double Swl_d, Swl_n, Swl_ssi, zwl_d, zwl_n, zwl_ssi; // Temporary weak layer markers
	double Swl_Sk38, zwl_Sk38;       // Temporary weak layer markers
	double Pk;                       // Penetration depth
	double cos_sl=cos(Xdata.SlopeAngle); // Cosine of slope angle

	StabilityData  STpar;        // Stability parameters
	SN_ELEM_DATA   *EMS;         // Avoids dereferencing the element data pointer

	// Dereference the element pointer containing micro-structure data
	EMS = &Xdata.Edata[0]; nE = Xdata.getNumberOfElements();
	vector<SN_NODE_DATA>& NDS = Xdata.Ndata; nN = Xdata.getNumberOfNodes();

	vector<InstabilityData> SIdata = vector<InstabilityData>(nN); // Parameters for structural instabilities

	// Return if bare soil or PLASTIC
	if ( (nE < Xdata.SoilNode+1) || plastic ) {
		Xdata.S_d = Xdata.S_n = Xdata.S_s = Xdata.S_4 = Stability::max_stability;
		Xdata.z_S_d = Xdata.z_S_n = Xdata.z_S_s = Xdata.z_S_4 = 0.;
		Xdata.S_class1 = Xdata.S_class2 = -1;
		return;
	}

	initStability(Stability::psi_ref, STpar, Xdata, SIdata);
	Pk = st_PenetrationDepth(Xdata);
	for (e = Xdata.SoilNode; e < nE; e++) {
		EMS[e].hard = CALL_MEMBER_FN(*this, mapHandHardness[hardness_model])(EMS[e]);
	}
	EMS[nE-1].s_strength = 100.;
	for (e=nE-2; e >= Xdata.SoilNode; e--) {
		calcReducedStresses(EMS[e+1].C, cos_sl, STpar);
		STpar.strength_up = EMS[e+1].s_strength;

		if ( !(CALL_MEMBER_FN(*this, mapShearStrength[strength_model])(Xdata.cH, cos_sl, Mdata.date,
		                                                               EMS[e], NDS[e+1], STpar))) {
			prn_msg(__FILE__, __LINE__, "msg-", -1., "Node %03d of %03d", e+1, nE+1);
		}
	}

	// Now find the weakest point in the stability profiles for natural and skier indices
	// Initialize
	Swl_lemon = 0;
	Swl_d = Swl_n = Swl_ssi = Swl_Sk38 = INIT_STABILITY;
	zwl_d = zwl_n = zwl_ssi = zwl_Sk38 = Xdata.cH;

	// Natural and "deformation rate" Stability Index
	// Discard Stability::minimum_slab (in m) at surface
	e = nE-1;
	while ( (((Xdata.cH - (NDS[e+1].z + NDS[e+1].u))/cos_sl) < Stability::minimum_slab) && (e > Xdata.SoilNode) ) {
		e--;
	}
	if ( e > Xdata.SoilNode ) {
		// Slab must be thicker than Stability::ground_rough (m)  for an avalanche to release.
		while ( (NDS[e+1].z + NDS[e+1].u)/cos_sl > Stability::ground_rough && (e >= Xdata.SoilNode) ) {
			// "deformation rate" Stability Index: find minimum ...
			EMS[e].S_dr = st_DeformationRateIndex(EMS[e]);
			if ( Swl_d > EMS[e].S_dr ) {
				Swl_d = EMS[e].S_dr;
				zwl_d = (NDS[e].z + NDS[e+1].z + NDS[e].u + NDS[e+1].u)/2.;
			}
			// Natural Stability Index: find minimum ...
			NDS[e+1].S_n = st_NaturalStabilityIndex(STpar);
			if ( Swl_n > NDS[e+1].S_n ) {
				Swl_n = NDS[e+1].S_n;
				zwl_n = NDS[e+1].z + NDS[e+1].u;
			}
			e--;
		}
		// Assign minimum to stability indices
		Xdata.S_d = Swl_d;    Xdata.z_S_d = zwl_d - Xdata.Ground;
		Xdata.S_n = Swl_n;    Xdata.z_S_n = zwl_n - Xdata.Ground;
	} else {
		// Assign bottom values to stability indices
		Xdata.S_d = EMS[Xdata.SoilNode].S_dr;  Xdata.z_S_d = EMS[Xdata.SoilNode].L;
		Xdata.S_n = NDS[Xdata.SoilNode+1].S_n; Xdata.z_S_n = EMS[Xdata.SoilNode].L;
	}

	// Skier Stability Index
	//   Snow depth must be larger than Stability::ground_rough (m) and at least Stability::min_depth_ssi (m)
	//   snow must be left after discarding Pk for a SSI value to be searched.
	if ( (Xdata.cH/cos_sl > Stability::ground_rough) && ((Xdata.cH/cos_sl - Pk) > Stability::min_depth_ssi) ) {
		// Discard penetration depth Pk (in m) at surface
		e = nE-1;
		while ( (((Xdata.cH - (NDS[e+1].z + NDS[e+1].u))/cos_sl) < Pk) && (e > Xdata.SoilNode) ) {
			e--;
		}
		if ( e > Xdata.SoilNode ) {
			// Only down to Pk + Stability::skier_depth (m)
			while ( (((Xdata.cH - (NDS[e+1].z + NDS[e+1].u))/cos_sl) < (Pk + Stability::skier_depth)) && ((NDS[e+1].z + NDS[e+1].u)/cos_sl > Stability::ground_rough) && (e >= Xdata.SoilNode) ) {
				// Skier Stability Index: find minimum OR consider number of structural instabilities in case of near equalities
				const double depth_lay = (Xdata.cH - (NDS[e+1].z + NDS[e+1].u))/cos_sl - Pk;
				NDS[e+1].S_s = st_SkierStabilityIndex(depth_lay, STpar);
				setStructuralStabilityIndex(EMS[e], EMS[e+1], NDS[e+1].S_s, SIdata[e+1]);
				if ( (Swl_ssi > SIdata[e+1].ssi) || ((fabs(Swl_ssi - SIdata[e+1].ssi) < 0.09) && (SIdata[e+1].n_lemon > Swl_lemon)) ) {
					Swl_ssi = SIdata[e+1].ssi;
					zwl_ssi = NDS[e+1].z + NDS[e+1].u ;
					Swl_lemon = SIdata[e+1].n_lemon;
					Swl_Sk38 = NDS[e+1].S_s;
					zwl_Sk38 = NDS[e+1].z + NDS[e+1].u;
				}
				e--;
			}
			// Assign minimum to stability indices
			Xdata.S_s = Swl_Sk38; Xdata.z_S_s = zwl_Sk38 - Xdata.Ground;
			Xdata.S_4 = Swl_ssi;  Xdata.z_S_4 = zwl_ssi - Xdata.Ground;
		} else {
			// Assign bottom values to stability indices
			Xdata.S_s = NDS[Xdata.SoilNode+1].S_s; Xdata.z_S_s = EMS[Xdata.SoilNode].L;
			Xdata.S_4 = SIdata[Xdata.SoilNode+1].ssi; Xdata.z_S_4 = EMS[Xdata.SoilNode].L;
		}
	} else {
		// Assign top values to stability indices
		Xdata.S_s = Stability::max_stability; Xdata.z_S_s = Xdata.cH;
		Xdata.S_4 = SIdata[nN-1].ssi; Xdata.z_S_4 = Xdata.cH;
	}

	switch ( Stability::prof_classi ) {
		case 0:
			// Classify in poor, fair and good based on master thesis of S. Bellaire (September 2005)
			if ( (Swl_ssi > 0.) && (Swl_ssi < 100.) ) {
				if ( Swl_ssi >= 1.55 ) {
					Xdata.S_class2 = 5;
				} else if ( Swl_ssi >= 1.25 ) {
					Xdata.S_class2 = 3;
				} else {
					Xdata.S_class2 = 1;
				}
			} else {
				Xdata.S_class2 = -1;
			}
			break;
		case 1:
		// Classify in poor, fair and good based on re-analysis by Schweizer/Bellaire (paper CRST 46 (2006) 52-59)
			if ( (Swl_ssi > 0.) && (Swl_ssi < 100.) ) {
				if ( Swl_Sk38 >= 0.45 ) {
					Xdata.S_class2 = 5;
				} else if ( Swl_ssi >= 1.32 ) {
						Xdata.S_class2 = 3;
					} else {
						Xdata.S_class2 = 1;
					}
			} else {
				Xdata.S_class2 = -1;
			}
			break;
		case 2:
			// Classify in poor, fair and good based on re-analysis after recalibration of settling,
			// Nov 2007(see rev 250/251)
			if ( (Swl_ssi > 0.) && (Swl_ssi < 100.) ) {
				if ( Swl_lemon >= 2 ) {
					Xdata.S_class2 = 1;
				} else if ( Swl_lemon == 1 ) {
					if ( Swl_Sk38 < 0.48 ) {
						Xdata.S_class2 = 1;
					} else {
						if ( Swl_Sk38 < 0.71 ) {
							Xdata.S_class2 = 3;
						} else {
							Xdata.S_class2 = 5;
						}
					}
				} else {
					Xdata.S_class2 = 3;
				}
			} else {
				Xdata.S_class2 = -1;
			}
			break;
		case 3:
			// Classify in 5 classes based on ideas from Schweizer & Wiesinger
			if ( !classifyProfileStability(Xdata) ) {
				prn_msg( __FILE__, __LINE__, "wrn", Mdata.date.getJulianDate(),
					    "Profile classification failed! (classifyProfileStability)");
			}
			break;
	}

	if ( !ALPINE3D ) {
		// Profile type based on "pattern recognition"; N types out of 10
		// We assume that we don't need it in Alpine3D
		if ( !recognizeProfileType(Mdata.date, Xdata) ) {
			prn_msg( __FILE__, __LINE__, "wrn", Mdata.date.getJulianDate(), "Profile not classifiable! (recognizeProfileType)");
		}
	}
} // End checkStability
