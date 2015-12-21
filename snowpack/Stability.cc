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

#include <snowpack/Stability.h>
#include <snowpack/StabilityAlgorithms.h>
#include <snowpack/Laws_sn.h>

#include <assert.h>

using namespace mio;
using namespace std;

const double Stability::psi_ref = 38.0; ///< Reference slope angle
const double Stability::max_stability = 6.0; ///< Upper stability limit
const double Stability::minimum_slab = 0.1; ///< Minimum slab thickness for natural and deformation stability index (m)
const double Stability::ground_rough = 0.2; ///< The first GROUND_ROUGH m of snow will not be unstable due to ground roughness
const double Stability::min_depth_ssi = 0.1; ///< MIN_DEPTH_SSI m of snow must be left after discarding penetration depth
const double Stability::skier_depth = 1.0; ///< Skiers will not trigger failures SKIER_DEPTH m below penetration depth
const double Stability::min_thick_crust = 0.03; ///< Minimum thickness for a supporting melt-freeze crust (perp to slope, in m)
const size_t Stability::nmax_lemon = 2; ///< Maximum number of structural instabilities looked at ("lemons")
const int Stability::sh_mod = 2; ///< Defines regression model for surface hoar shear strength

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
	mapHandHardness["MONTI"]    = &StabilityAlgorithms::getHandHardnessMONTI;
	mapHandHardness["BELLAIRE"]  = &StabilityAlgorithms::getHandHardnessBELLAIRE;
	mapHandHardness["ASARC"]    = &StabilityAlgorithms::getHandHardnessASARC;

	mapShearStrength["DEFAULT"] = &StabilityAlgorithms::setShearStrengthDEFAULT;
	mapShearStrength["NIED"]    = &StabilityAlgorithms::setShearStrength_NIED;

	return true;
}

/************************************************************
 * non-static section                                       *
 ************************************************************/

Stability::Stability(const SnowpackConfig& cfg, const bool& i_classify_profile)
           : strength_model(), hardness_parameterization(), hoar_density_buried(IOUtils::nodata), plastic(false),
             classify_profile(i_classify_profile)
{
	cfg.getValue("STRENGTH_MODEL", "SnowpackAdvanced", strength_model);
	cfg.getValue("HARDNESS_PARAMETERIZATION", "SnowpackAdvanced", hardness_parameterization);

	const map<string, StabMemFn>::const_iterator it1 = mapHandHardness.find(hardness_parameterization);
	if (it1 == mapHandHardness.end()) throw InvalidArgumentException("Unknown hardness parameterization: "+hardness_parameterization, AT);

	const map<string, StabFnShearStrength>::const_iterator it2 = mapShearStrength.find(strength_model);
	if (it2 == mapShearStrength.end()) throw InvalidArgumentException("Unknown strength model: "+strength_model, AT);

	cfg.getValue("PLASTIC", "SnowpackAdvanced", plastic); //To build a sandwich with a non-snow layer (plastic or wood chips) on top;

	// Density of BURIED surface hoar (kg m-3), default: 125./ Antarctica: 200.
	cfg.getValue("HOAR_DENSITY_BURIED", "SnowpackAdvanced", hoar_density_buried);
}

/**
 * @brief Initializes stability parameters
 * @param Xdata
 */
void Stability::initStability(SnowStation& Xdata)
{
	const size_t nN = Xdata.getNumberOfNodes();
	for(size_t n=Xdata.SoilNode; n<nN; n++) {
		Xdata.Ndata[n].S_n = Stability::max_stability;
		Xdata.Ndata[n].S_s = Stability::max_stability;
		if (n < nN-1) Xdata.Edata[n].S_dr = Stability::max_stability;
	}

	Xdata.S_d = Xdata.S_n = Xdata.S_s = Xdata.S_4 = Xdata.S_5 = Stability::max_stability;
	Xdata.z_S_d = Xdata.z_S_n = Xdata.z_S_s = Xdata.z_S_4 = Xdata.z_S_5 = 0.;
	Xdata.S_class1 = Xdata.S_class2 = -1;
}

/**
 * @brief Returns the structural stability index SSI
 * Adds one lemon to Sk for each structural instability found, presently hardness and grain size differences
 * above a given threshold
 * @param Edata_lower Xdata->Edata[e]
 * @param Edata_upper Xdata->Edata[e+1]
 * @param Sk Skier stability index Sk (Xdata->Ndata[e+1].S_s)
 * @param SIdata [e+1]
 * @return SIdata.ssi [e+1]
 */
double Stability::setStructuralStabilityIndex(const ElementData& Edata_lower, const ElementData& Edata_upper,
                                              const double& Sk, InstabilityData& SIdata)
{
	const double thresh_dhard=1.5, thresh_dgsz=0.5; // Thresholds for structural instabilities
	//const int nmax_lemon = 2; //Maximum number of structural instabilities looked at ("lemons")

	SIdata.n_lemon = 0;
	SIdata.dhard = fabs(Edata_lower.hard - Edata_upper.hard);
	if ( SIdata.dhard > thresh_dhard ) {
		SIdata.n_lemon++;
	}
	SIdata.dgsz = 2.*fabs(Edata_lower.rg - Edata_upper.rg);
	//double ref_gs= MIN (Edata_lower.rg,Edata_upper.rg);
	//SIdata.dgsz = (fabs(Edata_lower.rg - Edata_upper.rg))/(ref_gs);
	if ( SIdata.dgsz > thresh_dgsz ) {
		SIdata.n_lemon++;
	}
	// Skier Stability Index (SSI)
	SIdata.ssi = static_cast<double>(Stability::nmax_lemon - SIdata.n_lemon) + Sk;
	// Limit stability index to range {0.05, Stability::max_stability}
	SIdata.ssi = MAX(0.05, MIN (SIdata.ssi, Stability::max_stability));

	return SIdata.ssi;
}


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
 * @param Mdata CurrentMeteo
 * @param Xdata Profile
 */
void Stability::checkStability(const CurrentMeteo& Mdata, SnowStation& Xdata)
{
	const double cos_sl = Xdata.cos_sl; // Cosine of slope angle
	// Dereference the element pointer containing micro-structure data
	const size_t nN = Xdata.getNumberOfNodes();
	const size_t nE = nN-1;
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;

	vector<InstabilityData> SIdata(nN); // Parameters for structural instabilities
	StabilityData  STpar(Stability::psi_ref*mio::Cst::to_rad);        // Stability parameters

	initStability(Xdata);
	if ( (nE < Xdata.SoilNode+1) || plastic ) return; // Return if bare soil or PLASTIC

	const double Pk = StabilityAlgorithms::compPenetrationDepth(Xdata); // Skier penetration depth
	size_t e = nE; // Counter
	while (e-- > Xdata.SoilNode) {
		EMS[e].hard = (mapHandHardness[hardness_parameterization])(EMS[e], hoar_density_buried);
		EMS[e].S_dr = StabilityAlgorithms::setDeformationRateIndex(EMS[e]);
		StabilityAlgorithms::compReducedStresses(EMS[e].C, cos_sl, STpar);

		if ( !(mapShearStrength[strength_model])(Xdata.cH, cos_sl, Mdata.date,
		                                                               EMS[e], NDS[e+1], STpar)) {
			prn_msg(__FILE__, __LINE__, "msg-", Date(), "Node %03d of %03d", e+1, nN);
		}
		
		const double depth_lay = (Xdata.cH - (NDS[e+1].z + NDS[e+1].u))/cos_sl - Pk; // corrected for skier penetration depth Pk.
		NDS[e+1].S_n = StabilityAlgorithms::getNaturalStability(STpar);
		NDS[e+1].S_s = StabilityAlgorithms::getLayerSkierStability(depth_lay, STpar);
		if (e < nE-1)
			NDS[e+1].ssi = setStructuralStabilityIndex(EMS[e], EMS[e+1], NDS[e+1].S_s, SIdata[e+1]);
		else
			NDS[nN-1].ssi = SIdata[nN-1].ssi = Stability::max_stability;
	}

	// Now find the weakest point in the stability profiles for natural and skier indices
	// Initialize
	size_t Swl_lemon = 0; // Lemon counter
	double Swl_d, Swl_n, Swl_ssi, zwl_d, zwl_n, zwl_ssi; // Temporary weak layer markers
	double Swl_Sk38, zwl_Sk38;       // Temporary weak layer markers
	Swl_d = Swl_n = Swl_ssi = Swl_Sk38 = INIT_STABILITY;
	zwl_d = zwl_n = zwl_ssi = zwl_Sk38 = Xdata.cH;

	// Natural and "deformation rate" Stability Index
	// Discard Stability::minimum_slab (in m) at surface
	e = nE;
	while ((e-- > Xdata.SoilNode) && (((Xdata.cH - (NDS[e+1].z + NDS[e+1].u))/cos_sl) < Stability::minimum_slab)) {};
	if (e==static_cast<size_t>(-1)) e=0; //HACK: this is ugly: e got corrupted if SoilNode==0
	
	if ((e > Xdata.SoilNode) && (e != IOUtils::unodata)) {
		// Slab must be thicker than Stability::ground_rough (m)  for an avalanche to release.
		while ((e-- > Xdata.SoilNode) && ((NDS[e+1].z + NDS[e+1].u)/cos_sl > Stability::ground_rough)) {
			// "deformation rate" Stability Index: find minimum ...
			if (Swl_d > EMS[e].S_dr) {
				Swl_d = EMS[e].S_dr;
				zwl_d = (NDS[e].z + NDS[e+1].z + NDS[e].u + NDS[e+1].u)/2.;
			}
			// Natural Stability Index: find minimum ...
			if ( Swl_n > NDS[e+1].S_n ) {
				Swl_n = NDS[e+1].S_n;
				zwl_n = NDS[e+1].z + NDS[e+1].u;
			}
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
	if ((Xdata.cH/cos_sl > Stability::ground_rough) && ((Xdata.cH/cos_sl - Pk) > Stability::min_depth_ssi)) {
		// Discard penetration depth Pk (in m) at surface
		e = nE;
		while ((e-- > Xdata.SoilNode) && (((Xdata.cH - (NDS[e+1].z + NDS[e+1].u))/cos_sl) < Pk)) {};
		if (e==static_cast<size_t>(-1)) e=0; //HACK: this is ugly: e got corrupted if SoilNode==0
		
		if ((e > Xdata.SoilNode) && (e != IOUtils::unodata)) {
			// Only down to Pk + Stability::skier_depth (m)

			while ((e-- > Xdata.SoilNode) && (((Xdata.cH - (NDS[e+1].z + NDS[e+1].u))/cos_sl) < (Pk + Stability::skier_depth)) && ((NDS[e+1].z + NDS[e+1].u)/cos_sl > Stability::ground_rough)) {
				// Skier Stability Index: find minimum OR consider number of structural instabilities in case of near equalities

				if ( (Swl_ssi > SIdata[e+1].ssi) || ((fabs(Swl_ssi - SIdata[e+1].ssi) < 0.09) && (SIdata[e+1].n_lemon > Swl_lemon)) ) {
					Swl_ssi = SIdata[e+1].ssi;
					zwl_ssi = NDS[e+1].z + NDS[e+1].u ;
					Swl_lemon = SIdata[e+1].n_lemon;
					Swl_Sk38 = NDS[e+1].S_s;
					zwl_Sk38 = NDS[e+1].z + NDS[e+1].u;
				}
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

	switch (Stability::prof_classi) {
		case 0:
			StabilityAlgorithms::classifyStability_Bellaire(Swl_ssi, Xdata);
			break;
		case 1:
			StabilityAlgorithms::classifyStability_SchweizerBellaire(Swl_ssi, Swl_Sk38, Xdata);
			break;
		case 2:
			StabilityAlgorithms::classifyStability_SchweizerBellaire2(Swl_ssi, Swl_lemon, Swl_Sk38, Xdata);
			break;
		case 3:
			// Classify in 5 classes based on ideas from Schweizer & Wiesinger
			if (!StabilityAlgorithms::classifyStability_SchweizerWiesinger(Xdata)) {
				prn_msg( __FILE__, __LINE__, "wrn", Mdata.date,
					    "Profile classification failed! (classifyStability_SchweizerWiesinger)");
			}
			break;
	}

	if (classify_profile) {
		// Profile type based on "pattern recognition"; N types out of 10
		// We assume that we don't need it in Alpine3D
		if (!StabilityAlgorithms::classifyType_SchweizerLuetschg(Xdata)) {
			prn_msg( __FILE__, __LINE__, "wrn", Mdata.date, "Profile not classifiable! (classifyType_SchweizerLuetschg)");
		}
	}
} // End checkStability
