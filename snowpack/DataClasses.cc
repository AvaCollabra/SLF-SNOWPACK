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
 * @file DataClasses.cc
 * @version 11.03
 * @bug     -
 * @brief This module contains the definitions of data classes
 */

#include <snowpack/DataClasses.h>
#include <snowpack/Utils.h>
#include <snowpack/Canopy.h>
#include <snowpack/Metamorphism.h>

using namespace mio;
using namespace std;

/// Number of top elements left untouched by the join functions
const unsigned int SnowStation::number_top_elements = 5;
unsigned int SnowStation::number_of_solutes = 0;

/// Both elements must be smaller than JOIN_THRESH_L (m) for an action to be taken
const double SnowStation::join_thresh_l = 0.015;
/// Volumetric ice content (1), i.e., about 46 kg m-3
const double SnowStation::join_thresh_ice = 0.05;
const double SnowStation::join_thresh_water = 0.01; ///< Water content (1)
const double SnowStation::join_thresh_dd = 0.2;     ///< Dendricity (1)
const double SnowStation::join_thresh_sp = 0.05;    ///< Sphericity (1)
const double SnowStation::join_thresh_rg = 0.125;   ///< Grain radius (mm)

ZwischenData::ZwischenData()
{
	reset();
}

void ZwischenData::reset()
{
	hoar24  = std::vector<double>(48, 0.0);
	drift24 = std::vector<double>(48, 0.0);
	hn3     = std::vector<double>(144, 0.0);
	hn24    = std::vector<double>(144, 0.0);
}

/**
 * @brief Determines the averaged quantities of the current layer with another layer
 * @param Lp1 Thickness (weight) of layer Pdata
 * @param Lp0 Thickness (weight) of current layer
 * @param profile_layer to average with
 */
void SnowProfileLayer::average(const double& Lp0, const double& Lp1, const SnowProfileLayer& profile_layer)
{
	const double layerThickness = Lp0 + Lp1;

	height += Lp1;
	if (Lp1 > Lp0) {
		layerDate = profile_layer.layerDate;
	}
	rho         = (Lp1*profile_layer.rho + Lp0*rho) / layerThickness;
	T           = profile_layer.T;
	gradT       = (Lp1*profile_layer.gradT + Lp0*gradT) / layerThickness;
	strain_rate = (Lp1*profile_layer.strain_rate + Lp0*strain_rate) / layerThickness;
	theta_w     = (Lp1*profile_layer.theta_w + Lp0*theta_w) / layerThickness;
	theta_i     = (Lp1*profile_layer.theta_i + Lp0*theta_i) / layerThickness;
	dendricity  = (Lp1*profile_layer.dendricity + Lp0*dendricity) / layerThickness;
	sphericity  = (Lp1*profile_layer.sphericity + Lp0*sphericity) / layerThickness;
	coordin_num = (Lp1*profile_layer.coordin_num + Lp0*coordin_num) / layerThickness;
	grain_size  = (Lp1*profile_layer.grain_size + Lp0*grain_size) / layerThickness;
	bond_size   = (Lp1*profile_layer.bond_size + Lp0*bond_size) / layerThickness;
	hard        = (Lp1*profile_layer.hard + Lp0*hard) / layerThickness;
	marker      = MAX(profile_layer.marker, marker);
}

SurfaceFluxes::SurfaceFluxes()
  : lw_in(0.), lw_out(0.), lw_net(0.), qs(0.), ql(0.), hoar(0.), qr(0.), qg(0.), qg0(0.), sw_hor(0.),
    sw_in(0.), sw_out(0.), qw(0.), sw_dir(0.), sw_diff(0.), cA(0.), mA(0.), dIntEnergy(0.), drift(0.),
    dhs_corr(0.), cRho_hn(Constants::undefined), mRho_hn(Constants::undefined)
{
	mass.resize(N_MASS_CHANGES);
	load.resize(SnowStation::number_of_solutes);
}

void SurfaceFluxes::reset(const bool& cumsum_mass)
{
	if (cumsum_mass) { // Do not reset cumulated mass balance
		lw_in   = 0.;
		lw_out  = 0.;
		lw_net  = 0.;
		qs      = 0.;
		ql      = 0.;
		qr      = 0.;
		qg      = 0.;
		qg0     = 0.;
		sw_hor  = 0.;
		sw_in   = 0.;
		sw_out  = 0.;
		qw      = 0.;
		sw_dir  = 0.;
		sw_diff = 0.;
		cA      = 0.;
		mA      = 0.;
		dIntEnergy = 0.;
		mass[MS_HNW] = 0.;
		mass[MS_RAIN] = 0.;
	} else {
		*this = SurfaceFluxes(); //reset everything
	}
}

/**
 * @brief Assign surface data from SnowStation to SurfaceFluxes.
 * @param Sdata
 * @param Bdata
 * @param Xdata
 * @param Mdata
 * @param useSoilLayers then the model also deals with soil layers
 * @param soil_flux strength of soil heat flux at greater depth (where it can be considered to be constant)
 */
void SurfaceFluxes::CollectSurfaceFluxes(SurfaceFluxes& Sdata, const BoundCond& Bdata,
                                         SnowStation& Xdata, const CurrentMeteo& Mdata,
                                         const bool& useSoilLayers, const bool& soil_flux)
{
	// 1) Short wave fluxes and Albedo.
	//     Depending on settings (sw_mode) and conditions,
	//     sw_in and sw_out may differ slightly from the original input
	Sdata.sw_in  += Mdata.iswr;
	Sdata.sw_out += Mdata.rswr;
	Sdata.qw     += Mdata.iswr - Mdata.rswr;

	Sdata.cA += Xdata.cAlbedo;
	if (Xdata.mAlbedo != Constants::undefined)
		Sdata.mA += Xdata.mAlbedo;
	else
		Sdata.mA = Constants::undefined;

	// 2) Long wave and heat fluxes.
	Sdata.qs += Bdata.qs;
	//Sdata->ql += Bdata->ql; //HACK needed because latent heat ql not linearized w/ respect to Tss!!!
	Sdata.qr += Bdata.qr;
	Sdata.lw_out += Bdata.lw_out;
	Sdata.lw_net += Bdata.lw_net;
	Sdata.lw_in  += (Bdata.lw_net + Bdata.lw_out);

	// 3) Ground heat fluxes
	if (Xdata.getNumberOfElements() > 0) {
		// 3a) qg0: heat flux at soil/snow boundary
		if (Xdata.getNumberOfElements() > Xdata.SoilNode) {
			if ((Xdata.getNumberOfElements() < 3)
			        && (Xdata.Edata[0].theta[WATER] >= 0.9 * Xdata.Edata[0].res_wat_cont))
				Sdata.qg += 0.;
			else
				Sdata.qg0 += -Xdata.Edata[Xdata.SoilNode].k[TEMPERATURE]
				                 * Xdata.Edata[Xdata.SoilNode].gradT;
		} else if (Xdata.SoilNode > 0) {
			Sdata.qg0 += -Xdata.Edata[Xdata.SoilNode-1].k[TEMPERATURE]
			                 * Xdata.Edata[Xdata.SoilNode-1].gradT;
		} else {
			Sdata.qg0 += Constants::undefined;
		}
		// 3b) qg : geothermal heat flux or heat flux at bottom of snow in case of no soil
		if (useSoilLayers && soil_flux) {
			Sdata.qg += Bdata.qg;
		} else if ((Xdata.getNumberOfElements() < 3)
		               && (Xdata.Edata[0].theta[WATER] >= 0.9 * Xdata.Edata[0].res_wat_cont)) {
			Sdata.qg += 0.;
		} else {
			Sdata.qg += -Xdata.Edata[0].k[TEMPERATURE] * Xdata.Edata[0].gradT;
		}
	} else {
		Sdata.qg0 += Constants::undefined;
		Sdata.qg  += Constants::undefined;
	}

	// 4) Change of internal energy
	if (Xdata.getNumberOfElements() > Xdata.SoilNode)
		Sdata.dIntEnergy += Xdata.dIntEnergy;

	// 5) Compute total masses of snowpack
	Sdata.mass[MS_TOTALMASS] = Sdata.mass[MS_SWE]= Sdata.mass[MS_WATER] = 0.;
	for (size_t e = Xdata.SoilNode; e < Xdata.getNumberOfElements(); e++) {
		Sdata.mass[MS_TOTALMASS] += Xdata.Edata[e].M;
		Sdata.mass[MS_SWE] += Xdata.Edata[e].L * Xdata.Edata[e].Rho;
		Sdata.mass[MS_WATER] += Xdata.Edata[e].L * (Xdata.Edata[e].theta[WATER]
		                            * Constants::density_water);
	}
}

void CanopyData::reset(const bool& cumsum_mass)
{
	if (cumsum_mass) { // Do not reset cumulated mass balance
		// radiation
		rswrac=0.;
		iswrac=0.;
		rswrbc=0.;
		iswrbc=0.;
		ilwrac=0.;
		rlwrac=0.;
		ilwrbc=0.;
		rlwrbc=0.;
		rsnet=0.;
		rlnet=0.;
		// turbulent heat fluxes
		sensible=0.0;
		latent=0.0;
		latentcorr=0.0;
		// auxiliaries
		canopyalb=0.0;
		totalalb=0.0;
		intcapacity=0.0;
	} else {
		initializeSurfaceExchangeData();
	}
}

/**
 * @brief Function called to initialize the canopy "Surface" exchange
 * data (to enable accumulated mass and energy fluxes)
 */
void CanopyData::initializeSurfaceExchangeData()
{
	// radiation
	rswrac=0.; // upward shortwave above canopy
	iswrac=0.; // downward shortwave radiation above canopy
	rswrbc=0.; // upward shortwave below canopy
	iswrbc=0.; // downward shortwave radiation below canopy
	ilwrac=0.; // downward longwave radiation ABOVE canopy
	rlwrac=0.; // upward longwave radiation ABOVE canopy
	ilwrbc=0.; // downward longwave radiation BELOW canopy
	rlwrbc=0.; // upward longwave radiation BELOW canopy
	rsnet=0.;  // net shortwave radiation absorbed by canopy
	rlnet=0.;  // net longwave radiation absorbed by canopy
	// turbulent heat fluxes
	sensible = 0.0;
	latent = 0.0;
	latentcorr = 0.0;
	// mass fluxes
	transp = 0.0;
	intevap = 0.0;
	interception = 0.0;
	throughfall = 0.0;
	snowunload = 0.0;
	// auxiliaries
	canopyalb = 0.0;
	totalalb = 0.0;
	intcapacity = 0.0;
}

ElementData::ElementData() : depositionDate(), L0(0.), L(0.),
                             Te(0.), gradT(0.), Rho(0.), M(0.), sw_abs(0.), rg(0.), dd(0.), sp(0.), rb(0.), ps2rb(0.),
                             N3(0.), mk(0), type(0), metamo(0.), dth_w(0.), Qmf(0.),
                             dE(0.), E(0.), Ee(0.), Ev(0.), EDot(0.), EvDot(0.),
                             S(0.), C(0.), CDot(0.), S_dr(0.), s_strength(0.), hard(0.), dhf(0.)
{
	theta.resize(N_COMPONENTS);
	conc.resize(N_COMPONENTS, SnowStation::number_of_solutes);
	k.resize(N_SN_FIELDS);
	c.resize(N_SN_FIELDS);
	soil.resize(N_SOIL_FIELDS);
}

/**
 * @brief Check volumetric content
 * @version 11.01
 * @return sum of volumetric contents (1)
 */
bool ElementData::checkVolContent()
{
	bool ret = true;
	double sum = 0.;
	for (unsigned int i = 0; i < N_COMPONENTS; i++) {
		sum += theta[i];
	}
	if (sum <= 0.99 || sum >= 1.01) {
		prn_msg(__FILE__, __LINE__, "wrn", Date(), "SUM of volumetric contents = %1.4f", sum);
		ret = false;
	}
	if ((theta[SOIL] < -Constants::eps) || (theta[ICE] < -Constants::eps)
	        || (theta[WATER] < -Constants::eps) || (theta[AIR] < -Constants::eps) ) {
		prn_msg(__FILE__, __LINE__, "wrn", Date(), "Negative volumetric content!");
		ret = false;
	}
	return ret;
}

/**
 * @brief Computes heat capacity of an element based on volumetric contents
 * @version 11.01
 * @return Effective heat capacity (J kg-1 K-1)
 */
void ElementData::heatCapacity()
{
	double c_p;

	c_p  = Constants::density_air * theta[AIR] * Constants::specific_heat_air;
	c_p += Constants::density_ice * theta[ICE] * Constants::specific_heat_ice;
	c_p += Constants::density_water * theta[WATER] * Constants::specific_heat_water;
	c_p += soil[SOIL_RHO] * theta[SOIL] * soil[SOIL_C];
	c_p /= Rho;
	c[TEMPERATURE] = c_p;
}

/**
 * @brief Computes cold content of an element
 * @version 10.08
 * @return Cold content (J m-2)
 */
double ElementData::coldContent()
{
	return (Rho * c[TEMPERATURE] * (Te - Constants::melting_tk) * L);
}

/**
 * @brief Density dependent extinction coefficient -> Michi's magic trick... out of his magic hat
 * @version 9Y.mm
 * @return Density dependent extinction coefficient (m-1)
 */
double ElementData::extinction()
{
	return(Rho/10. + 30.);
	//return(Edata->Rho/10. + 30.);
	//return(Edata->Rho/7.  + 70.);
	//return(Edata->Rho/7.  + 75. - 0.0*Edata->theta[WATER]);
}

/**
 * @brief Estimate the residual water content RWC by Vol \n
 * From work by Coleou and Lesaffre, 1998, Ann. Glaciol., 26, 64-68. \n
 * Experimental range:
 * - density unsoaked: 235 to 580
 * - density soaked: 328 to 589 kg m-3
 * - RWC by Mass 0.049 to 0.029
 * @note That function will limit range to 0.0264 to 0.08 RWC by Vol
 * @version 11.01
 * @return residual water content of snow element (1)
 */
void ElementData::snowResidualWaterContent()
{
	res_wat_cont = snowResidualWaterContent(theta[ICE]);
}

double ElementData::snowResidualWaterContent(const double theta_i)
{
	double res_wat_cont;

	if (theta_i > 0.23) {
		res_wat_cont = 0.0264 + 0.0099 * (1. - theta_i) / theta_i;
	} else {
		res_wat_cont = 0.08 - 0.1023 * (theta_i - 0.03);
	}
	return MIN(res_wat_cont, 0.08); //NOTE: MIN() only needed in case of theta_i < 0.03
}

/**
 * @brief Field Capacity Soil is dependent on grain properties.
 * This nice formulation is based on some tedious curve fitting by
 * Martina Luetschg. The data stems from Scheffer und Schachtschabel
 * (Richtwerte Baugrund), which once more proves that "nomen est omen".
 * If my name was "Schachtschabel", I would never ever be dealing with
 * wet soils and Baugrund.
 * @author Michael Lehning
 * @version 9Y.mm
 * @return Soil field capacity (?)
 */
double ElementData::soilFieldCapacity()
{
	double fc;
	if (!(rg > 0.)) {
		fc = MIN(SnLaws::field_capacity_soil, (1. - theta[SOIL]) * 0.1);
	} else {
		//Follow implementation by Tobias Hipp master thesis.
		//Note that the value of 0.0976114 is more precise and the value of 60.8057 is
		//slightly different from what is mentioned in thesis, to make the function continuous over rg.
		if(rg<17.0) {
			fc = MIN(0.95, 0.32 / sqrt(rg) + 0.02);
		} else {
			if(rg<60.8057) {
				fc=0.0976114-0.002*(rg-17.0);
			} else {
				fc=0.01;
			}
		}
	}
	return fc;
}

/**
 * @brief SNOW ELASTICITY  :  This important routine was programmed by Marc Christen, who took it directly
 * from Mellor's famous 1975 paper on SNOW MECHANICS in the GRINDLEWALD symposium. Dimensions
 * are in [Pa]. (Presently, it is NOT temperature dependent.)
 * @version 9Y.mm
 * @return Module of elasticity (Pa)
 */
double ElementData::snowElasticity()
{
	double g, h;

	if (Rho > 1000.)
		return Constants::big;

	if (Rho >= 70.) {
		g = ((Rho / 1000.0) * 8.235) - 0.47;
	} else {
		g = ((70. / 1000.0) * 8.235 ) - 0.47;
	}
	h = pow(10.0, g);
	return (h * 100000.0);
}

/**
 * @brief Computes the enhancement of hydrostatically applied stress (overburden) in the necks (or bonds)
 * @version 11.01
 * @return Enhancement factor for neck stress (1))
 */
double ElementData::neckStressEnhancement()
{
	const double stressEnhance = (4. / (N3 * theta[ICE])) * (rg * rg) / (rb * rb);
	return stressEnhance;
}

/**
 * @brief A non-generic function to compute the concave neck radius (mm). \n
 * It is assumed that the neck is bound by a sphere fitting between the two grains side by side
 * such as the sphere goes to RB from the axis (this is a quick and dirty approximation)
 * @author Mathias Bavay
 * @version 9.mm
 * @return Concave neck radius (mm)
 */
double ElementData::concaveNeckRadius()
{
	if ( (rg - rb) < Constants::eps ) {
		prn_msg(__FILE__, __LINE__, "wrn", Date(), "Infinite radius of curvature, rg(%lf) = rb(%lf); return Constants::big!", rg, rb);
		return (Constants::big);
	} else {
		return (rb*rb / (2. * (rg - rb)));
	}
}

/**
 * @brief Computes the neck (or bond) length (mm)
 * @version 11.01
 * @return Neck length (mm)
 */
double ElementData::neckLength()
{
	double rc = concaveNeckRadius();
	return ((2. * rg * rc) / (rg + rc));
}

/**
 * @brief Relates the neck strain to the global volumetric strain
 * @version 11.01
 * @return Macro factor (1)
 */
double ElementData::neck2VolumetricStrain()
{
	const double Ln = neckLength();
	return (Ln / (2. * rg + Ln));
}
/**
 * @brief Determine the type of snow \n
 * First revisited by Fierz and Bellaire 2006 and 2007
 * TODO needs to be adapted to international classification
 * @version 11.11
 * @return snow type code according to old-fashioned Swiss tradition
 */

void ElementData::snowType()
{
	type = snowType(dd, sp, 2.*rg, mk%100, theta[WATER], res_wat_cont);
}

int ElementData::snowType(const double dendricity, const double sphericity,
                          const double grain_size, const int marker, const double theta_w, const double res_wat_cont)
{
	int a=-1,b=-1,c=0;
	int sw2;

	// Dry snow
	if (dendricity > 0.) { // Dry dendritic (new) snow: dendricity and sphericityhericity determine the class
		sw2 = (int)(sphericity*10.);
		if (dendricity > 0.80 ) { // ori 0.90, 27 Nov 2007 sb
			a = 1; b = 1; c = 0;
		} else if (dendricity > 0.70) { // ori 0.85, 27 Nov 2007 sb
			a = 1; b = 2; c = 1;
		} else if (dendricity > 0.65) { // ori 0.75, 27 Nov 2007 sb
			a = 2; b = 1; c = 0;
		} else if (dendricity > 0.60) { // ori 0.70, 27 Nov 2007 sb
			a = 2; b = 1; c = 1;
		} else if (dendricity > 0.30) {
			a = 2; b = 2; c = 0;
		} else if (dendricity > 0.05) {
			a = 2;
			switch (sw2) {
				case 0: case 1: case 2:
					b = 4; c = 0; break;
				case 3: case 4:
					b = 4; c = 1; break;
				case 5: case 6:
					b = 3; c = 1; break;
				default:
					b = 3; c = 0;
			}
		} else {
			switch (sw2) {
				case 0: case 1:
					a = 4; b = 4; c = 0; break;
				case 2: case 3: case 4:
					a = 4; b = 2; c = 1; break;
				case 5: case 6: case 7:
					a = 3;  b = 2; c = 1; break;
				default:
					a = 3; b = 3; c = 0;
			}
		}
	} else if (marker <= 2) {
		// Dry non-dendritic snow
		// Sphericity is most important for "a", while the marker is most important for "b","c"
		if (grain_size < 0.7) {
			sw2 = (int)(sphericity*10.);
			switch (sw2) {
				case 0: case 1:
					a = 4; b = 4; c = 0; break;
				case 2: case 3:
					a = 4; b = 3; c = 1; break;
				case 4: case 5:
					a = 4;  b = 3; c = 0; break;
				case 6: case 7:
					a = 3;  b = 4; c = 1; break;
				default:
					a = 3; b = 3; c = 0;
			}
		} else if (grain_size < 1.1) {
			if (sphericity < 0.2) {
				a = 4; b = 4; c = 0;
			} else if (sphericity < 0.4) {
				a = 4; b = 9; c = 0;
			} else { // sphericityhericity limited to sphericity_max=0.5 in Metamorphism.c
				a = 9; b = 9 ; c = 0;
			}
		} else if (grain_size < 1.5) {
			if (sphericity < 0.2) {
				a = 4; b = 5; c = 0;
			} else if (sphericity < 0.4) {
				a = 4; b = 9; c = 1;
			} else { // sphericityhericity limited to sphericity_max=0.5 in Metamorphism.c
				a = 9; b = 9 ; c = 0;
			}
		} else {
			if (sphericity < 0.2) {
				a = 5; b = 5; c = 0;
			} else if (sphericity < 0.4) {
				a = 5; b = 9; c = 1;
			} else { // sphericityhericity limited to sphericity_max=0.5 in Metamorphism.c
				a = 9; b = 5 ; c = 1;
			}
		}
	}
	// Melt forms
	if (marker >= 10) {
		if (dendricity > 0.) { // Wet dendritic snow
			if (sphericity > 0.7) {
				b = a; a = 7; c = 0;
			} else {
				b = 7 ; c = 1;
			}
		} else { // Wet non-dendritic snow
			b = 7; c = 0;
			if (sphericity > 0.75) {
				a = 7;
			} else if( sphericity > 0.4 ) {
				if (grain_size <= 0.7) {
					b = a; a = 7;
				} else if (marker != 13 ) {
					if (grain_size <= 1.5) {
						a = 7; b = 9; c = 1;
					} else {
						a = 7; b = 5; c = 1;
					}
				} else {
					a = 7; b = 6;  c = 1;
				}
			} else {
				if (grain_size <= 1.5) {
					a = 4;
				} else {
					a = 5;
				}
				if (sphericity <= 0.2) {
					c = 1;
				}
			}
		}
	}
	// Now treat a couple of exceptions - note that the order is important
	if (b < 0) b = a;
	if ((marker >= 20) && (theta_w < 0.1*res_wat_cont)) { // MFcr Melt-Freeze
		c = 2;
	}
	switch (marker) {
	case 3: // SH   Surface Hoar
		a = 6; b = 6; c = 0;
		break;
	case 4: // PPgp Graupel
		a = 0; b = 0; c = 0;
		break;
	case 7: case 8: // Glacier ice & IFil, that is, ice layers within the snowpack
		a = 8; b = 8; c = 0;
		break;
	}

	return (a*100 + b*10 + c);
}

std::ostream& operator<<(std::ostream& os, const ElementData& data)
{
	os << std::fixed << std::showpoint;
	os << data.depositionDate << " " << setprecision(4) << data.L << " - " << data.type << " (" << setprecision(2) << data.rg << "/" << data.rb << ") - " << data.mk;
	os << setprecision(2) << "(" << data.theta[SOIL] << "-" << data.theta[ICE] << "-" << data.theta[WATER] << "-" << data.theta[AIR] << ")";
	os << std::endl;

	return os;
}

SnowStation::SnowStation(const bool& i_useCanopyModel, const bool& i_useSoilLayers) :
	meta(), cAlbedo(0.), mAlbedo(0.), SoilAlb(0.), BareSoil_z0(0.), SoilNode(0), cH(0.),
	mH(0.), Ground(0.), hn(0.), rho_hn(0.), windward(false), ErosionLevel(0), ErosionMass(0.),
	S_class1(0), S_class2(0), S_d(0.), z_S_d(0.), S_n(0.), z_S_n(0.), S_s(0.), z_S_s(0.), S_4(0.),
	z_S_4(0.), S_5(0.), z_S_5(0.), Kt(NULL), ColdContent(0.), dIntEnergy(0.),
	SubSurfaceMelt('x'), SubSurfaceFrze('x'), Cdata(), tag_low(0),
	useCanopyModel(i_useCanopyModel), useSoilLayers(i_useSoilLayers),
	nNodes(0), nElems(0)
{
	Edata = vector<ElementData>();
	Ndata = vector<NodeData>();
}

SnowStation::~SnowStation()
{
	MYTYPE* pMat = (MYTYPE*) Kt;

	if (pMat != NULL) {
		if ( pMat->State == ConMatrix ){
			ReleaseConMatrix(&pMat->Mat.Con);
		} else if ( pMat->State == BlockMatrix  ){
			ReleaseBlockMatrix(&pMat->Mat.Block);
		}
		//free(pMat);
	}
}

/**
 * @brief Computes the internal energy change of the snowpack during one computation time step (J m-2)
 * @version 11.01
 */
void SnowStation::compSnowpackInternalEnergyChange(const double sn_dt)
{
	size_t e = SoilNode;
	double cold_content_in = ColdContent;
	double melt_freeze_energy = 0.;

	if (nElems > SoilNode) {
		ColdContent = 0.;
		for (; e<nElems; e++) {
			melt_freeze_energy -= Edata[e].Qmf * Edata[e].L * sn_dt;
			ColdContent += Edata[e].coldContent();
		}
		dIntEnergy = ColdContent - cold_content_in + melt_freeze_energy;
	} else {
		dIntEnergy = 0.;
	}
}

/**
 * @brief Returns modelled internal snow or/and soil temperature (instantaneous value; degC),
 *        at a given position z perpendicular to slope (m) \n
 *        z must be less than computed height (Xdata->cH), otherwise modeled temperature is set to Constants::undefined
 * @version 11.03
 * @param z Sensor position perpendicular to slope (m)
 */
double SnowStation::getModelledTemperature(const double& z) const
{
	int n_up;           // Upper node number
	double z_up, z_low; // Upper and lower nodes around position z of sensor

	if ( (z == Constants::nodata) || !((getNumberOfNodes() > 1) && (z < cH)) ) {
		return Constants::nodata;
	} else {
		n_up = findUpperNode(z, Ndata, getNumberOfNodes());
		z_low = (Ndata[n_up-1].z + Ndata[n_up-1].u);
		z_up = (Ndata[n_up].z + Ndata[n_up].u);
		return (K_TO_C(Ndata[n_up-1].T + (z - z_low)*(Ndata[n_up].T-Ndata[n_up-1].T)/(z_up-z_low)));
	}
}

/**
 * @brief Reallocate element and node data \n
 * Xdata->Edata, Xdata->Ndata and Xdata->nElems, Xdata->nNodes are reallocated or reset, respectively.
 * In case of augmenting the element number, the new elements are initialized to 0 (memset)
 * @param number_of_elements The new number of elements
 */
void SnowStation::resize(const unsigned int& number_of_elements)
{
	try {
		Edata.resize(number_of_elements, ElementData());
		Ndata.resize(number_of_elements + 1);
	}catch(const exception& e){
		throw IOException(e.what(), AT); //this will catch all allocation exceptions
	}

	nElems = (int)Edata.size();
	nNodes = (int)Ndata.size();
}

size_t SnowStation::getNumberOfElements() const
{
	return nElems;
}

size_t SnowStation::getNumberOfNodes() const
{
	return nNodes;
}

/**
 * @brief If more than NUMBER_TOP_ELEMENTS snow elements exist, attempt to reduce their number in the FEM mesh,
 * leaving NUMBER_TOP_ELEMENTS surface elements untouched \n
 * Pairs of elements within the snow cover satisfying the join conditions of joinCondition() are merged
 * by placing everything in the lower element, setting the density of upper element to Constants::undefined,
 * and getting rid of node in between. \n
 * The elements being very similar and thus the microstructure parameters being approximately equal
 * as defined in joinCondition(), simply average the microstructure properties \n
 * NOTE that the condense element check is placed at the end of a time step, allowing elements do develop on their own.
 * @param number_top_elements The number of surface elements to be left untouched
 */
void SnowStation::joinElements(const unsigned int& number_top_elements)
{
	size_t e0, e1;  // Lower (e0) and upper (e1) element index
	size_t rnE;     // Reduced number of elements
	size_t nJoin=0; // Number of elements to be removed

	if (nElems - SoilNode < number_top_elements+1) {
		return;
	}
	for (e0 = SoilNode, e1 = SoilNode+1; e0 < nElems-number_top_elements; e0++, e1++) {
	  if (joinCondition(Edata[e0], Edata[e1])) {
			mergeElements(Edata[e0], Edata[e1], true);
			nJoin++;
			Edata[e1].Rho = Constants::undefined;
			e0++; e1++;
		}
	}
	if (nJoin > 0) {
		rnE = nElems - nJoin;
		reduceNumberOfElements(rnE);
	}
}

/**
 * @brief Remove the upper "marked" element of two (snow only) \n
 * -# Joining two elements:
 *  - density is undefined
 *  - take the uppermost node of both
 * -# Removing melted or thin elements
 *  - density is undefined AND length negative (*= -1.) as the latter will be used!
 *  - keep upper node of lowest element
 * @param rnE Reduced number of elements
 */
void SnowStation::reduceNumberOfElements(const unsigned int& rnE)
{
	unsigned int e0;                    // Lower element index
	unsigned int eNew;                  // New element index
	double cH_old, dL=0.;

	for (e0 = SoilNode, eNew = SoilNode; e0 < nElems; e0++) {
		if (Edata[e0].Rho == Constants::undefined) {
			if (Edata[e0].L > 0.0) { // Joining elements
				Ndata[eNew] = Ndata[e0+1];
				Ndata[eNew].z = Ndata[e0+1].z + Ndata[e0+1].u + dL;
				Ndata[eNew].u = Ndata[e0].udot = 0.;
				Ndata[eNew].ssi = Ndata[e0+1].ssi;
				Ndata[eNew].S_s = Ndata[e0+1].S_s;
				Ndata[eNew].S_n = Ndata[e0+1].S_n;
			} else { // Removing elements for negative length L
				dL += Edata[e0].L;
			}
		} else {
			if (eNew < e0) {
				Edata[eNew] = Edata[e0];
				Ndata[eNew+1] = Ndata[e0+1];
			}
			Ndata[eNew+1].z = Ndata[e0+1].z + Ndata[e0+1].u + dL;
			Ndata[eNew+1].u = Ndata[e0+1].udot = 0.;
			Ndata[eNew+1].ssi = Ndata[e0+1].ssi;
			Ndata[eNew+1].S_s = Ndata[e0+1].S_s;
			Ndata[eNew+1].S_n = Ndata[e0+1].S_n;
			eNew++;
		}
	}

	resize(rnE);

	cH_old = cH;
	cH = Ndata[nNodes-1].z + Ndata[nNodes-1].u;
	mH -= (cH_old - cH);
	ErosionLevel = MAX(SoilNode, MIN(ErosionLevel, rnE-1));
}

/**
 * @brief This routine initializes the snow cover structure which contains all information about a station
 * including element, nodal and canopy data \n
 * Because you are working with layers, the first node is a special case; the rest can be generated in a loop ....
 * The bottom temperature at the beginning of the computation is given by the temperature at the top of the
 * lowest soil or snow layer \n
 * IMPORTANT: it is very important for Alpine3D that Cdata.height is initialized even if CANOPY = 0,
 * otherwise SnowInterface will not recognize the canopy grids (David 2007-06-25).
 * @version 10.02
 * @param SSdata
 * @param i_sector defines the exposition sector of the slope (width 360./number_slopes)
 */
void SnowStation::initialize(const SN_SNOWSOIL_DATA& SSdata, const unsigned int i_sector)
{
	size_t ll, le, e, n;      //  Counters for layers, layer elements, elements, and nodes
	int real_soil_no_sandwich = 1;  // Switch to count real soil layers

	cAlbedo = SSdata.Albedo;
	SoilAlb = SSdata.SoilAlb;
	BareSoil_z0 = SSdata.BareSoil_z0;

	meta = SSdata.meta;
	sector = i_sector;

	mH = cH = SSdata.Height;

	nNodes = SSdata.nN;
	nElems = SSdata.nN-1;
	resize(nElems);

	SoilNode = 0;
	Ground = 0.0;
	Ndata[0].z = 0.;
	if (SSdata.nLayers > 0) {
		Ndata[0].T = SSdata.Ldata[0].tl;
	} else {
		Ndata[0].T = C_TO_K(0.0);
	}
	Ndata[0].u = 0.;
	Ndata[0].f = 0.;
	Ndata[0].udot = 0.;
	Ndata[0].hoar = 0.;
	Ndata[0].S_n=6.;   // Interface static natural stability index
	Ndata[0].S_s=6.;   // Interface stability index Sk38 (skier)
	for (ll = 0, n = 1; ll < SSdata.nLayers; ll++) {
		double dT;
		// Update ground heigth and SoilNode number
		if (SSdata.Ldata[ll].phiSoil > 0.0 && real_soil_no_sandwich) {
			Ground += SSdata.Ldata[ll].hl;
			SoilNode += SSdata.Ldata[ll].ne;
		} else {
			real_soil_no_sandwich = 0;
		}
		if (ll == 0) {
			dT = 0.;
		} else {
			dT = (SSdata.Ldata[ll].tl - SSdata.Ldata[ll-1].tl) / (SSdata.Ldata[ll].ne);
		}
		for (le = 0; le < SSdata.Ldata[ll].ne; le++, n++ ) {
			Ndata[n].z = Ndata[n-1].z + SSdata.Ldata[ll].hl / SSdata.Ldata[ll].ne;
			Ndata[n].T = Ndata[n-1].T + dT;
			Ndata[n].u = 0.;
			Ndata[n].f = 0.;
			Ndata[n].udot = 0.;
			Ndata[n].S_n = INIT_STABILITY;   // Static natural stability index
			Ndata[n].S_s = INIT_STABILITY;   // Alternative Stability Index (skier stability)
		}
	}

	if (SoilNode == 0 && useSoilLayers) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "SNP_SOIL set but no soil layers given");
		throw IOException("Snowpack Initialization failed", AT);
	}

	// INITIALIZE THE ELEMENT DATA
	for (ll = 0, e = 0; ll<SSdata.nLayers; ll++) {
		for (le = 0; le < SSdata.Ldata[ll].ne; le++, e++) {
			// Element's JulianQ Date
			Edata[e].depositionDate = SSdata.Ldata[ll].layerDate;
			// Temperature data
			Edata[e].Te = (Ndata[e].T+Ndata[e+1].T) / 2.;
			Edata[e].L0 = Edata[e].L = (Ndata[e+1].z - Ndata[e].z);
			Edata[e].gradT = (Ndata[e+1].T-Ndata[e].T) / Edata[e].L;
			// Creep data
			Edata[e].E = Edata[e].S = Edata[e].EDot=0.0;
			Edata[e].Ev = Edata[e].Ee = Edata[e].EvDot=0.0;
			// Very important to initialize the increments in length and strain
			Edata[e].dE = 0.0;
			// Volumetric Components
			Edata[e].theta[SOIL]  = SSdata.Ldata[ll].phiSoil;
			Edata[e].theta[AIR]   = SSdata.Ldata[ll].phiVoids;
			Edata[e].theta[ICE]   = SSdata.Ldata[ll].phiIce;
			Edata[e].theta[WATER] = SSdata.Ldata[ll].phiWater;
			Edata[e].soil[SOIL_RHO] = SSdata.Ldata[ll].SoilRho;
			Edata[e].soil[SOIL_K]   = SSdata.Ldata[ll].SoilK;
			Edata[e].soil[SOIL_C]   = SSdata.Ldata[ll].SoilC;
			for (int ii = 0; ii < signed(SnowStation::number_of_solutes); ii++) {
				Edata[e].conc[SOIL][ii]  = SSdata.Ldata[ll].cSoil[ii];
				Edata[e].conc[ICE][ii]  = SSdata.Ldata[ll].cIce[ii];
				Edata[e].conc[WATER][ii] = SSdata.Ldata[ll].cWater[ii];
				Edata[e].conc[AIR][ii]  = SSdata.Ldata[ll].cVoids[ii];
			}
			Edata[e].Rho = Edata[e].theta[ICE]*Constants::density_ice +
				Edata[e].theta[WATER]*Constants::density_water + Edata[e].theta[SOIL]*Edata[e].soil[SOIL_RHO];
			// conductivities, specific heat and moisture content
			Edata[e].k[TEMPERATURE] = Edata[e].k[SEEPAGE] = Edata[e].k[SETTLEMENT] = 0.;
			Edata[e].heatCapacity();
			Edata[e].c[SEEPAGE] = Edata[e].c[SETTLEMENT] = 0.;
			Edata[e].snowResidualWaterContent();
			// Set the initial short wave radiation to zero
			Edata[e].sw_abs = 0.;
			// Phase change variables
			Edata[e].Qmf = 0.;
			Edata[e].dth_w = 0.;
			// Micro-structure data
			Edata[e].dd = SSdata.Ldata[ll].dd;
			Edata[e].sp = SSdata.Ldata[ll].sp;
			Edata[e].rg = SSdata.Ldata[ll].rg;
			Edata[e].rb = SSdata.Ldata[ll].rb;
			Edata[e].N3 = Metamorphism::getCoordinationNumberN3(Edata[e].Rho);
			Edata[e].mk = SSdata.Ldata[ll].mk;
			Edata[e].snowType();
			Ndata[e+1].hoar = SSdata.Ldata[ll].hr;
			// Memories, memories
			Edata[e].CDot = SSdata.Ldata[ll].CDot;
			Edata[e].metamo = SSdata.Ldata[ll].metamo;
			Edata[e].S_dr = INIT_STABILITY;
			Edata[e].hard = 0.;
			Edata[e].M = Edata[e].Rho * Edata[e].L0;
		} // end of element layer for
	} // end of layer for
	if (SSdata.ErosionLevel > 0) {
		ErosionLevel = SSdata.ErosionLevel;
	} else {
		ErosionLevel = MAX(SoilNode, nElems-1);
	}
	// Find the real Cauchy stresses
	e = nElems; double SigC = 0.0;
	while (e-- > 0) {
		if (e < nElems-1) SigC -= (Edata[e+1].M / 2.) * Constants::g * cos(DEG_TO_RAD(meta.getSlopeAngle()));
		SigC -= (Edata[e].M / 2.) * Constants::g * cos(DEG_TO_RAD(meta.getSlopeAngle()));
		Edata[e].C = SigC;
	}
	// Cold content
	compSnowpackInternalEnergyChange(900.); // Time (900 s) will not matter if Qmf == 0. for all layers

	// INITIALIZE CANOPY DATA
	if (useCanopyModel) {
		Cdata.height=SSdata.Canopy_Height;
		Cdata.storage=0.0;           // intercepted water (kg m-2 or mm Water Equivalent)
		Cdata.temp=273.15;	          // temperature (K)
		Cdata.canopyalb=Canopy::can_alb_dry; // albedo [-], which is a function of the dry canopy albedo and intercepted snow
		Cdata.wetfraction=0.0;
		Cdata.lai=SSdata.Canopy_LAI;

		Cdata.sigf=1.-exp(-Canopy::krnt_lai * (Cdata.lai)); // radiation transmissivity (-)
		Cdata.ec=1.0;               //longwave emissivity

		Cdata.z0m=Cdata.height*0.1;
		Cdata.z0h=Cdata.z0m*0.1;
		Cdata.zdispl=Cdata.height*0.66;
		Cdata.direct_throughfall=SSdata.Canopy_Direct_Throughfall;
		if (!(SSdata.Canopy_Direct_Throughfall >= 0. && SSdata.Canopy_Direct_Throughfall <= 1.)) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "Given Canopy Throughfall (*.sno file) = %lf but Canopy is set",
			        SSdata.Canopy_Direct_Throughfall);
			throw IOException("Snowpack Initialization failed", AT);
		}
		Cdata.ra=0.0;
		Cdata.rc=0.0;
		Cdata.rs=0.0;
		Cdata.rstransp=0.0;
	} else {
		Cdata.height=0.0;
		Cdata.storage=0.0;           // intercepted water (kg m-2 or mm Water Equivalent)
		Cdata.temp=273.15;	          // temperature (K)
		Cdata.canopyalb=Canopy::can_alb_dry; // albedo [-], which is a function of the dry canopy albedo and intercepted snow
		Cdata.wetfraction=0.0;
		Cdata.intcapacity=0.0;
		Cdata.lai=0.0;
		Cdata.sigf=1.0;              // radiation transmissivity (-)
		Cdata.ec=1.0;                //longwave emissivity

		Cdata.z0m=0.0;
		Cdata.z0h=0.0;
		Cdata.zdispl=0.0;
		Cdata.direct_throughfall=1.0;
		Cdata.ra=0.0;
		Cdata.rc=0.0;
		Cdata.rs=0.0;
		Cdata.rstransp=0.0;
	}
}

/**
 * @brief Boolean routine to check whether two snow elements can be joined
 * -# NO ACTION will be taken if one of the two elements is
 * 	- a soil element
 * 	- larger than join_thresh_l
 * 	- tagged
 * 	- dry surface hoar (mk=3),
 * 	- dendritic but not both
 * -# OTHERWISE we use criteria for dendricity, sphericity, volumetric ice or water content, grain size and marker
 * - NOTE July 2006: whatever type of thin elements are treated in wt_ElementRemoval()
 * @param Edata0 Lower element
 * @param Edata1 Upper element
 * @return true if the two elements should be joined, false otherwise
 */
bool SnowStation::joinCondition(const ElementData& Edata0, const ElementData& Edata1)
{
	if ( (Edata0.L > join_thresh_l) || (Edata1.L > join_thresh_l) )
		return false;

	if ( Edata0.mk%100 != Edata1.mk%100 )
		return false;

	if ( fabs(Edata0.sp - Edata1.sp) > join_thresh_sp )
		return false;

	if ( Edata0.theta[SOIL] > 0. || Edata1.theta[SOIL] > 0. )
		return false;

	if ( (Edata0.mk >= 100) || (Edata1.mk >= 100) )
		return false;

	if ( (Edata0.mk%100 == 3) || (Edata1.mk%100 == 3) )
		return false;

	if ( (Edata0.dd > join_thresh_dd || Edata1.dd > join_thresh_dd) &&
		!(Edata0.dd > join_thresh_dd && Edata1.dd > join_thresh_dd) ) {
		return false;
	} else if ( fabs(Edata0.dd - Edata1.dd) > join_thresh_dd ) {
		return false;
	}

	if ( fabs(Edata0.theta[ICE] - Edata1.theta[ICE]) > join_thresh_ice )
		return false;

	if ( fabs(Edata0.theta[WATER] - Edata1.theta[WATER]) > join_thresh_water )
		return false;

	if ( fabs(Edata0.rg - Edata1.rg) > join_thresh_rg )
		return false;

	return true;
}

/**
 * @brief Merging two elements
 * - Joining:
 * 	- Keep the lower element, that is, the lowest snow element always survives!
 * 	- The properties of the upper and lower elements are (depth) averaged in an appropriate way.
 * 	- The new length is the sum of both
 * 	- Keep the birthday of the upper element
 * 	- Keep the tag of the upper element only if the lower is untagged (mk >= 100)
 * 	- @note Joining two elements may cause the tag (marker >= 100) to "jump" abruptly
 * - Removing:
 * 	- Remaining ice, liquid water, solutes, etc. are added to the lower element
 * 	- The length of the lower element is kept
 * 	- Keep the birthday of the lower element
 * @param *Edata0 Properties of lower element
 * @param *Edata1 Properties of upper element
 * @param join Tells whether upper element is either joined with lower one or simply removed
 */
void SnowStation::mergeElements(ElementData& Edata0, const ElementData& Edata1, const bool& join)
{
	double L0, L1, LNew; // Lengths of lower (e0), upper (e1) elements, and "new" element, respectively

	L1 = Edata1.L;
	L0 = Edata0.L;
	LNew = L0;

	if (join) {
		LNew += L1;
		Edata0.depositionDate = Edata1.depositionDate;
		Edata0.rg = 0.5 * ( Edata0.rg + Edata1.rg );
		Edata0.dd = 0.5 * ( Edata0.dd + Edata1.dd );
		Edata0.sp = 0.5 * ( Edata0.sp + Edata1.sp );
		Edata0.rb = 0.5 * ( Edata0.rb + Edata1.rb );
		Edata0.CDot = 0.5 * ( Edata0.CDot + Edata1.CDot );
		Edata0.E = Edata0.Ev;
		Edata0.Ee = 0.0; // TODO (very old) Check whether not simply add the elastic
		                 //                 and viscous strains of the elements and average the stress?
	} else {
		Edata0.Ee = Edata0.E = Edata0.Ev = Edata0.dE = 0.0;
	}

	Edata0.L0 = Edata0.L = LNew;
	Edata0.M += Edata1.M;
	Edata0.theta[ICE] = (L1*Edata1.theta[ICE] + L0*Edata0.theta[ICE]) / LNew;
	Edata0.snowResidualWaterContent();
	// TODO Check whether we have space to accomodate all water percolating into Edata0
	Edata0.theta[WATER] = (L1*Edata1.theta[WATER] + L0*Edata0.theta[WATER]) / LNew;
	Edata0.theta[AIR] = 1.0 - Edata0.theta[WATER] - Edata0.theta[ICE];
	Edata0.Rho = (Edata0.theta[ICE]*Constants::density_ice) + (Edata0.theta[WATER]*Constants::density_water);

	for (int ii = 0; ii < signed(SnowStation::number_of_solutes); ii++) {
		for (unsigned int kk = 0; kk < N_COMPONENTS; kk++) {
			Edata0.conc[kk][ii] = (L1*Edata1.conc(kk,ii) + L0*Edata0.conc[kk][ii]) / LNew;
		}
	}
	Edata0.dth_w = (L1*Edata1.dth_w + L0*Edata0.dth_w) / LNew;
	Edata0.Qmf += Edata1.Qmf;
	Edata0.sw_abs += Edata1.sw_abs;
	if ((Edata1.mk >= 100) && (Edata0.mk < 100)) {
		Edata0.mk += (Edata1.mk/100)*100;
	}
}

CurrentMeteo::CurrentMeteo(const size_t& i_max_number_of_sensors)
	: n(0), date(), ta(0.), rh(0.), rh_avg(0.), vw(0.), vw_avg(0.), vw_max(0.), dw(0.),
	  vw_drift(0.), dw_drift(0.), ustar(0.), z0(0.), psi_s(0.),
	  iswr(0.), rswr(0.), diff(0.), elev(0.), ea(0.), tss(0.), tss_a12h(0.), tss_a24h(0.), ts0(0.),
	  hnw(0.), hs(0.), hs_a3h(0.), hs_rate(0.),
	  rho_hn(0.),
	  max_number_of_sensors(i_max_number_of_sensors)
{
	ts    = vector<double>(max_number_of_sensors, Constants::nodata);
	zv_ts = vector<double>(max_number_of_sensors, Constants::nodata);
	conc  = vector<double>(SnowStation::number_of_solutes, 0.);
}

void CurrentMeteo::reset()
{
	*this = CurrentMeteo(max_number_of_sensors);
}

std::ostream& operator<<(std::ostream &os, const CurrentMeteo& mdata)
{
	os << "<InterpolatedMeteoData>" << endl;

	os << "n:        " << mdata.n << endl;
	os << "Date:     " << mdata.date.toString(Date::ISO) << endl;
	os << "TA:       " << mdata.ta << endl;
	os << "RH:       " << mdata.rh << endl;
	os << "rh_avg:   " << mdata.rh_avg << endl;
	os << "VW:       " << mdata.vw << endl;
	os << "VW_MAX:   " << mdata.vw_max << endl;
	os << "DW:       " << mdata.dw << endl;
	os << "vw_avg:   " << mdata.vw_avg << endl;
	os << "USTAR:    " << mdata.ustar << endl;
	os << "z0:       " << mdata.z0 << endl;
	os << "psi_s:    " << mdata.psi_s << endl;
	os << "ISWR:     " << mdata.iswr << endl;
	os << "RSWR:     " << mdata.rswr << endl;
	os << "diff:     " << mdata.diff << endl;
	os << "ELEV:     " << mdata.elev << endl;
	os << "EA:       " << mdata.ea << endl;
	os << "TSS:      " << mdata.tss << endl;
	os << "TSG:      " << mdata.ts0 << endl;
	os << "HNW:      " << mdata.hnw << endl;
	os << "HS:       " << mdata.hs << endl;
	stringstream ss;
	for (unsigned int ii=1; ii<mdata.ts.size(); ii++) {
		ss << ""; ss << ii;
		os << "ts[" << ss.str() << "]:    " << mdata.ts[ii] << endl;
		os << "zv_ts[" << ss.str() << "]: " << mdata.zv_ts[ii] << endl;
	}
	for (unsigned int ii=1; ii<mdata.conc.size(); ii++) {
		ss << ""; ss << ii;
		os << "conc[" << ss.str() << "]: " << mdata.conc[ii] << endl;
	}
	os << "vw_drift: " << mdata.vw_drift << endl;
	os << "rho_hn:   " << mdata.rho_hn << endl;

	os << "</InterpolatedMeteoData>" << endl;
	return os;
}


LayerData::LayerData() : layerDate(), hl(0.), ne(0), tl(0.),
                     phiIce(0.), phiWater(0.), phiVoids(0.), phiSoil(0.), SoilRho(0.), SoilK(0.), SoilC(0.),
                     rg(0.), sp(0.), dd(0.), rb(0.), mk(0), hr(0.), CDot(0.), metamo(0.)
{
	cIce.resize(SnowStation::number_of_solutes);
	cWater.resize(SnowStation::number_of_solutes);
	cVoids.resize(SnowStation::number_of_solutes);
	cSoil.resize(SnowStation::number_of_solutes);
}
