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
 * @file Snowpack.c
 * @version 10.02
 * @bug     -
 * @brief This module contains the driving routines for the 1d snowpack model
 */

#include <snowpack/DataClasses.h>
#include <snowpack/Utils.h>
#include <snowpack/Canopy.h>
#include <snowpack/Metamorphism.h>
#include <snowpack/Aggregate.h>

using namespace mio;
using namespace std;

SN_SURFACE_DATA::SN_SURFACE_DATA(const unsigned int& i_max_number_of_solutes) 
  : dIntEnergy(0.), lw_in(0.), lw_out(0.), lw_net(0.), qs(0.), ql(0.), hoar(0.), qr(0.), qg(0.), qg0(0.), sw_hor(0.),
    sw_in(0.), sw_out(0.), qw(0.), sw_dir(0.), sw_diff(0.), cA(0.), mA(0.), drift(0.), dhs_corr(0.), 
    max_number_of_solutes(i_max_number_of_solutes)
{
	mass.resize(N_MASS_CHANGES);
	load.resize(max_number_of_solutes);
}

void SN_SURFACE_DATA::reset(const bool& cumsum_mass)
{
	if (cumsum_mass) { // Do not reset cumulated mass balance
		dIntEnergy = 0.;
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
		mass[MS_PRECIP] = 0.;
		mass[MS_RAIN] = 0.;
	} else {
		*this = SN_SURFACE_DATA(max_number_of_solutes); //reset everything
	}
}

void SN_CANOPY_DATA::reset(const bool& cumsum_mass)
{
	if ( cumsum_mass ) { // Do not reset cumulated mass balance
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
void SN_CANOPY_DATA::initializeSurfaceExchangeData()
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

SN_STATION_DATA::SN_STATION_DATA(const bool& i_useCanopyModel, const bool& i_useSnowLayers) : 
	Lat(0.), Lon(0.), Alt(0.), 
	SlopeAzi(0.), SlopeAngle(0.), Albedo(0.), SoilAlb(0.), BareSoil_z0(0.), SoilNode(0), cH(0.),
	mH(0.), Ground(0.), hn_slope(0.), rho_slope(0.), windward(0), ErosionLevel(0), ErosionMass(0.), 
	S_class1(0), S_class2(0), S_d(0.), z_S_d(0.), S_n(0.), z_S_n(0.), S_s(0.), z_S_s(0.), S_4(0.), 
	z_S_4(0.), S_5(0.), z_S_5(0.), Kt(NULL), Ks(NULL), ColdContent(0.), 
	SubSurfaceMelt('x'), SubSurfaceFrze('x'), Cdata(), tag_low(0),
	useCanopyModel(i_useCanopyModel), useSnowLayers(i_useSnowLayers), nNodes(0), nElems(0)
{
	Edata = vector<SN_ELEM_DATA>();
	Ndata = vector<SN_NODE_DATA>();
}

/**
 * @brief Reallocate element and node data \n
 * Xdata->Edata, Xdata->Ndata and Xdata->nElems, Xdata->nNodes are reallocated or reset, respectively.
 * In case of augmenting the element number, the new elements are initialized to 0 (memset)
 * @param **Xdata
 * @param NewnE New number of elements
 * @return ERROR if any, NO_ERROR if not
 */
void SN_STATION_DATA::resize(const int& number_of_elements)
{
	int oldNE = getNumberOfElements();
	int dnE = number_of_elements - oldNE;

	try {
		Edata.resize(number_of_elements);
		Ndata.resize(number_of_elements + 1);
	}catch(exception& e){
		throw IOException(e.what(), AT); //this will catch all allocation exceptions
	}

	if (dnE > 0) {
		memset(&Edata[oldNE], 0, (dnE)*sizeof(SN_ELEM_DATA));
	}

	nElems = (int)Edata.size();
	nNodes = (int)Ndata.size();
}

int SN_STATION_DATA::getNumberOfElements() const
{
	return nElems;
}

int SN_STATION_DATA::getNumberOfNodes() const
{
	return nNodes;
}


/**
 * @brief If more than NUMBER_TOP_ELEMENTS snow elements exist, attempt to reduce their number in the FEM mesh,
 * leaving NUMBER_TOP_ELEMENTS surface elements untouched \n
 * Pairs of elements within the snow cover satisfying the join conditions of sn_JoinCondition() are merged
 * by placing everything in the lower element, setting the density of upper element to NODATA,
 * and getting rid of node in between. \n
 * The elements being very similar and thus the microstructure parameters being approximately equal
 * as defined in sn_JoinCondition(), simply average the microstructure properties \n
 * NOTE that the condense element check is placed at the end of a time step, allowing elements do develop on their own.
 * @param *Xdata
 * @return ERROR if any, NO_ERROR if not
 */
void SN_STATION_DATA::joinElements(const int& number_top_elements)
{
	int e0, e1;  // Lower (e0) and upper (e1) element index
	int nE, rnE; // Original and reduced number of elements
	int nJoin=0; // Number of elements to be removed

	nE = nElems;
	if ( nE - SoilNode < number_top_elements+1 ) {
		return;
	}
	for (e0 = SoilNode, e1 = SoilNode+1; e0 < nE-number_top_elements; e0++, e1++) {
	  if ( SN_STATION_DATA::sn_JoinCondition(Edata[e0], Edata[e1]) ) {
	    SN_STATION_DATA::mergeElements(Edata[e0], Edata[e1], true);
			nJoin++;
			Edata[e1].Rho *= -1.;
			e0++; e1++;
		}
	}
	if ( nJoin > 0 ) {
		rnE = nE - nJoin;
		reduceNumberOfElements(rnE);
	}
}

/**
 * @brief Remove the upper "marked" element of two (snow only) \n
 * -# Joining two elements:
 *  - make density negative (*= -1.)
 *  - take the uppermost node of both
 * -# Removing melted or thin elements
 *  - make both density AND length negative (*= -1.) as the latter will be used!
 *  - keep upper node of lowest element
 * @param rnE Reduced number of elements
 */
void SN_STATION_DATA::reduceNumberOfElements(const int& rnE)
{
	int e0;                    // Lower element index
	int eNew;                  // New element index
	double cH_old, dL=0.;
	
	for (e0 = SoilNode, eNew = SoilNode; e0 < nElems; e0++) {
		if ( Edata[e0].Rho < 0.0 ) {
			if ( Edata[e0].L > 0.0 ) { // Joining elements
				Ndata[eNew] = Ndata[e0+1];
				Ndata[eNew].z = Ndata[e0+1].z + Ndata[e0+1].u + dL;
				Ndata[eNew].u = Ndata[e0].udot = 0.;
			} else { // Removing elements
				dL += Edata[e0].L;
			}
		} else {
			if ( eNew < e0 ) {
				Edata[eNew] = Edata[e0];
				Ndata[eNew+1] = Ndata[e0+1];
			}
			Ndata[eNew+1].z = Ndata[e0+1].z + Ndata[e0+1].u + dL;
			Ndata[eNew+1].u = Ndata[e0+1].udot = 0.;
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
 * The bottom temperature at the beginning of the calculation is given by the temperature at the top of the
 * lowest soil or snow layer \n
 * IMPORTANT: it is very important for Alpine3D that Cdata.height is initialized even if CANOPY = 0,
 * otherwise SnowInterface will not recognize the canopy grids (David 2007-06-25).
 * @author Perry Bartelt \n Michael Lehning \n Charles Fierz
 * @version 10.02
 * @param SSdata
 */
void SN_STATION_DATA::initialize(const SN_SNOWSOIL_DATA& SSdata)
{
	int nE;               //  Number elements
	int e, i, l, n;       //  Counters: element, layer and node
	int le;               //  Number of elements per layer
	int real_soil_no_sandwich = 1;  // Switch to count real soil layers
	double dT, SigC;      //  Change in temperature and Cauchy stress

	Albedo = SSdata.Albedo;
	SoilAlb = SSdata.SoilAlb;
	BareSoil_z0 = SSdata.BareSoil_z0;
	SlopeAngle = DEG_TO_RAD(SSdata.Angle);
	SlopeAzi   = DEG_TO_RAD(SSdata.Azi);
	Alt = SSdata.Alt;
	Lat = SSdata.Lat;
	Lon = SSdata.Lon;
	mH = cH = SSdata.Height;

	nElems=0;

	resize(SSdata.nN-1);

	nE = nElems;
	SoilNode = 0;
	Ground = 0.0;
	Ndata[0].z = 0.;
	if ( SSdata.nLayers>0 ) {
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
	for (l = 0, n = 1; l < SSdata.nLayers; l++) {
		// Update ground heigth and SoilNode number
		if (SSdata.Ldata[l].phiSoil > 0.0 && real_soil_no_sandwich) {
			Ground += SSdata.Ldata[l].hl;
			SoilNode += SSdata.Ldata[l].ne;
		} else {
			real_soil_no_sandwich = 0;
		}
		if (l == 0) {
			dT = 0.;
		} else {
			dT = (SSdata.Ldata[l].tl - SSdata.Ldata[l-1].tl) / (SSdata.Ldata[l].ne);
		}
		for (i = 0;i < SSdata.Ldata[l].ne; i++, n++ ) {
			Ndata[n].z = Ndata[n-1].z + SSdata.Ldata[l].hl / SSdata.Ldata[l].ne;
			Ndata[n].T = Ndata[n-1].T + dT;
			Ndata[n].u = 0.;
			Ndata[n].f = 0.;
			Ndata[n].udot = 0.;
			Ndata[n].S_n = INIT_STABILITY;   // Static natural stability index
			Ndata[n].S_s = INIT_STABILITY;   // Alternative Stability Index (skier stability)
		}
	}

	if (SoilNode == 0 && useSnowLayers) {
		prn_msg(__FILE__, __LINE__, "err", -1., "SNP_SOIL set but no soil layers given");
		throw IOException("Snowpack Initialization failed", AT);
	}

	// INITIALIZE THE ELEMENT DATA
	for (l = 0, e = 0; l<SSdata.nLayers; l++) {
		for (le = 0; le < SSdata.Ldata[l].ne; le++, e++) {
			// Element's JulianQ Date
			Edata[e].date = SSdata.Ldata[l].date;
			// Temperature data
			Edata[e].Te=(Ndata[e].T+Ndata[e+1].T)/2.;
			Edata[e].L0=Edata[e].L=(Ndata[e+1].z-Ndata[e].z);
			Edata[e].gradT=(Ndata[e+1].T-Ndata[e].T)/Edata[e].L;
			// Creep data
			Edata[e].E=Edata[e].S=Edata[e].EDot=0.0;
			Edata[e].Ev=Edata[e].Ee=Edata[e].EvDot=0.0;
			// Very important to initialize the increments in length and strain
			Edata[e].dE = 0.0;
			// Volumetric Components  (If the density is over 1000 assume PERMAFROST)
			Edata[e].theta[AIR] = SSdata.Ldata[l].phiVoids;
			Edata[e].theta[ICE] = SSdata.Ldata[l].phiIce;
			Edata[e].theta[WATER] = SSdata.Ldata[l].phiWater;
			for (i = 0; i < N_SOLUTES; i++) {
				Edata[e].conc[ICE][i]  = SSdata.Ldata[l].cIce[i];
				Edata[e].conc[WATER][i] = SSdata.Ldata[l].cWater[i];
				Edata[e].conc[AIR][i]  = SSdata.Ldata[l].cVoids[i];
			}
			Edata[e].theta[SOIL]    = SSdata.Ldata[l].phiSoil;
			Edata[e].soil[SOIL_RHO] = SSdata.Ldata[l].SoilRho;
			Edata[e].soil[SOIL_K]   = SSdata.Ldata[l].SoilK;
			Edata[e].soil[SOIL_C]   = SSdata.Ldata[l].SoilC;
			Edata[e].Rho = Edata[e].theta[ICE]*Constants::density_ice + 
				Edata[e].theta[WATER]*Constants::density_water + Edata[e].theta[SOIL]*Edata[e].soil[SOIL_RHO];
			for (i = 0; i < N_SOLUTES; i++) {
				Edata[e].conc[SOIL][i]  = SSdata.Ldata[l].cSoil[i];
			}
			// Zero out conductivities, specific heat and moisture content
			Edata[e].k[0] = Edata[e].k[1] = Edata[e].k[2] = 0.0;
			Edata[e].c[0] = Edata[e].c[1] = Edata[e].c[2] = 0.0;
			// Set the initial short wave radiation to zero
			Edata[e].sw_abs = 0.0;
			// Phase change variables
			Edata[e].Qmf=0.0;
			Edata[e].dth_w=0.0;
			// Micro-structure data
			Edata[e].dd = SSdata.Ldata[l].dd;
			Edata[e].sp = SSdata.Ldata[l].sp;
			Edata[e].rg = SSdata.Ldata[l].rg;
			Edata[e].rb = SSdata.Ldata[l].rb;
			Edata[e].N3 = Metamorphism::getCoordinationNumberN3(Edata[e].Rho);
			Edata[e].mk = SSdata.Ldata[l].mk;
			Edata[e].type = ml_ag_Classify(SSdata.Ldata[l].dd, SSdata.Ldata[l].sp, 2.*SSdata.Ldata[l].rg, SSdata.Ldata[l].mk%100, SSdata.Ldata[l].phiWater, SSdata.Ldata[l].phiIce);
			Ndata[e+1].hoar = SSdata.Ldata[l].hr;
			// Initialize the Layer Stability Index  (deformation rate)
			Edata[e].S_dr = INIT_STABILITY;
			Edata[e].hard = 0.0;
			// Finally the INITIAL MASS of the ELEMENT which is used to make sure that the mass balance of the snowpack is CORRECT.
			Edata[e].M = Edata[e].Rho * Edata[e].L0;
		} // end of element layer for
	} // end of layer for
	if ( SSdata.ErosionLevel > 0 ) {
		ErosionLevel = SSdata.ErosionLevel;
	} else {
		ErosionLevel = MAX(SoilNode, nE-1);
	}
	// Find the real Cauchy stresses
	for (e = nE-1, SigC = 0.0; e >=0; e--) {
		SigC += -(Edata[e].M)*Constants::g*cos(SlopeAngle);
		Edata[e].C = SigC;
	}

	// INITIALIZE CANOPY DATA
	if ( useCanopyModel ) {
		Cdata.height=SSdata.Canopy_Height;
		Cdata.storage=0.0;           // intercepted water (kg m-2 or mm Water Equivalent)
		Cdata.temp=273.15;	          // temperature (K)
		Cdata.canopyalb=CAN_ALB_DRY; // albedo [-], which is a function of the dry canopy albedo and intercepted snow
		Cdata.wetfraction=0.0;
		Cdata.lai=SSdata.Canopy_LAI;

		Cdata.sigf=1.-exp(-KRNT_LAI*(Cdata.lai)); // radiation transmissivity (-)
		Cdata.ec=1.0;               //longwave emissivity

		Cdata.z0m=Cdata.height*0.1;
		Cdata.z0h=Cdata.z0m*0.1;
		Cdata.zdispl=Cdata.height*0.66;
		Cdata.direct_throughfall=SSdata.Canopy_Direct_Throughfall;
		if ( !(SSdata.Canopy_Direct_Throughfall >= 0. && SSdata.Canopy_Direct_Throughfall <= 1.) ) {
			prn_msg(__FILE__, __LINE__, "err", -1., "Given Canopy Throughfall (*.sno file) = %lf but Canopy is set", 
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
		Cdata.canopyalb=CAN_ALB_DRY; // albedo [-], which is a function of the dry canopy albedo and intercepted snow
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
 * 	- larger than JOIN_THRESH_L
 * 	- tagged
 * 	- dry surface hoar (mk=3),
 * 	- dendritic but not both
 * -# OTHERWISE we use criteria for dendricity, sphericity, volumetric ice or water content, grain size and marker
 * - NOTE July 2006: whatever type of thin elements are treated in wt_ElementRemoval()
 * @param Edata0 Lower element
 * @param Edata1 Upper element
 * @return true if the two elements should be joined, false otherwise
 */
bool SN_STATION_DATA::sn_JoinCondition(const SN_ELEM_DATA& Edata0, const SN_ELEM_DATA& Edata1)
{
	if ( (Edata0.L > JOIN_THRESH_L) || (Edata1.L > JOIN_THRESH_L) )
		return false;
	
	if ( Edata0.mk%100 != Edata1.mk%100 )
		return false;

	if ( fabs(Edata0.sp - Edata1.sp) > JOIN_THRESH_SP )
		return false;

	if ( Edata0.theta[SOIL] > 0. || Edata1.theta[SOIL] > 0. )
		return false;

	if ( (Edata0.mk >= 100) || (Edata1.mk >= 100) )
		return false;

	if ( (Edata0.mk%100 == 3) || (Edata1.mk%100 == 3) )
		return false;

	if ( (Edata0.dd > JOIN_THRESH_DD || Edata1.dd > JOIN_THRESH_DD) && 
		!(Edata0.dd > JOIN_THRESH_DD && Edata1.dd > JOIN_THRESH_DD) ) {
		return false;
	} else if ( fabs(Edata0.dd - Edata1.dd) > JOIN_THRESH_DD ) {
		return false;
	}

	if ( fabs(Edata0.theta[ICE] - Edata1.theta[ICE]) > JOIN_THRESH_ICE )
		return false;

	if ( fabs(Edata0.theta[WATER] - Edata1.theta[WATER]) > JOIN_THRESH_WATER )
		return false;

	if ( fabs(Edata0.rg - Edata1.rg) > JOIN_THRESH_RG )
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
 * 	- NOTE Joining two elements may cause the tags (marker >= 100) to "jump" upwards
 * - Removing:
 * 	- Remaining ice, liquid water, solutes, etc. are added to the lower element
 * 	- The length of the lower element is kept
 * 	- Keep the birthday of the lower element
 * @param *Edata0 Properties of lower element
 * @param *Edata1 Properties of upper element
 * @param join Tells whether upper element is either joined with lower one or simply removed
 */
void SN_STATION_DATA::mergeElements(SN_ELEM_DATA& Edata0, const SN_ELEM_DATA& Edata1, const bool& join)
{
	int i, k;
	double L0, L1, LNew; // Length of lower (e0) and upper (e1) elements, respectively

	L1 = Edata1.L;
	L0 = Edata0.L;
	LNew = L0;

	if ( join ) {
		LNew += L1;
		Edata0.date = Edata1.date;
		Edata0.rg = 0.5 * ( Edata0.rg + Edata1.rg );
		Edata0.dd = 0.5 * ( Edata0.dd + Edata1.dd );
		Edata0.sp = 0.5 * ( Edata0.sp + Edata1.sp );
		Edata0.rb = 0.5 * ( Edata0.rb + Edata1.rb );
		Edata0.E = Edata0.Ev;
		Edata0.Ee = 0.0; // TODO Check whether not simply add the elastic and viscous strains of the elements and average the stress?
	} else {
		Edata0.Ee = Edata0.E = Edata0.Ev = Edata0.dE = 0.0;
	}

	Edata0.L0 = Edata0.L = LNew;
	Edata0.M += Edata1.M;
	Edata0.theta[ICE] = (L1*Edata1.theta[ICE] + L0*Edata0.theta[ICE]) / LNew;
	Edata0.theta[WATER] = (L1*Edata1.theta[WATER] + L0*Edata0.theta[WATER]) / LNew;
	Edata0.theta[AIR] = 1.0 - Edata0.theta[WATER] - Edata0.theta[ICE];
	Edata0.Rho = (Edata0.theta[ICE]*Constants::density_ice) + (Edata0.theta[WATER]*Constants::density_water);
	
	for (i = 0; i < N_SOLUTES; i++) {
		for (k = 0; k < N_COMPONENTS; k++) {
			Edata0.conc[k][i] = (L1*Edata1.conc[k][i] + L0*Edata0.conc[k][i]) / LNew;
		}
	}
	Edata0.dth_w = (L1*Edata1.dth_w + L0*Edata0.dth_w) / LNew;
	Edata0.Qmf += Edata1.Qmf;
	Edata0.sw_abs += Edata1.sw_abs;
	if ( (Edata1.mk >= 100) && (Edata0.mk < 100) ) {
		Edata0.mk += (Edata1.mk/100)*100;
	}
}


SN_MET_DATA::SN_MET_DATA(const unsigned int& i_max_number_of_sensors, const unsigned int& i_max_number_of_solutes) 
  : n(0), date(), ta(0.), rh(0.), rh_ave(0.), vw(0.), vw_ave(0.), vw_max(0.), dw(0.), ustar(0.), z0(0.), 
    psi_s(0.), iswr(0.), rswr(0.), diff(0.), elev(0.), ea(0.), tss(0.), ts0(0.), hnw(0.), hs1(0.), rho_hn(0.), 
    max_number_of_sensors(i_max_number_of_sensors), max_number_of_solutes(i_max_number_of_solutes)
{
	ts    = vector<double>(max_number_of_sensors, 0.);
	zv_ts = vector<double>(max_number_of_sensors, 0.);
	conc  = vector<double>(max_number_of_solutes, 0.);
}

void SN_MET_DATA::reset()
{
	*this = SN_MET_DATA(max_number_of_sensors, max_number_of_solutes);
}


std::ostream& operator<<(std::ostream &os, const SN_MET_DATA& mdata)
{
	os << "<InterpolatedMeteoData>" << endl;

	os << "n:        " << mdata.n << endl;
	os << "Date:     " << mdata.date.toString(Date::ISO) << endl;
	os << "TA:       " << mdata.ta << endl;
	os << "RH:       " << mdata.rh << endl;
	os << "RH_AVE:   " << mdata.rh_ave << endl;
	os << "VW:       " << mdata.vw << endl;
	os << "VW_AVE:   " << mdata.vw_ave << endl;
	os << "VW_MAX:   " << mdata.vw_max << endl;
	os << "DW:       " << mdata.dw << endl;
	os << "USTAR:    " << mdata.ustar << endl;
	os << "z0:       " << mdata.z0 << endl;
	os << "psi_s:    " << mdata.psi_s << endl;
	os << "ISWR:     " << mdata.iswr << endl;
	os << "RSWR:     " << mdata.rswr << endl;
	os << "diff:     " << mdata.diff << endl;
	os << "ELEV:     " << mdata.elev << endl;
	os << "EA:       " << mdata.ea << endl;
	os << "TSS:      " << mdata.tss << endl;
	os << "TS0:      " << mdata.ts0 << endl;
	os << "HNW:      " << mdata.hnw << endl;
	os << "HS1:      " << mdata.hs1 << endl;
	os << "ts[0]:    " << mdata.ts[0] << endl;
	os << "zv_ts[0]: " << mdata.zv_ts[0] << endl;
	os << "rho_hn:   " << mdata.rho_hn << endl;

	os << "</InterpolatedMeteoData>" << endl;
	return os;
}


SN_LAYER_DATA::SN_LAYER_DATA(const unsigned int& i_max_number_of_solutes) : date(0.), hl(0.), ne(0), tl(0.), phiIce(0.), 
					    phiWater(0.), phiVoids(0.), phiSoil(0.), SoilRho(0.), SoilK(0.), SoilC(0.), rg(0.), sp(0.), 
                             dd(0.), rb(0.), mk(0), hr(0.), max_number_of_solutes(i_max_number_of_solutes)
{
	cIce.resize(max_number_of_solutes);
	cWater.resize(max_number_of_solutes);
	cVoids.resize(max_number_of_solutes);
	cSoil.resize(max_number_of_solutes);
}
