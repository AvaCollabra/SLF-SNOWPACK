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
 * @file Snowpack.h
 * @version 10.02
 * This header file contains all the data structures needed for the 1d snowpack model
 */

#ifndef __DATACLASSES_H__
#define __DATACLASSES_H__

#include <snowpack/Constants.h>
#include <meteoio/MeteoIO.h>
#include <vector>
#include <sstream>

/// @brief The 3 different phases in the matrix
enum {
	SOLID,  ///< Solid
	LIQUID, ///< Liquid
	VOID,   ///< Gas
	NUMBER_OF_PHASES
};


/// @brief Optical Properties of snow
struct SNOW_OPTIC {
	double ggg;     ///< Asymmetry Parameter
	double exteff;  ///< Extinction Efficiency
	double ssa;     ///< Single Scattering Albedo
};

/// @brief Spectrum of incoming solar radiation
struct WL_STRUCT {
	double nm;      ///< Wavelength
	double perc;    ///< Percentage of Energy
};

/**
 * SurfaceFluxes contains all surface exchange data \n
 * Some of the most important results of the simulation are contained in this data structure.
 */
class SurfaceFluxes {
	public:
		/**
		 * @brief The different types of Mass fluxes
		 * Mass fluxes in kg m-2, computed and output in qr_WriteTimeSeries() \n
		 * Rates in kg m-2 h-1 (MS_HNW, MS_RAIN and MS_WIND)
		 */
		enum SN_MASS_CHANGES {
			MS_TOTALMASS,      ///<  This of course is the total mass of the snowpack at the present time
			MS_SWE,            ///<  This too, of course, but summing rho*L
			MS_WATER,          ///<  The total amount of water in the snowpack at the present time
			MS_HNW,            ///<  Solid precipitation rate
			MS_RAIN,           ///<  Rain rate
			MS_WIND,           ///<  Mass loss rate due to wind erosion
			MS_EVAPORATION,    ///<  The mass loss or gain of the top element due to water evaporating
			MS_SUBLIMATION,    ///<  The mass loss or gain of the top element due to snow (ice) sublimating
			MS_RUNOFF,         ///<  The total mass loss due to surface runoff; used to check mass balance of snowpack and also to compute cummulative discharge
			MS_SOIL_RUNOFF,    ///<  Equivalent to MS_RUNOFF but at bottom soil node
			MS_CORRECTION,     ///< Mass correction from either squeezing (neg) or blowing up (pos)
			N_MASS_CHANGES     ///< Total number of different mass change types
		};

		SurfaceFluxes();

		void reset(const bool& cumsum_mass);

		double dIntEnergy; ///< change of internal energy
		double lw_in;    ///< incoming longwave radiation
		double lw_out;   ///< outgoing longwave radiation
		double lw_net;   ///< net longwave radiation
		double qs;       ///< sensible heat
		double ql;       ///< latent heat
		double hoar;     ///< mass of surface hoar formed or sublimated
		double qr;       ///< rain energy
		double qg;       ///< geothermal heat flux or heat flux at bottom of combined soil-snow pack
		double qg0;      ///< ground heat flux at soil-snow interface
		double sw_hor;   ///< incoming global shortwave radiation in W m-2 on horizontal surface
		double sw_in;    ///< incoming global shortwave radiation in W m-2; on slopes projected
		double sw_out;   ///< reflected shortwave radiation in W m-2
		double qw;       ///< absorbed shortwave radiation in W m-2 on surface (net SW)
		double sw_dir;   ///< incoming direct shortwave radiation in W/m2; on slopes: projected
		double sw_diff;  ///< incoming diffuse shortwave radiation in W m-2
		double cA;       ///< computed Albedo (USED only for OUTPUT)
		double mA;       ///< measured Albedo (USED only for OUTPUT)
		double drift;    ///< the surface flux of drifting snow in kg m-1 s-1
		std::vector<double> mass; ///< Total mass of snowpack PLUS different amounts of total mass change, sublimation, runoff, erosion, etc. Basically the mass which crosses the surface
		std::vector<double> load; ///< Total load (kg m-2) in water runoff from solutes like nitrate
		double dhs_corr; ///< Snow depth correction in case of squezzing or blow-up
};

/**
 * @brief sn_Zdata contains "memory" information \n
 * It is used to prepare some parameters of qr_Hdata. This data is read from and written to *.sno files
 */
struct SN_ZWISCHEN_DATA {
	double hoar24[48];  ///< Twenty-four hour hoar index every half-hour over one day
	double drift24[48]; ///< Twenty-four hour drift index every half-hour over one day
	double hn3[144];    ///< Three hour new snow heights every half-hour over three days
	double hn24[144];   ///< Twenty-four hour snow heights every half-hour over three days
};

/**
 * @brief CurrentMeteo is the class of interpolated meteo data for the current calculation time step \n
 * It contains some additional and very important derived parameters such as the roughness length or running mean values.
 */
class CurrentMeteo {

	public:
		CurrentMeteo(const unsigned int& i_max_number_of_sensors);
		void reset();

		friend std::ostream& operator<<(std::ostream& os, const CurrentMeteo& mdata);

		int    n;      ///< record number of basic meteo input data

		mio::Date date;        ///< Date of current meteo data
		double ta;     ///< Air temperature (K)
		double rh;     ///< Relative humidity (% or 1)
		double rh_ave; ///< Running mean of relative humidity (1)
		double vw;     ///< Wind velocity at snow station (m s-1)
		double vw_ave; ///< Running mean of wind velocity at snow station (m s-1)
		double vw_max; ///< Maximum wind velocity at snow station (m s-1)
		double dw;     ///< Wind direction at snow station (deg)
		double vw_drift; ///< Wind velocity for blowing and drifting snow (operational: wind ridge station)
		double dw_drift; ///< Wind direction of blowing and drifting snow (operational: wind ridge station)
		double ustar;  ///< The friction velocity (m s-1) computed in mt_MicroMet() and also used later for the MeteoHeat fluxes
		double z0;     ///< The roughness length computed in SnowDrift and also used later for the MeteoHeat fluxes (m)
		double psi_s;  ///< Stability correction for scalar heat fluxes
		double iswr;   ///< Incoming SHORTWAVE radiation (W m-2)
		double rswr;   ///< Reflected SHORTWAVE radiation (W m-2) divide this value by the ALBEDO to get iswr
		double diff;   ///< Diffuse radiation from the sky (W m-2)
		double elev;   ///< Solar elevation to be used in Canopy.c (rad) => see also
		double ea;     ///< Atmospheric emissivity (1)
		double tss;    ///< Snow surface temperature (K)
		double ts0;    ///< Bottom temperatures of snow/soil pack (K)
		double hnw;    ///< The water equivalent of snowfall in mm w.e. (kg m-2) per CALCULATION_STEP_LENGTH
		double hs1;    ///< The measured height of snow (m), corrected for spikes etc. in ml_co_Control()

		std::vector<double> ts;    ///< Snowpack and/or Soil temperatures (K)
		std::vector<double> zv_ts; ///< Depth of temperature sensors (m)
		std::vector<double> conc;  ///< Solute concentrations in precipitation
		double rho_hn;             ///< Measured new sno density (kg m-3)

	private:
		unsigned int max_number_of_sensors;
};

///@brief BoundCond is a subset of surface energy exchange data required to set Neumann boundary conditions
class BoundCond {
	public:
		BoundCond() : lw_out(0.), lw_net(0.), qs(0.), ql(0.), qr(0.), qg(0.) {};

		double lw_out;  ///< outgoing longwave radiation
		double lw_net;  ///< net longwave radiation
		double qs;      ///< sensible heat
		double ql;      ///< latent heat
		double qr;      ///< rain energy
		double qg;      ///< geothermal heat flux or heat flux at bottom of combined soil-snow pack
};

/// @brief The 3 mathematical fields that can be solved
enum SN_FIELD{
	TEMPERATURE, ///< Temperature (phase change and metamorphism)
	SEEPAGE,     ///< Water transport
	SETTLEMENT,  ///< Creep displacements
	N_SN_FIELDS
};

/// @brief The 4 different components of the porous matrix
enum {
	SOIL,  ///< Soil
	ICE,   ///< Ice
	WATER, ///< Water
	AIR,   ///< Air
	N_COMPONENTS
};

/// @brief The different soil fields (PERMAFROST)
enum SN_SOIL_DATA{
	SOIL_RHO, ///< Bulk density of dry soil component (without accounting for porosity)
	SOIL_K,   ///< Bulk conductivity of dry soil component
	SOIL_C,   ///< Specific heat of dry soil component
	N_SOIL_FIELDS
};

/**
 * @brief Parameters of the different layers of the snowpack \n
 * It is used as a pointer (array) within the sn_SSdata (profile) data structure.
 */
class LayerData {

	public:
		LayerData(); 

		mio::Date layerDate;        ///< Date of deposition
		double hl;                  ///< The height of the layer in m
		unsigned int ne;            ///< Number of finite elements in the the layer (hl/ne defines elm. size)
		double tl;                  ///< Temperature at the top of the layer in K or degC
		double phiIce;              ///< Volumetric ice content in %
		double phiWater;            ///< Volumetric water content in %
		double phiVoids;            ///< Volumetric void content in %
		double phiSoil;             ///< Volumetric soil content in %
		std::vector<double> cIce;   ///< Solute concentrations in Ice
		std::vector<double> cWater; ///< Solute concentrations in Water
		std::vector<double> cVoids; ///< Solute concentrations in Air
		std::vector<double> cSoil;  ///< Solute concentrations in Soil
		double SoilRho;             ///< Density of soil in kg m-3
		double SoilK;               ///< Conductivity of soil
		double SoilC;               ///< Heat Capacity of soil
		double rg;                  ///< Micro-structure : Grainsize in mm
		double sp;                  ///< Micro-structure : Sphericity
		double dd;                  ///< Micro-structure : Dendricity
		double rb;                  ///< Micro-structure : Bond Radius in mm
		int mk;                     ///< Micro-structure : Marker
		double hr;                  ///< Surface hoar Mass in kg m-2
		double CDot;                ///< Stress rate (Pa s-1), that is the LAST overload change rate
		double metamo;              ///< keep track of metamorphism
};

/**
 * @brief sn_SSdata includes all important station parameters as well as LayerData \n
 * This data structure will have to be replaced by something a little more complicated soon ???
 * For now it is simply an efficient way of creating a snowpack to investigate.
 */
class SN_SNOWSOIL_DATA {

	public:

		SN_SNOWSOIL_DATA() : meta(), profileDate(0., 0.), nN(0), Height(0.),
                     nLayers(0), HS_last(0.), Albedo(0.), SoilAlb(0.), BareSoil_z0(0.),
                     Canopy_Height(0.), Canopy_LAI(0.), Canopy_Direct_Throughfall(0.),
                     WindScalingFactor(1.), ErosionLevel(0), TimeCountDeltaHS(0.)
		{
			Ldata.clear();
		}

		mio::StationData meta;            ///< Station meta data
		mio::Date profileDate;            ///< Date of profile
		unsigned int nN;                  ///< Total number of FE nodes
		double Height;                    ///< Total height of snowpack in m (sum of the layer heights)
		unsigned int nLayers;             ///< Total number of snowpack layers
		std::vector<LayerData> Ldata;     ///< contains all the information required to construct the Xdata
		double HS_last;                   ///< Last checked measured Snow Height
		double Albedo;                    ///< Snow albedo
		double SoilAlb;                   ///< Soil albedo; default 0.2
		double BareSoil_z0;               ///< Bare soil roughness in m; default 0.02 m
		double Canopy_Height;             ///< Canopy Height in m
		double Canopy_LAI;                ///< Canopy Leaf Area Index in m2 m-2
		double Canopy_Direct_Throughfall; ///< Direct throughfall [fraction of precipitation]
		double WindScalingFactor;         ///< Local scaling factor for wind at drift station
		int    ErosionLevel;              ///< Erosion Level in operational mode (flat field virtual erosion)
		double TimeCountDeltaHS;          ///< Time counter tracking erroneous settlement in operational mode
};

/**
 * @brief ELEMENT DATA used as a pointer in the SnowStation structure
 * NOTE on M below: this is the mass of an element that is neither changed by phase changes nor densification. \n
 * It is set in the data initialization and used to compute the stress field.
 * It can ONLY be changed by the WATER TRANSPORT or SURFACE SUBLIMATION or WIND TRANSPORT routines.
 */
class ElementData {
	
	public:
		ElementData();

		bool checkVolContent();
		void heatCapacity();
		double coldContent();
		double extinction();
		double snowResidualWaterContent();
		static double snowResidualWaterContent(const double theta_i);
		double soilFieldCapacity();

		double snowElasticity();
		double neckStressEnhancement();
		double concaveNeckRadius();
		double neckLength();
		double neck2VolumetricStrain();

		void snowType();
		static int snowType(const double dendricity, const double sphericity, const double grain_dia, const int marker,
                        const double theta_w, const double res_wat_cont);

		mio::Date depositionDate;  ///< Date of deposition
		double L0, L;              ///< Original and present element length (m)
		double Te;                 ///< mean element temperature (K)
		double gradT;              ///< temperature gradient over element (K m-1)
		std::vector<double> theta; ///< volumetric contents: SOIL, ICE, WATER, AIR (1)
		double Rho;                ///< mean element density (or BULK density; kg m-3), that is, rho=M/V=sum( theta(i)*rho(i) )
		mio::Array2D<double> conc; ///< Concentration for chemical constituents in (kg m-3)
		double M;                  ///< the total mass of the element (kg m-2)
		std::vector<double> k;     ///< For example, heat conductivity of TEMPERATURE field (W m-1 K-1)
		//   Stored in order to visualize constitutive laws
		//   Will be used for creep field hydraulic conductivity in m3 s kg-1
		std::vector<double> c;     ///< For example, specific heat of TEMPERATURE field (J kg)
		//   Will also be used for creep specific snow water capacity  in m3 J-1
		std::vector<double> soil;  ///< Contains the heat conductivity, capacity and dry density of the soil (solid, non-ice)  component phase
		double sw_abs;             ///< total absorbed shortwave radiation by the element (W m-2)
		// Snow Metamorphism Data
		double rg;                 ///< grain radius (mm)
		double dd;                 ///< snow dendricity: 0 = none, 1 = newsnow
		double sp;                 ///< sphericity: 1 = round, 0 = angular
		double rb;                 ///< grain bond radius (mm)
		double ps2rb;              ///< proportion of grain bond growth due to pressure sintering (1)
		double N3;                 ///< grain Coordination number (1)
		int    mk;                 ///< grain marker (history dependent)
		int    type;               ///< grain class
		double metamo;             ///< keep track of metamorphism
		double dth_w;              ///< Subsurface Melting & Freezing Data: change of water content
		double Qmf;                ///< Subsurface Melting & Freezing Data: change of energy due to phase changes (melt-freeze)
		double dE, E, Ee, Ev;      ///< Total element strain (GREEN'S strains -- TOTAL LAGRANGIAN FORMULATION.
		double EDot, EvDot;        ///< Total Strain Rate (s-1) (Simply, E/sn_dt)
		double S;                  ///< Total Element Stress (Pa), S being the energy conjugate stress
		double C;                  ///< Total Element Stress (Pa), C being the real or the Cauchy stress, which is output
		double CDot;               ///< Stress rate (Pa s-1), that is the overload change rate
		double S_dr;               ///< Stability Index based on deformation rate (Direct Action Avalanching)
		double s_strength;         ///< Parameterized snow shear strength (kPa)
		double hard;               ///< Parameterized hand hardness (1)
		//NIED (H. Hirashima)
		double dhf;
};

/// @brief NODAL DATA used as a pointer in the SnowStation structure
class NodeData {
	public:
		NodeData() : z(0.), u(0.), f(0.), udot(0.), T(0.), S_n(0.), S_s(0.), hoar(0.),
		             dhf(0.), S_dhf(0.), Sigdhf(0.) {}

		double z;    ///< nodal height from ground in m
		double u;    ///< creep displacements in m
		double f;    ///< reaction or unbalanced forces (CREEP)
		double udot; ///< downward creep velocity in m s-1
		double T;    ///< nodal temperature in K
		double S_n;  ///< Stability Index for natural avalanches
		double S_s;  ///< Stability Index for skier triggered avalanches
		double hoar; ///< Mass of surface hoar collected while node was exposed to surface

		//NIED (H. Hirashima)
		double dhf;
		double S_dhf;
		double Sigdhf;
};

/**
 * @brief Canopy data used as a pointer in the SnowStation structure
 * -# INSTANTANEOUS VARIABLES
 * 	-# Canopy "state" variables, and some auxiliaries
 * 	-# Properties which could be given here or as a parameter field
 * 	-# Aerodynamic resistances
 * -# CUMULATED/AVERAGE VARIABLES \n
 *    cumulated between time series output timesteps - these variables can be moved or mirrored in a canopy "surface" data structure
 * 	-# Albedo and similar auxiliaries
 * 	-# Radiation fluxes
 * 	-# Canopy turbulent heat fluxes
 * 	-# Canopy evaporative fluxes
 * 	-# Canopy mass fluxes
 */
class CanopyData {
	public:
		CanopyData() : storage(0.), temp(0.), sigf(0.), ec(0.), lai(0.), z0m(0.), z0h(0.), zdispl(0.),
		     height(0.), direct_throughfall(0.), ra(0.), rc(0.), rs(0.), rstransp(0.), canopyalb(0.),
		     totalalb(0.), wetfraction(0.), intcapacity(0.), rswrac(0.), iswrac(0.), rswrbc(0.),
		     iswrbc(0.), ilwrac(0.), rlwrac(0.), ilwrbc(0.), rlwrbc(0.), rsnet(0.), rlnet(0.),
		     sensible(0.), latent(0.), latentcorr(0.), transp(0.), intevap(0.), 
		     interception(0.), throughfall(0.), snowunload(0.) {}

		void reset(const bool& cumsum_mass);
		void initializeSurfaceExchangeData();

		// Aa
		double storage;  ///< intercepted water (mm or kg m-2)
		double temp;	   ///< temperature (K)
		double sigf;	   ///< radiation transmissivity (1)
		double ec;       ///< longwave emissivity (1)
		// Ab
		double lai;
		double z0m;
		double z0h;
		double zdispl;
		double height;
		double direct_throughfall;
		// Ac
		double ra;          ///< from canopy air to reference height
		double rc;          ///< from canopy to canopy air
		double rs;          ///< from subsurface to canpopy air
		double rstransp;    ///< stomatal surface resistance for transpiration
		// Ba
		double canopyalb;   ///< canopy albedo [-]
		double totalalb;    ///< total albedo above canopy and snow/soil surface [-]
		double wetfraction; ///< fraction of canopy covered by interception [-]
		double intcapacity; ///< maximum interception storage [mm]
		// Bb
		double rswrac;      ///< upward shortwave above canopy
		double iswrac;	    ///< downward shortwave radiation above canopy
		double rswrbc;      ///< upward shortwave below canopy
		double iswrbc;      ///< downward shortwave radiation below canopy
		double ilwrac;      ///< downward longwave radiation ABOVE canopy
		double rlwrac;      ///< upward longwave radiation ABOVE canopy
		double ilwrbc;      ///< downward longwave radiation BELOW canopy
		double rlwrbc;      ///< upward longwave radiation BELOW canopy
		double rsnet;       ///< net shortwave radiation
		double rlnet;       ///< net longwave radiation
		// Bc
		double sensible;
		double latent;
		double latentcorr;
		// Bd
		double transp;
		double intevap;
		// Be
		double interception;
		double throughfall;
		double snowunload;
};

/**
 * @brief Station data including all information on snowpack layers (elements and nodes) and on canopy \n
 * This is the PRIMARY data structure of the SNOWPACK program \n
 * It is used extensively not only during the finite element solution but also to control
 * the post-processing writes. It is initialized from SN_SNOWSOIL_DATA (at present).
 */
class SnowStation {

	public:
		SnowStation(const bool& i_useCanopyModel, const bool& i_useSoilLayers);

		void initialize(const SN_SNOWSOIL_DATA& SSdata);
		void resize(const unsigned int& number_of_elements);

		void reduceNumberOfElements(const unsigned int& rnE);
		void joinElements(const unsigned int& number_top_elements);
		static bool sn_joinCondition(const ElementData& Edata0, const ElementData& Edata1);
		static void mergeElements(ElementData& Edata0, const ElementData& Edata1, const bool& join);

		double compSnowpackInternalEnergyChange(const double sn_dt);
		double getModelledTemperature(const double& z) const;

		unsigned int getNumberOfElements() const;
		unsigned int getNumberOfNodes() const;

		mio::StationData meta;      ///< Station meta data
		double Albedo;              ///< Snow albedo
		double SoilAlb;             ///< Soil albedo
		double BareSoil_z0;         ///< Bare soil roughness in m
		unsigned int SoilNode;      ///< The top soil node, 0 in case of SNP_SOIL == 0
		double cH;                  ///< The CALCULATED snowpack height, including soil depth if SNP_SOIL == 1
		double mH;                  ///< The MEASURED snowpack height, including soil depth if SNP_SOIL == 1
		double Ground;              ///< The ground height -- meaning the height of the top soil node
		double hn;                  ///< Depth of new snow to be used on slopes
		double rho_hn;              ///< Density of new snow to be used on slopes
		bool windward;              ///< True for windward (luv) slope
		unsigned int ErosionLevel;  ///< Element where snow erosion stopped previously for the drift index
		double ErosionMass;         ///< Eroded mass either real or virtually (storage if less than one element)
		int S_class1;               ///< Stability class based on hand hardness, grain class ...
		int S_class2;               ///< Stability class based on hand hardness, grain class ...
		double S_d;                 ///< Minimum Direct Action Stability Index  ...
		double z_S_d;               ///< Depth of Minimum Direct Action Stability
		double S_n;                 ///< Minimum Natural Stability Index
		double z_S_n;               ///< Depth of Minimum Natural Stability
		double S_s;                 ///< Minimum Alternative Stability Index (ASI; skier stability)
		double z_S_s;               ///< Depth of Minimum ASI
		double S_4;                 ///< stab_index4
		double z_S_4;               ///< Depth of stab_index4
		double S_5;                 ///< stab_index5
		double z_S_5;               ///< Depth of stab_index5
		std::vector<NodeData> Ndata;    ///< pointer to nodal data array (e.g. T, z, u, etc..)
		std::vector<ElementData> Edata; ///< pointer to element data array (e.g. Te, L, Rho, etc..)
		void *Kt, *Ks;              ///< Pointer to pseudo-conductivity and stiffnes matrix
		double ColdContent;         ///< Cold content of snowpack (J m-2)
		char SubSurfaceMelt;        ///< Subsurface melting flag ( yes/no ) for exposition
		char SubSurfaceFrze;        ///< Subsurface refreezing flag ( yes/no ) for exposition
		CanopyData Cdata;           ///< Pointer to canopy data
		int tag_low;                ///< Lowest tag to dump, 0 means no tags at all

		static const double join_thresh_l, join_thresh_ice, join_thresh_water;
		static const double join_thresh_dd, join_thresh_sp, join_thresh_rg;
		static const unsigned int number_top_elements;
		static unsigned int number_of_solutes;  ///< The model treats that number of solutes

	private:
		bool useCanopyModel, useSoilLayers; ///< The model includes soil layers
		unsigned int nNodes;                ///< Actual number of nodes; different for each exposition
		unsigned int nElems;                ///< Actual number of elements (nElems=nNodes-1)
};

/// @brief Defines structure for snow profile layers
class SnowProfileLayer {

	public:
		void average(const double& w1, const double& w2, const SnowProfileLayer& Pdata);

		std::string stationname;
		int  loc_for_snow;
		int  loc_for_wind;

		// Version used, date, user, ...
		std::string sn_version;
		std::string sn_computation_date;
		double sn_jul_computation_date;
		std::string sn_user;

		mio::Date profileDate; ///< Date of profile

		mio::Date layerDate; ///< Date of deposition
		double height;       ///< 0 to 1000      (cm)
		double rho;          ///< 0 to 1000      (kg m-3)
		double T;            ///< -50 to 50, snow temperature at top of layer (degC)
		double gradT   ;     ///< -1000 to 1000, temperature gradient across layer (K m-1)
		double strain_rate;  ///< 0 to 1e-5      (s-1)
		double theta_w;      ///< 0 to 100       (%)
		double theta_i;      ///< 0 to 100       (%)
		double dendricity;   ///< 0 to 1         (-)
		double sphericity;   ///< 0 to 1         (-)
		double coordin_num;  ///< 0 to 10        (-)
		double grain_size;   ///< 0 to 100       (mm)
		double bond_size;    ///< 0 to 100       (mm)
		double hard;         ///< 0. to 5.       (-)
		int    marker;       ///< 0 to 100       (-)
		int    type;         ///< 0 to 100       (-)
};


/// Structure of double values for output to SDB
struct ProcessDat {
	mio::Date date; ///< Process date
	char stat_abbrev[16];
	int  loc_for_snow;
	int  loc_for_wind;
	// Version used, date, user, ...
	char sn_version[MAX_STRING_LENGTH];          ///< SNOWPACK version
	char sn_computation_date[MAX_STRING_LENGTH]; ///< Date of computation
	double sn_jul_computation_date;
	char sn_compilation_date[MAX_STRING_LENGTH]; ///< Date of compilation
	char sn_user[MAX_STRING_LENGTH];             ///< SNOWPACK user
	int nHz;               ///< Number of hazard steps

	// Snow depth
	double ch;             ///< cm
	// SWE, total liquid water content and runoff
	double swe;
	double tot_lwc;
	double runoff;
	// Surface hoar index
	double dewpt_def;      ///< degC
	double hoar_size;      ///< mm
	double hoar_ind6;      ///< kg m-2
	double hoar_ind24;     ///< kg m-2
	// Drift index
	double wind_trans;     ///< cm
	double wind_trans24;   ///< cm
	// New snow depths
	double hn_half_hour;   ///< cm
	double hn3;            ///< cm
	double hn6;            ///< cm
	double hn12;           ///< cm
	double hn24;           ///< cm
	double hn72;           ///< cm
	double hn72_24;        ///< cm;
	// New snow water equivalents
	double hnw_half_hour;  ///< kg m-2
	double hnw3;           ///< kg m-2
	double hnw6;           ///< kg m-2
	double hnw12;          ///< kg m-2
	double hnw24;          ///< kg m-2
	double hnw72;          ///< kg m-2
	// Stability indices
	int stab_class1;       ///< Stability class 1,3,5
	int stab_class2;       ///< Profile type 0..10
	double stab_index1;
	double stab_index2;
	double stab_index3;
	double stab_index4;
	double stab_index5;
	double stab_height1;   ///< cm
	double stab_height2;   ///< cm
	double stab_height3;   ///< cm
	double stab_height4;   ///< cm
	double stab_height5;   ///< cm
	// Special parameters
	double crust;          ///< cm
	double en_bal;         ///< kJ m-2
	double sw_net;         ///< kJ m-2
	double t_top1;         ///< degC
	double t_top2;         ///< degC
};

struct ProcessInd {
	short stat_abbrev;
	short loc_for_snow;
	short loc_for_wind;
	// Data
	short ch;

	short swe;
	short tot_lwc;
	short runoff;

	short dewpt_def;
	short hoar_size;
	short hoar_ind6;
	short hoar_ind24;

	short wind_trans;
	short wind_trans24;

	short hn3;
	short hn6;
	short hn12;
	short hn24;
	short hn72;
	short hn72_24;

	short hnw3;
	short hnw6;
	short hnw12;
	short hnw24;
	short hnw72;

	short stab_class1;
	short stab_class2;
	short stab_index1;
	short stab_index2;
	short stab_index3;
	short stab_index4;
	short stab_index5;
	short stab_height1;
	short stab_height2;
	short stab_height3;
	short stab_height4;
	short stab_height5;

	short crust;
	short en_bal;
	short sw_net;
	short t_top1;
	short t_top2;
};

#endif
