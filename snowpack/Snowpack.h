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

#ifndef __SNOWPACK_H__
#define __SNOWPACK_H__

#include <snowpack/DataClasses.h>
#include <snowpack/Constants.h>
#include <snowpack/Solver.h>
#include <snowpack/Hazard.h>
#include <snowpack/Utils.h>
#include <snowpack/Laws_sn.h>
#include <snowpack/Laws.h>
#include <snowpack/Metamorphism.h>
#include <snowpack/Aggregate.h>
#include <snowpack/PhaseChange.h>
#include <snowpack/SnowDrift.h>
#include <snowpack/Stability.h>
#include <snowpack/WaterTransport.h>

#include <meteoio/MeteoIO.h>
#include <string>
#include <sstream>
#include <errno.h>


/**
 * @name Element macro definitions used in Snowpack.cc only
 * TODO to be replaced by proper functions some day
 */
//@{
/// @brief The number of element incidences
#define N_OF_INCIDENCES 2
/// @brief Defines the assembly macro
#define EL_INCID(e, Ie) { \
	Ie[0] = e; Ie[1] = e+1; }
/// @brief Define the node to element temperature macro
#define EL_TEMP( Ie, Te0, Tei, T0, Ti ) { \
	Te0[ 0 ] = T0[ Ie[ 0 ] ].T; \
	Te0[ 1 ] = T0[ Ie[ 1 ] ].T; \
	Tei[ 0 ] = Ti[ Ie[ 0 ] ]; \
	Tei[ 1 ] = Ti[ Ie[ 1 ] ]; }
/// @brief Define the node to element coordinates macro
#define EL_COORD( Ze, Ie, NODES ) { \
	Ze[ 0 ] = NODES[ ( Ie[ 0 ] ) ].z; \
	Ze[ 1 ] = NODES[ ( Ie[ 1 ] ) ].z;  }
/// @brief Define the node to element coordinates macro
#define EL_U( Ue, Ie, NODES ) { \
	Ue[ 0 ] = NODES[ ( Ie[ 0 ] ) ].u; \
	Ue[ 1 ] = NODES[ ( Ie[ 1 ] ) ].u;  }
/// @brief Element right-hand side macro
#define EL_RGT_ASSEM(F, Ie, Fe) { \
	F[Ie[0]] += Fe[0]; F[Ie[1]] += Fe[1]; }
//@}

class Snowpack {

	public:
		///@brief New snow density models
		enum NewSnowDensityModel {
			ZWART,       ///< Costijn Zwart's model (elaborated 2006; in use since 4 Dec 2007
			LEHNING_OLD, ///< First model by M. Lehning
			LEHNING_NEW, ///< Improved model by M. Lehning, incl. ad-hoc wind & temperature effects (used until 06/07)
			BELLAIRE,    ///< Sascha Bellaire's model (elaborated 2007; used summer/fall 2007)
			PAHAUT,      ///< Edmond Pahaut's model, introduced Sep 1995 in CROCUS by G. Giraud
			EVENT,       ///< Introduced by Christine Groot Zwaaftink 2009 (Antarctica)
			MEASURED,    ///< Use measured new snow density read from meteo input
			FIXED,       ///< Fixed new snow density => FIXED_HN_DENSITY (SnowMIP NSD experiment, Antarctica)
			N_HNDM
		};

		Snowpack(const mio::Config& i_cfg);

		void runSnowpackModel(SN_MET_DATA& Mdata, SN_STATION_DATA& Xdata, double& cumu_hnw, 
						  SN_BOUNDARY_DATA& Bdata, SN_SURFACE_DATA& Sdata);

		static double sn_NewSnowDensityHendrikx(const double ta, const double tss, const double rh, const double vw);

		static double calculateNewSnowDensity(const SN_MET_DATA& Mdata, const SN_STATION_DATA& Xdata, 
									   const double& tss, const double& hnw, const NewSnowDensityModel& model);

		const static double new_snow_albedo;
		static NewSnowDensityModel hn_density_model; ///<New snow density model to be used

	private:
		/**
		 * @brief Specifies what kind of boundary condition is to be implemented at the top surface \n
		 * 0 : use surface fluxes (NEUMANN_BC)\n
		 * 1 : use prescribed surface temperature (DIRICHLET_BC)
		 */
		enum BoundaryCondition {
			NEUMANN_BC,
			DIRICHLET_BC
		};

		///@brief Types of events for computing new snow density
		enum EventType {
			EVENT_DEFAULT,  ///< Default
			EVENT_WIND,     ///< Wind driven deposition of snow
			N_EVENT_TYPE
		};

		static double sn_NewSnowDensityPara(double TA, double TSS, double RH, double VW, double HH);

		static double sn_NewSnowDensityEvent(const SN_MET_DATA& Mdata, const EventType& i_event_type);

		bool sn_SnowForces(SN_ELEM_DATA *Edata,  double dt, double cos_sl, double Zn[ N_OF_INCIDENCES ], 
					    double Un[ N_OF_INCIDENCES ], double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ], 
					    double Fc[ N_OF_INCIDENCES ], double Fi[ N_OF_INCIDENCES ], double Fe[ N_OF_INCIDENCES ]);

		void calcSnowCreep(const SN_MET_DATA& Mdata, SN_STATION_DATA& Xdata);
	
		bool sn_ElementKtMatrix(SN_ELEM_DATA *Edata, double dt, double dvdz, double T0[ N_OF_INCIDENCES ], 
						    double Tn[ N_OF_INCIDENCES ], double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ], 
						    double Fe[ N_OF_INCIDENCES ], char *SubSurfaceMelt, char *SubSurfaceFrze, 
						    double VaporEnhance);
		
		void updateMeteoHeatFluxes(const SN_MET_DATA& Mdata, const SN_STATION_DATA& Xdata, 
							  SN_BOUNDARY_DATA& Bdata, SN_SURFACE_DATA& Sdata);

		void sn_Neumann(const SN_MET_DATA& Mdata, SN_BOUNDARY_DATA& Bdata, const SN_STATION_DATA& Xdata, 
					const double& T_snow, const double& T_iter, 
					double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ], double Fe[ N_OF_INCIDENCES ]);

		void sn_NeumannSoil(const double& flux, const double& T_snow, 
						double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ], double Fe[ N_OF_INCIDENCES ]);

		void sn_SnowTemperature(SN_STATION_DATA& Xdata, SN_MET_DATA& Mdata, SN_BOUNDARY_DATA& Bdata, double& mAlb);
	
		void assignSomeFluxes(const SN_STATION_DATA& Xdata, const SN_MET_DATA& Mdata, const double& mAlb, 
						  SN_SURFACE_DATA& Sdata);
 
		void determineSnowFall(const SN_MET_DATA& Mdata, SN_STATION_DATA& Xdata, double& cumu_hnw);
		
		mio::Config cfg;
		BoundaryCondition surfaceCode;
		bool research_mode, useCanopyModel, enforce_measured_snow_heights, soil_flux, useSnowLayers;
		int sw_ref;
		double thresh_change_bc, geo_heat, height_of_meteo_values, meteo_step_length, hns_ne_height, thresh_rain, sn_dt;
		double t_crazy_min, t_crazy_max, thresh_rh;
		std::string viscosity_model, variant;
		bool multistream, join_elements, change_bc, meas_tss;
		double new_snow_dd, new_snow_sp, new_snow_dd_wind, new_snow_sp_wind, rh_lowlim, bond_factor_rh;
		double new_snow_grain_rad, new_snow_bond_rad;
		double density_hoar_surf, density_hoar_buried, min_size_hoar_buried;
		double minimum_l_element;
		bool vw_dendricity;

		//The following block is necessary, because calculateNewSnowDensity() is called from other modules
		//and shall remain static, therefore these variables also have to be static
		static void initStaticData(const std::string& variant);
		static double fixed_hn_density, max_hn_density, event_wind_lowlim, event_wind_highlim;
		static EventType event_type;

		const static bool jordy_new_snow, hydrometeor;
		const static double snowfall_warning;
		const static int new_snow_marker;
}; //end class Snowpack

#endif //End of Snowpack.h
