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
 * @todo to be replaced by proper functions some day
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
		Snowpack(const mio::Config& i_cfg);

		void runSnowpackModel(CurrentMeteo& Mdata, SnowStation& Xdata, double& cumu_hnw,
		                      BoundCond& Bdata, SurfaceFluxes& Sdata);

		double getThreshRain() const;

		const static double new_snow_albedo, min_ice_content;

	private:
		/**
		 * @brief Specifies what kind of boundary condition is to be implemented at the top surface \n
		 * - 0 : use surface fluxes (NEUMANN_BC)
		 * - 1 : use prescribed surface temperature (DIRICHLET_BC)
		 */
		enum BoundaryCondition {
			NEUMANN_BC,
			DIRICHLET_BC
		};

		bool compSnowForces(ElementData *Edata,  double dt, double cos_sl, double Zn[ N_OF_INCIDENCES ],
		                    double Un[ N_OF_INCIDENCES ], double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ],
		                    double Fc[ N_OF_INCIDENCES ], double Fi[ N_OF_INCIDENCES ],
		                    double Fe[ N_OF_INCIDENCES ]);

		void compSnowCreep(const CurrentMeteo& Mdata, SnowStation& Xdata);

		bool sn_ElementKtMatrix(ElementData *Edata, double dt, double dvdz, double T0[ N_OF_INCIDENCES ],
		                        double Tn[ N_OF_INCIDENCES ], double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ],
		                        double Fe[ N_OF_INCIDENCES ], char *SubSurfaceMelt, char *SubSurfaceFrze,
		                        double VaporEnhance);

		void updateMeteoHeatFluxes(const CurrentMeteo& Mdata, SnowStation& Xdata,
                                   BoundCond& Bdata);

		void neumannBoundaryConditions(const CurrentMeteo& Mdata, BoundCond& Bdata, const SnowStation& Xdata,
		                               const double& T_snow, const double& T_iter,
		                               double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ],
		                               double Fe[ N_OF_INCIDENCES ]);

		void neumannBoundaryConditionsSoil(const double& flux, const double& T_snow,
		                                   double Se[ N_OF_INCIDENCES ][ N_OF_INCIDENCES ],
		                                   double Fe[ N_OF_INCIDENCES ]);

		void compSnowTemperatures(SnowStation& Xdata, CurrentMeteo& Mdata, BoundCond& Bdata);

		void assignSomeFluxes(SnowStation& Xdata, const CurrentMeteo& Mdata, const double& mAlb,
		                      SurfaceFluxes& Sdata);

		void compSnowFall(const CurrentMeteo& Mdata, SnowStation& Xdata, double& cumu_hnw,
		                  SurfaceFluxes& Sdata);

		std::string hn_density, hn_density_model, viscosity_model, variant;
		const mio::Config& cfg;
		BoundaryCondition surfaceCode;
		bool research_mode, useCanopyModel, enforce_measured_snow_heights, soil_flux, useSoilLayers;
		bool multistream, join_elements, change_bc, meas_tss;
		bool vw_dendricity;
		int sw_mode;
		double thresh_change_bc, geo_heat, height_of_meteo_values, height_new_elem, thresh_rain, sn_dt;
		double t_crazy_min, t_crazy_max, thresh_rh;
		double new_snow_dd, new_snow_sp, new_snow_dd_wind, new_snow_sp_wind, rh_lowlim, bond_factor_rh;
		double new_snow_grain_rad, new_snow_bond_rad;
		double hoar_density_buried, hoar_density_surf, hoar_min_size_buried;
		double minimum_l_element;
		double fixed_albedo;
		double t_surf;

		const static bool hydrometeor;
		const static double snowfall_warning;
		const static unsigned int new_snow_marker;

		double meteo_step_length;
}; //end class Snowpack

#endif
