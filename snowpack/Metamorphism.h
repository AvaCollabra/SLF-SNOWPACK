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
 * @file Metamorphism.h
 */

#ifndef __METAMORPHISM_H__
#define __METAMORPHISM_H__

#include <snowpack/Constants.h>
#include <snowpack/DataClasses.h>
#include <snowpack/Utils.h>
#include <snowpack/Laws.h>
#include <snowpack/Laws_sn.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <fcntl.h>
#include <meteoio/MeteoIO.h>
#include <map>
#include <string>

class SN_STATION_DATA;
class Metamorphism;

typedef void (Metamorphism::*MetaModelFn)(const SN_MET_DATA&, SN_STATION_DATA&); 
typedef double (Metamorphism::*MetaSpRateFn)(const SN_ELEM_DATA&);

class Metamorphism {

	public:
		Metamorphism(const mio::Config& i_cfg);

		void runMetamorphismModel(const SN_MET_DATA& Mdata, SN_STATION_DATA& Xdata) throw();

		static double csPoreArea(const SN_ELEM_DATA& Edata);

		static double getCoordinationNumberN3(const double& Rho);

		static double ddRate(const SN_ELEM_DATA& Edata);

		static const double mm_tg_dpdz, ba_g_fudge, sa_g_fudge, max_grain_growth, bond_size_stop;
		static const double max_grain_bond_ratio, wind_slab_enhance, wind_slab_vw, wind_slab_depth;
		
	private:
		double TGBondRate(const SN_ELEM_DATA& Edata);
 
		double LatticeConstant0(const double& th_ice);

		double TGGrainRate(const SN_ELEM_DATA& Edata, const double& Tbot, const double& Ttop,
		                   const double& gradTSub, const double& gradTSup);

		double ETBondRate(const SN_ELEM_DATA& Edata);
		double ETGrainRate(const SN_ELEM_DATA& Edata);

		double PressureSintering(const SN_ELEM_DATA& Edata);

		void metamorphismDEFAULT(const SN_MET_DATA& Mdata, SN_STATION_DATA& Xdata);
		void metamorphismNIED(const SN_MET_DATA& Mdata, SN_STATION_DATA& Xdata);
		

		double spRateDEFAULT(const SN_ELEM_DATA& Edata);
		double spRateNIED(const SN_ELEM_DATA& Edata);

	private:
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static containers
		static std::map<std::string, MetaModelFn> mapMetamorphismModel;
		static std::map<std::string, MetaSpRateFn> mapSpRate;
		
		mio::Config cfg;
		std::string metamorphism_model;
		double sn_dt, new_snow_grain_rad;
};

#endif //End of Metamorphism.h
