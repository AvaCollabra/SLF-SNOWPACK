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
#ifndef __PHASE_CHANGE_H__
#define __PHASE_CHANGE_H__

#include <snowpack/Constants.h>
#include <snowpack/DataClasses.h>
#include <snowpack/Hazard.h>
#include <snowpack/Laws_sn.h>
#include <meteoio/MeteoIO.h>

/**
 * @class PhaseChange
 * @author Perry Bartelt \n Michael Lehning and others
 * @version 10.03
 * @brief This class contains the phase change routines for the 1d snowpack model
 *        It also updates the volumetric contents of each element.
 */
class PhaseChange {
 	public:
		PhaseChange(const mio::Config& i_cfg);

		void runPhaseChange(const SurfaceFluxes& Sdata, SnowStation& Xdata);

	private:
		void compSubSurfaceMelt(ElementData& Edata, const double& dt, double& ql_Rest);
		void compSubSurfaceFrze(ElementData& Edata, const double& dt);
		
		double sn_dt;
		const mio::Config& cfg;

		static const double theta_r; ///< Residual Water Content,  for now we say  0.0
		static const double theta_s; ///< Saturated Water Content, for now we say  1.0
};

#endif
