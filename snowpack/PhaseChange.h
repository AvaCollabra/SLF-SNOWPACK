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
		void initialize(SnowStation& Xdata);								//Call before first call to compPhaseChange in a time step
		void finalize(const SurfaceFluxes& Sdata, SnowStation& Xdata, const mio::Date& date_in);	//Call after last call to compPhaseChange in a time step
		void compPhaseChange(const SurfaceFluxes& Sdata, SnowStation& Xdata, const mio::Date& date_in, const bool& verbose=true);	//Call to do a phase change in a time step

		static const double theta_r; ///< Residual Water Content,  for now we say  0.0

	private:
		void compSubSurfaceMelt(ElementData& Edata, const unsigned int nSolutes, const double& dt,
		                        double& ql_Rest, const mio::Date& date_in);
		void compSubSurfaceFrze(ElementData& Edata, const unsigned int nSolutes, const double& dt,
		                        const mio::Date& date_in);

		double sn_dt; ///< The calculation_step_length in seconds

		double cold_content_in;		///< cold content before first PhaseChange call (for checking energy balance)
		double cold_content_out;	///< cold content after last PhaseChange call (for checking energy balance)

		static const double theta_s; ///< Saturated Water Content, for now we say  1.0
};

#endif
