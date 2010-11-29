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

#ifndef __AGGREGATE_H__
#define __AGGREGATE_H__

#include <snowpack/DataClasses.h>
#include <snowpack/Constants.h>
#include <snowpack/Utils.h>
#include <snowpack/Laws.h>
#include <vector>

/**
 * @class Aggregate
 * @version 7.03
 * @bug     -
 * @brief This module contains the routines to perform profile aggregation
 */
class Aggregate {

	public:
		static int aggregate(std::vector<Q_PROFILE_DAT>& Pdata);

	private:
		static void shift(const int& nE, std::vector<Q_PROFILE_DAT>& Pdata);
		static bool doJoin2(const int& e1, std::vector<Q_PROFILE_DAT>& Pdata);
		static bool doJoin(const int& e1, std::vector<Q_PROFILE_DAT>& Pdata);

		static const double limit_dry;     ///< Distinguishes between dry and wet snow layers (1)
		static const double diff_theta_w;  ///< Maximum water difference for aggregation (Vol %)
		static const double diff_jul;      ///< Maximum  age difference (days) for aggregation
		static const double diff_sp;       ///< Maximum  sphericity difference for aggregation
		static const double diff_dd;       ///< Maximum  dendricity difference for aggregation
		static const double diff_dg;       ///< Maximum  grain size difference for aggregation
		static const double diff_dg_rel;   ///< Maximum  relative grain size difference for aggregation
		static const double min_l_element; ///< Minimum length of element to be kept separate
};

#endif
