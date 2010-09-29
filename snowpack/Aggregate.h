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
 * @file Aggregate.h
 * @version 9.10
 */

#ifndef __AGGREGATE_H__
#define __AGGREGATE_H__

#include <snowpack/Constants.h>
#include <snowpack/Snowpack.h>
#include <snowpack/Hazard.h>
#include <snowpack/Utils.h>
#include <snowpack/Laws.h>

struct Q_PROFILE_DAT;

/*
 * CONSTANTS
 */
/**
 * @name Thresholds for aggregation
 */
//@{
/// Distinguishes between dry and wet snow layers (1)
#define LIMIT_DRY 0.001
/// Maximum water difference for aggregation (Vol %)
#define DIFF_THETA_W 0.7
/// Maximum  age difference (days) for aggregation
#define DIFF_JUL 1.0
/// Maximum  sphericity difference for aggregation
#define DIFF_SP 0.2
/// Maximum  dendricity difference for aggregation
#define DIFF_DD 0.2
/// Maximum  grain size difference for aggregation
#define DIFF_DG 0.25
/// Maximum  relative grain size difference for aggregation
#define DIFF_DG_REL 0.2
/// Minimum length of element to be kept separate
#define MIN_L_ELEMENT 0.3
//@}

/*
 * FUNCTION PROTOTYPES
 */
void ml_ag_Average(const int e, const int l, const double w1, const double w2, Q_PROFILE_DAT *Pdata);

void ml_ag_Shift(const int nE, Q_PROFILE_DAT *Pdata);

int ml_ag_Join2(const int e1, Q_PROFILE_DAT *Pdata);

int ml_ag_Join(const int e1, Q_PROFILE_DAT *Pdata);

int ml_ag_Classify(const double dendricity, const double sphericity, const double grain_dia, const int marker, const double theta_w, const double theta_i);

int ml_ag_Aggregate(int nE, Q_PROFILE_DAT *Pdata);

#endif //End of Template.h
