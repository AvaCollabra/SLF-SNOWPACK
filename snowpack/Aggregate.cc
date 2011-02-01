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

#include <snowpack/Aggregate.h>

/************************************************************
 * static section                                           *
 ************************************************************/

const double Aggregate::limit_dry     = 0.001; ///< Distinguishes between dry and wet snow layers (1)
const double Aggregate::diff_theta_w  = 0.7;   ///< Maximum water difference for aggregation (Vol %)
const double Aggregate::diff_jul      = 1.0;   ///< Maximum  age difference (days) for aggregation
const double Aggregate::diff_sp       = 0.2;   ///< Maximum  sphericity difference for aggregation
const double Aggregate::diff_dd       = 0.2;   ///< Maximum  dendricity difference for aggregation
const double Aggregate::diff_dg       = 0.25;  ///< Maximum  grain size difference for aggregation
const double Aggregate::diff_dg_rel   = 0.2;   ///< Maximum  relative grain size difference for aggregation
const double Aggregate::min_l_element = 0.3;   ///< Minimum length of element to be kept separate

/**
 * @brief Eliminate the "empty" layers shift the remaining layers to form a compact snowpack
 * @param nE (int)
 * @param Pdata A vector of Q_PROFILE_DAT
 */
void Aggregate::shift(const int& nE, std::vector<Q_PROFILE_DAT>& Pdata)
{
	int l = 0;
	for (int e=1; e<nE; e++) {
		if (Pdata[e].height != Constants::undefined) {
			l++;
			Pdata[l].height = Pdata[e].height;
			Pdata[l].layer_date = Pdata[e].layer_date;
			Pdata[l].rho = Pdata[e].rho;
			Pdata[l].tem = Pdata[e].tem;
			Pdata[l].tem_grad = Pdata[e].tem_grad;
			Pdata[l].strain_rate = Pdata[e].strain_rate;
			Pdata[l].theta_w = Pdata[e].theta_w;
			Pdata[l].theta_i = Pdata[e].theta_i;
			Pdata[l].dendricity = Pdata[e].dendricity;
			Pdata[l].sphericity = Pdata[e].sphericity;
			Pdata[l].coordin_num = Pdata[e].coordin_num;
			Pdata[l].grain_size = Pdata[e].grain_size;
			Pdata[l].bond_size = Pdata[e].bond_size;
			Pdata[l].hard = Pdata[e].hard;
			Pdata[l].marker = Pdata[e].marker;
		}
	} // for all old elements
} // End of shift


/**
 * @brief Decide whether to codense two layers
 * @param e1 (int)
 * @param Pdata A vector of Q_PROFILE_DAT
 */
bool Aggregate::doJoin2(const int& e1, std::vector<Q_PROFILE_DAT>& Pdata)
{
	int e2 = e1+1;

	// if a dry layer is involved
	if ((Pdata[e1].theta_w < limit_dry) || (Pdata[e2].theta_w < limit_dry)){
		// do not combine dry with moist layers
		if (fabs(Pdata[e1].theta_w - Pdata[e2].theta_w) > limit_dry)
			return false;

		// do not combine layers which are of quite different age
		if (fabs(Pdata[e1].layer_date - Pdata[e2].layer_date) > 2 * diff_jul)
			return false;

		// do not combine layers with different grain classes
		if (Pdata[e1].marker != Pdata[e2].marker) {
			return false;
		}
		// do not combine layers which are of quite different hardness
		if (fabs(Pdata[e1].hard - Pdata[e2].hard) > 1.5) {
			return false;
		}
	}
	// for two wet layers
	else {
		if ((Pdata[e1].grain_size < 0.75) || (Pdata[e1].sphericity < 0.5) ||
			(Pdata[e2].grain_size < 0.75) || (Pdata[e2].sphericity < 0.5)) {
			if (Pdata[e1].marker != Pdata[e2].marker) {
				return false;
			}
		}
		else {
			if (!((((Pdata[e1].marker > 9) &&  (Pdata[e1].marker < 13)) ||
				((Pdata[e1].marker > 19) &&  (Pdata[e1].marker < 23)) ) &&
				( ((Pdata[e1].marker > 9) &&  (Pdata[e1].marker < 13)) ||
				((Pdata[e1].marker > 19) &&  (Pdata[e1].marker < 23))))) {
				
				if (Pdata[e1].marker != Pdata[e2].marker)
					return false;
			}
		}
	}
	return true;
} // End of doJoin2


/**
 * @brief Decide whether to codense two layers
 * @param e1 (int)
 * @param Pdata A vector of Q_PROFILE_DAT
 */
bool Aggregate::doJoin(const int& e1, std::vector<Q_PROFILE_DAT>& Pdata)
{
	int e2 = e1 - 1;

	if ((Pdata[e1].theta_w < limit_dry) || (Pdata[e2].theta_w < limit_dry)) {
		if (fabs(Pdata[e1].theta_w - Pdata[e2].theta_w) > limit_dry)
			return false;
	} else {
		if (fabs(Pdata[e1].theta_w - Pdata[e2].theta_w) > diff_theta_w)
			return false;
	}

	if (fabs(Pdata[e1].layer_date - Pdata[e2].layer_date) > diff_jul)
		return false;

	if (Pdata[e1].marker != Pdata[e2].marker)
		return false;

	// do not combine layers which are of quite different hardness 020917;Fz
	if (fabs(Pdata[e1].hard - Pdata[e2].hard) > 1.0)
		return false;

	if ((Pdata[e1].dendricity == 0) && (Pdata[e2].dendricity == 0)){
		if ( fabs(Pdata[e1].sphericity - Pdata[e2].sphericity) > diff_sp)
			return false;
		
		if (fabs(Pdata[e1].grain_size - Pdata[e2].grain_size) >	MAX(diff_dg, diff_dg_rel * Pdata[e1].grain_size))
			return false;
	} else {
		if (fabs(Pdata[e1].sphericity - Pdata[e2].sphericity) > diff_sp)
			return false;

		if (fabs(Pdata[e1].dendricity - Pdata[e2].dendricity) > diff_dd)
			return false;
	}

	return true;
} // End of doJoin


/**
 * @brief Aggregate layers and compute the grain class
 * @param Pdata A vector of Q_PROFILE_DAT
 */
int Aggregate::aggregate(std::vector<Q_PROFILE_DAT>& Pdata)
{
	bool flag = false;
	int e, nL, l;
	double l1, l2;

	int nE = (int)Pdata.size();

	// Initialize number of elements  and aggregate only if more than 5 layers
	nL = nE;
	if (nL > 5) {
		// First Run -  leave top element alone
		// keep track of the coordinates and length of elements
		l2 = (Pdata[nE-2].height -  Pdata[nE-3].height);
		for (e=nE-2; e>0; e--) {
			l = e-1;
			l1 = l2;
			if (l>0) {
				l2 = (Pdata[l].height -  Pdata[l-1].height);
			} else {
				l2 = Pdata[l].height;
			}
			// if two layers are similar combine them
			if (doJoin(e, Pdata) && (Pdata[l].marker != 3) && (Pdata[e].marker != 3)) {
				nL--;
				Pdata[l].average(l1, l2, Pdata[e]);
				Pdata[e].height = Constants::undefined;
				l2 += l1;
			}
		}  // for all elements

		shift(nE, Pdata);
		nE = nL;

		// Second Run - aggregate remaining very thin layers
		if (nE > 2) {
			l2 = (Pdata[nE-2].height -  Pdata[nE-3].height);
			for(e=nE-2; e>0; e--) {
				l = e-1;
				l1 = l2;
				if (l>0) {
					l2 = (Pdata[l].height -  Pdata[l-1].height);
				} else {
					l2 = Pdata[l].height;
				}
				if ((Pdata[l].marker != 3) && (Pdata[e].marker != 3)) {
					// trick to try to join with upper or lower level -> use flag to mark thin layer
					if (flag || (l2 < (sqrt(Pdata[nE-1].height-Pdata[l].height)/4.)) || (l2 < min_l_element)) {
						// if two layers are similar or one layer is very very small combine them
						if (doJoin2(e, Pdata) || (l2 < min_l_element) || (l1 < min_l_element)){
							nL--;
							Pdata[l].average(l1, l2, Pdata[e]);
							Pdata[e].height = Constants::undefined;
							l2 += l1;
							flag = false;
						} else {
							flag = true;
						}
					} else {
						flag = false;
					}
				} // if not surface hoar layer
				else {
					flag = false;
				}
			}  // for all elements
			shift(nE, Pdata);
		} // if nE > 2
	} // if more than 5 layers

	// Update snow type
	for(e=0; e<nL; e++) {
		Pdata[e].type = ElementData::snowType(Pdata[e].dendricity, Pdata[e].sphericity,
                                          Pdata[e].grain_size, Pdata[e].marker, Pdata[e].theta_w/100.,
                                          ElementData::snowResidualWaterContent(Pdata[e].theta_i/100.));
	}
	return (nL);
} // End of aggregate
