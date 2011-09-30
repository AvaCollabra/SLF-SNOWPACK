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
 * @param Pdata A vector of SnowProfileLayer
 */
void Aggregate::shift(const unsigned int& nE, std::vector<SnowProfileLayer>& Pdata)
{
	for (unsigned int e=1, ll=1; e<nE; e++) {
		if (Pdata[e].height != Constants::undefined) {
			if (ll != e) {
				Pdata[ll].height = Pdata[e].height;
				Pdata[ll].layerDate = Pdata[e].layerDate;
				Pdata[ll].rho = Pdata[e].rho;
				Pdata[ll].T = Pdata[e].T;
				Pdata[ll].gradT = Pdata[e].gradT;
				Pdata[ll].strain_rate = Pdata[e].strain_rate;
				Pdata[ll].theta_w = Pdata[e].theta_w;
				Pdata[ll].theta_i = Pdata[e].theta_i;
				Pdata[ll].dendricity = Pdata[e].dendricity;
				Pdata[ll].sphericity = Pdata[e].sphericity;
				Pdata[ll].coordin_num = Pdata[e].coordin_num;
				Pdata[ll].grain_size = Pdata[e].grain_size;
				Pdata[ll].bond_size = Pdata[e].bond_size;
				Pdata[ll].hard = Pdata[e].hard;
				Pdata[ll].marker = Pdata[e].marker;
			}
			ll++;
		}
	}
}

/**
 * @brief Decide whether to codense two layers
 * @param e1 (int)
 * @param Pdata A vector of SnowProfileLayer
 */
bool Aggregate::doJoin2(const unsigned int& e1, std::vector<SnowProfileLayer>& Pdata)
{
	unsigned int e2 = e1+1;

	// if a dry layer is involved
	if ((Pdata[e1].theta_w < limit_dry) || (Pdata[e2].theta_w < limit_dry)){
		// do not combine dry with moist layers
		if (fabs(Pdata[e1].theta_w - Pdata[e2].theta_w) > limit_dry)
			return false;

		// do not combine layers which are of quite different age
		if (fabs(Pdata[e1].layerDate.getJulianDate() - Pdata[e2].layerDate.getJulianDate()) > 2 * diff_jul)
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
}

/**
 * @brief Decide whether to codense two layers
 * @param e1 (int)
 * @param Pdata A vector of SnowProfileLayer
 */
bool Aggregate::doJoin(const unsigned int& e1, std::vector<SnowProfileLayer>& Pdata)
{
	unsigned int e2 = e1 - 1;

	if ((Pdata[e1].theta_w < limit_dry) || (Pdata[e2].theta_w < limit_dry)) {
		if (fabs(Pdata[e1].theta_w - Pdata[e2].theta_w) > limit_dry)
			return false;
	} else {
		if (fabs(Pdata[e1].theta_w - Pdata[e2].theta_w) > diff_theta_w)
			return false;
	}

	if (fabs(Pdata[e1].layerDate.getJulianDate() - Pdata[e2].layerDate.getJulianDate()) > diff_jul)
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
}

/**
 * @brief Aggregate snow profile layers and compute the grain class
 * @param Pdata A vector of SnowProfileLayer
 */
unsigned int Aggregate::aggregate(std::vector<SnowProfileLayer>& Pdata)
{
	bool flag = false;
	unsigned int nL, nE = (unsigned)Pdata.size();
	double Lp0, Lp1; // lower and upper layer

	// Initialize number of layers and aggregate only if more than 5 layers
	nL = nE;
	if (nL > 5) {
		// First Run - do not touch top element
		// keep track of the coordinates and length of elements
		Lp0 = (Pdata[nE-2].height -  Pdata[nE-3].height);
		for (unsigned int l1=nE-2; l1 > 0; l1--) {
			int l0 = l1-1;
			Lp1 = Lp0;
			if (l0 > 0) {
				Lp0 = (Pdata[l0].height -  Pdata[l0-1].height);
			} else {
				Lp0 = Pdata[l0].height;
			}
			// if two layers are similar combine them
			if (doJoin(l1, Pdata) && (Pdata[l0].marker != 3) && (Pdata[l1].marker != 3)) {
				nL--;
				Pdata[l0].average(Lp0, Lp1, Pdata[l1]);
				Pdata[l1].height = Constants::undefined;
				Lp0 += Lp1;
			}
		}  // for all elements

		shift(nE, Pdata);
		nE = nL;

		// Second Run - aggregate remaining very thin layers
		if (nE > 2) {
			Lp0 = (Pdata[nE-2].height -  Pdata[nE-3].height);
			for(unsigned int l1 = nE-2; l1 > 0; l1--) {
				int l0 = signed(l1)-1;
				Lp1 = Lp0;
				if (l0 > 0) {
					Lp0 = (Pdata[l0].height -  Pdata[l0-1].height);
				} else {
					Lp0 = Pdata[l0].height;
				}
				if ((Pdata[l0].marker != 3) && (Pdata[l1].marker != 3)) {
					// trick to try to join with upper or lower level -> use flag to mark thin layer
					if (flag || (Lp0 < (sqrt(Pdata[nE-1].height-Pdata[l0].height)/4.)) || (Lp0 < min_l_element)) {
						// if two layers are similar or one layer is very very small combine them
						if (doJoin2(l1, Pdata) || (Lp0 < min_l_element) || (Lp1 < min_l_element)){
							nL--;
							Pdata[l0].average(Lp0, Lp1, Pdata[l1]);
							Pdata[l1].height = Constants::undefined;
							Lp0 += Lp1;
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
	for(unsigned int ll=0; ll<nL; ll++) {
		Pdata[ll].type = ElementData::snowType(Pdata[ll].dendricity, Pdata[ll].sphericity,
                                           Pdata[ll].grain_size, Pdata[ll].marker, Pdata[ll].theta_w/100.,
                                           ElementData::snowResidualWaterContent(Pdata[ll].theta_i/100.));
	}
	return (nL);
}
