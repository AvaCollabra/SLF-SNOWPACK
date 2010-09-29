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
 * @file Aggregate.c
 * @version 7.03
 * @bug     -
 * @brief This module contains the routines to perform profile aggregation
 */

#include <snowpack/Aggregate.h>


/**
 * @brief Determine the average quantities of the new layer
 * @param e Index of former layer
 * @param l Index of former/new layer
 * @param w1 Weight of former layer e
 * @param w2 Weight of former layer l
 * @param *Pdata (Q_PROFILE_DAT)
 */
void ml_ag_Average(const int e, const int l, const double w1, const double w2, Q_PROFILE_DAT *Pdata)
{
	Pdata[l].height += w1;
	Pdata[l].layer_date = (w1*Pdata[e].layer_date + w2*Pdata[l].layer_date)/(w1+w2) ;
	Pdata[l].rho = (w1*Pdata[e].rho + w2*Pdata[l].rho)/(w1+w2) ;
	Pdata[l].tem = (w1*Pdata[e].tem + w2*Pdata[l].tem)/(w1+w2) ;
	Pdata[l].tem_grad = (w1*Pdata[e].tem_grad + w2*Pdata[l].tem_grad)/(w1+w2) ;
	Pdata[l].strain_rate = (w1*Pdata[e].strain_rate + w2*Pdata[l].strain_rate)/(w1+w2) ;
	Pdata[l].theta_w = (w1*Pdata[e].theta_w + w2*Pdata[l].theta_w)/(w1+w2) ;
	Pdata[l].theta_i = (w1*Pdata[e].theta_i + w2*Pdata[l].theta_i)/(w1+w2) ;
	Pdata[l].dendricity = (w1*Pdata[e].dendricity + w2*Pdata[l].dendricity)/(w1+w2) ;
	Pdata[l].sphericity = (w1*Pdata[e].sphericity + w2*Pdata[l].sphericity)/(w1+w2) ;
	Pdata[l].coordin_num = (w1*Pdata[e].coordin_num + w2*Pdata[l].coordin_num)/(w1+w2) ;
	Pdata[l].grain_dia = (w1*Pdata[e].grain_dia + w2*Pdata[l].grain_dia)/(w1+w2) ;
	Pdata[l].bond_dia = (w1*Pdata[e].bond_dia + w2*Pdata[l].bond_dia)/(w1+w2) ;
	Pdata[l].hard = (w1*Pdata[e].hard + w2*Pdata[l].hard)/(w1+w2) ;
	Pdata[l].marker = MAX(Pdata[e].marker, Pdata[l].marker) ;

	Pdata[e].height = SNOWPACK_UNDEFINED;

	return;

} // End of ml_ag_Average


/**
 * @brief Eliminate the "empty" layers shift the remaining layers to form a compact snowpack
 * @param nE (int)
 * @param *Pdata (Q_PROFILE_DAT)
 */
void ml_ag_Shift(const int nE, Q_PROFILE_DAT *Pdata)
{
	int e, l=0;
	for (e=1; e<nE; e++) {
		if ( Pdata[e].height != SNOWPACK_UNDEFINED ) {
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
			Pdata[l].grain_dia = Pdata[e].grain_dia;
			Pdata[l].bond_dia = Pdata[e].bond_dia;
			Pdata[l].hard = Pdata[e].hard;
			Pdata[l].marker = Pdata[e].marker;
		}
	} // for all old elements

	return;

} // End of ml_ag_Shift


/**
 * @brief Decide whether to codense two layers
 * @param e1 (int)
 * @param *Pdata (Q_PROFILE_DAT)
 */
int ml_ag_Join2(const int e1, Q_PROFILE_DAT *Pdata)
{
	int e2, aggreg = ON;

	e2 = e1+1;

	// if a dry layer is involved
	if ( (Pdata[e1].theta_w < LIMIT_DRY) || (Pdata[e2].theta_w < LIMIT_DRY) ) {
		// do not combine dry with moist layers
		if ( fabs(Pdata[e1].theta_w - Pdata[e2].theta_w) > LIMIT_DRY ) {
			aggreg = OFF;
		}
		// do not combine layers which are of quite different age
		if ( fabs(Pdata[e1].layer_date - Pdata[e2].layer_date) > 2*DIFF_JUL ) {
			aggreg = OFF;
		}
		// do not combine layers with different grain classes
		if (Pdata[e1].marker != Pdata[e2].marker) {
			aggreg = OFF;
		}
		// do not combine layers which are of quite different hardness
		if ( fabs(Pdata[e1].hard - Pdata[e2].hard) > 1.5 ) {
			aggreg = OFF;
		}
	}
	// for two wet layers
	else {
		if ( (Pdata[e1].grain_dia < 0.75) || (Pdata[e1].sphericity < 0.5) ||
			(Pdata[e2].grain_dia < 0.75) || (Pdata[e2].sphericity < 0.5) ) {
			if (Pdata[e1].marker != Pdata[e2].marker) {
				aggreg = OFF;
			}
		}
		else {
			if ( !( ( ((Pdata[e1].marker > 9) &&  (Pdata[e1].marker < 13)) ||
				((Pdata[e1].marker > 19) &&  (Pdata[e1].marker < 23)) ) &&
				( ((Pdata[e1].marker > 9) &&  (Pdata[e1].marker < 13)) ||
				((Pdata[e1].marker > 19) &&  (Pdata[e1].marker < 23)) ) ) ) {
				if (Pdata[e1].marker != Pdata[e2].marker) {
					aggreg = OFF;
				}
			}
		}
	}

	return(aggreg);
} // End of ml_ag_Join2


/**
 * @brief Decide whether to codense two layers
 * @param e1 (int)
 * @param *Pdata (Q_PROFILE_DAT)
 */
int ml_ag_Join(const int e1, Q_PROFILE_DAT *Pdata)
{
	int e2, aggreg = ON;

	e2 = e1-1;

	if ( (Pdata[e1].theta_w < LIMIT_DRY) || (Pdata[e2].theta_w < LIMIT_DRY) ) {
		 if ( fabs(Pdata[e1].theta_w - Pdata[e2].theta_w) > LIMIT_DRY) {
			 aggreg = OFF;
		 }
	} else {
		if ( fabs(Pdata[e1].theta_w - Pdata[e2].theta_w) > DIFF_THETA_W) {
			aggreg = OFF;
		}
	}

	if ( fabs(Pdata[e1].layer_date - Pdata[e2].layer_date) > DIFF_JUL) {
		aggreg = OFF;
	}

	if (Pdata[e1].marker != Pdata[e2].marker) {
		aggreg = OFF;
	}

	// do not combine layers which are of quite different hardness 020917;Fz
	if ( fabs(Pdata[e1].hard - Pdata[e2].hard) > 1.0 ) {
		aggreg = OFF;
	}

	if ( (Pdata[e1].dendricity == 0) && (Pdata[e2].dendricity == 0) ) {
		if ( fabs(Pdata[e1].sphericity - Pdata[e2].sphericity) > DIFF_SP) {
			aggreg = OFF;
		}
		if ( fabs(Pdata[e1].grain_dia - Pdata[e2].grain_dia) >
		MAX(DIFF_DG, DIFF_DG_REL*Pdata[e1].grain_dia) ) {
			aggreg = OFF;
		}
	} else {
		if ( fabs(Pdata[e1].sphericity - Pdata[e2].sphericity) > DIFF_SP) {
			aggreg = OFF;
		}
		if ( fabs(Pdata[e1].dendricity - Pdata[e2].dendricity) > DIFF_DD) {
			aggreg = OFF;
		}
	}

	return(aggreg);
} // End of ml_ag_Join

/**
 * @brief Determine the grain class
 * Revisited by Fierz and Bellaire fall 2006
 * @param dendricity
 * @param sphericity
 * @param grain_dia
 * @param marker
 * @param theta_w
 * @param theta_i
 */
int ml_ag_Classify(const double dendricity, const double sphericity, const double grain_dia, const int marker, const double theta_w, const double theta_i)
{
	int a=-1,b=-1,c=0;
	int sw2;
	double res_wat_cont;

	res_wat_cont = lw_SnowResidualWaterContent(theta_i);

	// Dry snow
	if( dendricity > 0. ) {
		// Dry dendritic (new) snow: dendricity and sphericity determine the class
		sw2 = (int)(sphericity*10.);
		if( dendricity > 0.80 ) {
			// ori 0.90, 27 Nov 2007 sb
			a = 1; b = 1; c = 0;
		} else if( dendricity > 0.70 ) {
			// ori 0.85, 27 Nov 2007 sb
			a = 1; b = 2; c = 1;
		} else if( dendricity > 0.65 ) {
			// ori 0.75, 27 Nov 2007 sb
			a = 2; b = 1; c = 0;
		} else if( dendricity > 0.60 ) {
			// ori 0.70, 27 Nov 2007 sb
			a = 2; b = 1; c = 1;
		} else if( dendricity > 0.30 ) {
			a = 2; b = 2; c = 0;
		} else if( dendricity > 0.05 ) {
			a = 2;
			switch(sw2) {
				case 0: case 1: case 2:
				b = 4; c = 0; break;
				case 3: case 4:
				b = 4; c = 1; break;
				case 5: case 6:
				b = 3; c = 1; break;
				default:
				b = 3; c = 0;
			}
		} else {
			switch(sw2) {
				case 0: case 1:
				a = 4; b = 4; c = 0; break;
				case 2: case 3: case 4:
				a = 4; b = 2; c = 1; break;
				case 5: case 6: case 7:
				a = 3;  b = 2; c = 1; break;
				default:
				a = 3; b = 3; c = 0;
			}
		}
	} else if( marker <= 2) {
		/*
		 * Dry non-dendritic snow
		 * Sphericity is most important for "a", while the marker is most important for "b","c"
		 */
		if( grain_dia < 0.7 ) {
			sw2 = (int)(sphericity*10.);
			switch(sw2) {
				case 0: case 1:
				a = 4; b = 4; c = 0; break;
				case 2: case 3:
				a = 4; b = 3; c = 1; break;
				case 4: case 5:
				a = 4;  b = 3; c = 0; break;
				case 6: case 7:
				a = 3;  b = 4; c = 1; break;
				default:
				a = 3; b = 3; c = 0;
			}
		} else if( grain_dia < 1.1 ) {
			if( sphericity < 0.2 ) {
				a = 4; b = 4; c = 0;
			} else if( sphericity < 0.4 ) {
				a = 4; b = 9; c = 0;
			} else {
				// sphericity limited to sp_max=0.5 in Metamorphism.c
				a = 9; b = 9 ; c = 0;
			}
		} else if( grain_dia < 1.5 ) {
			if( sphericity < 0.2 ) {
				a = 4; b = 5; c = 0;
			} else if( sphericity < 0.4 ) {
				a = 4; b = 9; c = 1;
			} else {
				// sphericity limited to sp_max=0.5 in Metamorphism.c
				a = 9; b = 9 ; c = 0;
			}
		} else {
			if( sphericity < 0.2 ) {
				a = 5; b = 5; c = 0;
			} else if( sphericity < 0.4 ) {
				a = 5; b = 9; c = 1;
			} else {
				// sphericity limited to sp_max=0.5 in Metamorphism.c
				a = 9; b = 5 ; c = 1;
			}
		}
	} // end dry snow

	// Snow getting wet
	if( marker >= 10 ) {
		if( dendricity > 0.0 ) {
			// Wet dendritic snow
			if( sphericity > 0.7 ) {
				b = a; a = 7; c = 0;
			} else {
				b = 7 ; c = 1;
			}
		} else {
			// Wet non-dendritic snow
			b = 7; c = 0;
			if( sphericity > 0.75) {
				a = 7;
			} else if( sphericity > 0.4 ) {
				if( grain_dia <= 0.7 ) {
					b = a; a = 7;
				} else if( marker != 13 ) {
					if( grain_dia <= 1.5 ) {
						a = 7; b = 9; c = 1;
					} else {
						a = 7; b = 5; c = 1;
					}
				} else {
					a = 7; b = 6;  c = 1;
				}
			} else {
				if( grain_dia <= 1.5 ) {
					a = 4;
				} else {
					a = 5;
				}
				if( sphericity <= 0.2 ) {
					c = 1;
				}
			}
		}
	}

	// Now treat a couple of exceptions - note that the order is important
	if( b < 0 ) {
		b = a;
	}
	// Melt-Freeze
	if( (marker >= 20) && (theta_w < 0.1*res_wat_cont) ) {
		c = 2;
	}
	// Surface Hoar
	if( marker == 3 ) {
		 a = 6;
		 b = 6;
		 c = 0;
	}
	// Graupel
	if( marker == 4 ) {
		 a = 0;
		 b = 0;
		 c = 0;
	}
	// Ice Layer
	if( marker % 10 == 8 ) {
		 a = 8;
		 b = 8;
		 c = 0;
	}

	return (a*100 + b*10 + c);

} // End of ml_ag_Classify


/**
 * @brief Aggregate the Layers and calculate the grain class
 * @param nE (int)
 * @param *Pdata (Q_PROFILE_DAT)
 */
int ml_ag_Aggregate(int nE, Q_PROFILE_DAT *Pdata)
{
	int e, nL, l, flag=OFF;
	double l1, l2;

	// Initialize number of elements  and aggregate only if more than 5 layers
	nL = nE;
	if (nL > 5) {
		// First Run -  leave top element alone

		// keep track of the coordinates and length of elements
		l2 = (Pdata[nE-2].height -  Pdata[nE-3].height);
		for(e=nE-2; e>0; e--) {
			l = e-1;
			l1 = l2;
			if (l>0) {
				l2 = (Pdata[l].height -  Pdata[l-1].height);
			} else {
				l2 = Pdata[l].height;
			}
			// if two layers are similar combine them
			if ( (ml_ag_Join(e, Pdata) == ON) && (Pdata[l].marker != 3) && (Pdata[e].marker != 3) ) {
				nL--;
				ml_ag_Average(e, l, l1, l2, Pdata);
				l2 += l1;
			}

		}  // for all elements

		ml_ag_Shift(nE, Pdata);
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
				if ( (Pdata[l].marker != 3) && (Pdata[e].marker != 3) ) {
					// trick to try to join with upper or lower level -> use flag to mark thin layer
					if ( (flag == ON) || (l2 < (sqrt(Pdata[nE-1].height-Pdata[l].height)/4.)) || (l2 < MIN_L_ELEMENT) ) {
						// if two layers are similar or one layer is very very small combine them
						if ( (ml_ag_Join2(e, Pdata) == ON) || (l2 < MIN_L_ELEMENT) ||
							(l1 < MIN_L_ELEMENT)  ) {
							nL--;
							ml_ag_Average(e, l, l1, l2, Pdata);
							l2 += l1;
							flag = OFF;
						} else {
							flag = ON;
						}
					} else {
						flag = OFF;
					}
				} // if not surface hoar layer
				else {
					flag = OFF;
				}

			}  // for all elements

			ml_ag_Shift(nE, Pdata);
		} // if nE > 2

	} // if more than 5 layers

	// Now Calculate the grain class. The grain class is coded according to the Matt&Sommer profile visualization
	for(e=0; e<nL; e++) {
		Pdata[e].grain_class = ml_ag_Classify(Pdata[e].dendricity, Pdata[e].sphericity,
						Pdata[e].grain_dia, Pdata[e].marker,
						Pdata[e].theta_w/100., Pdata[e].theta_i/100.);
	}

	return (nL);

} // End of ml_ag_Aggregate

/*
 * End of Aggregate.c
 */
