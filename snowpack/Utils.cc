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
 * @file Utils.cc
 * @version 10.02
 * @brief This module contains all-purpose functions
 */

#include <snowpack/Utils.h>

using namespace std;
using namespace mio;

/**
 * @brief Print a message to screen (see p.351-352[,469] in C: The Complete Reference) \n
 * If JulianDate < 0., no running date will be written \n
 * If JulianDate >= 0., compute current date and write it (\<t>) \n
 * The output format depends on message type:
 * - "err"  : [E] [\<t>] [\<file>:\<line>] \<on JulianDate> \<msg>
 * - "wrn"  : [W] [\<t>] [\<file>:\<line>] \<on JulianDate> \<msg>
 * - "msg+" : [I] [\<t>] [\<file>:\<line>] \<on JulianDate> \<msg>
 * - "msg"  : [i] [\<t>] ---> \<msg>
 * - "msg-" : [i] []      \<msg>
 * @author Charles Fierz \n Mathias Bavay
 * @version 9.mm
 * @param *theFile
 * @param theLine
 * @param *msg_type See above
 * @param JulianDate Set to -1. if JulianDate is not available.
 * @param *format Format for message
 * @param ... Variable number of parameters to format
 */
void prn_msg(const char *theFile, int theLine, const char *msg_type, double JulianDate, const char *format, ...)
{
	va_list argptr; // get an arg ptr

	int msg_ok = 0;

	// Initialize argptr to point to the first argument after the format string
	va_start(argptr, format);

	//calculate time stamp
	string currentdate = Date(time(NULL), 1.).toString(Date::ISO); //default HACK: we should read TZ from io.ini
	if ( JulianDate > 0. ) currentdate = Date(JulianDate).toString(Date::ISO);
   	if ( JulianDate < 0. ) currentdate = "";


	//print message
	//printf("¬"); //if we need multiline output, use a special char as bloc delimiter
	if ( strcmp(msg_type, "err") == 0 ) {
		fprintf(stdout, "[E] [%s] [%s:%d] ",currentdate.c_str(), theFile, theLine);
		msg_ok=1;
	}
	if ( strcmp(msg_type, "wrn") == 0 ) {
		fprintf(stdout, "[W] [%s] [%s:%d] ",currentdate.c_str(), theFile, theLine);
		msg_ok=1;
	}
	if ( strcmp(msg_type, "msg+") == 0 ) {
		fprintf(stdout, "[I] [%s] [%s:%d] ",currentdate.c_str(), theFile, theLine);
		msg_ok=1;
	}
	if ( strcmp(msg_type, "msg-") == 0 ) {
		fprintf(stdout, "[i] []      ");
		msg_ok=1;
	}
	if ( strcmp(msg_type, "msg") == 0 ) {
		fprintf(stdout, "[i] [%s] ---> ",currentdate.c_str());
		msg_ok=1;
	}

	if ( msg_ok ) {
		vfprintf(stdout, format, argptr);
	} else {
		fprintf(stdout, "[W] [%s] [%s:%d] Message type '%s' unknown!", currentdate.c_str(), theFile, theLine, msg_type);
	}

	fprintf(stdout, "\n");

	// Clear ptr
	va_end(argptr);
}

/**
 * @brief Determines whether a certain time has been reached (e.g. to dump output)
 * The function returns TRUE whenever JulianDate >= start && JulianDate == start + n*days_between, n=0,1,2,...
 * - Start must be a Julian Date with origin 1900-01-01T00:00
 * - In case regular dumps are requested throughout the day, it is best to set start to 0.0.
 *   This is the preferred setting in operational mode with *_START = 0.0
 * - This version is the result of various efforts (Lehning, Kowalski, Fierz, Loewe)
 * @author Michael Lehning \n Julia Kowalski \n Charles Fierz \n Mathias Bavay
 * @version 9.mm
 * @param JulianDate Julian Date
 * @param days_between (double) number of days between two outputs
 * @param start (const double) start date as Julian Date
 * @return int
 */
int qr_BooleanTime(const double& JulianDate, double days_between, 
			    const double& start, const double& calculation_step_length)
{
	int ret;
	double jul_frc;
	const double step = M_TO_D(calculation_step_length);	//step length in days (converted from minutes)

	if ( JulianDate < (start - 0.5*step) ) {
		return 0;
	}

	// NOTE that days_between must be known to a high degree of accuracy to make the test below work
	//computing the output interval in units of "step" (ie: number of "step")
	//days_between = round(days_between/step) * step; //only in C99
	days_between = floor(days_between / step + 0.5) * step;//how to implement a replacement to round() using only floor()!
	if (days_between == 0.) {
		prn_msg(__FILE__, __LINE__, "err", -1., "Days_between is zero. Please consider changing data output intervals!");
		return 0;
	}
	jul_frc = (JulianDate - start) / days_between - floor((JulianDate - start) / days_between);
	ret = ( (jul_frc > (days_between - 0.5 * step) / days_between) || (jul_frc < 0.5 * step / days_between) );
	return ret;
}

/**
 * @brief Delete outdir.ext_n[ext] from outdir
 * @author Michael Lehning
 * @version 9Y.mm
 * @param outdir Output dir
 */
void deleteOldOutputFiles(const std::string& outdir, const std::string& experiment, 
					 const std::string& station, const int& number_expo)
{
	vector<string> vecExtension;
	vecExtension.push_back("sno"); //Snow-cover profile file (I/O)
	vecExtension.push_back("met"); //Meteo data input
	vecExtension.push_back("pro"); //Time series of modeled profile-type data
	vecExtension.push_back("ini"); //Record of run configuration

	char fp[MAX_STRING_LENGTH]="\000", exp[MAX_STRING_LENGTH]="\000";
	int j, n_files;

	if (experiment != "none"){
		snprintf(exp, MAX_STRING_LENGTH-2, "%s_%s", station.c_str(), experiment.c_str());
	}
	prn_msg(__FILE__, __LINE__, "msg+", -1., "Erasing old result file(s) %s*%s*", outdir.c_str(), exp);
	for (unsigned int ii=0; ii<vecExtension.size(); ii++){
		const string& ext = vecExtension[ii];

		cout << "Deleting " << ext << " files" << endl;

		n_files = 0;

		if ((ext == "sno") || (ext == "met") || (ext == "pro")){
			for (j = 0; j < number_expo; j++) {
				if ( j ) {
					snprintf(fp, MAX_STRING_LENGTH-1, "%s%d%s.%s", outdir.c_str(), j, exp, ext.c_str());
				} else {
					snprintf(fp, MAX_STRING_LENGTH-1, "%s%s.%s", outdir.c_str(), exp, ext.c_str());
				}
				if ( remove(fp) == 0 ) {
					n_files++;
				}
			}
			if ( n_files > 0 ) {
				prn_msg(__FILE__, __LINE__, "msg-", -1., "Erased %d *.%s file(s)", n_files, ext.c_str());
			} else {
				prn_msg(__FILE__, __LINE__, "msg-", -1., "No *.%s file(s) to erase", ext.c_str());
			}
		} else if (ext == "ini"){
			if (station != "IMIS"){
				snprintf(fp, MAX_STRING_LENGTH-2, "%s%s.%s", outdir.c_str(), exp, ext.c_str());
				if ( number_expo > 1 ) {
					snprintf(fp, MAX_STRING_LENGTH-3, "%s%s-%d.%s", outdir.c_str(), exp, number_expo-1, ext.c_str());
				}
				if ( remove(fp) == 0 ) {
					prn_msg(__FILE__, __LINE__, "msg-", -1., "Erased %s", fp);
				} else {
					prn_msg(__FILE__, __LINE__, "msg-", -1., "No file %s to erase", fp);
				}
			}
		}
	}
}

/**
 * @brief Returns modelled internal snow or/and soil temperature (instantaneous value; degC),
 *        at a given position z perpendicular to slope (m) \n
 *        z must be less than calculated height (Xdata->cH), otherwise modeled temperature is set to NODATA
 * @author Charles Fierz
 * @version 10.02
 * @param z Sensor position perpendicular to slope (m)
 * @param *Xdata
 */
double getModelledTemperature(const double& z, const SN_STATION_DATA& Xdata)
{
	int n_up;           // Upper node number
	double z_up, z_low; // Upper and lower nodes around position z of sensor

	const vector<SN_NODE_DATA>& NDS = Xdata.Ndata;
	if ( (z == NODATA) || !((Xdata.getNumberOfNodes() > 1) && (z < Xdata.cH)) ) {
		return NODATA;
	} else {
		n_up = findUpperNode(z, NDS, Xdata.getNumberOfNodes());
		z_low = (NDS[n_up-1].z + NDS[n_up-1].u);
		z_up = (NDS[n_up].z + NDS[n_up].u);
		return (K_TO_C(NDS[n_up-1].T + (z - z_low)*(NDS[n_up].T-NDS[n_up-1].T)/(z_up-z_low)));
	}
}


/**
 * @brief Returns number of lowest node above a given position z perpendicular to slope (m) \n
 * @author Charles Fierz
 * @version 10.02
 * @param z Position perpendicular to slope (m)
 * @param *Ndata
 * @param nN Number of nodes
 * @return Upper node number
 */
int findUpperNode(const double& z, const vector<SN_NODE_DATA>& Ndata, const int& nN)
{
	int n_up = nN-2;
	double z_low = Ndata[n_up].z + Ndata[n_up].u;
	while ( z < z_low &&  n_up > 0 ) {
		n_up--;
		z_low = Ndata[n_up].z + Ndata[n_up].u;
	}
	return ++n_up;
}

/**
 * @brief Fill the snowpack version number, date of computation, user, ...
 * @author Mathias Bavay
 * @version 8.mm
 * @param *version
 * @param *computation_date iso-format
 * @param *jul_computation_date Julian Date
 * @param *user
 */
void qr_VersionUserRuntime(char *version, char *computation_date, double *jul_computation_date, 
					  char *user, mio::Date& date)
{
	char *logname;

	time_t rawtime; // time in seconds since EPOC
	struct tm * timeinfo; // local time structure

	time (&rawtime);
	timeinfo = localtime (&rawtime);
	time_t ltime = mktime(timeinfo);
	Date localdate(ltime, 1.); //HACK: use TZ from io.ini
	
	date = localdate;

	// version and computation time
	snprintf(version, MAX_STRING_LENGTH-1, "%s", SN_VERSION);
	snprintf(computation_date, MAX_STRING_LENGTH-1, "%s", localdate.toString(Date::ISO).c_str());

	*jul_computation_date = localdate.getJulianDate();
	//logname=getlogin(); //other options possible, see man
	logname = getenv("LOGNAME");
	if ( logname == NULL ) {
		snprintf(user, MAX_STRING_LENGTH, "N/A");
	} else {
		snprintf(user, MAX_STRING_LENGTH, "%s", logname);
	}
}

/**
 * @brief Averages energy fluxes
 * @author Charles Fierz
 * @version 9Y.mm
 * @param Sdata
 * @param Xdata
 * @param n_steps Number of calculation time steps since last output
*/
void qr_AverageFluxTimeSeries(const int& n_steps, const bool& useCanopyModel, 
						SN_SURFACE_DATA& Sdata, SN_STATION_DATA& Xdata)
{
	// Mean energy input (J m-2)
	Sdata.dIntEnergy  /= n_steps;
	// Mean energy fluxes (W m-2), including albedo
	Sdata.lw_in   /= n_steps;
	Sdata.lw_out  /= n_steps;
	Sdata.lw_net  /= n_steps;
	Sdata.qs      /= n_steps;
	Sdata.ql      /= n_steps;
	Sdata.qr      /= n_steps;
	Sdata.qg      /= n_steps;
	Sdata.qg0     /= n_steps;
	Sdata.sw_hor  /= n_steps;
	Sdata.sw_in   /= n_steps;
	Sdata.sw_out  /= n_steps;
	Sdata.qw      /= n_steps;
	Sdata.sw_dir  /= n_steps;
	Sdata.sw_diff /= n_steps;
	Sdata.cA      /= n_steps;
	Sdata.mA      /= n_steps;

	if (useCanopyModel) {
		// *radiation
		Xdata.Cdata.rswrac /= n_steps;
		Xdata.Cdata.iswrac /= n_steps;
		Xdata.Cdata.rswrbc /= n_steps;
		Xdata.Cdata.iswrbc /= n_steps;
		Xdata.Cdata.ilwrac /= n_steps;
		Xdata.Cdata.rlwrac /= n_steps;
		Xdata.Cdata.ilwrbc /= n_steps;
		Xdata.Cdata.rlwrbc /= n_steps;
		Xdata.Cdata.rsnet /= n_steps;
		Xdata.Cdata.rlnet /= n_steps;
		// turbulent heat fluxes
		Xdata.Cdata.sensible /= n_steps;
		Xdata.Cdata.latent /= n_steps;
		Xdata.Cdata.latentcorr /= n_steps;
		// auxiliaries
		Xdata.Cdata.canopyalb /= n_steps;
		Xdata.Cdata.totalalb /= n_steps;
		Xdata.Cdata.intcapacity /= n_steps;
	}
}

/**
 * @brief Decompose type in its constituents
 * At present decomposition into Swiss numerical code \n
 * TODO Work with new international code
 * @author Charles Fierz
 * @version 9.12
 * @param F1 Majority grain shape
 * @param F2 Minority grain shape
 * @param F3 2 indicates a melt-freeze crust
 * @param type aggregated shape information
 */
void qr_TypeToCode(int *F1, int *F2, int *F3, int type)
{
	*F1   = (int)(floor(type/100.));
	type -= (int)((*F1)*100);
	*F2   = (int)(floor(type/10.));
	*F3   = (int)(type - (*F2)*10);
}

/**
 * @brief Performs mass balance check either before or after a calculation time step
 * - NOTE: AVGSUM_TIME_SERIES should not be set
 * @author Charles Fierz
 * @version 10.03
 * @param *Xdata
 * @param *Sdata
 * @param *tot_mass_in Total mass after last time step (kg m-2)
 * @return bool if mass error occured, putting balance terms on screen
 */
bool qr_MassBalanceCheck(const SN_STATION_DATA& Xdata, const SN_SURFACE_DATA& Sdata, double& tot_mass_in)
{
	bool mass_error = true;
	int e;
	double tot_mass=0., tot_swe=0., dmassE=0.;
	double hnw = Xdata.hn_slope*Xdata.rho_slope;
	double mass_change = hnw - Sdata.mass[SN_SURFACE_DATA::MS_RUNOFF] + Sdata.mass[SN_SURFACE_DATA::MS_RAIN] + Sdata.mass[SN_SURFACE_DATA::MS_SUBLIMATION] + Sdata.mass[SN_SURFACE_DATA::MS_EVAPORATION] - MAX(0., Xdata.ErosionMass);

	// Actual mass of snowpack
	for (e=Xdata.SoilNode; e<Xdata.getNumberOfElements(); e++) {
		tot_mass += Xdata.Edata[e].M;
		tot_swe  += Xdata.Edata[e].L * Xdata.Edata[e].Rho;
		dmassE = Xdata.Edata[e].M - (Xdata.Edata[e].L * Xdata.Edata[e].Rho);
		if ( fabs(dmassE) > Constants::eps ) {
			prn_msg(__FILE__, __LINE__, "msg", -1., "Mass error at element e=%d (nE=%d): mass(now)=%lf swe(now)=%lf dmassE=%lf", e, Xdata.getNumberOfElements(), tot_mass, tot_swe, dmassE);
			mass_error = false;
		}
	}
	// Mass balance check
	if ( tot_mass > Constants::eps ) {
		if ( (fabs(tot_swe/tot_mass) - 1.) > 0.5e-2 ) {
			prn_msg(__FILE__, __LINE__, "msg", -1., "Mass balance (theta): mass(now)=%lf swe(now)=%lf swe/mass=%lf mass-swe=%lf", tot_mass, tot_swe, tot_swe/tot_mass, tot_mass - tot_swe);
			mass_error = false;
		}
		if ( tot_mass_in > Constants::eps ) {
			if ( (fabs(tot_mass - (tot_mass_in + mass_change))/tot_mass) > 0.5e-4 ) {
				prn_msg(__FILE__, __LINE__, "msg", -1., "Mass balance: mass_err(now)=%lf", tot_mass - (tot_mass_in + mass_change));
				prn_msg(__FILE__, __LINE__, "msg", -1., "tot_mass_in=%lf tot_mass_now=%lf mass_change=%lf", tot_mass_in, tot_mass, mass_change);
				mass_error = false;
			}
		} else {
			tot_mass_in = tot_mass;
		}
	} else if ( tot_mass_in > Constants::eps ) {
		if ( fabs(tot_mass_in + mass_change) > 1.0e-3 ) {
			prn_msg(__FILE__, __LINE__, "msg", -1., "Mass balance error: tot_mass_in=%lf tot_mass_now=%lf mass_change=%lf", tot_mass_in, tot_mass, mass_change);
			mass_error = false;
		}
	}

	return mass_error;
}

/**
 * @brief Forced erosion of a missed event
 * @author Michael Lehning \n Mathis Bavay
 * @date 2008-03-07
 * @param hs1 Snow depth to correct
 * @param *Xdata
 * @return Eroded mass (kg m-2)
 */
double qro_ForcedErosion(const double hs1, SN_STATION_DATA *Xdata)
{
	int    nErode=0;        // Counters
	double massErode=0.;    // Eroded mass (kg m-2)

	massErode=0.;
	while ( (Xdata->getNumberOfElements() > Xdata->SoilNode) && (hs1 + 0.01) < (Xdata->cH - Xdata->Ground) ) {
		massErode += Xdata->Edata[Xdata->getNumberOfElements()-1].M; 
		Xdata->cH -= Xdata->Edata[Xdata->getNumberOfElements()-1].L; 
		Xdata->resize(Xdata->getNumberOfElements() - 1);
		nErode++;
	}
	Xdata->ErosionLevel = MIN(Xdata->getNumberOfElements()-1, Xdata->ErosionLevel);

	return(massErode);
}


/**
 * @brief Michi is very unhappy that the warning service wants to inflate or deflate the
 * snow if there is a consistent under- oder overestimation of settling respectively. \n
 * But since the will of the warning service is law at the SLF, we have no choice
 * and need to implement this additional terrible non mass-conserving and cheating
 * feature. But it will make the operational users happy, I hope. \n
 * Implemented on 2 Feb 2008 (and 7 Mar 2008: back to erosion) by Mathias Bavay
 * and Michi, who should be home with Leo being sick ...
 * @param *Mdata
 * @param *Xdata
 * @param *time_counter
 * @param *dhs_corr Correction on snow depth (m)
 * @param *mass_corr Mass correction (kg m-2)
 */
void qro_DeflateInflate(const SN_MET_DATA *Mdata, SN_STATION_DATA *Xdata, double *dhs_corr, double *mass_corr)
{
	int    e, nE, nSoil;                         // Element counter
	double factor_corr, sum_total_correction=0.; // Correction factor
	double ddL, dL=0.;                           // Length changes
	double cH, cH_old;                           // Snow depth

	// Dereference a few values
	vector<SN_NODE_DATA>& NDS = Xdata->Ndata;
	vector<SN_ELEM_DATA>& EMS = Xdata->Edata;
	nE = Xdata->getNumberOfElements(); nSoil = Xdata->SoilNode;
	cH = Xdata->cH - Xdata->Ground;
	/*
	 * First try to find erosion events, which have not been captured by the drift module
	 * (Maybe the wind sensor did not measure correctly due to riming, or s.th. else went
	 * wrong with the model. For now assume erosion if more than 3 cm are missing
	 */
	if ( (Mdata->hs1 + 0.03) < cH ) {
		*dhs_corr = Mdata->hs1 - cH;
		*mass_corr = qro_ForcedErosion(Mdata->hs1, Xdata);
		if ( 0 ) {
			prn_msg(__FILE__, __LINE__, "msg+", Mdata->date.getJulianDate(), "Missed erosion event detected");
			prn_msg(__FILE__, __LINE__, "msg-", -1., "Measured Snow Depth:%lf   Calculated Snow Depth:%lf", Mdata->hs1, cH);
		}
	} else {
		// assume settling error
		*dhs_corr = Mdata->hs1 - cH;
		*mass_corr = 0.;

		//Test whether normalization quantity does not lead to an arithmetic exception
		//This is a work around for weird cases in which the whole snowpack appears at once
		if (EMS[nE-1].date.getJulianDate() <= EMS[nSoil].date.getJulianDate())
			return;

		if ( 0 ) {
			prn_msg(__FILE__, __LINE__, "msg+", Mdata->date.getJulianDate(), 
				   "Small correction due to assumed settling error\n");
			prn_msg(__FILE__, __LINE__, "msg-", -1., 
				   "Measured Snow Depth:%lf   Calculated Snow Depth:%lf", Mdata->hs1, cH);
		}
		// Second find the normalization quantity, which we choose to be the age of the layer.
		for (e = nSoil; e < nE; e++) {
			if ( (!(EMS[e].mk > 20 || EMS[e].mk == 3))&&(Mdata->date.getJulianDate() > EMS[e].date.getJulianDate())){
				sum_total_correction += EMS[e].L
					* (1. - sqrt((EMS[nE-1].date.getJulianDate() - EMS[e].date.getJulianDate()) 
							   / (EMS[nE-1].date.getJulianDate() - EMS[nSoil].date.getJulianDate())));
			}
		}
		if ( sum_total_correction > 0. ) {
			factor_corr = (Mdata->hs1 - cH) / sum_total_correction;
		} else {
			*dhs_corr = 0.;
			return;
		}
		// ... above marked element (translation only) ...
		// Squeeze or blow-up
		for (e = nSoil; e < nE; e++) {
			if ( (!(EMS[e].mk > 20 || EMS[e].mk == 3)) && (Mdata->date.getJulianDate() > EMS[e].date.getJulianDate())){
				ddL = EMS[e].L 
					* MAX(-0.9, MIN(0.9, factor_corr 
								 * (1. - sqrt((EMS[nE-1].date.getJulianDate() - EMS[e].date.getJulianDate()) 
										    / (EMS[nE-1].date.getJulianDate() - EMS[nSoil].date.getJulianDate())))));
			} else {
				ddL = 0.;
			}
			dL += ddL;
			*mass_corr += ddL*EMS[e].Rho;
			NDS[e+1].z += dL + NDS[e+1].u;
			NDS[e+1].u  = 0.0;
			EMS[e].M += ddL*EMS[e].Rho;
			EMS[e].L0 = EMS[e].L += ddL;
			EMS[e].E  = EMS[e].dE = EMS[e].Ee = EMS[e].Ev = EMS[e].S = 0.0;
		}
		// Update the overall height
		cH_old     = Xdata->cH;
		Xdata->cH  = NDS[nE].z + NDS[nE].u;
		Xdata->mH -= (cH_old - Xdata->cH);
	}
} // End qro_ErodeDeflateInflate

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
int ml_ag_Classify(const double& dendricity, const double& sphericity, const double& grain_dia, 
			    const int& marker, const double& theta_w, const double& theta_i)
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
