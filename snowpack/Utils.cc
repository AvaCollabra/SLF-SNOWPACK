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
 * @file Utils.cc
 * @version 10.02
 * @brief This module contains all-purpose functions
 */

#include <snowpack/Utils.h>

using namespace std;
using namespace mio;

namespace snowpack {
std::string getLibVersion() {
	std::stringstream ss;
	ss << _VERSION << " compiled on " << __DATE__ << " " << __TIME__;
	return ss.str();
}
}

/**
 * @brief Print a message to screen (see p.351-352[,469] in C: The Complete Reference) \n
 * If date_in is Date(), no running date will be written, \n
 * otherwise write date_in to (\<t>) \n
 * The output format depends on message type:
 * - "err"  : [E] [\<t>] [\<file>:\<line>] \<on date_in> \<msg> \\n
 * - "wrn"  : [W] [\<t>] [\<file>:\<line>] \<on date_in> \<msg> \\n
 * - "msg+" : [I] [\<t>] [\<file>:\<line>] \<on date_in> \<msg> \\n
 * - "msg"  : [i] [\<t>] ---> \<msg> \<n>
 * - "msg-" : [i] []      \<msg> \<n>
 * @author Charles Fierz \n Mathias Bavay
 * @version 11.02
 * @param *theFile
 * @param theLine
 * @param *msg_type See above
 * @param date_in Use Date() if date_in is not available.
 * @param *format Format for message
 * @param ... Variable number of parameters to format
 */
void prn_msg(const char *theFile, const int theLine, const char *msg_type, const mio::Date& date_in, const char *format, ...)
{
	va_list argptr; // get an arg ptr

	int msg_ok = 0;

	// Initialize argptr to point to the first argument after the format string
	va_start(argptr, format);

	//compute time stamp
	string currentdate;
	if (date_in.isUndef()) {
		Date msg_date;
		msg_date.setFromSys(); //HACK: we don't put the timezone...
		currentdate = msg_date.toString(Date::ISO);
	} else {
		currentdate = date_in.toString(Date::ISO);
	}

	//print message
	//printf("¬"); //if we need multiline output, use a special char as bloc delimiter
	if (strcmp(msg_type, "err") == 0) {
		fprintf(stdout, "[E] [%s] [%s:%d] ", currentdate.c_str(), theFile, theLine);
		msg_ok=1;
	}
	if (strcmp(msg_type, "wrn") == 0) {
		fprintf(stdout, "[W] [%s] [%s:%d] ", currentdate.c_str(), theFile, theLine);
		msg_ok=1;
	}
	if (strcmp(msg_type, "msg+") == 0) {
		fprintf(stdout, "[I] [%s] [%s:%d] ", currentdate.c_str(), theFile, theLine);
		msg_ok=1;
	}
	if (strcmp(msg_type, "msg-") == 0) {
		fprintf(stdout, "[i] []                 ");
		msg_ok=1;
	}
	if (strcmp(msg_type, "msg") == 0) {
		fprintf(stdout, "[i] [%s] ---> ", currentdate.c_str());
		msg_ok=1;
	}

	if (msg_ok) {
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
 * @param days_between number of days between two outputs
 * @param start start date as Julian Date
 * @param calculation_step_length (min) //HACK: it should be in SI!!
 * @return int
 */
bool booleanTime(const double& JulianDate, double days_between,
                 const double& start, const double& calculation_step_length)
{
	int ret;
	double jul_frc;
	const double step = M_TO_D(calculation_step_length);	//step length in days (converted from minutes)

	if ( JulianDate < (start - 0.5*step) ) {
		return false;
	}

	// NOTE that days_between must be known to a high degree of accuracy to make the test below work
	//computing the output interval in units of "step" (ie: number of "step")
	//days_between = round(days_between/step) * step; //only in C99
	days_between = floor(days_between / step + 0.5) * step;//how to implement a replacement to round() using only floor()!
	if (days_between == 0.) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Days_between is zero. Please consider changing data output intervals!");
		return 0;
	}
	jul_frc = (JulianDate - start) / days_between - floor((JulianDate - start) / days_between);
	ret = ( (jul_frc > (days_between - 0.5 * step) / days_between) || (jul_frc < 0.5 * step / days_between) );
	return (ret != 0);
}

/**
 * @brief Delete old output files (*.sno, *.ini) from outdir
 * @version 11.03
 * @param outdir Output dir
 * @param experiment Name of ongoin experiment
 * @param stationID
 * @param nSlopes Number of slopes treated
 */
void deleteOldOutputFiles(const std::string& outdir, const std::string& experiment,
                          const std::string& stationID, const unsigned int& nSlopes)
{
	vector<string> vecExtension;
	vecExtension.push_back("met"); //Meteo data input
	vecExtension.push_back("pro"); //Time series of modeled profile-type data
	vecExtension.push_back("ini"); //Record of run configuration
	vecExtension.push_back("sno"); //Snow-cover profile file (I/O)

	char ftrunc[MAX_STRING_LENGTH]="\000", fname[MAX_STRING_LENGTH]="\000", exp[MAX_STRING_LENGTH]="\000";
	unsigned int n_files;

	if (experiment != "NO_EXP") {
		snprintf(exp, MAX_STRING_LENGTH-2, "%s_%s", stationID.c_str(), experiment.c_str());
	}
	for (unsigned int ii=0; ii<vecExtension.size(); ii++){
		const string& ext = vecExtension[ii];
		n_files = 0;

		snprintf(ftrunc, MAX_STRING_LENGTH-1, "%s/%s", outdir.c_str(), exp);
		if (ext == "ini") {
			if (stationID != "IMIS") {
				snprintf(fname, MAX_STRING_LENGTH-2, "%s.%s", ftrunc, ext.c_str());
				if (nSlopes > 1) {
					snprintf(fname, MAX_STRING_LENGTH-3, "%s-%d.%s", ftrunc, nSlopes-1, ext.c_str());
				}
				if (remove(fname) == 0) {
					prn_msg(__FILE__, __LINE__, "msg-", Date(), "Erased %s", fname);
				} else {
					prn_msg(__FILE__, __LINE__, "msg-", Date(), "No file %s to erase", fname);
				}
			}
		} else {
			for (unsigned int j = 0; j < nSlopes; j++) {
				if ( j ) {
					snprintf(fname, MAX_STRING_LENGTH-1, "%s%d.%s", ftrunc, j, ext.c_str());
				} else {
					snprintf(fname, MAX_STRING_LENGTH-1, "%s.%s", ftrunc, ext.c_str());
				}
				if (ext == "sno") {
					if (remove(fname) == 0) {
						n_files++;
					}
					if (j == nSlopes-1) {
						if (n_files > 0) {
							prn_msg(__FILE__, __LINE__, "msg-", Date(), "Erased %d *.%s file(s)", n_files, ext.c_str());
						} else {
							prn_msg(__FILE__, __LINE__, "msg-", Date(), "No *.%s file(s) to erase", ext.c_str());
						}
					}
				} else if ((j == 0) && IOUtils::fileExists(fname)) {
					prn_msg(__FILE__, __LINE__, "msg-", Date(), "Data in *.%s file(s) may be overwritten", ext.c_str());
				}
			}
		}
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
int findUpperNode(const double& z, const vector<NodeData>& Ndata, const int& nN)
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
 * @param time_zone
 * @param *version
 * @param *computation_date iso-format
 * @param *jul_computation_date Julian Date
 * @param *compilation_date Julian Date
 * @param *user
 */
void versionUserRuntime(const double& time_zone, char *version, char *computation_date, double *jul_computation_date,
                        char *compilation_date, char *user)
{
	Date localdate;
	localdate.setFromSys();
	localdate.setTimeZone(time_zone);

	// version as well as computation and compilation time
	strncpy(version, SN_VERSION, MAX_STRING_LENGTH-1);
	snprintf(computation_date, MAX_STRING_LENGTH-1, "%s", localdate.toString(Date::ISO).c_str());
	snprintf(compilation_date, MAX_STRING_LENGTH-1, "%s, %s", __DATE__, __TIME__);

	*jul_computation_date = localdate.getJulianDate();
	//logname=getlogin(); //other options possible, see man
	const std::string logname = IOUtils::getLogName();
	snprintf(user, MAX_STRING_LENGTH, "%s", logname.c_str());
}

/**
 * @brief Averages energy fluxes
 * @version 11.03
 * @param n_steps Number of calculation time steps since last output
 * @param useCanopyModel
 * @param Sdata
 * @param Xdata
*/
void averageFluxTimeSeries(const int& n_steps, const bool& useCanopyModel, SurfaceFluxes& Sdata, SnowStation& Xdata)
{
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
	if (Sdata.mA != Constants::undefined)
		Sdata.mA  /= n_steps;

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
 * TODO Adapt to new international code
 * @author Charles Fierz
 * @version 9.12
 * @param F1 Majority grain shape
 * @param F2 Minority grain shape
 * @param F3 2 indicates a melt-freeze crust
 * @param type aggregated shape information
 */
void typeToCode(int *F1, int *F2, int *F3, int type)
{
	*F1   = int (floor(type/100.));
	type -= int ((*F1)*100);
	*F2   = int (floor(type/10.));
	*F3   = int (type - (*F2)*10);
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
bool massBalanceCheck(const SnowStation& Xdata, const SurfaceFluxes& Sdata, double& tot_mass_in)
{
	bool mass_error = true;
	double tot_mass=0., tot_swe=0., dmassE=0.;
	const double hnw = Xdata.hn*Xdata.rho_hn;
	double mass_change = hnw - Sdata.mass[SurfaceFluxes::MS_RUNOFF] + Sdata.mass[SurfaceFluxes::MS_RAIN] + Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] + Sdata.mass[SurfaceFluxes::MS_EVAPORATION] - MAX(0., Xdata.ErosionMass);

	// Actual mass of snowpack
	for (unsigned int e=Xdata.SoilNode; e<Xdata.getNumberOfElements(); e++) {
		tot_mass += Xdata.Edata[e].M;
		tot_swe  += Xdata.Edata[e].L * Xdata.Edata[e].Rho;
		dmassE = Xdata.Edata[e].M - (Xdata.Edata[e].L * Xdata.Edata[e].Rho);
		if ( fabs(dmassE) > Constants::eps ) {
			prn_msg(__FILE__, __LINE__, "msg", Date(), "Mass error at element e=%d (nE=%d): mass(now)=%lf swe(now)=%lf dmassE=%lf", e, Xdata.getNumberOfElements(), tot_mass, tot_swe, dmassE);
			mass_error = false;
		}
	}
	// Mass balance check
	if ( tot_mass > Constants::eps ) {
		if ( (fabs(tot_swe/tot_mass) - 1.) > 0.5e-2 ) {
			prn_msg(__FILE__, __LINE__, "msg", Date(), "Mass balance (theta): mass(now)=%lf swe(now)=%lf swe/mass=%lf mass-swe=%lf", tot_mass, tot_swe, tot_swe/tot_mass, tot_mass - tot_swe);
			mass_error = false;
		}
		if ( tot_mass_in > Constants::eps ) {
			if ( (fabs(tot_mass - (tot_mass_in + mass_change))/tot_mass) > 0.5e-4 ) {
				prn_msg(__FILE__, __LINE__, "msg", Date(), "Mass balance: mass_err(now)=%lf", tot_mass - (tot_mass_in + mass_change));
				prn_msg(__FILE__, __LINE__, "msg", Date(), "tot_mass_in=%lf tot_mass_now=%lf mass_change=%lf", tot_mass_in, tot_mass, mass_change);
				mass_error = false;
			}
		} else {
			tot_mass_in = tot_mass;
		}
	} else if ( tot_mass_in > Constants::eps ) {
		if ( fabs(tot_mass_in + mass_change) > 1.0e-3 ) {
			prn_msg(__FILE__, __LINE__, "msg", Date(), "Mass balance error: tot_mass_in=%lf tot_mass_now=%lf mass_change=%lf", tot_mass_in, tot_mass, mass_change);
			mass_error = false;
		}
	}

	return mass_error;
}

/**
 * @brief Forced erosion of a missed event
 * @author Michael Lehning \n Mathis Bavay
 * @date 2008-03-07
 * @param hs Snow depth to assimilate
 * @param Xdata
 * @return Eroded mass (kg m-2)
 */
double forcedErosion(const double hs, SnowStation& Xdata)
{
	int    nErode=0;        // Counters
	double massErode=0.;    // Eroded mass (kg m-2)

	massErode=0.;
	while ( (Xdata.getNumberOfElements() > Xdata.SoilNode) && (hs + 0.01) < (Xdata.cH - Xdata.Ground) ) {
		massErode += Xdata.Edata[Xdata.getNumberOfElements()-1].M;
		Xdata.cH -= Xdata.Edata[Xdata.getNumberOfElements()-1].L;
		Xdata.resize(Xdata.getNumberOfElements() - 1);
		nErode++;
	}
	Xdata.ErosionLevel = MIN(Xdata.getNumberOfElements()-1, Xdata.ErosionLevel);

	return(massErode);
}


/**
 * @brief Deflates or inflates the snowpack for warning service purposes
 * - Forced erosion: dhs_corr < 0. & mass_corr > 0.
 * - Deflate       : dhs_corr < 0. & mass_corr < 0.
 * - Inlate        : dhs_corr > 0. & mass_corr > 0.
 * @note Michi is very unhappy that the warning service wants to deflate or inflate the
 *       snow if there is a consistent under- oder overestimation of settling, respectively.
 *       But since the will of the warning service is law at the SLF, we have no choice
 *       and need to implement this additional terrible non mass-conserving and cheating
 *       feature. But it will make the operational users happy, I hope. \n
 * Implemented on 2 Feb 2008 (and 7 Mar 2008: back to erosion) by Mathias Bavay
 * and Michi, who should be home with Leo being sick ...
 * @param Mdata
 * @param Xdata
 * @param dhs_corr Correction on calculated snow depth (m)
 * @param mass_corr Mass correction (kg m-2)
 */
void deflateInflate(const CurrentMeteo& Mdata, SnowStation& Xdata, double& dhs_corr, double& mass_corr)
{
	const size_t nE = Xdata.getNumberOfElements(), soil_node = Xdata.SoilNode;
	size_t e;
	double factor_corr, sum_total_correction=0.;  // Correction factors
	double ddL, dL = 0.;                          // Length changes
	const double cH = Xdata.cH - Xdata.Ground;    // Calculated snow depth
	const double mH = Xdata.mH - Xdata.Ground;    // Calculated snow depth
	double cH_old;                                // Temporary snow depth
	bool prn_CK = false;

	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;

	/*
	 * First try to find erosion events, which have not been captured by the drift module
	 * (Maybe the wind sensor did not measure correctly due to riming, or s.th. else went
	 * wrong with the model. For now assume erosion if more than 3 cm are missing
	 */
	if ((mH + 0.03) < cH) {
		dhs_corr = mH - cH;
		mass_corr = forcedErosion(mH, Xdata);
		if (prn_CK) { //HACK
			prn_msg(__FILE__, __LINE__, "msg+", Mdata.date, "Missed erosion event detected");
			prn_msg(__FILE__, __LINE__, "msg-", Date(), "Measured Snow Depth:%lf   Computed Snow Depth:%lf",
			        mH, cH);
		}
	} else if (cH > Constants::eps){ // assume settling error

		//Test whether normalization quantity does not lead to an arithmetic exception
		//This is a work around for weird cases in which the whole snowpack appears at once
		if (EMS[nE-1].depositionDate.getJulianDate() <= EMS[soil_node].depositionDate.getJulianDate())
			return;

		if (prn_CK) { //HACK
			prn_msg(__FILE__, __LINE__, "msg+", Mdata.date,
			          "Small correction due to assumed settling error");
			prn_msg(__FILE__, __LINE__, "msg-", Date(),
			          "True Measured Snow Depth:%lf   Computed Snow Depth:%lf", mH, cH);
		}
		// Second find the normalization quantity, which we choose to be the age of the layer.
		dhs_corr = mH - cH;
		for (e = soil_node; e < nE; e++) {
			if ((!(EMS[e].mk > 20 || EMS[e].mk == 3))
			        && (Mdata.date.getJulianDate() > EMS[e].depositionDate.getJulianDate())) {
				const double surf_date = EMS[nE-1].depositionDate.getJulianDate();
				const double current_layer_date = EMS[e].depositionDate.getJulianDate();
				const double first_snow_date = EMS[soil_node].depositionDate.getJulianDate();
				double age_fraction = (surf_date - current_layer_date) / (surf_date - first_snow_date);
				// Rounding errors could produce very small negative numbers ...
				if (age_fraction < 0.) {
					age_fraction = 0.;
				}
				sum_total_correction += EMS[e].L * (1. - sqrt(age_fraction) );
			}
		}
		if (sum_total_correction > 0.) {
			factor_corr = dhs_corr / sum_total_correction;
		} else {
			dhs_corr = 0.;
			return;
		}
		// ... above marked element (translation only) ...
		// Squeeze or blow-up
		for (e = soil_node; e < nE; e++) {
			if ((!(EMS[e].mk > 20 || EMS[e].mk == 3))
			        && (Mdata.date.getJulianDate() > EMS[e].depositionDate.getJulianDate())) {
				const double surf_date = EMS[nE-1].depositionDate.getJulianDate();
				const double current_layer_date = EMS[e].depositionDate.getJulianDate();
				const double first_snow_date = EMS[soil_node].depositionDate.getJulianDate();
				double age_fraction = (surf_date - current_layer_date) / (surf_date - first_snow_date);
				// Rounding errors could produce very small negative numbers ...
				if (age_fraction < 0.) {
					age_fraction = 0.;
				}
				ddL = EMS[e].L
				        * MAX(-0.9,
			                 MIN(0.9, factor_corr * (1. - sqrt(age_fraction))));
			} else {
				ddL = 0.;
			}
			dL += ddL;
			NDS[e+1].z += dL + NDS[e+1].u;
			NDS[e+1].u  = 0.0;
			mass_corr += ddL * EMS[e].Rho;
			EMS[e].M += ddL * EMS[e].Rho;
			EMS[e].L0 = EMS[e].L += ddL;
			EMS[e].E  = EMS[e].dE = EMS[e].Ee = EMS[e].Ev = EMS[e].S = 0.0;
		}
		// Update the overall height
		cH_old    = Xdata.cH;
		Xdata.cH  = NDS[nE].z + NDS[nE].u;
		//Xdata.mH -= (cH_old - Xdata.cH);
	} else {
		return;
	}
}

/**
 * @brief Logistic function
 * @version 10.11
 * @param input value
 * @param threshold Threshold of logistic function
 * @param width Width of logistic function
 */
double logisticFunction(const double input, const double threshold, const double width)
{
	const double x = (input - threshold) / width;
	return exp(x)/(1. + exp(x));
}

/**
 * @brief Cumulate if and only if value is defined
 * @version 11.01
 * @param accu cumulator
 * @param value to be added to accu
 */
void cumulate(double& accu, const double value)
{
	if (value != Constants::undefined) {
		accu += value;
	} else {
		accu = Constants::undefined;
	}

}
