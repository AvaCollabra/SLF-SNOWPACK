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
#include <snowpack/libsnowpack.h>
#include <meteoio/MeteoIO.h>

#include <iostream>
#include <string>
#include <sstream>

#ifdef _MSC_VER
	/*
	This software contains code under BSD license (namely, getopt for Visual C++).
	Therefore, this product includes software developed by the University of
	California, Berkeley and its contributors when compiling with Visual C++.
	*/
	#include "getopt.h"
#else
	//#include <unistd.h> //for getopt
	#include <getopt.h> //for getopt_long
#endif

using namespace std;
using namespace mio;

#ifdef DEBUG_ARITHM
	#ifndef _GNU_SOURCE
		#define _GNU_SOURCE
	#endif
	#ifndef __USE_GNU
		#define __USE_GNU
	#endif
	#include <fenv.h>
#endif


/**
 * @class Slope a C Fierz class ;-)
 * @brief . Peut-on l'être plus?
 * @version 11.01
 */
class Slope {

 public:
	 Slope(const mio::Config& i_cfg);

	 bool snow_redistribution;
	 bool virtual_slopes;
	 unsigned int nSlopes;
	 unsigned int station;      ///< main station, flat field or slope
	 unsigned int sector;       ///< current slope sector of width 360./MAX(1, nSlopes-1)
	 unsigned int first;        ///< first sector in computing sequence
	 unsigned int south;
	 unsigned int luv;
	 unsigned int lee;
	 double prevailing_wind_dir;

	 int getSectorDir(const double& dir_or_expo) const;
	 void setSlope(const int slope_sequence, vector<SnowStation>& vecXdata, double& wind_dir);

 private:
	 const mio::Config& cfg;
	 double sector_width;
};

/************************************************************
 * static section                                           *
 ************************************************************/

//Global variables in this file:
string cfgfile = "io.ini";
string mode = "RESEARCH";
mio::Date dateEnd;
vector<string> vecStationIDs;

/// @brief Main control parameters
struct MainControl
{
	double TimeN;      ///< Time of present computation (s)
	int    nStep;      ///< Time step number
	double Duration;   ///< Duration of run (s)
	int    TsDump;     ///< Flag for time series dump
	int    nAvg;       ///< Number of calculation time steps to average fluxes etc.
	int    HzStep;     ///< Hazard step number (should be half of nStep in operational mode)
	int    HzDump;     ///< Calculation of hazard information will be performed
	int    PrDump;     ///< Flag for profile dump
	int    XdataDump;  ///< Backup of Xdata will be performed
	int    TaglayDump; ///< Flag for tagged layer series dump
	int    PrTabDump;  ///< Flag for tabular profile dump
};

/************************************************************
 * non-static section                                       *
 ************************************************************/

Slope::Slope(const mio::Config& i_cfg) : cfg(i_cfg)
{
	cfg.getValue("NUMBER_SLOPES", "Snowpack", nSlopes);
	cfg.getValue("SNOW_REDISTRIBUTION", "Snowpack", snow_redistribution);
	virtual_slopes = (snow_redistribution && (nSlopes > 1) && (nSlopes % 2 == 1));

	station = 0;
	first = 1;
	cfg.getValue("PREVAILING_WIND_DIR", "SnowpackAdvanced", prevailing_wind_dir, mio::Config::nothrow);
	sector_width = 360./MAX(1, nSlopes-1);
	south = getSectorDir(180.);
}

/**
 * @brief Determine sector of either blowing wind or slope exposition\n
 *        NOTE that sector slope.first always corresponds to the prevailing wind direction
 * @param dir_or_expo direction of wind or exposition
 **/
int Slope::getSectorDir(const double& dir_or_expo) const
{
	double dir = dir_or_expo;
	if (dir > 360.) dir -= 360.;
	else if (dir < 0.) dir += 360.;
	unsigned int sectorDir = int (floor((dir + 0.5*sector_width)/sector_width));
	sectorDir++;
	if (sectorDir >= nSlopes) return 1;
	else return sectorDir;
}

/**
 * @brief Set slope variables
 * @param slope_sequence computation sequence for slopes
 * @param vecXdata
 * @param wind_dir direction of wind
 **/
void Slope::setSlope(const int slope_sequence, vector<SnowStation>& vecXdata, double& wind_dir)
{
	switch (slope_sequence) {
	case 0:
		for (unsigned int kk=0; kk<nSlopes; kk++) {
			vecXdata[kk].windward = false;
			vecXdata[kk].rho_hn   = 0.;
			vecXdata[kk].hn       = 0.;
		}
		if (nSlopes > 1) {
			luv = getSectorDir(wind_dir - prevailing_wind_dir);
			vecXdata[luv].windward = true;
			lee = (luv + nSlopes/2) % (nSlopes-1);
			if (lee == 0) lee = nSlopes - 1;
		} else {
			//requesting slope 0 of 0 expositions
			luv = lee = 0;
		}
		sector = station;
		break;
	case 1:
		sector = luv;
		break;
	default:
		sector++;
		if (sector == nSlopes) sector = 1;
	}
}

void Usage(const string& programname)
{
	cout << "Snowpack version " << _VERSION << " compiled on " << __DATE__ << " " << __TIME__ << endl
		<< "\tLibsnowpack " << snowpack::getLibVersion() << endl
		<< "\tMeteoIO " << mio::getLibVersion() << endl;
#ifdef _MSC_VER
	cout << "This version of Snowpack uses a BSD-licensed port of getopt for Visual C++. " << endl
		<< "It therefore includes software developed by the University of "
		<< "California, Berkeley and its contributors." << endl;
#endif
	cout << "Usage: " << programname << endl
		<< "\t-e, --enddate=YYYY-MM-DDTHH:MM (e.g.:2008-08-11T09:00)" << endl
		<< "\t[-c, --config=<ini file> (e.g. io.ini)]" << endl
		<< "\t[-m, --mode=[operational or research] (default: research)]" << endl
		<< "\t[-s, --stations=[comma delimited stationnames] (e.g. DAV2,WFJ2)] (NOTE: ONLY in operational mode)" << endl
		<< "\t[-h, --help] Print help message and version information]" << endl << endl;

	cout << "Example: " << programname << " -c io.ini -m research -e 1996-06-17T00:00" << endl << endl;
	exit(1);
}

void parseCmdLine(int argc, char **argv, string& end_date_str)
{
	int longindex=0, opt=-1;
	bool setEnd = false;

	struct option long_options[] =
	{
		{"enddate", required_argument, 0, 'e'},
		{"mode", required_argument, 0, 'm'},
		{"config", required_argument, 0, 'c'},
		{"stations", required_argument, 0, 's'},
		{"help", no_argument, 0, 'h'},
		{0, 0, 0, 0}
	};

	while ((opt=getopt_long( argc, argv, ":e:m:c:s:h", long_options, &longindex)) != -1) {
		switch (opt) {
		case 0:
			break;
		case 'e': {
			end_date_str=string(optarg); //we don't know yet the time zone, conversion will be done later
			setEnd = true;
			break;
		}
		case 'm':
			mode = string(optarg);
			mio::IOUtils::toUpper(mode);
			break;
		case 'c':
			cfgfile = string(optarg);
			break;
		case 's':
			mio::IOUtils::readLineToVec(string(optarg), vecStationIDs, ',');
			break;
		case ':': //operand missing
			cerr << endl << "ERROR: Command line parameter '-" << char(optopt) << "' requires an operand" << endl;
			Usage(string(argv[0]));
			break;
		case 'h':
			Usage(string(argv[0]));
			break;
		case '?':
			cerr << endl << "ERROR: Unknown argument detected" << endl;
			Usage(string(argv[0]));
			break;
		default:
			cerr << endl << "ERROR: getopt returned character code " <<  opt << endl;
			Usage(string(argv[0]));
		}
	}

	if (!setEnd) {
		cerr << endl << "ERROR: You must specify an enddate for the simulation!" << endl;
		Usage(string(argv[0]));
	}
}

void editMeteoData(mio::MeteoData& md, const string& variant)
{
	// Since we cannot deal with precipitation nodata, we set it to zero (HACK)
	if (md(MeteoData::HNW) == mio::IOUtils::nodata)
		md(MeteoData::HNW) = 0.0;

	if (md(MeteoData::VW) == mio::IOUtils::nodata)
		md(MeteoData::VW) = 1.0; // if no wind measurement exists assume 1 m/s; ori: 3 m/s

	if (md(MeteoData::DW) == mio::IOUtils::nodata)
		md(MeteoData::DW) = 0.;

	if (md(MeteoData::TSG) == mio::IOUtils::nodata)
		md(MeteoData::TSG) = 273.15;

	//Add the atmospheric emissivity as a parameter
	if (!md.param_exists("EA")) md.addParameter("EA");
	if (variant != "ANTARCTICA") {
		md("EA") = lw_AirEmissivity(md(MeteoData::ILWR), md(MeteoData::TA), md(MeteoData::RH));
	} else {
		//change min_air_emissivity
		md("EA") = lw_AirEmissivity(md(MeteoData::ILWR), md(MeteoData::TA), md(MeteoData::RH), 0.31);
	}

	// Snow stations without separate wind station use their own wind for drifting and blowing snow
	if (!md.param_exists("VW_DRIFT")) {
		md.addParameter("VW_DRIFT");
		md("VW_DRIFT") = md(MeteoData::VW);
	}
	if (!md.param_exists("DW_DRIFT")) {
		md.addParameter("DW_DRIFT");
		md("DW_DRIFT") = md(MeteoData::DW);
	}
}

// Return true if snowpack can compute the next timestep, else false
bool validMeteoData(const mio::MeteoData& md, const string& StationName, const string& variant)
{
	bool miss_ta=false, miss_rh=false, miss_precip=false, miss_rad=false;
	bool miss_ea=false;

	if (md(MeteoData::TA) == mio::IOUtils::nodata)
		miss_ta=true;
	if (md(MeteoData::RH) == mio::IOUtils::nodata)
		miss_rh=true;
	if ((variant != "ANTARCTICA")
	        && ((md(MeteoData::ISWR) == mio::IOUtils::nodata) && (md(MeteoData::RSWR) == mio::IOUtils::nodata)))
		miss_rad=true;
	if ((md(MeteoData::HNW) == mio::IOUtils::nodata) && (md(MeteoData::HS) == mio::IOUtils::nodata))
		miss_precip=true;
	if (md("EA") == mio::IOUtils::nodata)
		miss_ea=true;

	if(miss_ta || miss_rh || miss_rad || miss_precip || miss_ea) {
		mio::Date now;
		now.setFromSys();
		cout << "[E] [" << now.toString(mio::Date::ISO) << "] ";
		cout << StationName << " missing { ";
		if(miss_ta) cout << "TA ";
		if(miss_rh) cout << "RH ";
		if(miss_rad) cout << "radiation ";
		if(miss_precip) cout << "precipitations ";
		if(miss_ea) cout << "ea ";
		cout << "} on [" << md.date.toString(mio::Date::ISO) << "]\n";
		return false;
	}
	return true;
}

void copyMeteoData(const mio::MeteoData& md, CurrentMeteo& Mdata, const double prevailing_wind_dir, const double wind_scaling_factor)
{
	Mdata.date = md.date;
	Mdata.ta   = md(MeteoData::TA);
	Mdata.rh   = md(MeteoData::RH);
	Mdata.vw   = md(MeteoData::VW);
	Mdata.dw   = md(MeteoData::DW);
	Mdata.vw_max  = md(MeteoData::VW_MAX);

	Mdata.vw_drift = md("VW_DRIFT");
	if (Mdata.vw_drift != mio::IOUtils::nodata) Mdata.vw_drift *= wind_scaling_factor;
	Mdata.dw_drift = md("DW_DRIFT");
	if (Mdata.dw_drift == mio::IOUtils::nodata) Mdata.dw_drift = prevailing_wind_dir;

	if (md(MeteoData::ISWR) == mio::IOUtils::nodata) {
		Mdata.iswr   = 0.0;
	} else {
		Mdata.iswr   = md(MeteoData::ISWR);
	}

	if (md(MeteoData::RSWR) == mio::IOUtils::nodata){
		Mdata.rswr   = 0.0;
	} else {
		Mdata.rswr   = md(MeteoData::RSWR);
	}

	Mdata.ea  = md("EA");
	Mdata.tss = md(MeteoData::TSS);
	Mdata.ts0 = md(MeteoData::TSG);
	Mdata.hnw = md(MeteoData::HNW);
	Mdata.hs1 = md(MeteoData::HS);

	// Add measured new snow density if available
	if (md.param_exists("RHO_HN"))
		Mdata.rho_hn = md("RHO_HN");
}

void editSensorDepths(const mio::MeteoData& md, vector<double>& vecHTS)
{
	unsigned int jj;
	for (jj = 0; jj < vecHTS.size(); jj++) {
		stringstream ss;
		ss << "HTS" << jj+1;
		if (md.param_exists(ss.str()) && (md(ss.str()) != mio::IOUtils::nodata)) {
			vecHTS[jj] = md(ss.str());
		} else {
			if (jj > 0)
				vecHTS[jj] = vecHTS[jj-1] + 0.5;
			break;
		}
	}
}

void copySnowTemperatures(const mio::MeteoData& md, CurrentMeteo& Mdata, vector<double>& vecHTS, const int current_slope)
{
	for (unsigned int jj=0; jj < vecHTS.size(); jj++) {
		Mdata.zv_ts[jj] = vecHTS[jj];
		Mdata.ts[jj] = mio::IOUtils::nodata;
		if (current_slope == 0) {
			stringstream ss;
			ss << "TS" << jj+1;
			if (md.param_exists(ss.str()) && (md(ss.str()) != mio::IOUtils::nodata)) {
				Mdata.ts[jj] = md(ss.str());
			}
		}
	}
}

void copySolutes(const mio::MeteoData& md, CurrentMeteo& Mdata, const unsigned int& i_number_of_solutes)
{
	if (i_number_of_solutes > 0) {
		for (unsigned int jj=0; jj < i_number_of_solutes; jj++) {
			Mdata.conc[jj] = mio::IOUtils::nodata;
			stringstream ss;
			ss << "CONC" << jj;
			Mdata.conc[jj] = md(ss.str());
		}
	} else {
		return;
	}
}

//for a given config (that can be altered) and original meteo data, prepare the snowpack data structures
//This means that all tweaking of config MUST be reflected in the config object
void dataForCurrentTimeStep(CurrentMeteo& Mdata, SurfaceFluxes& surfFluxes, vector<SnowStation>& vecXdata,
                            const Slope& slope, SnowpackConfig& cfg,
                            PositionSun& Psolar, RadiationData& Rdata,
                            double& cumu_hnw, const double& lw_in, double& iswr_forced, double& tot_mass_in)
{
	/**
	 * @brief Switch defining whether input data (solar radiation, snow depth and precipitation rates)
	 * is from EB\@ALPINE3D \n
	 * Set to false in stand-alone applications.
	 */
	const bool ebalance_switch = false;

	const bool useCanopyModel = cfg.get("CANOPY", "Snowpack");
	const bool enforce_measured_snow_heights = cfg.get("ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack");
	const bool incoming_longwave = cfg.get("INCOMING_LONGWAVE", "Snowpack");
	int sw_mode = 0;
	cfg.getValue("SW_MODE", "Snowpack", sw_mode);
	sw_mode %= 10;
	if (Mdata.tss==mio::IOUtils::nodata) { //degraded, but better than nothing if no TSS!
		cfg.addKey("CHANGE_BC", "Snowpack", "0");
		cfg.addKey("MEAS_TSS", "Snowpack", "0");
		Mdata.tss=Constants::melting_tk; //TODO: it would be better to test for NEUMANN instead of relying on tss=melting_tk in the model...
	}

	const bool sw_mode_change = cfg.get("SW_MODE_CHANGE", "SnowpackAdvanced"); //Adjust for correct radiation input if ground is effectively bare. It has to be set to true in operational mode.
	const bool mass_balance = cfg.get("MASS_BALANCE", "SnowpackAdvanced");
	const bool perp_to_slope = cfg.get("PERP_TO_SLOPE", "SnowpackAdvanced");

	const bool avgsum_time_series = cfg.get("AVGSUM_TIME_SERIES", "Output");
	const bool cumsum_mass = cfg.get("CUMSUM_MASS", "Output");

	Radiation radiation(cfg);

	// Reset Surface and Canopy Data to zero if you seek current values
	if (!avgsum_time_series) {
		surfFluxes.reset(cumsum_mass);
		if (useCanopyModel)
			vecXdata[slope.sector].Cdata.reset(cumsum_mass);

		if (mass_balance) { // Do an initial mass balance check
			tot_mass_in = 0.;
			if (!massBalanceCheck(vecXdata[slope.sector], surfFluxes, tot_mass_in))
				prn_msg(__FILE__, __LINE__, "msg+", Mdata.date, "Mass error during initial check!");
		}
	}
	// Reset surfFluxes.drift and surfFluxes.mass[MS_WIND] anyway since you use these to transport eroded snow
	surfFluxes.drift = 0.;
	surfFluxes.mass[SurfaceFluxes::MS_WIND] = 0.;

	if (slope.sector == slope.station) {
		//split flat field radiation and compute potential radiation
		radiation.flatFieldRadiation(vecXdata[slope.station], Mdata, Psolar, Rdata);
		iswr_forced = -1.0; //initialize on station field
		if (sw_mode_change) {
			/*
			 * Sometimes, there is no snow left on the ground at the station (-> rswr is small)
			 * but there is still some snow left in the simulation, which then is hard to melt
			 * if we find this is such a situation, we set iswr to the potential radiation.
			 * Such a correction is only needed for flat field, the others will inherit it
			*/
			double hs, iswr_factor;
			// What snow depth should be used?
			if (enforce_measured_snow_heights &&
				(vecXdata[slope.station].meta.getSlopeAngle() < Constants::min_slope_angle)) {
				hs = Mdata.hs1;
			} else {
				hs = vecXdata[slope.station].cH;
			}
			iswr_factor=(Mdata.rswr+0.01)/(Rdata.pot_dir+Rdata.pot_diffsky+0.00001); //avoiding "0/0"
			if (hs<0.1 && Mdata.rh<0.7 && iswr_factor<0.3) {
				Rdata.dir_hor = Rdata.pot_dir;
				Rdata.diffsky = Rdata.pot_diffsky;
				Rdata.global_hor = Rdata.dir_hor + Rdata.diffsky;

				iswr_forced = Rdata.dir_hor + Rdata.diffsky;
				cfg.addKey("SW_MODE", "Snowpack", "2");  // as both Mdata.iswr and Mdata.rswr were reset
				sw_mode = 2;
			}
		}
		// Compute running mean over nHours window
		Mdata.vw_ave = Mdata.vw; //TODO running mean over +-50 hours required
		Mdata.rh_ave = Mdata.rh; //TODO running mean over +-50 hours required

	} else { //virtual slope
		cfg.addKey("CHANGE_BC", "Snowpack", "0");
		cfg.addKey("ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack", "1");
		cfg.addKey("MEAS_TSS", "Snowpack", "0");
		cfg.addKey("NUMBER_MEAS_TEMPERATURES", "Input", "0");
		Mdata.tss=Constants::melting_tk; //HACK: should we write tss=melting_tk here?
	}

	if (iswr_forced >= 0.) {
		//iswr_forced=-1 when starting on flat field. If iswr was recomputed, then iswr_forced>=0
		Mdata.iswr = iswr_forced;
		if ((Mdata.rswr/Mdata.iswr) < (2.0*vecXdata[slope.sector].SoilAlb))
			Mdata.rswr = Mdata.iswr*2.0 * vecXdata[slope.sector].SoilAlb;
	}

	Meteo meteo(cfg);

	// Project irradiance on slope; take care of measured snow depth and/or precipitations too
	if (!(ebalance_switch || perp_to_slope)) {
		radiation.radiationOnSlope(vecXdata[slope.sector], Mdata, surfFluxes, Psolar, Rdata);

		if (((sw_mode == 1) || (sw_mode == 2)) && (vecXdata[slope.sector].meta.getSlopeAngle() > Constants::min_slope_angle)) {
			cfg.addKey("SW_MODE", "Snowpack", "0"); // as Mdata.iswr is the sum of dir_slope and diff
		}

		//HACK A3D: resynchronize Mdata with the solar and radiation parameters
		Mdata.elev = Psolar.elev;
		Mdata.diff = Rdata.diffsky;
		meteo.projectPrecipitations(vecXdata[slope.sector].meta.getSlopeAngle(), Mdata.hnw, Mdata.hs1);
	}

	// Find the Wind Profile Parameters, w/ or w/o canopy; take care of canopy
	meteo.compMeteo(&Mdata, &vecXdata[slope.sector]);

	// Update precipitation memory
	if (slope.sector == slope.station) {
		cumu_hnw += Mdata.hnw;
		if (enforce_measured_snow_heights && (Mdata.hs1 != mio::IOUtils::nodata)) {
			vecXdata[slope.station].mH = Mdata.hs1 + vecXdata[slope.station].Ground;
		}
	}

	if (slope.sector != slope.station) { // Meteo data for virtual slope simulations
		double hn_slope = 0., rho_hn_slope = SnLaws::min_hn_density;
		// Compute new snow water equivalent and density
		if (vecXdata[slope.station].hn > 0.) {
			// Assign new snow depth and density from station field (usually flat)
			hn_slope = vecXdata[slope.station].hn * cos(DEG_TO_RAD(vecXdata[slope.sector].meta.getSlopeAngle()));
			rho_hn_slope = vecXdata[slope.station].rho_hn;
		}

		/*
		 * Snow redistribution on slopes: Add windward eroded snow to lee slope
		 * These are very important lines: Note that deposition is treated here (lee)
		 * while erosion is treated in SnowDrift.c (windward).
		*/
		if (slope.snow_redistribution && (slope.sector == slope.lee)) {
			// If it is not snowing, use surface snow density on windward slope
			if (!(hn_slope > 0.)) {
				rho_hn_slope = vecXdata[slope.luv].rho_hn;
			}
			// Add eroded mass from windward slope
			if (rho_hn_slope != 0.) {
				hn_slope += vecXdata[slope.luv].ErosionMass / rho_hn_slope;
			}
			vecXdata[slope.luv].ErosionMass = 0.;
		}

		// Check whether snow depth on slope needs to be changed
		if ((hn_slope > 0.) && (vecXdata[slope.station].cH > 0.01)) {
			// Increase snow depth
			vecXdata[slope.sector].hn  = hn_slope;
			vecXdata[slope.sector].rho_hn = rho_hn_slope;
		}

		// Check whether to use incoming longwave as estimated from station field
		if (!incoming_longwave) {
			double ea = 0.;
			if ((lw_in > 50.) && (lw_in < 300.)) {
				ea = lw_in/Mdata.ta/Mdata.ta/Mdata.ta/Mdata.ta;
				ea /= Constants::stefan_boltzmann;
			}
			if ((ea > 0.55) && (ea < 1.0)) {
				Mdata.ea = ea;
			}
		}
	}
}

/**
 * @brief Get 6 and 24 hour difference in snow depth
 * @version 11.02
 * @param delta_hs6 Snow depth difference over 6 hours
 * @param delta_hs24 Snow depth difference over 24 hours
 * @param io
 * @param md meteo data
 * @param i_stn station number
 * @param current_date
 */
void getDhs6Dhs24(double& delta_hs6, double& delta_hs24, mio::IOManager& io, const vector<mio::MeteoData>& md, const int i_stn, const mio::Date current_date)
{
	vector<mio::MeteoData> md6, md24; // meteo data 6 and 24 hours ago, repectively

	delta_hs6  = Constants::undefined;
	delta_hs24 = Constants::undefined;

	try {
		io.getMeteoData(current_date - 1.0, md24);
		io.getMeteoData(current_date - 6.0/24.0, md6);
	} catch (...) {
		cerr << "[E] failed to read meteo data for previous 6 and 24 hours\n";
		throw;
	}

	if ((md6.size() > 0) && (md[i_stn](MeteoData::HS) != mio::IOUtils::nodata) && (md6[i_stn](MeteoData::HS) != mio::IOUtils::nodata)) {
		delta_hs6 = MAX(0., md[i_stn](MeteoData::HS) - md6[i_stn](MeteoData::HS));
	}
	if ((md24.size() > 0) && (md[i_stn](MeteoData::HS) != mio::IOUtils::nodata) && (md24[i_stn](MeteoData::HS) != mio::IOUtils::nodata)) {
		delta_hs24 = MAX(0., md[i_stn](MeteoData::HS) - md24[i_stn](MeteoData::HS));
	}
}
/**
 * @brief determine which outputs need to be done for the current time step
 * @param step current time integration step
 * @param sno_step current step in the sno files (current sno profile)
 * @param mn_ctrl timestep control structure
 **/

void getOutputControl(MainControl& mn_ctrl, const mio::Date& step, const mio::Date& sno_step, const double& calculation_step_length,
					  const double& tsstart, const double& tsdaysbetween,
	   const double& profstart, const double& profdaysbetween,
	const double& first_backup, const double& backup_days_between)
{
//HACK: put all tsstart, tsdaysbetween, etc in MainControl as well as current timestep
	const double Dstep = step.getJulianDate();
	const double Dsno_step = sno_step.getJulianDate();
	if ( mn_ctrl.nStep ) {
		// Hazard data, every half-hour
		mn_ctrl.HzDump = booleanTime(Dstep, 0.5/24., 0.0, calculation_step_length);
		// Time series (*.met)
		double bool_start = H_TO_D(tsstart);
		if ( bool_start > 0. )
			bool_start += Dsno_step;
		mn_ctrl.TsDump = booleanTime(Dstep, tsdaysbetween, bool_start, calculation_step_length);
		// Profile (*.pro)
		bool_start = H_TO_D(profstart);
		if ( bool_start > 0. )
			bool_start += Dsno_step;
		mn_ctrl.PrDump = booleanTime(Dstep, profdaysbetween, bool_start, calculation_step_length);

	} else {
		mn_ctrl.HzDump = 0;
		mn_ctrl.TsDump = 1;
		mn_ctrl.PrDump = 1;
	}

	// Additional Xdata backup (*.<JulianDate>sno)
	const double bool_start = Dsno_step + first_backup;
	mn_ctrl.XdataDump = booleanTime(Dstep, backup_days_between, bool_start, calculation_step_length);
}

// SNOWPACK MAIN **************************************************************
int main (int argc, char *argv[])
{
	if (argc==1) Usage(string(argv[0]));

#ifdef DEBUG_ARITHM
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); //for halting the process at arithmetic exceptions
	//feenableexcept(FE_ALL_EXCEPT);
#endif
	mio::Timer prebuffering_timer;
	mio::Timer meteoRead_timer;
	mio::Timer run_timer;
	run_timer.start();
	time_t nowSRT = time(NULL);
	MainControl mn_ctrl; //Time step control parameters

	string end_date_str;
	parseCmdLine(argc, argv, end_date_str); //parse the command line arguments

	SnowpackConfig cfg(cfgfile);

	const double i_time_zone = cfg.get("TIME_ZONE", "Input"); //get user provided input time_zone
	if (end_date_str == "NOW") { //interpret user provided end date
		dateEnd.setFromSys();
		dateEnd.setTimeZone(i_time_zone);
		dateEnd.rnd(1800., mio::Date::DOWN);
	} else {
		mio::IOUtils::convertString(dateEnd, end_date_str, i_time_zone);
	}

	size_t nStations = 0;
	string outpath(""), experiment("");
	stringstream ss;

	if (vecStationIDs.size() == 0) {
		cfg.getValue("NROFSTATIONS", "Input", nStations);
		vecStationIDs.clear();
		for (size_t i_stn=0; i_stn<nStations; i_stn++) {
			string stationID("");
			ss.str("");
			ss << "STATION" << i_stn+1;
			cfg.getValue(ss.str(), "Input", stationID);
			vecStationIDs.push_back(stationID);
		}
	} else {
		nStations = vecStationIDs.size();
		ss.str("");
		ss << nStations;
		cfg.addKey("NROFSTATIONS", "Input", ss.str());
		for (size_t i_stn=0; i_stn<nStations; i_stn++) {
			ss.str("");
			ss << "STATION" << i_stn+1;
			cfg.addKey(ss.str(), "Input", vecStationIDs[i_stn]);
		}
	}

	for (size_t i_stn=0; i_stn<vecStationIDs.size(); i_stn++) {
		string meteosrc; cfg.getValue("METEO", "Input", meteosrc);
		ss.str("");
		ss << "METEOFILE" << i_stn+1;
		string meteofile; cfg.getValue(ss.str(), "Input", meteofile, mio::Config::nothrow);
		if (meteofile == "") {
			if (meteosrc == "SMET") {
				cfg.addKey(ss.str(), "Input", vecStationIDs[i_stn]+".smet");
			} else if (meteosrc == "SNOWPACK") {
				cfg.addKey(ss.str(), "Input", vecStationIDs[i_stn]+".inp");
			}
		}
	}

	string variant; cfg.getValue("VARIANT", "SnowpackAdvanced", variant, mio::Config::nothrow);
	const bool useSoilLayers = cfg.get("SNP_SOIL", "Snowpack");
	bool soil_flux = false;
	if (useSoilLayers) {
		cfg.getValue("SOIL_FLUX", "Snowpack", soil_flux);
		prn_msg(__FILE__, __LINE__, "msg",  mio::Date(), "Start SNOWPACK w/ soil layers in %s mode", mode.c_str());
	} else {
		prn_msg(__FILE__, __LINE__, "msg",  mio::Date(), "Start SNOWPACK in %s mode", mode.c_str());
	}
	if (variant != "DEFAULT") {
		prn_msg(__FILE__, __LINE__, "msg",  mio::Date(), "Variant is '%s'", variant.c_str());
	}
	prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "%s compiled on %s at %s", string(argv[0]).c_str(), __DATE__, __TIME__);

	const bool useCanopyModel = cfg.get("CANOPY", "Snowpack");
	if (mode == "OPERATIONAL") {
		cfg.addKey("RESEARCH", "SnowpackAdvanced", "0");
		cfg.addKey("AVGSUM_TIME_SERIES", "Output", "0");
		if (useCanopyModel) {
			throw mio::IOException("Please don't set CANOPY to 1 in OPERATIONAL mode", AT);
		}
	} else {
		cfg.getValue("EXPERIMENT", "Output", experiment, mio::Config::nothrow);
		cfg.getValue("METEOPATH", "Output", outpath);
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Experiment : %s", experiment.c_str());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Output dir : %s", outpath.c_str());
	}

	const double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	const double sn_dt = M_TO_S(calculation_step_length); //Calculation time step in seconds

	const unsigned int max_number_sensors = cfg.get("MAX_NUMBER_SENSORS", "SnowpackAdvanced", mio::Config::nothrow);
	const unsigned int number_meas_temperatures = cfg.get("NUMBER_MEAS_TEMPERATURES", "Input", mio::Config::nothrow);
	vector<double> fixed_sensor_depths;
	if (number_meas_temperatures > 0) {
		cfg.getValue("FIXED_SENSOR_DEPTHS", "Input", fixed_sensor_depths);
	} else {
		cfg.getValue("FIXED_SENSOR_DEPTHS", "Output", fixed_sensor_depths);
	}

	int nSolutes = Constants::inodata;
	cfg.getValue("NUMBER_OF_SOLUTES", "Input", nSolutes, mio::Config::nothrow);
	if (nSolutes > 0) SnowStation::number_of_solutes = unsigned(nSolutes);

	//Interval between profile backups (*.sno\<JulianDate\>) (d)
	double backup_days_between = 400.;
	cfg.getValue("BACKUP_DAYS_BETWEEN", "Output", backup_days_between, mio::Config::nothrow);
	//First additional profile backup (*.sno\<JulianDate\>) since start of simulation (d)
	double first_backup = 0.;
	cfg.getValue("FIRST_BACKUP", "Output", first_backup, mio::Config::nothrow);

	//To check mass balance if AVGSUM_TIME_SERIES is not set (screen output only)
	bool mass_balance = false;
	cfg.getValue("MASS_BALANCE", "Output", mass_balance, mio::Config::nothrow);

	const bool profwrite = cfg.get("PROF_WRITE", "Output");
	const double profstart = cfg.get("PROF_START", "Output");
	const double profdaysbetween = cfg.get("PROF_DAYS_BETWEEN", "Output");
	const bool tswrite = cfg.get("TS_WRITE", "Output");
	const double tsstart = cfg.get("TS_START", "Output");
	const double tsdaysbetween = cfg.get("TS_DAYS_BETWEEN", "Output");
	const int hazard_steps_between = cfg.get("HAZARD_STEPS_BETWEEN", "Output");

	const bool precip_rates = cfg.get("PRECIP_RATES", "Output", mio::Config::nothrow);
	const bool avgsum_time_series = cfg.get("AVGSUM_TIME_SERIES", "Output", mio::Config::nothrow);
	const bool cumsum_mass = cfg.get("CUMSUM_MASS", "Output", mio::Config::nothrow);

	SnowpackIO snowpackio(cfg);

	mio::IOManager io(cfg);
	io.setMinBufferRequirements(IOUtils::nodata, 1.1); //we require the buffer to contain at least 1.1 day before the current point

	/* START LOOP OVER ALL STATIONS */
	for (size_t i_stn=0; i_stn<vecStationIDs.size(); i_stn++) {
		cout << endl;
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Run on station %s", vecStationIDs[i_stn].c_str());
		run_timer.reset();
		meteoRead_timer.reset();
		prebuffering_timer.reset();

		Slope slope(cfg);
		if (slope.snow_redistribution && !((slope.nSlopes != 1) || (slope.nSlopes != 5) || (slope.nSlopes != 9)))
			throw mio::IOException("Please set NUMBER_SLOPES to 1, 5 or 9", AT);

		double lw_in = Constants::undefined;    // Storage for LWin from flat field energy balance
		double iswr_forced = -1; //if iswr was wrongly evaluated, this contains the incoming pot. radiation

		// Used to scale wind for blowing and drifting snowpack (from statistical analysis)
		double wind_scaling_factor = cfg.get("WIND_SCALING_FACTOR", "SnowpackAdvanced");

		// Control of time window: used for adapting diverging snow depth in operational mode
		double time_count_deltaHS = 0.;

		// To sum precips and erosion mass if !AVGSUM_TIME_SERIES
		double cumu_hnw = 0., cumsum_drift = 0., cumsum_snow = 0., cumsum_runoff = 0., cumsum_rain = 0.;
		vector<double> cumsum_erosion(slope.nSlopes, 0.); // Cumulated eroded mass; dumped to file as rate

		// Snowpack data (input/output)
		vector<SN_SNOWSOIL_DATA> vecSSdata(slope.nSlopes, SN_SNOWSOIL_DATA(/*number_of_solutes*/));
		vector<SnowStation> vecXdata(slope.nSlopes, SnowStation(useCanopyModel, useSoilLayers/*, number_of_solutes*/));
		SN_ZWISCHEN_DATA sn_Zdata;   // "Memory"-data, required for every operational station

		// Meteo data for the current time step (interpolated!)
		CurrentMeteo sn_MdataT(max_number_sensors/*, number_of_solutes*/);
		// To collect surface exchange data for output
		SurfaceFluxes surfFluxes/*(number_of_solutes)*/;
		// Boundary condition (fluxes)
		BoundCond sn_Bdata;

		if (mode == "OPERATIONAL") cfg.addKey("PERP_TO_SLOPE", "SnowpackAdvanced", "0");

		mio::Date current_date;
		string snowfile("");
		ss.str("");
		ss << "SNOWFILE" << i_stn+1;
		cfg.getValue(ss.str(), "Input", snowfile, mio::Config::nothrow);
		for (size_t sector=slope.station; sector<slope.nSlopes; sector++) { //Read SSdata for every sector
			try {
				if (sector == slope.station) {
					if (snowfile == "") {
						snowfile = vecStationIDs[i_stn];
					} else {
						const size_t pos_dot = snowfile.rfind(".");
						const size_t pos_slash = snowfile.rfind("/");
						if( pos_dot!= string::npos && pos_dot>pos_slash) //so that the dot is not in a directory name
							snowfile.erase(pos_dot, snowfile.size()-pos_dot);
					}
					snowpackio.readSnowCover(snowfile, vecStationIDs[i_stn], vecSSdata[slope.station], sn_Zdata);
					prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Reading snow cover data for station %s",
					            vecStationIDs[i_stn].c_str());
					// NOTE (Is it a HACK?) Reading station meta data provided in meteo data and prebuffering those data
					vector<mio::MeteoData> vectmpmd;
					current_date = vecSSdata[slope.station].profileDate;
					prebuffering_timer.start();
					io.getMeteoData(current_date, vectmpmd);
					prebuffering_timer.stop();
					if (vectmpmd.size() == 0)
						throw mio::IOException("No data found for station " + vecStationIDs[i_stn] + " on "
						                           + current_date.toString(mio::Date::ISO), AT);
					editSensorDepths(vectmpmd[i_stn], fixed_sensor_depths);
					vecSSdata[slope.station].meta = mio::StationData::merge(vectmpmd[i_stn].meta,
					                                    vecSSdata[slope.station].meta);
				} else {
					stringstream sec_snowfile;
					sec_snowfile << snowfile << sector;
					ss.str("");
					ss << vecSSdata[slope.station].meta.getStationID() << sector;
					snowpackio.readSnowCover(sec_snowfile.str(), ss.str(), vecSSdata[sector], sn_Zdata);
					vecSSdata[sector].meta.position = vecSSdata[slope.station].meta.getPosition();
					vecSSdata[sector].meta.stationName = vecSSdata[slope.station].meta.getStationName();
				}
				vecXdata[sector].initialize(vecSSdata[sector], sector); // Generate the corresponding Xdata
			} catch (const exception& e){
				if (sector == slope.first){
					prn_msg(__FILE__, __LINE__, "msg-", mio::Date(),
					            "No virtual slopes! Computation for main station %s only!",
					                vecStationIDs[i_stn].c_str());
					slope.nSlopes = 1;
					if ((mode == "OPERATIONAL")
					        && (vecSSdata[slope.station].meta.getSlopeAngle() > Constants::min_slope_angle)) {
						cfg.addKey("PERP_TO_SLOPE", "SnowpackAdvanced", "1");
					}
					break;
				} else {
					cout << e.what();
					throw;
				}
			}

			// Operational mode ONLY: Pass ... wind_factor, snow depth discrepancy time counter
			if ((sector == slope.station) && (mode == "OPERATIONAL")) {
				wind_scaling_factor = vecSSdata[slope.station].WindScalingFactor;
				time_count_deltaHS = vecSSdata[slope.station].TimeCountDeltaHS;
			}
		}
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Finished Initializing station %s", vecStationIDs[i_stn].c_str());

		//CHECK date inconsistencies between sno files
		bool dates_consistent(true);
		for (size_t sector=slope.first; sector<slope.nSlopes; sector++) {
			if (vecSSdata[sector].profileDate != vecSSdata[slope.station].profileDate) {
				prn_msg(__FILE__, __LINE__, "err", mio::Date(),
				        "Date virtual slope %d inconsistent with flat field %s", sector, vecStationIDs[i_stn].c_str());
				dates_consistent = false;

			}
		}
		if (!dates_consistent) continue; //go to next station

		// Do not go ahead if starting time is larger than maxtime!
		if (vecSSdata[slope.station].profileDate > dateEnd) {
			prn_msg(__FILE__, __LINE__, "err", mio::Date(), "Starting time (%.5lf) larger than end time(%.5lf), station %s!",
			        vecSSdata[slope.station].profileDate.getJulianDate(), dateEnd.getJulianDate(), vecStationIDs[i_stn].c_str());
			continue; //goto next station
		}

		memset(&mn_ctrl, 0, sizeof(MainControl));
		if (mode == "RESEARCH") {
			deleteOldOutputFiles(outpath, experiment, vecStationIDs[i_stn], slope.nSlopes);
			cfg.write(outpath + "/" + vecStationIDs[i_stn] + "_" + experiment + ".ini"); //output config
			current_date -= calculation_step_length/1440;
			mn_ctrl.nStep = -1;
		}

		mn_ctrl.Duration = (dateEnd.getJulianDate() - vecSSdata[slope.station].profileDate.getJulianDate() + 0.5/24)*24*3600;
		vector<ProcessDat> qr_Hdata;     //Hazard data for t=0...tn
		vector<ProcessInd> qr_Hdata_ind; //Hazard data Index for t=0...tn
		Hazard hazard(cfg, mn_ctrl.Duration);
		hazard.initializeHazard(sn_Zdata.drift24, vecXdata.at(0).meta.getSlopeAngle(), qr_Hdata, qr_Hdata_ind);

		prn_msg(__FILE__, __LINE__, "msg", mio::Date(), "Start simulation for %s on %s (%f) station time",
		        vecStationIDs[i_stn].c_str(), vecSSdata[slope.station].profileDate.toString(mio::Date::ISO).c_str(),
		        vecSSdata[slope.station].profileDate.getJulianDate());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "End date specified by user: %s", dateEnd.toString(mio::Date::ISO).c_str());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Integration step length: %f min", calculation_step_length);

		bool computed_one_timestep = false;

		// START TIME INTEGRATION LOOP
		do {
			current_date += calculation_step_length/1440;
			mn_ctrl.nStep++;
			mn_ctrl.nAvg++;

			// Get meteo data
			vector<mio::MeteoData> vecMyMeteo;
			meteoRead_timer.start();
			io.getMeteoData(current_date, vecMyMeteo);
			meteoRead_timer.stop();
			editMeteoData(vecMyMeteo[i_stn], variant);
			if (!validMeteoData(vecMyMeteo[i_stn], vecStationIDs[i_stn], variant)) {
				prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "No valid data for station %s on [%s]",
				        vecStationIDs[i_stn].c_str(), current_date.toString(mio::Date::ISO).c_str());
				break;
			}

			//determine which outputs will have to be done
			getOutputControl(mn_ctrl, current_date, vecSSdata[slope.station].profileDate, calculation_step_length,
			                 tsstart, tsdaysbetween, profstart, profdaysbetween,
			                 first_backup, backup_days_between);

			//Radiation data
			RadiationData Rdata;   // Radiation splitting data
			PositionSun   Psolar;  // Parameters to determine the position of the sun

			//Calculate average TSS and HS for first snow fall detection algorithm
			Meteo::compTSSavgHSrate(sn_MdataT, vecXdata[i_stn], io, current_date);

			// START LOOP OVER ASPECTS
			for (size_t slope_sequence=0; slope_sequence<slope.nSlopes; slope_sequence++) {
				double tot_mass_in = 0.; //Check mass balance over one CALCULATION_STEP_LENGTH if MASS_BALANCE is set
				SnowpackConfig tmpcfg(cfg);
				copyMeteoData(vecMyMeteo[i_stn], sn_MdataT, slope.prevailing_wind_dir, wind_scaling_factor);
				copySnowTemperatures(vecMyMeteo[i_stn], sn_MdataT, fixed_sensor_depths, slope_sequence);
				copySolutes(vecMyMeteo[i_stn], sn_MdataT, SnowStation::number_of_solutes);
				slope.setSlope(slope_sequence, vecXdata, sn_MdataT.dw_drift);
				dataForCurrentTimeStep(sn_MdataT, surfFluxes, vecXdata, slope, tmpcfg,
                               Psolar, Rdata, cumu_hnw, lw_in, iswr_forced, tot_mass_in);

				// Notify user every fifteen days of date being processed
				const double notify_start = floor(vecSSdata[slope.station].profileDate.getJulianDate()) + 15.5;
				if ((mode == "RESEARCH") && (slope.sector == slope.station)
				        && booleanTime(current_date.getJulianDate(false), 15., notify_start, calculation_step_length)) {
					prn_msg(__FILE__, __LINE__, "msg", mio::Date(),
					            "Station %s (%d slope(s)): advanced to %s (%f) station time",
					                vecSSdata[slope.station].meta.stationID.c_str(), slope.nSlopes,
					                    current_date.toString(mio::Date::DIN).c_str(), current_date.getJulianDate());
				}

				// SNOWPACK model (Temperature and Settlement computations)
				Snowpack snowpack(tmpcfg); //the snowpack model to use
				Stability stability(tmpcfg);
				snowpack.runSnowpackModel(sn_MdataT, vecXdata[slope.sector], cumu_hnw, sn_Bdata, surfFluxes);
				stability.checkStability(sn_MdataT, vecXdata[slope.sector]);

				// Update blowing snow flux
				if (slope.virtual_slopes) {
					if (slope.sector == slope.luv) { // Make sure the transported snow comes from the windward slope
						cumulate(cumsum_drift, surfFluxes.drift);
						if (mn_ctrl.HzDump) {
							int i_hz = mn_ctrl.HzStep-1;
							hazard.getDriftIndex(qr_Hdata[i_hz], qr_Hdata_ind[i_hz], &(sn_Zdata.drift24[0]), cumsum_drift,
							                     cos(DEG_TO_RAD(vecXdata[slope.sector].meta.getSlopeAngle())));
							cumsum_drift = 0.;
						}
					}
				} else if (slope.sector == slope.station) { // Consider erosion at the station only
					cumulate(cumsum_drift, surfFluxes.drift);
				}

				/***** OUTPUT SECTION *****/
				surfFluxes.CollectSurfaceFluxes(surfFluxes, sn_Bdata, vecXdata[slope.sector], sn_MdataT,
				                                useSoilLayers, soil_flux);
				if (slope.sector == slope.station) { // station field only (usualy flat)
					int i_hz = mn_ctrl.HzStep;
					// Calculate consistent lw_in for virtual slopes
					if ( vecXdata[slope.station].getNumberOfElements() > 0 ) {
						double k_eff, gradT;
						k_eff =
						    vecXdata[slope.station].Edata[vecXdata[slope.station].getNumberOfElements()-1].k[TEMPERATURE];
						gradT =
						    vecXdata[slope.station].Edata[vecXdata[slope.station].getNumberOfElements()-1].gradT;
						lw_in = k_eff*gradT + sn_Bdata.lw_out - sn_Bdata.qs - sn_Bdata.ql - sn_Bdata.qr;
					} else {
						lw_in = Constants::undefined;
					}
					// Deal with new snow densities
					if (vecXdata[slope.station].hn > 0.) {
						surfFluxes.cRho_hn = vecXdata[slope.station].rho_hn;
						surfFluxes.mRho_hn = sn_MdataT.rho_hn;
					}
					if (mode == "OPERATIONAL") {
						if (!cumsum_mass) {
							// Cumulate flat field runoff in operational mode
							qr_Hdata[i_hz].runoff += surfFluxes.mass[SurfaceFluxes::MS_RUNOFF];
							cumsum_runoff += surfFluxes.mass[SurfaceFluxes::MS_RUNOFF];
						}
						/*
						 * Snow depth and mass corrections (deflate-inflate):
						 *   Monitor snow depth discrepancy assumed to be due to ...
						 *   ... wrong settling, which in turn is assumed to be due to a wrong estimation ...
						 *   of fresh snow mass because Michi spent many painful days calibrating the settling ...
						 *   and therefore it can't be wrong, dixunt Michi and Charles.
						*/
						// Either missed erosion or settling not strong enough
						if ((time_count_deltaHS >= 0.) && (sn_MdataT.hs1 != Constants::nodata)) {
							if ((sn_MdataT.hs1 + 0.01) < (vecXdata[slope.station].cH
							                                 - vecXdata[slope.station].Ground)) {
								time_count_deltaHS += S_TO_D(sn_dt);
							} else {
								time_count_deltaHS = 0.;
							}
						}
						// Settling too strong
						if (time_count_deltaHS <= 0.) {
							if ((sn_MdataT.hs1 - 0.01) > (vecXdata[slope.station].cH - vecXdata[slope.station].Ground)) {
								time_count_deltaHS -= S_TO_D(sn_dt);
							} else {
								time_count_deltaHS = 0.;
							}
						}
						// There is a persistent error for at least one day => apply correction
						if (fabs(1. - time_count_deltaHS) < 0.5 * M_TO_D(calculation_step_length)) {
							deflateInflate(sn_MdataT, vecXdata[slope.station],
							               qr_Hdata[i_hz].dhs_corr, qr_Hdata[i_hz].mass_corr);
							time_count_deltaHS = 0.;
						}
					}

					// SAVE HAZARD DATA ...
					if (mn_ctrl.HzDump) {
						if (cumsum_drift == Constants::undefined) {
							surfFluxes.drift = Constants::undefined;
						} else if ((slope.nSlopes == 1) || !slope.snow_redistribution){
							// Drift Index from flat field
							surfFluxes.drift = cumsum_drift / hazard_steps_between;
							cumsum_drift = 0.;
						} else if (mn_ctrl.nStep == hazard_steps_between) {
							// Zeroth drift index to be written into *.met file
							surfFluxes.drift = cumsum_drift / MAX(1, (hazard_steps_between - 1));
						}

						strncpy(qr_Hdata[i_hz].stat_abbrev, vecStationIDs[i_stn].c_str(), 15);
						if (mode == "OPERATIONAL") {
							qr_Hdata[i_hz].loc_for_snow = (int)vecStationIDs[i_stn][vecStationIDs[i_stn].length()-1];
							//TODO: WHAT SHOULD WE SET HERE? wstat_abk (not existing yet in DB) and wstao_nr, of course;-)
							qr_Hdata_ind[i_hz].loc_for_wind = -1;
						} else {
							qr_Hdata[i_hz].loc_for_snow = 2;
							qr_Hdata[i_hz].loc_for_wind = 1;
						}

						double delta_hs6 = 0.0, delta_hs24 = 0.0;
						getDhs6Dhs24(delta_hs6, delta_hs24, io, vecMyMeteo, i_stn, sn_MdataT.date);

						hazard.getHazardData(qr_Hdata[i_hz], qr_Hdata_ind[i_hz], delta_hs6, delta_hs24,
						                     sn_MdataT, surfFluxes, sn_Zdata, vecXdata[slope.station],
						                     vecXdata[slope.south], slope.nSlopes, slope.virtual_slopes);
						mn_ctrl.HzStep++;
						surfFluxes.hoar = 0.;
					}

					// New snow water equivalent (kg m-2), rain was dealt with in Watertransport.cc
					surfFluxes.mass[SurfaceFluxes::MS_HNW] += vecXdata[slope.station].hn
					                                              * vecXdata[slope.station].rho_hn;
					if (!avgsum_time_series) { // Sum up precipitations
						cumsum_rain += surfFluxes.mass[SurfaceFluxes::MS_RAIN];
						cumsum_snow += surfFluxes.mass[SurfaceFluxes::MS_HNW];
					}
					// In case of erosion ...
					if (vecXdata[slope.station].ErosionMass > 0.) { // Real erosion
						cumsum_erosion[slope.station] = fabs(cumsum_erosion[slope.station])
						                                    + vecXdata[slope.station].ErosionMass;
					} else if (cumsum_erosion[slope.station] > 0.) { // Real and virtual erosion
						cumsum_erosion[slope.station] += surfFluxes.mass[SurfaceFluxes::MS_WIND];
					} else { // Virtual erosion only
						cumsum_erosion[slope.station] -= surfFluxes.mass[SurfaceFluxes::MS_WIND];
					}
				} else { // On slopes take care of cumulative eroded mass on virtual slopes
					cumsum_erosion[slope.sector] += vecXdata[slope.sector].ErosionMass;
				}

				// TIME SERIES (*.met)
				if (tswrite && mn_ctrl.TsDump) {
					// Average fluxes
					if (avgsum_time_series) {
						averageFluxTimeSeries(mn_ctrl.nAvg, useCanopyModel, surfFluxes, vecXdata[slope.sector]);
					} else {
						surfFluxes.mass[SurfaceFluxes::MS_RAIN] = cumsum_rain;
						surfFluxes.mass[SurfaceFluxes::MS_HNW] = cumsum_snow;
						if (slope.sector == slope.lee) {
							surfFluxes.mass[SurfaceFluxes::MS_HNW] += cumsum_erosion[slope.luv]
							        / cos(DEG_TO_RAD(vecXdata[slope.luv].meta.getSlopeAngle()));
						}
					}

					if (precip_rates) { // Precip rates in kg m-2 h-1
						surfFluxes.mass[SurfaceFluxes::MS_RAIN] /= mn_ctrl.nAvg*M_TO_H(calculation_step_length);
						surfFluxes.mass[SurfaceFluxes::MS_HNW] /= mn_ctrl.nAvg*M_TO_H(calculation_step_length);
						if ((mode == "OPERATIONAL") && (!cumsum_mass)) {
							surfFluxes.mass[SurfaceFluxes::MS_RUNOFF] = cumsum_runoff;
							surfFluxes.mass[SurfaceFluxes::MS_RUNOFF] /= mn_ctrl.nAvg*M_TO_H(calculation_step_length);
							cumsum_runoff = 0.;
						}
					}

					// Erosion mass rate in kg m-2 h-1
					surfFluxes.mass[SurfaceFluxes::MS_WIND] = cumsum_erosion[slope.sector];
					surfFluxes.mass[SurfaceFluxes::MS_WIND] /= mn_ctrl.nAvg*M_TO_H(calculation_step_length);

					// Dump
					const unsigned int i_hz = MAX(0, mn_ctrl.HzStep-1);
					double wind_trans24 = qr_Hdata[i_hz].wind_trans24;
					if (i_hz) wind_trans24 = qr_Hdata[i_hz-1].wind_trans24;
					ss.str("");
					ss << vecStationIDs[i_stn];
					if (slope.sector != slope.station) {
						ss << slope.sector;
					}
					snowpackio.writeTimeSeries(vecXdata[slope.sector], surfFluxes, sn_MdataT,
						                       qr_Hdata[i_hz], wind_trans24);

					if (avgsum_time_series) {
						surfFluxes.reset(cumsum_mass);
						if (useCanopyModel) vecXdata[slope.sector].Cdata.reset(cumsum_mass);
					}
					surfFluxes.mRho_hn = surfFluxes.cRho_hn = Constants::undefined;
					// reset cumulative variables
					if (slope_sequence == slope.nSlopes-1) {
						cumsum_erosion.assign(cumsum_erosion.size(), 0.);
						cumsum_rain = cumsum_snow = 0.;
						mn_ctrl.nAvg = 0;
					}
				}

				// SNOW PROFILES ...
				// ... for visualization (*.pro)
				if (profwrite && mn_ctrl.PrDump)
					snowpackio.writeProfile(current_date, vecXdata[slope.sector], qr_Hdata[0]);

				// ... backup Xdata (*.sno<JulianDate>)
				if (mn_ctrl.XdataDump) {
					ss.str("");
					ss << vecStationIDs[i_stn];
					if (slope.sector != slope.station) ss << slope.sector;
					snowpackio.writeSnowCover(current_date, vecXdata[slope.sector],
					                          vecSSdata[slope.sector], sn_Zdata, true);
					prn_msg(__FILE__, __LINE__, "msg", current_date,
					        "Backup Xdata dumped for station %s [%.2f days, step %d]", ss.str().c_str(),
					        (current_date.getJulianDate()
					            - (vecSSdata[slope.station].profileDate.getJulianDate() + 0.5/24)),
					        mn_ctrl.nStep);
				}

				// Mass balance check at end of time step
				if (mass_balance) {
					if (massBalanceCheck(vecXdata[slope.sector], surfFluxes, tot_mass_in) == false)
						prn_msg(__FILE__, __LINE__, "msg+", current_date, "Mass error at end of time step!");
				}
			} //end loop on sectors
			computed_one_timestep = true;
		} while ((dateEnd - current_date).getJulianDate(true) > calculation_step_length/(2.*1440));
		//end loop on timesteps

		if (computed_one_timestep) {
			// Dump the PROFILEs (Xdata) for all available sectors at the end of the Meteo data
			for (size_t sector=slope.station; sector<slope.nSlopes; sector++) {
				if ((mode == "OPERATIONAL") && (sector == slope.station)) {
					// Operational mode ONLY: dump snow depth discrepancy time counter
					vecSSdata[slope.station].TimeCountDeltaHS = time_count_deltaHS;
				}
				snowpackio.writeSnowCover(current_date, vecXdata[sector], vecSSdata[sector], sn_Zdata);
				if (sector == slope.station) {
					prn_msg(__FILE__, __LINE__, "msg", current_date,
					        "Writing data to sno file(s) for %s (station %s)",
					        vecSSdata[slope.station].meta.getStationName().c_str(),
					        vecStationIDs[i_stn].c_str());
				}
			}
			// Dump results to SDB
			if (mode == "OPERATIONAL") {
				mio::Timer sdbDump_timer;
				sdbDump_timer.reset();
				sdbDump_timer.start();
				if (snowpackio.writeHazardData(vecStationIDs[i_stn], qr_Hdata, qr_Hdata_ind, mn_ctrl.HzStep)) {
					sdbDump_timer.stop();
					prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Finished writing Hdata to SDB for station %s (%lf s)",
					        vecStationIDs[i_stn].c_str(), sdbDump_timer.getElapsed());
				}
			}
		}
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Prebuffering meteo data       : %lf s",
		        prebuffering_timer.getElapsed());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Total time to read meteo data : %lf s",
		        meteoRead_timer.getElapsed() + prebuffering_timer.getElapsed());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Runtime for station %s: %lf s",
		        vecStationIDs[i_stn].c_str(), run_timer.getElapsed());
	}

	time_t nowEND=time(NULL);
	cout << endl;
	cout << "[i] []                 STARTED  running SLF " << mode << " Snowpack Model on " << ctime(&nowSRT);
	cout << "                       ===========================================================================" << endl;
	cout << "                       FINISHED running SLF " << mode << " Snowpack Model on " << ctime(&nowEND) << endl;

	return 0;
}
