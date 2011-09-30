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

#include <snowpack/SnowpackConfig.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/
map<string,string> SnowpackConfig::advancedConfig;
map<string,string> SnowpackConfig::inputConfig;
map<string,string> SnowpackConfig::outputConfig;

const bool SnowpackConfig::__init = SnowpackConfig::initStaticData();

bool SnowpackConfig::initStaticData()
{
	//[SnowpackAdvanced] section
	advancedConfig["ALPINE3D"] = "false";
	advancedConfig["DOORSCHOT"] = "false";
	advancedConfig["FIXED_ALBEDO"] = "-999.";
	advancedConfig["FORCE_RH_WATER"] = "true";
	advancedConfig["HARDNESS_MODEL"] = "DEFAULT";
	advancedConfig["HEIGHT_NEW_ELEM"] = "0.02";
	advancedConfig["HN_DENSITY"] = "PARAMETERIZED";
	advancedConfig["HN_DENSITY_MODEL"] = "LEHNING_NEW";
	advancedConfig["HOAR_DENSITY_BURIED"] = "125.";
	advancedConfig["HOAR_DENSITY_SURF"] = "100.";
	advancedConfig["HOAR_MIN_SIZE_BURIED"] = "2.";
	advancedConfig["HOAR_MIN_SIZE_SURF"] = "0.5";
	advancedConfig["HOAR_THRESH_RH"] = "0.97";
	advancedConfig["HOAR_THRESH_VW"] = "3.5";
	advancedConfig["JAM"] = "false";
	advancedConfig["JOIN_ELEMENTS"] = "true";
	advancedConfig["MASS_BALANCE"] = "false";
	advancedConfig["MAX_NUMBER_SENSORS"] = "5";
	advancedConfig["METAMORPHISM_MODEL"] = "DEFAULT";
	advancedConfig["MIN_DEPTH_SUBSURF"] = "0.07";
	advancedConfig["MULTISTREAM"] = "true";
	advancedConfig["NEW_SNOW_GRAIN_RAD"] = "0.15";
	advancedConfig["NUMBER_FIXED_HEIGHTS"] = "5";
	advancedConfig["NUMBER_FIXED_RATES"] = "0";
	advancedConfig["PERP_TO_SLOPE"] = "false";
	advancedConfig["PLASTIC"] = "false";
	advancedConfig["PREVAILING_WIND_DIR"] = "0.";
	advancedConfig["RESEARCH"] = "true";
	advancedConfig["STATION_NAME"] = "station";
	advancedConfig["STRENGTH_MODEL"] = "DEFAULT";
	advancedConfig["SURFACECODE"] = "NEUMANN_BC";
	advancedConfig["SW_MODE_CHANGE"] = "false";
	advancedConfig["THRESH_RAIN"] = "1.2";
	advancedConfig["THRESH_RH"] = "0.5";
	advancedConfig["T_CRAZY_MAX"] = "340.";
	advancedConfig["T_CRAZY_MIN"] = "210.";
	advancedConfig["VARIANT"] = "DEFAULT";
	advancedConfig["VISCOSITY_MODEL"] = "DEFAULT";
	advancedConfig["WET_LAYER"] = "true";
	advancedConfig["WIND_SCALING_FACTOR"] = "1.0";

	//[Input] section
	inputConfig["FIXED_SENSOR_DEPTHS"] = "0.25 0.50 1.0 1.5 -0.1";
	inputConfig["METEOPATH"] = "./DATA/input";
	inputConfig["NUMBER_OF_SOLUTES"] = "0";
	inputConfig["NUMBER_MEAS_TEMPERATURES"] = "0";
	inputConfig["RHO_HN"] = "false";
	inputConfig["SOLUTE_NAMES"] = "NITRATE";
	inputConfig["USEANETZ"] = "false"; // Operational (ImisIO) only
	inputConfig["VW_DRIFT"] = "false";

	//[Output] section
	outputConfig["AVGSUM_TIME_SERIES"] = "true";
	outputConfig["BACKUP_DAYS_BETWEEN"] = "365.";
	outputConfig["CUMSUM_MASS"] = "false";
	outputConfig["EXPERIMENT"] = "NO_EXP";
	outputConfig["FIRST_BACKUP"] = "400.";
	outputConfig["METEOPATH"] = "./DATA";
	outputConfig["OUT_CANOPY"] = "false";
	outputConfig["OUT_HAZ"] = "true";
	outputConfig["OUT_HEAT"] = "true";
	outputConfig["OUT_LOAD"] = "false";
	outputConfig["OUT_LW"] = "true";
	outputConfig["OUT_MASS"] = "true";
	outputConfig["OUT_METEO"] = "true";
	outputConfig["OUT_STAB"] = "true";
	outputConfig["OUT_SW"] = "true";
	outputConfig["OUT_T"] = "true";
	outputConfig["PRECIP_RATES"] = "true";

	return true;
}

/************************************************************
 * non-static section                                       *
 ************************************************************/

SnowpackConfig::~SnowpackConfig() {}

SnowpackConfig::SnowpackConfig(const std::string& i_filename) : Config(i_filename)
{
	string variant; getValue("VARIANT", "SnowpackAdvanced", variant, Config::nothrow);

	int enforce_measured_snow_heights = get("ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack");

	addKey("MINIMUM_L_ELEMENT", "SnowpackAdvanced", "0.0025"); //Minimum element length (m)
	double minimum_l_element = get("MINIMUM_L_ELEMENT", "SnowpackAdvanced");
	if (enforce_measured_snow_heights) {
		addKey("HEIGHT_NEW_ELEM", "SnowpackAdvanced", "0.02");
	} else {
		stringstream ss;
		double tmp = 2.0 * minimum_l_element;
		ss << tmp;
		addKey("HEIGHT_NEW_ELEM", "SnowpackAdvanced", ss.str());
	}

	string hn_density;  getValue("HN_DENSITY", "SnowpackAdvanced", hn_density, Config::nothrow);
	if (hn_density == "MEASURED") {
		bool rho_hn = false;
		getValue("RHO_HN", "Input", rho_hn, Config::nothrow);
		if (!rho_hn)
			throw InvalidArgumentException("HN_DENSITY = " + hn_density + " while RHO_HN = false", AT);
	}
	string hn_density_model; getValue("HN_DENSITY_MODEL", "SnowpackAdvanced", hn_density_model, Config::nothrow);
	string metamorphism_model; getValue("METAMORPHISM_MODEL", "SnowpackAdvanced", metamorphism_model, Config::nothrow);
	string strength_model; getValue("STRENGTH_MODEL", "SnowpackAdvanced", strength_model, Config::nothrow);
	string viscosity_model; getValue("VISCOSITY_MODEL", "SnowpackAdvanced", viscosity_model, Config::nothrow);

	if ((variant == "") || (variant == "DEFAULT")) {

		// Use default settings

	} else if (variant == "JAPAN") {

		if (metamorphism_model == "") addKey("METAMORPHISM_MODEL", "SnowpackAdvanced", "NIED");
		if (strength_model == "") addKey("STRENGTH_MODEL", "SnowpackAdvanced", "NIED");
		if (viscosity_model == "") addKey("VISCOSITY_MODEL", "SnowpackAdvanced", "KOJIMA");

	} else if (variant == "ANTARCTICA") {

		if (hn_density == "") addKey("HN_DENSITY", "SnowpackAdvanced", "EVENT");

		addKey("MINIMUM_L_ELEMENT", "SnowpackAdvanced", "0.0001"); //Minimum element length (m)
		minimum_l_element = get("MINIMUM_L_ELEMENT", "SnowpackAdvanced");

		string hoar_density_buried; getValue("HOAR_DENSITY_BURIED", "SnowpackAdvanced", hoar_density_buried, Config::nothrow);
		if (hoar_density_buried == "") addKey("HOAR_DENSITY_BURIED", "SnowpackAdvanced", "200.0");

		string force_rh_water; getValue("FORCE_RH_WATER", "SnowpackAdvanced", force_rh_water, Config::nothrow);
		if (force_rh_water == "") addKey("FORCE_RH_WATER", "SnowpackAdvanced", "false");

		string thresh_rh; getValue("THRESH_RH", "SnowpackAdvanced", thresh_rh, Config::nothrow);
		if (thresh_rh == "") addKey("THRESH_RH", "SnowpackAdvanced", "0.7");

		if ( !enforce_measured_snow_heights) {
			stringstream ss;
			double tmp = 1.1 * minimum_l_element;
			ss << tmp;
			addKey("HEIGHT_NEW_ELEM", "SnowpackAdvanced", ss.str());
		}
		addKey("FIRST_BACKUP", "Output", "1500.");
		addKey("NUMBER_FIXED_HEIGHTS", "SnowpackAdvanced", "7");
		addKey("FIXED_RATES", "SnowpackAdvanced", "false");
		addKey("NUMBER_FIXED_RATES", "SnowpackAdvanced", "0");
		addKey("MAX_NUMBER_SENSORS", "SnowpackAdvanced", "7");
		addKey("MIN_DEPTH_SUBSURF", "SnowpackAdvanced", "0.");
		addKey("T_CRAZY_MIN", "SnowpackAdvanced", "165.");
		addKey("T_CRAZY_MAX", "SnowpackAdvanced", "300.");
		addKey("NEW_SNOW_GRAIN_RAD", "SnowpackAdvanced", "0.1");

	} else if (variant == "CALIBRATION") {

		if (hn_density_model == "") addKey("HN_DENSITY_MODEL", "SnowpackAdvanced", "ZWART");
		if (viscosity_model == "") addKey("VISCOSITY_MODEL", "SnowpackAdvanced", "CALIBRATION");

		string number_fixed_heights; getValue("NUMBER_FIXED_HEIGHTS", "SnowpackAdvanced", number_fixed_heights, Config::nothrow);
		if (number_fixed_heights == "") addKey("NUMBER_FIXED_HEIGHTS", "SnowpackAdvanced", "5");
		string number_fixed_rates; getValue("NUMBER_FIXED_RATES", "SnowpackAdvanced", number_fixed_rates, Config::nothrow);
		if (number_fixed_rates == "") addKey("NUMBER_FIXED_RATES", "SnowpackAdvanced", "0");
		string max_number_sensors; getValue("MAX_NUMBER_SENSORS", "SnowpackAdvanced", max_number_sensors, Config::nothrow);
		if (max_number_sensors == "") addKey("MAX_NUMBER_SENSORS", "SnowpackAdvanced", "5");
		string min_depth_subsurf; getValue("MIN_DEPTH_SUBSURF", "SnowpackAdvanced", min_depth_subsurf, Config::nothrow);
		if (min_depth_subsurf == "") addKey("MIN_DEPTH_SUBSURF", "SnowpackAdvanced", "0.0");

	} else {
		throw UnknownValueException("Unknown variant " + variant, AT);
	}

	/* For all parameters not set by the user or by the initialization above, the default values apply
	 * That is, loop through advancedConfig (then inputConfig & outputConfig) and check whether user has set
	 * the parameter in the corresponding section, if not add default value
	 */
	for(map<string,string>::const_iterator it = advancedConfig.begin(); it != advancedConfig.end(); it++){
		//[SnowpackAdvanced] section
		string value; getValue(it->first, "SnowpackAdvanced", value, Config::nothrow);
		if (value == "") addKey(it->first, "SnowpackAdvanced", it->second);
	}

	for(map<string,string>::const_iterator it = inputConfig.begin(); it != inputConfig.end(); it++){
		//[Input] section
		string value; getValue(it->first, "Input", value, Config::nothrow);
		if (value == "") addKey(it->first, "Input", it->second);
	}

	for(map<string,string>::const_iterator it = outputConfig.begin(); it != outputConfig.end(); it++){
		//[Output] section
		string value; getValue(it->first, "Output", value, Config::nothrow);
		if (value == "") addKey(it->first, "Output", it->second);
	}

	/**
	 * @brief Defines how energy and mass balance are output \n
	 * - AVGSUM_TIME_SERIES == 1 \n
	 *   Energy and mass fluxes are averaged and cumulated over TS_DAYS_BETWEEN, respectively. \n
	 *   Otherwise, instantaneous energy fluxes and mass fluxes cumulated over the last computation
	 *   time step (CALCULATION_STEP_LENGTH) are dumped. \n
	 *   @note Precipitations and Erosion are always given in rates per TS_DAYS_BETWEEN interval (h-1)
	 * - CUMSUM_MASS == 1 \n
	 *   Mass fluxes are cumulated over whole run period.
	 * - WARNING: In operational mode and if NUMBER_SLOPES > 1, the above two values are always unset!
	 */
	unsigned int nSlopes = get("NUMBER_SLOPES", "Snowpack");
	if (nSlopes > 1) {
		addKey("AVGSUM_TIME_SERIES", "Output", "false");
		addKey("CUMSUM_MASS", "Output", "false");
	}

	/**
	 * @brief Hazard data interval in units of CALCULATION_STEP_LENGTH \n
	 * WARNING: In operational mode, this has to result in a 30 min interval!
	 * It is a matter of consitency. If you change this, a big mess will result!!!
	 */
	double calculation_step_length = get("CALCULATION_STEP_LENGTH", "Snowpack");

	string hazard_steps_between; getValue("HAZARD_STEPS_BETWEEN", "Output", hazard_steps_between, Config::nothrow);
	if (hazard_steps_between == "") {
		stringstream ss;
		int tmp = (int)(30./calculation_step_length + 0.5);
		ss << tmp;
		addKey("HAZARD_STEPS_BETWEEN", "Output", ss.str());
	}

	/**
	 * @brief Defines how depths of snow temperature sensors are read in and output \n
	 * - If measured snow temperatures are available, default or user provided sensor depths
	 *     from the input section will be used for output
	 * - If no measured snow temperatures are available, default or user provided sensor depths
	 *     from the output section will be used for output
	 * @note default depths are provided for the input section only
	 */
	unsigned int number_meas_temperatures = get("NUMBER_MEAS_TEMPERATURES", "Input", Config::nothrow);
	if (number_meas_temperatures > 0) {
		string i_fixed_sensor_depths; getValue("FIXED_SENSOR_DEPTHS", "Input", i_fixed_sensor_depths, Config::nothrow);
		addKey("FIXED_SENSOR_DEPTHS", "Output", i_fixed_sensor_depths);
	} else {
		string o_fixed_sensor_depths; getValue("FIXED_SENSOR_DEPTHS", "Output", o_fixed_sensor_depths, Config::nothrow);
		if (o_fixed_sensor_depths == "")
			addKey("FIXED_SENSOR_DEPTHS", "Output", "0.25 0.50 1.0 1.5 -0.1");
		else
			addKey("FIXED_SENSOR_DEPTHS", "Output", o_fixed_sensor_depths);
		stringstream ss;
		ss << "-999.";
		addKey("FIXED_SENSOR_DEPTHS", "Input", ss.str());
	}
	unsigned int number_fixed_heights = get("NUMBER_FIXED_HEIGHTS", "SnowpackAdvanced", Config::nothrow);
	vector<double> fixed_sensor_depths = get("FIXED_SENSOR_DEPTHS", "Output");
	if (fixed_sensor_depths.size() > number_fixed_heights) {
		stringstream ss;
		ss << number_fixed_heights;
		throw InvalidArgumentException("At most NUMBER_FIXED_HEIGHTS="+ss.str()+" sensor depths allowed", AT);
	}
}

void checkUserConfiguration(mio::Config& /*cfg*/)
{
	/*
	if ( !((NUMBER_MEAS_TEMP > -1) && (NUMBER_MEAS_TEMPERATURES <= MAX_NUMBER_SENSORS)) ) {
	prn_msg(__FILE__, __LINE__, "err", Date(), "NUMBER_MEAS_TEMPERATURES=%d out of range (0, %d)", NR_MEAS_TEMP, MAX_NUMBER_SENSORS);
		return ERROR;
	}
	NUMBER_SENSORS = NUMBER_FIXED_HEIGHTS + NUMBER_FIXED_RATES;
	if ( !((NUMBER_SENSORS > -1) && (NUMBER_SENSORS <= MAX_NUMBER_SENSORS)) ) {
	prn_msg(__FILE__, __LINE__, "err", Date(), "%d NUMBER_FIXED_HEIGHTS + %d NUMBER_FIXED_RATES out of range (0, %d)", NUMBER_FIXED_HEIGHTS, NUMBER_FIXED_RATES, MAX_NUMBER_SENSORS);
		return ERROR;
	}


	// Print some infos to stdout: research mode
	if ( !MEAS_TSS && CHANGE_BC ){
	prn_msg(__FILE__, __LINE__, "wrn", Date(), "(!MEAS_TSS && CHANGE_BC) == 1");
	prn_msg(__FILE__, __LINE__, "msg", Date(), "Using Neumann boundary conditions because no measured TSS is available");
		CHANGE_BC = 0;
	}
	if ( !SNP_SOIL && CANOPY ) {
	prn_msg(__FILE__, __LINE__, "wrn", Date(), "Canopy Model is used WITHOUT soil data!");
	}


	// Check a few settings - operational mode
	if ( SNOW_REDISTRIBUTION && ((NUMBER_SLOPES > 5) || ((NUMBER_SLOPES-1)%4 != 0)) ) {
	prn_msg(__FILE__, __LINE__, "err", Date(), "NUMBER_SLOPES (%d) not compatible with SNOW_REDISTRIBUTION", NUMBER_SLOPES);
		exit(EXIT_FAILURE);
	}
	if ( CANOPY ) {
		CANOPY = 0;
	prn_msg(__FILE__, __LINE__, "wrn", Date(), "CANOPY was set! You may have run into troubles! Reset to 0");
	}
	*/
}

