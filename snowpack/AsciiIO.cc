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

#include <snowpack/AsciiIO.h>

using namespace std;
using namespace mio;

/************************************************************
 * static section                                           *
 ************************************************************/

/// @brief Defines whether hardness (R) is output either in N (Swiss scale) or steps
const bool AsciiIO::r_in_n = true;

/// @brief Defines whether surface temperature is included in comparison; If no
///        values are provided in Master-file, they will be extrapolated
const bool AsciiIO::t_srf = false;

/// @brief Defines whether snow/ground temperature is included in comparison; If no
///        values are provided in Master-file, they will be extrapolated
const bool AsciiIO::t_gnd = false;

/************************************************************
 * non-static section                                       *
 ************************************************************/

AsciiIO::AsciiIO(const mio::Config& cfg)
{
	/**
	 * @name Snow sensors
	 * @brief Defines the number of modelled and/or measured sensors that are monitored
	 */
	/**
	 * @brief description:
	 * - NUMBER_FIXED_HEIGHTS (default: 5) modelled and measured temperatures at fixed positions (m) can be monitored.
	 * 	- positive position values: heigth from ground surface (snow only)
	 * 	- negative position values: depth from either ground surface or snow surface if SNP_SOIL = 0
	 * - NOTE:
	 *  - A sensor must at least be covered by MIN_DEPTH_SUBSURF (m) snow for its temperature to be output
	 * 	- At most MAX_NUMBER_SENSORS can be monitored (maximum 10+44 columns are available).
	 *    If (CANOPY && OUT_CANOPY), however, only 10 values can be dumped to file, that is,
	 *    five modelled and five measured temperatures.
	 * 	- NUMBER_MEAS_TEMPERATURES (<= MAX_NUMBER_SENSORS) is the number of monitored PHYSICAL sensors,
	 *    typically thermometers (and their position). A mix of sensors at fixed positions (NUMBER_FIXED_HEIGHTS) and
	 *    variable positions is allowed.
	 * 		-# research mode: up to 5
	 * 		-# advanced mode: up to 5+22 = 27
	 * 			- Antarctica 7+7+9 = 23
	 * 			- Calibration 5+11 = 16
	 * 		-# operational mode: 3
	 * 	- Add data columns to the input file after the snow depth as follows:
	 * 		-# NUMBER_FIXED_HEIGHTS (\< NUMBER_MEAS_TEMPERATURES) columns of measured values at fixed positions
	 * 		-# NUMBER_FIXED_RATES (\< NUMBER_MEAS_TEMPERATURES-NUMBER_FIXED_HEIGHTS) double columns of measured values and positions given by
	 *       fixed negative settling rates (d-1). Enter the initial position (pos or neg, see above) at the time
	 *       the change is effective, followed by the rate. Additional such entries can be repeated when initial
	 *       position and/or rate change again. Otherwise fill the column with -999.9 [IOUtils::nodata]
	 * 		-# Up to NUMBER_MEAS_TEMPERATURES-NUMBER_FIXED_HEIGHTS-NUMBER_FIXED_RATES double columns of measured values and positions
	 */
	// Snowpack section
	cfg.getValue("CALCULATION_STEP_LENGTH", "Snowpack", calculation_step_length);
	cfg.getValue("CANOPY", "Snowpack", useCanopyModel);
	cfg.getValue("SNP_SOIL", "Snowpack", useSoilLayers);
	cfg.getValue("SW_MODE", "Snowpack", sw_mode);
	sw_mode %= 10;

	// SnowpackAdvanced section
	cfg.getValue("HOAR_DENSITY_SURF", "SnowpackAdvanced", hoar_density_surf); // Density of SH at surface node (kg m-3)
	cfg.getValue("HOAR_MIN_SIZE_SURF", "SnowpackAdvanced", hoar_min_size_surf); // Minimum size to show SH on surface (mm)
	cfg.getValue("MAX_NUMBER_SENSORS", "SnowpackAdvanced", max_number_sensors);
	cfg.getValue("MIN_DEPTH_SUBSURF", "SnowpackAdvanced", min_depth_subsurf);
	cfg.getValue("NUMBER_FIXED_HEIGHTS", "SnowpackAdvanced", number_fixed_heights);
	cfg.getValue("NUMBER_FIXED_RATES", "SnowpackAdvanced", number_fixed_rates);
	number_sensors = number_fixed_heights + number_fixed_rates;

	cfg.getValue("HN_DENSITY", "SnowpackAdvanced", hn_density);
	cfg.getValue("HN_DENSITY_MODEL", "SnowpackAdvanced", hn_density_model);

	cfg.getValue("PERP_TO_SLOPE", "SnowpackAdvanced", perp_to_slope);
	cfg.getValue("RESEARCH", "SnowpackAdvanced", research_mode);
	cfg.getValue("VARIANT", "SnowpackAdvanced", variant);

	// Input section
	cfg.getValue("NUMBER_MEAS_TEMPERATURES", "Input", number_meas_temperatures);
	if (number_meas_temperatures > 0)
		cfg.getValue("FIXED_SENSOR_DEPTHS", "Input", fixed_sensor_depths);
	cfg.getValue("METEOPATH", "Input", inpath, Config::nothrow);
	string snowpath("");
	cfg.getValue("SNOWPATH", "Input", snowpath, Config::nothrow);
	if (snowpath != "") {
		i_snopath = snowpath;
	} else {
		i_snopath = inpath;
	}
	cfg.getValue("STATION_NAME", "Input", station_name, Config::nothrow);
	cfg.getValue("TIME_ZONE", "Input", time_zone);

	// Output section
	if (number_meas_temperatures == 0)
		cfg.getValue("FIXED_SENSOR_DEPTHS", "Output", fixed_sensor_depths);
	cfg.getValue("AVGSUM_TIME_SERIES", "Output", avgsum_time_series, Config::nothrow);
	cfg.getValue("EXPERIMENT", "Output", experiment);
	cfg.getValue("HAZARD_STEPS_BETWEEN", "Output", hazard_steps_between);
	cfg.getValue("METEOPATH", "Output", outpath, Config::nothrow);
	cfg.getValue("OUT_CANOPY", "Output", out_canopy);
	cfg.getValue("OUT_HAZ", "Output", out_haz);
	cfg.getValue("OUT_HEAT", "Output", out_heat);
	cfg.getValue("OUT_LOAD", "Output", out_load);
	cfg.getValue("OUT_LW", "Output", out_lw);
	cfg.getValue("OUT_MASS", "Output", out_mass);
	cfg.getValue("OUT_METEO", "Output", out_meteo);
	cfg.getValue("OUT_STAB", "Output", out_stab);
	cfg.getValue("OUT_SW", "Output", out_sw);
	cfg.getValue("OUT_T", "Output", out_t);
	cfg.getValue("SNOWPATH", "Output", snowpath, Config::nothrow);
	if (snowpath != "") {
		o_snopath = snowpath;
	} else {
		o_snopath = outpath;
	}
	cfg.getValue("TS_DAYS_BETWEEN", "Output", ts_days_between);
}
/*
AsciiIO::~AsciiIO() throw()
{
	cleanup();
}

void AsciiIO::cleanup() throw()
{
	if (fout.is_open()) {//close fin if open
		fout.close();
	}
}
*/

/**
 * @brief This routine reads the status of the snow cover at program start
 * @version 10.02
 * @param i_snowfile file containing the initial state of the snowpack
 * @param stationID
 * @param SSdata
 * @param Zdata
 */
void AsciiIO::readSnowCover(const std::string& i_snowfile, const std::string& stationID,
                            SN_SNOWSOIL_DATA& SSdata, SN_ZWISCHEN_DATA& Zdata)
{
	unsigned int ii;
	string snofilename = getFilenamePrefix(i_snowfile, i_snopath, false);
	if (snofilename.rfind(".snoold") == string::npos) {
		snofilename += ".snoold";
	}

	FILE *fin=NULL;
	fin = fopen(snofilename.c_str(), "r");
	if (fin == NULL) {
		prn_msg(__FILE__, __LINE__, "msg+", Date(), "Cannot open profile INPUT file: %s", snofilename.c_str());
		throw IOException("Cannot generate Xdata from file "+snofilename, AT);
	}

	// Header, Station Name and Julian Date
	char station_name[MAX_STRING_LENGTH];
	if (fscanf(fin, " %*s") != 0) {
		fclose(fin);
		throw InvalidFormatException("Can not read header of file "+snofilename, AT);
	}
	if (fscanf(fin, "\nStationName= %s", station_name) != 1) {
		fclose(fin);
		throw InvalidFormatException("Can not read StationName in file "+snofilename, AT);
	}
	int YYYY, MM, DD, HH, MI, dum;
	if (fscanf(fin, "\nProfileDate= %4d %2d %2d %2d %2d", &YYYY, &MM, &DD, &HH, &MI) != 5) {
		fclose(fin);
		throw InvalidFormatException("Can not read ProfileDate in file "+snofilename, AT);
	}
	SSdata.profileDate = Date(YYYY, MM, DD, HH, MI, time_zone);

	// Last checked measured Snow Height used for data Control of next run
	if (fscanf(fin, "\nHS_Last=%lf", &SSdata.HS_last) != 1) {
		fclose(fin);
		throw InvalidFormatException("Can not read HS_Last in file "+snofilename, AT);
	}
	double latitude, longitude, altitude;
	if (fscanf(fin, "\nLatitude=%lf", &latitude) != 1) {
		fclose(fin);
		throw InvalidFormatException("Can not read Latitude in file "+snofilename, AT);
	}
	if (fscanf(fin, "\nLongitude=%lf", &longitude) != 1) {
		fclose(fin);
		throw InvalidFormatException("Can not read Longitude in file "+snofilename, AT);
	}
	if (fscanf(fin, "\nAltitude=%lf", &altitude) != 1) {
		fclose(fin);
		throw InvalidFormatException("Can not read Altitude in file "+snofilename, AT);
	}
	double slope_angle, azi;
	if (fscanf(fin, "\nSlopeAngle=%lf", &slope_angle) != 1) {
		fclose(fin);
		throw InvalidFormatException("Can not read SlopeAngle in file "+snofilename, AT);
	}
	if (fscanf(fin, "\nSlopeAzi=%lf", &azi) != 1) {
		fclose(fin);
		throw InvalidFormatException("Can not read SlopeAzi in file "+snofilename, AT);
	}

	mio::Coords tmppos;
	tmppos.setLatLon(latitude, longitude, altitude);
	SSdata.meta.setStationData(tmppos, stationID, station_name);
	SSdata.meta.setSlope(slope_angle, azi);

	// Check consistency with radiation switch
	if ((sw_mode == 2) && perp_to_slope && (SSdata.meta.getSlopeAngle() > Constants::min_slope_angle)) {
		prn_msg(__FILE__, __LINE__, "wrn", Date(),
		        "You want to use measured albedo in a slope steeper than 3 deg  with PERP_TO_SLOPE set!");
		throw IOException("Do not generate Xdata from file "+snofilename, AT);
	}

	// Check consistency of nXLayerData
	if (fscanf(fin, "\nnSoilLayerData=%d", &dum) != 1) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Missing 'nSoilLayerData'");
		throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
	}
	if (dum < 0) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "'nSoilLayerData' < 0 !!!");
		throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
	} else if (useSoilLayers && (dum < 1)) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "useSoilLayers set but 'nSoilLayerData' < 1 !!!");
		throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
	} else if (!useSoilLayers && (dum > 0)) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "useSoilLayers not set but 'nSoilLayerData' > 0 !!!");
		throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
	}
	SSdata.nLayers = dum;
	if (fscanf(fin, "\nnSnowLayerData=%d", &dum) != 1) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Missing 'nSnowLayerData'");
		throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
	}
	if (dum < 0) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "'nSnowLayerData' < 0  !!!");
		throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
	}
	SSdata.nLayers += dum;

	if (fscanf(fin, "\nSoilAlbedo=%lf", &SSdata.SoilAlb) != 1) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Missing 'SoilAlbedo'");
		throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
	}
	if (fscanf(fin, "\nBareSoil_z0=%lf", &SSdata.BareSoil_z0) != 1) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Missing 'BareSoil_z0'");
		throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
	}
	if (SSdata.BareSoil_z0==0.) {
		prn_msg(__FILE__, __LINE__, "wrn", Date(), "'BareSoil_z0'=0 from %s, reset to 0.02", snowfile.c_str());
		SSdata.BareSoil_z0=0.02;
	}
	if (SSdata.HS_last > 0.05) {
		SSdata.Albedo = 0.9;
	} else {
		SSdata.Albedo = SSdata.SoilAlb;
	}

	if (fscanf(fin, "\nCanopyHeight=%lf",&SSdata.Canopy_Height) != 1) {
		fclose(fin);
		throw InvalidFormatException("Can not read CanopyHeight in file "+snofilename, AT);
	}
	if (fscanf(fin, "\nCanopyLeafAreaIndex=%lf",&SSdata.Canopy_LAI) != 1) {
		fclose(fin);
		throw InvalidFormatException("Can not read CanopyLeafAreaIndex in file "+snofilename, AT);
	}
	if (fscanf(fin, "\nCanopyDirectThroughfall=%lf",&SSdata.Canopy_Direct_Throughfall) != 1) {
		fclose(fin);
		throw InvalidFormatException("Can not read CanopyDirectThroughfall in file "+snofilename, AT);
	}

	if (fscanf(fin, "\nWindScalingFactor=%lf",&SSdata.WindScalingFactor) != 1) {
		fclose(fin);
		throw InvalidFormatException("Can not read WindScalingFactor in file "+snofilename, AT);
	}
	if (fscanf(fin, "\nErosionLevel=%d",&SSdata.ErosionLevel) != 1) {
		fclose(fin);
		throw InvalidFormatException("Can not read ErosionLevel in file "+snofilename, AT);
	}
	if (fscanf(fin, "\nTimeCountDeltaHS=%lf",&SSdata.TimeCountDeltaHS) != 1) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Failed reading canopy or additional parameters");
		throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
	}

	if (fscanf(fin,"\nYYYY") < 0) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Failed reading layer header starting with 'YYYY'");
		throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
	}
	if( fscanf(fin, "%*[^\n]") != 0) {
		fclose(fin);
		throw InvalidFormatException("Can not read header end in file "+snofilename, AT);
	}

	int nFields = 0;
	if (SSdata.nLayers > 0)
		SSdata.Ldata.resize(SSdata.nLayers, LayerData());
	for (unsigned int ll = 0; ll < SSdata.nLayers; ll++) {
		if ((nFields = fscanf(fin, " %d %d %d %d %d", &YYYY, &MM, &DD, &HH, &MI)) != 5) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "Failed reading date: read %d fields", nFields);
			throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
		}
		SSdata.Ldata[ll].layerDate = Date(YYYY, MM, DD, HH, MI, time_zone);
		if (SSdata.Ldata[ll].layerDate > SSdata.profileDate) {
			prn_msg(__FILE__, __LINE__, "err", Date(),
			        "Layer %d from bottom is younger (%f) than ProfileDate (%f) !!!",
			        ll+1, SSdata.Ldata[ll].layerDate.getJulianDate(), SSdata.profileDate.getJulianDate());
			throw IOException("Cannot generate Xdata from file "+snofilename, AT);
		}
		if ((nFields = fscanf(fin, " %lf %lf %lf %lf %lf %lf",
		                      &SSdata.Ldata[ll].hl, &SSdata.Ldata[ll].tl, &SSdata.Ldata[ll].phiIce,
		                      &SSdata.Ldata[ll].phiWater, &SSdata.Ldata[ll].phiVoids, &SSdata.Ldata[ll].phiSoil)) != 6) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "Failed reading hl etc: read %d of 6 fields", nFields);
			throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
		}
		if (SSdata.Ldata[ll].tl < 100.) {
			SSdata.Ldata[ll].tl = C_TO_K(SSdata.Ldata[ll].tl);
		}
		if ((nFields = fscanf(fin, "%lf %lf %lf", &SSdata.Ldata[ll].SoilRho, &SSdata.Ldata[ll].SoilK,
		                      &SSdata.Ldata[ll].SoilC)) != 3) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "Failed reading SoilRho etc: read %d of 3 fields", nFields);
			throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
		}
		if ((nFields = fscanf(fin, "%lf %lf %lf %lf %u %lf %u", &SSdata.Ldata[ll].rg, &SSdata.Ldata[ll].rb,
		                      &SSdata.Ldata[ll].dd, &SSdata.Ldata[ll].sp, &SSdata.Ldata[ll].mk,
		                      &SSdata.Ldata[ll].hr, &SSdata.Ldata[ll].ne)) != 7) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "Failed reading rg etc: read %d of 7 fields", nFields);
			throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
		}
		if (SSdata.Ldata[ll].phiSoil==0. && (SSdata.Ldata[ll].rg<=0. || SSdata.Ldata[ll].rb<=0.)) { //Test only for snow layers
			std::stringstream ss;
			ss << "Invalid grain specification in layer " << ll+1 << " (from bottom) of file " << snofilename << ": ";
			ss << "grain radius = " << SSdata.Ldata[ll].rg << " bond radius = " << SSdata.Ldata[ll].rb;
			ss << " (they should be > 0).";
			throw InvalidArgumentException(ss.str(), AT);
		}
		if (SSdata.Ldata[ll].rg>0. && SSdata.Ldata[ll].rb >= SSdata.Ldata[ll].rg) {
			//HACK To avoid surprises in lwsn_ConcaveNeckRadius()
			SSdata.Ldata[ll].rb = Metamorphism::max_grain_bond_ratio * SSdata.Ldata[ll].rg;
			prn_msg(__FILE__, __LINE__, "wrn", Date(), "Layer %d from bottom: bond radius rb/rg larger than Metamorphism::max_grain_bond_ratio=%f (rb=%f mm, rg=%f mm)! Reset to Metamorphism::max_grain_bond_ratio",
			        ll+1, Metamorphism::max_grain_bond_ratio, SSdata.Ldata[ll].rb, SSdata.Ldata[ll].rg);
		}
		if ((nFields = fscanf(fin, "%lf %lf", &SSdata.Ldata[ll].CDot, &SSdata.Ldata[ll].metamo)) != 2) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "Failed reading CDot etc: read %d of 2 fields", nFields);
			throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
		}
		for (ii = 0; ii < SnowStation::number_of_solutes; ii++) {
			if ((nFields = fscanf(fin," %lf %lf %lf %lf ",
			                      &SSdata.Ldata[ll].cIce[ii], &SSdata.Ldata[ll].cWater[ii],
			                      &SSdata.Ldata[ll].cVoids[ii], &SSdata.Ldata[ll].cSoil[ii])) != 4) {
				prn_msg(__FILE__, __LINE__, "err", Date(),
				        "Failed reading impurity concentrations: read %d of 4 fields", nFields);
				throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
			}
		}
	}

	// Read the hoar, drift, and snowfall hazard data info (Zdata, needed for flat field only)
	if( fscanf(fin,"%*s ") != 0) {
		fclose(fin);
		throw InvalidFormatException("Can not read spacing in file "+snofilename, AT);
	}
	for (ii = 0; ii < 48; ii++) {
		if (fscanf(fin," %lf ", &Zdata.hoar24[ii]) != 1) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "While reading hoar data (48) !!!");
			throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
		}
	}
	if( fscanf(fin,"%*s ") != 0) {
		fclose(fin);
		throw InvalidFormatException("Can not read spacing in file "+snofilename, AT);
	}
	for (ii = 0; ii < 48; ii++) {
		if (fscanf(fin," %lf ", &Zdata.drift24[ii]) != 1) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "While reading drift data (48)  !!!");
			throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
		}
	}
	if( fscanf(fin,"%*s ") != 0) {
		fclose(fin);
		throw InvalidFormatException("Can not read spacing in file "+snofilename, AT);
	}
	for (ii = 0; ii < 144; ii++) {
		if (fscanf(fin," %lf ", &Zdata.hn3[ii]) != 1) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "While reading hn(3h) data (144) !!!");
			throw InvalidFormatException("While reading Zdata (hns3) !!!", AT);
		}
	}
	if( fscanf(fin,"%*s ") != 0) {
		fclose(fin);
		throw InvalidFormatException("Can not read spacing in file "+snofilename, AT);
	}
	for (ii = 0; ii < 144; ii++) {
		if (fscanf(fin," %lf ", &Zdata.hn24[ii]) != 1) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "While reading hn(24h) data (144)  !!!");
			throw InvalidFormatException("Cannot generate Xdata from file "+snofilename, AT);
		}
	}

	SSdata.nN = 1;
	SSdata.Height = 0.;
	for (unsigned int ll = 0; ll < SSdata.nLayers; ll++) {
		SSdata.nN    += SSdata.Ldata[ll].ne;
		SSdata.Height += SSdata.Ldata[ll].hl;
	}

	fclose(fin);
}

/**
 * @brief This routine writes the status of the snow cover at program termination and at specified backup times
 * @version 11.02
 * @param date current
 * @param Xdata
 * @param SSdata
 * @param Zdata
 * @param forbackup dump Xdata on the go
 */
void AsciiIO::writeSnowCover(const mio::Date& date, const SnowStation& Xdata, const SN_SNOWSOIL_DATA& SSdata,
                             const SN_ZWISCHEN_DATA& Zdata, const bool& forbackup)
{
	unsigned int ii, e;
	FILE *fout=NULL;
	string snofilename = getFilenamePrefix(Xdata.meta.getStationID().c_str(), o_snopath) + ".sno";

	if (forbackup){
		stringstream ss;
		ss << (int)(date.getJulianDate() + 0.5);
		snofilename += ss.str();
	}

	const vector<ElementData>& EMS = Xdata.Edata;
	fout = fopen(snofilename.c_str(), "w");
	if (fout == NULL) {
		prn_msg(__FILE__, __LINE__, "err", date,"Cannot open profile OUTPUT file: %s", snofilename.c_str());
		throw FileAccessException("Cannot dump final Xdata to file "+snofilename, AT);
	}

	// Header, Station Name and Julian Day
	fprintf(fout, "[SNOWPACK_INITIALIZATION]");
	fprintf(fout, "\nStationName= %s", Xdata.meta.getStationName().c_str());

	int yyyy,mm,dd,hh,mi;
	date.getDate(yyyy,mm,dd,hh,mi);

	fprintf(fout, "\nProfileDate= %04d %02d %02d %02d %02d", yyyy, mm, dd, hh, mi);

	// Last checked Snow Depth used for data Control of next run
	fprintf(fout, "\nHS_Last= %f", Xdata.cH - Xdata.Ground);

	// Latitude, Longitude, Altitude, Slope Angle, Slope Azimut
	fprintf(fout, "\nLatitude= %.4f",   Xdata.meta.position.getLat());
	fprintf(fout, "\nLongitude= %.4f",  Xdata.meta.position.getLon());
	fprintf(fout, "\nAltitude= %.0f",   Xdata.meta.position.getAltitude());
	fprintf(fout, "\nSlopeAngle= %.2f", Xdata.meta.getSlopeAngle());
	fprintf(fout, "\nSlopeAzi= %.2f",   Xdata.meta.getAzimuth());

	// Number of Soil Layer Data; in case of no soil used to store the erosion level
	fprintf(fout, "\nnSoilLayerData= %d", Xdata.SoilNode);
	// Number of Snow Layer Data
	fprintf(fout, "\nnSnowLayerData= %d", Xdata.getNumberOfElements() - Xdata.SoilNode);

	// Ground Characteristics (introduced June 2006)
	fprintf(fout, "\nSoilAlbedo= %.2f", Xdata.SoilAlb);
	fprintf(fout, "\nBareSoil_z0= %.3f", Xdata.BareSoil_z0);
	// Canopy Characteristics
	fprintf(fout, "\nCanopyHeight= %.2f",Xdata.Cdata.height);
	fprintf(fout, "\nCanopyLeafAreaIndex= %.6f",Xdata.Cdata.lai);
	fprintf(fout, "\nCanopyDirectThroughfall= %.2f",Xdata.Cdata.direct_throughfall);
	// Additional parameters
	fprintf(fout,"\nWindScalingFactor= %f",SSdata.WindScalingFactor);
	fprintf(fout,"\nErosionLevel= %d",Xdata.ErosionLevel);
	fprintf(fout,"\nTimeCountDeltaHS= %f",SSdata.TimeCountDeltaHS);

	// Layer Data
	fprintf(fout, "\nYYYY MM DD HH MI Layer_Thick  T  Vol_Frac_I  Vol_Frac_W  Vol_Frac_V ");
	fprintf(fout, " Vol_Frac_S Rho_S Conduc_S HeatCapac_S  rg  rb  dd  sp  mk mass_hoar ne");
	fprintf(fout, " CDot metamo");
	for (ii = 0; ii < Xdata.number_of_solutes; ii++) {
		fprintf(fout, " cIce cWater cAir  cSoil");
	}
	for (e = 0; e < Xdata.getNumberOfElements(); e++) {
		int YYYY, MM, DD, hh, mm;
		EMS[e].depositionDate.getDate(YYYY, MM, DD, hh, mm);
		fprintf(fout, "\n%04d %02d %02d %02d %02d", YYYY, MM, DD, hh, mm);
		fprintf(fout, " %7.5f %8.4f %7.5f %7.5f %7.5f",  EMS[e].L,
			Xdata.Ndata[e+1].T, EMS[e].theta[ICE], EMS[e].theta[WATER], EMS[e].theta[AIR]);
		fprintf(fout," %7.4f %6.1f %4.1f %6.1f %5.2f %5.2f %5.2f %5.2f %4d %7.5f 1",
			EMS[e].theta[SOIL], EMS[e].soil[SOIL_RHO], EMS[e].soil[SOIL_K], EMS[e].soil[SOIL_C],
			EMS[e].rg, EMS[e].rb, EMS[e].dd,  EMS[e].sp, EMS[e].mk, Xdata.Ndata[e+1].hoar);
		fprintf(fout," %10.3e %f", EMS[e].CDot, EMS[e].metamo);
		for (ii = 0; ii < Xdata.number_of_solutes; ii++) {
			fprintf(fout, "  %f %9.7f %9.7f %9.7f", EMS[e].conc(ICE,ii), EMS[e].conc(WATER,ii),
			        EMS[e].conc(AIR,ii), EMS[e].conc(SOIL,ii));
		}
	}

	// Print out the hoar hazard data info, contained in Zdata (needed for flat field only)
	fprintf(fout,"\nSurfaceHoarIndex\n");
	for (ii = 0; ii < 48; ii++) {
		fprintf(fout," %f ", Zdata.hoar24[ii]);
	}
	// Print out the drift hazard data info
	fprintf(fout,"\nDriftIndex\n");
	for (ii = 0; ii < 48; ii++) {
		fprintf(fout," %f ", Zdata.drift24[ii]);
	}
	// Print out the 3 hour new snowfall hazard data info
	fprintf(fout,"\nThreeHourNewSnow\n");
	for (ii = 0; ii < 144; ii++) {
		fprintf(fout," %f ", Zdata.hn3[ii]);
	}
	// Print out the 24 hour new snowfall hazard data info
	fprintf(fout,"\nTwentyFourHourNewSnow\n");
	for (ii = 0; ii < 144; ii++) {
		fprintf(fout," %f ", Zdata.hn24[ii]);
	}
	fprintf(fout, "\nEnd");

	fclose(fout);
}

std::string AsciiIO::getFilenamePrefix(const std::string& fnam, const std::string& path, const bool addexp)
{
	//TODO: read only once (in constructor)
	string filename_prefix = path + "/" + fnam;

	if (addexp && (experiment != "NO_EXP")) //in operational mode, nothing is appended
		filename_prefix += "_" + experiment; // complete filename_prefix

	return filename_prefix;
}

/**
 * @brief Write the Snow Profile Results, snow depth being taken VERTICALLY
 * Prepare Output File for JAVA visualization (SNOWPACK format, *.pro)
 * NOTE Parameters marked by an asterisk are available in RESEARCH visualisation only!
 * @version 11.03
 * @param i_date the current date
 * @param Xdata
 * @param Hdata
 */
void AsciiIO::writeProfile(const mio::Date& i_date, SnowStation& Xdata, const ProcessDat& Hdata)
{
	FILE *PFile=NULL;

	string stationname = Xdata.meta.getStationName();
	string filename = getFilenamePrefix(Xdata.meta.getStationID(), outpath) + ".pro";

	unsigned int e, nz;
	double cos_sl;

	const unsigned int nN = Xdata.getNumberOfNodes();
	const unsigned int nE = nN-1;
	const vector<ElementData>& EMS = Xdata.Edata;
	const vector<NodeData>& NDS = Xdata.Ndata;

	//Check whether file exists, if so check whether data can be appended
	//or file needs to be deleted
	if (IOUtils::fileExists(filename)) {
		bool append = appendFile(filename, i_date, "pro");
		if (!append && remove(filename.c_str()) != 0)
			prn_msg(__FILE__, __LINE__, "msg-", Date(), "Could not work on file %s", filename.c_str());
	}

	if (!checkHeader(filename.c_str(), "[STATION_PARAMETERS]", Hdata, "pro", &Xdata)) {
		prn_msg(__FILE__, __LINE__, "err", i_date,"Checking header in file %s", filename.c_str());
		throw IOException("Cannot dump profile " + filename + " for Java Visualisation", AT);
	}

	if (!(PFile = fopen(filename.c_str(), "a"))) {
		prn_msg(__FILE__, __LINE__, "err", i_date,
			   "Cannot open profile series file: %s", filename.c_str());
		throw IOException("Cannot dum profile " + filename + "for Java Visualisation", AT);
	}

	fprintf(PFile,"\n0500,%s", i_date.toString(Date::DIN).c_str());
	cos_sl = cos(DEG_TO_RAD(Xdata.meta.getSlopeAngle()));

	if (useSoilLayers)
		nz = nN;
	else
		nz = nE;
	//  501: height [> 0: top, < 0: bottom of elem.] (cm)
	fprintf(PFile,"\n0501,%d", nz);
	if (nz < 1) {
		// no soil and no snow
		fclose(PFile);
		return;
	}
	for (e = nN-nz; e < nN; e++)
		fprintf(PFile,",%.2f",M_TO_CM((NDS[e].z+NDS[e].u - NDS[Xdata.SoilNode].z)/cos_sl));
	//  502: element density (kg m-3)
	fprintf(PFile,"\n0502,%d", nE);
	for (e = 0; e < nE; e++)
		fprintf(PFile,",%.1f",EMS[e].Rho);
	//  503: element temperature (degC)
	fprintf(PFile,"\n0503,%d", nE);
	for (e = 0; e < nE; e++)
		fprintf(PFile,",%.2f",K_TO_C(EMS[e].Te));
	//  506: liquid water content by volume (%)
	fprintf(PFile,"\n0506,%d", nE);
	for (e = 0; e < nE; e++)
		fprintf(PFile,",%.1f",100.*EMS[e].theta[WATER]);
	// *508: dendricity (1)
	fprintf(PFile,"\n0508,%d", nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++)
		fprintf(PFile,",%.2f",EMS[e].dd);
	// *509: sphericity (1)
	fprintf(PFile,"\n0509,%d", nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++)
		fprintf(PFile,",%.2f",EMS[e].sp);
	// *510: coordination number (1)
	fprintf(PFile,"\n0510,%d", nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++)
		fprintf(PFile,",%.1f",EMS[e].N3);

	// *511: bond size (mm)
	fprintf(PFile,"\n0511,%d", nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++)
		fprintf(PFile,",%.2f",2.*EMS[e].rb);
	//  512: grain size (mm)
	fprintf(PFile,"\n0512,%d", nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++)
		fprintf(PFile,",%.2f",2.*EMS[e].rg);
	//  513: grain type (Swiss code F1F2F3)
	fprintf(PFile,"\n0513,%d", nE+1-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++)
		fprintf(PFile,",%03d",EMS[e].type);
	// surface hoar at surface? (depending on boundary conditions)
	if (M_TO_MM(NDS[nN-1].hoar/hoar_density_surf) > hoar_min_size_surf)
		fprintf(PFile,",660");
	else
		fprintf(PFile,",0");
	// *515: ice volume fraction (%)
	fprintf(PFile,"\n0515,%d", nE);
	for (e = 0; e < nE; e++)
		fprintf(PFile,",%.0f",100.*EMS[e].theta[ICE]);
	// *516: air volume fraction (%)
	fprintf(PFile,"\n0516,%d", nE);
	for (e = 0; e < nE; e++)
		fprintf(PFile,",%.0f",100.*EMS[e].theta[AIR]);
	// *517: stress (kPa)
	fprintf(PFile,"\n0517,%d", nE);
	for (e = 0; e < nE; e++)
		fprintf(PFile,",%.3e",1.e-3*EMS[e].C);
	// *518: viscosity (GPa s)
	fprintf(PFile,"\n0518,%d", nE);
	for (e = 0; e < nE; e++)
		fprintf(PFile,",%.3e",1.e-9*EMS[e].k[SETTLEMENT]);
	// *519: soil volume fraction (%)
	fprintf(PFile,"\n0519,%d", nE);
	for (e = 0; e < nE; e++)
		fprintf(PFile,",%.0f",100.*EMS[e].theta[SOIL]);
	// *520: temperature gradient (K m-1)
	fprintf(PFile,"\n0520,%d", nE);
	for (e = 0; e < nE; e++)
		fprintf(PFile,",%.3e",EMS[e].gradT);
	// *521: thermal conductivity (W K-1 m-1)
	fprintf(PFile,"\n0521,%d", nE);
	for (e = 0; e < nE; e++)
		fprintf(PFile,",%.3e",EMS[e].k[TEMPERATURE]);
	// *522: absorbed shortwave radiation (W m-2)
	fprintf(PFile,"\n0522,%d", nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++)
		fprintf(PFile,",%.1f",EMS[e].sw_abs);
	// *523: viscous deformation rate (1.e-6 s-1)
	fprintf(PFile,"\n0523,%d", nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++)
		fprintf(PFile,",%.1f",1.e6*EMS[e].EvDot);
	//  530: position (cm) and minimum stability indices
	fprintf(PFile,"\n0530,%d", 8);
	fprintf(PFile,",%d,%d,%.1f,%.2f,%.1f,%.2f,%.1f,%.2f", Xdata.S_class1, Xdata.S_class2, M_TO_CM(Xdata.z_S_d/cos_sl), Xdata.S_d, M_TO_CM(Xdata.z_S_n/cos_sl), Xdata.S_n, M_TO_CM(Xdata.z_S_s/cos_sl), Xdata.S_s);
	//  531: deformation rate stability index Sdef
	fprintf(PFile,"\n0531,%d" ,nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++)
		fprintf(PFile,",%.2f",EMS[e].S_dr);
	//  532: natural stability index Sn38
	fprintf(PFile,"\n0532,%d" ,nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode;  e < nE; e++)
		fprintf(PFile,",%.2f",NDS[e+1].S_n);
	//  533: stability index Sk38
	fprintf(PFile,"\n0533,%d" ,nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++)
		fprintf(PFile,",%.2f",NDS[e+1].S_s);
	//  534: hand hardness ...
	fprintf(PFile,"\n0534,%d" ,nE-Xdata.SoilNode);
	if (AsciiIO::r_in_n) { // ... either converted to newtons according to Swiss scale
		for (e = Xdata.SoilNode; e < nE; e++)
			fprintf(PFile,",%.1f",-1.*(19.472*pow(EMS[e].hard, 2.3607)));
	} else { // ... or in index steps (1)
		for (e = Xdata.SoilNode; e < nE; e++)
			fprintf(PFile,",%.1f", -EMS[e].hard);
	}
	//  535: inverse texture index ITI (Mg m-4)
	fprintf(PFile,"\n0535,%d" ,nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++) {
	 	if (EMS[e].dd < 0.005)
			fprintf(PFile,",%.1f",-1.*EMS[e].Rho/(2.*MM_TO_M(EMS[e].rg)));
		else
			fprintf(PFile,",%.1f",0.0);
	}

	if (variant == "CALIBRATION")
		writeFreeProfileCALIBRATION(Xdata, PFile);
	else
		writeFreeProfileDEFAULT(Xdata, PFile);

	fclose(PFile);
}

/**
 * @brief Default: dump special profiles to *.pro output file
 * @author Charles Fierz
 * @version 10.04
 * @param Xdata
 * @param *fout Output file
 */
void AsciiIO::writeFreeProfileDEFAULT(SnowStation& Xdata, FILE *fout)
{
	unsigned int e, ii, jj;
	const unsigned int nE = Xdata.getNumberOfElements();
	const vector<ElementData>& EMS = Xdata.Edata;
	const vector<NodeData>& NDS = Xdata.Ndata;

	if (out_load) {
		// *6nn: e.g. solute concentration
		for (jj = 2; jj < N_COMPONENTS-1; jj++) {
			for (ii = 0; ii < Xdata.number_of_solutes; ii++) {
				fprintf(fout,"\n06%02d,%d" , 10*jj + ii,nE-Xdata.SoilNode);
				for (e = Xdata.SoilNode; e < nE; e++) {
					fprintf(fout,",%.1f",EMS[e].conc(ii,jj));
				}
			}
		}
	} else {
		// 600-profile specials
		// *601: snow shear strength (kPa)
		fprintf(fout,"\n0601,%d" ,nE-Xdata.SoilNode);
		for (e = Xdata.SoilNode; e < nE; e++)
			fprintf(fout,",%.2f",EMS[e].s_strength);
		// *602: grain size difference (mm)
		fprintf(fout,"\n0602,%d" ,nE-Xdata.SoilNode);
		for (e = Xdata.SoilNode; e < nE-1; e++)
			fprintf(fout,",%.2f",2.*fabs(EMS[e].rg - EMS[e+1].rg));
		fprintf(fout,",0.");
		// *603: hardness difference (1)
		fprintf(fout,"\n0603,%d" ,nE-Xdata.SoilNode);
		for (e = Xdata.SoilNode; e < nE-1; e++)
			fprintf(fout,",%.2f",fabs(EMS[e].hard - EMS[e+1].hard));
		fprintf(fout,",0.");
		//  *604: ssi index
		fprintf(fout,"\n0604,%d" ,nE-Xdata.SoilNode);
		for (e = Xdata.SoilNode; e < nE-1; e++)
			fprintf(fout,",%.2f",NDS[e+1].ssi);
		fprintf(fout,",%.2f",Stability::max_stability);
	}
}

/**
 * @brief Calibration: dump special profiles to *.pro output file
 * @author Charles Fierz
 * @version 10.04
 * @param Xdata
 * @param *fout Output file
 */
void AsciiIO::writeFreeProfileCALIBRATION(SnowStation& Xdata, FILE *fout)
{
	unsigned int e;
	const unsigned int nE = Xdata.getNumberOfElements();

	const vector<ElementData>& EMS = Xdata.Edata;
	const vector<NodeData>& NDS = Xdata.Ndata;
	// 600-profile specials
	// *601: snow shear strength (kPa)
	fprintf(fout,"\n0601,%d",nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++)
		fprintf(fout,",%.2f",EMS[e].s_strength);
	// *602: grain size difference (mm)
	fprintf(fout,"\n0602,%d",nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE-1; e++)
		fprintf(fout,",%.2f",2.*fabs(EMS[e].rg - EMS[e+1].rg));
	fprintf(fout,",0.");
	// *603: hardness difference (1)
	fprintf(fout,"\n0603,%d",nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE-1; e++)
		fprintf(fout,",%.2f",fabs(EMS[e].hard - EMS[e+1].hard));
	fprintf(fout,",0.");

	// 700-profile specials for settling comparison
	// *701: SNOWPACK: settling rate due to metamorphism (sig0) (% h-1)
	fprintf(fout,"\n0701,%d",nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++)
		fprintf(fout, ",%.2f", -100.*H_TO_S(NDS[e].f));
	// *702: SNOWPACK: reaction to overload (% h-1) //ratio -Sig0 to load EMS[e].C (1)
	fprintf(fout,"\n0702,%d",nE-Xdata.SoilNode);
	for(e=Xdata.SoilNode; e<nE; e++)
		fprintf(fout,",%.2f", -100.*H_TO_S(EMS[e].EDot));
	// *703: SNOWPACK: settling rate due to load (% h-1)
	fprintf(fout,"\n0703,%d",nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++)
		fprintf(fout, ",%.2f", -100.*H_TO_S(NDS[e].udot));
	// *704: SNOWPACK: total settling rate (% h-1)
	fprintf(fout,"\n0704,%d",nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++)
		fprintf(fout,",%.2f", -100.*H_TO_S(EMS[e].EvDot));
	// *705: SNOWPACK: bond to grain ratio (1)
	fprintf(fout,"\n0705,%d",nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++)
		fprintf(fout,",%.4f", EMS[e].rb / EMS[e].rg);
	// *706: SNOWPACK: addLoad to load (%)
	fprintf(fout,"\n0706,%d",nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++)
		fprintf(fout,",%.4f", 100.*EMS[e].S);
	// SNTHERM.89
	// *891: SNTHERM: settling rate due to load (% h-1)
	fprintf(fout,"\n0891,%d" ,nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++) {
		const double eta_sntherm = (3.6e6*exp(0.08*(273.15-EMS[e].Te))*exp(0.021*EMS[e].Rho));
		fprintf(fout,",%.2f", -100.*H_TO_S(EMS[e].C/eta_sntherm));
	}
	// *892: SNTHERM: settling rate due to metamorphism (% h-1)
	fprintf(fout,"\n0892,%d" ,nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++) {
		double evdot = -2.778e-6*exp(-0.04*(273.15 - EMS[e].Te));
		if (EMS[e].Rho > 150.)
			evdot *= exp(-0.046*(EMS[e].Rho-150.));
		if( EMS[e].theta[WATER] > 0.01 )
			evdot *= 2.;
		fprintf(fout, ",%.2f", -100.*H_TO_S(evdot));
	}
	// *893: SNTHERM: viscosity (GPa s)
	fprintf(fout,"\n0893,%d" ,nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++) {
		const double eta_sntherm = (3.6e6*exp(0.08*(273.15-EMS[e].Te))*exp(0.021*EMS[e].Rho));
		fprintf(fout,",%.2f", 1.e-9*eta_sntherm);
	}
}


/**
 * @brief Dumps modelled (and measured) temperature at a given vertical height/depth z_vert (m) \n
 * Dumps also vertical height (cm) in case of fixed settling rate sensors
 * @author Charles Fierz
 * @version 10.05
 * @param *fout Output file
 * @param T Measured temperature (K)
 * @param z_vert Position of sensor measured vertically (m)
 * @param i Sensor number
 * @param *Xdata
 * @return Number of items dumped to file
 */
unsigned int AsciiIO::writeTemperatures(FILE *fout, const double& z_vert, const double& T,
                                        const unsigned int& ii, const SnowStation& Xdata)
{
	unsigned int jj=2;
	double perp_pos, temp;

	//HACK:
	/// @brief Initial height of snow needed to compute sensor position from ground if FIXED_RATES is set
	double INITIAL_HS=0;

	if (ii < number_fixed_heights) {
		perp_pos = compPerpPosition(z_vert, Xdata.cH, Xdata.Ground, Xdata.meta.getSlopeAngle());
	} else {
		if ((perp_pos = compPerpPosition(z_vert, INITIAL_HS, Xdata.Ground, Xdata.meta.getSlopeAngle()))
			    == Constants::nodata) {
			fprintf(fout, ",");
		} else {
			fprintf(fout, ",%.2f", M_TO_CM(perp_pos)/cos(DEG_TO_RAD(Xdata.meta.getSlopeAngle())));
		}
		jj++;
	}
	temp = Xdata.getModelledTemperature(perp_pos);
	fprintf(fout, ",%.2f", temp);
	if (ii < number_meas_temperatures) {
		temp = checkMeasuredTemperature(T, perp_pos, Xdata.mH);
		fprintf(fout,",%.2f", temp);
	} else {
		fprintf(fout, ",");
	}
	return jj;
}

/**
 * @brief Returns sensor position perpendicular to slope (m) \n
 * Negative vertical height indicates depth from either snow or ground surface \n
 * NOTE: Depth from snow surface cannot be used with SNP_SOIL set
 * @author Charles Fierz
 * @version 10.02
 * @param z_vert Vertical position of the sensor (m)
 * @param hs_ref Height of snow to refer to (m)
 * @param Ground Ground level (m)
 * @param slope_angle (deg)
 */
double AsciiIO::compPerpPosition(const double& z_vert, const double& hs_ref, const double& ground, const double& slope_angle)
{
	if (z_vert == Constants::nodata) {
		return Constants::nodata;
	} else if (!useSoilLayers && (z_vert < 0.)) {
		return (MAX(ground, hs_ref + z_vert * cos(DEG_TO_RAD(slope_angle))));
	} else {
		return (ground + z_vert * cos(DEG_TO_RAD(slope_angle)));
	}
}

/**
 * @brief Checks whether measured internal snow or/and soil temperature (instantaneous value) is valid \n
 * The temperature defaults to Constants::nodata if
 *  - the sensor is not covered by more than min_depth_subsurf snow (measured perpendicular to slope)
 * @author Charles Fierz
 * @version 10.01
 * @param T Measured temperature (K)
 * @param z Sensor position perpendicular to slope (m)
 * @param mH Measured snow height (m)
 * @return Measured temperature (degC) if OK, Constants::nodata else
 */
double AsciiIO::checkMeasuredTemperature(const double& T, const double& z, const double& mH)
{
	if ((z <= (mH - min_depth_subsurf)) && (T != Constants::nodata))
		return K_TO_C(T);
	else
		return Constants::nodata;
}

/**
 * @brief Find element with corresponding tag or return -1 if not found
 * @version 10.04
 * @param tag to look for
 * @param Xdata
 * @return Index of tagged element, -1 if not found
 */
int AsciiIO::findTaggedElement(const unsigned int& tag, const SnowStation& Xdata)
{
	for (unsigned int e=0; e<Xdata.getNumberOfElements(); e++) {
		if (Xdata.Edata[e].mk/100 == tag)
			return e;
	}
	return -1;
}

/**
 * @brief Dumps modelled and measured temperature for tag(ged layer)
 * @author Charles Fierz
 * @version 10.02
 * @param *fout Output file
 * @param tag Tag number;
 * @param *Mdata
 * @param *Xdata
 * @return Number of dumped values
 */
unsigned int AsciiIO::writeHeightTemperatureTag(FILE *fout, const unsigned int& tag,
                                                const CurrentMeteo& Mdata, const SnowStation& Xdata)
{
	unsigned int jj = 2;
	const unsigned int ii = number_fixed_heights + number_fixed_rates + (tag-1);
	int e;
	double perp_pos, temp;
	if ((e = findTaggedElement(tag, Xdata)) >= 0) {
		perp_pos = ((Xdata.Ndata[e].z + Xdata.Ndata[e].u + Xdata.Ndata[e+1].z
		                + Xdata.Ndata[e+1].u)/2. - Xdata.Ground);
		fprintf(fout,",%.2f,%.2f", M_TO_CM(perp_pos)
		            / cos(DEG_TO_RAD(Xdata.meta.getSlopeAngle())), K_TO_C(Xdata.Edata[e].Te));
	} else {
		fprintf(fout,",,%.2f", Constants::nodata);
	}
	if (ii < number_meas_temperatures) {
		if ((perp_pos = compPerpPosition(Mdata.zv_ts[ii], Xdata.cH, Xdata.Ground, Xdata.meta.getSlopeAngle()))
			    == Constants::nodata) {
			fprintf(fout,",,%.2f", Constants::nodata);
		} else {
			fprintf(fout,",%.2f", M_TO_CM(perp_pos)/cos(DEG_TO_RAD(Xdata.meta.getSlopeAngle())));
			temp = checkMeasuredTemperature(Mdata.ts[ii], perp_pos, Xdata.mH);
			fprintf(fout,",%.2f", temp);
		}
		jj += 2;
	}
	return jj;
}

/**
 * @brief Parse through a met file and read the last date for which meteo data has been written
 * @param eoln A char that represents the end of line character
 * @param start_date Holds the start date of this simulation
 * @param fin The file input stream to use
 * @param fout The output stream of a temporary file
 * @return TRUE if file may be appended, false if file needs to be overwritten
 */
bool AsciiIO::parseMetFile(const char& eoln, const mio::Date& start_date, std::istream& fin, std::ostream& ftmp)
{
	string tmpline = "";
	vector<string> vecTmp;
	Date current_date;

	bool append       = false; //true if file may be appended and false otherwise
	bool insert_endl  = false;
	bool data_started = false;

	do { //Loop going through the lines of the file
		getline(fin, tmpline, eoln); //read complete line

		if (data_started) {
			if (tmpline.length() > 20) {//the last line is without a carriage return
				IOUtils::trim(tmpline);
				IOUtils::readLineToVec(tmpline, vecTmp, ',');
				if ((vecTmp.size() >= 2) && (vecTmp[1].length() >= 16)) {
					string tmpdate = vecTmp[1].substr(6,4) + "-" + vecTmp[1].substr(3,2) + "-" + vecTmp[1].substr(0,2)
						+ "T" + vecTmp[1].substr(11,2) + ":" + vecTmp[1].substr(14,2);
					IOUtils::convertString(current_date, tmpdate, time_zone);

					if (current_date.getJulianDate() < (start_date.getJulianDate()-0.00001)) {
						append=true;
					} else {
						break; //the start date of the simulation is newer/equal than current_date
					}
				}
			}
		} else {
			IOUtils::trim(tmpline);
			if (tmpline == "[DATA]")
				data_started = true;
		}

		if (insert_endl)
			ftmp << endl;
		else
			insert_endl = true;

		ftmp << tmpline; //copy line to tmpfile

	} while(!fin.eof());

	return append;
}

/**
 * @brief Parse through a pro file and check whether it can be appended for current simulation
 * @param eoln A char that represents the end of line character
 * @param start_date Holds the start date of this simulation
 * @param fin The file input stream to use
 * @param fout The output stream of a temporary file
 * @return TRUE if file may be appended, false if file needs to be overwritten
 */
bool AsciiIO::parseProFile(const char& eoln, const mio::Date& start_date, std::istream& fin, std::ostream& ftmp)
{
	string tmpline = "";
	vector<string> vecTmp;
	Date current_date;

	bool append = false; //true if file may be appended and false otherwise
	bool insert_endl = false;
	bool data_started = false;

	do { //Loop going through the lines of the file
		getline(fin, tmpline, eoln); //read complete line
		IOUtils::readLineToVec(tmpline, vecTmp, ',');

		if (data_started) {
			if (vecTmp.size() >= 2) {
				if (vecTmp[0] == "0500"){ //The date tag
					if (vecTmp[1].length() >= 16) {
						string tmpdate = vecTmp[1].substr(6,4) + "-" + vecTmp[1].substr(3,2) + "-" + vecTmp[1].substr(0,2)
						                 + "T" + vecTmp[1].substr(11,2) + ":" + vecTmp[1].substr(14,2);
						IOUtils::convertString(current_date, tmpdate, time_zone);

						if (current_date.getJulianDate() < (start_date.getJulianDate()-0.00001)){
							append=true;
						} else {
							break; //the start date of the simulation is newer/equal than current_date
						}
					}
				}
			}
		} else {
			IOUtils::trim(tmpline);
			if (tmpline == "[DATA]") data_started = true;
		}

		if (insert_endl)
			ftmp << endl;
		else
			insert_endl = true;

		ftmp << tmpline; //copy line to tmpfile
	} while( !fin.eof() );

	return append;
}

/**
 * @brief Check whether data can be appended to a file, or whether file needs to be deleted and recreated
 *        The following logic is implemented if the file already contains data:
 *        - if the startdate lies before the data written in the file, overwrite the file
 *        - if the startdate lies within the data in the file then append from that date on, delete the rest
 *        - if the startdate is after the data in the file then simply append
 * @param filename The file to check (must exist)
 * @param startdate The start date of the data to be written
 * @param type A string representing the type of file, i.e. "pro" or "met"
 * @return A boolean, true if file can be appended, false otherwise
 */
bool AsciiIO::appendFile(const std::string& filename, const mio::Date& startdate, const std::string& ftype)
{
	//Check if file has already been checked
	set<string>::const_iterator it = setAppendableFiles.find(filename);
	if (it != setAppendableFiles.end()) //file was already checked
		return true;

	/* Go through file and parse meteo data date if current date
	 * is newer than the last one in the file, appending is possible
	 */
	ifstream fin;
	ofstream fout; //for the tmp file
	const string filename_tmp = filename + ".tmp";

	fin.open (filename.c_str());
	fout.open(filename_tmp.c_str());

	if (fin.fail()) throw FileAccessException(filename, AT);
	if (fout.fail()) throw FileAccessException(filename_tmp, AT);

	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	try {
		bool append_possible = false; //the temporary file will be copied

		if (ftype == "pro") {
			append_possible = parseProFile(eoln, startdate, fin, fout);
		} else if (ftype == "met") {
			append_possible = parseMetFile(eoln, startdate, fin, fout);
		}

		fin.close();
		fout.close();

		if (append_possible)
			IOUtils::copy_file(filename_tmp, filename);

		remove(filename_tmp.c_str()); //delete temporary file

		setAppendableFiles.insert(filename); //remember, that this file has been checked already
		return append_possible;
	} catch(...) {
		if (fin.is_open())  fin.close();
		if (fout.is_open()) fout.close();
		return false;
	}
}

/**
 * @brief Write all Time Series results (*.met)
 * All depths and water equivalents (mass) are taken VERTICALLY. \n
 * If AVGSUM_TIME_SERIES is set, mean fluxes and cumulated masses since last dump are written, \n
 * else current energy fluxes, cumulated masses over last computation_step_length (recommended setting in operational mode).
 * If CUMSUM_MASS is set, current value of cumulated masses since begin of run are dumped. \n
 * Precipitations (rain& snow, rain) are always dumped as rates (kg m-2 h-1). \n
 * NOTE:
 * 	-# neither AVGSUM_TIME_SERIES nor CUMSUM_MASS can be set if NUMBER_SLOPES > 1.
 * 	-# When running SNOW_REDISTRIBUTION on virtual slopes, eroded mass will be dumped
 *     to the windward *.met file and added to the solid precipitations of the lee *.met file!
 * \li DO NOT change the order of parameters below! Additional parameters may be dumped at pos.
 *     93[94] to 100 in writeFreeSeriesXXX()
 * @version 11.03
 * @param Xdata
 * @param Sdata
 * @param Mdata
 * @param Hdata
 * @param wind_trans24 eroded snow from flat field (present time step) or virtual windward slope (previous time step)
 */
void AsciiIO::writeTimeSeries(const SnowStation& Xdata, const SurfaceFluxes& Sdata,
                              const CurrentMeteo& Mdata, const ProcessDat& Hdata,
                              const double wind_trans24)
{
	FILE *TFile=NULL;

	string stationname = Xdata.meta.getStationName();
	string filename = getFilenamePrefix(Xdata.meta.getStationID(), outpath) + ".met";

	const vector<NodeData>& NDS = Xdata.Ndata;
	const int nN = Xdata.getNumberOfNodes();
	double cos_sl = cos(DEG_TO_RAD(Xdata.meta.getSlopeAngle()));

	//Check whether file exists, if so check whether data can be appended or file needs to be deleted
	if (IOUtils::fileExists(filename)) {
		bool append = appendFile(filename, Mdata.date, "met");
		if (!append && remove(filename.c_str()) != 0)
			prn_msg(__FILE__, __LINE__, "msg-", Date(), "Could not work on file %s", filename.c_str());
	}

	// Check file for header
	if (!checkHeader(filename.c_str(), "[STATION_PARAMETERS]", Hdata, "met", &Xdata)) {
		prn_msg(__FILE__, __LINE__, "err", Mdata.date, "Checking header in file %s", filename.c_str());
		throw InvalidFormatException("Writing Time Series data failed", AT);
	}

	if (!(TFile = fopen(filename.c_str(), "a"))) {
		prn_msg(__FILE__, __LINE__, "err", Mdata.date, "Cannot open time series file: %s", filename.c_str());
		throw FileAccessException(filename, AT);
	}
	// Print time stamp
	fprintf(TFile,"\n0203,%s", Mdata.date.toString(Date::DIN).c_str());
	if (out_heat)
		// 1-2: Turbulent fluxes (W m-2)
		fprintf(TFile,",%f,%f", Sdata.qs, Sdata.ql);
	else
		fprintf(TFile,",,");
	if (out_lw)
		// 3-5: Longwave radiation fluxes (W m-2)
		fprintf(TFile,",%f,%f,%f", Sdata.lw_out, Sdata.lw_in, Sdata.lw_net);
	else
		fprintf(TFile,",,,");
	if (out_sw)
		// 6-9: Shortwave radiation fluxes (W m-2) and computed albedo (1)
		fprintf(TFile,",%f,%f,%f,%f", Sdata.sw_out, Sdata.sw_in, Sdata.qw, Sdata.cA);
	else
		fprintf(TFile,",,,,");
	if (out_meteo)
		// 10-13: Air temperature, snow surface temperature (modeled and measured), temperature at bottom of snow/soil pack (degC)
		fprintf(TFile,",%f,%f,%f,%f", K_TO_C(Mdata.ta), K_TO_C(NDS[nN-1].T), K_TO_C(Mdata.tss), K_TO_C(NDS[0].T));
	else
		fprintf(TFile,",,,,");
	if (out_heat)
		// 14-17: Heat flux at lower boundary (W m-2), ground surface temperature (degC),
		//        Heat flux at gound surface (W m-2), rain energy (W m-2)
		fprintf(TFile,",%f,%f,%f,%f", Sdata.qg, K_TO_C(NDS[Xdata.SoilNode].T), Sdata.qg0, Sdata.qr);
	else
		fprintf(TFile,",,,,");
	if (out_sw)
		// 18-22: projected solar radiation (W m-2), meas. albedo (1)
		fprintf(TFile,",%f,%f,%f,%f,%f", Sdata.sw_hor, Sdata.sw_in, Sdata.sw_dir, Sdata.sw_diff, Sdata.mA);
	else
		fprintf(TFile,",,,,,");
	if (out_meteo)
		// 23-26: rH (%), wind (m s-1), wind_drift (m s-1), wind_dir (deg),
		// 27: solid precipitation rate (kg m-2 h-1),
		// 28-29: modeled and maesured vertical snow depth (cm)
		fprintf(TFile,",%f,%f,%f,%f,%f,%f,%f", 100.*Mdata.rh, Mdata.vw, Mdata.vw_drift, Mdata.dw,
		          Sdata.mass[SurfaceFluxes::MS_HNW], M_TO_CM((Xdata.cH - Xdata.Ground)/cos_sl),
				    M_TO_CM(Mdata.hs1/cos_sl));
	else
		fprintf(TFile,",,,,,,,");
	if (out_haz) {
		// 30-33: surface hoar size (mm), 24h drift index (cm), height of new snow HN (cm), 3d sum of daily new snow depths (cm)
		if (!perp_to_slope)
			fprintf(TFile,",%f,%f,%f,%f", Hdata.hoar_size, wind_trans24, Hdata.hn24, Hdata.hn72_24);
		else
			// dump vertical values if PERP_TO_SLOPE
			fprintf(TFile,",%f,%f,%f,%f", Hdata.hoar_size, wind_trans24, Hdata.hn24/cos_sl, Hdata.hn72_24/cos_sl);
	} else {
		fprintf(TFile,",,,,");
	}
	if (out_mass) {
		// 34-39: total mass, eroded mass, rain rate, runoff at bottom of snowpack, sublimation and evaporation, all in kg m-2 except rain as rate: kg m-2 h-1; see also 51-52 & 93
		fprintf(TFile,",%f,%f,%f,%f,%f,%f", Sdata.mass[SurfaceFluxes::MS_TOTALMASS]/cos_sl,
		        Sdata.mass[SurfaceFluxes::MS_WIND]/cos_sl, Sdata.mass[SurfaceFluxes::MS_RAIN],
		          Sdata.mass[SurfaceFluxes::MS_RUNOFF]/cos_sl, Sdata.mass[SurfaceFluxes::MS_SUBLIMATION]/cos_sl,
		            Sdata.mass[SurfaceFluxes::MS_EVAPORATION]/cos_sl);
	} else {
		fprintf(TFile,",,,,,,");
	}
	// 40-49: Internal Temperature Time Series at fixed heights, modeled and measured, all in degC
	if (out_t && (number_fixed_heights || number_fixed_rates)) {
		unsigned int jj = 0;
		for (unsigned int ii = 0; ii < MIN(5, number_fixed_heights); ii++)
			jj += writeTemperatures(TFile, Mdata.zv_ts[ii], Mdata.ts[ii], ii, Xdata);
		for (; jj < 10; jj++)
			fprintf(TFile,",");
	} else {
		fprintf(TFile,",,,,,,,,,,");
	}
	if (max_number_sensors == 5) {
		if (out_load)
		// 50: Solute load at ground surface
			fprintf(TFile,",%f",Sdata.load[0]);
		else
			fprintf(TFile,",");
		if (out_mass)
		// 51-52: SWE (for checks) and LWC (kg m-2); see also 34-39
			fprintf(TFile,",%f,%f", Sdata.mass[SurfaceFluxes::MS_SWE]/cos_sl,
			        Sdata.mass[SurfaceFluxes::MS_WATER]/cos_sl);
		else
			fprintf(TFile,",,");
		if (out_stab) {
			// 53-64: Stability Time Series, heights in cm
			fprintf(TFile,",%d,%d,%.1f,%.2f,%.1f,%.2f,%.1f,%.2f,%.1f,%.2f,%.1f,%.2f",
			        Xdata.S_class1, Xdata.S_class2, M_TO_CM(Xdata.z_S_d/cos_sl), Xdata.S_d,
			          M_TO_CM(Xdata.z_S_n/cos_sl), Xdata.S_n, M_TO_CM(Xdata.z_S_s/cos_sl), Xdata.S_s,
			            M_TO_CM(Xdata.z_S_4/cos_sl), Xdata.S_4, M_TO_CM(Xdata.z_S_5/cos_sl), Xdata.S_5);
		} else {
			fprintf(TFile,",,,,,,,,,,,,");
		}
		if (out_canopy && useCanopyModel)
			// 65-92 (28 columns)
			Canopy::cn_DumpCanopyData(TFile, &Xdata.Cdata, &Sdata, cos_sl);
		else
			fprintf(TFile,",,,,,,,,,,,,,,,,,,,,,,,,,,,,");
	} else if (out_t) {
		// 50-93 (44 columns)
		unsigned int ii, jj = 0;
		for (ii = MIN(5, number_fixed_heights); ii < number_fixed_heights+number_fixed_rates; ii++) {
			if ((jj += writeTemperatures(TFile, Mdata.zv_ts[ii], Mdata.ts[ii], ii, Xdata)) > 44) {
				prn_msg(__FILE__, __LINE__, "err", Mdata.date,
				        "There is not enough space to accomodate your temperature sensors: j=%d > 44!", jj);
				throw IOException("Writing Time Series data failed", AT);
			}
		}
		if (Xdata.tag_low) {
			unsigned int tag = Xdata.tag_low, j_lim;
			while ( (tag + ii) <= number_sensors ) {
				if ((tag + ii) <= number_meas_temperatures)
					j_lim = 41;
				else
					j_lim = 43;
				if (jj < j_lim) {
					jj += writeHeightTemperatureTag(TFile, tag, Mdata, Xdata);
					tag++;
				} else {
					break;
				}
			}
		}
		for (; jj < 44; jj++)
			fprintf(TFile,",");
	} else {
		fprintf(TFile,",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,");
	}
	// 93[94]-100 (8 or 7 free columns)
	unsigned int nCalcSteps = 1;
	double crust = 0., dhs_corr = 0., mass_corr = 0.;
	if (!avgsum_time_series)
		nCalcSteps = (int)(ts_days_between / M_TO_D(calculation_step_length) + 0.5);
	if (out_haz) {
		crust = Hdata.crust;
		dhs_corr = Hdata.dhs_corr;
		mass_corr = Hdata.mass_corr;
	}
	if (variant == "CALIBRATION") {
		writeFreeSeriesCALIBRATION(Xdata, Sdata, Mdata, crust, dhs_corr, mass_corr, nCalcSteps, TFile);
	} else if (variant == "ANTARCTICA") {
		writeFreeSeriesANTARCTICA(Xdata, Sdata, Mdata, crust, dhs_corr, mass_corr, nCalcSteps, TFile);
	} else {
		writeFreeSeriesDEFAULT(Xdata, Sdata, Mdata, crust, dhs_corr, mass_corr, nCalcSteps, TFile);
	}

	fclose (TFile);
}

/**
 * @brief Default: last 8 time series (columns 93 to 100) dumped to *.met output file
 * @version 11.05
 * @param Xdata
 * @param Sdata
 * @param Mdata
 * @param crust height
 * @param dhs_corr correction for height of snow in (operational mode only)
 * @param mass_corr mass correction due to dhs_corr (operational mode only)
 * @param nCalcSteps between outputs
 * @param *fout Output file
 */
void AsciiIO::writeFreeSeriesDEFAULT(const SnowStation& Xdata, const SurfaceFluxes& Sdata, const CurrentMeteo& Mdata, const double crust, const double dhs_corr, const double mass_corr, const unsigned int nCalcSteps, FILE *fout)
{
	double cos_sl = cos(DEG_TO_RAD(Xdata.meta.getSlopeAngle()));
	if (useSoilLayers)
		// 93: Soil Runoff (kg m-2); see also 34-39 & 51-52
		fprintf(fout,",%f",Sdata.mass[SurfaceFluxes::MS_SOIL_RUNOFF]/cos_sl);
	else
		fprintf(fout,",");
	if (out_heat) {
		// 94: change of internal energy (kJ m-2)
		if (Xdata.getNumberOfElements() > Xdata.SoilNode)
			fprintf(fout,",%.3f", ((Xdata.dIntEnergy * nCalcSteps)
		                             - (Sdata.qg0 * D_TO_S(ts_days_between))) / 1000.);
		else
			fprintf(fout, ",%f", Constants::undefined);
		// 95: sum of energy fluxes at surface (kJ m-2)
		fprintf(fout,",%f",
		            ((Sdata.qw + Sdata.lw_net + Sdata.qs + Sdata.ql + Sdata.qr)
		                * D_TO_S(ts_days_between)) / 1000.);
	} else {
		fprintf(fout,",,");
	}
	// 96-97: new snow densities, measured and in use (kg m-3)
	if(Sdata.cRho_hn > 0.) {
		fprintf(fout,",%.1f,%.1f", Sdata.mRho_hn, Sdata.cRho_hn);
	} else {
		double mRho_hn = Constants::undefined;
		if (Mdata.rho_hn != Constants::nodata)
			mRho_hn = -Mdata.rho_hn;
		fprintf(fout,",%.1f,%.1f", mRho_hn, Sdata.cRho_hn);
	}
	// 98: crust height (S-slope) (cm)
	fprintf(fout,",%f", crust);
	// 99-100:
	if (!research_mode)
		// snow depth (cm) and mass correction (kg m-2)
		fprintf(fout,",%f,%f", M_TO_CM(dhs_corr), mass_corr);
	else
		// for example, measured turbulent fluxes (W m-2); see also 1-2
		fprintf(fout,",,");
}

/**
 * @brief Antarctic: last [8]7 time series (columns [93]94 to 100) dumped to *.met output file
 * @version 11.05
 * @param Xdata
 * @param Sdata
 * @param Mdata
 * @param crust not available; replaced by potential erosion level
 * @param dhs_corr not available
 * @param mass_corr not available
 * @param nCalcSteps between outputs
 * @param *fout Output file
 */
void AsciiIO::writeFreeSeriesANTARCTICA(const SnowStation& Xdata, const SurfaceFluxes& Sdata,
                                        const CurrentMeteo& Mdata, const double crust,
                                        const double dhs_corr, const double mass_corr,
                                        const unsigned int nCalcSteps, FILE *fout)
{
	(void) crust; (void) dhs_corr; (void) mass_corr;
	if (max_number_sensors == 5)
		fprintf(fout, ",");
	// 94: change of internal energy (kJ m-2)
	if (out_heat) {
		// 94: change of internal energy (kJ m-2)
		if (Xdata.getNumberOfElements() > Xdata.SoilNode)
			fprintf(fout,",%.3f", ((Xdata.dIntEnergy * nCalcSteps)
		                             - (Sdata.qg0 * D_TO_S(ts_days_between))) / 1000.);
		else
			fprintf(fout, ",%f", Constants::undefined);
		// 95: sum of energy fluxes at surface (kJ m-2)
		fprintf(fout,",%f",
		            ((Sdata.qw + Sdata.lw_net + Sdata.qs + Sdata.ql + Sdata.qr)
		                * D_TO_S(ts_days_between)) / 1000.);
	} else {
		fprintf(fout,",,");
	}
	// 96-97: new snow densities, measured and in use (kg m-3)
	if(Sdata.cRho_hn > 0.) {
		fprintf(fout,",%.1f,%.1f", Sdata.mRho_hn, Sdata.cRho_hn);
	} else {
		double mRho_hn = Constants::undefined;
		if (Mdata.rho_hn != Constants::nodata)
			mRho_hn = -Mdata.rho_hn;
		fprintf(fout,",%.1f,%.1f", mRho_hn, Sdata.cRho_hn);
	}
	// 98: potential erosion level below surface (cm)
	fprintf(fout,",%f", M_TO_CM(Xdata.Ndata[Xdata.ErosionLevel+1].z - Xdata.cH));
	// 99-100
	if (out_meteo)
		// mean over 100 h of air humidity (%) and mean wind speed (m s-1)
		fprintf(fout,",%.2f,%.2f", 100. * Mdata.rh_avg, Mdata.vw_avg);
	else
		fprintf(fout,",,");
}

/**
 * @brief Calibration: last [8]7 time series (columns [93]94 to 100) dumped to *.met output file
 * @version 11.05
 * @param Xdata
 * @param Sdata not used
 * @param Mdata
 * @param crust height, not used
 * @param dhs_corr not available
 * @param mass_corr not available
 * @param nCalcSteps
 * @param *fout Output file
 */
void AsciiIO::writeFreeSeriesCALIBRATION(const SnowStation& Xdata, const SurfaceFluxes& Sdata, const CurrentMeteo& Mdata, const double crust, const double dhs_corr, const double mass_corr, const unsigned int nCalcSteps, FILE *fout)
{
	(void) crust; (void) dhs_corr; (void) mass_corr; (void) Sdata;
	double rho_hn;
	const double t_surf = MIN(C_TO_K(-0.1), Xdata.Ndata[Xdata.getNumberOfNodes()-1].T);
	if (max_number_sensors == 5)
		fprintf(fout,",");
	if (out_heat) {
		// 94: change of internal energy (kJ m-2)
		if (Xdata.getNumberOfElements() > Xdata.SoilNode)
			fprintf(fout,",%.3f", ((Xdata.dIntEnergy * nCalcSteps)
			                         - (Sdata.qg0 * D_TO_S(ts_days_between))) / 1000.);
		else
			fprintf(fout, ",%f", Constants::undefined);
		// 95: sum of energy fluxes at surface (kJ m-2)
		fprintf(fout,",%f",
		            ((Sdata.qw + Sdata.lw_net + Sdata.qs + Sdata.ql + Sdata.qr)
		                * D_TO_S(ts_days_between)) / 1000.);
	} else {
		fprintf(fout,",,");
	}
	// 96-100: new snow densities: measured, in use, newLe, bellaire, and crocus (kg m-3)
	if (Sdata.cRho_hn > 0.) {
		fprintf(fout,",%.1f,%.1f", Sdata.mRho_hn, Sdata.cRho_hn);
		rho_hn = SnLaws::compNewSnowDensity(hn_density, "LEHNING_NEW", Mdata, Xdata, t_surf, variant);
		fprintf(fout,",%.1f", rho_hn);
		rho_hn = SnLaws::compNewSnowDensity(hn_density, "BELLAIRE", Mdata, Xdata, t_surf, variant);
		fprintf(fout,",%.1f", rho_hn);
		rho_hn = SnLaws::compNewSnowDensity(hn_density, "CROCUS", Mdata, Xdata, t_surf, variant);
		fprintf(fout,",%.1f", rho_hn);
	} else {
		double mRho_hn = Constants::undefined;
		if (Mdata.rho_hn != Constants::nodata)
			mRho_hn = -Mdata.rho_hn;
		fprintf(fout,",%.1f,%.1f", mRho_hn, Sdata.cRho_hn);
		rho_hn = SnLaws::compNewSnowDensity(hn_density, "LEHNING_NEW", Mdata, Xdata, t_surf, variant);
		fprintf(fout,",%.1f", -rho_hn);
		rho_hn = SnLaws::compNewSnowDensity(hn_density, "BELLAIRE", Mdata, Xdata, t_surf, variant);
		fprintf(fout,",%.1f", -rho_hn);
		rho_hn = SnLaws::compNewSnowDensity(hn_density, "CROCUS", Mdata, Xdata, t_surf, variant);
		fprintf(fout,",%.1f", -rho_hn);
	}
}

/**
 * @brief This routine:
 * -# Checks for header in fnam by testing for first_string
 * -# If header is missing:
 *    - writes header in fnam according to file type (ext)
 *    - returns -1 (ext=="none")
 * @author Charles Fierz \n Mathias Bavay
 * @version 10.02
 * @param *fnam Filename
 * @param *first_string First string to be found in header
 * @param Hdata
 * @param *ext File extension
 * @return status
 */
bool AsciiIO::checkHeader(const char *fnam, const char *first_string, const ProcessDat& Hdata,
                          const char *ext, ...)
{
	FILE *fin=NULL;
	FILE *fout=NULL;
	va_list argptr; // get an arg ptr
	SnowStation *va_Xdata;
	char dummy[MAX_STRING_LENGTH]="\000", dummy_l[MAX_LINE_LENGTH]="\000";
	unsigned int  ii, jj;

	if ((fin = fopen(fnam, "r"))) {
		// Check header of existing file
		if( fgets(dummy_l, MAX_LINE_LENGTH, fin) == NULL) {
			fclose(fin);
			std::stringstream ss;
			ss << "Can not read header of file " << fnam;
			throw InvalidFormatException(ss.str(), AT);
		}
		sscanf(dummy_l, "%s", dummy);
		if ((strcmp(dummy, first_string) != 0)) {
			prn_msg(__FILE__, __LINE__, "err", Date(), "Header in %s should read %s, not %s", fnam, first_string, dummy);
			return false;
		}
		fclose(fin);
	} else if ((strcmp(ext, "none") == 0)) {
		// Check header only!
		return false;
	} else {
		if (!(fout = fopen(fnam, "w")))
			return false;
		// Initialize argptr to point to the first argument after the ext string
		va_start(argptr, ext);

		if ((strcmp(ext, "err") == 0)) {
			fprintf(fout, "[SNOWPACK_ERROR_LOG]");
			fprintf(fout, "\n          RUNTIME :  STN LOC LINE MSG [JULIAN]");
		} else if ((strcmp(ext, "met") == 0)) {
			va_Xdata = va_arg(argptr, SnowStation *);
			string stationname = va_Xdata->meta.getStationName();
			fprintf(fout, "[STATION_PARAMETERS]");
			fprintf(fout, "\nStationName= %s",   stationname.c_str());
			fprintf(fout, "\nLatitude= %.2f",   va_Xdata->meta.position.getLat());
			fprintf(fout, "\nLongitude= %.2f",  va_Xdata->meta.position.getLon());
			fprintf(fout, "\nAltitude= %.0f",   va_Xdata->meta.position.getAltitude());
			fprintf(fout, "\nSlopeAngle= %.2f", va_Xdata->meta.getSlopeAngle());
			fprintf(fout, "\nSlopeAzi= %.2f",   va_Xdata->meta.getAzimuth());
			fprintf(fout, "\nDepthTemp= %1d",    useSoilLayers);
			for (ii = 0; ii < number_fixed_heights; ii++)
				fprintf(fout, ",%.3f", fixed_sensor_depths[ii]);
			fprintf(fout, "\n\n[HEADER]");
			if (out_haz) // HACK To avoid troubles in A3D
				fprintf(fout, "\n#%s, Snowpack %s version %s run by \"%s\"", Hdata.sn_computation_date,
				        variant.c_str(), Hdata.sn_version, Hdata.sn_user);
			fprintf(fout, "\n,,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100");
			fprintf(fout, "\nID,Date,Sensible heat,Latent heat,Outgoing longwave radiation,Incoming longwave radiation,Net absorbed longwave radiation,Reflected shortwave radiation,Incoming shortwave radiation,Net absorbed shortwave radiation,Modelled surface albedo,Air temperature,Modeled surface temperature,Measured surface temperature,Temperature at bottom of snow or soil pack,Heat flux at bottom of snow or soil pack,Ground surface temperature,Heat flux at ground surface,Heat advected to the surface by liquid precipitation,Global solar radiation (horizontal)");
			fprintf(fout, ",Global solar radiation on slope,Direct solar radiation on slope,Diffuse solar radiation on slope,Measured surface albedo,Relative humidity,Wind speed,Max wind speed at snow station or wind speed at ridge station,Wind direction at snow station,Precipitation rate at surface (solid only),Modelled snow depth (vertical),Measured snow depth (vertical),Surface hoar size,24h Drift index (vertical),Height of new snow HN (24h vertical),3d sum of daily height of new snow (vertical),Total snowpack mass,Eroded mass,Rain rate,Surface runoff (without soil infiltration)");
			fprintf(fout, ",Sublimation,Evaporation,Temperature 1 (modelled),Temperature 1 (measured),Temperature 2 (modelled),Temperature 2 (measured),Temperature 3 (modelled),Temperature 3 (measured),Temperature 4 (modelled),Temperature 4 (measured),Temperature 5 (modelled),Temperature 5 (measured)");
			if (max_number_sensors == 5) {
				fprintf(fout, ",Solute load at soil surface,SWE (of snowpack),Liquid Water Content (of snowpack),Profile type,Stability class,z_Sdef,Deformation rate stability index Sdef,z_Sn38,Natural stability index Sn38,z_Sk38,Skier stability index Sk38,z_SSI,Structural Stability index SSI,z_S5,Stability index S5");
				if (useCanopyModel && out_canopy) {
					fprintf(fout, ",Interception storage,Canopy surface  temperature,Canopy albedo,Wet fraction,Interception capacity,Net shortwave radiation absorbed by canopy,Net longwave radiation absorbed by canopy,Net radiation canopy,Sensible heat flux into the canopy,Latent heat flux into the canopy,Transpiration of the canopy,Evaporation and sublimation of interception (liquid and frozen),Interception rate,Throughfall,Snow unload,Sensible heat flux to the canopy,Latent heat flux to the canopy,Longwave radiation up above canopy,Longwave radiation down above canopy");
					fprintf(fout, ",Net longwave radiation above canopy,Shortwave radiation up above canopy,Shortwave radiation down above canopy,Net shortwave radiation above canopy,Total land surface albedo,Total net radiation,Surface (radiative) temperature,Precipitation Above Canopy,Total Evapotranspiration");
				} else {
					fprintf(fout,",-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-");
				}
			} else if (out_t) {
				int i_prn;
				jj = 0;
				for (ii = MIN(5, number_fixed_heights); ii < number_fixed_heights+number_fixed_rates; ii++) {
					if (ii < number_fixed_heights) {
						i_prn = ii + 1;
						fprintf(fout, ",Temperature %d (modelled)", i_prn);
					} else {
						i_prn = (ii-number_fixed_heights)+1;
						fprintf(fout, ",Hfr %d", i_prn);
						fprintf(fout, ",Tfr %d (modelled)", i_prn);
						jj++;
					}
					if (ii < number_meas_temperatures) {
						if (ii < number_fixed_heights) {
							fprintf(fout, ",Temperature %d (measured)", i_prn);
						} else {
							fprintf(fout, ",Tfr %d (measured)", i_prn);
						}
					} else {
						fprintf(fout, ",");
					}
					jj += 2;
				}
				if (va_Xdata->tag_low) {
					unsigned int tag = va_Xdata->tag_low, j_lim;
					while ((tag + ii) <= number_sensors) {
						if ((tag + ii) <= number_meas_temperatures) {
							j_lim = 41;
						} else {
							j_lim = 43;
						}
						if (jj < j_lim) {
							fprintf(fout, ",H(tag%02d),T(tag%02d)", tag, tag);
							jj += 2;
							if (ii < number_meas_temperatures) {
								fprintf(fout, ",H(meas%02d),T(meas%02d)", tag, tag);
								jj += 2;
							}
							tag++;
						}
					}
				}
				for (; jj < 44; jj++) {
					fprintf(fout,",-");
				}
			} else {
				fprintf(fout, ",-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-");
			}

			if (variant == "ANTARCTICA") {
				if (max_number_sensors == 5)
					fprintf(fout, ",");
				fprintf(fout, ",Internal energy change,Surface input (sum fluxes),Measured new snow density,Modeled new snow density,Erosion level (from srf),Running mean relative humidity (100h),Running mean wind speed (100h)");
			} else if (variant == "CALIBRATION") {
				if (max_number_sensors == 5)
					fprintf(fout, ",");
				fprintf(fout, "Internal energy change,Surface input (sum fluxes),rho_hn(measured),rho_hn(Zwart),rho_hn(Lehning),rho_hn(Bellaire),rho_hn(CROCUS)");
			} else {
				fprintf(fout, ",Soil runoff,Internal energy change,Surface input (sum fluxes),Measured new snow density,Modeled new snow density,Crust thickness (S-slope)");
				if (!research_mode)
					fprintf(fout, ",Snow depth correction,Mass change");
				else
					fprintf(fout, ",Measured sensible heat,Measured latent heat");
			}

			fprintf(fout, "\n,,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,1,degC,degC,degC,degC,W m-2,degC,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,1,%%,m s-1,m s-1,deg,kg m-2 h-1,cm,cm,mm,cm,cm,cm,kg m-2,kg m-2 h-1,kg m-2 h-1,kg m-2,kg m-2,kg m-2,degC,degC,degC,degC,degC,degC,degC,degC,degC,degC");
			if (max_number_sensors == 5) {
				fprintf(fout, ",kg m-2,kg m-2,kg m-2,-,-,cm,1,cm,1,cm,1,cm,1,cm,1");
				if (out_canopy && useCanopyModel) {
					fprintf(fout, ",kg m-2,degC,-,-,kg m-2,W m-2,W m-2,W m-2,W m-2,W m-2,kg m-2 per timestep,kg m-2 per timestep,kg m-2,kg m-2,kg m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,degC,kg m-2,kg m-2 per timestep");
				} else {
					fprintf(fout,",,,,,,,,,,,,,,,,,,,,,,,,,,,,");
				}
			} else if (out_t) {
				jj = 0;
				for (ii = MIN(5, number_fixed_heights); ii < number_fixed_heights+number_fixed_rates; ii++) {
					if (ii >= number_fixed_heights) {
						fprintf(fout, ",cm");
						jj++;
					}
					fprintf(fout, ",degC");
					jj++;
					if (ii < number_meas_temperatures) {
						fprintf(fout, ",degC");
						jj++;
					}
				}
				if (va_Xdata->tag_low) {
					unsigned int tag = va_Xdata->tag_low, j_lim;
					while ((tag + ii) <= number_sensors) {
						if ((tag + ii) <= number_meas_temperatures)
							j_lim = 41;
						else
							j_lim = 43;
						if (jj < j_lim) {
							fprintf(fout, ",cm,degC");
							jj += 2;
							if (ii < number_meas_temperatures) {
								fprintf(fout, ",cm,degC");
								jj += 2;
							}
							tag++;
						}
					}
				}
				for (; jj < 44; jj++)
					fprintf(fout,",");
			} else {
				fprintf(fout, ",,,,,,,,,,,,,,,,,,,,,,,,,,,,");

			}
			if (variant == "ANTARCTICA") {
				if (max_number_sensors == 5)
					fprintf(fout, ",");
				fprintf(fout, ",kJ m-2,kJ m-2,kg m-3,kg m-3,cm,%%,m s-1");
			} else if (variant == "CALIBRATION") {
				if (max_number_sensors == 5)
					fprintf(fout, ",");
				fprintf(fout, ",kJ m-2,kJ m-2,kg m-3,kg m-3,kg m-3,kg m-3,kg m-3");
			} else {
				fprintf(fout, ",kg m-2,kJ m-2,kJ m-2,kg m-3,kg m-3,cm");
				if (!research_mode)
					fprintf(fout, ",cm,kg m-2");
				else
					fprintf(fout, ",W m-2,W m-2");
			}

			fprintf(fout, "\n\n[DATA]");
		} else if ((strcmp(ext, "pro") == 0)) {
			va_Xdata = va_arg(argptr, SnowStation *);
			string stationname = va_Xdata->meta.getStationName();
			fprintf(fout, "[STATION_PARAMETERS]");
			fprintf(fout, "\nStationName= %s",   stationname.c_str());
			fprintf(fout, "\nLatitude= %.2f",   va_Xdata->meta.position.getLat());
			fprintf(fout, "\nLongitude= %.2f",  va_Xdata->meta.position.getLon());
			fprintf(fout, "\nAltitude= %.0f",   va_Xdata->meta.position.getAltitude());
			fprintf(fout, "\nSlopeAngle= %.2f", va_Xdata->meta.getSlopeAngle());
			fprintf(fout, "\nSlopeAzi= %.2f",   va_Xdata->meta.getAzimuth());

			fprintf(fout, "\n\n[HEADER]");
			if (out_haz) // HACK To avoid troubles in A3D
				fprintf(fout, "\n#%s, Snowpack %s version %s run by \"%s\"",
					   Hdata.sn_computation_date, variant.c_str(), Hdata.sn_version, Hdata.sn_user);
			fprintf(fout, "\n0500,Date");
			fprintf(fout, "\n0501,nElems,height [> 0: top, < 0: bottom of elem.] (cm)");
			fprintf(fout, "\n0502,nElems,element density (kg m-3)");
			fprintf(fout, "\n0503,nElems,element temperature (degC)");
			fprintf(fout, "\n0506,nElems,liquid water content by volume (%%)");
			fprintf(fout, "\n0508,nElems,dendricity (1)");
			fprintf(fout, "\n0509,nElems,sphericity (1)");
			fprintf(fout, "\n0510,nElems,coordination number (1)");
			fprintf(fout, "\n0511,nElems,bond size (mm)");
			fprintf(fout, "\n0512,nElems,grain size (mm)");
			fprintf(fout, "\n0513,nElems,grain type (Swiss Code F1F2F3)");
			fprintf(fout, "\n0515,nElems,ice volume fraction (%%)");
			fprintf(fout, "\n0516,nElems,air volume fraction (%%)");
			fprintf(fout, "\n0517,nElems,stress in (kPa)");
			fprintf(fout, "\n0518,nElems,viscosity (GPa s)");
			fprintf(fout, "\n0519,nElems,soil volume fraction (%%)");
			fprintf(fout, "\n0520,nElems,temperature gradient (K m-1)");
			fprintf(fout, "\n0521,nElems,thermal conductivity (W K-1 m-1)");
			fprintf(fout, "\n0522,nElems,absorbed shortwave radiation (W m-2)");
			fprintf(fout, "\n0523,nElems,viscous deformation rate (1.e-6 s-1)");
			fprintf(fout, "\n0530,nElems,position (cm) and minimum stability indices:");
			fprintf(fout, "\n          profile type, stability class, z_Sdef, Sdef, z_Sn38, Sn38, z_Sk38, Sk38");
			fprintf(fout, "\n0531,nElems,deformation rate stability index Sdef");
			fprintf(fout, "\n0532,nElems,natural stability index Sn38");
			fprintf(fout, "\n0533,nElems,stability index Sk38");
			fprintf(fout, "\n0534,nElems,hand hardness either (N) or index steps (1)");
			fprintf(fout, "\n0535,nElems,inverse texture index ITI (Mg m-4)");
			fprintf(fout, "\n0601,nElems,snow shear strength (kPa)");
			fprintf(fout, "\n0602,nElems,grain size difference (mm)");
			fprintf(fout, "\n0603,nElems,hardness difference (1)");
			fprintf(fout, "\n0604,nElems,ssi");
			if (variant == "CALIBRATION") {
				fprintf(fout, "\n0701,nElems,SNOWPACK: total settling rate (%% h-1)");
				fprintf(fout, "\n0702,nElems,SNOWPACK: settling rate due to load (%% h-1)");
				fprintf(fout, "\n0703,nElems,SNOWPACK: settling rate due to metamorphism (sig0) (%% h-1)");
				fprintf(fout, "\n0704,nElems,SNOWPACK: ratio -Sig0 to load EMS[e].C (1)");
				fprintf(fout, "\n0705,nElems,SNOWPACK: bond to grain ratio (1)");
				fprintf(fout, "\n0891,nElems,SNTHERM: settling rate due to load (%% h-1)");
				fprintf(fout, "\n0892,nElems,SNTHERM: settling rate due to metamorphism (%% h-1)");
				fprintf(fout, "\n0893,nElems,SNTHERM: viscosity (GPa s)");
			}
			fprintf(fout, "\n\n[DATA]");
		} else {
			prn_msg(__FILE__, __LINE__, "wrn", Date(), "No header defined for files *.%s", ext);
		}
		va_end(argptr);
		fclose(fout);
	}

	return true;
}

bool AsciiIO::writeHazardData(const std::string& /*stationID*/, const std::vector<ProcessDat>& /*Hdata*/,
                              const std::vector<ProcessInd>& /*Hdata_ind*/, const int& /*num*/)
{
	/*
	fout.open(name.c_str());
	if (fout.fail()) throw FileAccessException(name, AT);

	try {
		// Print out the hoar hazard data info, contained in Zdata
		fout << "SurfaceHoarIndex" << endl;
		for(unsigned int e = 0; e < 48; e++) {
			if (e != 0) fout << " ";
			fout << Zdata.hoar24[e];
		}
		fout << endl;

		// Print out the drift hazard data info, contained in Zdata
		fout << "DriftIndex" << endl;
		for(unsigned int e = 0; e < 48; e++) {
			if (e != 0) fout << " ";
			fout << Zdata.drift24[e];
		}
		fout << endl;

		// Print out the 3 hour new snowfall hazard data info, contained in Zdata
		fout << "ThreeHourNewSnow" << endl;
		for(unsigned int e = 0; e < 144; e++) {
			if (e != 0) fout << " ";
			fout << Zdata.hns3[e];
		}
		fout << endl;

		// Print out the 24 hour new snowfall hazard data info, contained in Zdata
		fout << "TwentyFourHourNewSnow" << endl;
		for(unsigned int e = 0; e < 144; e++) {
			if (e != 0) fout << " ";
			fout << Zdata.hns24[e];
		}
		fout << endl;
	} catch (const exception&){
		cleanup();
		throw;
	}

	cleanup();
	*/
	return true;
}

