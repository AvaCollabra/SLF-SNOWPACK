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

#include <snowpack/AsciiIO.h>

using namespace std;
using namespace mio;

AsciiIO::AsciiIO(const mio::Config& i_cfg) : cfg(i_cfg)
{
	string tmp_variant = cfg.get("VARIANT", "Parameters"); 
	variant=tmp_variant;
	string tmp_experiment = cfg.get("EXPERIMENT", "Parameters");
	experiment = tmp_experiment;
	string tmp_outpath = cfg.get("OUTPATH", "Parameters");
	outpath = tmp_outpath;

	i_snopath = "input";
	string tmp_snopath = cfg.get("SNOWPATH", "INPUT", Config::nothrow);
	if (tmp_snopath != "") i_snopath = tmp_snopath;

	o_snopath = outpath;
	string tmp_snopath_out = cfg.get("SNOWPATH", "OUTPUT", Config::nothrow);
	if (tmp_snopath_out != "") o_snopath = tmp_snopath_out;

	NUMBER_EXPO = cfg.get("NUMBER_EXPO", "Parameters");
	OUT_HEAT = cfg.get("OUT_HEAT", "Parameters");
	OUT_LW = cfg.get("OUT_LW", "Parameters");
	OUT_SW = cfg.get("OUT_SW", "Parameters");
	OUT_METEO = cfg.get("OUT_METEO", "Parameters");
	OUT_HAZ = cfg.get("OUT_HAZ", "Parameters");
	OUT_MASS = cfg.get("OUT_MASS", "Parameters");
	OUT_T = cfg.get("OUT_T", "Parameters");
	OUT_LOAD = cfg.get("OUT_LOAD", "Parameters");
	OUT_STAB = cfg.get("OUT_STAB", "Parameters");
	OUT_CANOPY = cfg.get("OUT_CANOPY", "Parameters");
	PERP_TO_SLOPE = cfg.get("PERP_TO_SLOPE", "Parameters");	

	PERP_TO_SLOPE = cfg.get("PERP_TO_SLOPE", "Parameters");	

	sw_ref = cfg.get("SW_REF", "Parameters");
	hazard_steps_between = cfg.get("HAZARD_STEPS_BETWEEN", "Parameters");
	calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Parameters");
	useCanopyModel = cfg.get("CANOPY", "Parameters");
	useSnowLayers = cfg.get("SNP_SOIL", "Parameters");
	t_internal = cfg.get("T_INTERNAL", "Parameters");
	depth_of_sensors = cfg.get("DEPTH", "Parameters");	

	research_mode = cfg.get("RESEARCH", "Parameters");
	
	if (depth_of_sensors.size() != FIXED_HEIGHTS)
		throw InvalidArgumentException("FIXED_HEIGHTS and the number of values for key DEPTH must match", AT);
}

/**
 * @brief This routine reads the status of the snow cover at program start
 * @author Perry Bartelt \n Michael Lehning
 * @version 9Y.mm
 * @param *SSdata
 * @param *Zdata
 * @param *XFilename Input file containing info about snow-cover structure (layers)
 */
void AsciiIO::readSnowCover(const std::string& station, SN_SNOWSOIL_DATA& SSdata, SN_ZWISCHEN_DATA& Zdata)
{
	string filename = getFilenamePrefix(station, i_snopath) + ".sno";

	if (research_mode){ //research mode: the key SNOWFILE/section INPUT is a valid identifier for the sno file
		string tmpsno = cfg.get("SNOWFILE", "INPUT", Config::nothrow);
		if (tmpsno != "") filename = tmpsno;
	}

	FILE *fin=NULL;
	int l,i; // Counters
	int YYYY,MM,DD,HH,MI,dum;
	
	fin = fopen(filename.c_str(), "r");
	if (fin == NULL) {
		prn_msg(__FILE__, __LINE__, "err", -1., "Cannot open profile INPUT file: %s", filename.c_str());
		throw IOException("Cannot generate Xdata from file " + filename, AT);
	}

	// Header, Station Name and Julian Date
	fscanf(fin, " %*s");
	fscanf(fin, "\nStationName= %*s");
	fscanf(fin, "\nProfileDate= %4d %2d %2d %2d %2d", &YYYY, &MM, &DD, &HH, &MI);
	SSdata.date = Date(YYYY, MM, DD, HH, MI);

	// Last checked measured Snow Height used for data Control of next run
	fscanf(fin, "\nHS_Last=%lf", &SSdata.Hslast);
	// Latitude, Longitude, Altitude, Slope Angle, Slope Azimut
	fscanf(fin, "\nLatitude=%lf", &SSdata.Lat);
	fscanf(fin, "\nLongitude=%lf", &SSdata.Lon);
	fscanf(fin, "\nAltitude=%lf", &SSdata.Alt);
	fscanf(fin, "\nSlopeAngle=%lf", &SSdata.Angle);
	fscanf(fin, "\nSlopeAzi=%lf", &SSdata.Azi);
	if ( (sw_ref%10 == 2) && PERP_TO_SLOPE && !(SSdata.Angle <= 3.) ) {
		prn_msg(__FILE__, __LINE__, "wrn", -1., 
			   "You want to use measured albedo in a slope steeper than 3 deg  with PERP_TO_SLOPE set!");
	}
	// Number of Soil Layers
	if (fscanf(fin, "\nnSoilLayerData=%d", &SSdata.nLayers) != 1) {
		prn_msg(__FILE__, __LINE__, "err", -1., "Missing 'nSoilLayerData' in Snowfile (*.sno)");
		throw IOException("Cannot generate Xdata", AT);
	}
	// Check Consistency of Soil Data with Switch
	if ( useSnowLayers && (SSdata.nLayers < 1) ) {
		prn_msg(__FILE__, __LINE__, "err", -1., "Missing 'nSoilLayerData' in Snowfile (*.sno), but SNP_SOIL set");
		throw IOException("Cannot generate Xdata", AT);
	}
	// Use this to mask the Erosion Level in operational mode -- use negative sign
	if ( SSdata.nLayers < 0 ) {
		SSdata.ErosionLevel = -SSdata.nLayers;
		SSdata.nLayers = 0;
	} else {
		SSdata.ErosionLevel = 0;
	}
	// Number of Snow Layers
	if (fscanf(fin, "\nnSnowLayerData=%d", &dum) != 1) {
		prn_msg(__FILE__, __LINE__, "err", -1., "Missing 'nSnowLayerData' in Snowfile (*.sno)");
		throw IOException("Cannot generate Xdata", AT);
	}
	SSdata.nLayers += dum;
	// Ground Characteristics (June 2006: corresponds to former SOIL_ALBEDO and BARE_SOIL_z0)
	if ( fscanf(fin, "\nSoilAlbedo=%lf", &SSdata.SoilAlb) != 1 ) {
		prn_msg(__FILE__, __LINE__, "err", -1., "Missing 'SoilAlbedo' in Snowfile (*.sno)");
		throw IOException("Cannot generate Xdata", AT);
	}
	if ( fscanf(fin, "\nBareSoil_z0=%lf", &SSdata.BareSoil_z0) != 1 ) {
		//BUG HACK: we should prevent it from being 0!!
		prn_msg(__FILE__, __LINE__, "err", -1., "Missing 'BareSoil_z0' in Snowfile (*.sno)");
	}
	if ( SSdata.BareSoil_z0==0. ) {
		prn_msg(__FILE__, __LINE__, "wrn", -1., "'BareSoil_z0'=0 from Snowfile, setting it to 0.02");
		SSdata.BareSoil_z0=0.02;
	}
	if (SSdata.Hslast > 0.05) {
		SSdata.Albedo = 0.9;
	} else {
		SSdata.Albedo = SSdata.SoilAlb;
	}
	// Canopy Characteristics
	fscanf(fin, "\nCanopyHeight=%lf",&SSdata.Canopy_Height);
	fscanf(fin, "\nCanopyLeafAreaIndex=%lf",&SSdata.Canopy_LAI);
	fscanf(fin, "\nCanopyDirectThroughfall=%lf",&SSdata.Canopy_Direct_Throughfall);

	// Allocate Space for the Layers
	if (SSdata.nLayers > 0) {
		if ( (SSdata.Ldata = (SN_LAYER_DATA *) realloc(SSdata.Ldata, sizeof(SN_LAYER_DATA)*SSdata.nLayers)) == NULL ) {
			prn_msg(__FILE__, __LINE__, "err", -1., "Cannot ALLOCATE layer Data");
			throw IOException("Cannot generate Xdata", AT);
		}
	}
	// Read the layer data
	if ( fscanf(fin,"\nYYYY") < 0 ) {
		prn_msg(__FILE__, __LINE__, "err", -1., "Failed reading layer header in Snowfile (*.sno)");
		throw IOException("Cannot generate Xdata", AT);
	}
	fscanf(fin, "%*[^\n]");
	for (l = 0; l < SSdata.nLayers; l++) {
		fscanf(fin, " %d %d %d %d %d", &YYYY, &MM, &DD, &HH, &MI);
		SSdata.Ldata[l].date = Date(YYYY, MM, DD, HH, MI);

		if ( SSdata.Ldata[l].date > SSdata.date ) {
			prn_msg(__FILE__, __LINE__, "err", -1., "Layer %d from bottom is younger (%lf) than ProfileDate (%lf) in Snowfile (*.sno)", l+1, SSdata.Ldata[l].date.getJulianDate(), SSdata.date.getJulianDate());
			throw IOException("Cannot generate Xdata", AT);
		}
		fscanf(fin, " %lf %lf %lf %lf %lf %lf", &SSdata.Ldata[l].hl, &SSdata.Ldata[l].tl, &SSdata.Ldata[l].phiIce, &SSdata.Ldata[l].phiWater, &SSdata.Ldata[l].phiVoids, &SSdata.Ldata[l].phiSoil);
		if (SSdata.Ldata[l].tl < 100.) {
			SSdata.Ldata[l].tl = C_TO_K(SSdata.Ldata[l].tl);
		}
		fscanf(fin, "%lf %lf %lf", &SSdata.Ldata[l].SoilRho, &SSdata.Ldata[l].SoilK, &SSdata.Ldata[l].SoilC);
		fscanf(fin, "%lf %lf %lf %lf %d %lf %d", &SSdata.Ldata[l].rg, &SSdata.Ldata[l].rb, &SSdata.Ldata[l].dd, &SSdata.Ldata[l].sp, &SSdata.Ldata[l].mk, &SSdata.Ldata[l].hr, &SSdata.Ldata[l].ne);
		if ( SSdata.Ldata[l].rg>0. && SSdata.Ldata[l].rb >= SSdata.Ldata[l].rg ) {
			//HACK To avoid surprises in lwsn_ConcaveNeckRadius()
			SSdata.Ldata[l].rb = Metamorphism::max_grain_bond_ratio * SSdata.Ldata[l].rg;
			prn_msg(__FILE__, __LINE__, "wrn", -1., "Layer %d from bottom: bond radius rb/rg larger than Metamorphism::max_grain_bond_ratio=%lf (rb=%lf mm, rg=%lf mm)! Reset to Metamorphism::max_grain_bond_ratio", l+1, Metamorphism::max_grain_bond_ratio, SSdata.Ldata[l].rb, SSdata.Ldata[l].rg);
		}
		for (i = 0; i < N_SOLUTES; i++) {
			fscanf(fin," %lf %lf %lf %lf ", &SSdata.Ldata[l].cIce[i], &SSdata.Ldata[l].cWater[i], &SSdata.Ldata[l].cVoids[i], &SSdata.Ldata[l].cSoil[i]);
		}
	}

	// READ THE Zdata
	// Read the hoar hazard data info, contained in Zdata
	fscanf(fin,"%*s ");
	for (i = 0; i < 48; i++) {
		if (fscanf(fin," %lf ", &Zdata.hoar24[i]) != 1) {
			prn_msg(__FILE__, __LINE__, "err", -1., "Reading Zdata (hoar)");
			throw IOException("Cannot generate Xdata", AT);
		}
	}
	// Read the drift hazard data info, contained in Zdata
	fscanf(fin,"%*s");
	for (i = 0; i < 48; i++) {
		if (fscanf(fin," %lf ", &Zdata.drift24[i]) != 1) {
			prn_msg(__FILE__, __LINE__, "err", -1., "Reading Zdata (drift)");
			throw IOException("Cannot generate Xdata", AT);
		}
	}
	// Read the 3 hour new snowfall hazard data info, contained in Zdata
	fscanf(fin,"%*s");
	for (i = 0; i < 144; i++) {
		if (fscanf(fin," %lf ", &Zdata.hns3[i]) != 1) {
			prn_msg(__FILE__, __LINE__, "err", -1., "Reading Zdata (hns3)");
			throw IOException("Cannot generate Xdata", AT);
		}
	}
	// Read the 24 hour new snowfall hazard data info, contained in Zdata
	fscanf(fin,"%*s");
	for (i = 0; i < 144; i++) {
		if (fscanf(fin," %lf ", &Zdata.hns24[i]) != 1) {
			prn_msg(__FILE__, __LINE__, "err", -1., "Reading Zdata (hns24)");
			throw IOException("Cannot generate Xdata", AT);
		}
	}

	/*
	 * Loop over the number of layers to calculate the total height and total number of finite elements in the profile.
	 * NOTE: the number of nodes is one PLUS the number of elements.
	*/
	SSdata.nN = 1;
	SSdata.Height = 0.;
	for (l = 0; l < SSdata.nLayers;l++) {
		SSdata.nN    += SSdata.Ldata[l].ne;
		SSdata.Height += SSdata.Ldata[l].hl;
	}

	fclose(fin);
} // End readXData

/**
 * @brief This routine writes the status of the snow cover at program termination and at specified backup times
 * @author Perry Bartelt \n Michael Lehning
 * @version 9Y.mm
 * @param JulianDate double
 * @param *station char
 * @param *Xdata SN_STATION_DATA
 * @param *Zdata SN_ZWISCHEN_DATA
 * @param *XFilename char
 */
void AsciiIO::writeSnowCover(const mio::Date& date, const std::string& station, const SN_STATION_DATA& Xdata, 
					    const SN_ZWISCHEN_DATA& Zdata, const bool& forbackup)
{
	FILE *fout=NULL;
	int  e, i;

	string filename = getFilenamePrefix(station, o_snopath) + ".sno";

	if (forbackup){
		stringstream ss;
		ss << (int)(date.getJulianDate() + 0.5);
		filename += ss.str();
	}

	const vector<SN_ELEM_DATA>& EMS = Xdata.Edata;
	fout = fopen( filename.c_str(), "w");
	if (fout == NULL) {
		prn_msg(__FILE__, __LINE__, "err", date.getJulianDate(),"Cannot open profile OUTPUT file: %s", filename.c_str());
		throw FileAccessException("Cannot dump final Xdata to file " + filename, AT);
	}

	// Header, Station Name and Julian Day
	fprintf(fout, "[SNOWPACK_INITIALIZATION]");
	fprintf(fout, "\nStationName= %s", station.c_str());
	
	int yyyy,mm,dd,hh,mi;
	date.getDate(yyyy,mm,dd,hh,mi);

	fprintf(fout, "\nProfileDate= %04d %02d %02d %02d %02d", yyyy, mm, dd, hh, mi);

	// Last checked Snow Depth used for data Control of next run
	fprintf(fout, "\nHS_Last= %lf", Xdata.cH - Xdata.Ground);

	// Latitude, Longitude, Altitude, Slope Angle, Slope Azimut
	fprintf(fout, "\nLatitude= %.4lf",   Xdata.Lat);
	fprintf(fout, "\nLongitude= %.4lf",  Xdata.Lon);
	fprintf(fout, "\nAltitude= %.0lf",   Xdata.Alt);
	fprintf(fout, "\nSlopeAngle= %.2lf", RAD_TO_DEG(Xdata.SlopeAngle));
	fprintf(fout, "\nSlopeAzi= %.2lf",   RAD_TO_DEG(Xdata.SlopeAzi));

	// Number of Soil Layer Data; in case of no soil used to store the erosion level
	if ( Xdata.SoilNode > 0 ) {
		fprintf(fout, "\nnSoilLayerData= %d", Xdata.SoilNode);
	} else {
		fprintf(fout, "\nnSoilLayerData= %d", -Xdata.ErosionLevel);
	}

	// Number of Snow Layer Data
	fprintf(fout, "\nnSnowLayerData= %d", Xdata.getNumberOfElements() - Xdata.SoilNode);

	// Ground Characteristics (introduced June 2006)
	fprintf(fout, "\nSoilAlbedo= %.2lf", Xdata.SoilAlb);
	fprintf(fout, "\nBareSoil_z0= %.3lf", Xdata.BareSoil_z0);

	// Canopy Characteristics
	fprintf(fout, "\nCanopyHeight= %.2lf",Xdata.Cdata.height);
	fprintf(fout, "\nCanopyLeafAreaIndex= %.6lf",Xdata.Cdata.lai);
	fprintf(fout, "\nCanopyDirectThroughfall= %.2lf",Xdata.Cdata.direct_throughfall);

	// Layer Data
	fprintf(fout, "\nYYYY MM DD HH MI Layer_Thick  T  Vol_Frac_I  Vol_Frac_W  Vol_Frac_V ");
	fprintf(fout, " Vol_Frac_S Rho_S Conduc_S HeatCapac_S  rg  rb  dd  sp  mk mass_hoar ne");
	for (i = 0; i < N_SOLUTES; i++) {
		fprintf(fout, " cIce cWater cAir  cSoil");
	}
	for (e = 0; e < Xdata.getNumberOfElements(); e++) {
		int YYYY, MM, DD, hh, mm;
		EMS[e].date.getDate(YYYY, MM, DD, hh, mm);
		fprintf(fout, "\n%04d %02d %02d %02d %02d", YYYY, MM, DD, hh, mm); 
		fprintf(fout, " %7.5lf %8.4lf %7.5lf %7.5lf %7.5lf",  EMS[e].L,
			Xdata.Ndata[e+1].T, EMS[e].theta[ICE], EMS[e].theta[WATER], EMS[e].theta[AIR]);
		fprintf(fout," %7.4lf %6.1lf %4.1lf %6.1lf %5.2lf %5.2lf %5.2lf %5.2lf %4d %7.5lf 1",
			EMS[e].theta[SOIL], EMS[e].soil[SOIL_RHO], EMS[e].soil[SOIL_K], EMS[e].soil[SOIL_C],
			EMS[e].rg, EMS[e].rb, EMS[e].dd,  EMS[e].sp, EMS[e].mk, Xdata.Ndata[e+1].hoar);

		for (i = 0; i < N_SOLUTES; i++) {
			fprintf(fout, "  %lf %9.7lf %9.7lf %9.7lf", EMS[e].conc[ICE][i], EMS[e].conc[WATER][i],
					EMS[e].conc[AIR][i], EMS[e].conc[SOIL][i]);
		}
	}

	// Print out the hoar hazard data info, contained in Zdata
	fprintf(fout,"\nSurfaceHoarIndex\n");
	for(e = 0; e < 48; e++) {
		fprintf(fout," %lf ", Zdata.hoar24[e]);
	}

	// Print out the drift hazard data info, contained in Zdata
	fprintf(fout,"\nDriftIndex\n");
	for (e = 0; e < 48; e++) {
		fprintf(fout," %lf ", Zdata.drift24[e]);
	}

	// Print out the 3 hour new snowfall hazard data info, contained in Zdata
	fprintf(fout,"\nThreeHourNewSnow\n");
	for (e = 0; e < 144; e++) {
		fprintf(fout," %lf ", Zdata.hns3[e]);
	}

	// Print out the 24 hour new snowfall hazard data info, contained in Zdata
	fprintf(fout,"\nTwentyFourHourNewSnow\n");
	for (e = 0; e < 144; e++) {
		fprintf(fout," %lf ", Zdata.hns24[e]);
	}

	// End of Data
	fprintf(fout, "\nEnd");
	fclose(fout);
}

std::string AsciiIO::getFilenamePrefix(const std::string& stationname, const std::string& path)
{
	//TODO: read only once (in constructor)
	string filename_prefix = path + "/" + stationname;

	if ((experiment != "none") && (experiment != "opera")) //in operational mode, nothing is appended
		filename_prefix += "_" + experiment; // complete filename_prefix

	return filename_prefix;
}

/**
 * @brief Write the Snow Profile Results, snow depth being taken VERTICALLY
 * Prepare Output File for JAVA visualization (SNOWPACK format, *.pro)
 * NOTE Parameters marked by an asterisk are available in RESEARCH visualisation only!
 * @author Michael Lehning
 * @version 9Y.mm
 * @param filename string
 * @param i_date the date
 * @param Xdata SN_STATION_DATA
 * @param Hdata Q_PROCESS_DAT
 */
void AsciiIO::writeProfile(const mio::Date& i_date, const std::string& stationname, const unsigned int& expo,
					  const SN_STATION_DATA& Xdata, const Q_PROCESS_DAT& Hdata)
{
	stringstream ss;
	ss << stationname;
	if (expo != 0) ss << expo;

	const string station = ss.str();
	string filename = getFilenamePrefix(station, outpath) + ".pro";

	FILE *PFile=NULL;
	int e, nN, nE, nz;
	double cos_sl;

	const vector<SN_ELEM_DATA>& EMS = Xdata.Edata; nE = Xdata.getNumberOfElements();
	const vector<SN_NODE_DATA>& NDS = Xdata.Ndata; nN = Xdata.getNumberOfNodes();

	if ( !checkHeader(filename.c_str(), "[STATION_PARAMETERS]", &Hdata, station, "pro", &Xdata) ) {
		prn_msg(__FILE__, __LINE__, "err", i_date.getJulianDate(),"Checking header in file %s", filename.c_str());
		throw IOException("Cannot dump profile " + filename + " for Java Visualisation", AT);
	} else if ( !(PFile = fopen(filename.c_str(), "a")) ) {
		prn_msg(__FILE__, __LINE__, "err", i_date.getJulianDate(),
			   "Cannot open profile series file: %s", filename.c_str());
		throw IOException("Cannot dum profile " + filename + "for Java Visualisation", AT);
	}
	
	fprintf(PFile,"\n0500,%s", i_date.toString(Date::DIN).c_str());
	cos_sl = cos(Xdata.SlopeAngle);

	if (useSnowLayers){
		nz = nN;
	} else {
		nz = nE;
	}
	//  501: height [> 0: top, < 0: bottom of elem.] (cm)
	fprintf(PFile,"\n0501,%d", nz);
	if ( nz < 1 ) {
		// no soil and no snow
		fclose(PFile);
		return;
	}
	for (e = nN-nz; e < nN; e++) {
		fprintf(PFile,",%.2lf",M_TO_CM((NDS[e].z+NDS[e].u - NDS[Xdata.SoilNode].z)/cos_sl));
	}
	//  502: element density (kg m-3)
	fprintf(PFile,"\n0502,%d", nE);
	for (e = 0; e < nE; e++) {
		fprintf(PFile,",%.1lf",EMS[e].Rho);
	}
	//  503: element temperature (degC)
	fprintf(PFile,"\n0503,%d", nE);
	for (e = 0; e < nE; e++) {
		fprintf(PFile,",%.2lf",K_TO_C(EMS[e].Te));
	}
	//  506: liquid water content by volume (%)
	fprintf(PFile,"\n0506,%d", nE);
	for (e = 0; e < nE; e++) {
		fprintf(PFile,",%.1lf",100.*EMS[e].theta[WATER]);
	}
	// *508: dendricity (1)
	fprintf(PFile,"\n0508,%d", nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++) {
		fprintf(PFile,",%.2lf",EMS[e].dd);
	}
	// *509: sphericity (1)
	fprintf(PFile,"\n0509,%d", nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++) {
		fprintf(PFile,",%.2lf",EMS[e].sp);
	}
	// *510: coordination number (1)
	fprintf(PFile,"\n0510,%d", nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++) {
		fprintf(PFile,",%.1lf",EMS[e].N3);
	}
	// *511: bond size (mm)
	fprintf(PFile,"\n0511,%d", nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++) {
		fprintf(PFile,",%.2lf",2.*EMS[e].rb);
	}
	//  512: grain size (mm)
	fprintf(PFile,"\n0512,%d", nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++) {
		fprintf(PFile,",%.2lf",2.*EMS[e].rg);
	}
	//  513: grain type (Swiss code F1F2F3)
	fprintf(PFile,"\n0513,%d", nE+1-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++) {
		fprintf(PFile,",%03d",EMS[e].type);
	}
	// surface hoar at surface? (depending on boundary conditions)
	if ( M_TO_MM(NDS[nN-1].hoar/DENSITY_HOAR_SURF) > MIN_SIZE_HOAR_SURF ) {
		fprintf(PFile,",660");
	} else {
		fprintf(PFile,",0");
	}
	// *515: ice volume fraction (%)
	fprintf(PFile,"\n0515,%d", nE);
	for (e = 0; e < nE; e++) {
		fprintf(PFile,",%.0lf",100.*EMS[e].theta[ICE]);
	}
	// *516: air volume fraction (%)
	fprintf(PFile,"\n0516,%d", nE);
	for (e = 0; e < nE; e++) {
		fprintf(PFile,",%.0lf",100.*EMS[e].theta[AIR]);
	}
	// *517: stress (kPa)
	fprintf(PFile,"\n0517,%d", nE);
	for (e = 0; e < nE; e++) {
		fprintf(PFile,",%.3e",1.e-3*EMS[e].C);
	}
	// *518: viscosity (GPa s)
	fprintf(PFile,"\n0518,%d", nE);
	for (e = 0; e < nE; e++) {
		fprintf(PFile,",%.3e",1.e-9*EMS[e].k[SETTLEMENT]);
	}
	// *519: soil volume fraction (%)
	fprintf(PFile,"\n0519,%d", nE);
	for (e = 0; e < nE; e++) {
		fprintf(PFile,",%.0lf",100.*EMS[e].theta[SOIL]);
	}
	// *520: temperature gradient (K m-1)
	fprintf(PFile,"\n0520,%d", nE);
	for (e = 0; e < nE; e++) {
		fprintf(PFile,",%.3e",EMS[e].gradT);
	}
	// *521: thermal conductivity (W K-1 m-1)
	fprintf(PFile,"\n0521,%d", nE);
	for (e = 0; e < nE; e++) {
		fprintf(PFile,",%.3e",EMS[e].k[TEMPERATURE]);
	}
	// *522: absorbed shortwave radiation (W m-2)
	fprintf(PFile,"\n0522,%d", nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++) {
		fprintf(PFile,",%.1lf",EMS[e].sw_abs);
	}
	// *523: viscous deformation rate (1.e-6 s-1)
	fprintf(PFile,"\n0523,%d", nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++) {
		fprintf(PFile,",%.1lf",1.e6*EMS[e].EvDot);
	}
	//  530: position (cm) and minimum stability indices
	fprintf(PFile,"\n0530,%d", 8);
	fprintf(PFile,",%d,%d,%.1lf,%.2lf,%.1lf,%.2lf,%.1lf,%.2lf", Xdata.S_class1, Xdata.S_class2, M_TO_CM(Xdata.z_S_d/cos_sl), Xdata.S_d, M_TO_CM(Xdata.z_S_n/cos_sl), Xdata.S_n, M_TO_CM(Xdata.z_S_s/cos_sl), Xdata.S_s);
	//  531: deformation rate stability index Sdef
	fprintf(PFile,"\n0531,%d" ,nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++) {
		fprintf(PFile,",%.1lf",EMS[e].S_dr);
	}
	//  532: natural stability index Sn38
	fprintf(PFile,"\n0532,%d" ,nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode;  e < nE; e++) {
		fprintf(PFile,",%.1lf",NDS[e+1].S_n);
	}
	//  533: stability index Sk38
	fprintf(PFile,"\n0533,%d" ,nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++) {
		fprintf(PFile,",%.1lf",NDS[e+1].S_s);
	}
	//  534: hand hardness ...
	fprintf(PFile,"\n0534,%d" ,nE-Xdata.SoilNode);
	if ( R_IN_N ) { // ... either converted to newtons according to Swiss scale
		for (e = Xdata.SoilNode; e < nE; e++) {
			fprintf(PFile,",%.1lf",-1.*(19.472*pow(EMS[e].hard, 2.3607)));
		}
	} else { // ... or in index steps (1)
		for (e = Xdata.SoilNode; e < nE; e++) {
			fprintf(PFile,",%.1lf", -EMS[e].hard);
		}
	}
	//  535: inverse texture index ITI (Mg m-4)
	fprintf(PFile,"\n0535,%d" ,nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++) {
	 	if ( EMS[e].dd < 0.005 ) {
			fprintf(PFile,",%.1lf",-1.*EMS[e].Rho/(2.*MM_TO_M(EMS[e].rg)));
		} else {
			fprintf(PFile,",%.1lf",0.0);
		}
	}

	if (variant == "CALIBRATION"){
		writeFreeProfilesCALIBRATION(Xdata, PFile);
	} else {
		writeFreeProfilesDEFAULT(Xdata, PFile);
	}

	fclose(PFile);
}

/**
 * @brief Default: dump special profiles to *.pro output file
 * @author Charles Fierz
 * @version 10.04
 * @param *Xdata
 * @param *fout Output file
 */
void AsciiIO::writeFreeProfilesDEFAULT(const SN_STATION_DATA& Xdata, FILE *fout)
{
	int e;
	const int nE = Xdata.getNumberOfElements();

	const vector<SN_ELEM_DATA>& EMS = Xdata.Edata;
	if ( OUT_LOAD ) {
		// *6nn: e.g. solute concentration
		int i, j;
		for (j = 2; j < N_COMPONENTS-1; j++) {
			for (i = 0; i < N_SOLUTES; i++) {
				fprintf(fout,"\n06%02d,%d" , 10*j + i,nE-Xdata.SoilNode);
				for (e = Xdata.SoilNode; e < nE; e++) {
					fprintf(fout,",%.1lf",EMS[e].conc[i][j]);
				}
			}
		}
	} else {
		// 600-profile specials
		// *601: snow shear strength (kPa)
		fprintf(fout,"\n0601,%d" ,nE-Xdata.SoilNode);
		for (e = Xdata.SoilNode; e < nE; e++) {
			fprintf(fout,",%.2lf",EMS[e].s_strength);
		}
		// *602: grain size difference (mm)
		fprintf(fout,"\n0602,%d" ,nE-Xdata.SoilNode);
		for (e = Xdata.SoilNode; e < nE-1; e++) {
			fprintf(fout,",%.2lf",2.*fabs(EMS[e].rg - EMS[e+1].rg));
		}
		fprintf(fout,",0.");
		// *603: hardness difference (1)
		fprintf(fout,"\n0603,%d" ,nE-Xdata.SoilNode);
		for (e = Xdata.SoilNode; e < nE-1; e++) {
			fprintf(fout,",%.2lf",fabs(EMS[e].hard - EMS[e+1].hard));
		}
		fprintf(fout,",0.");
	}
}

/**
 * @brief Calibration: dump special profiles to *.pro output file
 * @author Charles Fierz
 * @version 10.04
 * @param *Xdata
 * @param *fout Output file
 */
void AsciiIO::writeFreeProfilesCALIBRATION(const SN_STATION_DATA& Xdata, FILE *fout)
{
	int e;
	const int nE = Xdata.getNumberOfElements();

	const vector<SN_ELEM_DATA>& EMS = Xdata.Edata;
	const vector<SN_NODE_DATA>& NDS = Xdata.Ndata;
	// 600-profile specials
	// *601: snow shear strength (kPa)
	fprintf(fout,"\n0601,%d",nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE; e++) {
		fprintf(fout,",%.2lf",EMS[e].s_strength);
	}
	// *602: grain size difference (mm)
	fprintf(fout,"\n0602,%d",nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE-1; e++) {
		fprintf(fout,",%.2lf",2.*fabs(EMS[e].rg - EMS[e+1].rg));
	}
	fprintf(fout,",0.");
	// *603: hardness difference (1)
	fprintf(fout,"\n0603,%d",nE-Xdata.SoilNode);
	for (e = Xdata.SoilNode; e < nE-1; e++) {
		fprintf(fout,",%.2lf",fabs(EMS[e].hard - EMS[e+1].hard));
	}
	fprintf(fout,",0.");

	// 700-profile specials for settling comparison
	// *701: SNOWPACK: total settling rate (% h-1)
	fprintf(fout,"\n0701,%d",nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++) {
		fprintf(fout,",%.2lf", -100.*H_TO_S(EMS[e].EvDot));
	}
	// *702: SNOWPACK: settling rate due to load (% h-1)
	fprintf(fout,"\n0702,%d",nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++) {
		fprintf(fout, ",%.2lf", -100.*H_TO_S(NDS[e].udot));
	}
	// *703: SNOWPACK: settling rate due to metamorphism (sig0) (% h-1)
	fprintf(fout,"\n0703,%d",nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++) {
		fprintf(fout, ",%.2lf", -100.*H_TO_S(NDS[e].f));
	}
	// *704: SNOWPACK: ratio -Sig0 to load EMS[e].C (1)
	fprintf(fout,"\n0704,%d",nE-Xdata.SoilNode);
	for(e=Xdata.SoilNode; e<nE; e++) {
		fprintf(fout,",%.2lf", EMS[e].EDot);
	}
	// *705: SNOWPACK: bond to grain ratio (1)
	fprintf(fout,"\n0705,%d",nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++) {
		fprintf(fout,",%.4lf", EMS[e].rb / EMS[e].rg);
	}
	// SNTHERM.89
	// *891: SNTHERM: settling rate due to load (% h-1)
	fprintf(fout,"\n0891,%d" ,nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++) {
		const double eta_sntherm = (3.6e6*exp(0.08*(273.15-EMS[e].Te))*exp(0.021*EMS[e].Rho));
		fprintf(fout,",%.2lf", -100.*H_TO_S(EMS[e].C/eta_sntherm));
	}
	// *892: SNTHERM: settling rate due to metamorphism (% h-1)
	fprintf(fout,"\n0892,%d" ,nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++) {
		double evdot = -2.778e-6*exp(-0.04*(273.15 - EMS[e].Te));
		if( EMS[e].Rho > 150. ){
			evdot *= exp(-0.046*(EMS[e].Rho-150.));
		}
		if( EMS[e].theta[WATER] > 0.01 ){
			evdot *= 2.;
		}
		fprintf(fout, ",%.2lf", -100.*H_TO_S(evdot));
	}
	// *893: SNTHERM: viscosity (GPa s)
	fprintf(fout,"\n0893,%d" ,nE-Xdata.SoilNode);
	for (e=Xdata.SoilNode; e<nE; e++){
		const double eta_sntherm = (3.6e6*exp(0.08*(273.15-EMS[e].Te))*exp(0.021*EMS[e].Rho));
		fprintf(fout,",%.2lf", 1.e-9*eta_sntherm);
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
int AsciiIO::writeTemperatures(FILE *fout, const double& z_vert, const double& T, 
						 const int& i, const SN_STATION_DATA& Xdata)
{
	int j=2;
	double perp_pos, temp;

	//HACK:
	/// @brief Initial height of snow needed to compute sensor position from ground if FIXED_RATES is set
	double INITIAL_HS=0;

	if ( i < FIXED_HEIGHTS ) {
		perp_pos = calcPerpPosition(z_vert, Xdata.cH, Xdata.Ground, Xdata.SlopeAngle);
	} else {
		if ( (perp_pos = calcPerpPosition(z_vert, INITIAL_HS, Xdata.Ground, Xdata.SlopeAngle)) == NODATA ) {
			fprintf(fout, ",");
		} else {
			fprintf(fout, ",%.2lf", M_TO_CM(perp_pos)/cos(Xdata.SlopeAngle));
		}
		j++;
	}
	temp = getModelledTemperature(perp_pos, Xdata);
	fprintf(fout, ",%.2lf", temp);
	if ( i < t_internal ) {
		temp = checkMeasuredTemperature(T, perp_pos, Xdata.mH);
		fprintf(fout,",%.2lf", temp);
	} else {
		fprintf(fout, ",");
	}
	return j;
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
 * @param SlopeAngle (rad)
 */
double AsciiIO::calcPerpPosition(const double& z_vert, const double& hs_ref, const double& ground, const double& angleSlope)
{
	if ( z_vert == NODATA ) {
		return NODATA;
	} else if ( !useSnowLayers && (z_vert < 0.) ) {
		return (MAX(ground, hs_ref + z_vert * cos(angleSlope)));
	} else {
		return (ground + z_vert * cos(angleSlope));
	}
}

/**
 * @brief Checks whether measured internal snow or/and soil temperature (instantaneous value) is valid \n
 * The temperature defaults to NODATA if
 *  - the sensor is not covered by more than MIN_DEPTH_SUBSURF snow (measured perpendicular to slope)
 * @author Charles Fierz
 * @version 10.01
 * @param T Measured temperature (K)
 * @param z Sensor position perpendicular to slope (m)
 * @param mH Measured snow height (m)
 * @return Measured temperature (degC) if OK, NODATA else
 */
double AsciiIO::checkMeasuredTemperature(const double& T, const double& z, const double& mH)
{
	if ( (z <= (mH - MIN_DEPTH_SUBSURF)) && (T != NODATA) ) {
		return K_TO_C(T);
	} else {
		return NODATA;
	}
}

/**
 * @brief Find element with corresponding tag or return -1 if not found
 * @author Charles Fierz
 * @version 10.04
 * @param tag Tag to look for
 * @param *Xdata
 * @return Index of tagged element, -1 if not found
 */
int AsciiIO::findTaggedElement(const int& tag, const SN_STATION_DATA& Xdata)
{
	int e;
	for ( e=0; e<Xdata.getNumberOfElements(); e++) {
		if ( Xdata.Edata[e].mk/100 == tag ) {
			return e;
		}
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
int AsciiIO::writeHeightTemperatureTag(FILE *fout, const int& tag, const SN_MET_DATA& Mdata, const SN_STATION_DATA& Xdata)
{
	int e, j=2;
	const int i = FIXED_HEIGHTS+FIXED_RATES + (tag-1);
	double perp_pos, temp;
	if ( (e = findTaggedElement(tag, Xdata)) >= 0 ) {
		perp_pos = ((Xdata.Ndata[e].z + Xdata.Ndata[e].u + Xdata.Ndata[e+1].z + Xdata.Ndata[e+1].u)/2. - Xdata.Ground);
		fprintf(fout,",%.2lf,%.2lf", M_TO_CM(perp_pos)/cos(Xdata.SlopeAngle), K_TO_C(Xdata.Edata[e].Te));
	} else {
		fprintf(fout,",,%.2lf", NODATA);
	}
	if (i < t_internal) {
		if ( (perp_pos = calcPerpPosition(Mdata.zv_ts[i], Xdata.cH, Xdata.Ground, Xdata.SlopeAngle)) == NODATA ) {
			fprintf(fout,",,%.2lf", NODATA);
		} else {
			fprintf(fout,",%.2lf", M_TO_CM(perp_pos)/cos(Xdata.SlopeAngle));
			temp = checkMeasuredTemperature(Mdata.ts[i], perp_pos, Xdata.mH);
			fprintf(fout,",%.2lf", temp);
		}
		j += 2;
	}
	return j;
}

/**
 * @brief Write all Time Series results (*.met)
 * All depths and water equivalents (mass) are taken VERTICALLY. \n
 * If AVGSUM_TIME_SERIES is set, mean fluxes and cumulated masses since last dump are written, \n
 * else current energy fluxes, cumulated masses over last calculation_step_length (recommended setting in operational mode).
 * If CUMSUM_MASS is set, current value of cumulated masses since begin of run are dumped. \n
 * Precipitations (rain& snow, rain) are always dumped as rates (kg m-2 h-1). \n
 * NOTE:
 * 	-# neither AVGSUM_TIME_SERIES nor CUMSUM_MASS can be set if NUMBER_EXPO > 1.
 * 	-# When running SNOW_REDISTRIBUTION on virtual slopes, the correct drift index will
 *     be contained in the windward and lee *.met files only!
 * \li DO NOT change the order of parameters below! Additional parameters may be dumped at pos.
 *     93[94] to 100 in writeFreeSeriesXXX()
 * @author Michael Lehning \n Charles Fierz
 * @version 9Y.mm
 * @param *Xdata
 * @param *Sdata
 * @param *Mdata
 * @param *Hdata
 * @param *TFilename Name of file holding time series
 */
void AsciiIO::writeTimeSeries(const std::string& station, const SN_STATION_DATA& Xdata, 
						const SN_SURFACE_DATA& Sdata, const SN_MET_DATA& Mdata, const Q_PROCESS_DAT& Hdata)
{
	FILE *TFile=NULL;
	int    nE, nN; // Number of elements and nodes
	int    i, j;
	double cos_sl;

	string filename = getFilenamePrefix(station, outpath) + ".met";

	const vector<SN_NODE_DATA>& NDS = Xdata.Ndata;
	nN=Xdata.getNumberOfNodes();
	nE=Xdata.getNumberOfElements();
	cos_sl = cos(Xdata.SlopeAngle);

	// Check file for header
	if ( !checkHeader(filename.c_str(), "[STATION_PARAMETERS]", &Hdata, station, "met", &Xdata) ) {
		prn_msg(__FILE__, __LINE__, "err", Mdata.date.getJulianDate(), "Checking header in file %s", filename.c_str());
		throw InvalidFormatException("Writing Time Series data failed", AT);
	} else if ( !(TFile = fopen(filename.c_str(), "a")) ) {
		prn_msg(__FILE__, __LINE__, "err", Mdata.date.getJulianDate(), 
			   "Cannot open time series file: %s", filename.c_str());
	     throw FileAccessException(filename, AT);
	}
	// Print time stamp
	fprintf(TFile,"\n0203,%s", Mdata.date.toString(Date::DIN).c_str());
	if ( OUT_HEAT ) {
		// 1-2: Turbulent fluxes (W m-2)
		fprintf(TFile,",%lf,%lf", Sdata.qs, Sdata.ql);
	} else {
		fprintf(TFile,",,");
	}
	if ( OUT_LW ) {
		// 3-5: Longwave radiation fluxes (W m-2)
		fprintf(TFile,",%lf,%lf,%lf", Sdata.lw_out, Sdata.lw_in, Sdata.lw_net);
	} else {
		fprintf(TFile,",,,");
	}
	if ( OUT_SW ) {
		// 6-9: Shortwave radiation fluxes (W m-2) and calc. albedo (1)
		fprintf(TFile,",%lf,%lf,%lf,%lf", Sdata.sw_out, Sdata.sw_in, Sdata.qw, Sdata.cA);
	} else {
		fprintf(TFile,",,,,");
	}
	if ( OUT_METEO ) {
		// 10-13: Air temperature, snow surface temperature (modeled and measured) (degC)
		fprintf(TFile,",%lf,%lf,%lf,%lf", K_TO_C(Mdata.ta), K_TO_C(NDS[nE].T), K_TO_C(Mdata.tss), K_TO_C(NDS[0].T));
	} else {
		fprintf(TFile,",,,,");
	}
	if ( OUT_HEAT ) {
		// 14-17: soil heat fluxes (W m-2), ground surface temperature (degC), rain energy (W m-2)
		fprintf(TFile,",%lf,%lf,%lf,%lf", Sdata.qg, K_TO_C(NDS[Xdata.SoilNode].T), Sdata.qg0, Sdata.qr);
	} else {
		fprintf(TFile,",,,,");
	}
	if ( OUT_SW ) {
		// 18-22: projected solar radiation (W m-2), meas. albedo (1)
		fprintf(TFile,",%lf,%lf,%lf,%lf,%lf", Sdata.sw_hor, Sdata.sw_in, Sdata.sw_dir, Sdata.sw_diff, Sdata.mA);
	} else {
		fprintf(TFile,",,,,,");
	}
	if ( OUT_METEO ) {
		// 23-26: rH (%), wind (m s-1), wind_max (m s-1), wind_dir (deg),
		// 27: solid precipitation rate (kg m-2 h-1),
		// 28-29: modeled and maesured vertical snow depth (cm)
		fprintf(TFile,",%lf,%lf,%lf,%lf,%lf,%lf,%lf", 100.*Mdata.rh, Mdata.vw, Mdata.vw_max, Mdata.dw, Sdata.mass[SN_SURFACE_DATA::MS_PRECIP]/cos_sl, M_TO_CM((Xdata.cH - Xdata.Ground)/cos(Xdata.SlopeAngle)), M_TO_CM(Mdata.hs1/cos(Xdata.SlopeAngle)));
	} else {
		fprintf(TFile,",,,,,,,");
	}
	if ( OUT_HAZ ) {
		// 30-33: surface hoar size (mm), 24h drift index (cm), height of new snow HN (cm), 3d sum of daily new snow depths (cm)
		if ( !PERP_TO_SLOPE ) {
			fprintf(TFile,",%lf,%lf,%lf,%lf", Hdata.hoar_size, Hdata.wind_trans24, Hdata.hns24, Hdata.hns72_24);
		} else {
			// dump vertical values if PERP_TO_SLOPE
			fprintf(TFile,",%lf,%lf,%lf,%lf", Hdata.hoar_size, Hdata.wind_trans24, Hdata.hns24/cos_sl, Hdata.hns72_24/cos_sl);
		}
	} else {
		fprintf(TFile,",,,,");
	}
	if ( OUT_MASS ) {
		// 34-39: total mass, eroded mass, rain rate, runoff at bottom of snowpack, sublimation and evaporation, all in kg m-2 except rain as rate: kg m-2 h-1; see also 51-52 & 93
		fprintf(TFile,",%lf,%lf,%lf,%lf,%lf,%lf", Sdata.mass[SN_SURFACE_DATA::MS_TOTALMASS]/cos_sl, Sdata.mass[SN_SURFACE_DATA::MS_WIND]/cos_sl, Sdata.mass[SN_SURFACE_DATA::MS_RAIN]/cos_sl, Sdata.mass[SN_SURFACE_DATA::MS_RUNOFF]/cos_sl, Sdata.mass[SN_SURFACE_DATA::MS_SUBLIMATION]/cos_sl, Sdata.mass[SN_SURFACE_DATA::MS_EVAPORATION]/cos_sl);
	} else{
		fprintf(TFile,",,,,,,");
	}
	// 40-49: Internal Temperature Time Series at fixed heights, modeled and measured, all in degC
	if ( OUT_T && (FIXED_HEIGHTS || FIXED_RATES) ) {
		j = 0;
		for (i = 0; i < MIN(5, FIXED_HEIGHTS); i++) {
			j += writeTemperatures(TFile, Mdata.zv_ts[i], Mdata.ts[i], i, Xdata);
		}
    for (; j < 10; j++) {
      fprintf(TFile,",");
    }
	} else {
		fprintf(TFile,",,,,,,,,,,");
	}
	if ( MAX_NUMBER_SENSORS == 5 ) {
		if ( OUT_LOAD ) {
		// 50: Solute load at ground surface
			fprintf(TFile,",%lf",Sdata.load[0]);
		} else {
			fprintf(TFile,",");
		}
		if ( OUT_MASS ) {
		// 51-52: SWE (for checks) and LWC (kg m-2); see also 34-39
			fprintf(TFile,",%lf,%lf", Sdata.mass[SN_SURFACE_DATA::MS_SWE]/cos_sl, Sdata.mass[SN_SURFACE_DATA::MS_WATER]/cos_sl);
		} else {
			fprintf(TFile,",,");
		}
		if ( OUT_STAB ) {
			// 53-64: Stability Time Series, heights in cm
			fprintf(TFile,",%d,%d,%.1lf,%.2lf,%.1lf,%.2lf,%.1lf,%.2lf,%.1lf,%.2lf,%.1lf,%.2lf", Xdata.S_class1, Xdata.S_class2, M_TO_CM(Xdata.z_S_d/cos_sl), Xdata.S_d, M_TO_CM(Xdata.z_S_n/cos_sl), Xdata.S_n, M_TO_CM(Xdata.z_S_s/cos_sl), Xdata.S_s, M_TO_CM(Xdata.z_S_4/cos_sl), Xdata.S_4, M_TO_CM(Xdata.z_S_5/cos_sl), Xdata.S_5);
		} else {
			fprintf(TFile,",,,,,,,,,,,,");
		}
		if ( OUT_CANOPY && useCanopyModel ) {
			// 65-92 (28 columns)
			Canopy::cn_DumpCanopyData(TFile, &Xdata.Cdata, &Sdata, cos_sl);
		} else {
			fprintf(TFile,",,,,,,,,,,,,,,,,,,,,,,,,,,,,");
		}
	} else if ( OUT_T ) {
		// 50-93 (44 columns)
		j = 0;
		for (i = MIN(5, FIXED_HEIGHTS); i < FIXED_HEIGHTS+FIXED_RATES; i++) {
			if ( (j += writeTemperatures(TFile, Mdata.zv_ts[i], Mdata.ts[i], i, Xdata)) > 44 ) {
				prn_msg(__FILE__, __LINE__, "err", Mdata.date.getJulianDate(), 
					   "There is not enough space to accomodate your temperature sensors: j=%d > 44!", j);
				throw IOException("Writing Time Series data failed", AT);
      }
		}
		if ( Xdata.tag_low ) {
			int tag = Xdata.tag_low, j_lim;
			while ( (tag + i) <= NUMBER_SENSORS ) {
			  if ( (tag + i) <= t_internal ) {
			    j_lim = 41;
			  } else {
          j_lim = 43;
        }
        if ( j < j_lim ) {
					j += writeHeightTemperatureTag(TFile, tag, Mdata, Xdata);
          tag++;
        } else {
          break;
        }
			}
		}
    for (; j < 44; j++) {
      fprintf(TFile,",");
    }
	} else {
		fprintf(TFile,",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,");
	}
	// 93[94]-100 (8 or 7 free columns)

	if (variant == "CALIBRATION"){
		writeFreeSeriesCALIBRATION(&Xdata, &Sdata, &Mdata, &Hdata, TFile);
	} else if (variant == "ANTARCTICA"){
		writeFreeSeriesANTARCTICA(&Xdata, &Sdata, &Mdata, &Hdata, TFile);
	} else {
		writeFreeSeriesDEFAULT(&Xdata, &Sdata, &Mdata, &Hdata, TFile);
	}

	fclose (TFile);
}


/**
 * @brief Default: last 8 time series (columns 93 to 100) dumped to *.met output file
 * @author Charles Fierz
 * @version 10.02
 * @param *Xdata
 * @param *Sdata
 * @param *Mdata
 * @param *Hdata
 * @param *fout Output file
 */
void AsciiIO::writeFreeSeriesDEFAULT(const SN_STATION_DATA *Xdata, const SN_SURFACE_DATA *Sdata, 
							  const SN_MET_DATA *Mdata, const Q_PROCESS_DAT *Hdata, FILE *fout)
{
	(void) *Mdata;
	double cos_sl = cos(Xdata->SlopeAngle);
	// 93-100 (8 columns)
	if ( useSnowLayers ) {
		// 93: Soil Runoff (kg m-2); see also 34-39 & 51-52
		fprintf(fout,",%lf",Sdata->mass[SN_SURFACE_DATA::MS_SOIL_RUNOFF]/cos_sl);
	} else {
		fprintf(fout,",");
	}
	if ( OUT_HEAT ) {
		// 94-95: Measured Turbulent Fluxes (W m-2); see also 1-2
		fprintf(fout,",%lf,%lf",0.0,0.0);
	} else {
		fprintf(fout,",,");
	}
	if ( OUT_HAZ ) {
		// 96: crust thickness (S-slope) (cm)
		fprintf(fout,",%lf", Hdata->crust); //lwsn_SnowpackEnergyInput(Xdata)); //
		// 97:
		if ( !research_mode ) { // energy input to surface (kJ m-2)
			fprintf(fout,",%lf", Hdata->en_bal);
		} else { // change of internal energy (kJ m-2)
			fprintf(fout,",%lf", Sdata->dIntEnergy / 1000.);
		}
	} else {
		fprintf(fout,",,");
	}
	// 98: sum of energy fluxes at surface (kJ m-2)
	if ( !research_mode ) {
		fprintf(fout,",%lf", ((Sdata->qw + Sdata->lw_net + Sdata->qs + Sdata->ql + Sdata->qr) * M_TO_S(calculation_step_length) * hazard_steps_between) / 1000.);
	} else {
		fprintf(fout,",%lf", ((Sdata->qw + Sdata->lw_net + Sdata->qs + Sdata->ql + Sdata->qr) * M_TO_S(calculation_step_length)) / 1000.);
	}
	if ( !research_mode ) {
		// 99-100: snow depth and mass correction
		fprintf(fout,",%lf,%lf", M_TO_CM(Sdata->dhs_corr), Sdata->mass[SN_SURFACE_DATA::MS_CORRECTION]);
	} else {
		fprintf(fout,",,");
	}
}

/**
 * @brief Antarctic: last 7 time series (columns 94 to 100) dumped to *.met output file
 * @author Charles Fierz
 * @version 10.02
 * @param *Xdata
 * @param *Sdata
 * @param *Mdata
 * @param *Hdata
 * @param *fout Output file
 */
void AsciiIO::writeFreeSeriesANTARCTICA(const SN_STATION_DATA *Xdata, const SN_SURFACE_DATA *Sdata,
								const SN_MET_DATA *Mdata, const Q_PROCESS_DAT *Hdata, FILE *fout)
{
	(void) *Hdata;
	if ( OUT_METEO ) {
		// 94-98: 6 Antarctic Special Time Series
		fprintf(fout,",%.2lf,%.2lf,%.1lf,%.1lf,%.1lf", 100.*Mdata->rh_ave, Mdata->vw_ave, Xdata->rho_slope, 
			   Mdata->rho_hn, M_TO_CM(Xdata->Ndata[Xdata->ErosionLevel+1].z - Xdata->cH));
	}
	if ( OUT_HEAT ) {
		// 99: change of internal energy (kJ m-2)
		fprintf(fout, ",%lf", Sdata->dIntEnergy / 1000.);
		// 100: sum of energy fluxes at surface (kJ m-2)
		fprintf(fout,",%lf", ((Sdata->qw + Sdata->lw_net + Sdata->qs + Sdata->ql + Sdata->qr) * 
						  M_TO_S(calculation_step_length)) / 1000.);
	}
}

/**
 * @brief Calibration: last 7 time series (columns 94 to 100) dumped to *.met output file
 * @author Charles Fierz
 * @version 10.04
 * @param *Xdata
 * @param *Sdata
 * @param *Mdata
 * @param *Hdata
 * @param *fout Output file
 */
void AsciiIO::writeFreeSeriesCALIBRATION(const SN_STATION_DATA *Xdata, const SN_SURFACE_DATA *Sdata, 
								 const SN_MET_DATA *Mdata, const Q_PROCESS_DAT *Hdata, FILE *fout)
{
	(void) *Hdata; (void) *Sdata;
	double rho;
	const double t_surf = MIN(C_TO_K(-0.1), Xdata->Ndata[Xdata->getNumberOfNodes()-1].T);
	// 94-97: 5 new snow densities zwart, newLe, bellaire, crocus, hendrikx
	if( Xdata->hn_slope > 0. ) {
		for(unsigned int i=0; i<5; i++) {
			if ( i != Snowpack::LEHNING_OLD ) {
				rho = Snowpack::calculateNewSnowDensity(*Mdata, *Xdata, t_surf, 0., Snowpack::NewSnowDensityModel(i));
				fprintf(fout,",%.1lf", rho);
			}
		}
		if ( Mdata->vw > 2.9 ) {
			rho = Snowpack::sn_NewSnowDensityHendrikx(K_TO_C(Mdata->ta), K_TO_C(t_surf), 100.*Mdata->rh, Mdata->vw);
			fprintf(fout,",%.1lf", rho);
		} else {
			fprintf(fout,",0.0");
		}
	} else {
		fprintf(fout,",0.0,0.0,0.0,0.0,0.0");
	}
	if ( OUT_HEAT ) {
		// 99: // change of internal energy (kJ m-2)
		fprintf(fout, ",%lf", Sdata->dIntEnergy / 1000.);
		// 100: sum of energy fluxes at surface (kJ m-2)
		fprintf(fout,",%lf", ((Sdata->qw + Sdata->lw_net + Sdata->qs + Sdata->ql + Sdata->qr) * 
						  M_TO_S(calculation_step_length)) / 1000.);
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
 * @param *ext File extension
 * @return status
 */
bool AsciiIO::checkHeader(const char *fnam, const char *first_string, const Q_PROCESS_DAT* Hdata, 
					    const std::string& station, const char *ext, ...)
{
	FILE *fin=NULL;
	FILE *fout=NULL;
	va_list argptr; // get an arg ptr
	char *va_char;
	SN_STATION_DATA *va_Xdata;
	char dummy[MAX_STRING_LENGTH]="\000", dummy_l[MAX_LINE_LENGTH]="\000";
	int  i, j;

	if ( (fin = fopen(fnam, "r")) ) {
		// Check header of existing file
		fgets(dummy_l, MAX_LINE_LENGTH, fin);
		sscanf(dummy_l, "%s", dummy);
		if ( (strcmp(dummy, first_string) != 0) ) {
			prn_msg(__FILE__, __LINE__, "err", -1., "Header in %s should read %s, not %s", fnam, first_string, dummy);
			return false;
		}
		fclose(fin);
	} else if ( (strcmp(ext, "none") == 0) ) {
		// Check header only!
		return -1;
	} else {
		if ( !(fout = fopen(fnam, "w")) ) {
			return false;
		}
		// Initialize argptr to point to the first argument after the ext string
		va_start(argptr, ext);

		if ( (strcmp(ext, "err") == 0) ) {
			fprintf(fout, "[SNOWPACK_ERROR_LOG]");
			fprintf(fout, "\n          RUNTIME :  STN LOC LINE MSG [JULIAN]");
		} else if ( (strcmp(ext, "met") == 0) ) {
			va_Xdata = va_arg(argptr, SN_STATION_DATA *);
			fprintf(fout, "[STATION_PARAMETERS]");
			fprintf(fout, "\nStationName= %s",   station.c_str());
			fprintf(fout, "\nLatitude= %.2lf",   va_Xdata->Lat);
			fprintf(fout, "\nLongitude= %.2lf",  va_Xdata->Lon);
			fprintf(fout, "\nAltitude= %.0lf",   va_Xdata->Alt);
			fprintf(fout, "\nSlopeAngle= %.2lf", RAD_TO_DEG(va_Xdata->SlopeAngle));
			fprintf(fout, "\nSlopeAzi= %.2lf",   RAD_TO_DEG(va_Xdata->SlopeAzi));
			fprintf(fout, "\nDepthTemp= %1d",    useSnowLayers);
			for (i = 0; i < FIXED_HEIGHTS; i++) {
				fprintf(fout, ",%.3lf", depth_of_sensors[i]);
			}
			fprintf(fout, "\n\n[HEADER]");
			if ( OUT_HAZ ) { // HACK To avoid troubles in A3D
				fprintf(fout, "\n#%s, Snowpack %s version %s run by \"%s\"", Hdata->sn_computation_date,
					   variant.c_str(),Hdata->sn_version,Hdata->sn_user);
			}
			fprintf(fout, "\n,,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100");
			fprintf(fout, "\nID,Date,Sensible heat,Latent heat,Outgoing longwave radiation,Incoming longwave radiation,Net absorbed longwave radiation,Reflected shortwave radiation,Incoming shortwave radiation,Net absorbed shortwave radiation,Modelled surface albedo,Air temperature,Modeled surface temperature,Measured surface temperature,Temperature at bottom of snow or soil pack,Heat flux at bottom of snow or soil pack,Ground surface temperature,Heat flux at ground surface,Heat advected to the surface by liquid precipitation,Global solar radiation (horizontal)");
			fprintf(fout, ",Global solar radiation on slope,Direct solar radiation on slope,Diffuse solar radiation on slope,Measured surface albedo,Relative humidity,Wind speed,Max wind speed at snow station or wind speed at ridge station,Wind direction at snow station,Precipitation rate at surface (solid only),Modelled snow depth (vertical),Measured snow depth (vertical),Surface hoar size,24h Drift index (vertical),Height of new snow HN (24h vertical),3d sum of daily height of new snow (vertical),Total snowpack mass,Eroded mass,Rain rate,Surface runoff (without soil infiltration)");
			fprintf(fout, ",Sublimation,Evaporation,Temperature 1 (modelled),Temperature 1 (measured),Temperature 2 (modelled),Temperature 2 (measured),Temperature 3 (modelled),Temperature 3 (measured),Temperature 4 (modelled),Temperature 4 (measured),Temperature 5 (modelled),Temperature 5 (measured)");
			if ( MAX_NUMBER_SENSORS == 5 ) {
				fprintf(fout, ",Solute load at soil surface,SWE (of snowpack),Liquid Water Content (of snowpack),Profile type,Stability class,z_Sdef,Deformation rate stability index Sdef,z_Sn38,Natural stability index Sn38,z_Sk38,Skier stability index Sk38,z_SSI,Structural Stability index SSI,z_S5,Stability index S5");
				if ( useCanopyModel && OUT_CANOPY ) {
					fprintf(fout, ",Interception storage,Canopy surface  temperature,Canopy albedo,Wet fraction,Interception capacity,Net shortwave radiation absorbed by canopy,Net longwave radiation absorbed by canopy,Net radiation canopy,Sensible heat flux into the canopy,Latent heat flux into the canopy,Transpiration of the canopy,Evaporation and sublimation of interception (liquid and frozen),Interception rate,Throughfall,Snow unload,Sensible heat flux to the canopy,Latent heat flux to the canopy,Longwave radiation up above canopy,Longwave radiation down above canopy");
					fprintf(fout, ",Net longwave radiation above canopy,Shortwave radiation up above canopy,Shortwave radiation down above canopy,Net shortwave radiation above canopy,Total land surface albedo,Total net radiation,Surface (radiative) temperature,Precipitation Above Canopy,Total Evapotranspiration");
				} else {
					fprintf(fout,",-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-");
				}
			} else if ( OUT_T ) {
				int i_prn;
				j = 0;
				for (i = MIN(5, FIXED_HEIGHTS); i < FIXED_HEIGHTS+FIXED_RATES; i++) {
					if ( i < FIXED_HEIGHTS ) {
						i_prn = i + 1;
						fprintf(fout, ",Temperature %d (modelled)", i_prn);
					} else {
						i_prn = (i-FIXED_HEIGHTS)+1;
						fprintf(fout, ",Hfr %d", i_prn);
						fprintf(fout, ",Tfr %d (modelled)", i_prn);
						j++;
					}
					if ( i < t_internal ) {
						if ( i < FIXED_HEIGHTS ) {
							fprintf(fout, ",Temperature %d (measured)", i_prn);
						} else {
							fprintf(fout, ",Tfr %d (measured)", i_prn);
						}
					} else {
						fprintf(fout, ",");
					}
					j += 2;
				}
				if ( va_Xdata->tag_low ) {
					int tag = va_Xdata->tag_low, j_lim;
					while ( (tag + i) <= NUMBER_SENSORS ) {
						if ( (tag + i) <= t_internal ) {
							j_lim = 41;
						} else {
							j_lim = 43;
						}
						if ( j < j_lim ) {
							fprintf(fout, ",H(tag%02d),T(tag%02d)", tag, tag);
							j += 2;
							if ( i < t_internal ) {
								fprintf(fout, ",H(meas%02d),T(meas%02d)", tag, tag);
								j += 2;
							}
							tag++;
						}
					}
				}
				for (; j < 44; j++) {
					fprintf(fout,",-");
				}
			} else {
				fprintf(fout, ",-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-");
			}

			if (variant == "ANTARCTICA"){
				fprintf(fout, ",Running mean relative humidity (100h),Running mean wind speed (100h),Modeled new snow density,Measured new snow density,Erosion level (from srf),Internal Energy change,Sum surface fluxes");
			} else if (variant == "CALIBRATION"){
				fprintf(fout, ",rho_hn(Zwart),rho_hn(Lehning),rho_hn(Bellaire),rho_hn(crocus),rho_hn(Hendrikx),Internal Energy change,Sum surface fluxes");
			} else {
				fprintf(fout, ",Soil runoff,Measured sensible heat,Measured latent heat,Crust thickness (S-slope),Internal energy change,Sum surface fluxes,free4,free5");
			}

			fprintf(fout, "\n,,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,1,degC,degC,degC,degC,W m-2,degC,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,1,%%,m s-1,m s-1,deg,kg m-2 h-1,cm,cm,mm,cm,cm,cm,kg m-2,kg m-2 h-1,kg m-2 h-1,kg m-2,kg m-2,kg m-2,degC,degC,degC,degC,degC,degC,degC,degC,degC,degC");
			if ( MAX_NUMBER_SENSORS == 5 ) {
				fprintf(fout, ",kg m-2,kg m-2,kg m-2,-,-,cm,1,cm,1,cm,1,cm,1,cm,1");
				if ( OUT_CANOPY && useCanopyModel ) {
					fprintf(fout, ",kg m-2,degC,-,-,kg m-2,W m-2,W m-2,W m-2,W m-2,W m-2,kg m-2 per timestep,kg m-2 per timestep,kg m-2,kg m-2,kg m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,degC,kg m-2,kg m-2 per timestep");
				} else {
					fprintf(fout,",,,,,,,,,,,,,,,,,,,,,,,,,,,,");
				}
			} else if ( OUT_T ) {
				j = 0;
				for (i = MIN(5, FIXED_HEIGHTS); i < FIXED_HEIGHTS+FIXED_RATES; i++) {
					if ( i >= FIXED_HEIGHTS ) {
						fprintf(fout, ",cm");
						j++;
					}
					fprintf(fout, ",degC");
					j++;
					if ( i < t_internal ) {
						fprintf(fout, ",degC");
						j++;
					}
				}
				if ( va_Xdata->tag_low ) {
					int tag = va_Xdata->tag_low, j_lim;
					while ( (tag + i) <= NUMBER_SENSORS ) {
						if ( (tag + i) <= t_internal ) {
							j_lim = 41;
						} else {
							j_lim = 43;
						}
						if ( j < j_lim ) {
							fprintf(fout, ",cm,degC");
							j += 2;
							if ( i < t_internal ) {
								fprintf(fout, ",cm,degC");
								j += 2;
							}
							tag++;
						}
					}
				}
				//while ( j < 43 ) {
				//fprintf(fout, ",cm,degC");
				//j += 2;
				//}
				for (; j < 44; j++) {
					fprintf(fout,",");
				}
			} else {
				fprintf(fout, ",,,,,,,,,,,,,,,,,,,,,,,,,,,,");

			}
			if (variant == "ANTARCTICA"){
				fprintf(fout, ",%%,m s-1,kg m-3,kg m-3,cm,kJ m-2,kJ m-2");
			} else if (variant == "CALIBRATION"){
				fprintf(fout, ",kg m-3,kg m-3,kg m-3,kg m-3,kg m-3,kJ m-2,kJ m-2");
			} else {
				fprintf(fout, ",kg m-2,W m-2,W m-2,cm,kJ m-2,kJ m-2,,");
			}

			fprintf(fout, "\n\n[DATA]");
		} else if ( (strcmp(ext, "pev") == 0) ) {
			va_char = va_arg(argptr, char *);
			fprintf(fout,"ProfEval @ station %s/%s for experiment %s\n",
				   station.c_str(), va_char, experiment.c_str());
			if ( T_SRF && T_GND ) {
				fprintf(fout, "(Tsrf and Tgnd BOTH included!)\n");
			} else if ( T_SRF ) {
				fprintf(fout, "(Tsrf included, Tgnd NOT included!)\n");
			} else if ( T_GND ) {
				fprintf(fout, "(Tsrf NOT included, Tgnd included!)\n");
			} else {
				fprintf(fout, "(NEITHER Tsrf NOR Tgnd included!)\n");
			}
			fprintf(fout,"Prf#\t     Termin     ");
			fprintf(fout,"\tshape\t size\twetness\thardness\t   hw\t  rho\t th_w\t    T\tT_50cm\toverall");
			fprintf(fout,"\tSLhsw/MAhsw\tSLhs/MAhs\tSL.hsw\tMA.hsw\t SL.hs\t MA.hs");
			fprintf(fout,"\tMA.nL\tMA.nHW\tMA.nRHO\tMA.nTH_W\tMA.nT\tMA.nTmD");
			fprintf(fout,"\tSL.nL\tSL.nHW\tSL.nRHO\tSL.nTH_W\tSL.nT\tSL.nTmD\n");
		} else if ( (strcmp(ext, "pro") == 0) ) {
			va_Xdata = va_arg(argptr, SN_STATION_DATA *);
			fprintf(fout, "[STATION_PARAMETERS]");
			fprintf(fout, "\nStationName= %s",   station.c_str());
			fprintf(fout, "\nLatitude= %.2lf",   va_Xdata->Lat);
			fprintf(fout, "\nLongitude= %.2lf",  va_Xdata->Lon);
			fprintf(fout, "\nAltitude= %.0lf",   va_Xdata->Alt);
			fprintf(fout, "\nSlopeAngle= %.2lf", RAD_TO_DEG(va_Xdata->SlopeAngle));
			fprintf(fout, "\nSlopeAzi= %.2lf",   RAD_TO_DEG(va_Xdata->SlopeAzi));

			fprintf(fout, "\n\n[HEADER]");
			if( OUT_HAZ ) { // HACK To avoid troubles in A3D
				fprintf(fout, "\n#%s, Snowpack %s version %s run by \"%s\"", 
					   Hdata->sn_computation_date, variant.c_str(), Hdata->sn_version, Hdata->sn_user);
			}
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
			fprintf(fout, "\n\n[DATA]");
		} else {
			prn_msg(__FILE__, __LINE__, "wrn", -1., "No header defined for files *.%s", ext);
		}
		va_end(argptr);
		fclose(fout);
	}
	
	return true;
}

void AsciiIO::writeHazardData(const std::string& station, const std::vector<Q_PROCESS_DAT>& Hdata, 
						const std::vector<Q_PROCESS_IND>& Hdata_ind, const int& num)
{
	throw IOException("Nothing implemented here!", AT);
}

