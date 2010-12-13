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

#include <snowpack/ImisDBIO.h>

using namespace std;
using namespace mio;
using namespace oracle;
using namespace oracle::occi;

const double ImisDBIO::in_tz = 1.; //All IMIS data is in gmt+1

bool ImisDBIO::research_mode = true;
double ImisDBIO::density_hoar_surf = 0.0;
double ImisDBIO::min_size_hoar_surf = 0.0;

const string ImisDBIO::sqlDeleteHdata = "DELETE FROM snowpack.ams_pmod WHERE stat_abk=:1 and stao_nr=:2 and wstao_nr = :3 and datum>=:4 and datum<=:5";

const string ImisDBIO::sqlInsertHdata = "INSERT INTO snowpack.ams_pmod(datum,stat_abk,stao_nr,wstao_nr,dewpt_def,hoar_ind6,hoar_ind24,wind_trans,hns3,hns6,hns12,hns24,hns72,hns72_24,wc3,wc6,wc12,wc24,wc72,hoar_size,wind_trans24,stab_class1,stab_class2,stab_index1,stab_height1,stab_index2,stab_height2,stab_index3,stab_height3,stab_index4,stab_height4,stab_index5,stab_height5,ch,crust,en_bal,sw_net,t_top1,t_top2,snowpack_version,calc_date,swe,tot_lwc,runoff) values (:1,:2,:3,:4,:5,:6,:7,:8,:9,:10,:11,:12,:13,:14,:15,:16,:17,:18,:19,:20,:21,:22,:23,:24,:25,:26,:27,:28,:29,:30,:31,:32,:33,:34,:35,:36,:37,:38,:39,:40,:41,:42,:43,:44)";

string ImisDBIO::oracleDB = "";
string ImisDBIO::oracleUser = "";
string ImisDBIO::oraclePassword = "";

const std::string ImisDBIO::profile_filename = "loaddata/pmodpro.dat";

ImisDBIO::ImisDBIO(const mio::Config& i_cfg) : cfg(i_cfg)
{
	cfg.getValue("DBNAME", "Output", oracleDB, Config::nothrow);
	cfg.getValue("DBUSER", "Output", oracleUser, Config::nothrow);
	cfg.getValue("DBPASS", "Output", oraclePassword, Config::nothrow);

	research_mode = cfg.get("RESEARCH", "Parameters");

	//Density of surface hoar (-> hoar index of surface node) (kg m-3)
	density_hoar_surf = cfg.get("DENSITY_HOAR_SURF", "Parameters");

	//Minimum size to show surface hoar on surface (mm)
	min_size_hoar_surf = cfg.get("MIN_SIZE_HOAR_SURF", "Parameters");
}

void ImisDBIO::readSnowCover(const std::string& /*station*/, SN_SNOWSOIL_DATA& /*SSdata*/, SN_ZWISCHEN_DATA& /*Zdata*/)
{
	throw IOException("Nothing implemented here!", AT);
}

void ImisDBIO::writeSnowCover(const mio::Date& /*date*/, const std::string& /*station*/, const SN_STATION_DATA& /*Xdata*/, 
                              const SN_ZWISCHEN_DATA& /*Zdata*/, const bool& /*forbackup*/)
{
	throw IOException("Nothing implemented here!", AT);
}
	
void ImisDBIO::writeTimeSeries(const std::string& /*station*/, const SN_STATION_DATA& /*Xdata*/, 
                               const SN_SURFACE_DATA& /*Sdata*/, const SN_MET_DATA& /*Mdata*/, const Q_PROCESS_DAT& /*Hdata*/)
{
	throw IOException("Nothing implemented here!", AT);
}

/**
 * @brief Dump snow profile to ASCII file for subsequent upload to SDBO
 */
void ImisDBIO::writeProfile(const mio::Date& date, const std::string& station, const unsigned int& expo,
                            const SN_STATION_DATA& Xdata, const Q_PROCESS_DAT& Hdata)
{

	if ((research_mode) || (expo != 0)) //write output only for the flat field (expo == 0)
		return;

	FILE *PFile=NULL;
	int n1, e, l = 0,  nL = 0, nE ;
	vector<Q_PROFILE_DAT> Pdata;

	const vector<SN_ELEM_DATA>& EMS = Xdata.Edata; nE = Xdata.getNumberOfElements();
	const vector<SN_NODE_DATA>& NDS = Xdata.Ndata; 

	// Allocate Memory, use double memory for possible surface hoar layers
	Pdata.resize(nE);

	// Generate the profile data from the element data (1 layer = 1 element)
	int snowloc = 1;
	string mystation = station;
	if (isdigit(station[station.length()-1])){
	    snowloc = station[station.length()-1] - '0';
	    if (mystation.length() > 2)
		    mystation = mystation.substr(0, mystation.length()-1);
	}
	

	for(e=0; e<nE; e++) {
		Pdata[l].date = date;
		Pdata[l].stationname = mystation;
		Pdata[l].loc_for_snow = snowloc;
		Pdata[l].loc_for_wind = 1;

		//write version and date
		Pdata[l].sn_version = Hdata.sn_version;
		Pdata[l].sn_computation_date = Hdata.sn_computation_date;
		Pdata[l].sn_jul_computation_date=Hdata.sn_jul_computation_date;
		n1=e+1;
		Pdata[l].height = M_TO_CM(NDS[n1].z + NDS[n1].u);
		Pdata[l].layer_date = EMS[e].date.getJulianDate();
		Pdata[l].rho = EMS[e].Rho;
		Pdata[l].tem = K_TO_C(EMS[e].Te);
		Pdata[l].tem_grad = EMS[e].gradT;
		Pdata[l].strain_rate = fabs(EMS[e].EDot);
		Pdata[l].theta_w = EMS[e].theta[WATER] * 100.;
		Pdata[l].theta_i = EMS[e].theta[ICE] * 100.;
		Pdata[l].dendricity = EMS[e].dd;
		Pdata[l].sphericity = EMS[e].sp;
		Pdata[l].coordin_num = EMS[e].N3;
		Pdata[l].grain_dia = 2. * EMS[e].rg;
		Pdata[l].bond_dia = 2. * EMS[e].rb;
		Pdata[l].marker = EMS[e].mk%100;
		Pdata[l].hard = EMS[e].hard;
		l++;
	}

	if ( (nL = Aggregate::aggregate(Pdata) ) < 0 ) {
		prn_msg(__FILE__, __LINE__, "err", date.getJulianDate(), "Cannot aggregate layers");
		throw IOException("Cannot aggregate layers", AT);
	}

	if ( !(PFile = fopen(profile_filename.c_str(), "a")) ) {
		prn_msg(__FILE__, __LINE__, "err", date.getJulianDate(), "Cannot open Profile file: %s", 
			   profile_filename.c_str());
		throw FileAccessException(profile_filename, AT);
	}

	for(e=0; e<nL; e++) {
		//HACK: these legacy offset should be removed.
		//This means specify a different import date format for the database and remove the offset here
		const double profile_date = Pdata[e].date.getJulianDate() - 2415021. + 0.5; //HACK
		const double layer_date = Pdata[e].layer_date - 2415021. + 0.5; //HACK

		fprintf(PFile,"%.5lf,%s,%d,%d,%.2lf,", profile_date, Pdata[e].stationname.c_str(),
			   Pdata[e].loc_for_snow, Pdata[e].loc_for_wind, Pdata[e].height);
		fprintf(PFile,"%.5lf,%.0lf,%.1lf,%.0lf,%.4e,%.0lf,%.0lf,%.2lf,%.2lf,%.1lf,%.1lf,%.2lf,%d\n", 
			   layer_date, Pdata[e].rho, Pdata[e].tem, Pdata[e].tem_grad, Pdata[e].strain_rate,
			   Pdata[e].theta_w, Pdata[e].theta_i, Pdata[e].dendricity, Pdata[e].sphericity, 
			   Pdata[e].coordin_num, Pdata[e].grain_dia, Pdata[e].bond_dia, Pdata[e].grain_class);
	}

	//HACK: The condition nL < nE added by Egger: the aggregation of the layers in Aggregate::aggregate, 
	//disregards the top layer only if nE > 5, otherwise nL might equal nE and the index e is invalid
	//It's a hack because if the heighest layer is actually hoar then the former loop should only loop until
	//nL - 1, that is if nL == nE
	if ((NDS[nE].hoar > MM_TO_M(min_size_hoar_surf) * density_hoar_surf) && (nL < nE)) {
		//HACK: these legacy offset should be removed.
		//This means specify a different import date format for the database and remove the offset here
		const double profile_date = Pdata.at(e).date.getJulianDate() - 2415021. + 0.5; //HACK
		const double layer_date = Pdata.at(e).layer_date - 2415021. + 0.5; //HACK

		double gsz_SH = M_TO_MM(NDS[nE].hoar / density_hoar_surf);
		e=nL-1;
		fprintf(PFile,"%.5lf,%s,%d,%d,%.2lf,", profile_date, Pdata[e].stationname.c_str(),
			   Pdata[e].loc_for_snow, Pdata[e].loc_for_wind, Pdata[e].height + MM_TO_CM(gsz_SH));
		fprintf(PFile,"%.5lf,%.0lf,%.1lf,%.0lf,%.4e,%.0lf,%.0lf,%.2lf,%.2lf,%.1lf,%.1lf,%.2lf,%d\n", 
			   layer_date, density_hoar_surf, Pdata[e].tem, Pdata[e].tem_grad, Pdata[e].strain_rate,
			   0., density_hoar_surf/Constants::density_ice, 0., 0., 2., gsz_SH, 0.6667*gsz_SH, 660);
	}

	fclose(PFile);
}

void ImisDBIO::writeHazardData(const std::string& station, const vector<Q_PROCESS_DAT>& Hdata, 
                               const vector<Q_PROCESS_IND>& Hdata_ind, const int& num)
{
	if ((num < 0) || (num >= (int)Hdata.size())){
		cout << "\tNo hazard data inserted or deleted from DB" << endl;
		return; //nothing to do
	}

	if ((oracleDB == "") || (oraclePassword == "") || (oracleUser == ""))
		throw IOException("You must set the output database, username and password", AT);

	string stationName="", stationNumber="";
	parseStationName(station, stationName, stationNumber);

	Environment *env = NULL;
	
	try {
		env = Environment::createEnvironment();// static OCCI function
		Connection *conn = NULL;

		conn = env->createConnection(oracleUser, oraclePassword, oracleDB);

		deleteHdata(stationName, stationNumber, Hdata[0].date, Hdata[num-1].date, env, conn);
		insertHdata(stationName, stationNumber, Hdata, Hdata_ind, num, env, conn);

		env->terminateConnection(conn);
		Environment::terminateEnvironment(env); // static OCCI function

	} catch (exception& e){
		Environment::terminateEnvironment(env); // static OCCI function
		throw IOException("Oracle Error: " + string(e.what()), AT); //Translation of OCCI exception to IOException
	}
}

/**
 * @brief This function breaks up the station name into two components (a string and a number e.g. KLO2 -> "KLO","2")
 * @param stationName The full name of the station (e.g. "KLO2")
 * @param stName      The string part of the name  (e.g. "KLO")
 * @param stNumber    The integer part of the name (e.g. "2")
 */
void ImisDBIO::parseStationName(const std::string& stationName, std::string& stName, std::string& stNumber)
{
	stName    = stationName.substr(0, stationName.length()-1); //The station name: e.g. KLO
	stNumber  = stationName.substr(stationName.length()-1, 1); //The station number: e.g. 2
	if(!std::isdigit(stNumber[0])) {
		//the station is one of these non-imis stations that don't contain a number...
		stName = stationName;
		stNumber = "0";
	}
}

void ImisDBIO::deleteHdata(const std::string& stationName, const std::string& stationNumber,
					  const mio::Date& dateStart, const mio::Date& dateEnd, 
					  oracle::occi::Environment*& env, oracle::occi::Connection*& conn)
{
	vector< vector<string> > vecResult;
	vector<int> datestart = vector<int>(5);
	vector<int> dateend   = vector<int>(5);

	//IMIS is in TZ=+1, so moving back to this timezone
	mio::Date dateS(dateStart), dateE(dateEnd);
	//dateS.setTimeZone(in_tz);
	//dateE.setTimeZone(in_tz);
	dateS.getDate(datestart[0], datestart[1], datestart[2], datestart[3], datestart[4]);
	dateE.getDate(dateend[0], dateend[1], dateend[2], dateend[3], dateend[4]);

	//Oracle can't deal with an integer for the hour of 24, hence the following workaround
	if (datestart[3] == 24){
		mio::Date tmpDate = dateStart + mio::Date(3.0/(60*60*24)); //add three seconds to omit 24 for 00 
		//tmpDate.setTimeZone(in_tz);
		tmpDate.getDate(datestart[0], datestart[1], datestart[2], datestart[3], datestart[4]);
	}
	if (dateend[3] == 24){
		mio::Date tmpDate = dateEnd + mio::Date(3.0/(60*60*24)); //add three seconds to omit 24 for 00 
		//tmpDate.setTimeZone(in_tz);
		tmpDate.getDate(dateend[0], dateend[1], dateend[2], dateend[3], dateend[4]);
	}

	Statement *stmt  = conn->createStatement(sqlDeleteHdata);

	// construct the oracle specific Date object: year, month, day, hour, minutes
	occi::Date begindate(env, datestart[0], datestart[1], datestart[2], datestart[3], datestart[4]);
	occi::Date enddate(env, dateend[0], dateend[1], dateend[2], dateend[3], dateend[4]);
	stmt->setString(1, stationName);   // set 1st variable's value (station name)
	stmt->setString(2, stationNumber); // set 2nd variable's value (station number)
	stmt->setString(3, stationNumber); // set 3rd variable's value (station number)
	stmt->setDate(4, begindate);       // set 4rd variable's value (begin date)
	stmt->setDate(5, enddate);         // set 5th variable's value (enddate)

	unsigned int rows_deleted = stmt->executeUpdate();
	conn->terminateStatement(stmt);

	cout << "\tDeleted " << rows_deleted << " rows in DB " << oracleDB << endl;
}

void ImisDBIO::insertHdata(const std::string& stationName, const std::string& stationNumber,
					  const std::vector<Q_PROCESS_DAT>& Hdata, const std::vector<Q_PROCESS_IND>& Hdata_ind, 
					  const int& num, oracle::occi::Environment*& env, oracle::occi::Connection*& conn)
{
	vector<int> sndate = vector<int>(5);
	unsigned int rows_inserted = 0;
	double sn_version;
	IOUtils::convertString(sn_version, Hdata[0].sn_version);

	mio::Date dateSn(Hdata[0].sn_compile_date);
	//dateSn.setTimeZone(in_tz); //IMIS is in TZ=+1, so moving back to this timezone
	dateSn.getDate(sndate[0], sndate[1], sndate[2], sndate[3], sndate[4]);
	if (sndate[3] == 24){
		mio::Date tmpDate = dateSn + mio::Date(3.0/(60*60*24)); //add three seconds to omit 24 for 00 
		//tmpDate.setTimeZone(in_tz);
		tmpDate.getDate(sndate[0], sndate[1], sndate[2], sndate[3], sndate[4]);
	}

	for (int i = 0; i<num; i++){ 		
		if (Hdata[i].date == mio::Date()) break; //catch the case that not all Hdata has been set properly
		
		vector<int> hzdate = vector<int>(5);
		mio::Date dateH(Hdata[i].date);
		//dateH.setTimeZone(in_tz); //IMIS is in TZ=+1, so moving back to this timezone
		dateH.getDate(hzdate[0], hzdate[1], hzdate[2], hzdate[3], hzdate[4]);

		//Oracle can't deal with an integer for the hour of 24, hence the following workaround
		if (hzdate[3] == 24){
			mio::Date tmpDate = dateH + mio::Date(3.0/(60*60*24)); //add three seconds to omit 24 for 00 
			//tmpDate.setTimeZone(in_tz);
			tmpDate.getDate(hzdate[0], hzdate[1], hzdate[2], hzdate[3], hzdate[4]);
		}

		Statement *stmt  = conn->createStatement(sqlInsertHdata);

		// construct the oracle specific Date object: year, month, day, hour, minutes
		occi::Date hazarddate(env, hzdate[0], hzdate[1], hzdate[2], hzdate[3], hzdate[4]);
		occi::Date computationdate(env, sndate[0], sndate[1], sndate[2], sndate[3], sndate[4]);

		int statNum = 0;
		IOUtils::convertString(statNum, stationNumber);

		stmt->setAutoCommit(true);

		unsigned int param = 1;

		stmt->setDate(param++, hazarddate);
		stmt->setString(param++, stationName);
		stmt->setNumber(param++, statNum);
		stmt->setNumber(param++, statNum);
		
		if (Hdata_ind[i].dewpt_def != -1)  stmt->setNumber(param++, Hdata[i].dewpt_def); 
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hoar_ind6 != -1)	stmt->setNumber(param++, Hdata[i].hoar_ind6);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hoar_ind24 != -1) stmt->setNumber(param++, Hdata[i].hoar_ind24);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].wind_trans != -1) stmt->setNumber(param++, Hdata[i].wind_trans);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hns3 != -1)       stmt->setNumber(param++, Hdata[i].hns3);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hns6 != -1)       stmt->setNumber(param++, Hdata[i].hns6);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hns12 != -1)      stmt->setNumber(param++, Hdata[i].hns12);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hns24 != -1)      stmt->setNumber(param++, Hdata[i].hns24);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hns72 != -1)      stmt->setNumber(param++, Hdata[i].hns72);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].hns72_24 != -1)   stmt->setNumber(param++, Hdata[i].hns72_24);
		else stmt->setNull(param++, occi::OCCINUMBER);
		
		if (Hdata_ind[i].wc3 != -1)        stmt->setNumber(param++, Hdata[i].wc3);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].wc6 != -1)        stmt->setNumber(param++, Hdata[i].wc6);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].wc12 != -1)       stmt->setNumber(param++, Hdata[i].wc12);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].wc24 != -1)       stmt->setNumber(param++, Hdata[i].wc24);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].wc72 != -1)       stmt->setNumber(param++, Hdata[i].wc72);
		else stmt->setNull(param++, occi::OCCINUMBER);

		if (Hdata_ind[i].hoar_size != -1)  stmt->setNumber(param++, Hdata[i].hoar_size);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].wind_trans24 != -1) stmt->setNumber(param++, Hdata[i].wind_trans24);
		else stmt->setNull(param++, occi::OCCINUMBER);

		if (Hdata_ind[i].stab_class1 != -1)  stmt->setNumber(param++, Hdata[i].stab_class1);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_class2 != -1)  stmt->setNumber(param++, Hdata[i].stab_class2);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_index1 != -1)  stmt->setNumber(param++, Hdata[i].stab_index1);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_height1 != -1) stmt->setNumber(param++, Hdata[i].stab_height1);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_index1 != -1)  stmt->setNumber(param++, Hdata[i].stab_index2);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_height1 != -1) stmt->setNumber(param++, Hdata[i].stab_height2);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_index1 != -1)  stmt->setNumber(param++, Hdata[i].stab_index3);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_height1 != -1) stmt->setNumber(param++, Hdata[i].stab_height3);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_index1 != -1)  stmt->setNumber(param++, Hdata[i].stab_index4);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_height1 != -1) stmt->setNumber(param++, Hdata[i].stab_height4);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_index1 != -1)  stmt->setNumber(param++, Hdata[i].stab_index5);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].stab_height1 != -1) stmt->setNumber(param++, Hdata[i].stab_height5);
		else stmt->setNull(param++, occi::OCCINUMBER);

		if (Hdata_ind[i].ch != -1)     stmt->setNumber(param++, Hdata[i].ch);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].crust != -1)  stmt->setNumber(param++, Hdata[i].crust);
		else stmt->setNull(param++, occi::OCCINUMBER);

		if (Hdata_ind[i].en_bal != -1) stmt->setNumber(param++, Hdata[i].en_bal);
		else stmt->setNull(param++, occi::OCCINUMBER);

		if (Hdata_ind[i].sw_net != -1)  stmt->setNumber(param++, Hdata[i].sw_net);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].t_top1 != -1)  stmt->setNumber(param++, Hdata[i].t_top1);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].t_top2 != -1)  stmt->setNumber(param++, Hdata[i].t_top2);
		else stmt->setNull(param++, occi::OCCINUMBER);

		stmt->setNumber(param++, sn_version);
		stmt->setDate(param++, computationdate);

		if (Hdata_ind[i].swe != -1)      stmt->setNumber(param++, Hdata[i].swe);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].tot_lwc != -1)  stmt->setNumber(param++, Hdata[i].tot_lwc);
		else stmt->setNull(param++, occi::OCCINUMBER);
		if (Hdata_ind[i].runoff != -1)   stmt->setNumber(param++, Hdata[i].runoff);
		else stmt->setNull(param++, occi::OCCINUMBER);

		rows_inserted += stmt->executeUpdate(); // execute the statement stmt
		conn->terminateStatement(stmt);
	}

	cout << "\tInserted " << rows_inserted << " rows in DB " << oracleDB << endl;
}
