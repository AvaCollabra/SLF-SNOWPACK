#include <meteoio/MeteoIO.h>
//#include <snowpack/snowpackCore/Snowpack.h>
#include "TestSnowpack.h"
#include "HeatEquationAnalytical.h"
#include <snowpack/libsnowpack.h>
#include <stdlib.h>
#include <fstream>
#include <cmath>

using namespace std;
using namespace mio;


// Forward declarations
size_t getNodeNumber(size_t elementNumber);

// PARAMETERS
const double tol_inf = 1e-6;  // Tolerance for the infinite-norm error
const double tol_2 = 1e-6;  // Tolerance for the 2-norm error




int main(int argc, char *argv[]) {

  /* Testing solver.h */
  if (argc < 1) {
    cout << "solverTest called without argument. Error.\n";
    cerr << "solverTest called without argument. Error.\n";
    exit(1);
  } else if (argc > 2) {
    cout << "solverTest called with too many arguments. Error.\n";
    cerr << "solverTest called with too many arguments. Error.\n";
    exit(1);
  }

  /* Parameters definition */
	const double totalThickness = 20.0;
  const double T0 = 273.15;
	const double T1 = 263.15;
	const double snowDensity = 100.0;
	const double thetaAir = 1.0 / Constants::conductivity_air;
  const size_t nTestLayers = 10000;
	const size_t timeStep = 1000;

  /* Needed variables */
  double layerThickness(0.0);

  /* Benchmark data from analytical solution */
  HeatEquationAnalytical he(
	    nTestLayers,
      totalThickness,
      T0,
      T1,
      Constants::density_air / Constants::conductivity_air
          * Constants::specific_heat_air / snowDensity,
      snowDensity, 1.0);


  /* Create static test data for comparison */
  // Config file
  static string cfgfile = "io.ini";
  SnowpackConfig cfg(cfgfile);
  cfg.addKey("METEO_STEP_LENGTH", "Snowpack", "1");

  // Define snowpack test object
  TestSnowpack testsnowpack(cfg);

  // Set surface bc
  testsnowpack.setBC(Snowpack::DIRICHLET_BC);

  // Create neutral weather
  CurrentMeteo Mdata(cfg);
  Mdata.iswr = Mdata.rswr = 0.0;
  Mdata.ts0 = T0;
  Mdata.ustar = 0.0;
  Mdata.z0 = 1.0;

  // Build static part of snow layer
  SnowStation Xdata(false, false);
  BoundCond Bdata;


  /* Create container for */
  size_t layersVectorTmp[] = { 5, 10, 50, 100, 500, 1000, 5000 };
	size_t timeVectorTmp[] = { 1 };

  vector<size_t> layersVector(
      layersVectorTmp,
      layersVectorTmp + sizeof(layersVectorTmp) / sizeof(size_t));

	vector<size_t> timeVector(
	    timeVectorTmp, timeVectorTmp + sizeof(timeVectorTmp) / sizeof(size_t));

	map<size_t, map<size_t, Error> > resultsMap;

	/* Loop over time step lengths */
	for (vector<size_t>::iterator it2 = timeVector.begin();
			it2 != timeVector.end(); it2++) {

		// Change time step
		testsnowpack.setSnDt(static_cast<double>(*it2));

		// Declare fresh first level map
		map<size_t, Error> tmpResultsMap;


		/* Loop over number of layers */
		for (vector<size_t>::iterator it = layersVector.begin();
		    it != layersVector.end(); it++) {

			// Clear element and node vector
			Xdata.Ndata.clear();
			Xdata.Edata.clear();

			// Resize SnowStation data to adequate size
			Xdata.resize(*it);

			// Actual layerThickness
			layerThickness = totalThickness / (double) *it;

			// Set all nodes to same temp
			for (size_t n = 0; n < Xdata.getNumberOfNodes(); n++) {
				Xdata.Ndata[n].T = T0;
			}

			Xdata.Ndata[Xdata.getNumberOfNodes() - 1].T = T1;

			// Set all elements to same property
			for (size_t n = 0; n < Xdata.getNumberOfElements(); n++) {

				Xdata.Edata[n].Rho = snowDensity;
				Xdata.Edata[n].theta[ICE] = 0.0;
				Xdata.Edata[n].theta[WATER] = 0.0;
				Xdata.Edata[n].theta[WATER_PREF] = 0.0;
				Xdata.Edata[n].L = layerThickness;
				Xdata.Edata[n].theta[AIR] = thetaAir;
				Xdata.Edata[n].res_wat_cont = 0.00001;
				Xdata.Edata[n].Te = T0;

			}

			/* Solve temp diffusion */
			int ii = 0;		// Counter for sub-timesteps to match one SNOWPACK time step
			bool LastTimeStep = false;	// Flag to indicate if it is the last sub-time step
			double p_dt = 0.;			// Cumulative progress of time steps

			// Computation of temp profile
			while (LastTimeStep == false) {
				if (testsnowpack.testCompTempProfile(Mdata, Xdata, Bdata, false)) {
					// Entered after convergence
					ii++;						// Update time step counter
					p_dt += testsnowpack.getSnDt();					// Update progress variable
					if (p_dt > timeStep - Constants::eps) {  // Check if it is the last sub-time step
						LastTimeStep = true;
					}
				}
			}

			// Store the error in first map for layer "it"
			Error actualError;
			he.computeErrors(timeStep, *it, Xdata.Ndata, actualError);
			tmpResultsMap.insert(
			    tmpResultsMap.end(),
			                     pair<size_t, Error>(
			        *it, actualError));

			}

		// Store error map in main map
		resultsMap.insert(resultsMap.end(),
		                  pair<size_t, map<size_t, Error> >(*it2, tmpResultsMap));

		}




  /* Output to file */
  // Open input stream
  string str = "myoutput.dat";
  std::ofstream ofs;
  ofs.open(str.c_str());
  if (!ofs.is_open()) {
    cerr << "\nCould not open data file " << argv[1] << "\n";
    exit(1);
  }


  // Retrieve results for 1000s step
	std::map<size_t, map<size_t, Error> >::iterator res;
	res = resultsMap.find(1);
	map<size_t, Error> tmpMap;
	if (res != resultsMap.end()) {

		tmpMap = res->second;

	} else {
		cerr << "big error";
	}



	for (map<size_t, Error>::iterator it = tmpMap.begin();
	    it != tmpMap.end(); it++) {

		for (vector<double>::iterator it2 = it->second.absError.begin();
		    it2 != it->second.absError.end();
		    it2++) {
			ofs << setprecision(12)
			    << *it2 << "\t";
		}

		ofs << "\n";

		cout << it->second.rmsError << "\t" << it->second.maxError << endl;

  }

	ofs.close();

  return 0;
	}

inline size_t getNodeNumber(size_t elementNumber) {
	return elementNumber + 1;
}


