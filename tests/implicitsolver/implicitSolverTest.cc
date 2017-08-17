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
	if (argc < 2) {
    cout << "solverTest called without argument. Error.\n";
    cerr << "solverTest called without argument. Error.\n";
    exit(1);
	} else if (argc > 3) {
    cout << "solverTest called with too many arguments. Error.\n";
    cerr << "solverTest called with too many arguments. Error.\n";
    exit(1);
  }

  /* Parameters definition */
	const double totalThickness = 10.0;
  const double T0 = 273.15;
	const double T1 = 263.15;
	const double snowDensity = 100.0;
	const double thetaAir = 1.0 / Constants::conductivity_air;
  const size_t nTestLayers = 10000;

	/* Read-in argument */
	size_t overallTimeStep = atoi(argv[1]);
	cout << "Running consistency checks for implicit integration scheme." << endl
	     << "Overall step length is " << overallTimeStep << "s." << endl;

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
	size_t layersVectorTmp[] = { 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000 };
	size_t timeVectorTmp[] = { 1000, 500, 200, 100, 50, 20, 10, 5, 2, 1 };
	//size_t timeVectorTmp[] = { 1 };

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
					if (p_dt > static_cast<double>(overallTimeStep) - Constants::eps) {  // Check if it is the last sub-time step
						LastTimeStep = true;
					}
				}
			}

			// Store the error in first map for layer "it"
			Error actualError;
			he.computeErrors(overallTimeStep, *it, Xdata.Ndata, actualError);
			tmpResultsMap.insert(
			    tmpResultsMap.end(),
			                     pair<size_t, Error>(
			        *it, actualError));

			}

		// Store error map in main map
		resultsMap.insert(resultsMap.begin(),
		                  pair<size_t, map<size_t, Error> >(*it2, tmpResultsMap));

		}




  /* Output to file */
  // Open input stream
  string str = "myoutput.dat";
  std::ofstream ofs;
  ofs.open(str.c_str());
  if (!ofs.is_open()) {
		cerr << "\nCould not open data file " << argv[0] << "\n";
    exit(1);
  }

	/* Physical dimension consistency check */
	cout << "Physical dimension consistency check" << endl;
	const size_t timeStep = 1;
	cout << "Integration time step is " << timeStep
	     << "s and overall time step is " << overallTimeStep << "s." << endl;
	cout << "| #layers [-] \t| RMS error [K] \t| MAX error [K] \t|" << endl;
  // Retrieve results for 1s step to check consistency
	std::map<size_t, map<size_t, Error> >::iterator res;
	res = resultsMap.find(timeStep);  // We choose the highest time resolution to ensure it is not the limiting factor
	map<size_t, Error> tmpMap;  // Map containing the results for 1s
	if (res != resultsMap.end()) {
		tmpMap = res->second;
	} else {
		cerr << "Could not retrieve result for integration time step " << timeStep
		     << "s. Please check timesteps vector." << endl;
		return 2;
	}


	// Huge start error
	double lastRmsError = 1e10;

	Matrix XTX(3, 3), XTY(3, 1);
	// Loop over grid sizes
	for (map<size_t, Error>::iterator it = tmpMap.begin();
	    it != tmpMap.end();
	    ++it) {

		// Output console
		cout << "| " << it->first << "\t| " << it->second.rmsError << "\t| "
		     << it->second.maxError << "\t|" << endl;

		// Check for consistency
		if (it->second.rmsError > lastRmsError) {
			cerr
			    << "CONSISTENCY ISSUE: RMS error does not monotically decrease with decreasing grid size.\nThe error for "
			    << (--it)->first
			    << " layers is smaller than for "
			    << (it)->first
			    << " layers. The time resolution might be not small enough.\n"
			    << "If this error appears again after a consequent time resolution improvement, then there is an actual consistency issue!"
			    << endl;
			it++;
			return 1;
		} else {
			lastRmsError = it->second.rmsError;
		}

		for (vector<double>::iterator it2 = it->second.absError.begin();
		    it2 != it->second.absError.end(); it2++) {
			ofs << setprecision(12) << *it2 << "\t";
		}

		ofs << "\n";



		// Fitting
		XTX(1, 1) += pow(it->first, 4.0);
		XTX(1, 2) += pow(it->first, 3.0);
		XTX(1, 3) += pow(it->first, 2.0);
		XTX(2, 1) += pow(it->first, 3.0);
		XTX(2, 2) += pow(it->first, 2.0);
		XTX(2, 3) += pow(it->first, 1.0);
		XTX(3, 1) += pow(it->first, 2.0);
		XTX(3, 2) += pow(it->first, 1.0);
		XTX(3, 3) += 1.0;

		XTY(1, 1) += pow(it->first, 2.0) * it->second.rmsError;
		XTY(2, 1) += it->first * it->second.rmsError;
		XTY(3, 1) += it->second.rmsError;


  }


	//Matrix coeffVec = Matrix::solve(XTX, XTY);

	ofs.close();

	/* END of physical dimension consistency check */

	/* Time dimension consistency check */
	cout << "Time dimension consistency check" << endl;

	// Retrieve results for finest grid resolution to check consistency
	const size_t nLayers = 2000;
	cout << "Number of layers is " << nLayers << " and overall time step is "
	     << overallTimeStep << "s." << endl;
	cout << "| step length [s] \t| RMS error [K] \t| MAX error [K] \t|" << endl;

	// Huge start error
	lastRmsError = 1e10;

	// Loop over step sizes, from lowest to finest resolution
	for (map<size_t, map<size_t, Error> >::reverse_iterator it =
	    resultsMap.rbegin(); it != resultsMap.rend(); it++) {

		// Find finest grid resolution for current time step and save Error
		map<size_t, Error>::iterator res;
		res = it->second.find(nLayers);
		Error tmpError;
		if (res != it->second.end()) {
			tmpError = res->second;
		} else {
			cerr << "Could not retrieve result for " << nLayers
			     << " layers. Please check layers vector." << endl;
			return 2;
		}

		// Output console
		cout << "| " << it->first << "\t| " << tmpError.rmsError << "\t| "
		     << tmpError.maxError << "\t|" << endl;

		// Check for consistency
		if (tmpError.rmsError > lastRmsError) {
				cerr
			    << "CONSISTENCY ISSUE: RMS error does not monotically decrease with decreasing time step.\nThe error for "
			    << (++it)->first
			    << "s per step is larger than for "
			    << (--it)->first
			    << "s per step. The grid resolution might be not small enough.\n"
			    << "If this error appears again after a consequent grid resolution improvement, then there is an actual consistency issue!"
			    << endl;
			return 1;
		} else {
			lastRmsError = tmpError.rmsError;
		}


	}









  return 0;
	}

inline size_t getNodeNumber(size_t elementNumber) {
	return elementNumber + 1;
}


