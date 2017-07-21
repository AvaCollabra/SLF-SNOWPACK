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

static string mode = "RESEARCH";

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
  const double totalThickness = 1.0;
  const double T0 = 273.15;
  const double T1 = 268.15;
  const double snowDensity = 100.0;
  const double thetaAir = 1.0 / Constants::conductivity_air;
  const size_t nTestLayers = 10000;
  const double layerTestThickness(totalThickness / nTestLayers);

  /* Needed variables */
  double layerThickness(0.0);

  /* Benchmark data from analytical solution */
  HeatEquationAnalytical he(
      totalThickness,
      T0,
      T1,
      Constants::density_air / Constants::conductivity_air
          * Constants::specific_heat_air / snowDensity,
      snowDensity, 1.0);

  vector<double> benchmarkTemp;
  for (size_t n = 0; n <= nTestLayers; n++) {

    benchmarkTemp.push_back(
        he.getTemperature((double) n * layerTestThickness, 900));

  }


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
  vector<size_t> layersVector(
      layersVectorTmp,
      layersVectorTmp + sizeof(layersVectorTmp) / sizeof(size_t));
  map<unsigned int, vector<NodeData> > nodesMap;

  /* Loop to fill container with start data and solve diffusion problem */
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

    // Compute temp profile
    testsnowpack.testCompTempProfile(Mdata, Xdata, Bdata, false);

    // Save result to map
    nodesMap.insert(nodesMap.end(),
                    pair<size_t, vector<NodeData> >(*it, Xdata.Ndata));

  }




  /* Output to file
  // Open input stream
  string str = "myoutput.dat";
  std::ofstream ofs;
  ofs.open(str.c_str());
  if (!ofs.is_open()) {
    cerr << "\nCould not open data file " << argv[1] << "\n";
    exit(1);
  }


  for (size_t n = 0; n < nN; n++) {
    ofs << setprecision(12) << NDS[n].T << "\t";
  }

  ofs << "\n";

  for (size_t n = 0; n < nN; n++) {
    ofs << setprecision(12)
        << he.getTemperature((double) n * layerThickness, 900)
        << "\t";
  }

   ofs << "\n";*/


  return 0;
}

inline size_t getNodeNumber(size_t elementNumber) {
  return elementNumber + 1;
}
