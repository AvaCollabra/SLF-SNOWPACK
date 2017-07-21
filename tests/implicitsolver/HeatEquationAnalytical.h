/*
 * HeatEquationAnalytical.h
 *
 *  Created on: Jul 7, 2017
 *      Author: julien
 */

#ifndef TESTS_IMPLICITSOLVER_HEATEQUATIONANALYTICAL_H_
#define TESTS_IMPLICITSOLVER_HEATEQUATIONANALYTICAL_H_

#include <math.h>
#include <iostream>
#include <stdlib.h>

class HeatEquationAnalytical {
 public:
  explicit HeatEquationAnalytical(const double & totalThickness,
                                  const double &T0,
                         const double &T1, const double &heatCapacity,
                                  const double &snowDensity = 400.0,
                                  const double &thermalConductivity = 1.0,
                                  const unsigned int &p = 100);

  virtual ~HeatEquationAnalytical();

  double getTemperature(const double &x, const double &t);

 private:

  inline void updateAlpha();

  double _T0, _T1;

  double _heatCapacity;
  double _snowDensity;
  double _thermalConductivity;


  unsigned int _p;

  double _alpha;
  double _totalThickness;

};

#endif /* TESTS_IMPLICITSOLVER_HEATEQUATIONANALYTICAL_H_ */
