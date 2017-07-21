/*
 * HeatEquationAnalytical.cpp
 *
 *  Created on: Jul 7, 2017
 *      Author: julien
 */

#include "HeatEquationAnalytical.h"

HeatEquationAnalytical::HeatEquationAnalytical(
    const double & totalThickness,
    const double &T0,
    const double &T1, const double &heatCapacity, const double &snowDensity,
    const double &thermalConductivity,
    const unsigned int &p)
    : _T0(T0),
      _T1(T1),
      _heatCapacity(heatCapacity),
      _snowDensity(snowDensity),
      _thermalConductivity(thermalConductivity),
      _p(p),
      _alpha(0.),
      _totalThickness(totalThickness) {

  updateAlpha();




}

HeatEquationAnalytical::~HeatEquationAnalytical() {
  // TODO Auto-generated destructor stub
}

double HeatEquationAnalytical::getTemperature(const double& x,
                                              const double& t) {

  if (t < 0) {
    std::cerr << "t must be greater than 0" << std::endl;
    exit(1);
  }

  if (x < 0 || x > _totalThickness) {
    std::cerr
        << "x shall not me smaller than 0 or bigger than the total thickness "
        << _totalThickness << std::endl;
    exit(1);
  }

  double result = 0.0;

  for (unsigned int p = _p; p > 0; p--) {  // Going backward for numerical reasons...

    double tmp = 2 / M_PI * (_T1 - _T0) / p * (p % 2 ? -1.0 : 1.0);

    tmp *= exp(
        -_alpha * M_PI * M_PI / _totalThickness / _totalThickness * p * p * t);

    tmp *= sin(p * M_PI / _totalThickness * x);

    result += tmp;
  }

  result += (_T1 - _T0) * x / _totalThickness + _T0;

  return result;

}

inline void HeatEquationAnalytical::updateAlpha() {
  _alpha = _thermalConductivity / (_snowDensity * _heatCapacity);
}

