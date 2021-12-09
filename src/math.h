#ifndef MATH_H
#define MATH_H

#include "raman_elastic_scattering.h"
#include <boost/math/special_functions/bessel.hpp>

using namespace std;
using namespace Eigen;

namespace Raman {
  template <class Real>
  Real mp_pi(void);

  template<class Real>
  Real mp_eps(void);

  template<class Real>
  complex<Real> mp_im_unit(void);

  template<class Real>
  RowArrayXr<Real> arr_bessel_j(RowArrayXr<Real>& nu, Real x);

  template<class Real>
  RowArrayXr<Real> arr_bessel_y(RowArrayXr<Real>& nu, Real x);

  const double PI = mp_pi<double>(), EPS = mp_eps<double>();
  const complex<double> I = mp_im_unit<double>();
}

#endif
