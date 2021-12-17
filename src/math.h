#ifndef MY_MATH_H
#define MY_MATH_H

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
  Tensor3c<Real> subtensor(Tensor3c<Real>& tensor,
      ArithmeticSequence<long int, long int, long int> slice_dim1,
      ArithmeticSequence<long int, long int, long int> slice_dim2,
      ArithmeticSequence<long int, long int, long int> slice_dim3);

  template<class Real>
  ArrayXr<Real> arr_bessel_j(ArrayXr<Real>& nu, Real x);

  template<class Real>
  ArrayXr<Real> arr_bessel_y(ArrayXr<Real>& nu, Real x);

  const double PI = mp_pi<double>(), EPS = mp_eps<double>();
  const complex<double> I = mp_im_unit<double>();
}

#endif
