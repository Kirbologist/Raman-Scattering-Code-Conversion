#ifndef RAMAN_ELASTC_SCATTERING_H
#define RAMAN_ELASTC_SCATTERING_H

#include <cmath>
#include <complex>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/CXX11/Tensor>

using namespace std;
using namespace Eigen;

namespace Raman {
  using RowArrayXi = Array<int, 1, Dynamic>;

  using RowArrayXb = Array<bool, 1, Dynamic>;

  template <class Real>
  using ArrayXXr = Array<Real, Dynamic, Dynamic>;

  template <class Real>
  using ArrayXr = Array<Real, Dynamic, 1>;

  template <class Real>
  using RowArrayXr = Array<Real, 1, Dynamic>;

  template <class Real>
  using ArrayXXc = Array<complex<Real>, Dynamic, Dynamic>;

  template <class Real>
  using ArrayXc = Array<complex<Real>, Dynamic, 1>;

  template <class Real>
  using RowArrayXc = Array<complex<Real>, 1, Dynamic>;

  template <class Real>
  using VectorXc = Matrix<complex<Real>, Dynamic, 1>;

  template <class Real>
  struct stParams {
    ArrayXr<Real> s;
    ArrayXr<Real> k1;
  };
}

#endif
