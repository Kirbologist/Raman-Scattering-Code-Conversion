#ifndef RAMAN_ELASTC_SCATTERING_H
#define RAMAN_ELASTC_SCATTERING_H

#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <memory>
#include <Eigen/Core>
#include <Eigen/CXX11/Tensor>
#include "TensorMatrixCast.h"

using namespace std;
using namespace Eigen;

namespace Raman {
  using RowArrayXi = Array<int, 1, Dynamic>;

  using ArrayXb = Array<bool, Dynamic, 1>;

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
  using MatrixXc = Matrix<complex<Real>, Dynamic, Dynamic>;

  template <class Real>
  using VectorXr = Matrix<Real, Dynamic, 1>;

  template <class Real>
  using VectorXc = Matrix<complex<Real>, Dynamic, 1>;

  template <class Real>
  using Tensor3c = Tensor<complex<Real>, 3>;

  template <class Real>
  using Tensor1c = Tensor<complex<Real>, 1>;

  enum sIncType {KxEz, KxEy, KyEx, KyEz, KzEx, KzEy, GENERAL};

  template<class Real>
  struct stIncPar {
    sIncType type;
    Real theta_p;
    Real phi_p;
    Real alpha_p;
    ArrayXi abs_m_vec;
  };

  template <class Real>
  struct stParams {
    Real a;
    Real c;
    ArrayXr<Real> k1;
    ArrayXr<Real> s;
    int N;
    int Nb_theta;
    sIncType inc_type;
    unique_ptr<stIncPar<Real>> inc_par;
    int Nb_theta_pst;
    Real lambda;
    Real epsilon1;
    complex<Real> epsilon2;
  };

  template <class Real>
  struct stOptions {
    bool get_R = false;
    int delta = 0;
    int NB = 0;
    ArrayXi abs_m_vec;
    bool get_symmetric_T = false;
    bool output = true;

    stOptions(int N) {
      this.abs_m_vec = ArrayXi::LinSpaced(N + 1, 0, N);
    }
  };
}

#endif
