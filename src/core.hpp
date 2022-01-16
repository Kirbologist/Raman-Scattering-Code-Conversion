#ifndef SMARTIES_CORE_H
#define SMARTIES_CORE_H

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <memory>
#include <Eigen/Core>
#include <Eigen/CXX11/Tensor>
#include "tensor_matrix_cast.hpp"

using namespace std;
using namespace Eigen;

namespace Smarties {
  using RowArrayXi = Array<int, 1, Dynamic>;

  using ArrayXb = Array<bool, Dynamic, 1>;

  using RowArrayXb = Array<bool, 1, Dynamic>;

  using RowArrayXd = Array<double, 1, Dynamic>;

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
  using MatrixXr = Matrix<Real, Dynamic, Dynamic>;

  template <class Real>
  using MatrixXc = Matrix<complex<Real>, Dynamic, Dynamic>;

  template <class Real>
  using VectorXr = Matrix<Real, Dynamic, 1>;

  template <class Real>
  using VectorXc = Matrix<complex<Real>, Dynamic, 1>;

  using Tensor3d = Tensor<double, 3>;

  template <class Real>
  using Tensor3r = Tensor<Real, 3>;

  template <class Real>
  using Tensor3c = Tensor<complex<Real>, 3>;

  template <class Real>
  using Tensor4c = Tensor<complex<Real>, 4>;

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

    stIncPar() {};

    stIncPar(const stIncPar<Real>& base) {
      this->type = base.type;
      this->theta_p = base.theta_p;
      this->phi_p = base.phi_p;
      this->alpha_p = base.alpha_p;
      this->abs_m_vec = base.abs_m_vec;
    }
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
    ArrayXr<Real> lambda;
    Real epsilon1;
    ArrayXr<Real> epsilon2;
  };

  struct stOptions {
    bool get_R = false;
    int delta = 0;
    int NB = 0;
    ArrayXi abs_m_vec;
    bool get_symmetric_T = false;
    bool output = true;

    stOptions() {};

    stOptions(int N) {
      this->abs_m_vec = ArrayXi::LinSpaced(N + 1, 0, N);
    }
  };
}

#endif
