/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Based on the SMARTIES MATLAB package by W.R.C. Somerville, B. Auguié, E.C. Le Ru
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.


This code contains all common or basic header files, typedefs, enums and structs used by the rest
of the SMARTIES code. Typedefs declared here are based on existing Eigen typedefs.
The -r suffix or 'Real' class usually denotes a scalar type that's used to represent real numbers.
The -c suffix denotes that each element is the complex extension of 'Real'.
*/

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

using namespace std;
using namespace Eigen;

namespace Smarties {
  using RowArrayXi = Array<int, 1, Dynamic>;

  using ArrayXb = Array<bool, Dynamic, 1>;

  using RowArrayXb = Array<bool, 1, Dynamic>;

  using RowArrayXf = Array<float, 1, Dynamic>;

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

  template <class Real>
  using Tensor3r = Tensor<Real, 3>;

  template <class Real>
  using Tensor3c = Tensor<complex<Real>, 3>;

  template <class Real>
  using Tensor4c = Tensor<complex<Real>, 4>;

  template <class Real>
  using Tensor1c = Tensor<complex<Real>, 1>;

  /*
  Used to describe the incident type of a photon in an electric field. The first x,y,z refers to
  wavevector direction and the second x,y,z (necessarily different to the first) is the direction of the
  linearly-polarised electric field. Otherwise, GENERAL if incident type can't be defined as such.
  */
  enum sIncType {Kx_Ez, Kx_Ey, Ky_Ex, Ky_Ez, Kz_Ex, Kz_Ey, GENERAL};

  /* Struct with parameters of an incident electric field */
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

  /*
  Various parameters used throughout SMARTIES functions to describe
  incident spheroids, electric fields and T matrices
  */
  template <class Real>
  struct stParams {
    Real a; // Spheroid semi-axis radius along x,y.
    Real c; // Spheroid semi-axis radius along z.
    ArrayXr<Real> k1; // Wavevector in embedding medium (possibly a wavelength-dependent vector)
    ArrayXr<Real> s; // Relative-refractive index
    int N; // The number of multipoles requested for the T (and R) matrix. The higher it is, the better the convergence
    int Nb_theta; // Number of angles used in Gaussian quadratures for the evaluation of P and Q matrix integrals
    sIncType inc_type;
    unique_ptr<stIncPar<Real>> inc_par;
    int Nb_theta_pst; // Number of angles for post-processing (should typically be larger than Nb_theta)
    ArrayXr<Real> lambda; // Wavelength (in free space. Same unit as a and c)
    Real epsilon1; // Relative dielectric constant of embedding medium (real positive)
    ArrayXr<Real> epsilon2; // Dielectric function of scatterer (may be complex, but that hasn't been implemented)
    bool output = true; // Originally not in the specification of User Guide, but some functions still use this anyway
  };

  /* Optional settings */
  struct stOptions {
    bool get_R = false; // Calculate R matrices and internal field coefficients
    int delta = 0; // Number of extra multipoles for P and Q matrices. If -1, sphEstimateDelta gets called
    int NB = 0; // Number of multipoles to compute the Bessel functions (NB >= N + delta)
    ArrayXi abs_m_vec; // Vector of values of abs(m) for which the T matrix is to be computed (0 <= abs(m) <= N)
    bool get_symmetric_T = false; // If true, T matrix is symmetrised
    bool output = true; // If true, prints some output to the terminal

    stOptions() {};

    stOptions(int N) {
      this->abs_m_vec = ArrayXi::LinSpaced(N + 1, 0, N); // Values of m in most cases of interest
    }
  };
}

#endif
