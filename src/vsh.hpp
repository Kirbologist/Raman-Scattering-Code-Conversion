/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Based on the SMARTIES MATLAB package by W.R.C. Somerville, B. Auguié, E.C. Le Ru
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.


This code contains all 'vsh' SMARTIES functions that are used in Raman scattering calculations,
i.e. mostly low-level functions handling the calculations of quantites related to vector spherical wave functions.
We often denote these vector spherical wave functions, or vector spherical harmonics as VSHs.
*/

#ifndef VSH_HPP
#define VSH_HPP

#include "core.hpp"
#include "misc.hpp"

namespace Smarties {

  enum sBessel {J, H1};

  /* Contains angular functions pi and tau as defined in Mishchenko 2002, and the Legendre polynomials d_{n0}*/
  template <class Real>
  struct stPinmTaunm {
    ArrayXXr<Real> pi_nm;
    ArrayXXr<Real> tau_nm;
    ArrayXXr<Real> p_n0; // Legendre polynomial with d_n0(theta) = P_n(cos(theta))
  };

  /* Contains the expansion coefficients of an incident plane wave */
  template <class Real>
  struct stIncEabnm {
    ArrayXc<Real> a_nm;
    ArrayXc<Real> b_nm;
  };

  /* Contains information about the three Z_n(rho) auxiliary functions for the radial dependence of VSHs */
  template <class Real>
  struct stZnAll {
    ArrayXXc<Real> Z0; // Z_n^0(rho) = z_n(rho)
    ArrayXXc<Real> Z1; // Z_n^1(rho) = z_n(rho)/rho
    ArrayXXc<Real> Z2; // Z_n^2(rho) = [rho*z_n(rho)]'/rho
  };

  template <class Real>
  struct stEAllPhi {
    RowArrayXr<Real> theta;
    RowArrayXr<Real> r_of_theta;
    vector<ArrayXXc<Real>> E_rm;
    vector<ArrayXXc<Real>> E_tm;
    vector<ArrayXXc<Real>> E_fm;
  };

  template <class Real>
  struct stEforPhi {
    ArrayXXc<Real> E_r;
    ArrayXXc<Real> E_t;
    ArrayXXc<Real> E_f;
    ArrayXr<Real> theta;
    Real phi0;
  };

  /*
  Create struct containing parameters of the incident electric field for plane wave excitation.
  Here, the parameters are specified using the enum type sIncType.
  The function name is a shortening of the function vshMakeIncidentParameters from the MATLAB code.
  Inputs:
    type - enumerated constant of the form `Ka_Eb`, where a and b are one of x, y or z and a != b.
           a refers to the wavevector direction, and b is the direction of the linearly-polarised electric field.
    N_max - the maximum order used in the expansions, required to determine the values of m required
  Output:
    Returns a unique pointer to a struct containing the incident parameters theta_p, phi_p, alpha_p and abs_m_vec
  Dependencies:
    mp_pi
  */
  template <class Real>
  unique_ptr<stIncPar<Real>> vshMakeIncidentParams(sIncType type, int N_max) {
    auto output = make_unique<stIncPar<Real>>();
    switch (type) {
      case Kx_Ez: {
        output->type = Kx_Ez;
        output->theta_p = mp_pi<Real>()/2;
        output->phi_p = 0;
        output->alpha_p = mp_pi<Real>();
        output->abs_m_vec = ArrayXi::LinSpaced(N_max + 1, 0, N_max);
        break;
      }

      case Kx_Ey: {
        output->type = Kx_Ey;
        output->theta_p = mp_pi<Real>()/2;
        output->phi_p = 0;
        output->alpha_p = mp_pi<Real>()/2;
        output->abs_m_vec = ArrayXi::LinSpaced(N_max + 1, 0, N_max);
        break;
      }

      case Ky_Ez: {
        output->type = Ky_Ez;
        output->theta_p = mp_pi<Real>()/2;
        output->phi_p = mp_pi<Real>()/2;
        output->alpha_p = mp_pi<Real>();
        output->abs_m_vec = ArrayXi::LinSpaced(N_max + 1, 0, N_max);
        break;
      }

      case Ky_Ex: {
        output->type = Ky_Ex;
        output->theta_p = mp_pi<Real>()/2;
        output->phi_p = mp_pi<Real>()/2;
        output->alpha_p = -mp_pi<Real>()/2;
        output->abs_m_vec = ArrayXi::LinSpaced(N_max + 1, 0, N_max);
        break;
      }

      case Kz_Ex: {
        output->type = Kz_Ex;
        output->theta_p = 0;
        output->phi_p = 0;
        output->alpha_p = 0;
        output->abs_m_vec = {1};
        break;
      }

      case Kz_Ey: {
        output->type = Kz_Ey;
        output->theta_p = 0;
        output->phi_p = 0;
        output->alpha_p = mp_pi<Real>()/2;
        output->abs_m_vec = {1};
        break;
      }

      default:
        output->type = Kz_Ex;
        output->theta_p = 0;
        output->phi_p = 0;
        output->alpha_p = 0;
        output->abs_m_vec = {1};
    }
    return output;
  }

  /*
  Create struct containing parameters of the incident electric field for plane wave excitation.
  Here, the parameters are the angles of the wavevector, and the electric field.
  The function name is a shortening of the function vshMakeIncidentParameters from the MATLAB code.
  Inputs:
    type - enumerated constant of the form `Ka_Eb`, where a and b are one of x, y or z and a != b.
           a refers to the wavevector direction, and b is the direction of the linearly-polarised electric field.
           If `type` is not of this form, the electric field parameters are retrieved from the other input parameters.
    N_max - the maximum order used in the expansions, required to determine the values of m required
    theta_p - this, along with phi_p, defines the direction of the incident wavevector k
    phi_p - this, along with theta_p, defines the direction of the incident wavevector k
    alpha_p - the defines the orientation of the electric field, in the plane orthogonal to wavevector k
  Output:
    Returns a unique pointer to a struct containing the incident parameters theta_p, phi_p, alpha_p and abs_m_vec
  */
  template <class Real>
  unique_ptr<stIncPar<Real>> vshMakeIncidentParams(
      sIncType type, int N_max, Real theta_p, Real phi_p, Real alpha_p) {
    unique_ptr<stIncPar<Real>> output;
    if (type == GENERAL) {
      output = make_unique<stIncPar<Real>>();
      output->type = GENERAL;
      output->theta_p = theta_p;
      output->phi_p = phi_p;
      output->alpha_p = alpha_p;
      output->abs_m_vec = ArrayXi::LinSpaced(N_max + 1, 0, N_max);
    } else {
      output = vshMakeIncidentParams<Real>(type, N_max);
    }
    return output;
  }

  /*
  Calculates angular functions pi_nm(theta) and tau_nm(theta) for n=0..N_max, |m|<=n, using recurrence relations.
  These functions are defined on p. 373 of Mishchenko 2002. It also returns the Legendre polynomials
  d_{n0}(theta) = P_n(cos(theta)), which are necessary for m = 0. Note that the case where n,m = 0 is not needed
  usually, but is included as padding to make the indexing simpler
  Inputs:
    N_max - maximum value of n
    theta - [T x 1] with theta angles (in radians). All thetas must be between 0 and pi.
  Output:
    Returns struct containing the functions pi_nm, tau_nm in an Eigen::Array of size [T x P] and
    the Legendre polynomials p_n0 of size [T x N_max] with d_n0(theta) = P_n(cos(theta)) where P = (N_max + 1)^2.
    The arrays pi_nm and tau_nm use the p-index, which stores the possible values of (n, m) in a linear array
    using the convention p = n*(n+1) + m.
  */
  template <class Real>
  unique_ptr<stPinmTaunm<Real>> vshPinmTaunm(int N_max, const ArrayXr<Real>& theta) {
    if ((theta < 0.0).any())
      cerr << "Warning: theta must be >= 0 in vshPinmTaunm..." << endl;
    int n_rows = size(theta);
    int P_max = (N_max + 1)*(N_max + 1); // maximum number of columns for pi_n and tau_n matrices
    auto output = make_unique<stPinmTaunm<Real>>();
    output->pi_nm = ArrayXXr<Real>::Zero(n_rows, P_max);
    output->tau_nm = ArrayXXr<Real>::Zero(n_rows, P_max);

    // Initialise Am for m = 0
    ArrayXr<Real> A_m_sin_mm1 = ArrayXr<Real>::Ones(n_rows); // [T x 1]
    ArrayXr<Real> mu_c = cos(theta); // [T x 1]
    ArrayXr<Real> mu_s = sin(theta); // [T x 1]

    // loop on m (m = 0 case is treated separately)
    for (int m = 1; m <= N_max; m++) {
      // Am * sin(theta)^(m - 1) is computed by recurrence on m
      A_m_sin_mm1 *= sqrt(static_cast<Real>(2*m - 1)/(2*m))*
          (m > 1 ? mu_s : ArrayXr<Real>::Ones(n_rows));
      int n_cols = N_max - m + 2;
      ArrayXXr<Real> pi_aux(n_rows, n_cols);
      pi_aux.col(0).setZero();
      pi_aux.col(1) = m*A_m_sin_mm1;
      // pi_aux contains pi_{m+j-2, m}, j = 1..(N_max-m+2), i.e. pi_{m-1, m}..pi_{N_max, m}

      // Get pi_{m+1, m} to pi_{N_max} by recurrence
      for (int j = 2, n = m + 1; j < n_cols; j++, n++)
        // pi_aux.col(j) is pi_{m+j-2, m}
        pi_aux.col(j) = (1/sqrt(static_cast<Real>((n - m) * (n + m)))) * ((2*n - 1)*mu_c*
            pi_aux.col(j - 1) - sqrt(static_cast<Real>((n - 1 - m)*(n - 1 + m))) * pi_aux.col(j - 2));

      ArrayXi n_vec = ArrayXi::LinSpaced(n_cols - 1, m, N_max); // [N2 x 1] with N2 = n_cols - 1
      ArrayXr<Real> n_vec_real = n_vec.template cast<Real>();
      ArrayXi p_vec = n_vec*(n_vec + 1) + m; // computes p for positive m
      ArrayXi p_vec_n = p_vec - 2*m; // computes p for negative m

      for (int n = 0; n < n_cols - 1; n++) {
        // fill in pi_nm matrix for positive or negative m
        (output->pi_nm).col(p_vec(n)) = pi_aux.col(n + 1);
        (output->pi_nm).col(p_vec_n(n)) = pow(-1, (m + 1) % 2) * pi_aux.col(n + 1);

        // fill in tau_nm matrix for positive or negative m
        (output->tau_nm).col(p_vec(n)) = pi_aux.col(n) *
            (-sqrt((n_vec_real(n) - m) * (n_vec_real(n) + m))/m)
            + (mu_c*(n_vec_real(n)/m))*pi_aux.col(n + 1);
        (output->tau_nm).col(p_vec_n(n)) = pow(-1, m % 2)*(output->tau_nm).col(p_vec(n));
      }
    }

    // Now do m = 0 case. Initialise reccurence p_0 = 1, p_1 = mu_c, t_0 = 0, t_1 = -mu_s
    // p_nm1 contains p_n, n = 1..N_max, same for t_nm1
    ArrayXXr<Real> p_nm1(n_rows, N_max + 1); // [T x (N + 1)]
    ArrayXXr<Real> t_nm1(n_rows, N_max + 1); // [T x (N + 1)]
    p_nm1.col(0).setOnes(); // [T x 1]
    p_nm1.col(1) = mu_c; // [T x 1]
    t_nm1.col(0).setZero(); // [T x 1]
    t_nm1.col(1) = -mu_s; // [T x 1]

    // Get p_2 to p_{N_max} and t_2 to t_{N_max} by recurrence.
    for (int n = 2; n <= N_max; n++) {
      p_nm1.col(n) = static_cast<Real>(2*n - 1)/n*mu_c * p_nm1.col(n - 1) -
          static_cast<Real>(n - 1)/n*p_nm1.col(n - 2);
      t_nm1.col(n) = mu_c*t_nm1.col(n - 1) - (n)*mu_s*p_nm1.col(n - 1);
    }

    // Return p_n matrix (except n = 0)
    output->p_n0 = p_nm1;
    // fill in tau_nm for m = 0 (pi_nm = 0 for m = 0 is already set)
    ArrayXi n_vec = ArrayXi::LinSpaced(N_max + 1, 0, N_max);
    ArrayXi p_vec = n_vec*(n_vec + 1);

    for (int n = 0; n <= N_max; n++) {
      output->pi_nm.col(p_vec(n)).setZero();
      output->tau_nm.col(p_vec(n)) = t_nm1.col(n);
    }

    return output;
  }

  /*
  Calculates expansion coefficients for a general linearly polarised incident plane wave (assuming E_0 = 1)
  for all m up to multipole N_max. The expressions are given in eqs. C.56-C.59 of [Mishchenko 2002].
  The name of this function is shortened from vshGetIncidentCoefficients from the MATLAB code.
  Inputs:
    N_max - maximum multipole order
    angles - unique pointer to struct with three angles defining the incident wave.
             theta_p and phi_p define the wavevector direction, while alpha_p defines the field polarisation.
  Output:
    A struct with member variables a_nm and b_nm, each [P x 1] with the expansion coefficients for all
    0<=n<=N_max, |m|<=n. While not usually not needed, the (n,m) = (0,0) case is included as padding to make
    the indexing much simpler. The indices (n,m) are condensed in a linear index (the p-index) with
    the convention p = n*(n+1) + m. Note that P = (N_max + 1)^2.
  Dependencies:
    mp_pi, mp_im_unit, vshPinmTaunm
  */
  template <class Real>
  unique_ptr<stIncEabnm<Real>> vshGetIncidentCoeffs(int N_max, const unique_ptr<stIncPar<Real>>& angles) {
    complex<Real> I = mp_im_unit<Real>();
    Real alpha_p = angles->alpha_p;
    Real phi_p = angles->phi_p;
    Array<Real, 1, 1> theta_p = {{angles->theta_p}};
    int P_max = (N_max + 1)*(N_max + 1); // for p-index

    // First get d_bar_nm (see SMARTIES user guide)
    ArrayXr<Real> n_vec = ArrayXr<Real>::LinSpaced(N_max + 1, 0, N_max);
    // may cause a divide-by-zero exception depending on hardware
    ArrayXc<Real> fact_n = pow(I, n_vec) * sqrt(4*mp_pi<Real>()*(2*n_vec + 1) / (n_vec*(n_vec+1))); // [N x 1]
    fact_n(0) = 0;
    ArrayXr<Real> m_vec = ArrayXr<Real>::LinSpaced(2*N_max + 1, -N_max, N_max);
    ArrayXc<Real> fact_m = -pow(-1, m_vec)*exp(-I*m_vec*phi_p); // [M x 1]

    ArrayXc<Real> d_bar_nm(P_max); // [P x 1]
    for (int n = 1; n <= N_max; n++) { // loop on n
      ArrayXi m = ArrayXi::LinSpaced(2*n + 1, -n, n); // all m at a time [M x 1]
      ArrayXi ind = n*(n + 1) + m; // [M x 1]
      for (int j = 0; j < 2*n + 1; j++)
        d_bar_nm(ind(j)) = fact_n(n) * fact_m(m(j) + N_max);
    }

    // Now get E.B and E.C part
    unique_ptr<stPinmTaunm<Real>> stPTp = vshPinmTaunm<Real>(N_max, theta_p);
    ArrayXc<Real> minus_EC_nm_star = static_cast<complex<Real>>(cos(alpha_p))*I*stPTp->pi_nm.row(0)
        + sin(alpha_p)*stPTp->tau_nm.row(0);
    ArrayXc<Real> i_EB_nm_star = static_cast<complex<Real>>(cos(alpha_p))*I*stPTp->tau_nm.row(0)
        + sin(alpha_p)*stPTp->pi_nm.row(0);

    auto output = make_unique<stIncEabnm<Real>>();
    // Multiply to get final result (eqs. C.57 and C.58 of [Mishchenko 2002])
    output->a_nm = d_bar_nm * minus_EC_nm_star;
    output->b_nm = d_bar_nm * i_EB_nm_star;
    return output;
  }

  /*
  Computes the three Z_n(rho) auxiliary functions for the radial dependence of VSHsfor n=1..N_max.
  Can be used for both regular VSHs (based on j(rho)) or irregular VSHs (based on h1(rho)).
  Inputs:
    N_max - maximum multipole order
    rho - [R x 1] arguments of the VSHs (no zero components allowed, even for regular VSH, for speed optimisation)
    type - enum type defining the Bessel function to be used, so either `J` or `H1`.
  Output:
    Returns a unique pointer to a struct containing the Z_n functions.
  Dependencies:
    mp_pi, mp_im_unit, ArrBesselJ, ArrBesselY
  */
  template <class Real>
  unique_ptr<stZnAll<Real>> vshGetZnAll(int N_max, const ArrayXr<Real>& rho, sBessel type) {
    if ((rho == 0).any())
      cerr << "Warning: rho = 0 arguments not allowed in vshZnAll..." << endl;

    ArrayXr<Real> nu = ArrayXr<Real>::LinSpaced(N_max + 1, 0.5, N_max + 0.5);
    ArrayXXc<Real> f(rho.size(), N_max + 1);

    for (int i = 0; i < rho.size(); i++) {
      f.row(i) = ArrBesselJ(nu, rho(i));
      if ((f.row(i) == static_cast<complex<Real>>(0)).any()) {
        cerr << "Warning: Bessel (j) calculation went beyond precision in vshGetZnAll()" << endl;
        cerr << "x = " << rho(i) << " N_max = " << N_max << endl;
      }
    }

    if (type == H1) {
      for (int i = 0; i < rho.size(); i++) {
        ArrayXr<Real> y = ArrBesselY(nu, rho(i));
        if ((f.row(i).isInf()).any()) {
          cerr << "Warning: Bessel (y) calculation went beyond precision in vshGetZnAll()" << endl;
          cerr << "x = " << rho(i) << ", N_max = " << N_max << endl;
        }
        f.row(i) += mp_im_unit<Real>()*y;
      }
    }
    // f is array [R x N_max + 1] of cylindrical Bessel Z_n{n+0.5}(rho), n = 0..N_max

    f.colwise() *= sqrt((mp_pi<Real>()/2) / rho);
    // f is now array of spherical Bessel z_n(rho), n = 0..N_max, or equivalently z_{n-1}(rho), n = 1..N_max+1

    RowArrayXc<Real> n = RowArrayXr<Real>::LinSpaced(N_max, 1, N_max);
    auto output = make_unique<stZnAll<Real>>();
    output->Z0 = f; // [R x N_max]
    output->Z1 = (output->Z0).colwise() / rho.template cast<complex<Real>>(); // [R x N_max]

    // Computes: Z2_n = z_{n-1} - n*Z1_n
    output->Z2 = ArrayXXc<Real>::Zero(rho.size(), N_max + 1);
    output->Z2.rightCols(N_max) = f(all, seq(0, last - 1)) - (output->Z1.rightCols(N_max)).rowwise() * n;
    return output;
  }

  /*
  Calculates VSH expansions for r(theta) and many lambda using series expansions. The VSHs used can either be of
  type 1 (`type` = `J`) or type 3 (`type` = `H1`). Use rt(theta) = Inf (with `H1`) to obtain far-field properties.
  The fields E_rm, E_tm and E_fm given in the results are discussed in the supplementary information of SMARTIES.
  Inputs:
    lambda - [L x 1] wavelengths
    epsilon - [1 x 1] or [L x 1] dielectric function where field is evaluated
    p_nm, q_nm - 2 arrays [L x P] where P = (N_max + 1)^2 containing the coefficients p_{n,m} and q_{n,m}
                 for the field expansion of M^(i)_{n,m} and N^(i)_{n,m} respectively
    rt - column array [T x 1] or possibly zero (if `type` = `J`) spherical coordinate r (in nm) of points,
         for each corresponding theta. If values are Inf, then these correspond to far-field values.
    theta - [T x 1] spherical coordinate theta of points
    type - enum type defining the Bessel function to be used. `J` for regular VSH or `H1` for irregular
    stPT - unique pointer to a struct containing the functions pi_n(theta) and tau_n(theta).
           If omitted, then the functions are computed from scratch. It is faster to pass in this struct
           if these functions have already been calculated.
  Output:
    Returns unique pointer to a struct containing the 3 field components E_r, E_t = E_theta and E_f = E_phi.
    Each are [2N + 1] vectors containing unique pointers to an array A_m for each m = -N..N.
    A_m is an array [L x T] for the wavelength and theta dependence such that A = sum_{m=-N..N} A_m exp(i*m*phi).
  Dependencies:
    mp_pi, mp_im_unit, vshGetZnAll, vshPinmTaunm
  */
  template <class Real>
  unique_ptr<stEAllPhi<Real>> vshEgenThetaAllPhi(
      const ArrayXr<Real>& lambda, const ArrayXr<Real>& epsilon, const ArrayXXc<Real>& p_nm,
      const ArrayXXc<Real>& q_nm, const RowArrayXr<Real>& rt, const RowArrayXr<Real>& theta,
      sBessel type, unique_ptr<stPinmTaunm<Real>> stPT = unique_ptr<stPinmTaunm<Real>>()) {
    int P_max = p_nm.cols();
    int N_max = static_cast<int>(round(sqrt(P_max) - 1));
    int Nb_lambda = lambda.size();
    if (rt.size() != theta.size() && rt(0) != 0 && !isinf(rt(0)))
      cerr << "vshEgenThetaAllPhi error: theta and rt must be the same size row arrays." << endl;
    int Nb_theta = theta.size();

    ArrayXr<Real> n = ArrayXr<Real>::LinSpaced(N_max + 1, 0, N_max); // [N_max x 1]
    // need n-dep coeff for series mu)n * sqrt(n*(n+1)) and mu_n/sqrt(n*(n+1))
    ArrayXr<Real> mu_n_times = sqrt((2*n + 1)*n*(n + 1)/(4*mp_pi<Real>())); // mu_n [N_max x 1] for E_t and E_f
    ArrayXr<Real> mu_n_divd_gen = mu_n_times/(n*(n + 1));
    mu_n_divd_gen(0) = 0;

    auto output = make_unique<stEAllPhi<Real>>();
    output->theta = theta;
    output->r_of_theta = rt;
    output->E_rm = vector<ArrayXXc<Real>>(2*N_max + 1);
    output->E_tm = vector<ArrayXXc<Real>>(2*N_max + 1);
    output->E_fm = vector<ArrayXXc<Real>>(2*N_max + 1);

    if (!rt(0)) { // if rt(0) = 0, then all of rt should be zeros
      for (int m = -N_max; m <= N_max; m++) {
        if (abs(m) > 1) {
          output->E_rm[m + N_max] = ArrayXXc<Real>::Zero(Nb_lambda, Nb_lambda);
          output->E_tm[m + N_max] = ArrayXXc<Real>::Zero(Nb_lambda, Nb_lambda);
          output->E_fm[m + N_max] = ArrayXXc<Real>::Zero(Nb_lambda, Nb_lambda);
        }
      }
      // special case where r0 = 0
      Real coeff1 = 1/sqrt(6*mp_pi<Real>());
      Real coeff2 = coeff1/sqrt(static_cast<Real>(2.0));

      // Results are all [L x 1] x [1 x T] = [L x T] matrices
      // m = 0
      output->E_rm[N_max] = (coeff1 * q_nm.col(2)).matrix() * cos(theta).matrix();
      output->E_tm[N_max] = (-coeff1 * q_nm.col(2)).matrix() * sin(theta).matrix();
      output->E_fm[N_max] = ArrayXc<Real>::Zero(Nb_lambda, Nb_theta);
      // m = 1
      output->E_rm[N_max + 1] = (coeff2 * q_nm.col(3)).matrix() * sin(theta).matrix();
      output->E_tm[N_max + 1] = (coeff2 * q_nm.col(3)).matrix() * cos(theta).matrix();
      output->E_fm[N_max + 1] = (-mp_im_unit<Real>()*coeff2 * q_nm.col(3)).replicate(1, Nb_theta);
      // m = -1
      output->E_rm[N_max - 1] = (coeff2 * q_nm.col(1)).matrix() * sin(theta).matrix();
      output->E_tm[N_max - 1] = (coeff2 * q_nm.col(1)).matrix() * cos(theta).matrix();
      output->E_fm[N_max - 1] = (-mp_im_unit<Real>()*coeff2 * q_nm.col(1)).replicate(1, Nb_theta);

      cout << "r0 = 0 in vshEgenThetaAllPhi" << endl;
      return output;
    }

    // r != 0 from here
    // get Z_n(rho) for radial dependence and derived functions
    ArrayXr<Real> rho_col;
    unique_ptr<stZnAll<Real>> st_zn_all_col;
    if (!isinf(rt(0))) {
      // Matrix product of [L x 1] by [1 x T]
      ArrayXXr<Real> kr = (2*mp_pi<Real>()*sqrt(epsilon)/lambda).matrix() * rt.matrix(); // [L x T]
      // Column array [LT x 1] for all [L x T] arguments
      rho_col = kr.transpose().reshaped();
      st_zn_all_col = vshGetZnAll(N_max, rho_col, type); // [LT x N_max]
    } else { // Special case for radiation profile
      st_zn_all_col = make_unique<stZnAll<Real>>();
      st_zn_all_col->Z0 = st_zn_all_col->Z1 = st_zn_all_col->Z2 =
          ArrayXXc<Real>::Ones(Nb_lambda*Nb_theta, N_max + 1);
    }

    // Get thet dependence if not provided
    if (!stPT)
      stPT = vshPinmTaunm<Real>(N_max, theta.transpose()); // [T x N_max]

    // Loop over lambda. At a fixed lambda, the sum over n for all theta can be carried out using
    // a matrix product of the theta-and-n-dependent [T x N_max] matrix by a n-dependent column vector [N_max x 1].
    // The result, a [T x 1] matrix, is then transposed to a [1 x T] line.
    ArrayXXc<Real> E_r_sum(Nb_lambda, Nb_theta);
    ArrayXXc<Real> E_t_sum(Nb_lambda, Nb_theta);
    ArrayXXc<Real> E_f_sum(Nb_lambda, Nb_theta);
    for (int m = -N_max; m <= N_max; m++) {
      ArrayXi n_vec = ArrayXi::LinSpaced(N_max - abs(m) + 1, abs(m), N_max); // [N2 x 1]
      ArrayXi p_vec = n_vec*(n_vec + 1) + m; // [N2 x 1]
      ArrayXXr<Real> pi_nm = stPT->pi_nm(all, p_vec); // [T x N2]
      ArrayXXr<Real> tau_nm = stPT->tau_nm(all, p_vec); // [T x N2]
      ArrayXXr<Real> d_nm = m ? pi_nm.colwise()*(sin(theta)/m).transpose() : stPT->p_n0; // [T x N2]

      ArrayXXc<Real> ip_nm_for_Z0, q_nm_for_Z1;
      ArrayXc<Real> mu_n_divd;
      if (isinf(rt(0))) {
        // for far-field radiation profile
        q_nm_for_Z1 = ArrayXXc<Real>::Zero(Nb_lambda, Nb_theta); // [L x 1]
        ip_nm_for_Z0 = p_nm(all, p_vec); // [L x N2]
        mu_n_divd = mu_n_divd_gen*pow(-mp_im_unit<Real>(), n + 1); // [L x N2]
      } else {
        q_nm_for_Z1 = q_nm(all, p_vec); // [L x N2]
        ip_nm_for_Z0 = mp_im_unit<Real>()*p_nm(all, p_vec); // [L x N2]
        mu_n_divd = mu_n_divd_gen; // [L x N2]
      }
      ArrayXXc<Real> q_nm_for_Z2 = q_nm(all, p_vec);

      for (int l = 0; l < Nb_lambda; l++) {
        ArrayXi ind_in_rho_col = ArrayXi::LinSpaced(Nb_theta, 0, Nb_theta - 1) + l*Nb_theta;
        // for E_r, vec_N = d_{n,1} * Z_n^1(rho) * mu_n * n*(n+1)
        VectorXc<Real> vec_N_dep = q_nm_for_Z1.row(l) * mu_n_times(n_vec).transpose();
        // E_r_sum = sum_n(pi_n(t) * vec_N_n). Do the sum as a matrix product
        E_r_sum.row(l) = ((d_nm*st_zn_all_col->Z1(ind_in_rho_col, n_vec)).matrix() * vec_N_dep).transpose(); // [1 x T]
        // for E_t and E_f
        vec_N_dep = (ip_nm_for_Z0.row(l) * mu_n_divd(n_vec).transpose()).transpose(); // [N2 x 1]
        VectorXc<Real> vec_N_dep2 = (q_nm_for_Z2.row(l) * mu_n_divd(n_vec).transpose()).transpose(); // [N2 x 1]
        // Do the sums as matrix products
        MatrixXc<Real> tmp1 = (pi_nm * st_zn_all_col->Z0(ind_in_rho_col, n_vec)).matrix() * vec_N_dep;
        MatrixXc<Real> tmp2 = (tau_nm * st_zn_all_col->Z2(ind_in_rho_col, n_vec)).matrix() * vec_N_dep2;
        E_t_sum.row(l) = (tmp1 + tmp2).transpose(); // [1 x T]

        tmp1 = (tau_nm * st_zn_all_col->Z0(ind_in_rho_col, n_vec)).matrix() * vec_N_dep;
        tmp2 = (pi_nm * st_zn_all_col->Z2(ind_in_rho_col, n_vec)).matrix() * vec_N_dep2;
        E_f_sum.row(l) = (tmp1 + tmp2).transpose(); // [1 x T]
      }

      output->E_rm[m + N_max] = ArrayXXc<Real>(pow(-1, m) * E_r_sum);
      output->E_tm[m + N_max] = ArrayXXc<Real>(pow(-1, m) * E_t_sum);
      output->E_fm[m + N_max] = ArrayXXc<Real>(mp_im_unit<Real>() * static_cast<complex<Real>>(pow(-1, m)) * E_f_sum);
    }

    return output;
  }

  /*
  Calculate the E field for all lambdas and thetas in a plane of some constant phi.
  Inputs:
    st_E_surf - unique pointer to struct containing field
    phi0 - value of phi at which the field is to be calculated
  Output:
    Returns unique pointer to struct containing the components of the E field, each [L x T],
    and containing the values of theta [T x 1] and phi0.
  Dependencies:
    mp_im_unit
  */
  template <class Real>
  unique_ptr<stEforPhi<Real>> vshEthetaForPhi(const unique_ptr<stEAllPhi<Real>>& st_E_surf, Real phi0) {
    auto output = make_unique<stEforPhi<Real>>();
    int N_max = (st_E_surf->E_rm.size() - 1) / 2;
    output->theta = st_E_surf->theta;
    output->phi0 = phi0;

    int Nb_lambda = st_E_surf->E_rm[0].rows();
    int Nb_theta = st_E_surf->theta.size();

    output->E_r = ArrayXXc<Real>::Zero(Nb_lambda, Nb_theta);
    output->E_t = ArrayXXc<Real>::Zero(Nb_lambda, Nb_theta);
    output->E_f = ArrayXXc<Real>::Zero(Nb_lambda, Nb_theta);

    for (int m = -N_max; m <= N_max; m++) {
      complex<Real> exp_phase = exp(mp_im_unit<Real>()*static_cast<Real>(m)*phi0);
      output->E_r += st_E_surf->E_rm[m + N_max] * exp_phase;
      output->E_t += st_E_surf->E_tm[m + N_max] * exp_phase;
      output->E_f += st_E_surf->E_fm[m + N_max] * exp_phase;
    }
    return output;
  }

  /*
  Calculates the Riccati-Bessel function chi_n(x) = x*y_n(x) for many n and x
  Inputs:
    n - [N X 1] array of values n for chi_n(x)
    x - [X x 1] array of values x for chi_n(x)
  Output:
    Returns chi(x) [X x N]
  Dependencies:
    mp_pi, ArrBesselY
  */
  template <class Real>
  ArrayXXc<Real> vshRBchi(ArrayXr<Real> n, const ArrayXr<Real>& x) {
    ArrayXXc<Real> chi_x(x.size(), n.size());
    n += 0.5;
    for (int i = 0; i < x.size(); i++) {
      ArrayXr<Real> yx = ArrBesselY(n, x(i));
      if ((yx.isInf()).any())
        cerr << "Warning: Bessel (y) calculation went beyond precision in vshRBchi()" << endl;
      chi_x.row(i) = sqrt(static_cast<complex<Real>>(x(i)*mp_pi<Real>()/2))*yx;
    }
    return chi_x;
  }

  /*
  Calculates the Riccati-Bessel function psi_n(x) = x*j_n(x) for many n and x
  Inputs:
    n - [N X 1] array of values n for psi_n(x)
    x - [X x 1] array of values x for psi_n(x)
  Output:
    Returns psi(x) [X x N]
  Dependnencies:
    mp_pi, ArrBesselJ
  */
  template <class Real>
  ArrayXXc<Real> vshRBpsi(ArrayXr<Real> n, const ArrayXr<Real>& x) {
    ArrayXXc<Real> psi_x(x.size(), n.size());
    n += 0.5;
    for (int i = 0; i < x.size(); i++) {
      ArrayXr<Real> jx = ArrBesselJ(n, x(i));
      if ((jx == 0.0).any())
        cerr << "Warning: Bessel (j) calculation went beyond precision in vshRBpsi()" << endl;
      psi_x.row(i) = sqrt(static_cast<complex<Real>>(x(i)*mp_pi<Real>()/2))*jx;
    }
    return psi_x;
  }
}

#endif
