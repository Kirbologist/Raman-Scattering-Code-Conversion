/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Based on the SMARTIES MATLAB package by W.R.C. Somerville, B. Auguié, E.C. Le Ru
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.


This code contains all 'sph' SMARTIES functions that are used in Raman scattering calculations,
i.e. T-matrix related functions specific to particles with mirror-reflection symmetry.
*/

#ifndef SPH_HPP
#define SPH_HPP

#include "core.hpp"
#include "smarties_aux.hpp"
#include "vsh.hpp"
#include "misc.hpp"
#include <stdexcept>

using namespace Eigen;
using namespace std;

namespace Smarties {

  template <class Real>
  struct stFpRow {
    ArrayXXr<Real> S;
    ArrayXXr<Real> loss_prec_S;
  };

  template <class Real>
  struct stFpovx {
    Tensor3c<Real> Fpovx; // The matrix F^+/x
    ArrayXXc<Real> rb_chi; // The Riccati-Bessel function chi_n(x)
    ArrayXXc<Real> rb_psi; // The Riccati-Bessel function psi_k(s*x)
  };

  /* Struct containing modified Beseel function products for spheroids */
  template <class Real>
  struct stBessel {
    Tensor3c<Real> xi_psi; // the contributing part of xi_n psi_k
    Tensor3c<Real> psi_psi; // the product psi_n*psi_k
    ArrayXXc<Real> chi_n; // chi_n(x)
    ArrayXXc<Real> psi_n; // psi_n(x)
    ArrayXXc<Real> psi_k; // chi_k(s*x)
  };

  template <class Real>
  struct stBesselPrimes {
    Tensor3c<Real> xi_psi; // the Bessel products xi_n*psi_k
    Tensor3c<Real> xi_prime_psi; // xi'_n*psi_k (for n + k odd or n = k)
    Tensor3c<Real> xi_psi_prime; // xi_n*psi'_k (for n + k odd or n = k)
    Tensor3c<Real> xi_prime_psi_prime; // xi'_n*psi'_k (for n + k even)
    Tensor3c<Real> xi_psi_over_sxx; // xi_n*psi_k/(s*x^2) (for n + k even)
    // xi'_n*psi'_k + n*(n+1)*xi_n*psi_k/(s*x^2) (for n + k even)
    Tensor3c<Real> xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx;
    // xi'_n*psi'_k + k*(k+1)*xi_n*psi_k/(s*x^2) (for n + k even)
    Tensor3c<Real> xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx;
    ArrayXXc<Real> for_diag_Lt1; // for L^tilde1 in diagonal
    ArrayXXc<Real> for_diag_Lt2; // for L^tilde2 in diagonal
    ArrayXXc<Real> for_diag_Lt3; // for L^tilde3 in diagonal
  };

  template <class Real>
  struct stBesselProducts {
    unique_ptr<stBesselPrimes<Real>> st_xi_psi_all; // all Bessel products that start with xi*psi (for Q-matrix)
    unique_ptr<stBesselPrimes<Real>> st_psi_psi_all; // same but with all the xi replaced by psi (for P-matrix)
  };

  /*
  Struct describes the blocks of the eo or oe matrix of a P, Q, T or R matrix for a spheroidal scatterer.
  The first number of a block is its row index in the matrix,
  and the second number is its column index in the matrix.
  */
  template <class Real>
  struct st4M {
    ArrayXXc<Real> M11; // The top-left block of the matrix
    ArrayXXc<Real> M12; // The top-right block of the matrix
    ArrayXXc<Real> M21; // The bottom-left block of the matrix
    ArrayXXc<Real> M22; // The bottom-right block of the matrix
    int m; // The index m of term in the multipole expansion of the P/Q/T/R matrix that the eo or oe matrix
           // corresponds to. m denotes the projected angular momentum of a vector spherical wave function.
    ArrayXi ind1; // ind1 and ind2 indicate which row nad column indices of the P/Q/T/R matrix
    ArrayXi ind2; // are contained in each block.
  };

  /*
  Struct describes the eo and oe matrices that make up a pair of matrices (either P, Q or T, R)
  of a spheroidal scatterer for a single value of m.
  */
  template <class Real>
  struct stMat {
    std::array<st4M<Real>, 4> st_4M_list; // The four eo/oe matrices that make up the matrix pair
    vector<string> mat_list; // Array specifying what matrices are defined in the struct

    stMat() {
      for (int i = 0; i < 4; i++)
        this->st_4M_list[i] = st4M<Real>();
    }

    stMat(const st4M<Real>& base) {
      for (int i = 0; i < 4; i++)
        this->st_4M_list[i] = base.st_4M_list[i];
      this->mat_list = base.mat_list;
    }
  };

  /* Specialisation of stMat for describing a P, Q matrix pair */
  template <class Real>
  struct stPQ : stMat<Real> {
    using stMat<Real>::stMat;
    inline st4M<Real>& st_4M_P_eo() { return this->st_4M_list[0]; }
    inline st4M<Real>& st_4M_P_oe() { return this->st_4M_list[1]; }
    inline st4M<Real>& st_4M_Q_eo() { return this->st_4M_list[2]; }
    inline st4M<Real>& st_4M_Q_oe() { return this->st_4M_list[3]; }
  };

  /*
  Calculates the matrix u, for one value of n, for the series implementation of the Bessel function product.
  This is valid for all k up to k = n.

  The matrix u is defined as [eq. B.10 of JQSRT 2013] u_{rb} = 2^b (d/dX)^b [X^(n - 1/2)*(1 - 1/X)^r] | X = 1
  which leads to the relations [eq. B.16 of JQSRT 2013]
    u_{0,0} = 0
    u_{r,0} = 1 (r > 0)
    u_{r,b+1} = (n - 1/2 - 2r)u_{rb} - (n - 1/2 - r)u_{r + 1,b}
  Inputs:
    n - the value of n required
  Output:
    u - [(B + 2) x (B + 1)] where B = floor(n/2), the matrix u as defined above
  */
  template <class Real>
  ArrayXXr<Real> sphGetUforFp(int n) {
    int b_max = n/2;
    ArrayXXr<Real> u(b_max + 2, b_max + 1); // [r x b], both from 0 to b_max
    u.setZero();
    u(0, 0) = 1;
    for (int b = 0; b < b_max; b++) {
      u(0, b + 1) = (2*n - 1)*u(0, b) - (2*n - 1)*u(1, b);
      for (int r = 1; r <= b + 1; r++) {
        u(r, b + 1) = (2*n - 1 - 4*r)*u(r, b) - (2*n - 1 - 2*r)*u(r + 1, b) + 2*r*u(r - 1, b);
      }
    }
    return u;
  }

  /*
  Calculate the product F^+_{nk} = P^+(x*chi_n(x)_psi_k(s*x)) for one n and all k with cancellations,
  i.e. all 0 <= k <= n - 4 and n + k even, for one value of s (so one wavelength), and possibly many x's.

  This uses a seris implementation to calculate the product, and should treat well both
  large and small (near 1) values of s, and large and small values of x.
  The result of this may be used along with a recursion scheme in order to calculate
  all the required terms for the EBCM integrals for spheroids.
  See appndix B of [JQSRT 2013] for further details
  Inputs:
    n - the value of n to use; typically we will do this for the last row, i.e. n = N + 1,
        where N is the highest value of n required in calculations.
    s - the relative refractive index of the particle
    x - [X x 1] the value of x to calculate for
  Output:
    Returns a unique pointer to a struct containing the following values:
    S - [(K = n - 3) x X] the values that are returned. This is for k = 0 to k = n - 4,
        and entries for odd n + k are zero.
    loss_prec_S - [K = n - 3 x X] potential loss of precision estimated as
                  abs(Max term in series)/S. Problem if loss_prec_S is large, i.e. 10^4 means 4-digit loss.
  Dependencies:
    sphGetUforFp
  */
  template <class Real>
  unique_ptr<stFpRow<Real>> sphGetFpRow(int n, Real s, const ArrayXr<Real>& x) {
    int num_x = x.size();
    RowArrayXr<Real> x_squared = x.pow(2).transpose(); // [1 x X] get x^2 as a row
    ArrayXXr<Real> u = sphGetUforFp<Real>(n); // [(B + 2) x (B + 1)] where B = floor(n/2)
    Real alpha_bar_k = 1; // \bar{alpha}_k [eq. B.23]
    int q_min;
    int q_int;

    ArrayXXr<Real> S = ArrayXXr<Real>::Zero(n - 3, num_x); // [K x X] this is the partial sum in the series
    // in eq. B.20 needed for F^+, we will add terms to it sequentially
    ArrayXXr<Real> max_term_S = ArrayXXr<Real>::Zero(n - 3, num_x); // [K x X] keeps track of maximum term
    // in the series for potential loss of precision
    ArrayXXr<Real> loss_prec_S = ArrayXXr<Real>::Zero(n - 3, num_x); // [K x X] keeps track of loss of precision

    // Computing the beta_{jq} [eq. B.14]
    // Note, the original size of beta in the MATLAB package is a bit too small for even numbered n.
    ArrayXXr<Real> beta = ArrayXXr<Real>::Zero(n + 1 - (n % 2), n - (n % 2)); // [j x q] 0<=j<=q, 1<=q<=q_int-1
    beta(0, 0) = 1; // j = 0, q = 2
    beta(1, 0) = beta(0, 0) * (s - 1) * (s + 1); // j = 1, q = 1
    // (s - 1) * (s + 1) is used instead of s^2 - 1 to minimise rounding error problems
    beta(0, 1) = 1; // j = 0, q = 2
    beta(1, 1) = beta(0, 1) * (s - 1) * (s + 1) * 2; // j = 1, q = 2
    beta(2, 1) = beta(1, 1) * (s - 1) * (s + 1) / 2; // j = 2, q = 2

    int num_to_test = 3; // this is the number of consecutive tests required to consider
                         // if convergence has been reached (a single term goes to zero by chance)

    for (int k = n - 4; k >= 0; k -= 2) {
      q_min = (n - k)/2 - 1; // expression is mathematically always an integer
      q_int = n - k; // "intermediate" value of q, where we switch which method applies
      alpha_bar_k /= -(n - k - 2); // [eq. B.24]
      // fill the last two beta columns here, as this is the first time that we need them
      beta(0, q_int - 2) = 1; // j = 0, q = q_int - 1
      for (int j = 1; j < q_int; j++)
        beta(j, q_int - 2) = beta(j - 1, q_int - 2)*(s - 1)*(s + 1)*(q_int - j)/j;
      beta(0, q_int - 1) = 1; // j = 0, q = q_int
      for (int j = 1; j <= q_int; j++)
        beta(j, q_int - 1) = beta(j - 1, q_int - 1)*(s - 1)*(s + 1)*(q_int - j + 1)/j;

      // more initialisations
      RowArrayXr<Real> alpha_q(num_x);
      alpha_q.fill(alpha_bar_k); // [1 x X] alpha_{qnk}(x) [eq. B.21] for q = q_min
      RowArrayXb x_not_converged = RowArrayXb::Ones(num_x);
      RowArrayXr<Real> current_term = RowArrayXr<Real>::Zero(num_x);
      RowArrayXi counter = RowArrayXi::Zero(num_x);
      RowArrayXb test = RowArrayXb::Zero(num_x);

      // Method 1 for q_min <= q <= n - k - 1
      for (int q = q_min; q < q_int; q++) {
        int b = n - k - q - 1;
        Real gamma_qnk = 0; // computes gamma_{qnk} from the sum in eq. B.13
        for (int j = max(0, 2 * q + k - n + 1); j <= q; j++)
          gamma_qnk += beta(j, q - 1) * u(q - j, b);

        // For the finite series of method 1, we only test convergence for |s| >= 2.
        // If s is close to 1, the sereis may appear to converge but the terms increase again,
        // so to be safe, we take all the terms
        if (abs(s) < 2) { // |s| < 2, no convergence test
          current_term = alpha_q * gamma_qnk;
          S.row(k) += current_term;
          max_term_S.row(k) = max_term_S.row(k).max(abs(current_term));
        } else { // test if convergence of the series S [eq. B.20] has been obtained
          test.fill(false); // reset the test flag to converged
          current_term = x_not_converged.select(alpha_q * gamma_qnk, current_term);
          if (x_not_converged.any())
            // converged if added term is less than epsilon times current value of the series.
            // test = false if this is the case, test = true for non-converged x-values
            test = x_not_converged.select(abs(current_term) > mp_eps<Real>() * abs(S.row(k)), test);
          if (test.any()) { // for the non-converged x values
            if (LogicalSlice(current_term, test).isInf().all()) // if all non-converged x values are infinite
              cout << "Problem (1) in sphGetFpovx..." << endl;
            // add the current term to the series
            S.row(k) += test.select(current_term, 0);
            max_term_S.row(k) = test.select(abs(current_term).max(max_term_S.row(k)), max_term_S.row(k));
          }

          // Full convergence is assumed if the convergence test is fulfilled for num_to_test consecutive terms
          counter += 1 - test.template cast<int>();
          counter *= 1 - test.template cast<int>();
          x_not_converged = (counter < num_to_test);
          if (!(x_not_converged.any()))
            break;
        }

        // Update alpha_q for next term (using eq. B.22)
        alpha_q *= -x_squared / (2*(q + 1));
      }

      // We now move on to method 2 for q >= q_int [sec. B.2]
      // Reinitialising coutners - convergence is now tested in all cases (otherwise the series would be infinite)
      x_not_converged.fill(true);
      counter.setZero();

      // Initialising c_{qqnk} [eq. B.19]
      Real c_qqnk = static_cast<Real>(1.0)/(2*k + 1); // this is for q = q_int = n - k
      int q = q_int;
      while (true) { // upper limit of series will be decided by convergence test
        Real c_iqnk = c_qqnk; // i = q term
        // initialise gamma_qnk = c_qqnk (this is i = q) and add the other term according to eq. B.3
        Real gamma_qnk = c_qqnk;
        for (int i = q - 1; i >= 0; i--) { // sum from i = 0 to i = q - 1
          // calculates c_iqnk by recurrence using eq. B.17
          c_iqnk *= pow(s, 2)*(i + 1)*(2*i + 1 - 2*n)/(q - i)/(2*k + 2*q - 2*i + 1);
          gamma_qnk += c_iqnk;
        }

        // convergence testing
        test.fill(false);
        for (int i = 0; i < num_x; i++) {
          if (x_not_converged(i)) {
            current_term(i) = alpha_q(i) * gamma_qnk;
            test(i) = abs(current_term(i)) > mp_eps<Real>() * abs(S(k, i));
          }
        }
        if (test.any()) {
          if ((1 - test + current_term.isInf()).all()) // if all non-converged x values are infinite
            cout << "Problem (2) in sphGetFpovx..." << endl;
          // Add the current term to the series
          for (int i = 0; i < num_x; i++) {
            if (test(i)) {
              S(k, i) += current_term(i);
              max_term_S(k, i) = max(abs(current_term(i)), max_term_S(k, i));
            }
          }
        }

        counter += 1 - test.template cast<int>();
        counter *= 1 - test.template cast<int>();
        x_not_converged = (counter < num_to_test);

        if (!(x_not_converged.any()))
          break;

        // update alpha_q for next term (using eq. B.22)
        for (int i = 0; i < num_x; i++) {
          if (x_not_converged(i))
            alpha_q(i) *= -x_squared(i) / (2*(q + 1));
        }
        // update startomg [point for c_qqnk using eq. B.18]
        c_qqnk /= (2*q + 1 - 2*n);
        q++;
      }

      // estimated loss of precision
      loss_prec_S.row(k) = abs(max_term_S.row(k) / S.row(k));
      // Finally, we scale everything with the missing factor (see eq. B.20)
      S.row(k) *= -pow(s, k + 1);
    }

    auto output = make_unique<stFpRow<Real>>();
    output->S = S;
    output->loss_prec_S = loss_prec_S;
    return output;
  }

  /*
  Calculate the matrix F^+/x (see eq. 46 of JQSRT 2013) where Fpovx = p^+(x*chi_n(x)*psi_k(s*x))/x
  for 0 <= (n,k) <= N_max, for one s (wavelength), and x can be a vector of length [T x 1]

  This function calculates F^+/x = F/x in regions where there are no cancellations, and calls sphGetFpRow
  to calculate the last row of F^+ for the region where there are cancellations. It then fills up the matrix
  using the "westward" recursion scheme (solving for n, k - 1, see fig. 3c). See sec. 4.2 for further details.
  Inputs:
    N_max - the maximum value of n that is desired
    s - the relative refractive index of the particle
    x - [T x 1] the values of x at which the function is evaluated
  Output:
    Returns a unique pointer to a stFpovx struct, where its Fpovx has size [(N + 1) x (N + 1) x T]
    and the Riccati-Bessel functions rb_chi and rb_psi have size [T x (N + 1)].
  Dependencies:
    TensorCast, sphGetFpRow, vshRBchi, vshRBpsi
  */
  template <class Real>
  unique_ptr<stFpovx<Real>> sphGetFpovx(int N_max, Real s, const ArrayXr<Real>& x) {
    int num_x = x.size();

    auto output = make_unique<stFpovx<Real>>();
    output->rb_chi = vshRBchi<Real>(ArrayXr<Real>::LinSpaced(N_max + 1, 0, N_max), x); // [T x (K + 1)]
    output->rb_psi = vshRBpsi<Real>(ArrayXr<Real>::LinSpaced(N_max + 1, 0, N_max), s*x); // [T x (N + 1)]
    output->Fpovx = Tensor3c<Real>(N_max + 1, N_max + 1, num_x); // [(N + 1) x (K + 1) x T]
    output->Fpovx.setZero();

    // To initisalise the westward recurrence, we need the last row (n = N), and the subdiagonal terms.
    // This first calculates the last row (n = N, all k) for Fpovx
    unique_ptr<stFpRow<Real>> FpRow = sphGetFpRow<Real>(N_max, s, x);
    ArrayXXc<Real> tmp;
    for (int i = 0; i < N_max - 3; i++) {
      tmp = (FpRow->S.row(i) / x.transpose()).template cast<complex<Real>>();
      output->Fpovx.chip(N_max, 0).chip(i, 0) = TensorCast(tmp, num_x);
    }

    // Do n, k - 1 (westward) recursion [solving for F^+_{n,k-1} in eq. 51 (dividing all terms by x)]
    // We only do n + k even
    Tensor1c<Real> col_shape(num_x);
    for (int k = N_max; k >= 0; k--) {
      // First fill in the matrix where there are no cancellations (n = 0..k+2)
      // N_min  = 0 (if k even) or 1 (if k odd), N_max = min(k + 2, N).
      // Here, F^+/x = F/x = chi_n(x) * psi_k(s*x). x-values are placed on the third dimension.
      for (int i = k % 2; i <= min(k + 2, N_max); i += 2)
        output->Fpovx.chip(i, 0).chip(k, 0) = TensorCast(output->rb_chi.col(i) * output->rb_psi.col(k), num_x);
      // Now calculate all terms in column k-1 for all relevant problematic n
      if (k) { // for all non-zero k
        for (int n = k + 3; n < N_max; n += 2) // n = k + 3 to N-1 with only n + k odd
          output->Fpovx.chip(n, 0).chip(k - 1, 0) = (output->Fpovx.chip(n + 1, 0).chip(k, 0) +
              output->Fpovx.chip(n - 1, 0).chip(k, 0)) *
              col_shape.constant(static_cast<Real>(2*k + 1)/(2*n + 1)/s) - output->Fpovx.chip(n, 0).chip(k + 1, 0);
      }
    }
    return output;
  }

  /*
  Calculates the products xi(n, x)*psi(k, s*x) and psi(n, x)*psi(k, s*x) for n, k from 0 to N_max + 1.
  For xi_psi, this calculates only the part that does not exhibit cancellations,
  i.e. F^+_{nk}/x (see eq. 46 of JQSRT 2013). NB >= N multipoles are used for the calculation
  (for improved precision).

  These are calculated for one s(wavelength), and multiple x, for (n, k) = 0..N_max + 1.
  The two extreme (0 and N_max + 1) are needed to compute the products involving derivatives.
  Only n + k even are calculated and returned.
  The single regular Bessel (not products) are also not returned.
  Inputs:
    N_max - the maximum value N required for the integrals
    s - the relative refractive index of the particle
    x - [T x 1] the values of x to calculate this at.
    NB - the number of n that are used to calculate the Bessel products (NB is required to be >= N)
  Output:
    Returns a unique pointer to a struct containing various Bessel functions and their products.
  Dependencies:
    sphGetFpovx, vshRBpsi
  */
  template <class Real>
  unique_ptr<stBessel<Real>> sphGetXiPsi(int N_max, Real s, const ArrayXr<Real>& x, int NB) {
    assert(NB >= N_max);
    int num_x = x.size();
    auto output = make_unique<stBessel<Real>>();
    // This call takes care of the chi*psi products
    unique_ptr<stFpovx<Real>> chi_psi_prod = sphGetFpovx<Real>(NB + 1, s, x);
    // chi_psi is [(NB + 2) x (NB + 2) x X]
    // chi_n , psi_k are [X x NB + 2]

    // The rest calculates the normal psi*psi products using standard C++/Boost library functions
    // to compute the Bessel functions (only up to N + 2)
    output->psi_psi = Tensor3c<Real>(N_max + 2, N_max + 2, num_x);
    output->psi_psi.setZero();
    output->psi_n = vshRBpsi<Real>(ArrayXr<Real>::LinSpaced(N_max + 2, 0, N_max + 1), x); // [X x (N + 2)]

    // To reduce computations, only the products we need are calculated
    // i.e. only n + k even (others are left as zero)
    for (int n = 0; n < N_max + 2; n++) {
      for (int k = n % 2; k < N_max + 2; k += 2)
        output->psi_psi.chip(n, 0).chip(k, 0) =
            TensorCast(chi_psi_prod->rb_psi.col(k) * output->psi_n.col(n), num_x); // [1 x K/2 x X]
    }

    std::array<int, 3> offsets = {0, 0, 0};
    std::array<int, 3> extents = {N_max + 2, N_max + 2, num_x};
    output->xi_psi = output->psi_psi + output->psi_psi.constant(mp_im_unit<Real>()) *
        chi_psi_prod->Fpovx.slice(offsets, extents);
    output->chi_n = chi_psi_prod->rb_chi.leftCols(N_max + 2);
    output->psi_k = chi_psi_prod->rb_psi.leftCols(N_max + 2);
    return output;
  }

  /*
  Calculates the derivatives of the bessel function products needed in the integrals,
  as well as other terms calculated in a similar manner.
  Inputs:
    prods - [(N + 2) x (N + 2) x X] Bessel function products (either xi_n(x)*psi_k(s*x) or psi_n(x)*psi_k(s*x)),
            for n = 0 to n = N + 1, for n + k even. These are obtained from sphGetXiPsi.
  Output:
    Returns unique pointer to struct containing tensors for the derivatives product.
    Each tensor is [(N + 1) x (N + 1) x X] (for n,k = 0..N). The cases where n or k = 0 is included just
    for padding, and it makes the indexing much simpler. Member variable names use "xi_psi" for convenience,
    but may also correspond to psipsi if prods = psi_psi.
  */
  template <class Real>
  unique_ptr<stBesselPrimes<Real>> sphGetBesselProductsPrimes(const Tensor3c<Real>& prods) {
    assert(prods.dimension(0) == prods.dimension(1));
    int N = prods.dimension(0) - 2; // [N x N]
    int X = prods.dimension(2);

    auto output = make_unique<stBesselPrimes<Real>>();
    std::array<long int, 3> offsets = {0, 0, 0};
    std::array<long int, 3> extents = {N + 1, N + 1, X};
    output->xi_psi = prods.slice(offsets, extents); // result already exists for xi_psi
    output->xi_prime_psi = Tensor3c<Real>(N + 1, N + 1, X);
    output->xi_psi_prime = Tensor3c<Real>(N + 1, N + 1, X);
    output->xi_prime_psi_prime = Tensor3c<Real>(N + 1, N + 1, X);
    output->xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx = Tensor3c<Real>(N + 1, N + 1, X);
    output->xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx = Tensor3c<Real>(N + 1, N + 1, X);
    output->xi_psi_over_sxx = Tensor3c<Real>(N + 1, N + 1, X);

    output->xi_prime_psi.setZero();
    output->xi_psi_prime.setZero();
    output->xi_prime_psi_prime.setZero();
    output->xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx.setZero();
    output->xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx.setZero();
    output->xi_psi_over_sxx.setZero();

    // Now calculate the products involving derivatives
    Tensor1c<Real> col_shape(X);
    for (int n = 1; n <= N; n++) {
      for (int k = 2 - n % 2; k <= N; k += 2) { // n + k even only
        int kkp1 = k*(k + 1);
        int nnp1 = n*(n + 1);

        // eq. 61
        output->xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx.chip(n, 0).chip(k, 0) = (
            col_shape.constant(static_cast<complex<Real>>((k + n + 1)*(k + 1))) * prods.chip(n - 1, 0).chip(k - 1, 0) +
            col_shape.constant(static_cast<complex<Real>>(kkp1 - k*(n + 1))) * prods.chip(n - 1, 0).chip(k + 1, 0) +
            col_shape.constant(static_cast<complex<Real>>(kkp1 - (k + 1)*n)) * prods.chip(n + 1, 0).chip(k - 1, 0) +
            col_shape.constant(static_cast<complex<Real>>(kkp1 + k*n)) * prods.chip(n + 1, 0).chip(k + 1, 0)) /
            col_shape.constant(static_cast<complex<Real>>((2*n + 1)*(2*k + 1)));

        // eq. 62
        output->xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx.chip(n, 0).chip(k, 0) = (
            col_shape.constant(static_cast<complex<Real>>(nnp1 + (k + 1)*(n + 1)))*prods.chip(n - 1, 0).chip(k - 1, 0) +
            col_shape.constant(static_cast<complex<Real>>((n - k)*(n + 1))) * prods.chip(n - 1, 0).chip(k + 1, 0) +
            col_shape.constant(static_cast<complex<Real>>((n - k)*n)) * prods.chip(n + 1, 0).chip(k - 1, 0) +
            col_shape.constant(static_cast<complex<Real>>(nnp1 + k*n)) * prods.chip(n + 1, 0).chip(k + 1, 0)) /
            col_shape.constant(static_cast<complex<Real>>((2*n + 1)*(2*k + 1)));
      }

      for (int k = 1 + n % 2; k <= N; k += 2) { // n + k odd only
        // eqs. 59-60
        output->xi_prime_psi.chip(n, 0).chip(k, 0) =
            (prods.chip(n - 1, 0).chip(k, 0) * col_shape.constant(static_cast<complex<Real>>(n + 1)) -
            col_shape.constant(static_cast<complex<Real>>(n))*prods.chip(n + 1, 0).chip(k, 0)) /
            col_shape.constant(static_cast<complex<Real>>(2*n + 1));
        output->xi_psi_prime.chip(n, 0).chip(k, 0) =
            (prods.chip(n, 0).chip(k - 1, 0) * col_shape.constant(static_cast<complex<Real>>(k + 1)) -
            col_shape.constant(static_cast<complex<Real>>(k))*prods.chip(n, 0).chip(k + 1, 0)) /
            col_shape.constant(static_cast<complex<Real>>(2*k + 1));
      }
    }
    return output;
  }

  /*
  Builds structs containing matrices of modified Bessel function products,
  suitable for calculating T-matrix integrals for spheroids. See JQSRT 123, 153 (2013)
  for further details. The modified Bessel products correspond to F^+_{nk}(s, x) / x (eq. 46)
  and equivalents as defined in sec. 4.3.

  The struct for psipsi also contains functions needed for the diagonals of P and Q matrices
  Inputs:
    N_max - maximum value of N required in the integrals
    s - relative refractive index of the particle
    x - [T X 1] the values of x at which the integrals are evaluated
    NB - the number of N that should be used when calculating the Bessel function products.
         Note it's required that NB >= N for proper functionality
  Output:
    Returns unique pointer to struct containing the matrices for the products xi*psi and psi*psi,
    as used in Q and P matrices. The member tensors all have size [(N + 1) X (N + 1) X T] or [(N + 1) X T].
    The cases where n or k = 0 are included as padding, as it makes the indexing simpler.
  Dependencies:
    sphGetBesselProductsPrimes, sphGetXiPsi
  */
  template <class Real>
  unique_ptr<stBesselProducts<Real>> sphGetModifiedBesselProducts(int N_max, Real s, const ArrayXr<Real>& x, int NB) {
    auto output = make_unique<stBesselProducts<Real>>();

    // start with F^+/x i.e. xi*psi and psi*psi
    unique_ptr<stBessel<Real>> prods = sphGetXiPsi<Real>(N_max, s, x, NB);
    // Then deduce from them the equivalent with derivatives
    output->st_xi_psi_all = sphGetBesselProductsPrimes<Real>(prods->xi_psi);
    output->st_psi_psi_all = sphGetBesselProductsPrimes<Real>(prods->psi_psi);

    // For diagonals of Q11 and P11 (for L^tilde1, eq. 25)
    // s*(xi'_n*psi_n)(s*x) - (xi_n*psi'_n)(s*x) = s*(xi_n*psi_{n+1})(s*x) - (xi_{n+1}*psi_n)(s*x) [eq. 65]
    // Same for P
    ArrayXXc<Real> psi_np1_psi_n =
        (prods->psi_n(all, seq(1, last)) * prods->psi_k(all, seq(0, last - 1))).transpose(); // [(N + 1) x X]
    ArrayXXc<Real> psi_n_psi_np1 =
        (prods->psi_n(all, seq(0, last - 1)) * prods->psi_k(all, seq(1, last))).transpose(); // [(N + 1) x X]
    ArrayXXc<Real> xi_np1_psi_n = psi_np1_psi_n + mp_im_unit<Real>() *
        (prods->chi_n(all, seq(1, last)) * prods->psi_k(all, seq(0, last - 1))).transpose(); // [(N + 1) x X]
    ArrayXXc<Real> xi_n_psi_np1 = psi_n_psi_np1 + mp_im_unit<Real>() *
        (prods->chi_n(all, seq(0, last - 1)) * prods->psi_k(all, seq(1, last))).transpose(); // [(N + 1) x X]
    output->st_psi_psi_all->for_diag_Lt1 = s*psi_n_psi_np1 - psi_np1_psi_n;
    output->st_xi_psi_all->for_diag_Lt1 = s*xi_n_psi_np1 - xi_np1_psi_n;
    // For diagonals of Q22 and P22 (for L^tilde2 and L^tilde3, eq. 26)
    ArrayXXc<Real> psi_n_psi_n =
        (prods->psi_n(all, seq(0, last - 1)) * prods->psi_k(all, seq(0, last - 1))).transpose(); // [(N + 1) x X]
    ArrayXXc<Real> xi_n_psi_n = psi_n_psi_n + mp_im_unit<Real>() *
        (prods->chi_n(all, seq(0, last - 1)) * prods->psi_k(all, seq(0, last - 1))).transpose(); // [(N + 1) x X]

    ArrayXc<Real> n_vec = ArrayXr<Real>::LinSpaced(N_max + 1, 1, N_max + 1);
    output->st_psi_psi_all->for_diag_Lt2 = psi_n_psi_np1 - s*psi_np1_psi_n + (s - 1)*(s + 1)/s *
        (n_vec.matrix() * (1/x.transpose()).matrix()).array() * psi_n_psi_n;
    output->st_xi_psi_all->for_diag_Lt2 = xi_n_psi_np1 - s*xi_np1_psi_n + (s - 1)*(s + 1)/s *
        (n_vec.matrix() * (1/x.transpose()).matrix()).array() *xi_n_psi_n;
    RowArrayXc<Real> tmp = s*x.pow(2).transpose();
    output->st_psi_psi_all->for_diag_Lt3 = psi_n_psi_n.rowwise() / tmp;
    output->st_xi_psi_all->for_diag_Lt3 = xi_n_psi_n.rowwise() / tmp;

    return output;
  }

  /*
  Calculates the geometry for spheroids and the node and weights of the Gaussian quadrature.
  By reflection symmetry, only theta from 0 to pi/2 are used; quadrature weights are
  doubled to compensate. If theta is specified, thet quadrature properties are not computed.
  Inputs:
    Nb_theta - the number of angles to use (in half of the shape)
    a - the semi-axis length for axes along x, y
    c - the semi-aixs length along z (axis of rotation)
    theta - theta values to calculate the geometry for. Else uses auxPrepareIntegrals to genereate quadrature
  Output:
    Returns pointer to struct containing geometry information. output->type = PTS if theta points were specified,
    otherwise output->type = GAUSS. Array members of the struct have size T = Nb_theta.
  Dependencies:
    auxPrepareIntegrals
  */
  template <class Real>
  unique_ptr<stRtfunc<Real>> sphMakeGeometry(int Nb_theta, Real a, Real c,
      const unique_ptr<ArrayXr<Real>> theta = unique_ptr<ArrayXr<Real>>()) {
    sInt type;
    auto output = make_unique<stRtfunc<Real>>();
    if (theta) { // points only
      output->w_theta = ArrayXr<Real>::Zero(theta->size());
      output->theta = *theta;
      output->Nb_theta = theta->size();
      type = PTS;
    } else {
      type = GAUSS;
      unique_ptr<stRtfunc<Real>> tmp = auxPrepareIntegrals<Real>(2*Nb_theta, type);
      output->theta = tmp->theta(seq(0, Nb_theta - 1));
      output->w_theta = tmp->w_theta(seq(0, Nb_theta - 1))*2;
      output->Nb_theta = Nb_theta;
    }

    output->a = a;
    output->c = c;
    output->type = type;

    // Defines geometry
    ArrayXr<Real> sin_t = sin(output->theta);
    ArrayXr<Real> cos_t = cos(output->theta);

    // r for spheroid with different semi-major and semi-minor axes at the origin
    output->r = a*c/sqrt(pow(c*sin_t, 2) + pow(a*cos_t, 2));
    output->dr_dt = (pow(a, 2) - pow(c, 2))/pow(a*c, 2)*sin_t * cos_t * output->r.pow(3);

    output->h = max(a, c)/min(a, c); // aspect ratio, >= 1
    output->r0 = pow(pow(a, 2)*c, static_cast<Real>(1.0)/3); // the radius of volume-equivalent sphere

    return output;
  }

  /*
  Determines required NB values for calculations of the Bessel function product for arguments `x` to be valid
  up to a given accuracy `acc`. This is asserted by checking the relative accuracy of
  the last row of the F^+/x matrix.
  Inputs:
    N_req - the minimum value of N that is required to be sufficiently accurate
    s - relative refractive index of the particle
    x - array containing the sizes to check, which can be representative of the whole particle
    acc - Required relative accuracy that the results are meant to converge within
    N_min - minimum value for NB (where the search starts)
  Output:
    Returns the minimum multipole order `NB` needed for convergence
  Dependencies:
    sphGetFpovx, TensorSlice
  */
  template <class Real>
  int sphCheckBesselConvergence(int N_req, Real s, const ArrayXr<Real>& x, Real acc, int N_min) {
    // Initial guess of N is the smallest possible size
    int NB_start = max(N_req, N_min);
    int NB = NB_start;
    unique_ptr<stFpovx<Real>> prod = sphGetFpovx<Real>(NB, s, x);
    bool to_continue = true;
    int max_N = 500; // The maximum value of N we are prepared to go to
    int NB_step = 16;

    int NB_next;
    unique_ptr<stFpovx<Real>> prod_new;
    Tensor<Real, 0> rel_acc_ee, rel_acc_oo;
    Real rel_acc;
    ArithmeticSequence<long int, long int, long int> seq1(0, N_req, 2);
    ArithmeticSequence<long int, long int, long int> seq2(1, N_req, 2);
    ArithmeticSequence<long int, long int, long int> seq3(0, x.size() - 1); // prod->Fpovx has no. of columns = x.size()
    long int dim1 = seq1.size();
    long int dim3 = seq3.size();
    Tensor3c<Real> tmp(dim1, dim1, dim3);
    while (to_continue && NB < max_N) {
      NB_next = NB + NB_step;
      prod_new = sphGetFpovx<Real>(NB_next, s, x);

      // Worst relative accuracy in all matrices up to n = N_req + 1 (+1 for derivatives to also be accurate)
      tmp = TensorSlice<Real>(prod->Fpovx, seq1, seq1, seq3) /
          TensorSlice<Real>(prod_new->Fpovx, seq1, seq1, seq3) - tmp.constant(static_cast<complex<Real>>(1));
      rel_acc_ee = tmp.abs().real().maximum();
      tmp = TensorSlice<Real>(prod->Fpovx, seq2, seq2, seq3) /
          TensorSlice<Real>(prod_new->Fpovx, seq2, seq2, seq3) - tmp.constant(static_cast<complex<Real>>(1));
      rel_acc_oo = tmp.abs().real().maximum();
      rel_acc = max(rel_acc_ee(0), rel_acc_oo(0));
      if (rel_acc < acc)
        to_continue = false; // we have found sufficiently high n, fine-tune NB by going back and changing step
      else {
        NB = NB_next;
        prod = move(prod_new);
      }
    }

    // Repeat with small steps to fine-tune NB
    if (NB > NB_start) {
      NB -= NB_step;
      NB_step = 1;
      to_continue = true;
      while (to_continue && NB < max_N) {
        NB += NB_step;
        prod = sphGetFpovx<Real>(NB, s, x);

        // We now use prod_new from the previous step as our reference
        // Worst relative accuracy in all matrices up to n = N_req + 1 (+1 for derivatives to also be accurate)
        tmp = TensorSlice<Real>(prod->Fpovx, seq1, seq1, seq3) /
            TensorSlice<Real>(prod_new->Fpovx, seq1, seq1, seq3) - tmp.constant(static_cast<complex<Real>>(1));
        rel_acc_ee = tmp.abs().maximum().real();
        tmp = TensorSlice<Real>(prod->Fpovx, seq2, seq2, seq3) /
            TensorSlice<Real>(prod_new->Fpovx, seq2, seq2, seq3) - tmp.constant(static_cast<complex<Real>>(1));
        rel_acc_oo = tmp.abs().maximum().real();
        rel_acc = max(rel_acc_ee(0), rel_acc_oo(0));
        if (rel_acc < acc)
          to_continue = false; // We have found sufficiently high n
      }
    }

    if (to_continue)
      throw(runtime_error("Problem in sphEstimateNB: convergence was not achieved"));
    return NB;
  }

  /*
  Finds the number of n (NB) required for accurate modified Bessel products
  Inputs:
    NQ - minimum number of multipoles required
    st_geometry - unique pointer to a struct containing geometric informaiton, as from sphMakeGeometry
    params - unique pointer to a struct containing simulation parameters k1 and s
    acc - Relative accuracy required, default is 1e-13
  Output:
    NB scalar with number of multipoles required
  Dependencies:
    sphCheckBesselConvergence, Subtensor2ArrMap
  */
  template <class Real>
  int sphEstimateNB(int NQ, const unique_ptr<stRtfunc<Real>>& st_geometry,
      const unique_ptr<stParams<Real>>& params, Real acc = 1e-13) {
    ArrayXr<Real> s = params->s;
    ArrayXr<Real> k1 = params->k1;

    // Find max size parameters
    ArrayXr<Real> x_max = {{k1.maxCoeff() * st_geometry->r.maxCoeff()}};
    // Find max and min relative refractive index
    typename ArrayXr<Real>::Index ind;
    abs(s).minCoeff(&ind);
    Real s_min = s(ind);
    abs(s).maxCoeff(&ind);
    Real s_max = s(ind);

    // Find required NB for each of the 2 extreme cases
    int N1 = sphCheckBesselConvergence(NQ, s_max, x_max, acc, NQ);
    int NB = sphCheckBesselConvergence(NQ, s_min, x_max, acc, N1);
    return NB;
  }

  /*
  Calculates P,Q matrices for a spheroid using the algorithm of [JQSRT 123, 153 (2013)] for all m in abs_m_vec.
  stRtfunc should contain the geometry, as obtained from sphMakeGeometry,
  and stParams should have fields k1 and s, both only containing 1 entry. The theta range for the geometry
  should either be [0, pi/2] or [0, pi].

  This makes use of modified Riccati-Bessel function products that are optimised for spheroids,
  i.e. they do not contain terms that integrate to zero, thus removing the loss of precision
  which affects traditional codes.
  Inputs:
    N_max - the maximum multipole order to return
    abs_m_vec - An array of size M containing values of m that the calculations should be carried out for.
                The order should be {0, 1, 2, ..., m_max}
    Rt_func - unique pointer to stRtfunc struct containing geometric information, as from sphMakeGeometry
    params - unique pointer to stParams struct containing simulation parameters
    NB - The number of multipoles that should be used to calculate the Bessel function products (NB >= N)
  Output:
    Returns a std::vector of size M containing unique pointers to stPQ structs, one for each m.
  Dependencies:
    sphGetModifiedBesselProducts, vshPinmTaunm
  */
  template <class Real>
  vector<unique_ptr<stPQ<Real>>> sphCalculatePQ(int N_max, const ArrayXi& abs_m_vec,
      const unique_ptr<stRtfunc<Real>>& Rt_func, const unique_ptr<stParams<Real>>& params, int NB = -1) {
    if (params->s.size() > 1 || params->k1.size() > 1)
      throw(runtime_error("params->s and params->k1 must be scalar"));
    if (NB < N_max)
      NB = N_max;
    int M = abs_m_vec.size();
    //if (params->output) (params shouldn't have an output member, perhaps use stOptions instead)
    cout << "sphCalculatePQ: Calculate P, Q for " << M << " m-values with N_Q = " <<
        N_max << ", NB = " << NB << ", N_Theta = " << Rt_func->Nb_theta << endl;

    vector<unique_ptr<stPQ<Real>>> output(M);
    for (int i = 0; i < M; i++)
      output[i] = make_unique<stPQ<Real>>();

    // The following are used for all m
    Real s = params->s(0);
    Real k1 = params->k1(0);
    int Nb_theta = Rt_func->Nb_theta;
    ArrayXr<Real> x = Rt_func->r * k1; // [T x 1] x(theta)
    ArrayXr<Real> x_theta = Rt_func->dr_dt * k1; // [T x 1] derivative of x(theta)

    unique_ptr<stPinmTaunm<Real>> stPT = vshPinmTaunm(N_max, Rt_func->theta); // Angular functions, pi and tau
    ArrayXr<Real> sin_t = sin(Rt_func->theta); // [T x 1]
    ArrayXr<Real> dx_dt_wt = x_theta * Rt_func->w_theta; // [T x 1]

    ArrayXr<Real> tmp1 = ArrayXr<Real>::LinSpaced(N_max + 1, 0, N_max); // [N x 1]
    // [N x 1] An (for prefactors)
    ArrayXr<Real> An_vec = sqrt((2*tmp1 + 1) / (2*tmp1 * (tmp1 + 1)));
    // [N x N] An_Ak matrix (for prefactors)
    ArrayXXr<Real> An_Ak = (An_vec.matrix() * An_vec.transpose().matrix()).array();

    // The hard part is to get the modified (non-problematic) radial functions.
    // Note that these do not depend on m.
    // The following function uses the algorithm of [JQSRT 2013] to achieve this
    unique_ptr<stBesselProducts<Real>> prods = sphGetModifiedBesselProducts(N_max, s, x, NB); // [N x K=N x T] tensors
    // pointers to the deepest members of prods are made for brevity of code
    Tensor3c<Real>* xi_prime_psi = &(prods->st_xi_psi_all->xi_prime_psi);
    Tensor3c<Real>* xi_psi_prime = &(prods->st_xi_psi_all->xi_psi_prime);
    Tensor3c<Real>* xi_psi = &(prods->st_xi_psi_all->xi_psi);
    Tensor3c<Real>* xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx =
        &(prods->st_xi_psi_all->xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx);
    Tensor3c<Real>* xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx =
        &(prods->st_xi_psi_all->xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx);
    ArrayXXc<Real>* for_Q_diag_Lt1 = &(prods->st_xi_psi_all->for_diag_Lt1);
    ArrayXXc<Real>* for_Q_diag_Lt2 = &(prods->st_xi_psi_all->for_diag_Lt2);
    ArrayXXc<Real>* for_Q_diag_Lt3 = &(prods->st_xi_psi_all->for_diag_Lt3);

    Tensor3c<Real>* psi_prime_psi = &(prods->st_psi_psi_all->xi_prime_psi);
    Tensor3c<Real>* psi_psi_prime = &(prods->st_psi_psi_all->xi_psi_prime);
    Tensor3c<Real>* psi_psi = &(prods->st_psi_psi_all->xi_psi);
    Tensor3c<Real>* psi_prime_psi_prime_plus_kkp1_psi_psi_over_sxx =
        &(prods->st_psi_psi_all->xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx);
    Tensor3c<Real>* psi_prime_psi_prime_plus_nnp1_psi_psi_over_sxx =
        &(prods->st_psi_psi_all->xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx);
    ArrayXXc<Real>* for_P_diag_Lt1 = &(prods->st_psi_psi_all->for_diag_Lt1);
    ArrayXXc<Real>* for_P_diag_Lt2 = &(prods->st_psi_psi_all->for_diag_Lt2);
    ArrayXXc<Real>* for_P_diag_Lt3 = &(prods->st_psi_psi_all->for_diag_Lt3);

    // The rest of the calculations are m-dependent, so we loop on m
    for (int i = 0; i < M; i++) {
      int m = abs_m_vec(i);
      int N_min = max(m, 1);
      int Nm = N_max - N_min + 1; // size of the matrices for a given m (since n,k >= m)
      ArrayXi n_vec = ArrayXi::LinSpaced(Nm, N_min, N_max); // [Nm X 1] a list of n values to use
      ArrayXi p_vec = n_vec*(n_vec + 1) + m; // [Nm X 1] p-index for positive m

      // angular functions
      ArrayXr<Real> n_vec_real = n_vec.template cast<Real>(); // [Nm X T]
      ArrayXXr<Real> pi_nm = stPT->pi_nm(all, p_vec).transpose(); // [Nm X T]
      ArrayXXr<Real> tau_nm = stPT->tau_nm(all, p_vec).transpose(); // [Nm X T]
      ArrayXXr<Real> d_n;
      if (m)
        d_n = pi_nm.rowwise() * (sin_t/m).transpose(); // [Nm X T]
      else
        d_n = stPT->p_n0(all, n_vec).transpose(); // [Nm X T]

      ArrayXXr<Real> d_n_times_nnp1 = d_n.colwise() * (n_vec_real*(n_vec_real + 1)); // [Nm X T]

      ArrayXXc<Real> K1, K2, L5, L6, K1P, K2P, L5P, L6P;
      K1 = K2 = L5 = L6 = K1P = K2P = L5P = L6P = ArrayXXc<Real>::Zero(Nm, Nm);

      // Note that the following code does not take into account the zeros in half of the terms
      // in the matrices and hence calculates twice as many terms as needed. This is because
      // in the original MATLAB code, this would've been as fast as more optimised code.
      // However, it may or may not speed up the code in C++ (it's possible that the zeros get optimised out
      // during compile time). Thus, for future speed ups, it may be worth trying to optimise this code.

      // The integrals are calculated as sums using Gaussian quadratures.
      // This is carried out for a given index k for all n (i.e. for a given column) by performing a matrix
      // multiplication of a [N X T] matrix by a [T X 1] vector. We therefore loop on k.
      for (int k = N_min; k <= N_max; k++) {
        // These parameters are used for truncating some Eigen::Tensors down to Eigen::Arrays
        int rows = Nm;
        int cols = Nb_theta;
        std::array<int, 3> offsets = {N_min, k, 0};
        std::array<int, 3> extents = {rows, 1, cols};

        int k_ind = k - N_min;
        ArrayXr<Real> d_k = d_n.row(k_ind).transpose();
        ArrayXr<Real> tau_k = tau_nm.row(k_ind).transpose();

        VectorXr<Real> dx_dt_tau_k_sin_t = dx_dt_wt * tau_k; // [T X 1]
        VectorXr<Real> dx_dt_d_k_sin_t = dx_dt_wt * d_k; // [T X 1]

        // For matrices K1 and K2, use eqs. 11/53 and 12/54 of JQSRT2013
        MatrixXc<Real> pi_n_xi_prime_psi = pi_nm *
            Subtensor2ArrMap(*xi_prime_psi, offsets, extents, rows, cols); // [Nm X T]
        MatrixXc<Real> pi_n_xi_psi_prime = pi_nm *
            Subtensor2ArrMap(*xi_psi_prime, offsets, extents, rows, cols); // [Nm X T]
        // The integrals are carried out as sums over theta's by taking a matrix product of
        // a [Nm X T] matrix by a [T X 1] vector
        K1.col(k_ind) = pi_n_xi_psi_prime * dx_dt_d_k_sin_t;
        K2.col(k_ind) = pi_n_xi_prime_psi * dx_dt_d_k_sin_t;

        // These are for Q, we also do the same for P. For efficiency, the same variable names are used
        // even if xi is replaced by psi
        pi_n_xi_prime_psi = pi_nm * Subtensor2ArrMap(*psi_prime_psi, offsets, extents, rows, cols);
        pi_n_xi_psi_prime = pi_nm * Subtensor2ArrMap(*psi_psi_prime, offsets, extents, rows, cols);
        K1P.col(k_ind) = pi_n_xi_psi_prime * dx_dt_d_k_sin_t;
        K2P.col(k_ind) = pi_n_xi_prime_psi * dx_dt_d_k_sin_t;

        // For L5, eqs. 18/52 is used
        MatrixXc<Real> d_n_xi_psi_nnp1 = d_n_times_nnp1 * Subtensor2ArrMap(*xi_psi, offsets, extents, rows, cols);
        MatrixXc<Real> tau_n_xi_psi = tau_nm * Subtensor2ArrMap(*xi_psi, offsets, extents, rows, cols);

        L5.col(k_ind) = d_n_xi_psi_nnp1 * dx_dt_tau_k_sin_t - tau_n_xi_psi * dx_dt_d_k_sin_t*k*(k + 1);

        // Same for P
        d_n_xi_psi_nnp1 = d_n_times_nnp1 * Subtensor2ArrMap(*psi_psi, offsets, extents, rows, cols);
        tau_n_xi_psi = tau_nm *  Subtensor2ArrMap(*psi_psi, offsets, extents, rows, cols);

        L5P.col(k_ind) = d_n_xi_psi_nnp1 * dx_dt_tau_k_sin_t - tau_n_xi_psi * dx_dt_d_k_sin_t*k*(k + 1);

        // For L6, eqs. 20-22/55-56 are used
        MatrixXc<Real> d_n_xi_prime_psi_prime_nnp1_plus_xi_psi_over_sxx_nnp1_kkp1 = d_n_times_nnp1 *
            Subtensor2ArrMap(*xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx, offsets, extents, rows, cols);
        MatrixXc<Real> tau_n_xi_prime_psi_prime_plus_xi_psi_over_sxx_nnp1 = tau_nm *
            Subtensor2ArrMap(*xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx, offsets, extents, rows, cols);

        L6.col(k_ind) = d_n_xi_prime_psi_prime_nnp1_plus_xi_psi_over_sxx_nnp1_kkp1 * dx_dt_tau_k_sin_t -
            tau_n_xi_prime_psi_prime_plus_xi_psi_over_sxx_nnp1 * dx_dt_d_k_sin_t * k*(k + 1);

        // Same for P
        d_n_xi_prime_psi_prime_nnp1_plus_xi_psi_over_sxx_nnp1_kkp1 = d_n_times_nnp1 *
            Subtensor2ArrMap(*psi_prime_psi_prime_plus_kkp1_psi_psi_over_sxx, offsets, extents, rows, cols);
        tau_n_xi_prime_psi_prime_plus_xi_psi_over_sxx_nnp1 = tau_nm *
            Subtensor2ArrMap(*psi_prime_psi_prime_plus_nnp1_psi_psi_over_sxx, offsets, extents, rows, cols);

        L6P.col(k_ind) = d_n_xi_prime_psi_prime_nnp1_plus_xi_psi_over_sxx_nnp1_kkp1 * dx_dt_tau_k_sin_t -
            tau_n_xi_prime_psi_prime_plus_xi_psi_over_sxx_nnp1 * dx_dt_d_k_sin_t * k*(k + 1);
      }

      // Prefactors
      ArrayXXc<Real> tmp2 = (n_vec_real * (n_vec_real + 1)).replicate(1, Nm).rowwise() -
      (n_vec_real * (n_vec_real + 1)).transpose();
      ArrayXXc<Real> prefactor1 = (s - 1)*(s + 1) / s * An_Ak(n_vec, n_vec); // [Nm X Nm] (s^2 - 1)/s * An*Ak
      ArrayXXc<Real> prefactor2 = mp_im_unit<Real>() * prefactor1 / tmp2; // [Nm X Nm] n(n + 1) - k(k + 1)

      // Get Q matrix with prefactors
      ArrayXXc<Real> Q12 = prefactor1 * K1; // Eq. 15
      ArrayXXc<Real> Q21 = -prefactor1 * K2; // Eq. 16
      ArrayXXc<Real> Q11 = prefactor2 * L5; // Eq. 17
      ArrayXXc<Real> Q22 = prefactor2 * L6; // Eq. 19

      // Get P matrix with prefactors
      ArrayXXc<Real> P12 = prefactor1 * K1P; // Eq. 15
      ArrayXXc<Real> P21 = -prefactor1 * K2P; // Eq. 16
      ArrayXXc<Real> P11 = prefactor2 * L5P; // Eq. 17
      ArrayXXc<Real> P22 = prefactor2 * L6P; // Eq. 19

      // To finish, we need to do the diagonal terms of Q11/P11 and Q22/P22
      // calculated differently to avoid problems
      ArrayXc<Real> prefact_diag1 = -mp_im_unit<Real>()/s * (2*n_vec_real + 1) / (2*n_vec_real*(n_vec_real + 1));
      ArrayXc<Real> prefact_diag2 = -mp_im_unit<Real>() *
          static_cast<complex<Real>>((s - 1)*(s + 1)/(2*s)) * (2*n_vec_real + 1);
      ArrayXXc<Real> pi2_p_tau2 = (pi_nm.pow(2) + tau_nm.pow(2));

      // Fill in the diagonals by calculating the integrals as matrix products of [Nm X T] by [T X 1] as before

      // Diagonal of Q11, from eqs. 23, 25, 65 (all Nm of them)
      ArrayXc<Real> L_tilde1 = ((pi2_p_tau2 * (*for_Q_diag_Lt1)(n_vec, all)).matrix()
          * Rt_func->w_theta.matrix()).array();
      // Diagonal of Q22, from eqs. 24, 26, 27 (all Nm of them)
      ArrayXc<Real> L_tilde2 = ((pi2_p_tau2 * (*for_Q_diag_Lt2)(n_vec, all)).matrix()
          * Rt_func->w_theta.matrix()).array();
      ArrayXc<Real> L_tilde3 = ((d_n * tau_nm * (*for_Q_diag_Lt3)(n_vec, all)).matrix() * dx_dt_wt.matrix()).array();

      // Fill in the diagonals of Q
      for (int j = 0; j < Nm; j++) {
        Q11(j, j) = prefact_diag1(j) * L_tilde1(j);
        Q22(j, j) = prefact_diag1(j) * L_tilde2(j) + prefact_diag2(j) * L_tilde3(j);
      }

      // Diagonal of P11, from eqs. 23, 25, 65 (all Nm of them)
      L_tilde1 = (pi2_p_tau2 * (*for_P_diag_Lt1)(n_vec, all)).matrix() * Rt_func->w_theta.matrix();
      // Diagonal of P22, from eqs. 24, 26, 27 (all Nm of them)
      L_tilde2 = (pi2_p_tau2 * (*for_P_diag_Lt2)(n_vec, all)).matrix() * Rt_func->w_theta.matrix();
      L_tilde3 = (d_n * tau_nm * (*for_P_diag_Lt3)(n_vec, all)).matrix() * dx_dt_wt.matrix();

      // Fill in the diagonals of P
      for (int j = 0; j < Nm; j++) {
        P11(j, j) = prefact_diag1(j) * L_tilde1(j);
        P22(j, j) = prefact_diag1(j) * L_tilde2(j) + prefact_diag2(j) * L_tilde3(j);
      }

      // Write results to output
      // if N_min is even, then the indices 0, 2, 4, .. are even. Otherwise, 1, 3, 5... are even.
      ArrayXi inde = Seq2Array(N_min % 2, Nm - 1, 2);
      ArrayXi indo = Seq2Array(1 - N_min % 2, Nm - 1, 2);

      output[i]->st_4M_Q_eo().M12 = Q12(inde, indo);
      output[i]->st_4M_Q_eo().M21 = Q21(indo, inde);
      output[i]->st_4M_Q_eo().M11 = Q11(inde, inde);
      output[i]->st_4M_Q_eo().M22 = Q22(indo, indo);
      output[i]->st_4M_Q_eo().m = m;
      output[i]->st_4M_Q_eo().ind1 = inde;
      output[i]->st_4M_Q_eo().ind2 = indo;

      output[i]->st_4M_Q_oe().M12 = Q12(indo, inde);
      output[i]->st_4M_Q_oe().M21 = Q21(inde, indo);
      output[i]->st_4M_Q_oe().M11 = Q11(indo, indo);
      output[i]->st_4M_Q_oe().M22 = Q22(inde, inde);
      output[i]->st_4M_Q_oe().m = m;
      output[i]->st_4M_Q_oe().ind1 = indo;
      output[i]->st_4M_Q_oe().ind2 = inde;

      output[i]->st_4M_P_eo().M12 = P12(inde, indo);
      output[i]->st_4M_P_eo().M21 = P21(indo, inde);
      output[i]->st_4M_P_eo().M11 = P11(inde, inde);
      output[i]->st_4M_P_eo().M22 = P22(indo, indo);
      output[i]->st_4M_P_eo().m = m;
      output[i]->st_4M_P_eo().ind1 = inde;
      output[i]->st_4M_P_eo().ind2 = indo;

      output[i]->st_4M_P_oe().M12 = P12(indo, inde);
      output[i]->st_4M_P_oe().M21 = P21(inde, indo);
      output[i]->st_4M_P_oe().M11 = P11(indo, indo);
      output[i]->st_4M_P_oe().M22 = P22(inde, inde);
      output[i]->st_4M_P_oe().m = m;
      output[i]->st_4M_P_oe().ind1 = indo;
      output[i]->st_4M_P_oe().ind2 = inde;

      output[i]->mat_list.push_back("st_4M_Q");
      output[i]->mat_list.push_back("st_4M_P");
    }
    return output;
  }
}

#endif
