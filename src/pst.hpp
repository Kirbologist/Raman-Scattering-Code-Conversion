/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Based on the SMARTIES MATLAB package by W.R.C. Somerville, B. Auguié, E.C. Le Ru
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.


This code contains all 'pst' SMARTIES functions that are used in Raman scattering calculations,
i.e. high-level functions used for post-processing. Much of the code here is based on the Fortran code in tmd.lp.f
written by Mishcehnko and co-workers publicly available at http://www.giss.nasa.gov/staff/mmishchenko/t_matrix.html
*/

#ifndef PST_HPP
#define PST_HPP

#include "core.hpp"
#include "rvh.hpp"
#include "vsh.hpp"
#include "misc.hpp"

using namespace Eigen;
using namespace std;

namespace Smarties {

  template <class Real>
  struct stRes {
    ArrayXc<Real> p_nm;
    ArrayXc<Real> q_nm;
    ArrayXc<Real> a_nm;
    ArrayXc<Real> b_nm;
    ArrayXc<Real> c_nm;
    ArrayXc<Real> d_nm;
    int N_max;
    ArrayXr<Real> lambda;
    Real epsilon1;
    ArrayXr<Real> epsilon2;
    unique_ptr<stIncPar<Real>> inc_par;
    Real a;
    Real c;

    stRes() {};

    stRes(const stAbcdnm<Real>& base) {
      this->p_nm = base.p_nm;
      this->q_nm = base.q_nm;
      this->a_nm = base.a_nm;
      this->b_nm = base.b_nm;
      this->c_nm = base.c_nm;
      this->d_nm = base.d_nm;
    }
  };

  /*
  L_max -
  ALF1n, ALF2n, ALF3n, ALF4n, BET1n, BET2n [N + 1 X 1] -
      expansion coefficients for the scattering matrix as defined in eqs. 4.75-4.80 and obtained from
      eqs. 4.109-4.114 of Mishchenko 2002. Note that those start at n = 0 and are zero for n > L_max
  asym_par - Asymmetry parameter (<cos theta>)
  theta, theta
  */
  template <class Real>
  struct stSM {
    int L_max; // maximum multipole order for expansion (called s_max in eqs. 4.75-4.80 of Mishchenko 2002)
    // ALFXn and BETXn are expansion coefficients for the scattering matrix as defined in eqs. 4.75-4.80 and
    // obtained from eqs. 4.109-4.114 of Mishchenko 2002. Note that those start at n = 0 and are zero for n > L_max.
    ArrayXr<Real> ALF1n;
    ArrayXr<Real> ALF2n;
    ArrayXr<Real> ALF3n;
    ArrayXr<Real> ALF4n;
    ArrayXr<Real> BET1n;
    ArrayXr<Real> BET2n;
    Real asym_par; // Asymmetry parameter (<cos theta>)
    ArrayXr<Real> theta; // Angles in radians for the scattering matrix calculation
    ArrayXr<Real> theta_deg; // Angles in degrees for the scattering matrix calculation
    // FXX are scattering matrix elements as defined in Mishchenko 2002.
    ArrayXr<Real> F11;
    ArrayXr<Real> F22;
    ArrayXr<Real> F33;
    ArrayXr<Real> F44;
    ArrayXr<Real> F12;
    ArrayXr<Real> F34;
    Array<Real, Dynamic, 7> all_AB; // Contains all the expansion coefficients
    Array<Real, Dynamic, 7> all_F; // Contains all the angles and scattering matrix elements
  };


  /*
  Creates struct with the necessary parameters to calculate surface fields. An alternate template definition
  is given where most of the input parameters are given through a single stParams struct.
  Inputs:
    st_abdcnm - The field expansion coefficients (from rvhGetFieldCoefficients).
                This may be for one or several wavelengths.
    N_max - The maximum multipole order N
    lambda - The values of wavelengths used.
    epsilon2 - The values of the dielectric functions in the particle
    epsilon1 - The dielectric function of the surrounding medium
    stIncPar - Struct defining the incident plane wave, (from vshMakeIncidentParams)
    a - radius of semi-axis of spheroid along x,y
    c - radius of semi-axis of spheroid along z

  Output:
    Returns an stRes struct which combines the members of st_abcdnm along with all other parameters.
  */
  template <class Real>
  unique_ptr<stRes<Real>> pstMakeStructForField(
      const unique_ptr<stAbcdnm<Real>>& st_abcdnm, int N_max, ArrayXr<Real> lambda,
      ArrayXr<Real> epsilon2, Real epsilon1, unique_ptr<stIncPar<Real>> inc_par,
      Real a = numeric_limits<Real>::quiet_NaN(), Real c = numeric_limits<Real>::quiet_NaN()) {
    if (!st_abcdnm->c_nm.size())
      cout << "Pb in  pstMakeStructForField: " <<
          "the structure stAbcdnm should contain c_nm and d_nm for internal fields" << endl;
    auto output = make_unique<stRes<Real>>(*st_abcdnm);
    output->N_max = N_max;
    output->lambda = lambda;
    output->epsilon1 = epsilon1;
    output->epsilon2 = epsilon2;
    output->inc_par = make_unique<stIncPar<Real>>(*inc_par);
    if (!isnan(a))
      output->a = a;
    if (!isnan(c))
      output->a = c;
    return output;
  }

  /*
  Creates struct with the necessary parameters to calculate surface fields. An alternate template definition
  is given where the relevant members of stParams are given directly as parameters.
  Inputs:
    st_abdcnm - The field expansion coefficients (from rvhGetFieldCoefficients).
                This may be for one or several wavelengths.
    params - Struct containing the rest of the relevant parameters.
             Only N, lambda, epsilon1, epsilon2, a, c need to be defined, with inc_par being optional
             and all others being unused.

  Output:
    Returns an stRes struct which combines the members of st_abcdnm along with all other parameters.
  Dependencies:
    vshMakeIncidentParameters
  */
  template <class Real>
  unique_ptr<stRes<Real>> pstMakeStructForField(
      const unique_ptr< stAbcdnm<Real>>& st_abcdnm, const unique_ptr<stParams<Real>>& params) {
    if (!st_abcdnm->c_nm.size())
      cout << "Pb in  pstMakeStructForField: " <<
          "the structure stAbcdnm should contain c_nm and d_nm for internal fields" << endl;
    auto output = make_unique<stRes<Real>>(*st_abcdnm);
    output->N_max = params->N;
    output->lambda = params->lambda;
    output->epsilon1 = params->epsilon1;
    output->epsilon2 = params->epsilon2;
    output->a = params->a;
    output->c = params->c;
    if (params->inc_par)
      output->inc_par = make_unique<stIncPar<Real>>(*(params->inc_par));
    else {
      sIncType inc_type = params->inc_type;
      output->inc_par = vshMakeIncidentParams<Real>(inc_type, params->N);
    }
    return output;
  }

  /*
  Subroutine used in pstScatteringMatrixOA.
  Calculates Clebcsh-Gordan coefficients for given n, n1, m, mm
  Based on the Fortran subroutine of the same name in tmd.lp.f
  */
  template <class Real>
  Real CCGIN(int n, int n1, int m, int mm, ArrayXr<Real> F, ArrayXi s_sign) {
    int m1 = mm - m;
    if (n < abs(m) || n1 < abs(m1) || abs(mm) > n + n1) {
      cout << "Error in CCGIN" << endl;
      return 0;
    }
    if (abs(mm) <= abs(n - n1)) {
      if (n1 > n) {
        swap(n, n1);
        swap(m, m1);
      }
      int N2 = 2*n;
      int N12 = 2*n1;
      return s_sign(n1 + m1) * exp(F(n + m) + F(n - m) + F(N12) + F(N2 - N12 + 1)
          - F(N2 + 1) - F(n1 + m1) - F(n1 - m1) - F(n - n1 + mm) - F(n - n1 - mm));
    }
    // abs(mm) > abs(n - n1)
    int A = 1;
    if (mm < 0) {
      mm = -mm;
      m = -m;
      m1 = -m1;
      A = s_sign(mm + n + n1);
    }
    return A * s_sign(n + m) * exp(F(2 * mm + 1) + F(n + n1 - mm) + F(n + m) + F(n1 + m1)
        - F(n + n1 + mm + 1) - F(n - n1 + mm) - F(-n + n1 + mm) - F(n - m) - F(n1 - m1));
  }

  /*
  Subroutine used in pstScatteringMatrixOA.
  Based on the Fortran subroutine DIRECT in tmd.lp.f
  */
  template <class Real>
  Real CGDIRECT(int n, int m, int n1, int m1, ArrayXr<Real> F) {
    Real C = F(2*n) + F(2*n1) + F(n + n1 + m + m1) + F(n + n1 - m - m1);
    C -= F(2*(n + n1)) + F(n + m) + F(n - m) + F(n1 + m1) + F(n1 - m1);
    return exp(C);
  }

  /*
  Subroutine used in pstScatteringMatrixOA.
  Calculates Clebsch-Gordan Coefficients for given n and n1, based on the Fortran subroutine
  of the same name in tmd.lp.f
  Dependencies:
    CCGIN, CGDIRECT
  */
  template <class Real>
  ArrayXXr<Real> CCG(int n, int n1, int N_max, int K1, int K2, ArrayXr<Real> log_fact, ArrayXi s_sign) {
    if (n1 < 0 || n1 > N_max + n || n < 0 || n > N_max) {
      cout << "Error in CCG" << endl;
      return ArrayXXr<Real>::Zero(1, 1);
    }
    ArrayXXr<Real> GG = ArrayXXr<Real>::Zero(2*N_max + 1, N_max + 1);
    ArrayXXr<Real> CD = ArrayXXr<Real>::Zero(2*N_max + 1, 2*N_max + 1);
    ArrayXXr<Real> CU = ArrayXXr<Real>::Zero(2*N_max + 1, 2*N_max + 1);

    int NPN6 = N_max;
    int NNF = min(n + n1, N_max);
    int m_init = (K1 == 1 && K2 == 0) ? NPN6 : NPN6 - n;
    int m_fin = NPN6 + n;

    for (int m_ind = m_init; m_ind <= m_fin; m_ind++) {
      int m = m_ind - NPN6;
      int mm = m * K1 + K2;
      int m1 = mm - m;

      int NNL = max(abs(mm), abs(n - n1));
      if (abs(m1) > n1 || NNL > NNF)
        continue;
      int NNU = n + n1;
      int NNM = (NNU + NNL) / 2;
      if (NNU == NNL)
        NNM = NNL;
      Real C = CCGIN(n, n1, m, mm, log_fact, s_sign);
      CU(NNL) = C;
      if (NNL != NNF) {
        Real C2 = 0;
        Real C1 = C;
        for (int nn = NNL + 1; nn <= min(NNM, NNF); nn++) {
          Real A = static_cast<Real>((nn + mm) * (nn - mm) * (n1 - n + nn));
          A *= (n - n1 + nn) * (n + n1 - nn + 1) * (n + n1 + nn + 1);
          A = (4 * nn * nn) / A;
          A *= (2*nn + 1) * (2*nn - 1);
          A = sqrt(A);
          Real B, D;
          if (nn == 1) {
            B = 0.5*(m - m1);
            D = 0;
          } else {
            B = 2 * nn * (nn - 1);
            B = ((2*m - mm) * nn * (nn - 1) - mm * n * (n + 1) + mm * n1 * (n1 + 1)) / B;
            D = 4 * (nn - 1) * (nn - 1);
            D *= (2*nn - 3) * (2*nn - 1);
            D = ((nn - mm - 1) * (nn + mm - 1) * (n1 - n + nn - 1)) / D;
            D *= (n - n1 + nn - 1) * (n + n1 - nn + 2) * (n + n1 + nn);
            D = sqrt(D);
          }
          C = A * (B * C1 - D * C2);
          C2 = C1;
          C1 = C;
          CU(nn) = C;
        }
        if (NNF > NNM) {
          C = CGDIRECT(n, m, n1, m1, log_fact);
          CD(NNU) = C;
          if (NNU != NNM + 1) {
            C2 = 0;
            C1 = C;
            for (int nn = (NNU - 1); nn > NNM; nn--) {
              Real A = static_cast<Real>((nn - mm + 1) * (nn + mm + 1) * (n1 - n + nn + 1));
              A *= ((n - n1 + nn + 1) * (n + n1 - nn) * (n + n1 + nn + 2));
              A = (4 * (nn + 1) * (nn + 1)) / A;
              A *= (2*nn + 1) * (2*nn + 3);
              A = sqrt(A);
              Real B = static_cast<Real>(2 * (nn + 2) * (nn + 1));
              B = ((2 * m - mm) * (nn + 2) * (nn + 1) - mm * n * (n + 1) + mm * n1 * (n1 + 1)) / B;
              Real D = static_cast<Real>(4 * (nn + 2) * (nn + 2));
              D *= (2 * nn + 5) * (2 * nn + 3);
              D = ((nn + mm + 2) * (nn - mm + 2) * (n1 - n + nn + 2)) / D;
              D *= (n - n1 + nn + 2) * ( n + n1 - nn - 1) * (n + n1 + nn + 3);
              D = sqrt(D);
              C = A * (B * C1 - D * C2);
              C2 = C1;
              C1 = C;
              CD(nn) = C;
            }
          }
        }
      }

      for (int nn = NNL; nn <= NNF; nn++) {
        if (nn <= NNM)
          GG(m_ind, nn) = CU(nn);
        else
          GG(m_ind, nn) = CD(nn);
      }
    }
    return GG;
  }

  /*
  Calculates the quantities related to the scattering matrix of an ensemble of
  randomly-oriented scatterers.
  Based on the Fortran code by Mischenko and co-workers. Variable names and loop structures are kept from there.
  The method consists in computing first a number of expansion coefficients
  alpha_n and beta_n up to some maximum n<=LMAX, from which the scattering
  matrix elements F11 etc... can be computed accurately at any angle in a second stage.
  Inputs:
    st_TR_list - vector containing unique pointers to TR structs, constituating the T-matrix for all m.
    lambda - wavelength (in the same unit as C_sca)
    C_sca - orientation-averaged scattering cross-section. Note that C_sca and lambda are only needed
            to get the correct scaling factor for the scattering matrix
    Nb_theta - number of angles for which the scattering matrix is computed
               (the angles are linearly spaced between 0 and pi).
               This parameter is the same as NPNA in the Fortran codes.
  Output:
    Returns a unique pointer to a stSM struct (which are similar to the results returned by the Fortran codes).
    ALFXn and BETXn are [(N + 1) X 1] arrays, theta and theta_deg are [Nb_theta X 1] arrays,
    all_AB is a [(N + 1) X 7] array and all_F is a [(Nb_theta) X 7] array.
  Dependencies:
    CCG
  */
  template <class Real>
  unique_ptr<stSM<Real>> pstScatteringMatrixOA(const vector<unique_ptr<stTR<Real>>>& st_TR_list,
      Real lambda, Real C_sca, int Nb_theta = 2) {

    // Calculate arrays B1 and B2
    int K1 = 1, K2 = 0, K3 = 0, K4 = 1, K5 = 1, K6 = 2;
    int N1 = st_TR_list.size();
    int N = N1 - 1; // assumes that all m of T-matrix have been calculated

    RowArrayXr<Real> n_vec = RowArrayXr<Real>::LinSpaced(N1, 0, N);
    ArrayXr<Real> k_vec = ArrayXr<Real>::LinSpaced(N1, 0, N);

    // Encodes SSI(n) and SSJ(n) from the Fortran codes
    ArrayXXr<Real> FF_kn = sqrt((2*n_vec + 1).replicate(N1, 1).colwise() / (2*k_vec + 1));
    // Factor for A
    ArrayXXc<Real> FA_kn = pow(mp_im_unit<Real>(), -n_vec).replicate(N1, 1).colwise() *
        (pow(mp_im_unit<Real>(), k_vec) / sqrt(2*k_vec + 1));
    // Log of factorials array
    ArrayXr<Real> log_fact = ArrayXr<Real>::Zero(4*N1);
    for (int i = 2; i < 4*N1; i++)
      log_fact(i) = log_fact(i - 1) + 0.5*log(i);

    // Sign array (-1)^n
    ArrayXi s_sign = pow(-1, ArrayXi::LinSpaced(4*N1, 0, 4*N1 - 1));

    ArrayXXc<Real> T1_mk = ArrayXXc<Real>::Zero(2*N + 1, N1);
    ArrayXXc<Real> T2_mk = ArrayXXc<Real>::Zero(2*N + 1, N1);
    // The following tensors originally initialised to different sizes in the MATLAB version.
    // Those sizes are actually too small; MATLAB was dynamically resizing the tensor to fit.
    // Cases where n = 0 is included as padding to make the indexing more convenient
    Tensor3c<Real> B1_n1_mn(2*N + 1, 2*N + 1, N1);
    B1_n1_mn.setZero();
    Tensor3c<Real> B2_n1_mn(2*N + 1, 2*N + 1, N1);
    B2_n1_mn.setZero();

    // Get T-matrix in normal form
    vector<ArrayXXc<Real>> CT11(N1);
    vector<ArrayXXc<Real>> CT12(N1);
    vector<ArrayXXc<Real>> CT21(N1);
    vector<ArrayXXc<Real>> CT22(N1);

    for (int m = 0; m <= N; m++) {
      int N_min_m1 = max(1, m);
      CT11[m] = ArrayXXc<Real>::Zero(N1, N1);
      CT12[m] = ArrayXXc<Real>::Zero(N1, N1);
      CT21[m] = ArrayXXc<Real>::Zero(N1, N1);
      CT22[m] = ArrayXXc<Real>::Zero(N1, N1);
      for (int i = 0; i < st_TR_list[m]->st_4M_T_eo().ind1.size(); i++) {
        for (int j = 0; j < st_TR_list[m]->st_4M_T_eo().ind1.size(); j++)
          CT11[m](N_min_m1 + st_TR_list[m]->st_4M_T_eo().ind1(i),
              N_min_m1 + st_TR_list[m]->st_4M_T_eo().ind1(j)) = st_TR_list[m]->st_4M_T_eo().M11(i, j);
        for (int j = 0; j < st_TR_list[m]->st_4M_T_eo().ind2.size(); j++)
          CT12[m](N_min_m1 + st_TR_list[m]->st_4M_T_eo().ind1(i),
              N_min_m1 + st_TR_list[m]->st_4M_T_eo().ind2(j)) = st_TR_list[m]->st_4M_T_eo().M12(i, j);
      }
      for (int i = 0; i < st_TR_list[m]->st_4M_T_eo().ind2.size(); i++) {
        for (int j = 0; j < st_TR_list[m]->st_4M_T_eo().ind1.size(); j++)
          CT21[m](N_min_m1 + st_TR_list[m]->st_4M_T_eo().ind2(i),
              N_min_m1 + st_TR_list[m]->st_4M_T_eo().ind1(j)) = st_TR_list[m]->st_4M_T_eo().M21(i, j);
        for (int j = 0; j < st_TR_list[m]->st_4M_T_eo().ind2.size(); j++)
          CT22[m](N_min_m1 + st_TR_list[m]->st_4M_T_eo().ind2(i),
              N_min_m1 + st_TR_list[m]->st_4M_T_eo().ind2(j)) = st_TR_list[m]->st_4M_T_eo().M22(i, j);
      }
      for (int i = 0; i < st_TR_list[m]->st_4M_T_oe().ind1.size(); i++) {
        for (int j = 0; j < st_TR_list[m]->st_4M_T_oe().ind1.size(); j++)
          CT11[m](N_min_m1 + st_TR_list[m]->st_4M_T_oe().ind1(i),
              N_min_m1 + st_TR_list[m]->st_4M_T_oe().ind1(j)) = st_TR_list[m]->st_4M_T_oe().M11(i, j);
        for (int j = 0; j < st_TR_list[m]->st_4M_T_oe().ind2.size(); j++)
          CT12[m](N_min_m1 + st_TR_list[m]->st_4M_T_oe().ind1(i),
              N_min_m1 + st_TR_list[m]->st_4M_T_oe().ind2(j)) = st_TR_list[m]->st_4M_T_oe().M12(i, j);
      }
      for (int i = 0; i < st_TR_list[m]->st_4M_T_oe().ind2.size(); i++) {
        for (int j = 0; j < st_TR_list[m]->st_4M_T_oe().ind1.size(); j++)
          CT21[m](N_min_m1 + st_TR_list[m]->st_4M_T_oe().ind2(i),
              N_min_m1 + st_TR_list[m]->st_4M_T_oe().ind1(j)) = st_TR_list[m]->st_4M_T_oe().M21(i, j);
        for (int j = 0; j < st_TR_list[m]->st_4M_T_oe().ind2.size(); j++)
          CT22[m](N_min_m1 + st_TR_list[m]->st_4M_T_oe().ind2(i),
              N_min_m1 + st_TR_list[m]->st_4M_T_oe().ind2(j)) = st_TR_list[m]->st_4M_T_oe().M22(i, j);
      }
    }

    // Calculation of B1 and B2 mnk
    for (int n = 1; n <= N; n++) {
      // Calculation of T1 and T2 arrays
      for (int nn = 1; nn <= N; nn++) {
        int m_max = min(n, nn);
        for (int m = 0; m <= m_max; m++) {
          int m_ind = N + m;
          // This makes it eo-oe compatible
          complex<Real> TT1 = CT11[m](n, nn);
          complex<Real> TT2 = CT12[m](n, nn);
          complex<Real> TT3 = CT21[m](n, nn);
          complex<Real> TT4 = CT22[m](n, nn);
          complex<Real> T1 = TT1 + TT2;
          complex<Real> T2 = TT3 + TT4;
          T1_mk(m_ind, nn) = T1 + T2; // [2N + 1 X N]
          T2_mk(m_ind, nn) = T1 - T2;
          if (m > 0) { // Calculate for m < 0
            T1 = TT1 - TT2;
            T2 = TT3 - TT4;
            m_ind = N - m;
            T1_mk(m_ind, nn) = T1 - T2;
            T2_mk(m_ind, nn) = T1 + T2;
          }
        }
      }
      // Finished calculating arrays T1 and T2

      int nn1_max = N + n;
      for (int n1 = 0; n1 <= nn1_max; n1++) { // loop over nn1 = n1 + 1 in B
        // Calculate arrays A1 and A2
        ArrayXXr<Real> G1 = CCG(n, n1, N, K1, K2, log_fact, s_sign); // [2N + 1 X N + 1]
        int nn_max = min(N, n1 + n);
        int nn_min = max(1, abs(n - n1));
        int kn = n + n1;
        ArrayXc<Real> A1k = ArrayXc<Real>::Zero(N1);
        ArrayXc<Real> A2k = ArrayXc<Real>::Zero(N1);
        for (int nn = nn_min; nn <= nn_max; nn++) { // Loop over n' for A
          int SIG = s_sign(kn + nn); // (-1)^(kn + nn + 1)
          int m_ind_max = min(n, nn) + N;
          complex<Real> AA1 = static_cast<complex<Real>>(0);
          complex<Real> AA2 = static_cast<complex<Real>>(0);
          for (int m_ind = N; m_ind <= m_ind_max; m_ind++) { // Loop on m from 0 to min(n, nn)
            int m = m_ind - N;
            complex<Real> SSS = G1(m_ind, nn);
            complex<Real> R1 = T1_mk(m_ind, nn);
            complex<Real> R2 = T2_mk(m_ind, nn);
            if (m > 0) { // Calculate for m < 0
              int m_ind_n = N - m;
              R1 += T1_mk(m_ind_n, nn) * static_cast<complex<Real>>(SIG);
              R2 += T2_mk(m_ind_n, nn) * static_cast<complex<Real>>(SIG);
            }
            AA1 += SSS * R1;
            AA2 += SSS * R2;
          }
          A1k(nn) = AA1 * FA_kn(nn, n);
          A2k(nn) = AA2 * FA_kn(nn, n);
        }
        // Finished calculating arrays A1 and A2

        // K3 = 0, K4 = 1
        ArrayXXr<Real> G2 = CCG(n, n1, N, K3, K4, log_fact, s_sign); // [2N + 1 X N + 1]
        int m_max = min(n1 + 1, n);
        int m_min = max(-n1 + 1, -n);
        for (int m = m_min; m <= m_max; m++) {
          int m_ind = m + N;
          complex<Real> BB1 = static_cast<complex<Real>>(0);
          complex<Real> BB2 = static_cast<complex<Real>>(0);
          for (int nn = nn_min; nn <= nn_max; nn++) {
            complex<Real> SSS = G2(m_ind, nn);
            BB1 += SSS * A1k(nn);
            BB2 += SSS * A2k(nn);
          }
          B1_n1_mn(n1, m_ind, n) = BB1;
          B2_n1_mn(n1, m_ind, n) = BB2;
        }
      }
    }
    // Finished calculating arrays B1 and B2

    Tensor3r<Real> D1_m_kn(2*N + 1, N1, N1);
    Tensor3r<Real> D2_m_kn(2*N + 1, N1, N1);
    Tensor3r<Real> D3_m_kn(2*N + 1, N1, N1);
    Tensor3r<Real> D4_m_kn(2*N + 1, N1, N1);
    Tensor3c<Real> D5_m_kn(2*N + 1, N1, N1);
    D1_m_kn.setZero();
    D2_m_kn.setZero();
    D3_m_kn.setZero();
    D4_m_kn.setZero();
    D5_m_kn.setZero();

    // Calculate arrays D1, D2, D3, D4 and D5
    for (int n = 1; n <= N; n++) {
      for (int nn = 1; nn <= N; nn++) {
        int m_ind = min(n, nn);
        int m_ind_max = N + m_ind;
        int m_ind_min = N - m_ind;
        int n1_max = m_ind_max;
        for (m_ind = m_ind_min; m_ind <= m_ind_max; m_ind++) {
          int m = m_ind - N;
          int n1_min = abs(m - 1);
          Real DD1 = 0;
          Real DD2 = 0;
          for (int n1 = n1_min; n1 <= n1_max; n1++) { // Loop over n1 for sums of BB in D
            int XX = 2*n1 + 1;
            DD1 += XX * real(B1_n1_mn(n1, m_ind, n) * conj(B1_n1_mn(n1, m_ind, nn)));
            DD2 += XX * real(B2_n1_mn(n1, m_ind, n) * conj(B2_n1_mn(n1, m_ind, nn)));
          }
          D1_m_kn(m_ind, nn, n) = DD1; // Note the real part is not present in eqs. 5.132-5.136 of the book
          D2_m_kn(m_ind, nn, n) = DD2; // But it would give the same final results
        }
        int m_max = min(n, nn + 2);
        int m_min = max(-n, -nn + 2);
        m_ind_max = N + m_max;
        m_ind_min = N + m_min;
        for (int m_ind = m_ind_min; m_ind <= m_ind_max; m_ind++) {
          int m = m_ind - N;
          int n1_min = abs(m - 1);
          // Same comment as above regarding the real part
          Real DD3 = 0;
          Real DD4 = 0;
          complex<Real> DD5 = static_cast<complex<Real>>(0);
          int m_ind2 = -m + 2 + N;
          for (int n1 = n1_min; n1 <= n1_max; n1++) {
            Real XX = 2*n1 + 1;
            DD3 += XX * real(B1_n1_mn(n1, m_ind, n) * conj(B1_n1_mn(n1, m_ind2, nn)));
            DD4 += XX * real(B2_n1_mn(n1, m_ind, n) * conj(B2_n1_mn(n1, m_ind2, nn)));
            DD5 += XX * B2_n1_mn(n1, m_ind, n) * conj(B1_n1_mn(n1, m_ind2, nn));
          }
          D3_m_kn(m_ind, nn, n) = DD3;
          D4_m_kn(m_ind, nn, n) = DD4;
          D5_m_kn(m_ind, nn, n) = DD5;
        }
      }
    }
    // Finished calculating D-arrays

    auto output = make_unique<stSM<Real>>();
    output->ALF1n = ArrayXr<Real>::Zero(2*N + 1);
    output->ALF2n = ArrayXr<Real>::Zero(2*N + 1);
    output->ALF3n = ArrayXr<Real>::Zero(2*N + 1);
    output->ALF4n = ArrayXr<Real>::Zero(2*N + 1);
    output->BET1n = ArrayXr<Real>::Zero(2*N + 1);
    output->BET2n = ArrayXr<Real>::Zero(2*N + 1);

    // Calculate the expansion coefficients
    int L_max = 0;
    Real DK = pow(lambda, 2) / (4*C_sca*mp_pi<Real>());
    for (int l = 0; l < 2*N + 1; l++) { // loop on s (l here)
      Real G1L = 0;
      Real G2L = 0;
      Real G3L = 0;
      Real G4L = 0;
      complex<Real> G5L = static_cast<complex<Real>>(0);
      Real SL = (2*l + 1)*DK;
      for (int n = 1; n <= N; n++) { // sum from n = 0
        int nn_min = max(1, abs(n - l));
        int nn_max = min(N, n + l);
        if (nn_min <= nn_max) { // only sum the terms that fit this condition
          ArrayXXr<Real> G1 = CCG(n, l, N, K1, K2, log_fact, s_sign); // K1 = 1, K2 = 0
          ArrayXXr<Real> G2;
          if (l >= 2)
            G2 = CCG(n, l, N, K5, K6, log_fact, s_sign); // K5 = 1, K6 = 2
          for (int nn = nn_min; nn <= nn_max; nn++) { // sum over i (nn here)
            int m_max = min(n, nn);
            int m_ind_min = N - m_max;
            int m_ind_max = N + m_max;
            int SI = s_sign(n + l + nn);
            Real DM1 = 0;
            Real DM2 = 0;
            for (int m_ind = m_ind_min; m_ind <= m_ind_max; m_ind++) { // sum over m
              int m = m_ind - N;
              Real SSS1;
              if (m >= 0)
                SSS1 = G1(m_ind, nn);
              else {
                int m_ind_n = N - m;
                SSS1 = G1(m_ind_n, nn) * SI;
              }
              DM1 += SSS1 * D1_m_kn(m_ind, nn, n);
              DM2 += SSS1 * D2_m_kn(m_ind, nn, n);
            }
            Real FFN = FF_kn(nn, n); // Prefactors
            Real SSS = G1(N + 1, nn) * FFN; // G1(m = 1, nn)
            G1L += SSS * DM1;
            G2L += SSS * DM2 * SI;
            if (l >= 2) {
              Real DM3 = 0;
              Real DM4 = 0;
              complex<Real> DM5 = static_cast<complex<Real>>(0);
              int m_max = min(n, nn + 2);
              int m_min = max(-n, -nn + 2);
              int m_ind_max = N + m_max;
              int m_ind_min = N + m_min;
              for (int m_ind = m_ind_min; m_ind <= m_ind_max; m_ind++) {
                int m = m_ind - N;
                int m_ind_n = N - m;
                Real SSS1 = G2(m_ind_n, nn);
                DM3 += SSS1 * D3_m_kn(m_ind, nn, n);
                DM4 += SSS1 * D4_m_kn(m_ind, nn, n);
                DM5 += SSS1 * D5_m_kn(m_ind, nn, n);
              }
              G5L -= SSS * DM5;
              SSS = G2(N - 1, nn) * FFN;
              G3L += SSS * DM3;
              G4L += SSS * DM4 * SI;
            } // end if (l >= 2)
          } // finished loop on nn
        } // end if (nn_min <= nn_max)
      } // finished loop on n
      G1L *= SL;
      G2L *= SL;
      G3L *= SL;
      G4L *= SL;
      G5L *= SL;
      output->ALF1n(l) = G1L + G2L;
      output->ALF2n(l) = G3L + G4L;
      output->ALF3n(l) = G3L - G4L;
      output->ALF4n(l) = G1L - G2L;
      output->BET1n(l) = real(G5L) * 2;
      output->BET2n(l) = imag(G5L) * 2;

      // Keep more orders than in the Fortran code for maximum precision
      L_max = l;
      if (abs(G1L) < 1e-15)
        break;
    }
    output->L_max = L_max;
    output->asym_par = output->ALF1n(1) / 3; // asymmetry parameter (<cos theta>)

    // Calculate the scattering matrix for given expansion coefficients
    ArrayXr<Real> theta = ArrayXr<Real>::LinSpaced(Nb_theta, 0, mp_pi<Real>());
    ArrayXr<Real> theta_deg = theta * 180 / mp_pi<Real>();
    output->theta = theta;
    output->theta_deg = theta_deg;
    output->F11 = ArrayXr<Real>::Zero(Nb_theta);
    output->F22 = ArrayXr<Real>::Zero(Nb_theta);
    output->F33 = ArrayXr<Real>::Zero(Nb_theta);
    output->F44 = ArrayXr<Real>::Zero(Nb_theta);
    output->F12 = ArrayXr<Real>::Zero(Nb_theta);
    output->F34 = ArrayXr<Real>::Zero(Nb_theta);

    Real DN = static_cast<Real>(1.0)/(Nb_theta - 1);
    Real DA = DN * mp_pi<Real>(); // theta
    Real DB = DN * 180; // theta_deg
    Real TB = -DB;
    Real TAA = -DA;
    Real D6 = sqrt(static_cast<Real>(6.0)) / 4;

    for (int I1 = 0; I1 < Nb_theta; I1++) { // Loop over scatering angles theta
      TAA += DA; // theta(I1)
      TB += DB; // theta_deg(I1)
      Real U = cos(TAA); // u = cos(theta)
      Real F11, F2, F3, F44, F12, F34, P1, P2, P3, P4;
      F11 = F2 = F3 = F44 = F12 = F34 = P1 = P2 = P3 = P4 = 0;
      Real PP1 = 1;
      Real PP2 = pow(1 + U, 2)/4;
      Real PP3 = pow(1 - U, 2)/4;
      Real PP4 = D6 * (pow(U, 2) - 1);
      for (int L = 0, L1 = 1; L <= L_max; L++, L1++) {
        F11 += output->ALF1n(L) * PP1;
        F44 += output->ALF4n(L) * PP1;
        int PL1;
        if (L != L_max) {
          PL1 = 2*L + 1;
          Real P = (PL1 * U * PP1 - L * P1) / L1;
          P1 = PP1;
          PP1 = P;
        }
        if (L >= 2) {
          F2 += (output->ALF2n(L) + output->ALF3n(L)) * PP2;
          F3 += (output->ALF2n(L) - output->ALF3n(L)) * PP3;
          F12 += output->BET1n(L) * PP4;
          F34 += output->BET2n(L) * PP4;
          if (L != L_max) {
            Real PL2 = L * L1 * U;
            int PL3 = L1 * (L*L - 4);
            Real PL4 = static_cast<Real>(1.0) / (L*(L1*L1 - 4));
            Real P = (PL1 * (PL2 - 4) * PP2 - PL3 * P2) * PL4;
            P2 = PP2;
            PP2 = P;
            P = (PL1 * (PL2 + 4) * PP3 - PL3 * P3) * PL4;
            P3 = PP3;
            PP3 = P;
            P = (PL1 * U * PP4 - sqrt(static_cast<Real>(L*L) - 4) * P4) / sqrt(static_cast<Real>(L1*L1) - 4);
            P4 = PP4;
            PP4 = P;
          }
        }
      } // finished looping over L
      Real F22 = (F2 + F3) / 2;
      Real F33 = (F2 - F3) / 2;
      output->F11(I1) = F11;
      output->F22(I1) = F22;
      output->F33(I1) = F33;
      output->F44(I1) = F44;
      output->F12(I1) = F12;
      output->F34(I1) = F34;
    } // finished looping over I1

    output->all_AB = Array<Real, Dynamic, 7>(2*N + 1, 7);
    output->all_AB << ArrayXr<Real>::LinSpaced(2*N + 1, 1, 2*N + 1),
        output->ALF1n, output->ALF2n, output->ALF3n, output->ALF4n, output->BET1n, output->BET2n;

    output->all_F = Array<Real, Dynamic, 7>(Nb_theta, 7);
    output->all_F << theta_deg, output->F11, output->F22, output->F33, output->F44,
        output->F12, output->F34;

    return output;
  }
}

#endif
