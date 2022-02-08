/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Based on the SMARTIES MATLAB package by W.R.C. Somerville, B. Auguié, E.C. Le Ru
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.


This code contains all 'rvh' SMARTIES functions that are used in Raman scattering calculations,
i.e. T-matrix related functions specific to particles with mirror-reflection symmetry.
Henceforth, such symmetry shall be called rvh symmetry.
*/

#ifndef RVH_HPP
#define RVH_HPP

#include "core.hpp"
#include "sph.hpp"
#include "vsh.hpp"
#include "misc.hpp"

using namespace Eigen;
using namespace std;

namespace Smarties {

  /* Specialisation of stMat for a describing a T, R matrix pair */
  template <class Real>
  struct stTR : stMat<Real> {
    using stMat<Real>::stMat;
    // Two non-zero submatrices of a T-matrix
    inline st4M<Real>& st_4M_T_eo() { return this->st_4M_list[0]; }
    inline st4M<Real>& st_4M_T_oe() { return this->st_4M_list[1]; }
    // Two non-zero submatrices of an R-matrix
    inline st4M<Real>& st_4M_R_eo() { return this->st_4M_list[2]; }
    inline st4M<Real>& st_4M_R_oe() { return this->st_4M_list[3]; }
  };

  /* Struct containing vectors of expansion coefficients of the multipole expansions of a set of electric fields */
  template <class Real>
  struct stAbcdnm {
    // Scattering field coefficients
    ArrayXc<Real> p_nm;
    ArrayXc<Real> q_nm;
    // Incident field coefficients
    ArrayXc<Real> a_nm;
    ArrayXc<Real> b_nm;
    // Internal field coefficients
    ArrayXc<Real> c_nm;
    ArrayXc<Real> d_nm;
  };

  /* Struct containing orientation-averaged cross-sections for various wavelengths*/
  template <class Real>
  struct stCrossSection {
    ArrayXr<Real> C_ext; // Wavelength-dependent extinction coefficients
    ArrayXr<Real> C_sca; // Wavelength-dependent scattering coefficients
    ArrayXr<Real> C_abs; // Wavelength-dependent absorption coefficients
  };


  /*
  Calculate T (and possibly R) matrices from P, Q matrices, for scatterers with a plane of symmetry.
  This makes use of block inversion. See sec. 4.5 and eq. 70 in JQSRT2013 for details.
  Inputs:
    st_PQ_list - std::vector of size M, each element of which contains a unique pointer to an stPQ struct
                 (one entry for each m in abs_m_vec), and which makes use of the rvh symmetry
    get_R - if true, R is computed and returned
  Output:
    A std::vector of size M of unique pointers to stTR structs.
  Dependencies:
    InvertLUcol
  */
  template <class Real>
  vector<unique_ptr<stTR<Real>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<Real>>>& st_PQ_list, bool get_R = false) {
    // Note that all matrix inversions uses LU decomposition with partial pivoting of the columns
    int num_entries = st_PQ_list.size();
    vector<unique_ptr<stTR<Real>>> output(num_entries);
    for (int i = 0; i < num_entries; i++) {
      output[i] = make_unique<stTR<Real>>();
      // Get blocks of P and Q matrices
      unique_ptr<stPQ<Real>>& st_PQ = st_PQ_list[i];
      unique_ptr<stTR<Real>>& st_TR = output[i];
      // Get blocks of P_eo and Q_eo matrices
      MatrixXc<Real> P11_ee = st_PQ->st_4M_P_eo().M11;
      MatrixXc<Real> P12_eo = st_PQ->st_4M_P_eo().M12;
      MatrixXc<Real> P21_oe = st_PQ->st_4M_P_eo().M21;
      MatrixXc<Real> P22_oo = st_PQ->st_4M_P_eo().M22;
      MatrixXc<Real> Q11_ee = st_PQ->st_4M_Q_eo().M11;
      MatrixXc<Real> Q12_eo = st_PQ->st_4M_Q_eo().M12;
      MatrixXc<Real> Q21_oe = st_PQ->st_4M_Q_eo().M21;
      MatrixXc<Real> Q22_oo = st_PQ->st_4M_Q_eo().M22;

      // Get blocks of P_oe and Q_oe matrices
      MatrixXc<Real> P11_oo = st_PQ->st_4M_P_oe().M11;
      MatrixXc<Real> P12_oe = st_PQ->st_4M_P_oe().M12;
      MatrixXc<Real> P21_eo = st_PQ->st_4M_P_oe().M21;
      MatrixXc<Real> P22_ee = st_PQ->st_4M_P_oe().M22;
      MatrixXc<Real> Q11_oo = st_PQ->st_4M_Q_oe().M11;
      MatrixXc<Real> Q12_oe = st_PQ->st_4M_Q_oe().M12;
      MatrixXc<Real> Q21_eo = st_PQ->st_4M_Q_oe().M21;
      MatrixXc<Real> Q22_ee = st_PQ->st_4M_Q_oe().M22;

      int num_even = Q11_ee.rows();
      int num_odd = Q11_oo.rows();
      ArrayXi ind1_eo = st_PQ->st_4M_Q_eo().ind1;
      ArrayXi ind2_eo = st_PQ->st_4M_Q_eo().ind2;
      ArrayXi ind1_oe = st_PQ->st_4M_Q_oe().ind1;
      ArrayXi ind2_oe = st_PQ->st_4M_Q_oe().ind2;

      int m = st_PQ->st_4M_P_eo().m;

      // m == 0 case is done separately, since M12 and M21 are empty.
      if (!m) {
        st_TR->st_4M_R_oe().M11 = InvertLUcol(Q11_oo);
        st_TR->st_4M_T_oe().M11 = -P11_oo * st_TR->st_4M_R_oe().M11.matrix();
        st_TR->st_4M_R_oe().M22 = InvertLUcol(Q22_ee);
        st_TR->st_4M_T_oe().M22 = -P22_ee * st_TR->st_4M_R_oe().M22.matrix();
        st_TR->st_4M_R_eo().M11 = InvertLUcol(Q11_ee);
        st_TR->st_4M_T_eo().M11 = -P11_ee * st_TR->st_4M_R_eo().M11.matrix();
        st_TR->st_4M_R_eo().M22 = InvertLUcol(Q22_oo);
        st_TR->st_4M_T_eo().M22 = -P22_oo * st_TR->st_4M_R_eo().M22.matrix();

        st_TR->st_4M_T_eo().M12 = ArrayXXc<Real>::Zero(num_even, num_odd);
        st_TR->st_4M_T_eo().M21 = ArrayXXc<Real>::Zero(num_odd, num_even);
        st_TR->st_4M_T_oe().M12 = ArrayXXc<Real>::Zero(num_odd, num_even);
        st_TR->st_4M_T_oe().M21 = ArrayXXc<Real>::Zero(num_even, num_odd);

        if (get_R) {
          st_TR->st_4M_R_eo().M12 = ArrayXXc<Real>::Zero(num_even, num_odd);
          st_TR->st_4M_R_eo().M21 = ArrayXXc<Real>::Zero(num_odd, num_even);
          st_TR->st_4M_R_oe().M12 = ArrayXXc<Real>::Zero(num_odd, num_even);
          st_TR->st_4M_R_oe().M21 = ArrayXXc<Real>::Zero(num_even, num_odd);
        }
      } else { // Do inversion using JQSRT2013 eq. 70

        // For M_eo matrices
        MatrixXc<Real> Q11_inv = InvertLUcol(Q11_ee);

        MatrixXc<Real> G1 = P11_ee * Q11_inv;
        MatrixXc<Real> G3 = P21_oe * Q11_inv;
        MatrixXc<Real> G5 = Q21_oe * Q11_inv;
        MatrixXc<Real> F2_m1 = Q22_oo - G5 * Q12_eo;
        MatrixXc<Real> F2 = InvertLUcol(F2_m1);

        MatrixXc<Real> G2 = P22_oo * F2;
        MatrixXc<Real> G4 = P12_eo * F2;
        MatrixXc<Real> G6 = Q12_eo * F2;

        st_TR->st_4M_T_eo().M12 = G1 * G6 - G4;
        st_TR->st_4M_T_eo().M22 = G3 * G6 - G2;
        st_TR->st_4M_T_eo().M11 = -G1 - st_TR->st_4M_T_eo().M12.matrix() * G5;
        st_TR->st_4M_T_eo().M21 = -G3 - st_TR->st_4M_T_eo().M22.matrix() * G5;

        if (get_R) {
          st_TR->st_4M_R_eo().M12 = -Q11_inv * G6;
          st_TR->st_4M_R_eo().M22 = F2;
          st_TR->st_4M_R_eo().M11 = Q11_inv - st_TR->st_4M_R_eo().M12.matrix() * G5;
          st_TR->st_4M_R_eo().M21 = -st_TR->st_4M_R_eo().M22.matrix() * G5;
        }

        // For M_oe matrices
        Q11_inv = InvertLUcol(Q11_oo);

        G1 = P11_oo * Q11_inv;
        G3 = P21_eo * Q11_inv;
        G5 = Q21_eo * Q11_inv;
        F2_m1 = Q22_ee - G5 * Q12_oe;
        F2 = InvertLUcol(F2_m1);

        G2 = P22_ee * F2;
        G4 = P12_oe * F2;
        G6 = Q12_oe * F2;

        st_TR->st_4M_T_oe().M12 = G1 * G6 - G4;
        st_TR->st_4M_T_oe().M22 = G3 * G6 - G2;
        st_TR->st_4M_T_oe().M11 = -G1 - st_TR->st_4M_T_oe().M12.matrix() * G5;
        st_TR->st_4M_T_oe().M21 = -G3 - st_TR->st_4M_T_oe().M22.matrix() * G5;

        if (get_R) {
          st_TR->st_4M_R_oe().M12 = -Q11_inv * G6;
          st_TR->st_4M_R_oe().M22 = F2;
          st_TR->st_4M_R_oe().M11 = Q11_inv - st_TR->st_4M_R_oe().M12.matrix() * G5;
          st_TR->st_4M_R_oe().M21 = -st_TR->st_4M_R_oe().M22.matrix() * G5;
        }
      }

      st_TR->st_4M_T_eo().m = st_TR->st_4M_T_oe().m = m;
      st_TR->st_4M_T_eo().ind1 = ind1_eo;
      st_TR->st_4M_T_eo().ind2 = ind2_eo;
      st_TR->st_4M_T_oe().ind1 = ind1_oe;
      st_TR->st_4M_T_oe().ind2 = ind2_oe;
      st_TR->mat_list.push_back("st_4M_T");

      if (get_R) {
        st_TR->st_4M_R_eo().m = st_TR->st_4M_R_oe().m = m;
        st_TR->st_4M_R_eo().ind1 = ind1_eo;
        st_TR->st_4M_R_eo().ind2 = ind2_eo;
        st_TR->st_4M_R_oe().ind1 = ind1_oe;
        st_TR->st_4M_R_oe().ind2 = ind2_oe;
        st_TR->mat_list.push_back("st_4M_R");
      }
    }
    return output;
  }

  /*
  Truncate T and R matrices in st_mat_list such that the maximum number of multipoles N is now N_max.
  This also removes values of m that are larger than N_max.
  This function works on st_mat_list which makes use of rvh symmetry.
  Note that this function has only been implemented for T and R matrices, whereas the original MATLAB
  code also allows truncation of P and Q matrices. Although this could supposedly be implemented in
  by overloading the function or adding another template argument.
  Inputs:
    st_mat_list - a vector of unique pointers to stTR structs, with one entry for each m of abs_m_vec.
    N_max - The maximum value of N (and m) that is desired in the output.
  Output:
    returns a deep copy of the same st_mat_list vector, but where all M11, M12, M21 and M22 members
    have been truncated to N = N_max and values of m > N_max are removed.
  Dependencies:
    LogicalIndices
  */
  template <class Real>
  vector<unique_ptr<stTR<Real>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<Real>>>& st_mat_list, int N_max) {
    vector<unique_ptr<stTR<Real>>> output;
    int num_entries = st_mat_list.size();
    for (int i = 0 ; i < num_entries; i++) {
      int num_st_4M = st_mat_list[i]->mat_list.size() * 2; // Count how many matrices in the struct are defined
      int m;
      if (num_st_4M > 0 && (m = st_mat_list[i]->st_4M_list[0].m) <= N_max) {
        auto output_st_TR = make_unique<stTR<Real>>();
        output_st_TR->mat_list = st_mat_list[i]->mat_list;
        for (int j = 0; j < num_st_4M; j++) {
          int new_size = N_max - max(1, m) + 1;
          ArrayXb ind1_valid = st_mat_list[i]->st_4M_list[j].ind1 <= new_size;
          ArrayXb ind2_valid = st_mat_list[i]->st_4M_list[j].ind2 <= new_size;

          ArrayXi new_ind1 = LogicalIndices(ind1_valid);
          ArrayXi new_ind2 = LogicalIndices(ind2_valid);

          output_st_TR->st_4M_list[j].ind1 = st_mat_list[i]->st_4M_list[j].ind1(new_ind1);
          output_st_TR->st_4M_list[j].ind2 = st_mat_list[i]->st_4M_list[j].ind2(new_ind2);
          output_st_TR->st_4M_list[j].m = m;

          ArrayXXc<Real> current_matrix = st_mat_list[i]->st_4M_list[j].M11;
          output_st_TR->st_4M_list[j].M11 = current_matrix(new_ind1, new_ind1);
          current_matrix = st_mat_list[i]->st_4M_list[j].M12;
          output_st_TR->st_4M_list[j].M12 = current_matrix(new_ind1, new_ind2);
          current_matrix = st_mat_list[i]->st_4M_list[j].M21;
          output_st_TR->st_4M_list[j].M21 = current_matrix(new_ind2, new_ind1);
          current_matrix = st_mat_list[i]->st_4M_list[j].M22;
          output_st_TR->st_4M_list[j].M22 = current_matrix(new_ind2, new_ind2);
        }
        output.push_back(move(output_st_TR));
      }
    }
    return output;
  }

  /*
  Symmetrises the matrices given in mat_list, and leaves other matrices present in st_mat_list unchanged.
  Uses the upper triangular matrices to deduce the lower triangular parts from the symmetry relations
  of the T-matrix obtained from eqs. 5.34 and 5.37 of [Mishchenko 2002], i.e.
  T11 = T11.transpose(), T22 = T22.transpose(), T12 = -T21.transpose(), T21 = -T12.transpose()
  This function assumes rvh symmetry. While the original MATLAB code could theoretically work with P and Q matrices,
  the implementation here only works with T and R matrices, although it's intended only for T and R matrices anyway.
  Inputs:
    st_mat_list - a std::vector containing unique pointers to stTR structs containing
                  T (and possibly R) matrices to symmetrise
    mat_list - a std::vector of strings specifying which structs are to be symmetrised
  Outputs:
    A deep copy of st_mat_list, but with all the desired matrices symmetrised.
  */
  template <class Real>
  vector<unique_ptr<stTR<Real>>> rvhGetSymmetricMat(const vector<unique_ptr<stTR<Real>>>& st_mat_list,
      vector<string> mat_list = {"st_4M_T"}) {
    enum parity {EO, OE, END};
    int num_entries = st_mat_list.size();
    vector<unique_ptr<stTR<Real>>> output(num_entries);

    for (int i = 0; i < num_entries; i++) {
      output[i] = make_unique<stTR<Real>>(*(st_mat_list[i])); // Make a deep copy at the ith entry of st_mat_list
      int num_st_4M = mat_list.size();
      for (int j = 0; j < num_st_4M; j++) {
        for (int k = EO; k != END; k++) { // Symmetrise eo matrix, then oe matrix
          st4M<Real>* st_4M;
          if (mat_list[j] == "st_4M_T") {
            if (k == EO)
              st_4M = &(output[i]->st_4M_T_eo());
            else if (k == OE)
              st_4M = &(output[i]->st_4M_T_oe());
          } else if (mat_list[j] == "st_4M_R") {
            if (k == EO)
              st_4M = &(output[i]->st_4M_R_eo());
            else if (k == OE)
              st_4M = &(output[i]->st_4M_R_oe());
          } else
            continue;
          ArrayXi ind1 = st_4M->ind1;
          ArrayXi ind2 = st_4M->ind2;

          if (ind1.size() && ind2.size()) { // Check that neither two are empty
            int offset12 = (ind2(0) - ind1(0) + 1)/2;
            // int offset21 = 1 - offset12;

            // Upper triangular without diagonal
            MatrixXc<Real> upper = st_4M->M11.matrix().template triangularView<StrictlyUpper>();
            // Diagonal only
            MatrixXc<Real> diag = static_cast<VectorXc<Real>>(st_4M->M11.matrix().diagonal()).asDiagonal();
            st_4M->M11 = (upper + diag + upper.transpose()).array(); // Symmtrised matrix

            upper = st_4M->M22.matrix().template triangularView<StrictlyUpper>();
            diag = static_cast<VectorXc<Real>>(st_4M->M22.matrix().diagonal()).asDiagonal();
            st_4M->M22 = (upper + diag + upper.transpose()).array();

            // In either case, the transpose of upper2 complements upper1
            MatrixXc<Real> upper1, upper2;
            if (offset12) {
              upper1 = st_4M->M12.matrix().template triangularView<StrictlyUpper>();
              upper2 = st_4M->M21.matrix().template triangularView<Upper>();
            } else {
              upper1 = st_4M->M12.matrix().template triangularView<Upper>();
              upper2 = st_4M->M21.matrix().template triangularView<StrictlyUpper>();
            }
            st_4M->M12 = (upper1 - upper2.transpose()).array(); // Symmetrised matrix
            st_4M->M21 = -st_4M->M12.transpose(); // Symmetrised matrix
          }
        }
      }
    }

    return output;
  }

  /*
  Calculates the field expansion coefficients from the T/R-matrices for a given incident plane wave
  (for one wavelength only). If R is not defined in st_TR_list, then the internal fields are not computed.
  If stIncEabnm is given, the expansion coefficients for the incident wave are not recalculated.
  This method is valid for scatterers with rvh symmetry.

  If st_TR_list contains m values that are not required, then they are ignored.
  Inputs:
    N_max - the maximum multipole order to use in the expansions
            (should be less than or equal to the T-matrix size)
    st_TR_list - std::vector of size M containing unique pointers to stTR structs for each m.
                 This should contain all m present in st_inc_par.abs_m_vec, i.e. all m where |m| <= N,
                 and should make use of the rvh symmetry.
    st_inc_par - A struct containing information about the incident field, as from vshMakeIncidentParams.
    st_inc_E_abmm - A struct containing the expansion coefficients a_nm and b_nm of the incident wave.
                    These are recalculated if not provided.
  Output:
    returns a stAbcdnm struct containing all the calculated expansion coefficients.
    The case where n = 0 and m = 0 is included for padding, so every vector in the struct has size
    P = (N_max + 1)^2. The padding exists to simplify the indexing.
    c_nm and d_nm are only defined if R is defined in st_TR_list.
  Dependencies:
    Seq2Array, LogicalSlice, vshGetIncidentCoefficients
  */
  template <class Real>
  unique_ptr<stAbcdnm<Real>> rvhGetFieldCoefficients(int N_max,
      const vector<unique_ptr<stTR<Real>>>& st_TR_list, const unique_ptr<stIncPar<Real>>& st_inc_par,
      unique_ptr<stIncEabnm<Real>> st_inc_E_abnm = unique_ptr<stIncEabnm<Real>>()) {

    // Coeff of incident wave
    if (!st_inc_E_abnm)
      st_inc_E_abnm = vshGetIncidentCoeffs(N_max, st_inc_par);

    bool get_R = false;
    if (find(st_TR_list[0]->mat_list.begin(), st_TR_list[0]->mat_list.end(), "st_4M_R")
        != st_TR_list[0]->mat_list.end())
      get_R = true;

    // Truncate to relevant m's
    ArrayXb abs_m_vec_valid = st_inc_par->abs_m_vec <= N_max;
    ArrayXi abs_m_vec = LogicalSlice(st_inc_par->abs_m_vec, abs_m_vec_valid);

    int num_M = st_TR_list.size(); // number of m-values in T, not all of which may be needed
    int M = abs_m_vec.size();
    int P_max = (N_max + 1)*(N_max + 1); // number of (n,m) coupled indices (with (0,0) included as padding)

    ArrayXc<Real> a_nm = st_inc_E_abnm->a_nm; // [P X 1]
    ArrayXc<Real> b_nm = st_inc_E_abnm->b_nm; // [P X 1]
    ArrayXc<Real> p_nm = ArrayXc<Real>::Zero(P_max); // [P X 1]
    ArrayXc<Real> q_nm = ArrayXc<Real>::Zero(P_max); // [P X 1]
    ArrayXc<Real> c_nm = ArrayXc<Real>::Zero(get_R ? P_max : 1); // if get_R, [P X 1]
    ArrayXc<Real> d_nm = ArrayXc<Real>::Zero(get_R ? P_max : 1); // if get_R, [P X 1]
    ArrayXi m_ind_for_T(M);

    // Find the corresponding m-indices in st_TR_list (in case they are not in standard order)
    for (int m_ind11 = 0; m_ind11 < M; m_ind11++) {
      for (int m_ind = 0; m_ind < num_M; m_ind++) {
        if (st_TR_list[m_ind]->st_4M_T_eo().m == abs_m_vec(m_ind11)) {
          m_ind_for_T(m_ind11) = m_ind;
          break;
        }
      }
    }

    for (int m_ind11 = 0; m_ind11 < M; m_ind11++) {
      int m_ind = m_ind_for_T(m_ind11);
      int m = abs_m_vec(m_ind11);
      int N_min = max(m, 1);

      ArrayXi n_vec_e = Seq2Array(N_min + (N_min % 2), N_max, 2); // indices for even n [N_e X 1]
      ArrayXi n_vec_o = Seq2Array(N_min + 1 - (N_min % 2), N_max, 2); // indices for odd n [N_o X 1]

      ArrayXi p_vec_e = n_vec_e * (n_vec_e + 1) + m; // [N_e X 1]
      ArrayXi p_vec_o = n_vec_o * (n_vec_o + 1) + m; // [N_o X 1]
      ArrayXi n_vec_T_e = ArrayXi::LinSpaced(n_vec_e.size(), 0, n_vec_e.size() - 1);
      ArrayXi n_vec_T_o = ArrayXi::LinSpaced(n_vec_o.size(), 0, n_vec_o.size() - 1);

      VectorXc<Real> p_nm_e = st_TR_list[m_ind]->st_4M_T_eo().M11(n_vec_T_e, n_vec_T_e).matrix() *
          a_nm(p_vec_e).matrix() + st_TR_list[m_ind]->st_4M_T_eo().M12(n_vec_T_e, n_vec_T_o).matrix() *
          b_nm(p_vec_o).matrix();
      VectorXc<Real> q_nm_o = st_TR_list[m_ind]->st_4M_T_eo().M21(n_vec_T_o, n_vec_T_e).matrix() *
          a_nm(p_vec_e).matrix() + st_TR_list[m_ind]->st_4M_T_eo().M22(n_vec_T_o, n_vec_T_o).matrix() *
          b_nm(p_vec_o).matrix();
      VectorXc<Real> p_nm_o = st_TR_list[m_ind]->st_4M_T_oe().M11(n_vec_T_o, n_vec_T_o).matrix() *
          a_nm(p_vec_o).matrix() + st_TR_list[m_ind]->st_4M_T_oe().M12(n_vec_T_o, n_vec_T_e).matrix() *
          b_nm(p_vec_e).matrix();
      VectorXc<Real> q_nm_e = st_TR_list[m_ind]->st_4M_T_oe().M21(n_vec_T_e, n_vec_T_o).matrix() *
          a_nm(p_vec_o).matrix() + st_TR_list[m_ind]->st_4M_T_oe().M22(n_vec_T_e, n_vec_T_e).matrix() *
          b_nm(p_vec_e).matrix();

      for (int i = 0; i < n_vec_e.size(); i++) {
        p_nm(p_vec_e(i)) = p_nm_e(i);
        q_nm(p_vec_e(i)) = q_nm_e(i);
      }
      for (int i = 0; i < n_vec_o.size(); i++) {
        p_nm(p_vec_o(i)) = p_nm_o(i);
        q_nm(p_vec_o(i)) = q_nm_o(i);
      }

      if (get_R) {
        VectorXc<Real> c_nm_e = st_TR_list[m_ind]->st_4M_R_eo().M11(n_vec_T_e, n_vec_T_e).matrix() *
            a_nm(p_vec_e).matrix() + st_TR_list[m_ind]->st_4M_R_eo().M12(n_vec_T_e, n_vec_T_o).matrix() *
            b_nm(p_vec_o).matrix();
        VectorXc<Real> d_nm_o = st_TR_list[m_ind]->st_4M_R_eo().M21(n_vec_T_o, n_vec_T_e).matrix() *
            a_nm(p_vec_e).matrix() + st_TR_list[m_ind]->st_4M_R_eo().M22(n_vec_T_o, n_vec_T_o).matrix() *
            b_nm(p_vec_o).matrix();
        VectorXc<Real> c_nm_o = st_TR_list[m_ind]->st_4M_R_oe().M11(n_vec_T_o, n_vec_T_o).matrix() *
            a_nm(p_vec_o).matrix() + st_TR_list[m_ind]->st_4M_R_oe().M12(n_vec_T_o, n_vec_T_e).matrix() *
            b_nm(p_vec_e).matrix();
        VectorXc<Real> d_nm_e = st_TR_list[m_ind]->st_4M_R_oe().M21(n_vec_T_e, n_vec_T_o).matrix() *
            a_nm(p_vec_o).matrix() + st_TR_list[m_ind]->st_4M_R_oe().M22(n_vec_T_e, n_vec_T_e).matrix() *
            b_nm(p_vec_e).matrix();

        for (int i = 0; i < n_vec_e.size(); i++) {
          c_nm(p_vec_e(i)) = c_nm_e(i);
          d_nm(p_vec_e(i)) = d_nm_e(i);
        }
        for (int i = 0; i < n_vec_o.size(); i++) {
          c_nm(p_vec_o(i)) = c_nm_o(i);
          d_nm(p_vec_o(i)) = d_nm_o(i);
        }
      }

      // Calculate for m < 0
      if (m) {
        ArrayXi p_vec_n_e = n_vec_e * (n_vec_e + 1) - m; // [N_e X 1]
        ArrayXi p_vec_n_o = n_vec_o * (n_vec_o + 1) - m; // [N_e X 1]

        // For negative m, using eq. 5.37 of [Mishchenko 2002]
        VectorXc<Real> p_nm_n_e = st_TR_list[m_ind]->st_4M_T_eo().M11(n_vec_T_e, n_vec_T_e).matrix() *
            a_nm(p_vec_n_e).matrix() - st_TR_list[m_ind]->st_4M_T_eo().M12(n_vec_T_e, n_vec_T_o).matrix() *
            b_nm(p_vec_n_o).matrix();
        VectorXc<Real> q_nm_n_o = -st_TR_list[m_ind]->st_4M_T_eo().M21(n_vec_T_o, n_vec_T_e).matrix() *
            a_nm(p_vec_n_e).matrix() + st_TR_list[m_ind]->st_4M_T_eo().M22(n_vec_T_o, n_vec_T_o).matrix() *
            b_nm(p_vec_n_o).matrix();
        VectorXc<Real> p_nm_n_o = st_TR_list[m_ind]->st_4M_T_oe().M11(n_vec_T_o, n_vec_T_o).matrix() *
            a_nm(p_vec_n_o).matrix() - st_TR_list[m_ind]->st_4M_T_oe().M12(n_vec_T_o, n_vec_T_e).matrix() *
            b_nm(p_vec_n_e).matrix();
        VectorXc<Real> q_nm_n_e = -st_TR_list[m_ind]->st_4M_T_oe().M21(n_vec_T_e, n_vec_T_o).matrix() *
            a_nm(p_vec_n_o).matrix() + st_TR_list[m_ind]->st_4M_T_oe().M22(n_vec_T_e, n_vec_T_e).matrix() *
            b_nm(p_vec_n_e).matrix();

        for (int i = 0; i < n_vec_e.size(); i++) {
          p_nm(p_vec_n_e(i)) = p_nm_n_e(i);
          q_nm(p_vec_n_e(i)) = q_nm_n_e(i);
        }
        for (int i = 0; i < n_vec_o.size(); i++) {
          p_nm(p_vec_n_o(i)) = p_nm_n_o(i);
          q_nm(p_vec_n_o(i)) = q_nm_n_o(i);
        }

        if (get_R) {
          VectorXc<Real> c_nm_n_e = st_TR_list[m_ind]->st_4M_R_eo().M11(n_vec_T_e, n_vec_T_e).matrix() *
              a_nm(p_vec_n_e).matrix() - st_TR_list[m_ind]->st_4M_R_eo().M12(n_vec_T_e, n_vec_T_o).matrix() *
              b_nm(p_vec_n_o).matrix();
          VectorXc<Real> d_nm_n_o = -st_TR_list[m_ind]->st_4M_R_eo().M21(n_vec_T_o, n_vec_T_e).matrix() *
              a_nm(p_vec_n_e).matrix() + st_TR_list[m_ind]->st_4M_R_eo().M22(n_vec_T_o, n_vec_T_o).matrix() *
              b_nm(p_vec_n_o).matrix();
          VectorXc<Real> c_nm_n_o = st_TR_list[m_ind]->st_4M_R_oe().M11(n_vec_T_o, n_vec_T_o).matrix() *
              a_nm(p_vec_n_o).matrix() - st_TR_list[m_ind]->st_4M_R_oe().M12(n_vec_T_o, n_vec_T_e).matrix() *
              b_nm(p_vec_n_e).matrix();
          VectorXc<Real> d_nm_n_e = -st_TR_list[m_ind]->st_4M_R_oe().M21(n_vec_T_e, n_vec_T_o).matrix() *
              a_nm(p_vec_n_o).matrix() + st_TR_list[m_ind]->st_4M_R_oe().M22(n_vec_T_e, n_vec_T_e).matrix() *
              b_nm(p_vec_n_e).matrix();

          for (int i = 0; i < n_vec_e.size(); i++) {
            c_nm(p_vec_n_e(i)) = c_nm_n_e(i);
            d_nm(p_vec_n_e(i)) = d_nm_n_e(i);
          }
          for (int i = 0; i < n_vec_o.size(); i++) {
            c_nm(p_vec_n_o(i)) = c_nm_n_o(i);
            d_nm(p_vec_n_o(i)) = d_nm_n_o(i);
          }
        }
      }
    }

    auto output = make_unique<stAbcdnm<Real>>();
    output->p_nm = p_nm; // [P X 1]
    output->q_nm = q_nm; // [P X 1]
    output->a_nm = a_nm; // [P X 1]
    output->b_nm = b_nm; // [P X 1]
    output->c_nm = c_nm; // if get_R, [P X 1]
    output->d_nm = d_nm; // if get_R, [P X 1]
    return output;
  }

  /*
  Calculate absorption, scattering and extinction cross sections, for the case of orientationally-averaged incidence.
  This is for matrices which make use of rvh symmetry.
  Inputs:
    k1 - wavevector [L X 1]
    st_TR_list - std::vector containing unique pointers to stTR structs (in rvh-block form).
                 Should be [L X M] and should contain all m where 0 <= m <= N
  Output:
    Returns a stCrossSection struct, where each member is of size [L X 1] and
    contains information for each wavelength.
  Dependencies:
    mp_pi
  */
  template <class Real>
  unique_ptr<stCrossSection<Real>> rvhGetAverageCrossSections(
      const ArrayXr<Real>& k1, const vector<vector<unique_ptr<stTR<Real>>>>& st_TR_list) {
    int L = st_TR_list.size();
    int M = st_TR_list[0].size();

    ArrayXc<Real> ext_sum = ArrayXc<Real>::Zero(L);
    ArrayXc<Real> sca_sum = ArrayXc<Real>::Zero(L);

    for (int l = 0; l < L; l++) { // Loop on lambda
      if (st_TR_list[l][0]->st_4M_T_eo().ind1.size() + st_TR_list[l][0]->st_4M_T_eo().ind2.size() + 1 != M)
        cout << "Warning in rvhGetAverageCrossSections: CstTRa does not seem to contain T-matrices for all m" << endl;
      for (int m = 0; m < M; m++) {
        if (st_TR_list[l][m]->st_4M_T_eo().m != m || st_TR_list[l][m]->st_4M_T_oe().m != m)
          cout << "Warning in rvhGetAverageCrossSections: CstTRa does not seem to contain T-matrices for all m" << endl;
        Real m_factor = m ? 1 : 0.5;

        // From eq. 5.107 of Mishchenko 2002
        ext_sum(l) += m_factor * (
            st_TR_list[l][m]->st_4M_T_eo().M11.matrix().diagonal().sum() +
            st_TR_list[l][m]->st_4M_T_oe().M11.matrix().diagonal().sum() +
            st_TR_list[l][m]->st_4M_T_eo().M22.matrix().diagonal().sum() +
            st_TR_list[l][m]->st_4M_T_oe().M22.matrix().diagonal().sum());

        // From eq. 5.141 of Mishchenko 2002
        sca_sum(l) += m_factor * (
            st_TR_list[l][m]->st_4M_T_eo().M11.abs().pow(2).sum() +
            st_TR_list[l][m]->st_4M_T_eo().M12.abs().pow(2).sum() +
            st_TR_list[l][m]->st_4M_T_eo().M21.abs().pow(2).sum() +
            st_TR_list[l][m]->st_4M_T_eo().M22.abs().pow(2).sum() +
            st_TR_list[l][m]->st_4M_T_oe().M11.abs().pow(2).sum() +
            st_TR_list[l][m]->st_4M_T_oe().M12.abs().pow(2).sum() +
            st_TR_list[l][m]->st_4M_T_oe().M21.abs().pow(2).sum() +
            st_TR_list[l][m]->st_4M_T_oe().M22.abs().pow(2).sum());
      }
    }

    auto output = make_unique<stCrossSection<Real>>();
    output->C_ext = -4*mp_pi<Real>()/k1.pow(2) * ext_sum.real(); // [L X 1]
    output->C_sca = 4*mp_pi<Real>()/k1.pow(2) * sca_sum.real(); // [L X 1]
    output->C_abs = -4*mp_pi<Real>()/k1.pow(2) * (ext_sum.real() + sca_sum.real()); // [L X 1]

    return output;
  }
}

#endif
