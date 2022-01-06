#include "rvh.h"
#include "math.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  vector<unique_ptr<stTR<Real>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<Real>>>& st_PQ_list, bool get_R) {
    size_t num_entries = st_PQ_list.size();
    vector<unique_ptr<stTR<Real>>> output(num_entries);
    for (size_t i = 0; i < num_entries; i++) {
      output[i] = make_unique<stTR<Real>>();
      unique_ptr<stPQ<Real>>& st_PQ = st_PQ_list[i];
      unique_ptr<stTR<Real>>& st_TR = output[i];
      MatrixXc<Real> P11_ee = st_PQ->st_4M_P_eo().M11;
      MatrixXc<Real> P12_eo = st_PQ->st_4M_P_eo().M12;
      MatrixXc<Real> P21_oe = st_PQ->st_4M_P_eo().M21;
      MatrixXc<Real> P22_oo = st_PQ->st_4M_P_eo().M22;
      MatrixXc<Real> Q11_ee = st_PQ->st_4M_Q_eo().M11;
      MatrixXc<Real> Q12_eo = st_PQ->st_4M_Q_eo().M12;
      MatrixXc<Real> Q21_oe = st_PQ->st_4M_Q_eo().M21;
      MatrixXc<Real> Q22_oo = st_PQ->st_4M_Q_eo().M22;

      MatrixXc<Real> P11_oo = st_PQ->st_4M_P_oe().M11;
      MatrixXc<Real> P12_oe = st_PQ->st_4M_P_oe().M12;
      MatrixXc<Real> P21_eo = st_PQ->st_4M_P_oe().M21;
      MatrixXc<Real> P22_ee = st_PQ->st_4M_P_oe().M22;
      MatrixXc<Real> Q11_oo = st_PQ->st_4M_Q_oe().M11;
      MatrixXc<Real> Q12_oe = st_PQ->st_4M_Q_oe().M12;
      MatrixXc<Real> Q21_eo = st_PQ->st_4M_Q_oe().M21;
      MatrixXc<Real> Q22_ee = st_PQ->st_4M_Q_oe().M22;

      size_t num_even = Q11_ee.rows();
      size_t num_odd = Q11_oo.rows();
      ArrayXi ind1_eo = st_PQ->st_4M_Q_eo().ind1;
      ArrayXi ind2_eo = st_PQ->st_4M_Q_eo().ind2;
      ArrayXi ind1_oe = st_PQ->st_4M_Q_oe().ind1;
      ArrayXi ind2_oe = st_PQ->st_4M_Q_oe().ind2;

      int m = st_PQ->st_4M_P_eo().m;

      if (!m) {
        st_TR->st_4M_R_oe().M11 = invertLUcol(Q11_oo);
        st_TR->st_4M_T_oe().M11 = -P11_oo * st_TR->st_4M_R_oe().M11.matrix();
        st_TR->st_4M_R_oe().M22 = invertLUcol(Q22_ee);
        st_TR->st_4M_T_oe().M22 = -P22_ee * st_TR->st_4M_R_oe().M22.matrix();
        st_TR->st_4M_R_eo().M11 = invertLUcol(Q11_ee);
        st_TR->st_4M_T_eo().M11 = -P11_ee * st_TR->st_4M_R_eo().M11.matrix();
        st_TR->st_4M_R_eo().M22 = invertLUcol(Q22_oo);
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
      } else {
        MatrixXc<Real> Q11_inv = invertLUcol(Q11_ee);

        MatrixXc<Real> G1 = P11_ee * Q11_inv;
        MatrixXc<Real> G3 = P21_oe * Q11_inv;
        MatrixXc<Real> G5 = Q21_oe * Q11_inv;
        MatrixXc<Real> F2_m1 = Q22_oo - G5 * Q12_eo;
        MatrixXc<Real> F2 = invertLUcol(F2_m1);

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

        Q11_inv = invertLUcol(Q11_oo);

        G1 = P11_oo * Q11_inv;
        G3 = P21_eo * Q11_inv;
        G5 = Q21_eo * Q11_inv;
        F2_m1 = Q22_ee - G5 * Q12_oe;
        F2 = invertLUcol(F2_m1);

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

  template <class Real>
  vector<unique_ptr<stTR<Real>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<Real>>>& st_mat_list, int N_max) {
    vector<unique_ptr<stTR<Real>>> output;

    for (size_t i = 0 ; i < st_mat_list.size(); i++) {
      size_t num_st_4M = st_mat_list[i]->mat_list.size() * 2;
      int m;
      if (num_st_4M > 0 && (m = st_mat_list[i]->st_4M_list[0].m) <= N_max) {
        unique_ptr<stTR<Real>> output_st_TR = make_unique<stTR<Real>>();
        output_st_TR->mat_list = st_mat_list[i]->mat_list;
        for (size_t j = 0; j < num_st_4M; j++) {
          int new_size = N_max - max(1, m) + 1;
          ArrayXb ind1_valid = st_mat_list[i]->st_4M_list[j].ind1 <= new_size;
          ArrayXb ind2_valid = st_mat_list[i]->st_4M_list[j].ind2 <= new_size;

          ArrayXi new_ind1 = logicalIndices(ind1_valid);
          ArrayXi new_ind2 = logicalIndices(ind2_valid);

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

  template <class Real>
  vector<unique_ptr<stTR<Real>>> rvhGetSymmetricMat(const vector<unique_ptr<stTR<Real>>>& st_mat_list, vector<string> mat_list) {
    enum parity {EO, OE, END};
    size_t num_entries = st_mat_list.size();
    vector<unique_ptr<stTR<Real>>> output(num_entries);

    for (size_t i = 0; i < num_entries; i++) {
      output[i] = make_unique<stTR<Real>>(*(st_mat_list[i])); //Makes a deep copy of st_mat_list[i]
      for (size_t j = 0; j < 2*mat_list.size(); j++) {
        for (int k = EO; k != END; k++) {
          st4M<Real>* st_4M;
          if (mat_list[j] == "st_4M_T") {
            if (k == EO)
              st_4M = &(output[i]->st_4M_T_eo());
            else if (k == OE)
              st_4M = &(output[i]->st_4M_T_oe());
          }
          else if (mat_list[j] == "st_4M_R") {
            if (k == EO)
              st_4M = &(output[i]->st_4M_T_eo());
            else if (k == OE)
              st_4M = &(output[i]->st_4M_T_oe());
          }
          else
            continue;
          ArrayXi ind1 = st_4M->ind1;
          ArrayXi ind2 = st_4M->ind2;

          if (!ind1.size() && !ind2.size()) {
            int offset12 = (ind1(0) - ind2(0) + 1)/2;
            // int offset21 = 1 - offset12;

            MatrixXc<Real> upper = st_4M->M11.matrix().template triangularView<StrictlyUpper>();
            MatrixXc<Real> diag = static_cast<VectorXc<Real>>(st_4M->M11.matrix().diagonal()).asDiagonal();
            st_4M->M11 = (upper + diag + upper.transpose()).array();

            upper = st_4M->M22.matrix().template triangularView<StrictlyUpper>();
            diag = static_cast<VectorXc<Real>>(st_4M->M22.matrix().diagonal()).asDiagonal();
            st_4M->M22 = (upper + diag + upper.transpose()).array();

            MatrixXc<Real> upper1, upper2;
            if (offset12) {
              upper1 = st_4M->M12.matrix().template triangularView<StrictlyUpper>();
              upper2 = st_4M->M21.matrix().template triangularView<Upper>();
            } else {
              upper1 = st_4M->M12.matrix().template triangularView<Upper>();
              upper2 = st_4M->M21.matrix().template triangularView<StrictlyUpper>();
            }
            st_4M->M12 = (upper1 - upper2.transpose()).array();
            st_4M->M21 = st_4M->M12.transpose();
          }
        }
      }
    }

    return output;
  }

  template <class Real>
  unique_ptr<stAbcdnm<Real>> rvhGetFieldCoefficients(int N_max, const vector<unique_ptr<stTR<Real>>>& st_TR_list,
      const unique_ptr<stIncPar<Real>>& st_inc_par, unique_ptr<stIncEabnm<Real>> st_inc_E_abnm) {
    if (!st_inc_E_abnm)
      st_inc_E_abnm = vshGetIncidentCoeffs(N_max, st_inc_par);

    bool get_R = false;
    if (find(st_TR_list[0]->mat_list.begin(), st_TR_list[0]->mat_list.end(), "st_4M_R") != st_TR_list[0]->mat_list.end())
      get_R = true;

    ArrayXb abs_m_vec_valid = st_inc_par->abs_m_vec <= N_max;
    ArrayXi abs_m_vec = logicalSlice(st_inc_par->abs_m_vec, abs_m_vec_valid);

    size_t num_M = st_TR_list.size();
    size_t M = abs_m_vec.size();
    int P_max = (N_max + 1)*(N_max + 1);

    ArrayXc<Real> a_nm = st_inc_E_abnm->a_nm;
    ArrayXc<Real> b_nm = st_inc_E_abnm->b_nm;
    ArrayXc<Real> p_nm(P_max);
    ArrayXc<Real> q_nm(P_max);
    ArrayXc<Real> c_nm = get_R ? ArrayXc<Real>(P_max) : ArrayXc<Real>::Zero(1);
    ArrayXc<Real> d_nm = get_R ? ArrayXc<Real>(P_max) : ArrayXc<Real>::Zero(1);
    ArrayXi m_for_T(M);

    for (size_t m_ind11 = 0; m_ind11 < M; m_ind11++) {
      for (size_t m_ind = 0; m_ind < num_M; m_ind++) {
        if (st_TR_list[m_ind]->st_4M_T_eo().m == abs_m_vec(m_ind11)) {
          m_for_T(m_ind11) = m_ind;
          break;
        }
      }
    }

    for (size_t m_ind11 = 0; m_ind11 < M; m_ind11++) {
      int m_ind = m_for_T(m_ind11);
      int m = abs_m_vec(m_ind11);

      ArrayXi n_vec_e = seq2Array(m + m % 2, N_max, 2);
      ArrayXi n_vec_o = seq2Array(m + 1 - m % 2, N_max, 2);

      ArrayXi p_vec_e = n_vec_e * (n_vec_e + 1) + m;
      ArrayXi p_vec_o = n_vec_o * (n_vec_o + 1) + m;
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
        VectorXc<Real> c_nm_e = st_TR_list[m_ind]->st_4M_T_eo().M11(n_vec_T_e, n_vec_T_e).matrix() *
            a_nm(p_vec_e).matrix() + st_TR_list[m_ind]->st_4M_T_eo().M12(n_vec_T_e, n_vec_T_o).matrix() *
            b_nm(p_vec_o).matrix();
        VectorXc<Real> d_nm_o = st_TR_list[m_ind]->st_4M_T_eo().M21(n_vec_T_o, n_vec_T_e).matrix() *
            a_nm(p_vec_e).matrix() + st_TR_list[m_ind]->st_4M_T_eo().M22(n_vec_T_o, n_vec_T_o).matrix() *
            b_nm(p_vec_o).matrix();
        VectorXc<Real> c_nm_o = st_TR_list[m_ind]->st_4M_T_oe().M11(n_vec_T_o, n_vec_T_o).matrix() *
            a_nm(p_vec_o).matrix() + st_TR_list[m_ind]->st_4M_T_oe().M12(n_vec_T_o, n_vec_T_e).matrix() *
            b_nm(p_vec_e).matrix();
        VectorXc<Real> d_nm_e = st_TR_list[m_ind]->st_4M_T_oe().M21(n_vec_T_e, n_vec_T_o).matrix() *
            a_nm(p_vec_o).matrix() + st_TR_list[m_ind]->st_4M_T_oe().M22(n_vec_T_e, n_vec_T_e).matrix() *
            b_nm(p_vec_e).matrix();

        for (int i = 0; i < n_vec_e.size(); i++) {
          c_nm(p_vec_e(i)) = p_nm_e(i);
          d_nm(p_vec_e(i)) = q_nm_e(i);
        }
        for (int i = 0; i < n_vec_o.size(); i++) {
          c_nm(p_vec_o(i)) = p_nm_o(i);
          d_nm(p_vec_o(i)) = q_nm_o(i);
        }
      }

      if (m) {
        ArrayXi p_vec_n_e = n_vec_e * (n_vec_e + 1) - m;
        ArrayXi p_vec_n_o = n_vec_o * (n_vec_o + 1) - m;

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
          p_nm(p_vec_n_e(i)) = p_nm_e(i);
          q_nm(p_vec_n_e(i)) = q_nm_e(i);
        }
        for (int i = 0; i < n_vec_o.size(); i++) {
          p_nm(p_vec_n_o(i)) = p_nm_o(i);
          q_nm(p_vec_n_o(i)) = q_nm_o(i);
        }

        if (get_R) {
          VectorXc<Real> c_nm_n_e = st_TR_list[m_ind]->st_4M_T_eo().M11(n_vec_T_e, n_vec_T_e).matrix() *
              a_nm(p_vec_n_e).matrix() - st_TR_list[m_ind]->st_4M_T_eo().M12(n_vec_T_e, n_vec_T_o).matrix() *
              b_nm(p_vec_n_o).matrix();
          VectorXc<Real> d_nm_n_o = -st_TR_list[m_ind]->st_4M_T_eo().M21(n_vec_T_o, n_vec_T_e).matrix() *
              a_nm(p_vec_n_e).matrix() + st_TR_list[m_ind]->st_4M_T_eo().M22(n_vec_T_o, n_vec_T_o).matrix() *
              b_nm(p_vec_n_o).matrix();
          VectorXc<Real> c_nm_n_o = st_TR_list[m_ind]->st_4M_T_oe().M11(n_vec_T_o, n_vec_T_o).matrix() *
              a_nm(p_vec_n_o).matrix() - st_TR_list[m_ind]->st_4M_T_oe().M12(n_vec_T_o, n_vec_T_e).matrix() *
              b_nm(p_vec_n_e).matrix();
          VectorXc<Real> d_nm_n_e = -st_TR_list[m_ind]->st_4M_T_oe().M21(n_vec_T_e, n_vec_T_o).matrix() *
              a_nm(p_vec_n_o).matrix() + st_TR_list[m_ind]->st_4M_T_oe().M22(n_vec_T_e, n_vec_T_e).matrix() *
              b_nm(p_vec_n_e).matrix();

          for (int i = 0; i < n_vec_e.size(); i++) {
            c_nm(p_vec_n_e(i)) = p_nm_e(i);
            d_nm(p_vec_n_e(i)) = q_nm_e(i);
          }
          for (int i = 0; i < n_vec_o.size(); i++) {
            c_nm(p_vec_n_o(i)) = p_nm_o(i);
            d_nm(p_vec_n_o(i)) = q_nm_o(i);
          }
        }
      }
    }

    unique_ptr<stAbcdnm<Real>> output = make_unique<stAbcdnm<Real>>();
    output->p_nm = p_nm;
    output->q_nm = q_nm;
    output->a_nm = a_nm;
    output->b_nm = b_nm;
    output->c_nm = c_nm;
    output->d_nm = d_nm;
    return output;
  }

  template <class Real>
  unique_ptr<stCrossSection<Real>> rvhGetAverageCrossSections(
      const ArrayXr<Real>& k1, const vector<vector<unique_ptr<stTR<Real>>>>& st_TR_list) {
    int L = st_TR_list.size();
    int M = st_TR_list[0].size();

    ArrayXc<Real> ext_sum = ArrayXc<Real>::Zero(L);
    ArrayXc<Real> sca_sum = ArrayXc<Real>::Zero(L);

    for (int l = 0; l < L; l++) {
      if (st_TR_list[l][0]->st_4M_T_eo().ind1.size() + st_TR_list[l][0]->st_4M_T_eo().ind2.size() != M)
        cout << "Warning in rvhGetAverageCrossSections: CstTRa does not seem to contain T-matrices for all m" << endl;
      for (int m = 0; m < M; m++) {
        if (st_TR_list[l][m]->st_4M_T_eo().m != m || st_TR_list[l][m]->st_4M_T_oe().m != m)
          cout << "Warning in rvhGetAverageCrossSections: CstTRa does not seem to contain T-matrices for all m" << endl;
        Real m_factor = m ? 1 : 0.5;

        ext_sum(l) += m_factor * (
            st_TR_list[l][m]->st_4M_T_eo().M11.matrix().diagonal().sum() +
            st_TR_list[l][m]->st_4M_T_oe().M11.matrix().diagonal().sum() +
            st_TR_list[l][m]->st_4M_T_eo().M22.matrix().diagonal().sum() +
            st_TR_list[l][m]->st_4M_T_oe().M22.matrix().diagonal().sum());

        sca_sum(l) += m_factor * (
            st_TR_list[l][m]->st_4M_T_eo().M11.pow(2).sum() +
            st_TR_list[l][m]->st_4M_T_eo().M12.pow(2).sum() +
            st_TR_list[l][m]->st_4M_T_eo().M21.pow(2).sum() +
            st_TR_list[l][m]->st_4M_T_eo().M22.pow(2).sum() +
            st_TR_list[l][m]->st_4M_T_oe().M11.pow(2).sum() +
            st_TR_list[l][m]->st_4M_T_oe().M12.pow(2).sum() +
            st_TR_list[l][m]->st_4M_T_oe().M21.pow(2).sum() +
            st_TR_list[l][m]->st_4M_T_oe().M22.pow(2).sum());
      }
    }

    unique_ptr<stCrossSection<Real>> output = make_unique<stCrossSection<Real>>();
    output->ext = -4*mp_pi<Real>()/k1.pow(2) * ext_sum.real();
    output->sca = 4*mp_pi<Real>()/k1.pow(2) * sca_sum.real();
    output->abs = -4*mp_pi<Real>()/k1.pow(2) * (ext_sum.real() + sca_sum.real());

    return output;
  }

  template vector<unique_ptr<stTR<double>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<double>>>&, bool);
  template vector<unique_ptr<stTR<double>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<double>>>&, int);
  template vector<unique_ptr<stTR<double>>> rvhGetSymmetricMat(const vector<unique_ptr<stTR<double>>>&, vector<string>);
  template unique_ptr<stAbcdnm<double>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<double>>>&,
      const unique_ptr<stIncPar<double>>&, unique_ptr<stIncEabnm<double>>);
  template unique_ptr<stCrossSection<double>> rvhGetAverageCrossSections(
      const ArrayXr<double>&, const vector<vector<unique_ptr<stTR<double>>>>&);
}
