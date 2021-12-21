#include "rvh.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  stTR<Real>* rvhGetTRfromPQ(const vector<stPQ<Real>>* st_PQ_list, bool get_R) {
    size_t num_entries = st_PQ_list->size();
    vector<stTR<Real>>* output = new vector<stTR<Real>>(num_entries, stTR<Real>());
    for (size_t i = 0; i < num_entries; i++) {
      stPQ<Real>* st_PQ = &(st_PQ_list->at(i));
      stPQ<Real>* st_TR = &(output->at(i));
      MatrixXc<Real> P11_ee = st_PQ->st_4M_P_eo->M11;
      MatrixXc<Real> P12_eo = st_PQ->st_4M_P_eo->M12;
      MatrixXc<Real> P21_oe = st_PQ->st_4M_P_eo->M21;
      MatrixXc<Real> P22_oo = st_PQ->st_4M_P_eo->M22;
      MatrixXc<Real> Q11_ee = st_PQ->st_4M_Q_eo->M11;
      MatrixXc<Real> Q12_eo = st_PQ->st_4M_Q_eo->M12;
      MatrixXc<Real> Q21_oe = st_PQ->st_4M_Q_eo->M21;
      MatrixXc<Real> Q22_oo = st_PQ->st_4M_Q_eo->M22;

      MatrixXc<Real> P11_oo = st_PQ->st_4M_P_oe->M11;
      MatrixXc<Real> P12_oe = st_PQ->st_4M_P_oe->M12;
      MatrixXc<Real> P21_eo = st_PQ->st_4M_P_oe->M21;
      MatrixXc<Real> P22_ee = st_PQ->st_4M_P_oe->M22;
      MatrixXc<Real> Q11_oo = st_PQ->st_4M_Q_oe->M11;
      MatrixXc<Real> Q12_oe = st_PQ->st_4M_Q_oe->M12;
      MatrixXc<Real> Q21_eo = st_PQ->st_4M_Q_oe->M21;
      MatrixXc<Real> Q22_ee = st_PQ->st_4M_Q_oe->M22;

      size_t num_even = Q11_ee.rows();
      size_t num_odd = Q11_oo.rows();
      ArrayXi ind1_eo = st_PQ->st_4M_Q_eo->ind1;
      ArrayXi ind2_eo = st_PQ->st_4M_Q_eo->ind2;
      ArrayXi ind1_oe = st_PQ->st_4M_Q_oe->ind1;
      ArrayXi ind2_oe = st_PQ->st_4M_Q_oe->ind2;

      int m = st_PQ->st_4M_P_eo->m;

      if (!m) {
        st_TR->st_4M_R_oe->M11 = invertLUcol(Q11_oo);
        st_TR->st_4M_T_oe->M11 = -P11_oo * st_TR->st_4M_R_oe->M11;
        st_TR->st_4M_R_oe->M22 = invertLUcol(Q22_ee);
        st_TR->st_4M_T_oe->M22 = -P22_ee * st_TR->st_4M_R_oe->M22;
        st_TR->st_4M_R_eo->M11 = invertLUcol(Q11_ee);
        st_TR->st_4M_T_eo->M11 = -P11_ee * st_TR->st_4M_R_eo->M11;
        st_TR->st_4M_R_eo->M22 = invertLUcol(Q22_oo);
        st_TR->st_4M_T_eo->M22 = -P22_oo * st_TR->st_4M_R_eo->M22;

        st_TR->st_4M_T_eo->M12 = ArrayXXc<Real>::Zero(num_even, num_odd);
        st_TR->st_4M_T_eo->M21 = ArrayXXc<Real>::Zero(num_odd, num_even);
        st_TR->st_4M_T_oe->M12 = ArrayXXc<Real>::Zero(num_odd, num_even);
        st_TR->st_4M_T_oe->M21 = ArrayXXc<Real>::Zero(num_even, num_odd);

        if (get_R) {
          st_TR->st_4M_R_eo->M12 = ArrayXXc<Real>::Zero(num_even, num_odd);
          st_TR->st_4M_R_eo->M21 = ArrayXXc<Real>::Zero(num_odd, num_even);
          st_TR->st_4M_R_oe->M12 = ArrayXXc<Real>::Zero(num_odd, num_even);
          st_TR->st_4M_R_oe->M21 = ArrayXXc<Real>::Zero(num_even, num_odd);
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

        st_TR->st_4M_T_eo->M12 = G1 * G6 - G4;
        st_TR->st_4M_T_eo->M22 = G3 * G6 - G2;
        st_TR->st_4M_T_eo->M11 = -G1 - st_TR->st_4M_T_eo->M12 * G5;
        st_TR->st_4M_T_eo->M21 = -G3 - st_TR->st_4M_T_eo->M22 * G5;

        if (get_R) {
          st_TR->st_4M_R_eo->M12 = -Q11_inv * G6;
          st_TR->st_4M_R_eo->M22 = F2;
          st_TR->st_4M_R_eo->M11 = Q11_inv - st_TR->st_4M_R_eo->M12 * G5;
          st_TR->st_4M_R_eo->M21 = -st_TR->st_4M_R_eo->M22 * G5;
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

        st_TR->st_4M_T_oe->M12 = G1 * G6 - G4;
        st_TR->st_4M_T_oe->M22 = G3 * G6 - G2;
        st_TR->st_4M_T_oe->M11 = -G1 - st_TR->st_4M_T_oe->M12 * G5;
        st_TR->st_4M_T_oe->M21 = -G3 - st_TR->st_4M_T_oe->M22 * G5;

        if (get_R) {
          st_TR->st_4M_R_oe->M12 = -Q11_inv * G6;
          st_TR->st_4M_R_oe->M22 = F2;
          st_TR->st_4M_R_oe->M11 = Q11_inv - st_TR->st_4M_R_oe->M12 * G5;
          st_TR->st_4M_R_oe->M21 = -st_TR->st_4M_R_oe->M22 * G5;
        }
      }

      st_TR->st_4M_T_eo->m = st_TR->st_4M_T_oe->m = m;
      st_TR->st_4M_T_eo->ind1 = ind1_eo;
      st_TR->st_4M_T_eo->ind2 = ind2_eo;
      st_TR->st_4M_T_oe->ind1 = ind1_oe;
      st_TR->st_4M_T_oe->ind2 = ind2_oe;
      st_TR->mat_list.push_back("st_4M_T");

      if (get_R) {
        st_TR->st_4M_R_eo->m = st_TR->st_4M_R_oe->m = m;
        st_TR->st_4M_R_eo->ind1 = ind1_eo;
        st_TR->st_4M_R_eo->ind2 = ind2_eo;
        st_TR->st_4M_R_oe->ind1 = ind1_oe;
        st_TR->st_4M_R_oe->ind2 = ind2_oe;
        st_TR->mat_list.push_back("st_4M_R");
      }
    }
    return output;
  }
}
