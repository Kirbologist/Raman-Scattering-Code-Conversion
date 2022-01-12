#include "pst.h"
#include "vsh.h"
#include "math.h"

using namespace Eigen;
using namespace std;

namespace Smarties {
  template <class Real>
  unique_ptr<stRes<Real>> pstMakeStructForField(const unique_ptr<stAbcdnm<Real>>& st_abcdnm,
    int N_max, ArrayXr<Real> lambda, ArrayXr<Real> epsilon2, Real epsilon1,
    unique_ptr<stIncPar<Real>> inc_par, Real a, Real c) {
    if (!st_abcdnm->c_nm.size())
      cout << "Pb in  pstMakeStructForField: the structure stAbcdnm should contain c_nm and d_nm for internal fields" << endl;
    unique_ptr<stRes<Real>> output = make_unique<stRes<Real>>(*st_abcdnm);
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

  template <class Real>
  unique_ptr<stRes<Real>> pstMakeStructForField(
      const unique_ptr< stAbcdnm<Real>>& st_abcdnm, const unique_ptr<stParams<Real>>& params) {
    if (!st_abcdnm->c_nm.size())
      cout << "Pb in  pstMakeStructForField: the structure stAbcdnm should contain c_nm and d_nm for internal fields" << endl;
    unique_ptr<stRes<Real>> output = make_unique<stRes<Real>>(*st_abcdnm);
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

  template <class Real>
  Real CCGIN(int n, int n1, int m, int mm, ArrayXr<Real> F, ArrayXi s_sign) {
    int m1 = mm - m;
    if (n < abs(m) || n1 < abs(m1) || abs(mm) > n + n1) {
      cout << "ERROR IN SUBROUTINE CCGIN" << endl;
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

  template <class Real>
  Real CGDIRECT(int n, int m, int n1, int m1, ArrayXr<Real> F) {
    Real C = F(2*n) + F(2*n1) + F(n + n1 + m + m1) + F(n + n1 - m - m1);
    C -= F(2*(n + n1)) + F(n + m) + F(n - m) + F(n1 + m1) + F(n1 - m1);
    return exp(C);
  }

  template <class Real>
  ArrayXXr<Real> CCG(int n, int n1, int N_max, int K1, int K2, ArrayXr<Real> log_fact, ArrayXi s_sign) {
    if (n1 < 0 || n1 > N_max + n || n < 0 || n > N_max) {
      cout << "ERROR IN CCG" << endl;
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

  template <class Real>
  unique_ptr<stSM<Real>> pstScatteringMatrixOA(const vector<unique_ptr<stTR<Real>>>& st_TR_list,
      Real lambda, Real sca, int Nb_theta) {
    int K1 = 1, K2 = 0, K3 = 0, K4 = 1, K5 = 1, K6 = 2;
    int N1 = st_TR_list.size();
    int N = N1 - 1;

    RowArrayXr<Real> n_vec = RowArrayXr<Real>::LinSpaced(N1, 0, N);
    ArrayXr<Real> k_vec = ArrayXr<Real>::LinSpaced(N1, 0, N);
    ArrayXXr<Real> FF_kn = sqrt((2*n_vec + 1).replicate(N1, 1).colwise() / (2*k_vec + 1));
    ArrayXXc<Real> FA_kn = pow(mp_im_unit<Real>(), -n_vec).replicate(N1, 1).colwise() *
        (pow(mp_im_unit<Real>(), k_vec) / sqrt(2*k_vec + 1));
    ArrayXr<Real> log_fact = ArrayXr<Real>::Zero(4*N1);
    for (int i = 2; i < 4*N1; i++)
      log_fact(i) = log_fact(i - 1) + 0.5*log(i);

    ArrayXi s_sign = pow(-1, ArrayXi::LinSpaced(4*N1, 0, 4*N1 - 1));

    ArrayXXc<Real> T1_mk = ArrayXXc<Real>::Zero(2*N + 1, N1);
    ArrayXXc<Real> T2_mk = ArrayXXc<Real>::Zero(2*N + 1, N1);
    Tensor3c<Real> B1_n1_mn(2*N + 1, 2*N + 1, N1);
    B1_n1_mn.setZero();
    Tensor3c<Real> B2_n1_mn(2*N + 1, 2*N + 1, N1);
    B2_n1_mn.setZero();

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

    for (int n = 0; n <= N; n++) {
      for (int nn = 0; nn <= N; nn++) {
        int m_max = min(n, nn);
        for (int m = 0; m <= m_max; m++) {
          int m_ind = N + m;
          complex<Real> TT1 = CT11[m](n, nn);
          complex<Real> TT2 = CT12[m](n, nn);
          complex<Real> TT3 = CT21[m](n, nn);
          complex<Real> TT4 = CT22[m](n, nn);
          complex<Real> T1 = TT1 + TT2;
          complex<Real> T2 = TT3 + TT4;
          T1_mk(m_ind, nn) = T1 + T2;
          T2_mk(m_ind, nn) = T1 - T2;
          if (m > 0) {
            T1 = TT1 - TT2;
            T2 = TT3 - TT4;
            m_ind = N - m;
            T1_mk(m_ind, nn) = T1 - T2;
            T2_mk(m_ind, nn) = T1 + T2;
          }
        }
      }

      int nn1_max = N + n;
      for (int n1 = 0; n1 <= nn1_max; n1++) {
        ArrayXXr<Real> G1 = CCG(n, n1, N, K1, K2, log_fact, s_sign);
        int nn_max = min(N, n1 + n);
        int nn_min = max(1, abs(n - n1));
        int kn = n + n1;
        ArrayXc<Real> A1k = ArrayXc<Real>::Zero(N1);
        ArrayXc<Real> A2k = ArrayXc<Real>::Zero(N1);
        for (int nn = nn_min; nn <= nn_max; nn++) {
          int SIG = s_sign(kn + nn);
          int m_ind_max = min(n, nn) + N;
          complex<Real> AA1 = 0;
          complex<Real> AA2 = 0;
          for (int m_ind = N; m_ind <= m_ind_max; m_ind++) {
            int m = m_ind - N;
            complex<Real> SSS = G1(m_ind, nn);
            complex<Real> R1 = T1_mk(m_ind, nn);
            complex<Real> R2 = T2_mk(m_ind, nn);
            if (m > 0) {
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

        ArrayXXr<Real> G2 = CCG(n, n1, N, K3, K4, log_fact, s_sign);
        int m_max = min(n1 + 1, n);
        int m_min = max(-n1 + 1, -n);
        for (int m = m_min; m <= m_max; m++) {
          int m_ind = m + N;
          complex<Real> BB1 = 0;
          complex<Real> BB2 = 0;
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

    for (int n = 0; n <= N; n++) {
      for (int nn = 0; nn <= N; nn++) {
        int m_ind = min(n, nn);
        int m_ind_max = N + m_ind;
        int m_ind_min = N - m_ind;
        int n1_max = m_ind_max;
        for (m_ind = m_ind_min; m_ind <= m_ind_max; m_ind++) {
          int m = m_ind - N;
          int n1_min = abs(m - 1);
          Real DD1 = 0;
          Real DD2 = 0;
          for (int n1 = n1_min; n1 <= n1_max; n1++) {
            int XX = 2*n1 + 1;
            DD1 += XX * real(B1_n1_mn(n1, m_ind, n) * conj(B1_n1_mn(n1, m_ind, nn)));
            DD2 += XX * real(B2_n1_mn(n1, m_ind, n) * conj(B2_n1_mn(n1, m_ind, nn)));
          }
          D1_m_kn(m_ind, nn, n) = DD1;
          D2_m_kn(m_ind, nn, n) = DD2;
        }
        int m_max = min(n, nn + 2);
        int m_min = max(-n, -nn + 2);
        m_ind_max = N + m_max;
        m_ind_min = N + m_min;
        for (int m_ind = m_ind_min; m_ind <= m_ind_max; m_ind++) {
          int m = m_ind - N;
          int n1_min = abs(m - 1);
          Real DD3 = 0;
          Real DD4 = 0;
          complex<Real> DD5 = 0;
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

    unique_ptr<stSM<Real>> output = make_unique<stSM<Real>>();
    output->ALF1n = ArrayXr<Real>::Zero(2*N + 1);
    output->ALF2n = ArrayXr<Real>::Zero(2*N + 1);
    output->ALF3n = ArrayXr<Real>::Zero(2*N + 1);
    output->ALF4n = ArrayXr<Real>::Zero(2*N + 1);
    output->BET1n = ArrayXr<Real>::Zero(2*N + 1);
    output->BET2n = ArrayXr<Real>::Zero(2*N + 1);

    int L_max = 0;
    Real DK = pow(lambda, 2) / (4*sca*mp_pi<Real>());
    for (int l = 0; l < 2*N + 1; l++) {
      Real G1L = 0;
      Real G2L = 0;
      Real G3L = 0;
      Real G4L = 0;
      complex<Real> G5L = 0;
      Real SL = (2*l + 1)*DK;
      for (int n = 0; n <= N; n++) {
        int nn_min = max(1, abs(n - l));
        int nn_max = min(N, n + l);
        if (nn_min <= nn_max) {
          ArrayXXr<Real> G1 = CCG(n, l, N, K1, K2, log_fact, s_sign);
          ArrayXXr<Real> G2;
          if (l >= 2)
            G2 = CCG(n, l, N, K5, K6, log_fact, s_sign);
          for (int nn = nn_min; nn <= nn_max; nn++) {
            int m_max = min(n, nn);
            int m_ind_min = N - m_max;
            int m_ind_max = N + m_max;
            int SI = s_sign(n + l + nn);
            Real DM1 = 0;
            Real DM2 = 0;
            for (int m_ind = m_ind_min; m_ind <= m_ind_max; m_ind++) {
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
            Real FFN = FF_kn(nn, n);
            Real SSS = G1(N + 1, nn) * FFN;
            G1L += SSS * DM1;
            G2L += SSS * DM2 * SI;
            if (l >= 2) {
              Real DM3 = 0;
              Real DM4 = 0;
              complex<Real> DM5 = 0;
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
            }
          }
        }
      }
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

      L_max = l;
      if (abs(G1L) < mp_eps<Real>())
        break;
    }
    output->L_max = L_max;
    output->asym_par = output->ALF1n(1) / 3;

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
    Real DA = DN * mp_pi<Real>();
    Real DB = DN * 180;
    Real TB = -DB;
    Real TAA = -DA;
    Real D6 = sqrt(static_cast<Real>(6.0)) / 4;

    for (int I1 = 0; I1 < Nb_theta; I1++) {
      TAA += DA;
      TB += DB;
      Real U = cos(TAA);
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
      }
      Real F22 = (F2 + F3) / 2;
      Real F33 = (F2 - F3) / 2;
      output->F11(I1) = F11;
      output->F22(I1) = F22;
      output->F33(I1) = F33;
      output->F44(I1) = F44;
      output->F12(I1) = F12;
      output->F34(I1) = F34;
    }

    output->all_AB = Array<Real, Dynamic, 7>(2*N + 1, 7);
    output->all_AB << ArrayXr<Real>::LinSpaced(2*N + 1, 1, 2*N + 1),
        output->ALF1n, output->ALF2n, output->ALF3n, output->ALF4n, output->BET1n, output->BET2n;

    output->all_F = Array<Real, Dynamic, 7>(Nb_theta, 7);
    output->all_F << theta_deg, output->F11, output->F22, output->F33, output->F44,
        output->F12, output->F34;

    return output;
  }

  template unique_ptr<stRes<double>> pstMakeStructForField(const unique_ptr<stAbcdnm<double>>&,
      int, ArrayXd, ArrayXd, double, unique_ptr<stIncPar<double>>, double, double);
  template unique_ptr<stRes<double>> pstMakeStructForField(
      const unique_ptr<stAbcdnm<double>>&, const unique_ptr<stParams<double>>&);
  template unique_ptr<stSM<double>> pstScatteringMatrixOA(const vector<unique_ptr<stTR<double>>>&,
      double, double, int);
}
