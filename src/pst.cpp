#include "pst.h"
#include "vsh.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  unique_ptr<stRes<Real>> pstMakeStructForField(unique_ptr<stAbcdnm<Real>> st_abcdnm,
    int N_max, ArrayXr<Real> lambda, ArrayXr<Real> epsilon2, Real epsilon1,
    unique_ptr<stIncPar<Real>> inc_par, Real a, Real c) {
    if (!st_abcdnm->c_nm.size())
      cout << "Pb in  pstMakeStructForField: the structure stAbcdnm should contain c_nm and d_nm for internal fields" << endl;
    unique_ptr<stRes<Real>> output = make_unique<stRes<Real>>();
    output->N_max = N_max;
    output->lambda = lambda;
    output->epsilon1 = epsilon1;
    output->epsilon2 = epsilon2;
    output->inc_par = move(inc_par);
    if (!isnan(a))
      output->a = a;
    if (!isnan(c))
      output->a = c;
    return output;
  }

  template <class Real>
  unique_ptr<stRes<Real>> pstMakeStructForField(
      unique_ptr<stAbcdnm<Real>> st_abcdnm, unique_ptr<stParams<Real>> params) {
    if (!st_abcdnm->c_nm.size())
      cout << "Pb in  pstMakeStructForField: the structure stAbcdnm should contain c_nm and d_nm for internal fields" << endl;
    unique_ptr<stRes<Real>> output = make_unique<stRes<Real>>();
    output->N_max = params->N;
    output->lambda = params->lambda;
    output->epsilon1 = params->epsilon1;
    output->epsilon2 = params->epsilon2;
    output->a = params->a;
    output->c = params->c;
    if (params->inc_par)
      output->inc_par = move(params->inc_par);
    else {
      sIncType inc_type = params->inc_type;
      output->inc_par = vshMakeIncidentParams<Real>(inc_type, params->N);
    }
    return output;
  }

  template <class Real>
  Real CCGIN(int n, int n1, int m, int mm, ArrayXr<Real> F, ArrayXr<Real> sign) {
    int m1 = mm - m;
    if (n < abs(m) || n1 < abs(m1) || abs(mm) > n + n1) {
      cout << "ERROR IN SUBROUTINE CCGIN" << endl;
      return ArrayXXr<Real>::Zero(1, 1);
    }
    if (abs(mm) > abs(n - n1)) {
      if (n1 > n) {
        swap(n, n1);
        swap(m, m1);
      }
      int N2 = 2*n;
      int N12 = 2*n1;
      return sign(n1 + m1) * exp(F(n + m) + F(n - m) + F(N12) + F(N2 - N12 + 1)
          - F(N2 + 1) - F(n1 + m1) - F(n1 - m1) - F(n - n1 + mm) - F(n - n1 - mm));
    }
    int A = 1;
    if (mm < 0) {
      mm = -mm;
      m = -m;
      m1 = -m1;
      A = sign(mm + n + n1);
    }
    return A * sign(n1 + m1) * exp(F(2 * mm + 1) + F(n + n1 - mm) + F(n + m) + F(n1 + m1)
        - F(n + n1 + mm + 1) - F(n - n1 + mm) - F(-n + n1 + mm) - F(n - m) - F(n1 - m1));
  }

  template <class Real>
  Real CGDIRECT(int n, int m, int n1, int m1, ArrayXr<Real> F) {
    Real C = F(2*n) + F(2*n1) + F(n + n1 + m + m1) + F(n + n1 - m - m1);
    C -= F(2*(n + n1)) + F(n + m) + F(n - m) + F(n1 + m1) - F(n1 - m1);
    return exp(C);
  }

  template <class Real>
  ArrayXXr<Real> CCG(int n, int n1, int N_max, int K1, int K2, ArrayXr<Real> log_fact, ArrayXr<Real> sign) {
    if (n1 < 0 || n1 > N_max + n || n < 1 || n > N_max) {
      cout << "ERROR IN CCG" << endl;
      return ArrayXXr<Real>::Zero(1, 1);
    }
    ArrayXXr<Real> GG = ArrayXXr<Real>::Zero(2*N_max + 1, N_max + 1);
    ArrayXXr<Real> CD = ArrayXXr<Real>::Zero(2*N_max + 1, 2*N_max + 1);
    ArrayXXr<Real> CU = ArrayXXr<Real>::Zero(2*N_max + 1, 2*N_max + 1);

    int NPN6 = N_max;
    int NNF = min(n + n1, N_max);
    int m_init = (K1 == 1 && K2 == 0) ? NPN6 - n : NPN6;
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
      Real C = CCGIN(n, n1, m, mm, log_fact, sign);
      CU(NNL) = C;
      if (NNL != NNF) {
        Real C2 = 0;
        Real C1 = C;
        for (int nn = NNL + 1; nn <= min(NNM, NNF); nn++) {
          Real A = static_cast<Real>((nn + mm) * (nn - mm) * (n1 - n + nn));
          A *= (n - n1 + nn) * (n + n1 - nn + 1) * (n + n1 + nn + 1);
          A = (4 * nn * nn) / A;
          A *= (2*nn + 1) * (2*nn + 1);
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
              B = static_cast<Real>(2 * (nn + 2) * (nn + 1));
              B = ((2 * m - mm) * (nn + 2) * (nn + 1) - mm * n * (n + 1) + mm * n1 * (n1 + 1)) / B
              D = static_cast<Real>(4 * (nn + 2) * (nn + 2));
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

      for (int nn = NNL; nn <= NNF; nn ++) {
        if (nn <= NNM)
          GG(m_ind, nn) = CU(nn);
        else
          GG(m_ind, nn) = CD(nn);
      }
    }
    return GG;
  }

  template <class Real>
  unique_ptr<stSM<Real>> pstScatteringMatrixOA(vector<unique_ptr<stTR>> st_TR_list,
      Real lambda, const ArrayXr<Real>& sca, int Nb_theta) {
    int K1 = 0, K2 = 0, K3 = 0, K4 = 1, K5 = 1, K6 = 2;
    int N = st_TR_list.size() - 1;

    RowArrayXr<Real> n_vec = RowArrayXi::LinSpaced(N, 1, N);
    ArrayXr<Real> k_vec = ArrayXi::LinSpaced(N, 1, N);
    ArrayXXr<Real> FF_kn = sqrt((2*n_vec + 1).replicate(N, 1).colwise() / (2*k_vec + 1));
    ArrayXXc<Real> FA_kn = (pow(I, -n_vec).replicate(N, 1).colwise() * (pow(I, k_vec) / sqrt(2*k_vec + 1)));
    ArrayXr<Real> log_fact = ArrayXXr<Real>::Zero(4*(N + 1));
    for (int i = 2; i < 4*(N + 1); i++)
      log_fact(i) = log_fact(i - 1) + 0.5*log(i - 1);

    ArrayXi sign = pow(-1, ArrayXi::LinSpaced(4*(N + 1), 0, 4*N + 3));

    ArrayXXc<Real> T1_mk = ArrayXXc<Real>::Zero(2*N + 1, N);
    ArrayXXc<Real> T2_mk = ArrayXXc<Real>::Zero(2*N + 1, N);
    Tensor3c<Real> B1_n1_mn(N + 1, 2*N + 1, N);
    B2_n1_mn.setZero();
    Tensor3c<Real> B1_n1_mn(N + 1, 2*N + 1, N);
    B2_n1_mn.setZero();

    vector<ArrayXXc<Real>> CT11(N + 1, ArrayXXc<Real>::Zero(N, N));
    vector<ArrayXXc<Real>> CT12(N + 1, ArrayXXc<Real>::Zero(N, N));
    vector<ArrayXXc<Real>> CT21(N + 1, ArrayXXc<Real>::Zero(N, N));
    vector<ArrayXXc<Real>> CT22(N + 1, ArrayXXc<Real>::Zero(N, N));

    /**
    for (int m = 0; m <= N; m++) {
      CT11[m]()
    }
    */

    return unique_ptr<stSM<Real>>();
  }

  template unique_ptr<stRes<double>> pstMakeStructForField(unique_ptr<stAbcdnm<double>>,
      int, ArrayXr<double>, ArrayXr<double>, double, unique_ptr<stIncPar<double>>, double, double);
  template unique_ptr<stRes<double>> pstMakeStructForField(
      unique_ptr<stAbcdnm<double>>, unique_ptr<stParams<double>>);
}
