#include "sph.h"
#include "vsh.h"
#include "math.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  ArrayXXr<Real>* sphGetUforFp(int n) {
    int b_max = n/2;
    ArrayXXr<Real>* u = new ArrayXXr<Real>(b_max + 2, b_max + 1);
    (*u).setZero();
    (*u)(0, 0) = 1;
    for (int b = 0; b < b_max; b++) {
      (*u)(0, b + 1) = (2*n - 1)*(*u)(0, b) - (2*n - 1)*(*u)(1, b);
      for (int r = 1; r <= b + 1; r++) {
        (*u)(r, b + 1) = (2*n - 1 - 4*r)*(*u)(r, b) - (2*n - 1 - 2*r)*(*u)(r + 1, b) + 2*r*(*u)(r - 1, b);
      }
    }
    return u;
  }

  template <class Real>
  stFprow<Real>* sphGetFpRow(int n, Real s, const ArrayXr<Real>& x) {
    int num_x = x.size();
    RowArrayXr<Real> x_squared = x.pow(2).transpose();
    ArrayXXr<Real>* u = sphGetUforFp<Real>(n);
    Real alpha_bar_k = 1;
    int q_min;
    int q_int;

    ArrayXXr<Real> S = ArrayXXr<Real>::Zero(n - 3, num_x);
    ArrayXXr<Real> max_term_S = ArrayXXr<Real>::Zero(n - 3, num_x);
    ArrayXXr<Real> loss_prec_S = ArrayXXr<Real>::Zero(n - 3, num_x);

    // Original code dynamically increased the array size as needed
    ArrayXXr<Real> beta = ArrayXXr<Real>::Zero(n + 1 - (n % 2), n - (n % 2));
    beta(0, 0) = 1;
    beta(1, 0) = beta(0, 0) * (s - 1) * (s + 1);
    beta(0, 1) = 1;
    beta(1, 1) = beta(0, 1) * (s - 1) * (s + 1) * 2;
    beta(2, 1) = beta(1, 1) * (s - 1) * (s + 1) / 2;

    RowArrayXr<Real> alpha_q(num_x), current_term(num_x);
    RowArrayXi counter(num_x);
    RowArrayXb x_not_converged(num_x), test(num_x);
    int b, num_to_test = 3;
    Real gamma_qnk;

    for (int k = n - 4; k >= 0; k -= 2) {
      q_min = (n - k)/2 - 1; // expression is mathematically always an integer
      q_int = n - k;
      alpha_bar_k = -alpha_bar_k/(n - k - 2);
      beta(0, q_int - 2) = 1;
      for (int j = 1; j < q_int; j++)
        beta(j, q_int - 2) = beta(j - 1, q_int - 2)*(s - 1)*(s + 1)*(q_int - j)/j;
      beta(0, q_int - 1) = 1;
      for (int j = 1; j <= q_int; j++)
        beta(j, q_int - 1) = beta(j - 1, q_int - 1)*(s - 1)*(s + 1)*(q_int - j + 1)/j;

      alpha_q.fill(alpha_bar_k);
      x_not_converged.fill(true);
      current_term.setZero();
      counter.setZero();
      test.fill(false);

      for (int q = q_min; q < q_int; q++) {
        b = n - k - q - 1;
        gamma_qnk = 0;
        for (int j = max(0, 2 * q + k - n + 1); j <= q; j++)
          gamma_qnk += beta(j, q - 1) * (*u)(q - j, b);

        if (abs(s) < 2) {
          current_term = alpha_q * gamma_qnk;
          S.row(k) += current_term;
          max_term_S.row(k) = max_term_S.row(k).max(abs(current_term));
        } else {
          test.fill(false);
          for (int i = 0; i < num_x; i++) {
            if (x_not_converged(i)) {
              current_term(i) = alpha_q(i) * gamma_qnk;
              test(i) = abs(current_term(i)) > EPS * abs(S(k, i));
            }
          }
          if (test.any()) {
            if ((1 - test + current_term.isInf()).all()) // if all non-converged x values are infinite
              cout << "Problem (1) in sphGetFpovx..." << endl;
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
        }

        alpha_q *= -x_squared / (2*(q + 1));
      }

      x_not_converged.fill(true);
      counter.setZero();

      Real c_qqnk = static_cast<Real>(1)/(2*k + 1), c_iqnk;
      int q = q_int;
      while (true) {
        gamma_qnk = c_iqnk = c_qqnk;
        for (int i = q - 1; i >= 0; i--) {
          c_iqnk *= pow(s, 2)*(i + 1)*(2*i + 1 - 2*n)/(q - i)/(2*k + 2*q - 2*i + 1);
          gamma_qnk += c_iqnk;
        }

        test.fill(false);
        for (int i = 0; i < num_x; i++) {
          if (x_not_converged(i)) {
            current_term(i) = alpha_q(i) * gamma_qnk;
            test(i) = abs(current_term(i)) > EPS * abs(S(k, i));
          }
        }
        if (test.any()) {
          if ((1 - test + current_term.isInf()).all()) // if all non-converged x values are infinite
            cout << "Problem (2) in sphGetFpovx..." << endl;
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

        for (int i = 0; i < num_x; i++) {
          if (x_not_converged(i))
            alpha_q(i) *= -x_squared(i) / (2*(q + 1));
        }

        c_qqnk /= (2*q + 1 - 2*n);
        q++;
      }

      loss_prec_S.row(k) = abs(max_term_S.row(k) / S.row(k));
      S.row(k) *= -pow(s, k + 1);
    }
    delete u;

    stFprow<Real>* output = new stFprow<Real>();
    output->S = S;
    output->loss_prec_S = loss_prec_S;
    return output;
  }

  template <class Real>
  stFpovx<Real>* sphGetFpovx(int n_n_max, Real s, const ArrayXr<Real>& x) {
    int num_x = x.size();

    stFpovx<Real>* output = new stFpovx<Real>();
    output->rb_chi = vshRBchi<Real>(ArrayXr<double>::LinSpaced(n_n_max + 1, 0, n_n_max), x);
    output->rb_psi = vshRBpsi<Real>(ArrayXr<double>::LinSpaced(n_n_max + 1, 0, n_n_max), s*x);
    output->Fpovx = new ArrayXXc<Real>[n_n_max + 1]();
    for (int i = 0; i <= n_n_max; i++)
      output->Fpovx[i] = ArrayXXc<Real>::Zero(n_n_max + 1, num_x);

    stFprow<Real>* FpRow = sphGetFpRow<Real>(n_n_max, s, x);
    output->Fpovx[n_n_max].topRows(n_n_max - 3) = (FpRow->S).rowwise() / x.transpose();
    delete FpRow;

    for (int k = n_n_max; k >= 0; k--) {
      for (int i = k % 2; i <= min(k + 2, n_n_max); i += 2)
        (output->Fpovx[i]).row(k) = (output->rb_chi.col(i) * output->rb_psi.col(k)).transpose();
      if (k > 0) {
        for (int n = k + 3; n < n_n_max; n += 2)
          (output->Fpovx[n]).row(k - 1) = ((output->Fpovx[n + 1]).row(k) + (output->Fpovx[n - 1]).row(k)) *
              (2*k + 1) / (2*n + 1) / s - (output->Fpovx[n]).row(k + 1);
      }
    }
    return output;
  }

  template <class Real>
  stBessel<Real>* sphGetXiPsi(int n_n_max, Real s, const ArrayXr<Real>& x, int N_B) {
    stBessel<Real>* output = new stBessel<Real>();
    stFpovx<Real>* chi_psi_prod = sphGetFpovx<Real>(N_B + 1, s, x);
    output->psi_psi = new ArrayXXc<Real>[n_n_max + 2]();
    for (int i = 0; i < n_n_max + 2; i++)
      output->psi_psi[i] = ArrayXXc<Real>::Zero(n_n_max + 2, x.size());
    output->psi_n = vshRBpsi<Real>(ArrayXr<Real>::LinSpaced(n_n_max + 2, 0, n_n_max + 1), x);

    for (int n = 0; n < n_n_max + 2; n++) {
      for (int k = n % 2; k < n_n_max + 2; k += 2)
        output->psi_psi[n].row(k) = (chi_psi_prod->rb_psi.col(k) * output->psi_n.col(n)).transpose();
    }

    output->xi_psi = new ArrayXXc<Real>[n_n_max + 2]();
    for (int i = 0; i < n_n_max + 2; i++)
      output->xi_psi[i] = output->psi_psi[i] + I*chi_psi_prod->Fpovx[i](seq(0, n_n_max + 1), all);
    output->chi_n = chi_psi_prod->rb_chi.leftCols(n_n_max + 2);
    output->psi_k = chi_psi_prod->rb_psi.leftCols(n_n_max + 2);
    delete[] chi_psi_prod->Fpovx;
    delete chi_psi_prod;
    return output;
  }

  template <class Real>
  stBesselPrimes<Real>* sphGetBesselProductsPrimes(ArrayXXc<Real>* prods, int N) {
    int X = prods[0].cols();

    stBesselPrimes<Real>* output = new stBesselPrimes<Real>();
    output->xi_psi = new ArrayXXc<Real>[N]();
    output->xi_prime_psi = new ArrayXXc<Real>[N]();
    output->xi_psi_prime = new ArrayXXc<Real>[N]();
    output->xi_prime_psi_prime = new ArrayXXc<Real>[N]();
    output->xi_prime_psi_prime_plus_nnp1_xi_psi_over_ssx = new ArrayXXc<Real>[N]();
    output->xi_prime_psi_prime_plus_kkp1_xi_psi_over_ssx = new ArrayXXc<Real>[N]();
    output->xi_psi_over_sxx = new ArrayXXc<Real>[N]();
    for (int i = 0; i < N; i++) {
      output->xi_psi[i] = prods[i + 1](seq(1, last - 1), all);
      output->xi_prime_psi[i] = ArrayXXc<Real>::Zero(N, X);
      output->xi_psi_prime[i] = ArrayXXc<Real>::Zero(N, X);
      output->xi_prime_psi_prime[i] = ArrayXXc<Real>::Zero(N, X);
      output->xi_prime_psi_prime_plus_nnp1_xi_psi_over_ssx[i] = ArrayXXc<Real>::Zero(N, X);
      output->xi_prime_psi_prime_plus_kkp1_xi_psi_over_ssx[i] = ArrayXXc<Real>::Zero(N, X);
      output->xi_psi_over_sxx[i] = ArrayXXc<Real>::Zero(N, X);
    }

    int k_p1, n_p1;
    for (int n = 1; n <= N; n++) {
      for (int k = 2 - n % 2; k <= N; k += 2) {
        k_p1 = k*(k + 1);
        n_p1 = n*(n + 1);

        output->xi_prime_psi_prime_plus_kkp1_xi_psi_over_ssx[n - 1].row(k - 1) = (
            (k + n + 1)*(k + 1)*prods[n - 1].row(k - 1) +
            (k_p1 - k*(n + 1))*prods[n - 1].row(k + 1) +
            (k_p1 - (k + 1)*n)*prods[n + 1].row(k - 1) +
            (k_p1 + k*n)*prods[n + 1].row(k + 1)) / ((2*n + 1)*(2*k + 1));

        output->xi_prime_psi_prime_plus_nnp1_xi_psi_over_ssx[n - 1].row(k - 1) = (
            (n_p1 + (k + 1)*(n + 1))*prods[n - 1].row(k - 1) +
            (n - k)*(n + 1)*prods[n - 1].row(k + 1) +
            (n - k)*n*prods[n + 1].row(k - 1) +
            (n_p1 + k*n)*prods[n + 1].row(k + 1)) / ((2*n + 1)*(2*k + 1));
      }
      for (int k = 1 + n % 2; k < N; k += 2) {
        output->xi_prime_psi[n - 1].row(k - 1) = (prods[n - 1].row(k)*(n + 1) - n*prods[n + 1].row(k)) / (2 * n + 1);
        output->xi_psi_prime[n - 1].row(k - 1) = (prods[n].row(k - 1)*(k + 1) - k*prods[n].row(k + 1)) / (2 * k + 1);
      }
    }
    return output;
  }

  template <class Real>
  stBesselProducts<Real>* sphGetModifiedBesselProducts(int n_n_max, Real s, const ArrayXr<Real>& x, int N_B) {
    stBesselProducts<Real>* output = new stBesselProducts<Real>();
    stBessel<Real>* prods = sphGetXiPsi<Real>(n_n_max, s, x, N_B);
    output->st_xi_psi_all = sphGetBesselProductsPrimes<Real>(prods->xi_psi, n_n_max);
    output->st_psi_psi_all = sphGetBesselProductsPrimes<Real>(prods->psi_psi, n_n_max);

    ArrayXXc<Real> psi_n_p1_psi_n = (prods->psi_n(all, seq(2, last)) * prods->psi_k(all, seq(1, last - 1))).transpose();
    ArrayXXc<Real> psi_n_psi_n_p1 = (prods->psi_n(all, seq(1, last - 1)) * prods->psi_k(all, seq(2, last))).transpose();
    ArrayXXc<Real> xi_n_p1_psi_n = psi_n_p1_psi_n + I*(prods->chi_n(all, seq(2, last)) * prods->psi_k(all, seq(1, last - 1))).transpose();
    ArrayXXc<Real> xi_n_psi_n_p1 = psi_n_psi_n_p1 + I*(prods->chi_n(all, seq(1, last - 1)) * prods->psi_k(all, seq(2, last))).transpose();
    output->st_psi_psi_all->for_diag_Lt1 = s*psi_n_psi_n_p1 - psi_n_p1_psi_n;
    output->st_xi_psi_all->for_diag_Lt1 = s*xi_n_psi_n_p1 - xi_n_p1_psi_n;
    ArrayXXc<Real> psi_n_psi_n = (prods->psi_n(all, seq(1, last - 1)) * prods->psi_k(all, seq(1, last - 1))).transpose();
    ArrayXXc<Real> xi_n_psi_n = psi_n_psi_n + I*(prods->chi_n(all, seq(1, last - 1)) * prods->psi_k(all, seq(1, last - 1))).transpose();

    delete[] prods->xi_psi;
    delete[] prods->psi_psi;
    delete prods;

    ArrayXXc<Real> n_vec = ArrayXc<Real>::LinSpaced(n_n_max, 2, n_n_max + 1);
    output->st_psi_psi_all->for_diag_Lt2 = psi_n_psi_n_p1 - s*psi_n_p1_psi_n + (s - 1)*(s + 1)/s *
        (n_vec.matrix() * (1/x.transpose()).matrix()).array() * psi_n_psi_n;
    output->st_xi_psi_all->for_diag_Lt2 = xi_n_psi_n_p1 - s*xi_n_p1_psi_n + (s - 1)*(s + 1)/s *
        (n_vec.matrix() * (1/x.transpose()).matrix()).array() *xi_n_psi_n;
    RowArrayXc<Real> tmp = s*x.pow(2).transpose();
    output->st_psi_psi_all->for_diag_Lt3 = psi_n_psi_n.rowwise() / tmp;
    output->st_xi_psi_all->for_diag_Lt3 = xi_n_psi_n.rowwise() / tmp;

    return output;
  }

  template <class Real>
  stRtfunc<Real>* sphMakeGeometry(size_t n_Nb_theta, Real a, Real c, ArrayXr<Real>* theta) {
    sInt type;
    stRtfunc<Real>* output = new stRtfunc<Real>();
    if (theta != nullptr) {
      output->w_theta = ArrayXr<Real>::Zero(theta->size());
      output->theta = *theta;
      output->n_Nb_theta = theta->size();
      type = PTS;
    } else {
      type = GAUSS;
      stRtfunc<Real>* tmp = auxPrepareIntegrals<Real>(2*n_Nb_theta, type);
      output->theta = tmp->theta(seq(0, n_Nb_theta - 1));
      output->w_theta = tmp->w_theta(seq(0, n_Nb_theta - 1))*2;
      output->n_Nb_theta = n_Nb_theta;
      delete tmp;
    }

    output->a = a;
    output->c = c;
    output->type = type;

    ArrayXr<Real> sin_t = sin(output->theta), cos_t = cos(output->theta);

    output->r = a*c/sqrt(pow(c*sin_t, 2) + pow(a*cos_t, 2));
    output->dr_dt = (pow(a, 2) - pow(c, 2))/pow(a*c, 2)*sin_t * cos_t * output->r.pow(3);

    output->h = max(a, c)/min(a, c);
    output->r0 = pow(pow(a, 2)*c, static_cast<Real>(1)/3);

    return output;
  }

  template ArrayXXr<double>* sphGetUforFp(int);
  template stFprow<double>* sphGetFpRow(int, double, const ArrayXr<double>&);
  template stFpovx<double>* sphGetFpovx(int, double, const ArrayXr<double>&);
  template stBessel<double>* sphGetXiPsi(int, double, const ArrayXr<double>&, int);
  template stBesselPrimes<double>* sphGetBesselProductsPrimes(ArrayXXc<double>*, int);
  template stBesselProducts<double>* sphGetModifiedBesselProducts(int, double, const ArrayXr<double>&, int);
  template stRtfunc<double>* sphMakeGeometry(size_t, double, double, ArrayXr<double>*);
}
