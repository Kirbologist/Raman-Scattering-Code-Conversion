#include "low_level.h"
#include "math.h"
#include "high_level.h"

using namespace std;
using namespace Eigen;

namespace Raman{
  template <class Real>
  typename LowLevel<Real>::stGLQuad* LowLevel<Real>::auxInitLegendreQuad(
      size_t N1, Real a, Real b) {
    int N = N1 - 1, N2 = N1 + 1;
    ArrayXr<Real> xu, y, y0(N1), Lp;
    ArrayXXr<Real> L(N1, N2);

    xu = ArrayXr<Real>::LinSpaced(N1, -1, 1);
    y = cos((2*ArrayXr<Real>::LinSpaced(N1, 0, N) + 1)*PI/(2*N + 2)) +
        (0.27/N1)*sin(PI*xu*N/N2);
    y0.fill(2);
    L.col(0).setOnes();

    int n_iter = 0;
    while ((n_iter < 15) && EPS <
        abs(acos(y) - acos(y0.template cast<complex<Real>>())).maxCoeff()) {
      n_iter++;
      L.col(1) = y;

      for (size_t i = 1; i < N1; i++)
        L.col(i + 1) = ((2*i + 1)*y * L.col(i) - i*L.col(i - 1))/(i + 1);

      Lp = N2*(L.col(N) - y*L.col(N1))/(1 - y.pow(2));

      y0 = y;
      y = y0 - L.col(N1)/Lp;
    }

    stGLQuad* output = new stGLQuad();

    output->x = (a*(1 - y) + b*(1 + y))/2;
    output->w = ArrayXr<Real>::Constant(N1, b - a) /
        ((1 - y.pow(2))*Lp.pow(2)) * pow(static_cast<Real>(N2) / N1, 2);
    return output;
  }

  // Missing functionality
  template <class Real>
  typename LowLevel<Real>::stRtfunc* LowLevel<Real>::auxPrepareIntegrals(
    size_t nNint, sInt type) {
    stRtfunc* output = new stRtfunc();
    output->nNbTheta = nNint;

    switch (type) {
      case GAUSS: {
        stGLQuad* weights = LowLevel<Real>::auxInitLegendreQuad(nNint);
        output->theta = acos(weights->x);
        output->wTheta = weights->w;
        delete weights;
        break;
      }

      case RECTANGLE: {
        output->theta = ArrayXr<Real>::LinSpaced(nNint, 0, PI);
        Real dTheta = PI/(nNint - 1);
        output->wTheta = dTheta*sin(output->theta);
        break;
      }

      default:
        cout << "Integration type not recognized" << endl;
    }

    return output;
  }

  template <class Real>
  typename LowLevel<Real>::stFpovx* LowLevel<Real>::sphGetFpovx(
      int n_n_max, Real s, ArrayXr<Real>& x) {
    int num_x = x.size();

    stFpovx* output = new stFpovx();
    output->rb_chi = vshRBchi(ArrayXr<double>::LinSpaced(n_n_max + 1, 0, n_n_max), x);
    output->rb_psi = vshRBpsi(ArrayXr<double>::LinSpaced(n_n_max + 1, 0, n_n_max), s*x);
    output->Fpovx = new ArrayXXc<Real>[n_n_max + 1]();
    for (int i = 0; i <= n_n_max; i++)
      output->Fpovx[i] = ArrayXXc<Real>::Zero(n_n_max + 1, num_x);

    stFprow* FpRow = sphGetFpRow(n_n_max, s, x);
    output->Fpovx[n_n_max].topRows(n_n_max - 3) = (FpRow->S).rowwise() / x.transpose();
    delete FpRow;

    for (int k = n_n_max; k >= 0; k--) {
      for (int i = k % 2; i <= min(k + 2, n_n_max); i += 2)
        (output->Fpovx[i]).row(k) = (output->rb_chi->col(i) * output->rb_psi->col(k)).transpose();
      if (k > 0) {
        for (int n = k + 3; n < n_n_max; n += 2)
          (output->Fpovx[n]).row(k - 1) = ((output->Fpovx[n + 1]).row(k) + (output->Fpovx[n - 1]).row(k)) *
              (2*k + 1) / (2*n + 1) / s - (output->Fpovx[n]).row(k + 1);
      }
    }
    return output;
  }

  template <class Real>
  typename LowLevel<Real>::stFprow* LowLevel<Real>::sphGetFpRow(
      int n, Real s, ArrayXr<Real>& x) {
    int num_x = x.size();
    RowArrayXr<Real> x_squared = pow(x, 2).transpose();
    ArrayXXr<Real>* u = sphGetUforFp(n);
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
            if ((1 - test + current_term.isInf()).all()) // if all non-converged x values are infinite...
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
          if ((1 - test + current_term.isInf()).all()) // if all non-converged x values are infinite...
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

    stFprow* output = new stFprow();
    output->S = S;
    output->loss_prec_S = loss_prec_S;
    return output;
  }

  template <class Real>
  ArrayXXr<Real>* LowLevel<Real>::sphGetUforFp(int n) {
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
  typename LowLevel<Real>::stBessel* LowLevel<Real>::sphGetXiPsi(
      int n_n_max, Real s, ArrayXr<Real>& x, int N_B) {
    stBessel* output = new stBessel();
    stFpovx* chi_psi_prod = sphGetFpovx(N_B + 1, s, x);
    output->psi_psi = new ArrayXXc<Real>[n_n_max + 2]();
    for (int i = 0; i < n_n_max + 2; i++)
      output->psi_psi[i] = ArrayXXc<Real>::Zero(n_n_max + 2, x.size());
    output->psi_n = vshRBpsi(ArrayXr<Real>::LinSpaced(n_n_max + 2, 0, n_n_max + 1), x);

    for (int n = 0; n < n_n_max + 2; n++) {
      for (int k = n % 2; k < n_n_max + 2; k += 2)
        output->psi_psi[n].row(k) = (chi_psi_prod->rb_psi->col(k) * output->psi_n->col(n)).transpose();
    }

    output->xi_psi = new ArrayXXc<Real>[n_n_max + 2]();
    for (int i = 0; i < n_n_max + 2; i++)
      output->xi_psi[i] = output->psi_psi[i] + I*chi_psi_prod->Fpovx[i](seq(0, n_n_max + 1), all);
    output->chi_n = new ArrayXXc<Real>();
    output->psi_k = new ArrayXXc<Real>();
    *(output->chi_n) = chi_psi_prod->rb_chi->leftCols(n_n_max + 2);
    *(output->psi_k) = chi_psi_prod->rb_psi->leftCols(n_n_max + 2);
    delete[] chi_psi_prod->Fpovx;
    delete chi_psi_prod->rb_chi;
    delete chi_psi_prod->rb_psi;
    return output;
  }

  // Many versions in original code
  template <class Real>
  typename LowLevel<Real>::stEAllPhi* LowLevel<Real>::vshEgenThetaAllPhi(
      ArrayXr<Real>& lambda, ArrayXr<Real>& epsilon, ArrayXXc<Real>& p_nm,
      ArrayXXc<Real>& q_nm, RowArrayXr<Real>& rt, RowArrayXr<Real>& theta,
      sBessel type, stPinmTaunm* stPT) {
    int n_p_max = p_nm.cols(), n_n_max = static_cast<int>(round(sqrt(n_p_max) - 1)),
        n_nb_lambda = lambda.size();
    if (rt.size() != theta.size() && rt(0) != 0 && !isinf(rt(0)))
      cout << "vshEgenThetaAllPhi error: theta and rt must be the same size row arrays." << endl;
    int n_nb_theta = theta.size();

    ArrayXr<Real> n = ArrayXr<Real>::LinSpaced(n_n_max + 1, 0, n_n_max);
    ArrayXr<Real> mu_n_times = sqrt((2*n + 1)*n*(n + 1)/(4*PI));
    ArrayXr<Real> mu_n_divd_gen = mu_n_times/(n*(n + 1));

    stEAllPhi* output = new stEAllPhi();
    output->theta = theta;
    output->r_of_theta = rt;
    ArrayXXc<Real>* CErm = new ArrayXXc<Real>[2*n_n_max + 1]();
    ArrayXXc<Real>* CEtm = new ArrayXXc<Real>[2*n_n_max + 1]();
    ArrayXXc<Real>* CEfm = new ArrayXXc<Real>[2*n_n_max + 1]();

    if (!rt(0)) {
      for (int m = -n_n_max; m <= n_n_max; m++) {
        if (abs(m) > 1) {
          CErm[m + n_n_max] = CEtm[m + n_n_max] = CEfm[m + n_n_max] =
              ArrayXc<Real>::Zero(n_nb_lambda, n_nb_lambda);
        }
      }
      Real coeff1 = 1/sqrt(6*PI), coeff2 = coeff1/sqrt(2);

      CErm[n_n_max] = ((coeff1 * q_nm.col(1)).matrix() * cos(theta).matrix()).array();
      CEtm[n_n_max] = ((-coeff1 * q_nm.col(1)).matrix() * sin(theta).matrix()).array();
      CEfm[n_n_max] = ArrayXc<Real>::Zero(n_nb_lambda, n_nb_theta);

      CErm[n_n_max + 1] = ((coeff2 * q_nm.col(2)).matrix() * sin(theta).matrix()).array();
      CEtm[n_n_max + 1] = ((coeff2 * q_nm.col(2)).matrix() * cos(theta).matrix()).array();
      CEfm[n_n_max + 1] = (-I*coeff2 * q_nm.col(2)).replicate(1, n_nb_theta);

      CErm[n_n_max - 1] = ((coeff2 * q_nm.col(0)).matrix() * sin(theta).matrix()).array();
      CEtm[n_n_max - 1] = ((coeff2 * q_nm.col(0)).matrix() * cos(theta).matrix()).array();
      CEfm[n_n_max - 1] = (-I*coeff2 * q_nm.col(0)).replicate(1, n_nb_theta);

      output->CErm = CErm;
      output->CEtm = CEtm;
      output->CEfm = CEfm;

      cout << "r0 = 0 in vshEgenThetaAllPhi" << endl;
      return output;
    }

    ArrayXXr<Real> kr;
    ArrayXr<Real> rho_col;
    stZnAll* st_zn_all_col = new stZnAll();
    if (!isinf(rt(0))) {
      kr = ((2*PI*sqrt(epsilon)/lambda).matrix() * rt.matrix()).array();
      rho_col = kr.transpose().reshaped();
      st_zn_all_col = vshGetZnAll(n_n_max, rho_col, type);
    } else {
      st_zn_all_col->Z0 = st_zn_all_col->Z1 = st_zn_all_col->Z2 =
          ArrayXXc<Real>::Ones(n_nb_lambda*n_nb_theta, n_n_max + 1);
    }

    if (stPT == nullptr)
      stPT = vshPinmTaunm(n_n_max, theta.transpose());

    ArrayXi n_vec, p_vec;
    ArrayXr<Real> pi_nm, tau_nm, d_nm;
    VectorXc<Real> vec_n_dep, vec_n_dep2, mu_n_divd;
    ArrayXXc<Real> Er_sum(n_nb_lambda, n_nb_theta);
    ArrayXXc<Real> Et_sum(n_nb_lambda, n_nb_theta);
    ArrayXXc<Real> Ef_sum(n_nb_lambda, n_nb_theta);
    ArrayXXc<Real> q_nm_for_Z1, ip_nm_for_Z0, q_nm_for_Z2, tmp1, tmp2;
    for (int m = -n_n_max; m <= n_n_max; m++) {
      n_vec = ArrayXi::LinSpaced(n_n_max - abs(m) + 1, abs(m), n_n_max);
      p_vec = n_vec*(n_vec + 1) + m;
      pi_nm = stPT->pi_nm(all, p_vec);
      tau_nm = stPT->tau_nm(all, p_vec);
      d_nm = m ? pi_nm.colwise()*(sin(theta)/m).transpose() : stPT->p_n0;

      if (isinf(rt(0))) {
        q_nm_for_Z1 = ArrayXXc<Real>::Zero(n_nb_lambda, n_nb_theta);
        ip_nm_for_Z0 = p_nm(all, p_vec);
        mu_n_divd = mu_n_divd_gen*pow(-I, n + 1);
      } else {
        q_nm_for_Z1 = q_nm(all, p_vec);
        ip_nm_for_Z0 = q_nm(all, p_vec);
        mu_n_divd = mu_n_divd_gen;
      }
      q_nm_for_Z2 = q_nm(all, p_vec);

      ArrayXi ind_in_rho_col;
      for (int l = 0; l < n_nb_lambda; l++) {
        ind_in_rho_col = ArrayXi::LinSpaced(n_nb_theta, 0, n_nb_theta - 1) + l*n_nb_theta;
        vec_n_dep = (q_nm_for_Z1.row(l) * mu_n_times(n_vec)).transpose().matrix();
        Er_sum.row(l) = ((d_nm*st_zn_all_col->Z1(ind_in_rho_col, n_vec)).matrix() * vec_n_dep).transpose().array();
        tmp1 = mu_n_divd(0, n_vec);
        vec_n_dep = (ip_nm_for_Z0.row(l) * tmp1).transpose().matrix();
        vec_n_dep2 = (q_nm_for_Z2.row(l) * tmp1).transpose().matrix();

        tmp1 = (pi_nm * st_zn_all_col->Z0(ind_in_rho_col, n_vec)).matrix() * vec_n_dep;
        tmp2 = (tau_nm * st_zn_all_col->Z2(ind_in_rho_col, n_vec)).matrix() * vec_n_dep2;
        Et_sum.row(l) = (tmp1 + tmp2).transpose().array();

        tmp1 = (tau_nm * st_zn_all_col->Z0(ind_in_rho_col, n_vec)).matrix() * vec_n_dep;
        tmp2 = (pi_nm * st_zn_all_col->Z2(ind_in_rho_col, n_vec)).matrix() * vec_n_dep2;
        Ef_sum.row(l) = (tmp1 + tmp2).transpose().array();
      }

      output->CErm[m + n_n_max] = pow(-1, m) * Er_sum;
      output->CEtm[m + n_n_max] = pow(-1, m) * Et_sum;
      output->CEfm[m + n_n_max] = pow(-1, m) * Ef_sum;
    }
    delete stPT;
    delete st_zn_all_col;

    return output;
  }

  // Many versions in original code
  template <class Real>
  typename LowLevel<Real>::stZnAll* LowLevel<Real>::vshGetZnAll(size_t n_n_max,
      ArrayXr<Real>& rho, sBessel type) {
    if ((rho == 0).any())
      cout << "Warning: rho = 0 arguments not allowed in vshZnAll..." << endl;

    ArrayXr<Real> nu = ArrayXr<Real>::LinSpaced(n_n_max + 1, 0.5, n_n_max + 0.5);
    ArrayXc<Real> f(rho.size(), n_n_max + 1);

    for (int i = 0; i < rho.size(); i++) {
      f.row(i) = arr_bessel_j(nu, rho(i));
      if ((f.row(i) == 0).any()) {
        cout << "Warning: Bessel (j) calculation went beyond precision in vshGetZnAll()" << endl;
        cout << "x = " << rho(i) << "n_max = " << n_n_max << endl;
      }
    }

    if (type == H1) {
      ArrayXr<Real> y;
      for (int i = 0; i < rho.size(); i++) {
        y = arr_bessel_y(nu, rho(i));
        if ((f.row(i).isInf()).any()) {
          cout << "Warning: Bessel (y) calculation went beyond precision in vshGetZnAll()" << endl;
          cout << "x = " << rho(i) << "n_max = " << n_n_max << endl;
        }
        f.row(i) += I*y;
      }
    }

    f.colwise() *= sqrt((PI/2) / rho);

    RowArrayXc<Real> n = RowArrayXc<Real>::LinSpaced(n_n_max, 1, n_n_max);
    stZnAll* output = new stZnAll();
    output->Z0 = f(all, seq(1, last));
    output->Z1 = (output->Z0).colwise() / rho.template cast<complex<Real>>();
    output->Z2 = f(all, seq(0, last - 1)) - (output->Z1).rowwise() * n;
    return output;
  }

  template <class Real>
  typename LowLevel<Real>::stIncEabnm* LowLevel<Real>::vshGetIncidentCoeffs(
      int n_max, typename HighLevel<Real>::stIncPar* angles) {
    Real alpha_p = angles->alpha_p, phi_p = angles->phi_p;
    Array<Real, 1, 1> theta_p = {{angles->theta_p}};
    int n_p_max = (n_max + 1)*(n_max + 1);

    ArrayXr<Real> n_vec = ArrayXr<Real>::LinSpaced(n_max + 1, 0, n_max);
    // may cause a divide-by-zero exception depending on hardware
    ArrayXc<Real> fact_n = pow(I, n_vec) * sqrt(4*PI*(2*n_vec + 1) / (n_vec*(n_vec+1)));
    fact_n(0) = 0;
    ArrayXr<Real> m_vec = ArrayXr<Real>::LinSpaced(2*n_max + 1, -n_max, n_max);
    ArrayXc<Real> fact_m = -pow(-1, m_vec)*exp(-I*m_vec*phi_p);

    ArrayXi m, ind;
    ArrayXc<Real> d_bar_nm(n_p_max);
    for (int n = 0; n <= n_max; n++) {
      m = ArrayXi::LinSpaced(2*n + 1, -n, n);
      ind = n*(n + 1) + m;
      for (int j = 0; j < 2*n + 1; j++)
        d_bar_nm(ind(j)) = fact_n(n) * fact_m(m(j) + n_max);
    }

    stPinmTaunm* stPTp = vshPinmTaunm(n_max, theta_p);
    ArrayXc<Real> minus_EC_nm_star = cos(alpha_p)*I*stPTp->pi_nm.row(0) + sin(alpha_p)*stPTp->tau_nm.row(0);
    ArrayXc<Real> i_EB_nm_star = I*cos(alpha_p)*stPTp->tau_nm.row(0) + sin(alpha_p)*stPTp->pi_nm.row(0);
    delete stPTp;

    stIncEabnm* output = new stIncEabnm();
    output->a_nm = d_bar_nm * minus_EC_nm_star;
    output->b_nm = d_bar_nm * i_EB_nm_star;
    return output;
  }

  // Many versions in orignial code
  template <class Real>
  typename LowLevel<Real>::stPinmTaunm* LowLevel<Real>::vshPinmTaunm(
      size_t n_max, const ArrayXr<Real>& theta) {
    if ((theta < 0.0).any())
      cout << "Warning: theta must be >= 0 in vshPinmTaunm..." << endl;
    size_t n_rows = size(theta), n_cols, n_p_max = (n_max + 1)*(n_max + 1);
    stPinmTaunm* output = new stPinmTaunm();
    output->pi_nm = ArrayXXr<Real>::Zero(n_rows, n_p_max);
    output->tau_nm = ArrayXXr<Real>::Zero(n_rows, n_p_max);
    ArrayXr<Real> A_m_sin_mm1 = ArrayXr<Real>::Ones(n_rows);
    ArrayXr<Real> mu_c = cos(theta), mu_s = sin(theta), n_vec_real;
    ArrayXi n_vec, p_vec, p_vec_n;

    for (size_t m = 1; m <= n_max; m++) {
      A_m_sin_mm1 *= sqrt(static_cast<Real>(2*m - 1)/(2*m))*
          (m > 1 ? mu_s : ArrayXr<Real>::Ones(n_rows));
      n_cols = n_max - m + 2;
      ArrayXXr<Real> pi_aux(n_rows, n_cols);
      pi_aux.col(0).setZero();
      pi_aux.col(1) = m*A_m_sin_mm1;

      for (size_t j = 2, n = m + 1; j < n_cols; j++, n++)
        pi_aux.col(j) = (1/sqrt((n - m) * (n + m))) * ((2*n - 1)*mu_c*
            pi_aux.col(j - 1) - sqrt((n - 1 - m)*(n - 1 + m)) * pi_aux.col(j - 2));

      n_vec = ArrayXi::LinSpaced(n_cols - 1, m, n_max);
      n_vec_real = n_vec.template cast<Real>();
      p_vec = n_vec*(n_vec + 1) + m;
      p_vec_n = p_vec - 2*m;

      for (size_t n = 0; n < n_cols - 1; n++) {
        (output->pi_nm).col(p_vec(n)) = pi_aux.col(n + 1);
        (output->pi_nm).col(p_vec_n(n)) = pow(-1, (m + 1) % 2) * pi_aux.col(n + 1);

        (output->tau_nm).col(p_vec(n)) = pi_aux.col(n) *
            (-sqrt((n_vec_real(n) - m) * (n_vec_real(n) + m))/m)
            + (mu_c*(n_vec_real(n)/m))*pi_aux.col(n + 1);
        (output->tau_nm).col(p_vec_n(n)) = pow(-1, m % 2)*(output->tau_nm).col(p_vec(n));
      }
    }

    ArrayXXr<Real> p_nm1(n_rows, n_max + 1), t_nm1(n_rows, n_max + 1);
    p_nm1.col(0).setOnes();
    p_nm1.col(1) = mu_c;
    t_nm1.col(0).setZero();
    t_nm1.col(1) = -mu_s;

    for (size_t n = 1; n < n_max; n++) {
      p_nm1.col(n + 1) = static_cast<Real>(2*n + 1)/(n + 1)*mu_c * p_nm1.col(n) -
          static_cast<Real>(n)/(n + 1)*p_nm1.col(n - 1);
      t_nm1.col(n + 1) = mu_c*t_nm1.col(n) - (n+1)*mu_s*p_nm1.col(n);
    }

    output->p_n0 = p_nm1;
    n_vec = ArrayXi::LinSpaced(n_max + 1, 0, n_max);
    p_vec = n_vec*(n_vec + 1);

    for (size_t n = 0; n <= n_max; n++) {
      output->pi_nm.col(p_vec(n)).setZero();
      output->tau_nm.col(p_vec(n)) = t_nm1.col(n);
    }

    return output;
  }

  // Many versions in original code
  template <class Real>
  ArrayXXc<Real>* LowLevel<Real>::vshRBchi(ArrayXr<Real> n, const ArrayXr<Real>& x) {
    ArrayXXc<Real>* chi_x = new ArrayXXc<Real>(x.size(), n.size());
    ArrayXr<Real> yx;
    n += 0.5;
    for (int i = 0; i < x.size(); i++) {
      yx = arr_bessel_y(n, x(i));
      if ((yx.isInf()).any())
        cout << "Warning: Bessel (y) calculation went beyond precision in vshRBchi()" << endl;
      (*chi_x).row(i) = sqrt(static_cast<complex<Real>>(x(i)*PI/2))*yx;
    }
    return chi_x;
  }

  // Many versions in original code
  template <class Real>
  ArrayXXc<Real>* LowLevel<Real>::vshRBpsi(ArrayXr<Real> n, const ArrayXr<Real>& x) {
    ArrayXXc<Real>* psi_x = new ArrayXXc<Real>(x.size(), n.size());
    ArrayXr<Real> jx;
    n += 0.5;
    for (int i = 0; i < x.size(); i++) {
      jx = arr_bessel_j(n, x(i));
      if ((jx == 0.0).any())
        cout << "Warning: Bessel (j) calculation went beyond precision in vshRBpsi()" << endl;
      (*psi_x).row(i) = sqrt(static_cast<complex<Real>>(x(i)*PI/2))*jx;
    }
    return psi_x;
  }

  template class LowLevel<double>;
}
