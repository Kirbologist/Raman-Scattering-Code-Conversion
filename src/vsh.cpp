#include "vsh.h"
#include "math.h"

namespace Raman {
  template <class Real>
  stIncPar<Real>* vshMakeIncidentParams(sIncType type, size_t n_max) {
    stIncPar<Real>* output = new stIncPar<Real>();
    switch (type) {
      case KxEz: {
        output->type = KxEz;
        output->theta_p = PI/2;
        output->phi_p = 0;
        output->alpha_p = PI;
        output->abs_m_vec = ArrayXi::LinSpaced(n_max + 1, 0, n_max);
        break;
      }

      case KxEy: {
        output->type = KxEy;
        output->theta_p = PI/2;
        output->phi_p = 0;
        output->alpha_p = PI/2;
        output->abs_m_vec = ArrayXi::LinSpaced(n_max + 1, 0, n_max);
        break;
      }

      case KyEz: {
        output->type = KyEz;
        output->theta_p = PI/2;
        output->phi_p = PI/2;
        output->alpha_p = PI;
        output->abs_m_vec = ArrayXi::LinSpaced(n_max + 1, 0, n_max);
        break;
      }

      case KyEx: {
        output->type = KyEx;
        output->theta_p = PI/2;
        output->phi_p = PI/2;
        output->alpha_p = -PI/2;
        output->abs_m_vec = ArrayXi::LinSpaced(n_max + 1, 0, n_max);
        break;
      }

      case KzEx: {
        output->type = KzEx;
        output->theta_p = 0;
        output->phi_p = 0;
        output->alpha_p = 0;
        output->abs_m_vec = {1};
        break;
      }

      case KzEy: {
        output->type = KzEy;
        output->theta_p = 0;
        output->phi_p = 0;
        output->alpha_p = PI/2;
        output->abs_m_vec = {1};
        break;
      }

      default:
        output->type = KzEx;
        output->theta_p = 0;
        output->phi_p = 0;
        output->alpha_p = 0;
        output->abs_m_vec = {1};
    }
    return output;
  }

  template <class Real>
  stIncPar<Real>* vshMakeIncidentParams(
      sIncType type, size_t n_max, Real theta_p, Real phi_p, Real alpha_p) {
    stIncPar<Real>* output = new stIncPar<Real>();
    if (type == GENERAL) {
      output->type = GENERAL;
      output->theta_p = theta_p;
      output->phi_p = phi_p;
      output->alpha_p = alpha_p;
      output->abs_m_vec = ArrayXi::LinSpaced(n_max + 1, 0, n_max);
    } else {
      output = vshMakeIncidentParams<Real>(type, n_max);
    }
    return output;
  }

  // Many versions in orignial code
  template <class Real>
  stPinmTaunm<Real>* vshPinmTaunm(size_t n_max, const ArrayXr<Real>& theta) {
    if ((theta < 0.0).any())
      cout << "Warning: theta must be >= 0 in vshPinmTaunm..." << endl;
    size_t n_rows = size(theta), n_cols, n_p_max = (n_max + 1)*(n_max + 1);
    stPinmTaunm<Real>* output = new stPinmTaunm<Real>();
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

  template <class Real>
  stIncEabnm<Real>* vshGetIncidentCoeffs(
      int n_max, stIncPar<Real>* angles) {
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

    stPinmTaunm<Real>* stPTp = vshPinmTaunm<Real>(n_max, theta_p);
    ArrayXc<Real> minus_EC_nm_star = cos(alpha_p)*I*stPTp->pi_nm.row(0) + sin(alpha_p)*stPTp->tau_nm.row(0);
    ArrayXc<Real> i_EB_nm_star = I*cos(alpha_p)*stPTp->tau_nm.row(0) + sin(alpha_p)*stPTp->pi_nm.row(0);
    delete stPTp;

    stIncEabnm<Real>* output = new stIncEabnm<Real>();
    output->a_nm = d_bar_nm * minus_EC_nm_star;
    output->b_nm = d_bar_nm * i_EB_nm_star;
    return output;
  }

  // Many versions in original code
  template <class Real>
  stZnAll<Real>* vshGetZnAll(size_t n_n_max, const ArrayXr<Real>& rho, sBessel type) {
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
    stZnAll<Real>* output = new stZnAll<Real>();
    output->Z0 = f(all, seq(1, last));
    output->Z1 = (output->Z0).colwise() / rho.template cast<complex<Real>>();
    output->Z2 = f(all, seq(0, last - 1)) - (output->Z1).rowwise() * n;
    return output;
  }

  // Many versions in original code
  template <class Real>
  stEAllPhi<Real>* vshEgenThetaAllPhi(
      const ArrayXr<Real>& lambda, const ArrayXr<Real>& epsilon, const ArrayXXc<Real>& p_nm,
      const ArrayXXc<Real>& q_nm, const RowArrayXr<Real>& rt, const RowArrayXr<Real>& theta,
      sBessel type, stPinmTaunm<Real>* stPT) {
    int n_p_max = p_nm.cols(), n_n_max = static_cast<int>(round(sqrt(n_p_max) - 1)),
        n_nb_lambda = lambda.size();
    if (rt.size() != theta.size() && rt(0) != 0 && !isinf(rt(0)))
      cout << "vshEgenThetaAllPhi error: theta and rt must be the same size row arrays." << endl;
    int n_Nb_theta = theta.size();

    ArrayXr<Real> n = ArrayXr<Real>::LinSpaced(n_n_max + 1, 0, n_n_max);
    ArrayXr<Real> mu_n_times = sqrt((2*n + 1)*n*(n + 1)/(4*PI));
    ArrayXr<Real> mu_n_divd_gen = mu_n_times/(n*(n + 1));

    stEAllPhi<Real>* output = new stEAllPhi<Real>();
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
      CEfm[n_n_max] = ArrayXc<Real>::Zero(n_nb_lambda, n_Nb_theta);

      CErm[n_n_max + 1] = ((coeff2 * q_nm.col(2)).matrix() * sin(theta).matrix()).array();
      CEtm[n_n_max + 1] = ((coeff2 * q_nm.col(2)).matrix() * cos(theta).matrix()).array();
      CEfm[n_n_max + 1] = (-I*coeff2 * q_nm.col(2)).replicate(1, n_Nb_theta);

      CErm[n_n_max - 1] = ((coeff2 * q_nm.col(0)).matrix() * sin(theta).matrix()).array();
      CEtm[n_n_max - 1] = ((coeff2 * q_nm.col(0)).matrix() * cos(theta).matrix()).array();
      CEfm[n_n_max - 1] = (-I*coeff2 * q_nm.col(0)).replicate(1, n_Nb_theta);

      output->CErm = CErm;
      output->CEtm = CEtm;
      output->CEfm = CEfm;

      cout << "r0 = 0 in vshEgenThetaAllPhi" << endl;
      return output;
    }

    ArrayXXr<Real> kr;
    ArrayXr<Real> rho_col;
    stZnAll<Real>* st_zn_all_col = new stZnAll<Real>();
    if (!isinf(rt(0))) {
      kr = ((2*PI*sqrt(epsilon)/lambda).matrix() * rt.matrix()).array();
      rho_col = kr.transpose().reshaped();
      st_zn_all_col = vshGetZnAll(n_n_max, rho_col, type);
    } else {
      st_zn_all_col->Z0 = st_zn_all_col->Z1 = st_zn_all_col->Z2 =
          ArrayXXc<Real>::Ones(n_nb_lambda*n_Nb_theta, n_n_max + 1);
    }

    if (stPT == nullptr)
      stPT = vshPinmTaunm<Real>(n_n_max, theta.transpose());

    ArrayXi n_vec, p_vec;
    ArrayXr<Real> pi_nm, tau_nm, d_nm;
    VectorXc<Real> vec_n_dep, vec_n_dep2, mu_n_divd;
    ArrayXXc<Real> Er_sum(n_nb_lambda, n_Nb_theta);
    ArrayXXc<Real> Et_sum(n_nb_lambda, n_Nb_theta);
    ArrayXXc<Real> Ef_sum(n_nb_lambda, n_Nb_theta);
    ArrayXXc<Real> q_nm_for_Z1, ip_nm_for_Z0, q_nm_for_Z2, tmp1, tmp2;
    for (int m = -n_n_max; m <= n_n_max; m++) {
      n_vec = ArrayXi::LinSpaced(n_n_max - abs(m) + 1, abs(m), n_n_max);
      p_vec = n_vec*(n_vec + 1) + m;
      pi_nm = stPT->pi_nm(all, p_vec);
      tau_nm = stPT->tau_nm(all, p_vec);
      d_nm = m ? pi_nm.colwise()*(sin(theta)/m).transpose() : stPT->p_n0;

      if (isinf(rt(0))) {
        q_nm_for_Z1 = ArrayXXc<Real>::Zero(n_nb_lambda, n_Nb_theta);
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
        ind_in_rho_col = ArrayXi::LinSpaced(n_Nb_theta, 0, n_Nb_theta - 1) + l*n_Nb_theta;
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
  ArrayXXc<Real>& vshRBchi(ArrayXr<Real> n, const ArrayXr<Real>& x) {
    ArrayXXc<Real>* chi_x = new ArrayXXc<Real>(x.size(), n.size());
    ArrayXr<Real> yx;
    n += 0.5;
    for (int i = 0; i < x.size(); i++) {
      yx = arr_bessel_y(n, x(i));
      if ((yx.isInf()).any())
        cout << "Warning: Bessel (y) calculation went beyond precision in vshRBchi()" << endl;
      (*chi_x).row(i) = sqrt(static_cast<complex<Real>>(x(i)*PI/2))*yx;
    }
    return *chi_x;
  }

  // Many versions in original code
  template <class Real>
  ArrayXXc<Real>& vshRBpsi(ArrayXr<Real> n, const ArrayXr<Real>& x) {
    ArrayXXc<Real>* psi_x = new ArrayXXc<Real>(x.size(), n.size());
    ArrayXr<Real> jx;
    n += 0.5;
    for (int i = 0; i < x.size(); i++) {
      jx = arr_bessel_j(n, x(i));
      if ((jx == 0.0).any())
        cout << "Warning: Bessel (j) calculation went beyond precision in vshRBpsi()" << endl;
      (*psi_x).row(i) = sqrt(static_cast<complex<Real>>(x(i)*PI/2))*jx;
    }
    return *psi_x;
  }

  template stIncPar<double>* vshMakeIncidentParams(sIncType, size_t);
  template stIncPar<double>* vshMakeIncidentParams(sIncType, size_t, double, double, double);
  template stPinmTaunm<double>* vshPinmTaunm(size_t, const ArrayXr<double>&);
  template stIncEabnm<double>* vshGetIncidentCoeffs(int, stIncPar<double>*);
  template stZnAll<double>* vshGetZnAll(size_t, const ArrayXr<double>&, sBessel);
  template stEAllPhi<double>* vshEgenThetaAllPhi(const ArrayXr<double>&,
      const ArrayXr<double>&, const ArrayXXc<double>&, const ArrayXXc<double>&,
      const RowArrayXr<double>&, const RowArrayXr<double>&, sBessel, stPinmTaunm<double>*);
  template ArrayXXc<double>& vshRBchi(ArrayXr<double>, const ArrayXr<double>&);
  template ArrayXXc<double>& vshRBpsi(ArrayXr<double>, const ArrayXr<double>&);
}