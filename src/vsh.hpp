#ifndef VSH_HPP
#define VSH_HPP

#include "core.hpp"
#include "math.hpp"

namespace Smarties {

  enum sBessel {J, H1};

  template <class Real>
  struct stPinmTaunm {
    ArrayXXr<Real> pi_nm;
    ArrayXXr<Real> tau_nm;
    ArrayXXr<Real> p_n0;
  };

  template <class Real>
  struct stIncEabnm {
    ArrayXc<Real> a_nm;
    ArrayXc<Real> b_nm;
  };

  template <class Real>
  struct stZnAll {
    ArrayXXc<Real> Z0;
    ArrayXXc<Real> Z1;
    ArrayXXc<Real> Z2;
  };

  template <class Real>
  struct stEAllPhi {
    RowArrayXr<Real> theta;
    RowArrayXr<Real> r_of_theta;
    vector<unique_ptr<ArrayXXc<Real>>> Erm;
    vector<unique_ptr<ArrayXXc<Real>>> Etm;
    vector<unique_ptr<ArrayXXc<Real>>> Efm;
  };

  template <class Real>
  struct stEforPhi {
    ArrayXXc<Real> Er;
    ArrayXXc<Real> Et;
    ArrayXXc<Real> Ef;
    ArrayXr<Real> theta;
    Real phi0;
  };

  template <class Real>
  unique_ptr<stIncPar<Real>> vshMakeIncidentParams(sIncType type, int N_max) {
    auto output = make_unique<stIncPar<Real>>();
    switch (type) {
      case KxEz: {
        output->type = KxEz;
        output->theta_p = mp_pi<Real>()/2;
        output->phi_p = 0;
        output->alpha_p = mp_pi<Real>();
        output->abs_m_vec = ArrayXi::LinSpaced(N_max + 1, 0, N_max);
        break;
      }

      case KxEy: {
        output->type = KxEy;
        output->theta_p = mp_pi<Real>()/2;
        output->phi_p = 0;
        output->alpha_p = mp_pi<Real>()/2;
        output->abs_m_vec = ArrayXi::LinSpaced(N_max + 1, 0, N_max);
        break;
      }

      case KyEz: {
        output->type = KyEz;
        output->theta_p = mp_pi<Real>()/2;
        output->phi_p = mp_pi<Real>()/2;
        output->alpha_p = mp_pi<Real>();
        output->abs_m_vec = ArrayXi::LinSpaced(N_max + 1, 0, N_max);
        break;
      }

      case KyEx: {
        output->type = KyEx;
        output->theta_p = mp_pi<Real>()/2;
        output->phi_p = mp_pi<Real>()/2;
        output->alpha_p = -mp_pi<Real>()/2;
        output->abs_m_vec = ArrayXi::LinSpaced(N_max + 1, 0, N_max);
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
        output->alpha_p = mp_pi<Real>()/2;
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
  unique_ptr<stIncPar<Real>> vshMakeIncidentParams(
      sIncType type, int N_max, Real theta_p, Real phi_p, Real alpha_p) {
    unique_ptr<stIncPar<Real>> output;
    if (type == GENERAL) {
      output = make_unique<stIncPar<Real>>();
      output->type = GENERAL;
      output->theta_p = theta_p;
      output->phi_p = phi_p;
      output->alpha_p = alpha_p;
      output->abs_m_vec = ArrayXi::LinSpaced(N_max + 1, 0, N_max);
    } else {
      output = vshMakeIncidentParams<Real>(type, N_max);
    }
    return output;
  }

  // Many versions in orignial code
  template <class Real>
  unique_ptr<stPinmTaunm<Real>> vshPinmTaunm(int N_max, const ArrayXr<Real>& theta) {
    if ((theta < 0.0).any())
      cerr << "Warning: theta must be >= 0 in vshPinmTaunm..." << endl;
    int n_rows = size(theta), n_cols, P_max = (N_max + 1)*(N_max + 1);
    auto output = make_unique<stPinmTaunm<Real>>();
    output->pi_nm = ArrayXXr<Real>::Zero(n_rows, P_max);
    output->tau_nm = ArrayXXr<Real>::Zero(n_rows, P_max);
    ArrayXr<Real> A_m_sin_mm1 = ArrayXr<Real>::Ones(n_rows);
    ArrayXr<Real> mu_c = cos(theta), mu_s = sin(theta), n_vec_real;
    ArrayXi n_vec, p_vec, p_vec_n;

    for (int m = 1; m <= N_max; m++) {
      A_m_sin_mm1 *= sqrt(static_cast<Real>(2*m - 1)/(2*m))*
          (m > 1 ? mu_s : ArrayXr<Real>::Ones(n_rows));
      n_cols = N_max - m + 2;
      ArrayXXr<Real> pi_aux(n_rows, n_cols);
      pi_aux.col(0).setZero();
      pi_aux.col(1) = m*A_m_sin_mm1;

      for (int j = 2, n = m + 1; j < n_cols; j++, n++)
        pi_aux.col(j) = (1/sqrt((n - m) * (n + m))) * ((2*n - 1)*mu_c*
            pi_aux.col(j - 1) - sqrt((n - 1 - m)*(n - 1 + m)) * pi_aux.col(j - 2));

      n_vec = ArrayXi::LinSpaced(n_cols - 1, m, N_max);
      n_vec_real = n_vec.template cast<Real>();
      p_vec = n_vec*(n_vec + 1) + m;
      p_vec_n = p_vec - 2*m;

      for (int n = 0; n < n_cols - 1; n++) {
        (output->pi_nm).col(p_vec(n)) = pi_aux.col(n + 1);
        (output->pi_nm).col(p_vec_n(n)) = pow(-1, (m + 1) % 2) * pi_aux.col(n + 1);

        (output->tau_nm).col(p_vec(n)) = pi_aux.col(n) *
            (-sqrt((n_vec_real(n) - m) * (n_vec_real(n) + m))/m)
            + (mu_c*(n_vec_real(n)/m))*pi_aux.col(n + 1);
        (output->tau_nm).col(p_vec_n(n)) = pow(-1, m % 2)*(output->tau_nm).col(p_vec(n));
      }
    }

    ArrayXXr<Real> p_nm1(n_rows, N_max + 1), t_nm1(n_rows, N_max + 1);
    p_nm1.col(0).setOnes();
    p_nm1.col(1) = mu_c;
    t_nm1.col(0).setZero();
    t_nm1.col(1) = -mu_s;

    for (int n = 1; n < N_max; n++) {
      p_nm1.col(n + 1) = static_cast<Real>(2*n + 1)/(n + 1)*mu_c * p_nm1.col(n) -
          static_cast<Real>(n)/(n + 1)*p_nm1.col(n - 1);
      t_nm1.col(n + 1) = mu_c*t_nm1.col(n) - (n+1)*mu_s*p_nm1.col(n);
    }

    output->p_n0 = p_nm1;
    n_vec = ArrayXi::LinSpaced(N_max + 1, 0, N_max);
    p_vec = n_vec*(n_vec + 1);

    for (int n = 0; n <= N_max; n++) {
      output->pi_nm.col(p_vec(n)).setZero();
      output->tau_nm.col(p_vec(n)) = t_nm1.col(n);
    }

    return output;
  }

  template <class Real>
  unique_ptr<stIncEabnm<Real>> vshGetIncidentCoeffs(int N_max, const unique_ptr<stIncPar<Real>>& angles) {
    complex<Real> I = mp_im_unit<Real>();
    Real alpha_p = angles->alpha_p, phi_p = angles->phi_p;
    Array<Real, 1, 1> theta_p = {{angles->theta_p}};
    int P_max = (N_max + 1)*(N_max + 1);

    ArrayXr<Real> n_vec = ArrayXr<Real>::LinSpaced(N_max + 1, 0, N_max);
    // may cause a divide-by-zero exception depending on hardware
    ArrayXc<Real> fact_n = pow(I, n_vec) * sqrt(4*mp_pi<Real>()*(2*n_vec + 1) / (n_vec*(n_vec+1)));
    fact_n(0) = 0;
    ArrayXr<Real> m_vec = ArrayXr<Real>::LinSpaced(2*N_max + 1, -N_max, N_max);
    ArrayXc<Real> fact_m = -pow(-1, m_vec)*exp(-I*m_vec*phi_p);

    ArrayXi m, ind;
    ArrayXc<Real> d_bar_nm(P_max);
    for (int n = 0; n <= N_max; n++) {
      m = ArrayXi::LinSpaced(2*n + 1, -n, n);
      ind = n*(n + 1) + m;
      for (int j = 0; j < 2*n + 1; j++)
        d_bar_nm(ind(j)) = fact_n(n) * fact_m(m(j) + N_max);
    }

    auto stPTp = vshPinmTaunm<Real>(N_max, theta_p);
    ArrayXc<Real> minus_EC_nm_star = static_cast<complex<Real>>(cos(alpha_p))*I*stPTp->pi_nm.row(0)
        + sin(alpha_p)*stPTp->tau_nm.row(0);
    ArrayXc<Real> i_EB_nm_star = static_cast<complex<Real>>(cos(alpha_p))*I*stPTp->tau_nm.row(0)
        + sin(alpha_p)*stPTp->pi_nm.row(0);

    auto output = make_unique<stIncEabnm<Real>>();
    output->a_nm = d_bar_nm * minus_EC_nm_star;
    output->b_nm = d_bar_nm * i_EB_nm_star;
    return output;
  }

  // Many versions in original code
  template <class Real>
  unique_ptr<stZnAll<Real>> vshGetZnAll(int N_max, const ArrayXr<Real>& rho, sBessel type) {
    if ((rho == 0).any())
      cerr << "Warning: rho = 0 arguments not allowed in vshZnAll..." << endl;

    ArrayXr<Real> nu = ArrayXr<Real>::LinSpaced(N_max + 1, 0.5, N_max + 0.5);
    ArrayXXc<Real> f(rho.size(), N_max + 1);

    for (int i = 0; i < rho.size(); i++) {
      f.row(i) = arr_bessel_j(nu, rho(i));
      if ((f.row(i) == static_cast<complex<Real>>(0)).any()) {
        cerr << "Warning: Bessel (j) calculation went beyond precision in vshGetZnAll()" << endl;
        cerr << "x = " << rho(i) << "N_max = " << N_max << endl;
      }
    }

    if (type == H1) {
      ArrayXr<Real> y;
      for (int i = 0; i < rho.size(); i++) {
        y = arr_bessel_y(nu, rho(i));
        if ((f.row(i).isInf()).any()) {
          cerr << "Warning: Bessel (y) calculation went beyond precision in vshGetZnAll()" << endl;
          cerr << "x = " << rho(i) << "N_max = " << N_max << endl;
        }
        f.row(i) += mp_im_unit<Real>()*y;
      }
    }

    f.colwise() *= sqrt((mp_pi<Real>()/2) / rho);

    RowArrayXc<Real> n = RowArrayXr<Real>::LinSpaced(N_max, 1, N_max);
    auto output = make_unique<stZnAll<Real>>();
    output->Z0 = f;
    output->Z1 = (output->Z0).colwise() / rho.template cast<complex<Real>>();
    output->Z2 = ArrayXXc<Real>::Zero(rho.size(), N_max + 1);
    output->Z2.rightCols(N_max) = f(all, seq(0, last - 1)) - (output->Z1.rightCols(N_max)).rowwise() * n;
    return output;
  }

  // Many versions in original code
  template <class Real>
  unique_ptr<stEAllPhi<Real>> vshEgenThetaAllPhi(
      const ArrayXr<Real>& lambda, const ArrayXr<Real>& epsilon, const ArrayXXc<Real>& p_nm,
      const ArrayXXc<Real>& q_nm, const RowArrayXr<Real>& rt, const RowArrayXr<Real>& theta,
      sBessel type, unique_ptr<stPinmTaunm<Real>> stPT = unique_ptr<stPinmTaunm<Real>>()) {
    int P_max = p_nm.cols(), N_max = static_cast<int>(round(sqrt(P_max) - 1)),
        Nb_lambda = lambda.size();
    if (rt.size() != theta.size() && rt(0) != 0 && !isinf(rt(0)))
      cout << "vshEgenThetaAllPhi error: theta and rt must be the same size row arrays." << endl;
    int Nb_theta = theta.size();

    ArrayXr<Real> n = ArrayXr<Real>::LinSpaced(N_max + 1, 0, N_max);
    ArrayXr<Real> mu_n_times = sqrt((2*n + 1)*n*(n + 1)/(4*mp_pi<Real>()));
    ArrayXr<Real> mu_n_divd_gen = mu_n_times/(n*(n + 1));
    mu_n_divd_gen(0) = 0;

    auto output = make_unique<stEAllPhi<Real>>();
    output->theta = theta;
    output->r_of_theta = rt;
    output->Erm = vector<unique_ptr<ArrayXXc<Real>>>(2*N_max + 1);
    output->Etm = vector<unique_ptr<ArrayXXc<Real>>>(2*N_max + 1);
    output->Efm = vector<unique_ptr<ArrayXXc<Real>>>(2*N_max + 1);

    if (!rt(0)) {
      for (int m = -N_max; m <= N_max; m++) {
        if (abs(m) > 1) {
          output->Erm[m + N_max] = make_unique<ArrayXXc<Real>>(Nb_lambda, Nb_lambda);
          output->Etm[m + N_max] = make_unique<ArrayXXc<Real>>(Nb_lambda, Nb_lambda);
          output->Efm[m + N_max] = make_unique<ArrayXXc<Real>>(Nb_lambda, Nb_lambda);
          output->Erm[m + N_max]->setZero();
          output->Etm[m + N_max]->setZero();
          output->Efm[m + N_max]->setZero();
        }
      }
      Real coeff1 = 1/sqrt(6*mp_pi<Real>()), coeff2 = coeff1/sqrt(2);

      output->Erm[N_max] = make_unique<ArrayXXc<Real>>(((coeff1 * q_nm.col(2)).matrix() * cos(theta).matrix()).array());
      output->Etm[N_max] = make_unique<ArrayXXc<Real>>(((-coeff1 * q_nm.col(2)).matrix() * sin(theta).matrix()).array());
      output->Efm[N_max] = make_unique<ArrayXXc<Real>>(ArrayXc<Real>::Zero(Nb_lambda, Nb_theta));

      output->Erm[N_max + 1] = make_unique<ArrayXXc<Real>>(((coeff2 * q_nm.col(3)).matrix() * sin(theta).matrix()).array());
      output->Etm[N_max + 1] = make_unique<ArrayXXc<Real>>(((coeff2 * q_nm.col(3)).matrix() * cos(theta).matrix()).array());
      output->Efm[N_max + 1] = make_unique<ArrayXXc<Real>>((-mp_im_unit<Real>()*coeff2 * q_nm.col(3)).replicate(1, Nb_theta));

      output->Erm[N_max - 1] = make_unique<ArrayXXc<Real>>(((coeff2 * q_nm.col(1)).matrix() * sin(theta).matrix()).array());
      output->Etm[N_max - 1] = make_unique<ArrayXXc<Real>>(((coeff2 * q_nm.col(1)).matrix() * cos(theta).matrix()).array());
      output->Efm[N_max - 1] = make_unique<ArrayXXc<Real>>((-mp_im_unit<Real>()*coeff2 * q_nm.col(1)).replicate(1, Nb_theta));

      cout << "r0 = 0 in vshEgenThetaAllPhi" << endl;
      return output;
    }

    ArrayXr<Real> rho_col;
    unique_ptr<stZnAll<Real>> st_zn_all_col;
    if (!isinf(rt(0))) {
      ArrayXXr<Real> kr = ((2*mp_pi<Real>()*sqrt(epsilon)/lambda).matrix() * rt.matrix()).array();
      rho_col = kr.transpose().reshaped();
      st_zn_all_col = vshGetZnAll(N_max, rho_col, type);
    } else {
      st_zn_all_col = make_unique<stZnAll<Real>>();
      st_zn_all_col->Z0 = st_zn_all_col->Z1 = st_zn_all_col->Z2 =
          ArrayXXc<Real>::Ones(Nb_lambda*Nb_theta, N_max + 1);
    }

    if (!stPT)
      stPT = vshPinmTaunm<Real>(N_max, theta.transpose());

    ArrayXi n_vec, p_vec;
    ArrayXXr<Real> pi_nm, tau_nm, d_nm;
    VectorXc<Real> vec_n_dep, vec_n_dep2, mu_n_divd;
    ArrayXXc<Real> Er_sum(Nb_lambda, Nb_theta);
    ArrayXXc<Real> Et_sum(Nb_lambda, Nb_theta);
    ArrayXXc<Real> Ef_sum(Nb_lambda, Nb_theta);
    ArrayXXc<Real> q_nm_for_Z1, ip_nm_for_Z0, q_nm_for_Z2, tmp1, tmp2;
    for (int m = -N_max; m <= N_max; m++) {
      n_vec = ArrayXi::LinSpaced(N_max - abs(m) + 1, abs(m), N_max);
      p_vec = n_vec*(n_vec + 1) + m;
      pi_nm = stPT->pi_nm(all, p_vec);
      tau_nm = stPT->tau_nm(all, p_vec);
      d_nm = m ? pi_nm.colwise()*(sin(theta)/m).transpose() : stPT->p_n0;

      if (isinf(rt(0))) {
        q_nm_for_Z1 = ArrayXXc<Real>::Zero(Nb_lambda, Nb_theta);
        ip_nm_for_Z0 = p_nm(all, p_vec);
        mu_n_divd = mu_n_divd_gen*pow(-mp_im_unit<Real>(), n + 1);
      } else {
        q_nm_for_Z1 = q_nm(all, p_vec);
        ip_nm_for_Z0 = mp_im_unit<Real>()*p_nm(all, p_vec);
        mu_n_divd = mu_n_divd_gen;
      }
      q_nm_for_Z2 = q_nm(all, p_vec);

      ArrayXi ind_in_rho_col;
      for (int l = 0; l < Nb_lambda; l++) {
        ind_in_rho_col = ArrayXi::LinSpaced(Nb_theta, 0, Nb_theta - 1) + l*Nb_theta;
        vec_n_dep = (q_nm_for_Z1.row(l) * mu_n_times(n_vec).transpose()).matrix();
        Er_sum.row(l) = ((d_nm*st_zn_all_col->Z1(ind_in_rho_col, n_vec)).matrix() * vec_n_dep).transpose().array();
        tmp1 = mu_n_divd(n_vec);
        vec_n_dep = (ip_nm_for_Z0.row(l) * tmp1.transpose()).transpose().matrix();
        vec_n_dep2 = (q_nm_for_Z2.row(l) * tmp1.transpose()).transpose().matrix();

        tmp1 = (pi_nm * st_zn_all_col->Z0(ind_in_rho_col, n_vec)).matrix() * vec_n_dep;
        tmp2 = (tau_nm * st_zn_all_col->Z2(ind_in_rho_col, n_vec)).matrix() * vec_n_dep2;
        Et_sum.row(l) = (tmp1 + tmp2).transpose().array();

        tmp1 = (tau_nm * st_zn_all_col->Z0(ind_in_rho_col, n_vec)).matrix() * vec_n_dep;
        tmp2 = (pi_nm * st_zn_all_col->Z2(ind_in_rho_col, n_vec)).matrix() * vec_n_dep2;
        Ef_sum.row(l) = (tmp1 + tmp2).transpose().array();
      }

      output->Erm[m + N_max] = make_unique<ArrayXXc<Real>>(pow(-1, m) * Er_sum);
      output->Etm[m + N_max] = make_unique<ArrayXXc<Real>>(pow(-1, m) * Et_sum);
      output->Efm[m + N_max] = make_unique<ArrayXXc<Real>>(mp_im_unit<Real>() *
          static_cast<complex<Real>>(pow(-1, m)) * Ef_sum);
    }

    return output;
  }

  template <class Real>
  unique_ptr<stEforPhi<Real>> vshEthetaForPhi(const unique_ptr<stEAllPhi<Real>>& stEsurf, Real phi0) {
    auto output = make_unique<stEforPhi<Real>>();
    int N_max = (stEsurf->Erm.size() - 1) / 2;
    output->theta = stEsurf->theta;
    output->phi0 = phi0;

    int Nb_lambda = stEsurf->Erm[0]->rows();
    int Nb_theta = stEsurf->theta.size();

    output->Er = ArrayXXc<Real>::Zero(Nb_lambda, Nb_theta);
    output->Et = ArrayXXc<Real>::Zero(Nb_lambda, Nb_theta);
    output->Ef = ArrayXXc<Real>::Zero(Nb_lambda, Nb_theta);

    complex<Real> exp_phase;
    for (int m = -N_max; m <= N_max; m++) {
      exp_phase = exp(mp_im_unit<Real>()*static_cast<Real>(m)*phi0);
      output->Er += *(stEsurf->Erm[m + N_max]) * exp_phase;
      output->Et += *(stEsurf->Etm[m + N_max]) * exp_phase;
      output->Ef += *(stEsurf->Efm[m + N_max]) * exp_phase;
    }
    return output;
  }

  // Many versions in original code
  template <class Real>
  ArrayXXc<Real> vshRBchi(ArrayXr<Real> n, const ArrayXr<Real>& x) {
    ArrayXXc<Real> chi_x(x.size(), n.size());
    ArrayXr<Real> yx;
    n += 0.5;
    for (int i = 0; i < x.size(); i++) {
      yx = arr_bessel_y(n, x(i));
      if ((yx.isInf()).any())
        cout << "Warning: Bessel (y) calculation went beyond precision in vshRBchi()" << endl;
      chi_x.row(i) = sqrt(static_cast<complex<Real>>(x(i)*mp_pi<Real>()/2))*yx;
    }
    return chi_x;
  }

  // Many versions in original code
  template <class Real>
  ArrayXXc<Real> vshRBpsi(ArrayXr<Real> n, const ArrayXr<Real>& x) {
    ArrayXXc<Real> psi_x(x.size(), n.size());
    ArrayXr<Real> jx;
    n += 0.5;
    for (int i = 0; i < x.size(); i++) {
      jx = arr_bessel_j(n, x(i));
      if ((jx == 0.0).any())
        cout << "Warning: Bessel (j) calculation went beyond precision in vshRBpsi()" << endl;
      psi_x.row(i) = sqrt(static_cast<complex<Real>>(x(i)*mp_pi<Real>()/2))*jx;
    }
    return psi_x;
  }
}

#endif
