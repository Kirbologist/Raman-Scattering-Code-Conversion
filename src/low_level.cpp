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
    ArrayXr<Real> xu(N1, 1), y(N1, 1), y0(N1, 1), Lp(N1, 1);
    ArrayXXr<Real> L(N1, N2);

    xu = ArrayXr<Real>::LinSpaced(N1, -1, 1);
    y = cos((2*ArrayXr<Real>::LinSpaced(N1, 0, N) + 1)*PI/(2*N + 2)) +
        (0.27/N1)*sin(PI*xu*N/N2);
    y0 = 2*ArrayXr<Real>::Ones(N1);
    L.col(0) = ArrayXr<Real>::Ones(N1);

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
    ArrayXXc<Real> CErm[2*n_n_max + 1], CEtm[2*n_n_max + 1], CEfm[2*n_n_max + 1];

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
    ArrayXXc<Real> Er_sum, Et_sum, Ef_sum, q_nm_for_Z1, ip_nm_for_Z0, q_nm_for_Z2,
        tmp1, tmp2;
    for (int m = -n_n_max; m <= n_n_max; m++) {
      n_vec = ArrayXi::LinSpaced(n_n_max - abs(m) + 1, abs(m), n_n_max);
      p_vec = n_vec*(n_vec + 1) + m;
      pi_nm = stPT->pi_nm(all, p_vec);
      tau_nm = stPT->tau_nm(all, p_vec);
      d_nm = m ? pi_nm*(sin(theta)/m).transpose().replicate(1, n_n_max - abs(m) + 1) : stPT->p_n0;

      Er_sum = Et_sum = Ef_sum = ArrayXXc<Real>::Zero(n_nb_lambda, n_nb_theta);

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

    return output;
  }

  template <class Real>
  typename LowLevel<Real>::stZnAll* LowLevel<Real>::vshGetZnAll(size_t n_n_max,
      ArrayXr<Real>& rho, sBessel type) {
    if ((rho == 0).any())
      cout << "Warning: rho = 0 arguments not allowed in vshZnAll..." << endl;

    //RowArrayXi wasn't supported for some reason
    RowArrayXr<Real> n = RowArrayXr<Real>::LinSpaced(n_n_max, 1, n_n_max);
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

    f *= sqrt((PI/2) / rho.replicate(1, n_n_max + 1));

    stZnAll* output = new stZnAll();
    output->Z0 = f(all, seq(1, last));
    output->Z1 = output->Z0 * 1/rho.replicate(1, n_n_max);
    output->Z2 = f(all, seq(0, last - 1)) - output->Z1 * n.replicate(rho.size(), 1);
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
      pi_aux.col(0) = ArrayXr<Real>::Zero(n_rows);
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
    p_nm1.col(0) = ArrayXr<Real>::Ones(n_rows);
    p_nm1.col(1) = mu_c;
    t_nm1.col(0) = ArrayXr<Real>::Zero(n_rows);
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
      output->pi_nm.col(p_vec(n)) = ArrayXr<Real>::Zero(n_rows);
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
