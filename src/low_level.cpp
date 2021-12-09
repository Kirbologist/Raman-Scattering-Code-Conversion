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
    size_t nNint, scheme scheme_type) {
    stRtfunc* output = new stRtfunc();
    output->nNbTheta = nNint;

    switch (scheme_type) {
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
  typename LowLevel<Real>::stIncEabnm* LowLevel<Real>::vshGetIncidentCoeffs(
      int n_max, typename HighLevel<Real>::stIncPar& angles) {
    Real alpha_p = angles.alpha_p, phi_p = angles.phi_p;
    Array<Real, 1, 1> theta_p = {{angles.theta_p}};
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

    cout << "fact n" << fact_n << endl << "fact_m" << fact_m << endl << "d_bar_nm" << d_bar_nm << endl;

    stPinmTaunm* stPTp = vshPinmTaunm(n_max, theta_p);
    ArrayXc<Real> minus_EC_nm_star = cos(alpha_p)*I*stPTp->pi_nm.row(0) + sin(alpha_p)*stPTp->tau_nm.row(0);
    ArrayXc<Real> i_EB_nm_star = I*cos(alpha_p)*stPTp->tau_nm.row(0) + sin(alpha_p)*stPTp->pi_nm.row(0);

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
  ArrayXXc<Real>* LowLevel<Real>::vshRBchi(RowArrayXr<Real> n, const ArrayXr<Real>& x) {
    ArrayXXc<Real>* chi_x = new ArrayXXc<Real>(x.size(), n.size());
    RowArrayXr<Real> yx;
    n += 0.5;
    for (int i = 0; i < size(x); i++) {
      yx = arr_bessel_y(n, x(i));
      if ((yx.isInf()).any())
        cout << "Warning: Bessel (y) calculation went beyond precision" << endl;
      (*chi_x).row(i) = sqrt(static_cast<complex<Real>>(x(i)*PI/2))*yx;
    }
    return chi_x;
  }

  // Many versions in original code
  template <class Real>
  ArrayXXc<Real>* LowLevel<Real>::vshRBpsi(RowArrayXr<Real> n, const ArrayXr<Real>& x) {
    ArrayXXc<Real>* psi_x = new ArrayXXc<Real>(x.size(), n.size());
    RowArrayXr<Real> jx;
    n += 0.5;
    for (int i = 0; i < size(x); i++) {
      jx = arr_bessel_j(n, x(i));
      if ((jx == 0.0).any())
        cout << "Warning: Bessel (j) calculation went beyond precision" << endl;
      (*psi_x).row(i) = sqrt(static_cast<complex<Real>>(x(i)*PI/2))*jx;
    }
    return psi_x;
  }

  template class LowLevel<double>;
}
