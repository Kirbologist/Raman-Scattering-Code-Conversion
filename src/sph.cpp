#include "sph.h"
#include "vsh.h"
#include "math.h"
#include <stdexcept>

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
  stFpovx<Real>* sphGetFpovx(int N_max, Real s, const ArrayXr<Real>& x) {
    int num_x = x.size();

    stFpovx<Real>* output = new stFpovx<Real>();
    output->rb_chi = vshRBchi<Real>(ArrayXr<double>::LinSpaced(N_max + 1, 0, N_max), x);
    output->rb_psi = vshRBpsi<Real>(ArrayXr<double>::LinSpaced(N_max + 1, 0, N_max), s*x);
    output->Fpovx = Tensor3c<Real>(N_max + 1, N_max + 1, num_x);
    output->Fpovx.setZero();

    stFprow<Real>* FpRow = sphGetFpRow<Real>(N_max, s, x);
    ArrayXXc<Real> tmp;
    for (int i = 0; i < N_max - 3; i++) {
      tmp = (FpRow->S.row(i) / x.transpose()).template cast<complex<Real>>();
      output->Fpovx.chip(N_max, 0).chip(i, 0) = TensorCast(tmp, num_x);
    }
    delete FpRow;

    Tensor1c<Real> col_shape(num_x);
    for (int k = N_max; k >= 0; k--) {
      for (int i = k % 2; i <= min(k + 2, N_max); i += 2)
        output->Fpovx.chip(i, 0).chip(k, 0) = TensorCast(output->rb_chi.col(i) * output->rb_psi.col(k), num_x);
      if (k > 0) {
        for (int n = k + 3; n < N_max; n += 2)
          output->Fpovx.chip(n, 0).chip(k - 1, 0) = (output->Fpovx.chip(n + 1, 0).chip(k, 0) +
              output->Fpovx.chip(n - 1, 0).chip(k, 0)) * col_shape.constant(static_cast<Real>(2*k + 1) / (2*n + 1) / s) -
              output->Fpovx.chip(n, 0).chip(k + 1, 0);
      }
    }
    return output;
  }

  // Expects: NB >= N_max
  template <class Real>
  stBessel<Real>* sphGetXiPsi(int N_max, Real s, const ArrayXr<Real>& x, int NB) {
    int num_x = x.size();
    stBessel<Real>* output = new stBessel<Real>();
    stFpovx<Real>* chi_psi_prod = sphGetFpovx<Real>(NB + 1, s, x);
    output->psi_psi = Tensor3c<Real>(N_max + 2, N_max + 2, num_x);
    output->psi_psi.setZero();
    output->psi_n = vshRBpsi<Real>(ArrayXr<Real>::LinSpaced(N_max + 2, 0, N_max + 1), x);

    for (int n = 0; n < N_max + 2; n++) {
      for (int k = n % 2; k < N_max + 2; k += 2)
        output->psi_psi.chip(n, 0).chip(k, 0) =
            TensorCast(chi_psi_prod->rb_psi.col(k) * output->psi_n.col(n), num_x);
    }

    std::array<int, 3> offsets = {0, 0, 0}, extents = {N_max + 2, N_max + 2, num_x};\
    output->xi_psi = output->psi_psi + output->psi_psi.constant(I) *
        chi_psi_prod->Fpovx.slice(offsets, extents);
    output->chi_n = chi_psi_prod->rb_chi.leftCols(N_max + 2);
    output->psi_k = chi_psi_prod->rb_psi.leftCols(N_max + 2);
    delete chi_psi_prod;
    return output;
  }

  // Expects: prods = [(N + 2) x (N + 2) x X] tensor
  template <class Real>
  stBesselPrimes<Real>* sphGetBesselProductsPrimes(Tensor3c<Real>& prods) {
    int N = prods.dimension(0) - 2;
    int X = prods.dimension(2);

    stBesselPrimes<Real>* output = new stBesselPrimes<Real>();
    std::array<int, 3> offsets = {1, 1, 0}, extents = {N, N, X};
    output->xi_psi = prods.slice(offsets, extents);
    output->xi_prime_psi = Tensor3c<Real>(N, N, X);
    output->xi_psi_prime = Tensor3c<Real>(N, N, X);
    output->xi_prime_psi_prime = Tensor3c<Real>(N, N, X);
    output->xi_prime_psi_prime_plus_nnp1_xi_psi_over_ssx = Tensor3c<Real>(N, N, X);
    output->xi_prime_psi_prime_plus_kkp1_xi_psi_over_ssx = Tensor3c<Real>(N, N, X);
    output->xi_psi_over_sxx = Tensor3c<Real>(N, N, X);

    output->xi_prime_psi.setZero();
    output->xi_psi_prime.setZero();
    output->xi_prime_psi_prime.setZero();
    output->xi_prime_psi_prime_plus_nnp1_xi_psi_over_ssx.setZero();
    output->xi_prime_psi_prime_plus_kkp1_xi_psi_over_ssx.setZero();
    output->xi_psi_over_sxx.setZero();

    int k_p1, n_p1;
    Tensor1c<Real> col_shape(X);
    for (int n = 1; n <= N; n++) {
      for (int k = 2 - n % 2; k <= N; k += 2) {
        k_p1 = k*(k + 1);
        n_p1 = n*(n + 1);

        output->xi_prime_psi_prime_plus_kkp1_xi_psi_over_ssx.chip(n - 1, 0).chip(k - 1, 0) = (
            col_shape.constant((k + n + 1)*(k + 1)) * prods.chip(n - 1, 0).chip(k - 1, 0) +
            col_shape.constant(k_p1 - k*(n + 1)) * prods.chip(n - 1, 0).chip(k + 1, 0) +
            col_shape.constant(k_p1 - (k + 1)*n) * prods.chip(n + 1, 0).chip(k - 1, 0) +
            col_shape.constant(k_p1 + k*n) * prods.chip(n + 1, 0).chip(k + 1, 0)) /
            col_shape.constant((2*n + 1)*(2*k + 1));

        output->xi_prime_psi_prime_plus_nnp1_xi_psi_over_ssx.chip(n - 1, 0).chip(k - 1, 0) = (
            col_shape.constant(n_p1 + (k + 1) * (n + 1))*prods.chip(n - 1, 0).chip(k - 1, 0) +
            col_shape.constant((n - k)*(n + 1)) * prods.chip(n - 1, 0).chip(k + 1, 0) +
            col_shape.constant((n - k)*n) * prods.chip(n + 1, 0).chip(k - 1, 0) +
            col_shape.constant(n_p1 + k*n) * prods.chip(n + 1, 0).chip(k + 1, 0)) /
            col_shape.constant((2*n + 1)*(2*k + 1));
      }

      for (int k = 1 + n % 2; k < N; k += 2) {
        output->xi_prime_psi.chip(n - 1, 0).chip(k - 1, 0) =
            (prods.chip(n - 1, 0).chip(k, 0) * col_shape.constant(n + 1) -
            col_shape.constant(n)*prods.chip(n + 1, 0).chip(k, 0)) / col_shape.constant(2*n + 1);
        output->xi_psi_prime.chip(n - 1, 0).chip(k - 1, 0) =
            (prods.chip(n, 0).chip(k - 1, 0) * col_shape.constant(k + 1) -
            col_shape.constant(k)*prods.chip(n, 0).chip(k + 1, 0)) / col_shape.constant(2*k + 1);
      }
    }
    return output;
  }

  // Expects: NB >= N_max
  template <class Real>
  stBesselProducts<Real>* sphGetModifiedBesselProducts(int N_max, Real s, const ArrayXr<Real>& x, int NB) {
    stBesselProducts<Real>* output = new stBesselProducts<Real>();
    stBessel<Real>* prods = sphGetXiPsi<Real>(N_max, s, x, NB);
    output->st_xi_psi_all = sphGetBesselProductsPrimes<Real>(prods->xi_psi);
    output->st_psi_psi_all = sphGetBesselProductsPrimes<Real>(prods->psi_psi);

    ArrayXXc<Real> psi_n_p1_psi_n = (prods->psi_n(all, seq(2, last)) * prods->psi_k(all, seq(1, last - 1))).transpose();
    ArrayXXc<Real> psi_n_psi_n_p1 = (prods->psi_n(all, seq(1, last - 1)) * prods->psi_k(all, seq(2, last))).transpose();
    ArrayXXc<Real> xi_n_p1_psi_n = psi_n_p1_psi_n + I*(prods->chi_n(all, seq(2, last)) * prods->psi_k(all, seq(1, last - 1))).transpose();
    ArrayXXc<Real> xi_n_psi_n_p1 = psi_n_psi_n_p1 + I*(prods->chi_n(all, seq(1, last - 1)) * prods->psi_k(all, seq(2, last))).transpose();
    output->st_psi_psi_all->for_diag_Lt1 = s*psi_n_psi_n_p1 - psi_n_p1_psi_n;
    output->st_xi_psi_all->for_diag_Lt1 = s*xi_n_psi_n_p1 - xi_n_p1_psi_n;
    ArrayXXc<Real> psi_n_psi_n = (prods->psi_n(all, seq(1, last - 1)) * prods->psi_k(all, seq(1, last - 1))).transpose();
    ArrayXXc<Real> xi_n_psi_n = psi_n_psi_n + I*(prods->chi_n(all, seq(1, last - 1)) * prods->psi_k(all, seq(1, last - 1))).transpose();

    delete prods;

    ArrayXXc<Real> n_vec = ArrayXc<Real>::LinSpaced(N_max, 2, N_max + 1);
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
  stRtfunc<Real>* sphMakeGeometry(size_t Nb_theta, Real a, Real c, ArrayXr<Real>* theta) {
    sInt type;
    stRtfunc<Real>* output = new stRtfunc<Real>();
    if (theta) {
      output->w_theta = ArrayXr<Real>::Zero(theta->size());
      output->theta = *theta;
      output->Nb_theta = theta->size();
      type = PTS;
    } else {
      type = GAUSS;
      stRtfunc<Real>* tmp = auxPrepareIntegrals<Real>(2*Nb_theta, type);
      output->theta = tmp->theta(seq(0, Nb_theta - 1));
      output->w_theta = tmp->w_theta(seq(0, Nb_theta - 1))*2;
      output->Nb_theta = Nb_theta;
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

  template <class Real>
  size_t sphCheckBesselConvergence(size_t N_req, Real s, ArrayXr<Real> x, Real acc, size_t N_min) {
    size_t NB_start = max(N_req, N_min);
    size_t NB = NB_start;
    stFpovx<Real>* prod = sphGetFpovx<Real>(NB, s, x);
    bool to_continue = true;
    size_t max_N = 500, NB_step = 16;

    int NB_next;
    stFpovx<Real>* prod_new;
    Tensor<Real, 0> rel_acc_ee, rel_acc_oo;
    Real rel_acc;
    ArithmeticSequence<long int, long int, long int>
        seq1(0, N_req, 2), seq2(1, N_req, 2), seq3(0, x.size() - 1); // prod->Fpovx has no. of columns = x.size()
    long int dim1 = seq1.size(), dim3 = seq3.size();
    Tensor3c<Real> tmp(dim1, dim1, dim3);
    while (to_continue && NB < max_N) {
      NB_next = NB + NB_step;
      prod_new = sphGetFpovx<Real>(NB_next, s, x);
      tmp = subtensor<Real>(prod->Fpovx, seq1, seq1, seq3) /
          subtensor<Real>(prod_new->Fpovx, seq1, seq1, seq3) - tmp.constant(1);
      rel_acc_ee = tmp.abs().real().maximum();
      tmp = subtensor<Real>(prod->Fpovx, seq2, seq2, seq3) /
          subtensor<Real>(prod_new->Fpovx, seq2, seq2, seq3) - tmp.constant(1);
      rel_acc_oo = tmp.abs().real().maximum();
      rel_acc = max(rel_acc_ee(0), rel_acc_oo(0));
      if (rel_acc < acc)
        to_continue = false;
      else {
        NB = NB_next;
        delete prod;
        prod = prod_new;
      }
    }
    delete prod_new;

    if (NB > NB_start) {
      NB -= NB_step;
      NB_step = 1;
      to_continue = true;
      while (to_continue && NB < max_N) {
        NB += NB_step;
        delete prod;
        prod = sphGetFpovx<Real>(NB, s, x);
        tmp = subtensor<Real>(prod->Fpovx, seq1, seq1, seq3) /
            subtensor<Real>(prod_new->Fpovx, seq1, seq1, seq3) - tmp.constant(1);
        rel_acc_ee = tmp.abs().maximum().real();
        tmp = subtensor<Real>(prod->Fpovx, seq2, seq2, seq3) /
            subtensor<Real>(prod_new->Fpovx, seq2, seq2, seq3) - tmp.constant(1);
        rel_acc_oo = tmp.abs().maximum().real();
        rel_acc = max(rel_acc_ee(0), rel_acc_oo(0));
        if (rel_acc < acc)
          to_continue = false;
      }
    }

    delete prod;
    if (to_continue)
      throw(runtime_error("Problem in sphEstimateNB: convergence was not achieved"));
    return NB;
  }

  template <class Real>
  size_t sphEstimateNB(size_t NQ, stRtfunc<Real>* stGeometry, stParams<Real>* params, Real acc) {
    ArrayXr<Real> s = params->s;
    ArrayXr<Real> k1 = params->k1;

    ArrayXr<Real> x_max = {{k1.maxCoeff() * stGeometry->r.maxCoeff()}};
    typename ArrayXr<Real>::Index ind;
    abs(s).minCoeff(&ind);
    Real s_min = s(ind);
    abs(s).maxCoeff(&ind);
    Real s_max = s(ind);

    size_t N1 = sphCheckBesselConvergence(NQ, s_max, x_max, acc, NQ);
    size_t NB = sphCheckBesselConvergence(NQ, s_min, x_max, acc, N1);
    return NB;
  }

  /*
  template <class Real>
  st4M<Real>* sphCalculatePQ(int N_max, ArrayXi abs_m_vec,
      stRtfunc<Real>* Rt_func, stParams<Real>* params, int NB) {
    if (NB < N_max)
      NB = N_max;
    size_t M = abs_m_vec.size();
    bool print_output = params->output;
    if (print_output)
      cout << "sphCalculatePQ: Calculate P, Q for " << M << "m-values with N_Q = " <<
          N_max << ", NB = " << NB << ", N_Theta = " << Rt_func->Nb_theta << endl;

    stPQa<Real>* output = new stPQa[M]();
    ArrayXr<Real> k1 = params->k1, s = params->s;
    int N_int = Rt_func->Nb_theta, T = N_int;
    ArrayXr<Real> x = Rt_func->r * k1, x_theta = Rt_func->dr_dt * k1

    stPinmTaunm<Real>* stPT = vshPinmTaunm(N_max, Rt_func->theta);
    ArrayXr<Real> sin_t = sin(Rt_func->theta);
    ArrayXr<Real> dx_dt_wt = x_theta * Rt_func->w_theta;
    ArrayXr<Real> n_vec = ArrayXr<Real>::LinSpaced(N_max, 1, N_max);
    ArrayXr<Real> An_vec = sqrt((2*n_vec + 1) / (2*n_vec * (n_vec + 1)));
    ArrayXXr<Real> An_Ak = (An_vec.matrix() * An_vec.transpose().matrix()).array();

    stBesselProducts<Real>* prods = sphGetModifiedBesselProducts(N_max, s, x, NB);
    ArrayXXc<Real>* xi_prime_psi = new ArrayXXc<Real>[N_max]();
    ArrayXXc<Real>* xi_psi_prime = new ArrayXXc<Real>[N_max]();
    ArrayXXc<Real>* xi_psi = new ArrayXXc<Real>[N_max]();
    ArrayXXc<Real>* xi_psi = new ArrayXXc<Real>[N_max]();
  }
  */

  template ArrayXXr<double>* sphGetUforFp(int);
  template stFprow<double>* sphGetFpRow(int, double, const ArrayXr<double>&);
  template stFpovx<double>* sphGetFpovx(int, double, const ArrayXr<double>&);
  template stBessel<double>* sphGetXiPsi(int, double, const ArrayXr<double>&, int);
  template stBesselPrimes<double>* sphGetBesselProductsPrimes(Tensor3c<double>&);
  template stBesselProducts<double>* sphGetModifiedBesselProducts(int, double, const ArrayXr<double>&, int);
  template stRtfunc<double>* sphMakeGeometry(size_t, double, double, ArrayXr<double>*);
  template size_t sphCheckBesselConvergence(size_t, double, ArrayXr<double>, double, size_t);
  template size_t sphEstimateNB(size_t, stRtfunc<double>*, stParams<double>*, double);
}
