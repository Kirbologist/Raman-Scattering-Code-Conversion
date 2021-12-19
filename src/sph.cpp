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
  stBesselPrimes<Real>* sphGetBesselProductsPrimes(const Tensor3c<Real>& prods) {
    int N = prods.dimension(0) - 2;
    int X = prods.dimension(2);

    stBesselPrimes<Real>* output = new stBesselPrimes<Real>();
    std::array<int, 3> offsets = {0, 0, 0}, extents = {N + 1, N + 1, X};
    output->xi_psi = prods.slice(offsets, extents);
    output->xi_prime_psi = Tensor3c<Real>(N + 1, N + 1, X);
    output->xi_psi_prime = Tensor3c<Real>(N + 1, N + 1, X);
    output->xi_prime_psi_prime = Tensor3c<Real>(N + 1, N + 1, X);
    output->xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx = Tensor3c<Real>(N + 1, N + 1, X);
    output->xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx = Tensor3c<Real>(N + 1, N + 1, X);
    output->xi_psi_over_sxx = Tensor3c<Real>(N + 1, N + 1, X);

    output->xi_prime_psi.setZero();
    output->xi_psi_prime.setZero();
    output->xi_prime_psi_prime.setZero();
    output->xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx.setZero();
    output->xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx.setZero();
    output->xi_psi_over_sxx.setZero();

    int kkp1, nnp1;
    Tensor1c<Real> col_shape(X);
    for (int n = 0; n <= N; n++) {
      for (int k = n % 2; k <= N; k += 2) {
        kkp1 = k*(k + 1);
        nnp1 = n*(n + 1);

        output->xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx.chip(n, 0).chip(k, 0) = (
            col_shape.constant((k + n + 1)*(k + 1)) * prods.chip(n - 1, 0).chip(k - 1, 0) +
            col_shape.constant(kkp1 - k*(n + 1)) * prods.chip(n - 1, 0).chip(k + 1, 0) +
            col_shape.constant(kkp1 - (k + 1)*n) * prods.chip(n + 1, 0).chip(k - 1, 0) +
            col_shape.constant(kkp1 + k*n) * prods.chip(n + 1, 0).chip(k + 1, 0)) /
            col_shape.constant((2*n + 1)*(2*k + 1));

        output->xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx.chip(n, 0).chip(k, 0) = (
            col_shape.constant(nnp1 + (k + 1) * (n + 1))*prods.chip(n - 1, 0).chip(k - 1, 0) +
            col_shape.constant((n - k)*(n + 1)) * prods.chip(n - 1, 0).chip(k + 1, 0) +
            col_shape.constant((n - k)*n) * prods.chip(n + 1, 0).chip(k - 1, 0) +
            col_shape.constant(nnp1 + k*n) * prods.chip(n + 1, 0).chip(k + 1, 0)) /
            col_shape.constant((2*n + 1)*(2*k + 1));
      }

      for (int k = 1 - n % 2; k <= N; k += 2) {
        output->xi_prime_psi.chip(n, 0).chip(k, 0) =
            (prods.chip(n - 1, 0).chip(k, 0) * col_shape.constant(n + 1) -
            col_shape.constant(n)*prods.chip(n + 1, 0).chip(k, 0)) / col_shape.constant(2*n + 1);
        output->xi_psi_prime.chip(n, 0).chip(k, 0) =
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

    ArrayXXc<Real> psi_np1_psi_n = (prods->psi_n(all, seq(1, last)) * prods->psi_k(all, seq(0, last - 1))).transpose();
    ArrayXXc<Real> psi_n_psi_np1 = (prods->psi_n(all, seq(0, last - 1)) * prods->psi_k(all, seq(1, last))).transpose();
    ArrayXXc<Real> xi_np1_psi_n = psi_np1_psi_n + I*(prods->chi_n(all, seq(1, last)) * prods->psi_k(all, seq(0, last - 1))).transpose();
    ArrayXXc<Real> xi_n_psi_np1 = psi_n_psi_np1 + I*(prods->chi_n(all, seq(0, last - 1)) * prods->psi_k(all, seq(1, last))).transpose();
    output->st_psi_psi_all->for_diag_Lt1 = s*psi_n_psi_np1 - psi_np1_psi_n;
    output->st_xi_psi_all->for_diag_Lt1 = s*xi_n_psi_np1 - xi_np1_psi_n;
    ArrayXXc<Real> psi_n_psi_n = (prods->psi_n(all, seq(0, last - 1)) * prods->psi_k(all, seq(0, last - 1))).transpose();
    ArrayXXc<Real> xi_n_psi_n = psi_n_psi_n + I*(prods->chi_n(all, seq(0, last - 1)) * prods->psi_k(all, seq(0, last - 1))).transpose();

    delete prods;

    ArrayXXc<Real> n_vec = ArrayXc<Real>::LinSpaced(N_max + 1, 1, N_max + 1);
    output->st_psi_psi_all->for_diag_Lt2 = psi_n_psi_np1 - s*psi_np1_psi_n + (s - 1)*(s + 1)/s *
        (n_vec.matrix() * (1/x.transpose()).matrix()).array() * psi_n_psi_n;
    output->st_xi_psi_all->for_diag_Lt2 = xi_n_psi_np1 - s*xi_np1_psi_n + (s - 1)*(s + 1)/s *
        (n_vec.matrix() * (1/x.transpose()).matrix()).array() *xi_n_psi_n;
    RowArrayXc<Real> tmp = s*x.pow(2).transpose();
    output->st_psi_psi_all->for_diag_Lt3 = psi_n_psi_n.rowwise() / tmp;
    output->st_xi_psi_all->for_diag_Lt3 = xi_n_psi_n.rowwise() / tmp;

    return output;
  }

  template <class Real>
  stRtfunc<Real>* sphMakeGeometry(size_t Nb_theta, Real a, Real c, const ArrayXr<Real>* theta) {
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
  size_t sphCheckBesselConvergence(size_t N_req, Real s, const ArrayXr<Real>& x, Real acc, size_t N_min) {
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
  size_t sphEstimateNB(size_t NQ, const stRtfunc<Real>* stGeometry, const stParams<Real>* params, Real acc) {
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

  template <class Real>
  stPQa<Real>* sphCalculatePQ(int N_max, const ArrayXi& abs_m_vec,
      const stRtfunc<Real>* Rt_func, const stParams<Real>* params, int NB) {
    if (params->s.size() > 1 || params->k1.size() > 1)
      throw(runtime_error("params->s and params->k1 must be scalar"));
    if (NB < N_max)
      NB = N_max;
    int M = abs_m_vec.size();
    //if (params->output) (params shouldn't have an output member, perhaps use stOptions instead)
    cout << "sphCalculatePQ: Calculate P, Q for " << M << "m-values with N_Q = " <<
        N_max << ", NB = " << NB << ", N_Theta = " << Rt_func->Nb_theta << endl;

    stPQa<Real>* output = new stPQa<Real>[M]();
    Real s = params->s(0);
    Real k1 = params->k1(0);
    int N_int = Rt_func->Nb_theta, T = N_int; // Number of theta's
    ArrayXr<Real> x = Rt_func->r * k1; // [T x 1] x(theta)
    ArrayXr<Real> x_theta = Rt_func->dr_dt * k1; // [T x 1] x'(theta)

    stPinmTaunm<Real>* stPT = vshPinmTaunm(N_max, Rt_func->theta);
    ArrayXr<Real> sin_t = sin(Rt_func->theta); // [T x 1]
    ArrayXr<Real> dx_dt_wt = x_theta * Rt_func->w_theta; // [T x 1]
    ArrayXr<Real> tmp1 = ArrayXr<Real>::LinSpaced(N_max + 1, 0, N_max); // [N x 1]
    ArrayXr<Real> An_vec = sqrt((2*tmp1 + 1) / (2*tmp1 * (tmp1 + 1))); // [N x 1]
    ArrayXXr<Real> An_Ak = (An_vec.matrix() * An_vec.transpose().matrix()).array(); // [N x N]

    stBesselProducts<Real>* prods = sphGetModifiedBesselProducts(N_max, s, x, NB); // prods contains [N x N x T] tensors
    // Converting these [N x N x T] tensors into [N x T x N] tensors
    Tensor3c<Real> xi_prime_psi(N_max + 1, T, N_max + 1);
    Tensor3c<Real> xi_psi_prime(N_max + 1, T, N_max + 1);
    Tensor3c<Real> xi_psi(N_max + 1, T, N_max + 1);
    Tensor3c<Real> xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx(N_max + 1, T, N_max + 1);
    Tensor3c<Real> xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx(N_max + 1, T, N_max + 1);
    for (int k = 0; k <= N_max; k++) {
      xi_psi.chip(k, 2) = prods->st_xi_psi_all->xi_psi.chip(k, 1);
      xi_prime_psi.chip(k, 2) = prods->st_xi_psi_all->xi_prime_psi.chip(k, 1);
      xi_psi_prime.chip(k, 2) = prods->st_xi_psi_all->xi_psi_prime.chip(k, 1);
      xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx.chip(k, 2) =
          prods->st_xi_psi_all->xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx.chip(k, 1);
      xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx.chip(k, 2) =
          prods->st_xi_psi_all->xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx.chip(k, 1);
    }
    ArrayXXc<Real> for_Q_diag_Lt1 = prods->st_xi_psi_all->for_diag_Lt1;
    ArrayXXc<Real> for_Q_diag_Lt2 = prods->st_xi_psi_all->for_diag_Lt2;
    ArrayXXc<Real> for_Q_diag_Lt3 = prods->st_xi_psi_all->for_diag_Lt3;
    delete prods->st_xi_psi_all;

    Tensor3c<Real> psi_prime_psi(N_max + 1, T, N_max + 1);
    Tensor3c<Real> psi_psi_prime(N_max + 1, T, N_max + 1);
    Tensor3c<Real> psi_psi(N_max + 1, T, N_max + 1);
    Tensor3c<Real> psi_prime_psi_prime_plus_kkp1_psi_psi_over_sxx(N_max + 1, T, N_max + 1);
    Tensor3c<Real> psi_prime_psi_prime_plus_nnp1_psi_psi_over_sxx(N_max + 1, T, N_max + 1);
    for (int k = 0; k <= N_max; k++) {
      psi_psi.chip(k, 2) = prods->st_psi_psi_all->xi_psi.chip(k, 1);
      psi_prime_psi.chip(k, 2) = prods->st_psi_psi_all->xi_prime_psi.chip(k, 1);
      psi_psi_prime.chip(k, 2) = prods->st_psi_psi_all->xi_psi_prime.chip(k, 1);
      psi_prime_psi_prime_plus_kkp1_psi_psi_over_sxx.chip(k, 2) =
          prods->st_psi_psi_all->xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx.chip(k, 1);
      psi_prime_psi_prime_plus_nnp1_psi_psi_over_sxx.chip(k, 2) =
          prods->st_psi_psi_all->xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx.chip(k, 1);
    }
    ArrayXXc<Real> for_P_diag_Lt1 = prods->st_psi_psi_all->for_diag_Lt1;
    ArrayXXc<Real> for_P_diag_Lt2 = prods->st_psi_psi_all->for_diag_Lt2;
    ArrayXXc<Real> for_P_diag_Lt3 = prods->st_psi_psi_all->for_diag_Lt3;
    delete prods->st_psi_psi_all;
    delete prods;

    for (int i = 0; i < M; i++) {
      int m = abs_m_vec(i);
      int Nm = N_max - m + 1;
      ArrayXi n_vec = ArrayXi::LinSpaced(Nm, m, N_max);
      ArrayXi p_vec = n_vec*(n_vec + 1) + m;
      ArrayXr<Real> n_vec_real = n_vec.template cast<Real>();
      ArrayXXr<Real> pi_nm = stPT->pi_nm(all, p_vec).transpose();
      ArrayXXr<Real> tau_nm = stPT->tau_nm(all, p_vec).transpose();
      ArrayXXr<Real> d_n;
      if (m)
        d_n = pi_nm.rowwise() * (sin_t/m).transpose();
      else
        d_n = stPT->p_n0(all, n_vec).transpose();

      ArrayXXr<Real> d_n_times_nnp1 = d_n.colwise() * (n_vec_real*(n_vec_real + 1));

      ArrayXXc<Real> K1, K2, L5, L6, K1P, K2P, L5P, L6P;
      K1 = K2 = L5 = L6 = K1P = K2P = L5P = L6P = ArrayXXc<Real>::Zero(Nm, Nm);

      for (int k = m; k <= N_max; k++) {
        int k_ind = k - m;
        ArrayXr<Real> d_k = d_n.row(k_ind).transpose();
        ArrayXr<Real> tau_k = tau_nm.row(k_ind).transpose();

        VectorXr<Real> dx_dt_tau_k_sin_t = dx_dt_wt * tau_k;
        VectorXr<Real> dx_dt_d_k_sin_t = dx_dt_wt * d_k;

        MatrixXc<Real> pi_n_xi_prime_psi = pi_nm * reduceAndSlice(xi_prime_psi, k, Nm);
        MatrixXc<Real> pi_n_xi_psi_prime = pi_nm * reduceAndSlice(xi_psi_prime, k, Nm);
        K1.col(k_ind) = pi_n_xi_psi_prime * dx_dt_d_k_sin_t;
        K2.col(k_ind) = pi_n_xi_prime_psi * dx_dt_d_k_sin_t;

        pi_n_xi_prime_psi = pi_nm * reduceAndSlice(psi_prime_psi, k, Nm);
        pi_n_xi_psi_prime = pi_nm * reduceAndSlice(psi_psi_prime, k, Nm);
        K1P.col(k_ind) = pi_n_xi_psi_prime * dx_dt_d_k_sin_t;
        K2P.col(k_ind) = pi_n_xi_prime_psi * dx_dt_d_k_sin_t;

        MatrixXc<Real> d_n_xi_psi_nnp1 = d_n_times_nnp1 * reduceAndSlice(xi_psi, k, Nm);
        MatrixXc<Real> tau_n_xi_psi = tau_nm * reduceAndSlice(xi_psi, k, Nm);

        L5.col(k_ind) = d_n_xi_psi_nnp1 * dx_dt_tau_k_sin_t - tau_n_xi_psi * dx_dt_d_k_sin_t*k*(k + 1);

        d_n_xi_psi_nnp1 = d_n_times_nnp1 * reduceAndSlice(psi_psi, k, Nm);
        tau_n_xi_psi = tau_nm *  reduceAndSlice(psi_psi, k, Nm);

        L5P.col(k_ind) = d_n_xi_psi_nnp1 * dx_dt_tau_k_sin_t - tau_n_xi_psi * dx_dt_d_k_sin_t*k*(k + 1);

        MatrixXc<Real> d_n_xi_prime_psi_prime_nnp1_plus_xi_psi_over_sxx_nnp1_kkp1 =
            d_n_times_nnp1 * reduceAndSlice(xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx, k, Nm);
        MatrixXc<Real> tau_n_xi_prime_psi_prime_plus_xi_psi_over_sxx_nnp1 =
            tau_nm * reduceAndSlice(xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx, k, Nm);

        L6.col(k_ind) = d_n_xi_prime_psi_prime_nnp1_plus_xi_psi_over_sxx_nnp1_kkp1 * dx_dt_tau_k_sin_t -
            tau_n_xi_prime_psi_prime_plus_xi_psi_over_sxx_nnp1 * dx_dt_d_k_sin_t * k*(k + 1);

        d_n_xi_prime_psi_prime_nnp1_plus_xi_psi_over_sxx_nnp1_kkp1 =
            d_n_times_nnp1 * reduceAndSlice(psi_prime_psi_prime_plus_kkp1_psi_psi_over_sxx, k, Nm);
        tau_n_xi_prime_psi_prime_plus_xi_psi_over_sxx_nnp1 =
            tau_nm * reduceAndSlice(psi_prime_psi_prime_plus_nnp1_psi_psi_over_sxx, k, Nm);

        L6P.col(k_ind) = d_n_xi_prime_psi_prime_nnp1_plus_xi_psi_over_sxx_nnp1_kkp1 * dx_dt_tau_k_sin_t -
            tau_n_xi_prime_psi_prime_plus_xi_psi_over_sxx_nnp1 * dx_dt_d_k_sin_t * k*(k + 1);
      }

      ArrayXXc<Real> tmp2 = (n_vec_real * (n_vec_real + 1)).replicate(1, Nm).rowwise() -
      (n_vec_real * (n_vec_real + 1)).transpose();
      ArrayXXc<Real> prefactor1 = (s - 1)*(s + 1) / s * An_Ak(n_vec, n_vec);
      ArrayXXc<Real> prefactor2 = I * prefactor1 / tmp2;

      ArrayXXc<Real> Q12 = prefactor1 * K1;
      ArrayXXc<Real> Q21 = -prefactor1 * K2;
      ArrayXXc<Real> Q11 = prefactor2 * L5;
      ArrayXXc<Real> Q22 = prefactor2 * L6;

      ArrayXXc<Real> P12 = prefactor1 * K1P;
      ArrayXXc<Real> P21 = -prefactor1 * K2P;
      ArrayXXc<Real> P11 = prefactor2 * L5P;
      ArrayXXc<Real> P22 = prefactor2 * L6P;

      ArrayXc<Real> prefact_diag1 = -I/s * (2*n_vec_real + 1) / (2*n_vec_real*(n_vec_real + 1));
      ArrayXc<Real> prefact_diag2 = -I*(s - 1)*(s + 1)/s/2.0 * (2*n_vec_real + 1);
      ArrayXXc<Real> pi2_p_tau2 = (pi_nm.pow(2) + tau_nm.pow(2));

      ArrayXc<Real> Ltilde1 = ((pi2_p_tau2 * for_Q_diag_Lt1(n_vec, all)).matrix() * Rt_func->w_theta.matrix()).array();
      ArrayXc<Real> Ltilde2 = ((pi2_p_tau2 * for_Q_diag_Lt2(n_vec, all)).matrix() * Rt_func->w_theta.matrix()).array();
      ArrayXc<Real> Ltilde3 = ((d_n * tau_nm * for_Q_diag_Lt3(n_vec, all)).matrix() * dx_dt_wt.matrix()).array();

      for (int j = 0; j < Nm; j++) {
        Q11(j, j) = prefact_diag1(j) * Ltilde1(j);
        Q22(j, j) = prefact_diag1(j) * Ltilde2(j) + prefact_diag2(j) * Ltilde3(j);
      }

      Ltilde1 = (pi2_p_tau2 * for_P_diag_Lt1(n_vec, all)).matrix() * Rt_func->w_theta.matrix();
      Ltilde2 = (pi2_p_tau2 * for_P_diag_Lt2(n_vec, all)).matrix() * Rt_func->w_theta.matrix();
      Ltilde3 = (d_n * tau_nm * for_P_diag_Lt3(n_vec, all)).matrix() * dx_dt_wt.matrix();

      for (int j = 0; j < Nm; j++) {
        P11(j, j) = prefact_diag1(j) * Ltilde1(j);
        P22(j, j) = prefact_diag1(j) * Ltilde2(j) + prefact_diag2(j) * Ltilde3(j);
      }

      ArrayXi inde = seq2Array(m % 2, Nm - 1, 2);
      ArrayXi indo = seq2Array(1 - m % 2, Nm - 1, 2);

      output[i].st4MQeo = new st4M<Real>();
      output[i].st4MQeo->M12 = Q12(inde, indo);
      output[i].st4MQeo->M21 = Q21(indo, inde);
      output[i].st4MQeo->M11 = Q11(inde, inde);
      output[i].st4MQeo->M22 = Q22(indo, indo);
      output[i].st4MQeo->m = m;
      output[i].st4MQeo->ind1 = inde;
      output[i].st4MQeo->ind2 = indo;

      output[i].st4MQoe = new st4M<Real>();
      output[i].st4MQoe->M12 = Q12(indo, inde);
      output[i].st4MQoe->M21 = Q21(inde, indo);
      output[i].st4MQoe->M11 = Q11(indo, indo);
      output[i].st4MQoe->M22 = Q22(inde, inde);
      output[i].st4MQoe->m = m;
      output[i].st4MQoe->ind1 = inde;
      output[i].st4MQoe->ind2 = indo;

      output[i].st4MPeo = new st4M<Real>();
      output[i].st4MPeo->M12 = P12(inde, indo);
      output[i].st4MPeo->M21 = P21(indo, inde);
      output[i].st4MPeo->M11 = P11(inde, inde);
      output[i].st4MPeo->M22 = P22(indo, indo);
      output[i].st4MPeo->m = m;
      output[i].st4MPeo->ind1 = inde;
      output[i].st4MPeo->ind2 = indo;

      output[i].st4MPoe = new st4M<Real>();
      output[i].st4MPoe->M12 = P12(indo, inde);
      output[i].st4MPoe->M21 = P21(inde, indo);
      output[i].st4MPoe->M11 = P11(indo, indo);
      output[i].st4MPoe->M22 = P22(inde, inde);
      output[i].st4MPoe->m = m;
      output[i].st4MPoe->ind1 = inde;
      output[i].st4MPoe->ind2 = indo;

      output[i].has_st4MQ = true;
      output[i].has_st4MP = true;
    }
    delete stPT;
    return output;
  }

  template ArrayXXr<double>* sphGetUforFp(int);
  template stFprow<double>* sphGetFpRow(int, double, const ArrayXr<double>&);
  template stFpovx<double>* sphGetFpovx(int, double, const ArrayXr<double>&);
  template stBessel<double>* sphGetXiPsi(int, double, const ArrayXr<double>&, int);
  template stBesselPrimes<double>* sphGetBesselProductsPrimes(const Tensor3c<double>&);
  template stBesselProducts<double>* sphGetModifiedBesselProducts(int, double, const ArrayXr<double>&, int);
  template stRtfunc<double>* sphMakeGeometry(size_t, double, double, const ArrayXr<double>*);
  template size_t sphCheckBesselConvergence(size_t, double, const ArrayXr<double>&, double, size_t);
  template size_t sphEstimateNB(size_t, const stRtfunc<double>*, const stParams<double>*, double);
  template stPQa<double>* sphCalculatePQ(int, const ArrayXi&, const stRtfunc<double>*, const stParams<double>*, int);
}
