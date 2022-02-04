#ifndef SPH_HPP
#define SPH_HPP

#include "core.hpp"
#include "smarties_aux.hpp"
#include "vsh.hpp"
#include "misc.hpp"
#include <stdexcept>

using namespace Eigen;
using namespace std;

namespace Smarties {

  template <class Real>
  struct stFprow {
    ArrayXXr<Real> S;
    ArrayXXr<Real> loss_prec_S;
  };

  template <class Real>
  struct stFpovx {
    Tensor3c<Real> Fpovx;
    ArrayXXc<Real> rb_chi;
    ArrayXXc<Real> rb_psi;
  };

  template <class Real>
  struct stBessel {
    Tensor3c<Real> xi_psi;
    Tensor3c<Real> psi_psi;
    ArrayXXc<Real> chi_n;
    ArrayXXc<Real> psi_n;
    ArrayXXc<Real> psi_k;
  };

  template <class Real>
  struct stBesselPrimes {
    Tensor3c<Real> xi_psi;
    Tensor3c<Real> xi_prime_psi;
    Tensor3c<Real> xi_psi_prime;
    Tensor3c<Real> xi_prime_psi_prime;
    Tensor3c<Real> xi_psi_over_sxx;
    Tensor3c<Real> xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx;
    Tensor3c<Real> xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx;
    ArrayXXc<Real> for_diag_Lt1;
    ArrayXXc<Real> for_diag_Lt2;
    ArrayXXc<Real> for_diag_Lt3;
  };

  template <class Real>
  struct stBesselProducts {
    unique_ptr<stBesselPrimes<Real>> st_xi_psi_all;
    unique_ptr<stBesselPrimes<Real>> st_psi_psi_all;
  };

  template <class Real>
  struct st4M {
    ArrayXXc<Real> M11;
    ArrayXXc<Real> M12;
    ArrayXXc<Real> M21;
    ArrayXXc<Real> M22;
    int m;
    ArrayXi ind1;
    ArrayXi ind2;
  };

  template <class Real>
  struct stMat {
    std::array<st4M<Real>, 4> st_4M_list;
    vector<string> mat_list;

    stMat() {
      for (int i = 0; i < 4; i++)
        this->st_4M_list[i] = st4M<Real>();
    }

    stMat(const st4M<Real>& base) {
      for (int i = 0; i < 4; i++)
        this->st_4M_list[i] = base.st_4M_list[i];
      this->mat_list = base.mat_list;
    }
  };

  template <class Real>
  struct stPQ : stMat<Real> {
    using stMat<Real>::stMat;
    inline st4M<Real>& st_4M_P_eo() { return this->st_4M_list[0]; }
    inline st4M<Real>& st_4M_P_oe() { return this->st_4M_list[1]; }
    inline st4M<Real>& st_4M_Q_eo() { return this->st_4M_list[2]; }
    inline st4M<Real>& st_4M_Q_oe() { return this->st_4M_list[3]; }
  };

  template <class Real>
  unique_ptr<ArrayXXr<Real>> sphGetUforFp(int n) {
    int b_max = n/2;
    auto u = make_unique<ArrayXXr<Real>>(b_max + 2, b_max + 1);
    u->setZero();
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
  unique_ptr<stFprow<Real>> sphGetFpRow(int n, Real s, const ArrayXr<Real>& x) {
    int num_x = x.size();
    RowArrayXr<Real> x_squared = x.pow(2).transpose();
    unique_ptr<ArrayXXr<Real>> u = sphGetUforFp<Real>(n);
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

    int num_to_test = 3;

    for (int k = n - 4; k >= 0; k -= 2) {
      q_min = (n - k)/2 - 1; // expression is mathematically always an integer
      q_int = n - k;
      alpha_bar_k /= -(n - k - 2);
      beta(0, q_int - 2) = 1;
      for (int j = 1; j < q_int; j++)
        beta(j, q_int - 2) = beta(j - 1, q_int - 2)*(s - 1)*(s + 1)*(q_int - j)/j;
      beta(0, q_int - 1) = 1;
      for (int j = 1; j <= q_int; j++)
        beta(j, q_int - 1) = beta(j - 1, q_int - 1)*(s - 1)*(s + 1)*(q_int - j + 1)/j;

      RowArrayXr<Real> alpha_q(num_x);
      alpha_q.fill(alpha_bar_k);
      RowArrayXb x_not_converged = RowArrayXb::Ones(num_x);
      RowArrayXr<Real> current_term = RowArrayXr<Real>::Zero(num_x);
      RowArrayXi counter = RowArrayXi::Zero(num_x);
      RowArrayXb test = RowArrayXb::Zero(num_x);

      for (int q = q_min; q < q_int; q++) {
        int b = n - k - q - 1;
        Real gamma_qnk = 0;
        for (int j = max(0, 2 * q + k - n + 1); j <= q; j++)
          gamma_qnk += beta(j, q - 1) * (*u)(q - j, b);

        if (abs(s) < 2) {
          current_term = alpha_q * gamma_qnk;
          S.row(k) += current_term;
          max_term_S.row(k) = max_term_S.row(k).max(abs(current_term));
        } else {
          test.fill(false);
          current_term = x_not_converged.select(alpha_q * gamma_qnk, current_term);
          if (x_not_converged.any())
            test = x_not_converged.select(abs(current_term) > mp_eps<Real>() * abs(S.row(k)), test);
          if (test.any()) {
            if (LogicalSlice(current_term, test).isInf().all()) // if all non-converged x values are infinite
              cout << "Problem (1) in sphGetFpovx..." << endl;
            S.row(k) += test.select(current_term, 0);
            max_term_S.row(k) = test.select(abs(current_term).max(max_term_S.row(k)), max_term_S.row(k));
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

      Real c_qqnk = static_cast<Real>(1.0)/(2*k + 1);
      int q = q_int;
      while (true) {
        Real c_iqnk = c_qqnk;
        Real gamma_qnk = c_iqnk;
        for (int i = q - 1; i >= 0; i--) {
          c_iqnk *= pow(s, 2)*(i + 1)*(2*i + 1 - 2*n)/(q - i)/(2*k + 2*q - 2*i + 1);
          gamma_qnk += c_iqnk;
        }

        test.fill(false);
        for (int i = 0; i < num_x; i++) {
          if (x_not_converged(i)) {
            current_term(i) = alpha_q(i) * gamma_qnk;
            test(i) = abs(current_term(i)) > mp_eps<Real>() * abs(S(k, i));
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

    auto output = make_unique<stFprow<Real>>();
    output->S = S;
    output->loss_prec_S = loss_prec_S;
    return output;
  }

  template <class Real>
  unique_ptr<stFpovx<Real>> sphGetFpovx(int N_max, Real s, const ArrayXr<Real>& x) {
    int num_x = x.size();

    auto output = make_unique<stFpovx<Real>>();
    output->rb_chi = vshRBchi<Real>(ArrayXr<Real>::LinSpaced(N_max + 1, 0, N_max), x);
    output->rb_psi = vshRBpsi<Real>(ArrayXr<Real>::LinSpaced(N_max + 1, 0, N_max), s*x);
    output->Fpovx = Tensor3c<Real>(N_max + 1, N_max + 1, num_x);
    output->Fpovx.setZero();

    unique_ptr<stFprow<Real>> FpRow = sphGetFpRow<Real>(N_max, s, x);
    ArrayXXc<Real> tmp;
    for (int i = 0; i < N_max - 3; i++) {
      tmp = (FpRow->S.row(i) / x.transpose()).template cast<complex<Real>>();
      output->Fpovx.chip(N_max, 0).chip(i, 0) = TensorCast(tmp, num_x);
    }

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
  unique_ptr<stBessel<Real>> sphGetXiPsi(int N_max, Real s, const ArrayXr<Real>& x, int NB) {
    int num_x = x.size();
    auto output = make_unique<stBessel<Real>>();
    unique_ptr<stFpovx<Real>> chi_psi_prod = sphGetFpovx<Real>(NB + 1, s, x);
    output->psi_psi = Tensor3c<Real>(N_max + 2, N_max + 2, num_x);
    output->psi_psi.setZero();
    output->psi_n = vshRBpsi<Real>(ArrayXr<Real>::LinSpaced(N_max + 2, 0, N_max + 1), x);

    for (int n = 0; n < N_max + 2; n++) {
      for (int k = n % 2; k < N_max + 2; k += 2)
        output->psi_psi.chip(n, 0).chip(k, 0) =
            TensorCast(chi_psi_prod->rb_psi.col(k) * output->psi_n.col(n), num_x);
    }

    std::array<int, 3> offsets = {0, 0, 0};
    std::array<int, 3> extents = {N_max + 2, N_max + 2, num_x};
    output->xi_psi = output->psi_psi + output->psi_psi.constant(mp_im_unit<Real>()) *
        chi_psi_prod->Fpovx.slice(offsets, extents);
    output->chi_n = chi_psi_prod->rb_chi.leftCols(N_max + 2);
    output->psi_k = chi_psi_prod->rb_psi.leftCols(N_max + 2);
    return output;
  }

  // Note: the products won't be correct at n = 0 or k = 0
  // Expects: prods = [(N + 2) x (N + 2) x X] tensor
  template <class Real>
  unique_ptr<stBesselPrimes<Real>> sphGetBesselProductsPrimes(const Tensor3c<Real>& prods) {
    int N = prods.dimension(0) - 2;
    int X = prods.dimension(2);

    auto output = make_unique<stBesselPrimes<Real>>();
    std::array<long int, 3> offsets = {0, 0, 0};
    std::array<long int, 3> extents = {N + 1, N + 1, X};
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

    Tensor1c<Real> col_shape(X);
    for (int n = 1; n <= N; n++) {
      for (int k = 2 - n % 2; k <= N; k += 2) {
        int kkp1 = k*(k + 1);
        int nnp1 = n*(n + 1);

        output->xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx.chip(n, 0).chip(k, 0) = (
            col_shape.constant(static_cast<complex<Real>>((k + n + 1)*(k + 1))) * prods.chip(n - 1, 0).chip(k - 1, 0) +
            col_shape.constant(static_cast<complex<Real>>(kkp1 - k*(n + 1))) * prods.chip(n - 1, 0).chip(k + 1, 0) +
            col_shape.constant(static_cast<complex<Real>>(kkp1 - (k + 1)*n)) * prods.chip(n + 1, 0).chip(k - 1, 0) +
            col_shape.constant(static_cast<complex<Real>>(kkp1 + k*n)) * prods.chip(n + 1, 0).chip(k + 1, 0)) /
            col_shape.constant(static_cast<complex<Real>>((2*n + 1)*(2*k + 1)));

        output->xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx.chip(n, 0).chip(k, 0) = (
            col_shape.constant(static_cast<complex<Real>>(nnp1 + (k + 1) * (n + 1)))*prods.chip(n - 1, 0).chip(k - 1, 0) +
            col_shape.constant(static_cast<complex<Real>>((n - k)*(n + 1))) * prods.chip(n - 1, 0).chip(k + 1, 0) +
            col_shape.constant(static_cast<complex<Real>>((n - k)*n)) * prods.chip(n + 1, 0).chip(k - 1, 0) +
            col_shape.constant(static_cast<complex<Real>>(nnp1 + k*n)) * prods.chip(n + 1, 0).chip(k + 1, 0)) /
            col_shape.constant(static_cast<complex<Real>>((2*n + 1)*(2*k + 1)));
      }

      for (int k = 1 + n % 2; k <= N; k += 2) {
        output->xi_prime_psi.chip(n, 0).chip(k, 0) =
            (prods.chip(n - 1, 0).chip(k, 0) * col_shape.constant(static_cast<complex<Real>>(n + 1)) -
            col_shape.constant(static_cast<complex<Real>>(n))*prods.chip(n + 1, 0).chip(k, 0)) /
            col_shape.constant(static_cast<complex<Real>>(2*n + 1));
        output->xi_psi_prime.chip(n, 0).chip(k, 0) =
            (prods.chip(n, 0).chip(k - 1, 0) * col_shape.constant(static_cast<complex<Real>>(k + 1)) -
            col_shape.constant(static_cast<complex<Real>>(k))*prods.chip(n, 0).chip(k + 1, 0)) /
            col_shape.constant(static_cast<complex<Real>>(2*k + 1));
      }
    }
    return output;
  }

  // Expects: NB >= N_max
  template <class Real>
  unique_ptr<stBesselProducts<Real>> sphGetModifiedBesselProducts(int N_max, Real s, const ArrayXr<Real>& x, int NB) {
    auto output = make_unique<stBesselProducts<Real>>();
    unique_ptr<stBessel<Real>> prods = sphGetXiPsi<Real>(N_max, s, x, NB);
    output->st_xi_psi_all = sphGetBesselProductsPrimes<Real>(prods->xi_psi);
    output->st_psi_psi_all = sphGetBesselProductsPrimes<Real>(prods->psi_psi);

    ArrayXXc<Real> psi_np1_psi_n = (prods->psi_n(all, seq(1, last)) * prods->psi_k(all, seq(0, last - 1))).transpose();
    ArrayXXc<Real> psi_n_psi_np1 = (prods->psi_n(all, seq(0, last - 1)) * prods->psi_k(all, seq(1, last))).transpose();
    ArrayXXc<Real> xi_np1_psi_n = psi_np1_psi_n + mp_im_unit<Real>() *
        (prods->chi_n(all, seq(1, last)) * prods->psi_k(all, seq(0, last - 1))).transpose();
    ArrayXXc<Real> xi_n_psi_np1 = psi_n_psi_np1 + mp_im_unit<Real>() *
        (prods->chi_n(all, seq(0, last - 1)) * prods->psi_k(all, seq(1, last))).transpose();
    output->st_psi_psi_all->for_diag_Lt1 = s*psi_n_psi_np1 - psi_np1_psi_n;
    output->st_xi_psi_all->for_diag_Lt1 = s*xi_n_psi_np1 - xi_np1_psi_n;
    ArrayXXc<Real> psi_n_psi_n = (prods->psi_n(all, seq(0, last - 1)) * prods->psi_k(all, seq(0, last - 1))).transpose();
    ArrayXXc<Real> xi_n_psi_n = psi_n_psi_n + mp_im_unit<Real>() *
        (prods->chi_n(all, seq(0, last - 1)) * prods->psi_k(all, seq(0, last - 1))).transpose();

    ArrayXc<Real> n_vec = ArrayXr<Real>::LinSpaced(N_max + 1, 1, N_max + 1);
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
  unique_ptr<stRtfunc<Real>> sphMakeGeometry(int Nb_theta, Real a, Real c,
      const unique_ptr<ArrayXr<Real>> theta = unique_ptr<ArrayXr<Real>>()) {
    sInt type;
    auto output = make_unique<stRtfunc<Real>>();
    if (theta) {
      output->w_theta = ArrayXr<Real>::Zero(theta->size());
      output->theta = *theta;
      output->Nb_theta = theta->size();
      type = PTS;
    } else {
      type = GAUSS;
      unique_ptr<stRtfunc<Real>> tmp = auxPrepareIntegrals<Real>(2*Nb_theta, type);
      output->theta = tmp->theta(seq(0, Nb_theta - 1));
      output->w_theta = tmp->w_theta(seq(0, Nb_theta - 1))*2;
      output->Nb_theta = Nb_theta;
    }

    output->a = a;
    output->c = c;
    output->type = type;

    ArrayXr<Real> sin_t = sin(output->theta);
    ArrayXr<Real> cos_t = cos(output->theta);

    output->r = a*c/sqrt(pow(c*sin_t, 2) + pow(a*cos_t, 2));
    output->dr_dt = (pow(a, 2) - pow(c, 2))/pow(a*c, 2)*sin_t * cos_t * output->r.pow(3);

    output->h = max(a, c)/min(a, c);
    output->r0 = pow(pow(a, 2)*c, static_cast<Real>(1.0)/3);

    return output;
  }

  template <class Real>
  int sphCheckBesselConvergence(int N_req, Real s, const ArrayXr<Real>& x, Real acc, int N_min) {
    int NB_start = max(N_req, N_min);
    int NB = NB_start;
    unique_ptr<stFpovx<Real>> prod = sphGetFpovx<Real>(NB, s, x);
    bool to_continue = true;
    int max_N = 500;
    int NB_step = 16;

    int NB_next;
    unique_ptr<stFpovx<Real>> prod_new;
    Tensor<Real, 0> rel_acc_ee, rel_acc_oo;
    Real rel_acc;
    ArithmeticSequence<long int, long int, long int> seq1(0, N_req, 2);
    ArithmeticSequence<long int, long int, long int> seq2(1, N_req, 2);
    ArithmeticSequence<long int, long int, long int> seq3(0, x.size() - 1); // prod->Fpovx has no. of columns = x.size()
    long int dim1 = seq1.size();
    long int dim3 = seq3.size();
    Tensor3c<Real> tmp(dim1, dim1, dim3);
    while (to_continue && NB < max_N) {
      NB_next = NB + NB_step;
      prod_new = sphGetFpovx<Real>(NB_next, s, x);
      tmp = TensorSlice<Real>(prod->Fpovx, seq1, seq1, seq3) /
          TensorSlice<Real>(prod_new->Fpovx, seq1, seq1, seq3) - tmp.constant(static_cast<complex<Real>>(1));
      rel_acc_ee = tmp.abs().real().maximum();
      tmp = TensorSlice<Real>(prod->Fpovx, seq2, seq2, seq3) /
          TensorSlice<Real>(prod_new->Fpovx, seq2, seq2, seq3) - tmp.constant(static_cast<complex<Real>>(1));
      rel_acc_oo = tmp.abs().real().maximum();
      rel_acc = max(rel_acc_ee(0), rel_acc_oo(0));
      if (rel_acc < acc)
        to_continue = false;
      else {
        NB = NB_next;
        prod = move(prod_new);
      }
    }

    if (NB > NB_start) {
      NB -= NB_step;
      NB_step = 1;
      to_continue = true;
      while (to_continue && NB < max_N) {
        NB += NB_step;
        prod = sphGetFpovx<Real>(NB, s, x);
        tmp = TensorSlice<Real>(prod->Fpovx, seq1, seq1, seq3) /
            TensorSlice<Real>(prod_new->Fpovx, seq1, seq1, seq3) - tmp.constant(static_cast<complex<Real>>(1));
        rel_acc_ee = tmp.abs().maximum().real();
        tmp = TensorSlice<Real>(prod->Fpovx, seq2, seq2, seq3) /
            TensorSlice<Real>(prod_new->Fpovx, seq2, seq2, seq3) - tmp.constant(static_cast<complex<Real>>(1));
        rel_acc_oo = tmp.abs().maximum().real();
        rel_acc = max(rel_acc_ee(0), rel_acc_oo(0));
        if (rel_acc < acc)
          to_continue = false;
      }
    }

    if (to_continue)
      throw(runtime_error("Problem in sphEstimateNB: convergence was not achieved"));
    return NB;
  }

  template <class Real>
  int sphEstimateNB(int NQ, const unique_ptr<stRtfunc<Real>>& stGeometry,
      const unique_ptr<stParams<Real>>& params, Real acc = 1e-13) {
    ArrayXr<Real> s = params->s;
    ArrayXr<Real> k1 = params->k1;

    ArrayXr<Real> x_max = {{k1.maxCoeff() * stGeometry->r.maxCoeff()}};
    typename ArrayXr<Real>::Index ind;
    abs(s).minCoeff(&ind);
    Real s_min = s(ind);
    abs(s).maxCoeff(&ind);
    Real s_max = s(ind);

    int N1 = sphCheckBesselConvergence(NQ, s_max, x_max, acc, NQ);
    int NB = sphCheckBesselConvergence(NQ, s_min, x_max, acc, N1);
    return NB;
  }

  template <class Real>
  vector<unique_ptr<stPQ<Real>>> sphCalculatePQ(int N_max, const ArrayXi& abs_m_vec,
      const unique_ptr<stRtfunc<Real>>& Rt_func, const unique_ptr<stParams<Real>>& params, int NB = -1) {
    if (params->s.size() > 1 || params->k1.size() > 1)
      throw(runtime_error("params->s and params->k1 must be scalar"));
    if (NB < N_max)
      NB = N_max;
    int M = abs_m_vec.size();
    //if (params->output) (params shouldn't have an output member, perhaps use stOptions instead)
    cout << "sphCalculatePQ: Calculate P, Q for " << M << " m-values with N_Q = " <<
        N_max << ", NB = " << NB << ", N_Theta = " << Rt_func->Nb_theta << endl;

    vector<unique_ptr<stPQ<Real>>> output(M);
    for (int i = 0; i < M; i++)
      output[i] = make_unique<stPQ<Real>>();
    Real s = params->s(0);
    Real k1 = params->k1(0);
    ArrayXr<Real> x = Rt_func->r * k1; // [T x 1] x(theta)
    ArrayXr<Real> x_theta = Rt_func->dr_dt * k1; // [T x 1] x'(theta)

    unique_ptr<stPinmTaunm<Real>> stPT = vshPinmTaunm(N_max, Rt_func->theta);
    ArrayXr<Real> sin_t = sin(Rt_func->theta); // [T x 1]
    ArrayXr<Real> dx_dt_wt = x_theta * Rt_func->w_theta; // [T x 1]
    ArrayXr<Real> tmp1 = ArrayXr<Real>::LinSpaced(N_max + 1, 0, N_max); // [N x 1]
    ArrayXr<Real> An_vec = sqrt((2*tmp1 + 1) / (2*tmp1 * (tmp1 + 1))); // [N x 1]
    ArrayXXr<Real> An_Ak = (An_vec.matrix() * An_vec.transpose().matrix()).array(); // [N x N]

    std::array<int, 3> shuffle = {0, 2, 1};
    unique_ptr<stBesselProducts<Real>> prods = sphGetModifiedBesselProducts(N_max, s, x, NB); // prods contains [N x N x T] tensors
    // Converting these [N x N x T] tensors into [N x T x N] tensors
    Tensor3c<Real> xi_psi = prods->st_xi_psi_all->xi_psi.shuffle(shuffle);
    Tensor3c<Real> xi_prime_psi = prods->st_xi_psi_all->xi_prime_psi.shuffle(shuffle);
    Tensor3c<Real> xi_psi_prime = prods->st_xi_psi_all->xi_psi_prime.shuffle(shuffle);
    Tensor3c<Real> xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx =
        prods->st_xi_psi_all->xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx.shuffle(shuffle);
    Tensor3c<Real> xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx =
        prods->st_xi_psi_all->xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx.shuffle(shuffle);
    ArrayXXc<Real> for_Q_diag_Lt1 = prods->st_xi_psi_all->for_diag_Lt1;
    ArrayXXc<Real> for_Q_diag_Lt2 = prods->st_xi_psi_all->for_diag_Lt2;
    ArrayXXc<Real> for_Q_diag_Lt3 = prods->st_xi_psi_all->for_diag_Lt3;

    Tensor3c<Real> psi_psi = prods->st_psi_psi_all->xi_psi.shuffle(shuffle);
    Tensor3c<Real> psi_prime_psi = prods->st_psi_psi_all->xi_prime_psi.shuffle(shuffle);
    Tensor3c<Real> psi_psi_prime = prods->st_psi_psi_all->xi_psi_prime.shuffle(shuffle);
    Tensor3c<Real> psi_prime_psi_prime_plus_kkp1_psi_psi_over_sxx =
        prods->st_psi_psi_all->xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx.shuffle(shuffle);
    Tensor3c<Real> psi_prime_psi_prime_plus_nnp1_psi_psi_over_sxx =
        prods->st_psi_psi_all->xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx.shuffle(shuffle);
    ArrayXXc<Real> for_P_diag_Lt1 = prods->st_psi_psi_all->for_diag_Lt1;
    ArrayXXc<Real> for_P_diag_Lt2 = prods->st_psi_psi_all->for_diag_Lt2;
    ArrayXXc<Real> for_P_diag_Lt3 = prods->st_psi_psi_all->for_diag_Lt3;

    for (int i = 0; i < M; i++) {
      int m = abs_m_vec(i);
      int N_min = max(m, 1);
      int Nm = N_max - N_min + 1;
      ArrayXi n_vec = ArrayXi::LinSpaced(Nm, N_min, N_max);
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

      for (int k = N_min; k <= N_max; k++) {
        int k_ind = k - N_min;
        ArrayXr<Real> d_k = d_n.row(k_ind).transpose();
        ArrayXr<Real> tau_k = tau_nm.row(k_ind).transpose();

        VectorXr<Real> dx_dt_tau_k_sin_t = dx_dt_wt * tau_k;
        VectorXr<Real> dx_dt_d_k_sin_t = dx_dt_wt * d_k;

        MatrixXc<Real> pi_n_xi_prime_psi = pi_nm * ReduceAndSlice(xi_prime_psi, k, Nm);
        MatrixXc<Real> pi_n_xi_psi_prime = pi_nm * ReduceAndSlice(xi_psi_prime, k, Nm);
        K1.col(k_ind) = pi_n_xi_psi_prime * dx_dt_d_k_sin_t;
        K2.col(k_ind) = pi_n_xi_prime_psi * dx_dt_d_k_sin_t;

        pi_n_xi_prime_psi = pi_nm * ReduceAndSlice(psi_prime_psi, k, Nm);
        pi_n_xi_psi_prime = pi_nm * ReduceAndSlice(psi_psi_prime, k, Nm);
        K1P.col(k_ind) = pi_n_xi_psi_prime * dx_dt_d_k_sin_t;
        K2P.col(k_ind) = pi_n_xi_prime_psi * dx_dt_d_k_sin_t;

        MatrixXc<Real> d_n_xi_psi_nnp1 = d_n_times_nnp1 * ReduceAndSlice(xi_psi, k, Nm);
        MatrixXc<Real> tau_n_xi_psi = tau_nm * ReduceAndSlice(xi_psi, k, Nm);

        L5.col(k_ind) = d_n_xi_psi_nnp1 * dx_dt_tau_k_sin_t - tau_n_xi_psi * dx_dt_d_k_sin_t*k*(k + 1);

        d_n_xi_psi_nnp1 = d_n_times_nnp1 * ReduceAndSlice(psi_psi, k, Nm);
        tau_n_xi_psi = tau_nm *  ReduceAndSlice(psi_psi, k, Nm);

        L5P.col(k_ind) = d_n_xi_psi_nnp1 * dx_dt_tau_k_sin_t - tau_n_xi_psi * dx_dt_d_k_sin_t*k*(k + 1);

        MatrixXc<Real> d_n_xi_prime_psi_prime_nnp1_plus_xi_psi_over_sxx_nnp1_kkp1 =
            d_n_times_nnp1 * ReduceAndSlice(xi_prime_psi_prime_plus_kkp1_xi_psi_over_sxx, k, Nm);
        MatrixXc<Real> tau_n_xi_prime_psi_prime_plus_xi_psi_over_sxx_nnp1 =
            tau_nm * ReduceAndSlice(xi_prime_psi_prime_plus_nnp1_xi_psi_over_sxx, k, Nm);

        L6.col(k_ind) = d_n_xi_prime_psi_prime_nnp1_plus_xi_psi_over_sxx_nnp1_kkp1 * dx_dt_tau_k_sin_t -
            tau_n_xi_prime_psi_prime_plus_xi_psi_over_sxx_nnp1 * dx_dt_d_k_sin_t * k*(k + 1);

        d_n_xi_prime_psi_prime_nnp1_plus_xi_psi_over_sxx_nnp1_kkp1 =
            d_n_times_nnp1 * ReduceAndSlice(psi_prime_psi_prime_plus_kkp1_psi_psi_over_sxx, k, Nm);
        tau_n_xi_prime_psi_prime_plus_xi_psi_over_sxx_nnp1 =
            tau_nm * ReduceAndSlice(psi_prime_psi_prime_plus_nnp1_psi_psi_over_sxx, k, Nm);

        L6P.col(k_ind) = d_n_xi_prime_psi_prime_nnp1_plus_xi_psi_over_sxx_nnp1_kkp1 * dx_dt_tau_k_sin_t -
            tau_n_xi_prime_psi_prime_plus_xi_psi_over_sxx_nnp1 * dx_dt_d_k_sin_t * k*(k + 1);
      }

      ArrayXXc<Real> tmp2 = (n_vec_real * (n_vec_real + 1)).replicate(1, Nm).rowwise() -
      (n_vec_real * (n_vec_real + 1)).transpose();
      ArrayXXc<Real> prefactor1 = (s - 1)*(s + 1) / s * An_Ak(n_vec, n_vec);
      ArrayXXc<Real> prefactor2 = mp_im_unit<Real>() * prefactor1 / tmp2;

      ArrayXXc<Real> Q12 = prefactor1 * K1;
      ArrayXXc<Real> Q21 = -prefactor1 * K2;
      ArrayXXc<Real> Q11 = prefactor2 * L5;
      ArrayXXc<Real> Q22 = prefactor2 * L6;

      ArrayXXc<Real> P12 = prefactor1 * K1P;
      ArrayXXc<Real> P21 = -prefactor1 * K2P;
      ArrayXXc<Real> P11 = prefactor2 * L5P;
      ArrayXXc<Real> P22 = prefactor2 * L6P;

      ArrayXc<Real> prefact_diag1 = -mp_im_unit<Real>()/s * (2*n_vec_real + 1) / (2*n_vec_real*(n_vec_real + 1));
      ArrayXc<Real> prefact_diag2 = -mp_im_unit<Real>()*static_cast<complex<Real>>((s - 1)*(s + 1)/(2*s)) * (2*n_vec_real + 1);
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

      ArrayXi inde = Seq2Array(N_min % 2, Nm - 1, 2);
      ArrayXi indo = Seq2Array(1 - N_min % 2, Nm - 1, 2);

      output[i]->st_4M_Q_eo().M12 = Q12(inde, indo);
      output[i]->st_4M_Q_eo().M21 = Q21(indo, inde);
      output[i]->st_4M_Q_eo().M11 = Q11(inde, inde);
      output[i]->st_4M_Q_eo().M22 = Q22(indo, indo);
      output[i]->st_4M_Q_eo().m = m;
      output[i]->st_4M_Q_eo().ind1 = inde;
      output[i]->st_4M_Q_eo().ind2 = indo;

      output[i]->st_4M_Q_oe().M12 = Q12(indo, inde);
      output[i]->st_4M_Q_oe().M21 = Q21(inde, indo);
      output[i]->st_4M_Q_oe().M11 = Q11(indo, indo);
      output[i]->st_4M_Q_oe().M22 = Q22(inde, inde);
      output[i]->st_4M_Q_oe().m = m;
      output[i]->st_4M_Q_oe().ind1 = indo;
      output[i]->st_4M_Q_oe().ind2 = inde;

      output[i]->st_4M_P_eo().M12 = P12(inde, indo);
      output[i]->st_4M_P_eo().M21 = P21(indo, inde);
      output[i]->st_4M_P_eo().M11 = P11(inde, inde);
      output[i]->st_4M_P_eo().M22 = P22(indo, indo);
      output[i]->st_4M_P_eo().m = m;
      output[i]->st_4M_P_eo().ind1 = inde;
      output[i]->st_4M_P_eo().ind2 = indo;

      output[i]->st_4M_P_oe().M12 = P12(indo, inde);
      output[i]->st_4M_P_oe().M21 = P21(inde, indo);
      output[i]->st_4M_P_oe().M11 = P11(indo, indo);
      output[i]->st_4M_P_oe().M22 = P22(inde, inde);
      output[i]->st_4M_P_oe().m = m;
      output[i]->st_4M_P_oe().ind1 = indo;
      output[i]->st_4M_P_oe().ind2 = inde;

      output[i]->mat_list.push_back("st_4M_Q");
      output[i]->mat_list.push_back("st_4M_P");
    }
    return output;
  }
}

#endif
