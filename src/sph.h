#ifndef SPH_H
#define SPH_H

#include "raman_elastic_scattering.h"
#include "raman_aux.h"

using namespace Eigen;

namespace Raman {

  // Might be worth re-testing
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
  unique_ptr<ArrayXXr<Real>> sphGetUforFp(int n);

  template <class Real>
  unique_ptr<stFprow<Real>> sphGetFpRow(int n, Real s, const ArrayXr<Real>& x);

  template <class Real>
  unique_ptr<stFpovx<Real>> sphGetFpovx(int N_max, Real s, const ArrayXr<Real>& x);

  template <class Real>
  unique_ptr<stBessel<Real>> sphGetXiPsi(int N_max, Real s, const ArrayXr<Real>& x, int NB);

  template <class Real>
  unique_ptr<stBesselPrimes<Real>> sphGetBesselProductsPrimes(const Tensor3c<Real>& prods);

  template <class Real>
  unique_ptr<stBesselProducts<Real>> sphGetModifiedBesselProducts(int N_max, Real s, const ArrayXr<Real>& x, int NB);

  template <class Real>
  unique_ptr<stRtfunc<Real>> sphMakeGeometry(size_t Nb_theta, Real a, Real c, const unique_ptr<ArrayXr<Real>> theta = unique_ptr<ArrayXr<Real>>());

  // Could use some more thorough testing
  template <class Real>
  size_t sphCheckBesselConvergence(size_t N_req, Real s, const ArrayXr<Real>& x, Real acc, size_t N_min);

  // Could use some more thorough testing
  template <class Real>
  size_t sphEstimateNB(size_t NQ, const unique_ptr<stRtfunc<Real>>& stGeometry, const unique_ptr<stParams<Real>>& params, Real acc = 1e-13);

  template <class Real>
  vector<unique_ptr<stPQ<Real>>> sphCalculatePQ(int N_max, const ArrayXi& abs_m_vec, const unique_ptr<stRtfunc<Real>>& Rt_func, const unique_ptr<stParams<Real>>& params, int NB = -1);
}

#endif
