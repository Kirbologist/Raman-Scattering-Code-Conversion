#ifndef SPH_H
#define SPH_H

#include "raman_elastic_scattering.h"
#include "raman_aux.h"

using namespace Eigen;

namespace Raman {

  template <class Real>
  struct stFprow {
    ArrayXXr<Real> S;
    ArrayXXr<Real> loss_prec_S;
  };

  template <class Real>
  struct stFpovx {
    ArrayXXc<Real>* Fpovx; // array of Arrays
    ArrayXXc<Real>* rb_chi;
    ArrayXXc<Real>* rb_psi;
  };

  template <class Real>
  struct stBessel {
    ArrayXXc<Real>* xi_psi; // array of Arrays
    ArrayXXc<Real>* psi_psi; // array of Arrays
    ArrayXXc<Real>* chi_n;
    ArrayXXc<Real>* psi_n;
    ArrayXXc<Real>* psi_k;
  };

  template <class Real>
  struct stBesselPrimes {
    ArrayXXc<Real>* xi_psi; // array of Arrays
    ArrayXXc<Real>* xi_prime_psi; // array of Arrays
    ArrayXXc<Real>* xi_psi_prime; // array of Arrays
    ArrayXXc<Real>* xi_prime_psi_prime; // array of Arrays
    ArrayXXc<Real>* xi_psi_over_sxx; // array of Arrays
    ArrayXXc<Real>* xi_prime_psi_prime_plus_nnp1_xi_psi_over_ssx; // array of Arrays
    ArrayXXc<Real>* xi_prime_psi_prime_plus_kkp1_xi_psi_over_ssx; // array of Arrays
    ArrayXXc<Real> for_diag_Lt1;
    ArrayXXc<Real> for_diag_Lt2;
    ArrayXXc<Real> for_diag_Lt3;
  };

  template <class Real>
  struct stBesselProducts {
    stBesselPrimes<Real>* st_xi_psi_all;
    stBesselPrimes<Real>* st_psi_psi_all;
  };

  template <class Real>
  void destructStBesselProducts(stBesselProducts<Real>* st) {
    delete[] st->st_xi_psi_all->xi_psi;
    delete[] st->st_xi_psi_all->xi_prime_psi;
    delete[] st->st_xi_psi_all->xi_psi_prime;
    delete[] st->st_xi_psi_all->xi_prime_psi_prime;
    delete[] st->st_xi_psi_all->xi_psi_over_sxx;
    delete[] st->st_xi_psi_all->xi_prime_psi_prime_plus_nnp1_xi_psi_over_ssx;
    delete[] st->st_xi_psi_all->xi_prime_psi_prime_plus_kkp1_xi_psi_over_ssx;
    delete st->st_xi_psi_all;

    delete[] st->st_psi_psi_all->xi_psi;
    delete[] st->st_psi_psi_all->xi_prime_psi;
    delete[] st->st_psi_psi_all->xi_psi_prime;
    delete[] st->st_psi_psi_all->xi_prime_psi_prime;
    delete[] st->st_psi_psi_all->xi_psi_over_sxx;
    delete[] st->st_psi_psi_all->xi_prime_psi_prime_plus_nnp1_xi_psi_over_ssx;
    delete[] st->st_psi_psi_all->xi_prime_psi_prime_plus_kkp1_xi_psi_over_ssx;
    delete st->st_psi_psi_all;

    delete st;
  }

  template <class Real>
  ArrayXXr<Real>* sphGetUforFp(int n);

  template <class Real>
  stFprow<Real>* sphGetFpRow(int n, Real s, ArrayXr<Real>& x);

  template <class Real>
  stFpovx<Real>* sphGetFpovx(int n_n_max, Real s, ArrayXr<Real>& x);

  template <class Real>
  stBessel<Real>* sphGetXiPsi(int n_n_max, Real s, ArrayXr<Real>& x, int N_B);

  template <class Real>
  stBesselPrimes<Real>* sphGetBesselProductsPrimes(ArrayXXc<Real>* prods, int N);

  template <class Real>
  stBesselProducts<Real>* sphGetModifiedBesselProducts(int n_n_max, Real s, ArrayXr<Real>& x, int N_B);

  template <class Real>
  stRtfunc<Real>* sphMakeGeometry(size_t n_Nb_theta, Real a, Real c, ArrayXr<Real>* theta = nullptr);
}

#endif
