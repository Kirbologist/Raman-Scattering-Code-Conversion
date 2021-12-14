#ifndef LOW_LEVEL_H
#define LOW_LEVEL_H

#include "raman_elastic_scattering.h"
#include "high_level.h"
#include <Eigen/CXX11/Tensor>

using namespace Eigen;
using namespace std;

namespace Raman {
  enum sInt {GAUSS, RECTANGLE, PTS};
  enum sBessel {J, H1};

  template <class Real>
  class LowLevel {
  public:
    struct stGLQuad {
      ArrayXr<Real> x;
      ArrayXr<Real> w;
    };

    struct stRtfunc {
      size_t nNbTheta;
      ArrayXr<Real> theta;
      ArrayXr<Real> wTheta;
    };

    struct stFpovx {
      ArrayXXc<Real>* Fpovx; // array of Arrays
      ArrayXXc<Real>* rb_chi;
      ArrayXXc<Real>* rb_psi;
    };

    struct stFprow {
      ArrayXXr<Real> S;
      ArrayXXr<Real> loss_prec_S;
    };

    struct stBesselPrimes {
      ArrayXXc<Real>* xi_psi; // array of Arrays
      ArrayXXc<Real>* xi_prime_psi; // array of Arrays
      ArrayXXc<Real>* xi_psi_prime; // array of Arrays
      ArrayXXc<Real>* xi_prime_psi_prime; // array of Arrays
      ArrayXXc<Real>* xi_psi_over_sxx; // array of Arrays
      ArrayXXc<Real>* xi_prime_psi_prime_plus_nnp1_xi_psi_over_ssx; // array of Arrays
      ArrayXXc<Real>* xi_prime_psi_prime_plus_kkp1_xi_psi_over_ssx; // array of Arrays
    };

    struct stXiPsiAll {
      stBesselPrimes* xi_primes;
      ArrayXXc<Real> for_diag_Lt1;
      ArrayXXc<Real> for_diag_Lt2;
      ArrayXXc<Real> for_diag_Lt3;
    };

    struct stPsiPsiAll {
      stBesselPrimes* psi_primes;
      ArrayXXc<Real> for_diag_Lt1;
      ArrayXXc<Real> for_diag_Lt2;
      ArrayXXc<Real> for_diag_Lt3;
    };

    struct stBesselProducts {
      stXiPsiAll* xi_struct;
      stPsiPsiAll* psi_struct;
    };

    struct stBessel {
      ArrayXXc<Real>* xi_psi; // array of Arrays
      ArrayXXc<Real>* psi_psi; // array of Arrays
      ArrayXXc<Real>* chi_n;
      ArrayXXc<Real>* psi_n;
      ArrayXXc<Real>* psi_k;
    };

    struct stEAllPhi {
      RowArrayXr<Real> theta;
      RowArrayXr<Real> r_of_theta;
      ArrayXXc<Real>* CErm; // array of Arrays
      ArrayXXc<Real>* CEtm; // array of Arrays
      ArrayXXc<Real>* CEfm; // array of Arrays
    };

    struct stZnAll {
      ArrayXXc<Real> Z0;
      ArrayXXc<Real> Z1;
      ArrayXXc<Real> Z2;
    };

    struct stIncEabnm {
      ArrayXc<Real> a_nm;
      ArrayXc<Real> b_nm;
    };

    struct stPinmTaunm {
      ArrayXXr<Real> pi_nm;
      ArrayXXr<Real> tau_nm;
      ArrayXXr<Real> p_n0;
    };

    static void destructStBesselProducts(stBesselProducts* st) {
      delete[] st->xi_struct->xi_primes->xi_psi;
      delete[] st->xi_struct->xi_primes->xi_prime_psi;
      delete[] st->xi_struct->xi_primes->xi_psi_prime;
      delete[] st->xi_struct->xi_primes->xi_prime_psi_prime;
      delete[] st->xi_struct->xi_primes->xi_psi_over_sxx;
      delete[] st->xi_struct->xi_primes->xi_prime_psi_prime_plus_nnp1_xi_psi_over_ssx;
      delete[] st->xi_struct->xi_primes->xi_prime_psi_prime_plus_kkp1_xi_psi_over_ssx;
      delete st->xi_struct->xi_primes;
      delete st->xi_struct;


      delete[] st->psi_struct->psi_primes->xi_psi;
      delete[] st->psi_struct->psi_primes->xi_prime_psi;
      delete[] st->psi_struct->psi_primes->xi_psi_prime;
      delete[] st->psi_struct->psi_primes->xi_prime_psi_prime;
      delete[] st->psi_struct->psi_primes->xi_psi_over_sxx;
      delete[] st->psi_struct->psi_primes->xi_prime_psi_prime_plus_nnp1_xi_psi_over_ssx;
      delete[] st->psi_struct->psi_primes->xi_prime_psi_prime_plus_kkp1_xi_psi_over_ssx;
      delete st->psi_struct->psi_primes;
      delete st->psi_struct;

      delete st;
    }

    static stGLQuad* auxInitLegendreQuad(size_t N1, Real a = -1.0, Real b = 1.0);

    static stRtfunc* auxPrepareIntegrals(size_t nNint, sInt type);

    static stFpovx* sphGetFpovx(int n_n_max, Real s, ArrayXr<Real>& x);

    static stFprow* sphGetFpRow(int n, Real s, ArrayXr<Real>& x);

    static stBesselProducts* sphGetModifiedBesselProducts(int n_n_max, Real s, ArrayXr<Real>& x, int N_B);

    static stBessel* sphGetXiPsi(int n_n_max, Real s, ArrayXr<Real>& x, int N_B);

    // Untested
    static stEAllPhi* vshEgenThetaAllPhi(ArrayXr<Real>& lambda,
        ArrayXr<Real>& epsilon, ArrayXXc<Real>& p_nm, ArrayXXc<Real>& q_nm,
        RowArrayXr<Real>& rt, RowArrayXr<Real>& theta, sBessel type, stPinmTaunm* stPT = nullptr);

    static stIncEabnm* vshGetIncidentCoeffs(int n_max, typename HighLevel<Real>::stIncPar* angles);

    static stPinmTaunm* vshPinmTaunm(size_t n_ax, const ArrayXr<Real>& theta);

    static ArrayXXc<Real>* vshRBchi(ArrayXr<Real> n, const ArrayXr<Real>& x);

    static ArrayXXc<Real>* vshRBpsi(ArrayXr<Real> n, const ArrayXr<Real>& x);

  private:
    LowLevel() {};

    static stBesselPrimes* sphGetBesselProductsPrimes(ArrayXXc<Real>* prods, int N);

    static ArrayXXr<Real>* sphGetUforFp(int n);

    // Untested
    static stZnAll* vshGetZnAll(size_t n_n_max, ArrayXr<Real>& rho, sBessel type);
  };
}

#endif
