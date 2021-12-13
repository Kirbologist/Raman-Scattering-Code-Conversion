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
      ArrayXXc<Real>* Fpovx;
      ArrayXXc<Real>* rb_chi;
      ArrayXXc<Real>* rb_psi;
    };

    struct stFprow {
      ArrayXXr<Real> S;
      ArrayXXr<Real> loss_prec_S;
    };

    struct stEAllPhi {
      RowArrayXr<Real> theta;
      RowArrayXr<Real> r_of_theta;
      ArrayXXc<Real>* CErm; // Why not just make this a tensor
      ArrayXXc<Real>* CEtm;
      ArrayXXc<Real>* CEfm;
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

    static stGLQuad* auxInitLegendreQuad(size_t N1, Real a = -1.0, Real b = 1.0);

    static stRtfunc* auxPrepareIntegrals(size_t nNint, sInt type);

    static stFpovx* sphGetFpovx(int n_n_max, Real s, ArrayXr<Real>& x);

    static stFprow* sphGetFpRow(int n, Real s, ArrayXr<Real>& x);

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

    // Untested
    static stZnAll* vshGetZnAll(size_t n_n_max, ArrayXr<Real>& rho, sBessel type);

    static ArrayXXr<Real>* sphGetUforFp(int n);
  };
}

#endif
