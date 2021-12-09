#ifndef LOW_LEVEL_H
#define LOW_LEVEL_H

#include "raman_elastic_scattering.h"
#include "high_level.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  enum scheme {GAUSS, RECTANGLE, PTS};

  template <typename Real>
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

    static stRtfunc* auxPrepareIntegrals(size_t nNint, scheme scheme_type);

    static stIncEabnm* vshGetIncidentCoeffs(int n_max, typename HighLevel<Real>::stIncPar& angles);

    static stPinmTaunm* vshPinmTaunm(size_t n_ax, const ArrayXr<Real>& theta);

    static ArrayXXc<Real>* vshRBchi(const RowArrayXr<Real> n, const ArrayXr<Real>& x);

    static ArrayXXc<Real>* vshRBpsi(const RowArrayXr<Real> n, const ArrayXr<Real>& x);

  private:
    LowLevel() {};
  };
}

#endif
