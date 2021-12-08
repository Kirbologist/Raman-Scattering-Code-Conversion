#ifndef LOW_LEVEL_H
#define LOW_LEVEL_H

#include "raman_elastic_scattering.h"
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

namespace Raman {
  enum scheme {GAUSS, RECTANGLE, PTS};

  template <typename Real>
  class LowLevel;

  template <typename Real>
  class LowLevel {
  public:
    struct legendreQuad {
      ArrayXr<Real> x;
      ArrayXr<Real> w;
    };

    struct stRtfunc {
      size_t nNbTheta;
      ArrayXr<Real> theta;
      ArrayXr<Real> wTheta;
    };

    struct stPinmTaunm {
      ArrayXXr<Real> pi_nm;
      ArrayXXr<Real> tau_nm;
      ArrayXXr<Real> p_n0;
    }

    static legendreQuad* auxInitLegendreQuad(size_t N1, Real a = -1.0, Real b = 1.0);

    static stRtfunc* auxPrepareIntegrals(size_t nNint, scheme scheme_type);

    static stPinmTaunm* vshPinmTaunm(size_t n_ax, const ArrayXr<Real>& theta);

    static ArrayXXc<Real>* vshRBchi(const RowArrayXr<Real> n, const ArrayXr<Real>& x);

    static ArrayXXc<Real>* vshRBpsi(const RowArrayXr<Real> n, const ArrayXr<Real>& x);

  private:
    LowLevel() {};
  };
};

#endif
