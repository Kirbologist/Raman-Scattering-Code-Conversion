#ifndef AUX_H
#define AUX_H

#include "raman_elastic_scattering.h"

using namespace Eigen;

namespace Raman {
  enum sInt {GAUSS, RECTANGLE, PTS};

  template <class Real>
  struct stGLQuad {
    ArrayXr<Real> x;
    ArrayXr<Real> w;
  };

  template <class Real>
  struct stRtfunc {
    size_t n_Nb_theta;
    ArrayXr<Real> theta;
    ArrayXr<Real> w_theta;
    ArrayXr<Real> r;
    ArrayXr<Real> dr_dt;
    Real a;
    Real c;
    Real h;
    Real r0;
    sInt type;
  };

  template <class Real>
  stGLQuad<Real>* auxInitLegendreQuad(size_t N1, Real a = -1.0, Real b = 1.0);

  template <class Real>
  stRtfunc<Real>* auxPrepareIntegrals(size_t nNint, sInt type);
}

#endif