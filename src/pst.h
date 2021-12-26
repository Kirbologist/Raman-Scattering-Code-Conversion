#ifndef PST_H
#define PST_H

#include "raman_elastic_scattering.h"
#include "rvh.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  struct stRes {
    unique_ptr<stAbcdnm<Real>> st_abcdnm;
    int N_max;
    ArrayXr<Real> lambda;
    Real epsilon1;
    ArrayXc<Real> epsilon2;
    unique_ptr<stIncPar<Real>> inc_par;
    Real a;
    Real c;
  };

  // Untested
  template <class Real>
  unique_ptr<stRes<Real>> pstMakeStructForField(
    unique_ptr<stAbcdnm<Real>> st_abcdnm, int N_max, ArrayXr<Real> lambda,
    ArrayXr<Real> epsilon2, Real epsilon1, unique_ptr<stIncPar<Real>> inc_par,
    Real a = numeric_limits<Real>::quiet_NaN(), Real c = numeric_limits<Real>::quiet_NaN());

  // Untested
  template <class Real>
  unique_ptr<stRes<Real>> pstMakeStructForField(
      unique_ptr<stAbcdnm<Real>> st_abcdnm, unique_ptr<stParams<Real>> params);
}

#endif
