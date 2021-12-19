#ifndef VSH_H
#define VSH_H

#include "raman_elastic_scattering.h"

using namespace Eigen;

namespace Raman {
  enum sBessel {J, H1};

  template <class Real>
  struct stPinmTaunm {
    ArrayXXr<Real> pi_nm;
    ArrayXXr<Real> tau_nm;
    ArrayXXr<Real> p_n0;
  };

  template <class Real>
  struct stIncEabnm {
    ArrayXc<Real> a_nm;
    ArrayXc<Real> b_nm;
  };

  template <class Real>
  struct stZnAll {
    ArrayXXc<Real> Z0;
    ArrayXXc<Real> Z1;
    ArrayXXc<Real> Z2;
  };

  template <class Real>
  struct stEAllPhi {
    RowArrayXr<Real> theta;
    RowArrayXr<Real> r_of_theta;
    ArrayXXc<Real>* Erm;
    ArrayXXc<Real>* Etm;
    ArrayXXc<Real>* Efm;
    int N_max;
  };

  template <class Real>
  struct stEforPhi {
    ArrayXXc<Real> Er;
    ArrayXXc<Real> Et;
    ArrayXXc<Real> Ef;
    ArrayXr<Real> theta;
    Real phi0;
  };

  // Untested
  template <class Real>
  stIncPar<Real>* vshMakeIncidentParams(sIncType type, size_t N_max);

  // Untested
  template <class Real>
  stIncPar<Real>* vshMakeIncidentParams(sIncType type, size_t N_max,
      Real theta_p, Real phi_p, Real alpha_p);

  template <class Real>
  stPinmTaunm<Real>* vshPinmTaunm(size_t N_max, const ArrayXr<Real>& theta);

  template <class Real>
  stIncEabnm<Real>* vshGetIncidentCoeffs(int N_max, stIncPar<Real>* angles);

  // Untested
  template <class Real>
  stZnAll<Real>* vshGetZnAll(size_t N_max, const ArrayXr<Real>& rho, sBessel type);

  // Untested
  template <class Real>
  stEAllPhi<Real>* vshEgenThetaAllPhi(
      const ArrayXr<Real>& lambda, const ArrayXr<Real>& epsilon, const ArrayXXc<Real>& p_nm,
      const ArrayXXc<Real>& q_nm, const RowArrayXr<Real>& rt, const RowArrayXr<Real>& theta,
      sBessel type, stPinmTaunm<Real>* stPT = nullptr);

  template <class Real>
  stEforPhi<Real>* vshEthetaForPhi(stEAllPhi<Real>* stEsurf, Real phi0);

  template <class Real>
  ArrayXXc<Real> vshRBchi(ArrayXr<Real> n, const ArrayXr<Real>& x);

  template <class Real>
  ArrayXXc<Real> vshRBpsi(ArrayXr<Real> n, const ArrayXr<Real>& x);
}

#endif
