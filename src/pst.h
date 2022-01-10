#ifndef PST_H
#define PST_H

#include "core.h"
#include "rvh.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  struct stRes {
    ArrayXc<Real> p_nm;
    ArrayXc<Real> q_nm;
    ArrayXc<Real> a_nm;
    ArrayXc<Real> b_nm;
    ArrayXc<Real> c_nm;
    ArrayXc<Real> d_nm;
    int N_max;
    ArrayXr<Real> lambda;
    Real epsilon1;
    ArrayXr<Real> epsilon2;
    unique_ptr<stIncPar<Real>> inc_par;
    Real a;
    Real c;

    stRes() {};

    stRes(const stAbcdnm<Real>& base) {
      this->p_nm = base.p_nm;
      this->q_nm = base.q_nm;
      this->a_nm = base.a_nm;
      this->b_nm = base.b_nm;
      this->c_nm = base.c_nm;
      this->d_nm = base.d_nm;
    }
  };

  template <class Real>
  struct stSM {
    int L_max;
    ArrayXr<Real> ALF1n;
    ArrayXr<Real> ALF2n;
    ArrayXr<Real> ALF3n;
    ArrayXr<Real> ALF4n;
    ArrayXr<Real> BET1n;
    ArrayXr<Real> BET2n;
    Real asym_par;
    ArrayXr<Real> theta;
    ArrayXr<Real> theta_deg;
    ArrayXr<Real> F11;
    ArrayXr<Real> F22;
    ArrayXr<Real> F33;
    ArrayXr<Real> F44;
    ArrayXr<Real> F12;
    ArrayXr<Real> F34;
    Array<Real, Dynamic, 7> all_AB;
    Array<Real, Dynamic, 7> all_F;
  };

  // Untested
  template <class Real>
  unique_ptr<stRes<Real>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<Real>>& st_abcdnm, int N_max, ArrayXr<Real> lambda,
    ArrayXr<Real> epsilon2, Real epsilon1, unique_ptr<stIncPar<Real>> inc_par,
    Real a = numeric_limits<Real>::quiet_NaN(), Real c = numeric_limits<Real>::quiet_NaN());

  // Untested
  template <class Real>
  unique_ptr<stRes<Real>> pstMakeStructForField(
      const unique_ptr<stAbcdnm<Real>>& st_abcdnm, const unique_ptr<stParams<Real>>& params);

  // Untested
  template <class Real>
  unique_ptr<stSM<Real>> pstScatteringMatrixOA(const vector<unique_ptr<stTR<Real>>>& st_TR_list,
      Real lambda, Real sca, int Nb_theta = 2);
}

#endif
