#ifndef RVH_H
#define RVH_H

#include "raman_elastic_scattering.h"
#include "sph.h"
#include "vsh.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  struct stTR : stMat<Real> {
    using stMat<Real>::stMat;
    inline st4M<Real>& st_4M_T_eo() { return this->st_4M_list[0]; }
    inline st4M<Real>& st_4M_T_oe() { return this->st_4M_list[1]; }
    inline st4M<Real>& st_4M_R_eo() { return this->st_4M_list[2]; }
    inline st4M<Real>& st_4M_R_oe() { return this->st_4M_list[3]; }
  };

  template <class Real>
  struct stAbcdnm {
    ArrayXc<Real> pnm;
    ArrayXc<Real> qnm;
    ArrayXc<Real> anm;
    ArrayXc<Real> bnm;
    ArrayXc<Real> cnm;
    ArrayXc<Real> dnm;
  };

  // Untested
  template <class Real>
  vector<unique_ptr<stTR<Real>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<Real>>>& st_PQ_list, bool get_R = false);

  template <class Real>
  vector<unique_ptr<stPQ<Real>>> rvhTruncateMatrices(const vector<unique_ptr<stPQ<Real>>>& st_mat_list, int N_max);

  // Untested
  template <class Real>
  vector<unique_ptr<stTR<Real>>> rvhGetSymmetricMat(const vector<unique_ptr<stTR<Real>>>& st_mat_list, vector<string> mat_list = {"st_4M_T"});

  template <class Real>
  unique_ptr<stAbcdnm<Real>> rvhGetFieldCoefficients(int N_max, const vector<unique_ptr<stTR<Real>>>& st_TR_list,
      const unique_ptr<stIncPar<Real>>& st_inc_par, stIncEabnm<Real>* st_inc_E_abnm = nullptr);
}

#endif
