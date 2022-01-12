#ifndef RVH_H
#define RVH_H

#include "core.h"
#include "sph.h"
#include "vsh.h"

using namespace Eigen;
using namespace std;

namespace Smarties {
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
    ArrayXc<Real> p_nm;
    ArrayXc<Real> q_nm;
    ArrayXc<Real> a_nm;
    ArrayXc<Real> b_nm;
    ArrayXc<Real> c_nm;
    ArrayXc<Real> d_nm;
  };

  template <class Real>
  struct stCrossSection {
    ArrayXr<Real> ext;
    ArrayXr<Real> sca;
    ArrayXr<Real> abs;
  };

  template <class Real>
  vector<unique_ptr<stTR<Real>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<Real>>>& st_PQ_list, bool get_R = false);

  template <class Real>
  vector<unique_ptr<stTR<Real>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<Real>>>& st_mat_list, int N_max);

  // Untested
  template <class Real>
  vector<unique_ptr<stTR<Real>>> rvhGetSymmetricMat(const vector<unique_ptr<stTR<Real>>>& st_mat_list, vector<string> mat_list = {"st_4M_T"});

  template <class Real>
  unique_ptr<stAbcdnm<Real>> rvhGetFieldCoefficients(int N_max, const vector<unique_ptr<stTR<Real>>>& st_TR_list,
      const unique_ptr<stIncPar<Real>>& st_inc_par, unique_ptr<stIncEabnm<Real>> st_inc_E_abnm = unique_ptr<stIncEabnm<Real>>());

  template <class Real>
  unique_ptr<stCrossSection<Real>> rvhGetAverageCrossSections(
      const ArrayXr<Real>& k1, const vector<vector<unique_ptr<stTR<Real>>>>& st_TR_list);
}

#endif
