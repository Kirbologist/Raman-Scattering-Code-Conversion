#ifndef RVH_H
#define RVH_H

#include "raman_elastic_scattering.h"
#include "sph.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  struct stTR : stMat <Real> {
    inline st4M<Real>& st_4M_T_eo() { return this->st_4M_list[0]; }
    inline st4M<Real>& st_4M_T_oe() { return this->st_4M_list[1]; }
    inline st4M<Real>& st_4M_R_eo() { return this->st_4M_list[2]; }
    inline st4M<Real>& st_4M_R_oe() { return this->st_4M_list[3]; }

    stTR() {
      this->st_4M_T_eo() = st4M<Real>();
      this->st_4M_T_oe() = st4M<Real>();
      this->st_4M_R_eo() = st4M<Real>();
      this->st_4M_R_oe() = st4M<Real>();
    }
  };

  // Untested
  template <class Real>
  vector<unique_ptr<stTR<Real>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<Real>>>& st_PQ_list, bool get_r = false);

  template <class Real>
  vector<unique_ptr<stPQ<Real>>>& rvhTruncateMatrices(vector<unique_ptr<stPQ<Real>>>& st_mat_list, int N_max);
}

#endif
