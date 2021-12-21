#ifndef RVH_H
#define RVH_H

#include "raman_elastic_scattering.h"
#include "sph.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  struct stTR {
    unique_ptr<st4M<Real>> st_4M_T_eo;
    unique_ptr<st4M<Real>> st_4M_T_oe;
    unique_ptr<st4M<Real>> st_4M_R_eo;
    unique_ptr<st4M<Real>> st_4M_R_oe;
    vector<string> mat_list;

    stTR() {
      this->st_4M_T_eo = make_unique<st4M<Real>>();
      this->st_4M_T_oe = make_unique<st4M<Real>>();
      this->st_4M_R_eo = make_unique<st4M<Real>>();
      this->st_4M_R_oe = make_unique<st4M<Real>>();
    }
  };

  // Untested
  template <class Real>
  unique_ptr<vector<stTR<Real>>> rvhGetTRfromPQ(const unique_ptr<vector<stPQ<Real>>>& st_PQ_list, bool get_r = false);
}

#endif
