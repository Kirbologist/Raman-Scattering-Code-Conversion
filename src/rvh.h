#ifndef RVH_H
#define RVH_H

#include "raman_elastic_scattering.h"
#include "sph.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  struct stTR {
    st4M<Real>* st_4M_T_eo;
    st4M<Real>* st_4M_T_oe;
    st4M<Real>* st_4M_R_eo;
    st4M<Real>* st_4M_R_oe;
    vector<string> mat_list;

    stTR() {
      this.st_4M_T_eo = new st4M<Real>();
      this.st_4M_T_oe = new st4M<Real>();
      this.st_4M_R_eo = new st4M<Real>();
      this.st_4M_R_oe = new st4M<Real>();
    }
  };

  template <class Real>
  vector<stTR<Real>>* rvhGetTRfromPQ(const vector<stPQ<Real>>*, bool get_r = false);
}

#endif
