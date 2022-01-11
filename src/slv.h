#ifndef SLV_H
#define SLV_H

#include "core.h"
#include "rvh.h"
#include "vsh.h"
#include "sph.h"

using namespace std;

namespace Raman {
  template <class Real>
  struct stTmatrix {
    unique_ptr<stCrossSection<Real>> st_coa;
    vector<unique_ptr<stTR<Real>>> st_TR_list;
  };

  template <class Real>
  unique_ptr<stTmatrix<Real>> slvForT(const unique_ptr<stParams<Real>>& params,
      const unique_ptr<stOptions<Real>>& options, unique_ptr<stRtfunc<Real>> stGeometry = unique_ptr<stRtfunc<Real>>());
}

#endif
