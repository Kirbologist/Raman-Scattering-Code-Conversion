#ifndef HIGH_LEVEL_H
#define HIGH_LEVEL_H

#include "raman_elastic_scattering.h"
#include "low_level.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  enum sIncType {KxEz, KxEy, KyEx, KyEz, KzEx, KzEy, GENERAL};

  template <class Real>
  class HighLevel {
  public:
    struct stIncPar {
      sIncType type;
      Real theta_p;
      Real phi_p;
      Real alpha_p;
      ArrayXi abs_m_vec;
    };

    static LowLevel<Real>::stRtfunc* sphMakeGeometry(size_t n_Nb_theta, Real a, Real c, ArrayXr<Real>* theta = nullptr);

    // Untested
    static stIncPar* vshMakeIncidentParams(sIncType type, size_t n_max);

    // Untested
    static stIncPar* vshMakeIncidentParams(sIncType type, size_t n_max,
        Real theta_p, Real phi_p, Real alpha_p);

  private:
    HighLevel() {};
  };
}

#endif
