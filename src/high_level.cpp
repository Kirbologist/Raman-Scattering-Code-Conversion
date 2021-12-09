#include "high_level.h"
#include "math.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  typename HighLevel<Real>::stIncPar* HighLevel<Real>::vshMakeIncidentParams(
      sIncType type, size_t n_max) {
    stIncPar* output = new stIncPar();
    switch (type) {
      case KxEz: {
        output->type = KxEz;
        output->theta_p = PI/2;
        output->phi_p = 0;
        output->alpha_p = PI;
        output->abs_m_vec = ArrayXi::LinSpaced(n_max + 1, 0, n_max);
        break;
      }

      case KxEy: {
        output->type = KxEy;
        output->theta_p = PI/2;
        output->phi_p = 0;
        output->alpha_p = PI/2;
        output->abs_m_vec = ArrayXi::LinSpaced(n_max + 1, 0, n_max);
        break;
      }

      case KyEz: {
        output->type = KyEz;
        output->theta_p = PI/2;
        output->phi_p = PI/2;
        output->alpha_p = PI;
        output->abs_m_vec = ArrayXi::LinSpaced(n_max + 1, 0, n_max);
        break;
      }

      case KyEx: {
        output->type = KyEx;
        output->theta_p = PI/2;
        output->phi_p = PI/2;
        output->alpha_p = -PI/2;
        output->abs_m_vec = ArrayXi::LinSpaced(n_max + 1, 0, n_max);
        break;
      }

      case KzEx: {
        output->type = KzEx;
        output->theta_p = 0;
        output->phi_p = 0;
        output->alpha_p = 0;
        output->abs_m_vec = {1};
        break;
      }

      case KzEy: {
        output->type = KzEy;
        output->theta_p = 0;
        output->phi_p = 0;
        output->alpha_p = PI/2;
        output->abs_m_vec = {1};
        break;
      }

      default:
        output->type = KzEx;
        output->theta_p = 0;
        output->phi_p = 0;
        output->alpha_p = 0;
        output->abs_m_vec = {1};
    }

    return output;
  }

  template <class Real>
  typename HighLevel<Real>::stIncPar* HighLevel<Real>::vshMakeIncidentParams(
      sIncType type, size_t n_max, Real theta_p, Real phi_p, Real alpha_p) {
    stIncPar* output = new stIncPar();
    if (type == GENERAL) {
      output->type = GENERAL;
      output->theta_p = theta_p;
      output->phi_p = phi_p;
      output->alpha_p = alpha_p;
      output->abs_m_vec = ArrayXi::LinSpaced(n_max + 1, 0, n_max);
    } else {
      output = HighLevel<Real>::vshMakeIncidentParams(type, n_max);
    }
    return output;
  }

  template class HighLevel<double>;
}
