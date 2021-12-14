#include "high_level.h"
#include "math.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  typename LowLevel<Real>::stRtfunc* HighLevel<Real>::sphMakeGeometry(
      size_t n_Nb_theta, Real a, Real c, ArrayXr<Real>* theta) {
    sInt type;
    stRtfunc* output;
    if (theta != nullptr) {
      output->w_theta = ArrayXr<Real>::Zero(theta->size());
      output->theta = *theta;
      output->n_Nb_theta = theta->size();
      type = PTS;
    } else {
      type = GAUSS;
      output = auxPrepareIntegrals(2*n_Nb_theta, type);
      output->theta = output->theta(seq(0, n_Nb_theta - 1));
      output->w_theta = output->w_theta(seq(0, n_Nb_theta - 1))*2;
      output->n_Nb_theta = n_Nb_theta;
    }

    output->a = a;
    output->c = c;
    output->type = type;

    ArrayXr<Real> sin_t = sin(output->theta), cos_t = cos(output->theta);

    output->r = a*c/sqrt(pow(c*sin_t, 2) + pow(a*cos_t, 2));
    output->dr_dt = (pow(a, 2) - pow(c, 2))/pow(a*c, 2)*sin_t * cos_t * output->r.pow(3);

    output->h = max(a, c)/min(a, c);
    output->r0 = pow(pow(a, 2)*c, static_cast<Real>(1)/3);

    return output;
  }

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
