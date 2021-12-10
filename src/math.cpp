#include "math.h"
#include "raman_elastic_scattering.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template<class Real>
  Real mp_pi(void) {
    return 3.14159265359;
  }

  template<class Real>
  Real mp_eps(void) {
    return numeric_limits<Real>::epsilon();
  }

  template<class Real>
  complex<Real> mp_im_unit(void) {
    complex<Real> i(0.0, 1.0);
    return i;
  }

  template<class Real>
  ArrayXr<Real> arr_bessel_j(ArrayXr<Real>& nu, Real x) {
    ArrayXr<Real> output(1, nu.size());
    for (int i = 0; i < nu.size(); i++)
      output(i) = cyl_bessel_j(nu(i), x);
    return output;
  }

  template<class Real>
  ArrayXr<Real> arr_bessel_y(ArrayXr<Real>& nu, Real x) {
    ArrayXr<Real> output(1, nu.size());
    for (int i = 0; i < nu.size(); i++)
      output(i) = cyl_neumann(nu(i), x);
    return output;
  }

  template ArrayXr<double> arr_bessel_j(ArrayXr<double>& nu, double x);
  template ArrayXr<double> arr_bessel_y(ArrayXr<double>& nu, double x);
}
