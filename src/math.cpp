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
  RowArrayXr<Real> arr_bessel_j(RowArrayXr<Real>& nu, Real x) {
    RowArrayXr<Real> output(1, nu.size());
    for (int i = 0; i < nu.size(); i++)
      output(i) = cyl_bessel_j(nu(i), x);
    return output;
  }

  template<class Real>
  RowArrayXr<Real> arr_bessel_y(RowArrayXr<Real>& nu, Real x) {
    RowArrayXr<Real> output(1, nu.size());
    for (int i = 0; i < nu.size(); i++)
      output(i) = cyl_neumann(nu(i), x);
    return output;
  }

  template RowArrayXr<double> arr_bessel_j(RowArrayXr<double>& nu, double x);
  template RowArrayXr<double> arr_bessel_y(RowArrayXr<double>& nu, double x);
}
