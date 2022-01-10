#include "smarties_aux.h"
#include "math.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  unique_ptr<stGLQuad<Real>> auxInitLegendreQuad(size_t N1, Real a, Real b) {
    int N = N1 - 1, N2 = N1 + 1;
    ArrayXr<Real> xu, y, y0(N1), Lp;
    ArrayXXr<Real> L(N1, N2);

    xu = ArrayXr<Real>::LinSpaced(N1, -1, 1);
    y = cos((2*ArrayXr<Real>::LinSpaced(N1, 0, N) + 1)*mp_pi<Real>()/(2*N + 2)) +
        (0.27/N1)*sin(mp_pi<Real>()*xu*N/N2);
    y0.fill(2);
    L.col(0).setOnes();

    int n_iter = 0;
    while ((n_iter < 15) && mp_eps<Real>() <
        abs(acos(y) - acos(y0.template cast<complex<Real>>())).maxCoeff()) {
      n_iter++;
      L.col(1) = y;

      for (size_t i = 1; i < N1; i++)
        L.col(i + 1) = ((2*i + 1)*y * L.col(i) - i*L.col(i - 1))/(i + 1);

      Lp = N2*(L.col(N) - y*L.col(N1))/(1 - y.pow(2));

      y0 = y;
      y = y0 - L.col(N1)/Lp;
    }

    unique_ptr<stGLQuad<Real>> output = make_unique<stGLQuad<Real>>();

    output->x = (a*(1 - y) + b*(1 + y))/2;
    output->w = ArrayXr<Real>::Constant(N1, b - a) /
        ((1 - y.pow(2))*Lp.pow(2)) * pow(static_cast<Real>(N2) / N1, 2);
    return output;
  }

  // Missing functionality
  template <class Real>
  unique_ptr<stRtfunc<Real>> auxPrepareIntegrals(size_t N_int, sInt type) {
    unique_ptr<stRtfunc<Real>> output = make_unique<stRtfunc<Real>>();
    output->Nb_theta = N_int;

    switch (type) {
      case GAUSS: {
        unique_ptr<stGLQuad<Real>> GL_quad = auxInitLegendreQuad<Real>(N_int);
        output->theta = acos(GL_quad->x);
        output->w_theta = GL_quad->w;
        break;
      }

      case RECTANGLE: {
        output->theta = ArrayXr<Real>::LinSpaced(N_int, 0, mp_pi<Real>());
        Real d_theta = mp_pi<Real>()/(N_int - 1);
        output->w_theta = d_theta*sin(output->theta);
        break;
      }

      default:
        cout << "Integration type not recognized" << endl;
    }

    return output;
  }

  template unique_ptr<stGLQuad<double>> auxInitLegendreQuad(size_t, double, double);
  template unique_ptr<stRtfunc<double>> auxPrepareIntegrals(size_t, sInt);
}
