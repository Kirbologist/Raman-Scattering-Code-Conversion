#include "raman_aux.h"
#include "math.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  stGLQuad<Real>* auxInitLegendreQuad(size_t N1, Real a, Real b) {
    int N = N1 - 1, N2 = N1 + 1;
    ArrayXr<Real> xu, y, y0(N1), Lp;
    ArrayXXr<Real> L(N1, N2);

    xu = ArrayXr<Real>::LinSpaced(N1, -1, 1);
    y = cos((2*ArrayXr<Real>::LinSpaced(N1, 0, N) + 1)*PI/(2*N + 2)) +
        (0.27/N1)*sin(PI*xu*N/N2);
    y0.fill(2);
    L.col(0).setOnes();

    int n_iter = 0;
    while ((n_iter < 15) && EPS <
        abs(acos(y) - acos(y0.template cast<complex<Real>>())).maxCoeff()) {
      n_iter++;
      L.col(1) = y;

      for (size_t i = 1; i < N1; i++)
        L.col(i + 1) = ((2*i + 1)*y * L.col(i) - i*L.col(i - 1))/(i + 1);

      Lp = N2*(L.col(N) - y*L.col(N1))/(1 - y.pow(2));

      y0 = y;
      y = y0 - L.col(N1)/Lp;
    }

    stGLQuad<Real>* output = new stGLQuad<Real>();

    output->x = (a*(1 - y) + b*(1 + y))/2;
    output->w = ArrayXr<Real>::Constant(N1, b - a) /
        ((1 - y.pow(2))*Lp.pow(2)) * pow(static_cast<Real>(N2) / N1, 2);
    return output;
  }

  // Missing functionality
  template <class Real>
  stRtfunc<Real>* auxPrepareIntegrals(size_t nNint, sInt type) {
    stRtfunc<Real>* output = new stRtfunc<Real>();
    output->n_Nb_theta = nNint;

    switch (type) {
      case GAUSS: {
        stGLQuad<Real>* weights = auxInitLegendreQuad<Real>(nNint);
        output->theta = acos(weights->x);
        output->w_theta = weights->w;
        delete weights;
        break;
      }

      case RECTANGLE: {
        output->theta = ArrayXr<Real>::LinSpaced(nNint, 0, PI);
        Real dTheta = PI/(nNint - 1);
        output->w_theta = dTheta*sin(output->theta);
        break;
      }

      default:
        cout << "Integration type not recognized" << endl;
    }

    return output;
  }

  template stGLQuad<double>* auxInitLegendreQuad(size_t, double, double);
  template stRtfunc<double>* auxPrepareIntegrals(size_t, sInt);
}
