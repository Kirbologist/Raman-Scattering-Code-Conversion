#include "raman_elastic_scattering.h"
#include "sph.h"

using namespace std;
using namespace Eigen;
using namespace Raman;

int main(int argc, char** argv) {
  //ArrayXr<double> x = ArrayXr<double>::LinSpaced(5, 1, 5);
  //stBesselProducts<double>* test = sphGetModifiedBesselProducts<double>(10, 0.5, x, 12);
  //destructStBesselProducts<double>(test);
  stRtfunc<double>* test = sphMakeGeometry<double>(5, 100, 150);
  stParams<double>* params = new stParams<double>();
  params->k1 = {{2*M_PI/355}};
  params->s = {{1.35}};
  size_t test2 = sphEstimateNB<double>(8, test, params);
  cout << test2 << endl;
  return 0;
}
