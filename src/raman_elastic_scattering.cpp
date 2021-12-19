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
  //size_t test2 = sphEstimateNB<double>(8, test, params);
  stPQa<double>* test2 = sphCalculatePQ(10, ArrayXi::LinSpaced(11, 0, 10), test, params, 12);
  cout << test2[1].st4MPoe->M11 << endl << test2[1].st4MPoe->M12 << endl << test2[1].st4MPoe->M21 << endl << test2[1].st4MPoe->M22 << endl;
  cout << test2[1].st4MPoe->m << endl << test2[1].st4MPoe->ind1 << endl << test2[1].st4MPoe->ind2 << endl;
  delete params;
  delete test;
  for (int i = 0; i <= 10; i++) {
    delete test2[i].st4MQeo;
    delete test2[i].st4MQoe;
    delete test2[i].st4MPeo;
    delete test2[i].st4MPoe;
  }
  delete[] test2;
  return 0;
}
