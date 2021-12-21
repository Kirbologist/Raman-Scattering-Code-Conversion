#include "raman_elastic_scattering.h"
#include "sph.h"

using namespace std;
using namespace Eigen;
using namespace Raman;

int main(int argc, char** argv) {
  //ArrayXr<double> x = ArrayXr<double>::LinSpaced(5, 1, 5);
  //stBesselProducts<double>* test = sphGetModifiedBesselProducts<double>(10, 0.5, x, 12);
  //destructStBesselProducts<double>(test);
  unique_ptr<stRtfunc<double>> test = sphMakeGeometry<double>(5, 100, 150);
  unique_ptr<stParams<double>> params = make_unique<stParams<double>>();
  params->k1 = {{2*M_PI/355}};
  params->s = {{1.35}};
  //size_t test2 = sphEstimateNB<double>(8, test, params);
  unique_ptr<vector<stPQ<double>>> test2 = sphCalculatePQ(10, ArrayXi::LinSpaced(11, 0, 10), test, params, 12);
  cout << test2->size() << endl;
  cout << "Now trying to access test2->at(0):" << endl;
  cout << test2->at(0).st_4M_P_oe->m << endl;
  cout << "Now trying to access test2->at(10):" << endl;
  cout << test2->at(10).st_4M_P_oe->m << endl;
  test2->at(10).st_4M_P_oe->m = 15;
  cout << "Now trying to access test2->at(0):" << endl;
  cout << test2->at(0).st_4M_P_oe->m << endl;
  cout << "Now trying to access test2->at(10):" << endl;
  cout << test2->at(10).st_4M_P_oe->m << endl;
  return 0;
}
