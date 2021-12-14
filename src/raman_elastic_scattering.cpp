#include "raman_elastic_scattering.h"
#include "low_level.h"

using namespace std;
using namespace Eigen;
using namespace Raman;

int main(int argc, char** argv) {
  ArrayXr<double> x = ArrayXr<double>::LinSpaced(5, 1, 5);
  LowLevel<double>::stBesselProducts* test = LowLevel<double>::sphGetModifiedBesselProducts(10, 0.5, x, 12);
  LowLevel<double>::destructStBesselProducts(test);
  return 0;
}
