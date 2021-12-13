#include "raman_elastic_scattering.h"
#include "low_level.h"

using namespace std;
using namespace Eigen;
using namespace Raman;

int main(int argc, char** argv) {
  ArrayXr<double> x = ArrayXr<double>::LinSpaced(5, 1, 5);
  LowLevel<double>::stFpovx* test = LowLevel<double>::sphGetFpovx(10, 0.5, x);
  delete[] test->Fpovx;
  delete test->rb_chi;
  delete test->rb_psi;
  delete test;
  return 0;
}
