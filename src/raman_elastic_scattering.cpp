#include "raman_elastic_scattering.h"
#include "low_level.h"

using namespace std;
using namespace Eigen;
using namespace Raman;

int main(int argc, char** argv) {
  ArrayXr<double> x = ArrayXr<double>::LinSpaced(5, 1, 5);
  LowLevel<double>::stBessel* test = LowLevel<double>::sphGetXiPsi(10, 0.5, x, 12);
  delete[] test->xi_psi;
  delete[] test->psi_psi;
  delete test->chi_n;
  delete test->psi_n;
  delete test->psi_k;
  delete test;
  return 0;
}
