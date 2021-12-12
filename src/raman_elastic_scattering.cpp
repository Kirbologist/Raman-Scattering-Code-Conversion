#include "raman_elastic_scattering.h"
#include "low_level.h"

using namespace std;
using namespace Eigen;
using namespace Raman;

int main(int argc, char** argv) {
  ArrayXr<double> x = ArrayXr<double>::LinSpaced(5, 1, 5);
  LowLevel<double>::stFprow* test = LowLevel<double>::sphGetFpRow(9, 0.5, x);
  cout << "S:" << test->S << endl << "loss prec S:" << test->loss_prec_S << endl;
  delete test;
  return 0;
}
