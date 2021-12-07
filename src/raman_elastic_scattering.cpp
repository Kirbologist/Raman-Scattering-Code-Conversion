#include "raman_elastic_scattering.h"
#include "low_level.h"

using namespace std;
using namespace Eigen;
using namespace Raman;

int main(int argc, char** argv) {
  Array<double, 5, 1> x;
  x = {1, 2, 3, 4, 5};
  Array<double, 1, 4> n;
  n = {6, 7, 8, 9};
  //Array<complex<double>, Dynamic, Dynamic> test = LowLevel<double>::vshRBchi(n, x);
  ArrayXXc<double>* test = LowLevel<double>::vshRBpsi(n, x);
  cout << *test << endl;
  return 0;
}
