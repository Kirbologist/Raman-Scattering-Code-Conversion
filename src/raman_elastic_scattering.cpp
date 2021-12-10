#include "raman_elastic_scattering.h"
#include "low_level.h"

using namespace std;
using namespace Eigen;
using namespace Raman;

int main(int argc, char** argv) {
  HighLevel<double>::stIncPar* test = HighLevel<double>::vshMakeIncidentParams(sIncType::KxEz, 3);
  LowLevel<double>::stIncEabnm* test2 =  LowLevel<double>::vshGetIncidentCoeffs(3, test);
  cout << "a_nm" << test2->a_nm << endl << "b_nm" << test2->b_nm << endl;
  delete test;
  delete test2;
  return 0;
}
