#include "src/raman_elastic_scattering.h"
#include <boost/multiprecision/mpfr.hpp>

int main(int argc, char** argv) {
  string calc_type = (argc > 7) ? argv[7] : "double";
  RamanElasticScattering<double>(argc, argv);
  return 0;
}
