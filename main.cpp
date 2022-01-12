#include "src/raman_elastic_scattering.hpp"

extern template RamanElasticScattering<double>(int, char**);

int main(int argc, char** argv) {
  RamanElasticScattering<double>(argc, argv);
  return 0;
}
