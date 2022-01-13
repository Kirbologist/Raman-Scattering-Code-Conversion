#include "src/raman_elastic_scattering.hpp"
#include "src/smarties.hpp"
#include "src_mp/core_mp.hpp"

using namespace boost::multiprecision;
using namespace Eigen;
using namespace Smarties;

extern template void RamanElasticScattering<double>(int, char**);
extern template void RamanElasticScattering<long double>(int, char**);
extern template void RamanElasticScattering<raman_float>(int, char**);
extern template void RamanElasticScatteringMpDouble<long double>(int, char**);
extern template void RamanElasticScatteringMpDouble<raman_float>(int, char**);

int main(int argc, char** argv) {
  string calc_type = (argc > 7) ? argv[7] : "double";
  if (calc_type == "mp")
    RamanElasticScattering<raman_float>(argc, argv);
  else if (calc_type == "mp-double")
    RamanElasticScatteringMpDouble<raman_float>(argc, argv);
  else if (calc_type == "quad")
    RamanElasticScattering<long double>(argc, argv);
  else if (calc_type == "quad-double")
    RamanElasticScatteringMpDouble<long double>(argc, argv);
  else
    RamanElasticScattering<double>(argc, argv);
  return 0;
}
