#include "src/raman_elastic_scattering.hpp"

extern template void RamanElasticScattering<double>(string, string);
extern template void RamanElasticScattering<long double>(string, string);
extern template void RamanElasticScatteringMpDouble<long double>(string, string);

const string in_file_name = "config.txt";
const string out_file_name = "output/results.txt";

int main(int argc, char** argv) {
  string calc_type = GetCalcType(in_file_name);
  if (calc_type == "quad")
    RamanElasticScattering<long double>(in_file_name, out_file_name);
  else if (calc_type == "quad-double")
    RamanElasticScatteringMpDouble<long double>(in_file_name, out_file_name);
  else
    RamanElasticScattering<double>(in_file_name, out_file_name);
  return 0;
}
