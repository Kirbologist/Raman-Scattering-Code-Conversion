#include "src/raman_elastic_scattering.hpp"
#include <omp.h>

using namespace std;

extern template void RamanElasticScattering<double>(string, string);
extern template void RamanElasticScattering<long double>(string, string);
extern template void RamanElasticScatteringMpDouble<long double>(string, string);

const string in_file_name = "config.txt";
const string log_file_name = "output/log.txt";

int main(int argc, char** argv) {
  omp_set_num_threads(THREADS);
  string calc_type = GetCalcType(in_file_name);
  if (calc_type == "quad")
    RamanElasticScattering<long double>(in_file_name, log_file_name);
  else if (calc_type == "quad-double")
    RamanElasticScatteringMpDouble<long double>(in_file_name, log_file_name);
  else {
    if (calc_type == "mp" || calc_type == "mp-double")
      cerr << "Warning: attempted to run with multiprecision. Please recompile with command 'make mp'. " <<
          "Running with double precision." << endl;
    else if (calc_type != "double")
      cerr << "Warning: unrecognised calculation type. Running with double precision." << endl;
    RamanElasticScattering<double>(in_file_name, log_file_name);
  }
  return 0;
}
