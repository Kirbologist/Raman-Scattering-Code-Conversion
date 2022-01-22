#include "src/raman_elastic_scattering.hpp"
#include "src_mp/core_mp.hpp"
#include <omp.h>

using namespace boost::multiprecision;
using namespace Eigen;
using namespace Smarties;

extern template void RamanElasticScattering<float, float>(string, string);
extern template void RamanElasticScattering<float, double>(string, string);
extern template void RamanElasticScattering<float, long double>(string, string);
extern template void RamanElasticScattering<float, raman_float>(string, string);
extern template void RamanElasticScattering<double, float>(string, string);
extern template void RamanElasticScattering<double, double>(string, string);
extern template void RamanElasticScattering<double, long double>(string, string);
extern template void RamanElasticScattering<double, raman_float>(string, string);
extern template void RamanElasticScattering<long double, float>(string, string);
extern template void RamanElasticScattering<long double, double>(string, string);
extern template void RamanElasticScattering<long double, long double>(string, string);
extern template void RamanElasticScattering<long double, raman_float>(string, string);
extern template void RamanElasticScattering<raman_float, float>(string, string);
extern template void RamanElasticScattering<raman_float, double>(string, string);
extern template void RamanElasticScattering<raman_float, long double>(string, string);
extern template void RamanElasticScattering<raman_float, raman_float>(string, string);

const string in_file_name = "config.txt";
const string log_file_name = "output/log.txt";

int main(int argc, char** argv) {
  int num_CPUs = GetNumCPUs("config.txt");
  if (num_CPUs > 0)
    omp_set_num_threads(num_CPUs);
  std::array<CalcType, 2> calc_types = GetCalcType(in_file_name);
  switch (calc_types[0]) {
    case CalcType::SINGLE : {
      switch (calc_types[1]) {
        case CalcType::SINGLE : {
          RamanElasticScattering<float, float>(in_file_name, log_file_name);
          break;
        }
        case CalcType::DOUBLE : {
          RamanElasticScattering<float, double>(in_file_name, log_file_name);
          break;
        }
        case CalcType::QUAD : {
          RamanElasticScattering<float, long double>(in_file_name, log_file_name);
          break;
        }
        case CalcType::CUSTOM : {
          RamanElasticScattering<float, raman_float>(in_file_name, log_file_name);
          break;
        }
        default : {
          cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, log_file_name);
        }
      }
      break;
    }
    case CalcType::DOUBLE : {
      switch (calc_types[1]) {
        case CalcType::SINGLE : {
          RamanElasticScattering<double, float>(in_file_name, log_file_name);
          break;
        }
        case CalcType::DOUBLE : {
          RamanElasticScattering<double, double>(in_file_name, log_file_name);
          break;
        }
        case CalcType::QUAD : {
          RamanElasticScattering<double, long double>(in_file_name, log_file_name);
          break;
        }
        case CalcType::CUSTOM : {
          RamanElasticScattering<double, raman_float>(in_file_name, log_file_name);
          break;
        }
        default : {
          cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, log_file_name);
        }
      }
      break;
    }
    case CalcType::QUAD : {
      switch (calc_types[1]) {
        case CalcType::SINGLE : {
          RamanElasticScattering<long double, float>(in_file_name, log_file_name);
          break;
        }
        case CalcType::DOUBLE : {
          RamanElasticScattering<long double, double>(in_file_name, log_file_name);
          break;
        }
        case CalcType::QUAD : {
          RamanElasticScattering<long double, long double>(in_file_name, log_file_name);
          break;
        }
        case CalcType::CUSTOM : {
          RamanElasticScattering<long double, raman_float>(in_file_name, log_file_name);
          break;
        }
        default : {
          cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, log_file_name);
        }
      }
      break;
    }
    case CalcType::CUSTOM : {
      switch (calc_types[1]) {
        case CalcType::SINGLE : {
          RamanElasticScattering<raman_float, float>(in_file_name, log_file_name);
          break;
        }
        case CalcType::DOUBLE : {
          RamanElasticScattering<raman_float, double>(in_file_name, log_file_name);
          break;
        }
        case CalcType::QUAD : {
          RamanElasticScattering<raman_float, long double>(in_file_name, log_file_name);
          break;
        }
        case CalcType::CUSTOM : {
          RamanElasticScattering<raman_float, raman_float>(in_file_name, log_file_name);
          break;
        }
        default : {
          cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, log_file_name);
        }
      }
      break;
    }
    default : {
      cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
      RamanElasticScattering<double, double>(in_file_name, log_file_name);
    }
  }

  return 0;
}
