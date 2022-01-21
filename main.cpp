#include "src/raman_elastic_scattering.hpp"
#include <omp.h>

using namespace std;

extern template void RamanElasticScattering<float, float>(string, string);
extern template void RamanElasticScattering<float, double>(string, string);
extern template void RamanElasticScattering<float, long double>(string, string);
extern template void RamanElasticScattering<double, float>(string, string);
extern template void RamanElasticScattering<double, double>(string, string);
extern template void RamanElasticScattering<double, long double>(string, string);
extern template void RamanElasticScattering<long double, float>(string, string);
extern template void RamanElasticScattering<long double, double>(string, string);
extern template void RamanElasticScattering<long double, long double>(string, string);

const string in_file_name = "config.txt";
const string log_file_name = "output/log.txt";

int main(int argc, char** argv) {
  omp_set_num_threads(THREADS);
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
          cerr << "Warning: attempted to run with multiprecision. Please recompile with command 'make mp'. " <<
              "Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, log_file_name);
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
          cerr << "Warning: attempted to run with multiprecision. Please recompile with command 'make mp'. " <<
              "Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, log_file_name);
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
          cerr << "Warning: attempted to run with multiprecision. Please recompile with command 'make mp'. " <<
              "Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, log_file_name);
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
      cerr << "Warning: attempted to run with multiprecision. Please recompile with command 'make mp'. " <<
          "Running with double-double precision." << endl;
      RamanElasticScattering<double, double>(in_file_name, log_file_name);
      break;
    }
    default : {
      cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
      RamanElasticScattering<double, double>(in_file_name, log_file_name);
    }
  }

  return 0;
}
