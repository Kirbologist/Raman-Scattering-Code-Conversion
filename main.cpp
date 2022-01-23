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

string in_file_name = "config.txt";
string out_file_name = "output/results.txt";
string log_file_name = "output/log.txt";

int main(int argc, char** argv) {
  vector<string> flags = {"--input=", "--log="};
  for (int i = 0; i < argc; i++) {
    string arg = argv[i];
    for (size_t j = 0; j < flags.size(); j++) {
      if (arg.find(flags[j]) == 0) {
        switch (j) {
          case 0 : {
            in_file_name = arg.substr(flags[j].size());
            break;
          }
          case 1 : {
            log_file_name = arg.substr(flags[j].size());
            break;
          }
        }
      }
    }
  }

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
