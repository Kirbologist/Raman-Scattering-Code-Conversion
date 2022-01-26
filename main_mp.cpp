#include "src/raman_elastic_scattering.hpp"
#include "src_mp/core_mp.hpp"
#include <omp.h>

using namespace boost::multiprecision;
using namespace Eigen;
using namespace Smarties;

extern template void RamanElasticScattering<float, float>(string, string, int);
extern template void RamanElasticScattering<float, double>(string, string, int);
extern template void RamanElasticScattering<float, long double>(string, string, int);
extern template void RamanElasticScattering<float, RamanFloat>(string, string, int);
extern template void RamanElasticScattering<double, float>(string, string, int);
extern template void RamanElasticScattering<double, double>(string, string, int);
extern template void RamanElasticScattering<double, long double>(string, string, int);
extern template void RamanElasticScattering<double, RamanFloat>(string, string, int);
extern template void RamanElasticScattering<long double, float>(string, string, int);
extern template void RamanElasticScattering<long double, double>(string, string, int);
extern template void RamanElasticScattering<long double, long double>(string, string, int);
extern template void RamanElasticScattering<long double, RamanFloat>(string, string, int);
extern template void RamanElasticScattering<RamanFloat, float>(string, string, int);
extern template void RamanElasticScattering<RamanFloat, double>(string, string, int);
extern template void RamanElasticScattering<RamanFloat, long double>(string, string, int);
extern template void RamanElasticScattering<RamanFloat, RamanFloat>(string, string, int);

string in_file_name = "config.txt";
string out_dir = "output";

int main(int argc, char** argv) {
  vector<string> flags = {"--input=", "--output-dir="};
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
            out_dir = arg.substr(flags[j].size());
            break;
          }
        }
      }
    }
  }

  int num_CPUs = GetNumCPUs(in_file_name);
  if (num_CPUs > 0)
    omp_set_num_threads(num_CPUs);
  bool can_write_output = CanWriteOutput(in_file_name);
  if (!can_write_output)
    out_dir = "";
  int num_particle_CPUs = GetNumParticleCPUs(in_file_name);

  std::array<CalcType, 2> calc_types = GetCalcType(in_file_name);
  switch (calc_types[0]) {
    case CalcType::SINGLE : {
      switch (calc_types[1]) {
        case CalcType::SINGLE : {
          RamanElasticScattering<float, float>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        case CalcType::DOUBLE : {
          RamanElasticScattering<float, double>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        case CalcType::QUAD : {
          RamanElasticScattering<float, long double>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        case CalcType::CUSTOM : {
          RamanElasticScattering<float, RamanFloat>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        default : {
          cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, out_dir, num_particle_CPUs);
        }
      }
      break;
    }
    case CalcType::DOUBLE : {
      switch (calc_types[1]) {
        case CalcType::SINGLE : {
          RamanElasticScattering<double, float>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        case CalcType::DOUBLE : {
          RamanElasticScattering<double, double>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        case CalcType::QUAD : {
          RamanElasticScattering<double, long double>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        case CalcType::CUSTOM : {
          RamanElasticScattering<double, RamanFloat>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        default : {
          cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, out_dir, num_particle_CPUs);
        }
      }
      break;
    }
    case CalcType::QUAD : {
      switch (calc_types[1]) {
        case CalcType::SINGLE : {
          RamanElasticScattering<long double, float>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        case CalcType::DOUBLE : {
          RamanElasticScattering<long double, double>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        case CalcType::QUAD : {
          RamanElasticScattering<long double, long double>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        case CalcType::CUSTOM : {
          RamanElasticScattering<long double, RamanFloat>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        default : {
          cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, out_dir, num_particle_CPUs);
        }
      }
      break;
    }
    case CalcType::CUSTOM : {
      switch (calc_types[1]) {
        case CalcType::SINGLE : {
          RamanElasticScattering<RamanFloat, float>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        case CalcType::DOUBLE : {
          RamanElasticScattering<RamanFloat, double>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        case CalcType::QUAD : {
          RamanElasticScattering<RamanFloat, long double>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        case CalcType::CUSTOM : {
          RamanElasticScattering<RamanFloat, RamanFloat>(in_file_name, out_dir, num_particle_CPUs);
          break;
        }
        default : {
          cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, out_dir, num_particle_CPUs);
        }
      }
      break;
    }
    default : {
      cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
      RamanElasticScattering<double, double>(in_file_name, out_dir, num_particle_CPUs);
    }
  }

  return 0;
}
