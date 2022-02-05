/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Based on the SMARTIES MATLAB package by W.R.C. Somerville, B. Auguié, E.C. Le Ru
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.
*/

#include "raman_elastic_scattering.hpp"
#include "core_mp.hpp"
#include <omp.h>

using namespace boost::multiprecision;
using namespace Eigen;
using namespace Smarties;

extern template void RamanElasticScattering<float, float>(string, string);
extern template void RamanElasticScattering<float, double>(string, string);
extern template void RamanElasticScattering<float, long double>(string, string);
extern template void RamanElasticScattering<float, RamanFloat>(string, string);
extern template void RamanElasticScattering<double, float>(string, string);
extern template void RamanElasticScattering<double, double>(string, string);
extern template void RamanElasticScattering<double, long double>(string, string);
extern template void RamanElasticScattering<double, RamanFloat>(string, string);
extern template void RamanElasticScattering<long double, float>(string, string);
extern template void RamanElasticScattering<long double, double>(string, string);
extern template void RamanElasticScattering<long double, long double>(string, string);
extern template void RamanElasticScattering<long double, RamanFloat>(string, string);
extern template void RamanElasticScattering<RamanFloat, float>(string, string);
extern template void RamanElasticScattering<RamanFloat, double>(string, string);
extern template void RamanElasticScattering<RamanFloat, long double>(string, string);
extern template void RamanElasticScattering<RamanFloat, RamanFloat>(string, string);

string in_file_name = "config.txt"; // input text file containing all the parameters
string out_dir = "output"; // directory for program to write outputs to

/*
Runs Raman Elastic Scattering calculations with fully custom multiprecision functionality
(requires multiprecision libraries to be installed)
*/
int main(int argc, char** argv) {
  // Check flags
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

  // Get the rest of the run parameters
  int num_CPUs = GetNumCPUs(in_file_name);
  if (num_CPUs > 0)
    omp_set_num_threads(num_CPUs);
  bool can_write_output = CanWriteOutput(in_file_name);
  if (!can_write_output)
    out_dir = "";
  std::array<CalcType, 2> calc_types = GetCalcType(in_file_name);

  // Run the calculations
  switch (calc_types[0]) {
    case CalcType::SINGLE : {
      switch (calc_types[1]) {
        case CalcType::SINGLE : {
          RamanElasticScattering<float, float>(in_file_name, out_dir);
          break;
        }
        case CalcType::DOUBLE : {
          RamanElasticScattering<float, double>(in_file_name, out_dir);
          break;
        }
        case CalcType::QUAD : {
          RamanElasticScattering<float, long double>(in_file_name, out_dir);
          break;
        }
        case CalcType::CUSTOM : {
          RamanElasticScattering<float, RamanFloat>(in_file_name, out_dir);
          break;
        }
        default : {
          cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, out_dir);
        }
      }
      break;
    }
    case CalcType::DOUBLE : {
      switch (calc_types[1]) {
        case CalcType::SINGLE : {
          RamanElasticScattering<double, float>(in_file_name, out_dir);
          break;
        }
        case CalcType::DOUBLE : {
          RamanElasticScattering<double, double>(in_file_name, out_dir);
          break;
        }
        case CalcType::QUAD : {
          RamanElasticScattering<double, long double>(in_file_name, out_dir);
          break;
        }
        case CalcType::CUSTOM : {
          RamanElasticScattering<double, RamanFloat>(in_file_name, out_dir);
          break;
        }
        default : {
          cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, out_dir);
        }
      }
      break;
    }
    case CalcType::QUAD : {
      switch (calc_types[1]) {
        case CalcType::SINGLE : {
          RamanElasticScattering<long double, float>(in_file_name, out_dir);
          break;
        }
        case CalcType::DOUBLE : {
          RamanElasticScattering<long double, double>(in_file_name, out_dir);
          break;
        }
        case CalcType::QUAD : {
          RamanElasticScattering<long double, long double>(in_file_name, out_dir);
          break;
        }
        case CalcType::CUSTOM : {
          RamanElasticScattering<long double, RamanFloat>(in_file_name, out_dir);
          break;
        }
        default : {
          cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, out_dir);
        }
      }
      break;
    }
    case CalcType::CUSTOM : {
      switch (calc_types[1]) {
        case CalcType::SINGLE : {
          RamanElasticScattering<RamanFloat, float>(in_file_name, out_dir);
          break;
        }
        case CalcType::DOUBLE : {
          RamanElasticScattering<RamanFloat, double>(in_file_name, out_dir);
          break;
        }
        case CalcType::QUAD : {
          RamanElasticScattering<RamanFloat, long double>(in_file_name, out_dir);
          break;
        }
        case CalcType::CUSTOM : {
          RamanElasticScattering<RamanFloat, RamanFloat>(in_file_name, out_dir);
          break;
        }
        default : {
          cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
          RamanElasticScattering<double, double>(in_file_name, out_dir);
        }
      }
      break;
    }
    default : {
      cerr << "Warning: unrecognised calculation type. Running with double-double precision." << endl;
      RamanElasticScattering<double, double>(in_file_name, out_dir);
    }
  }

  return 0;
}
