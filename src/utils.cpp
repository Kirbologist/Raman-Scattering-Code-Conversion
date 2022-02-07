/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Based on the SMARTIES MATLAB package by W.R.C. Somerville, B. Auguié, E.C. Le Ru
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.


This files runs the storeGLquadrature script, which was located in the utils folder in the original MATLAB
code. Note, this needs to be linked with the multiprecision libraries gmp and mpfr.
storeGLquadrature is ran using whichever one or two calculation types are written in the input file.
*/

#include "smarties_aux.hpp"
#include "raman_elastic_scattering.hpp"
#include "core_mp.hpp"

using namespace Eigen;
using namespace std;

string in_file_name = "config.txt";

int main(int argc, char** argv) {
  // Check flags
  vector<string> flags = {"--input="};
  for (int i = 0; i < argc; i++) {
    string arg = argv[i];
    for (size_t j = 0; j < flags.size(); j++) {
      if (arg.find(flags[j]) == 0) {
        switch (j) {
          case 0 : {
            in_file_name = arg.substr(flags[j].size());
            break;
          }
        }
      }
    }
  }
  
  std::array<CalcType, 2> calc_types = GetCalcType(in_file_name);
  for (CalcType calc_type : calc_types) {
    switch (calc_type) {
      case SINGLE : {
        StoreGLquadrature<float>();
        break;
      }
      case DOUBLE : {
        StoreGLquadrature<double>();
        break;
      }
      case QUAD : {
        StoreGLquadrature<long double>();
        break;
      }
      case CUSTOM : {
        StoreGLquadrature<RamanFloat>();
        break;
      }
      default :
        continue;
    }
  }
  return 0;
}
