/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.
*/

#include "core.hpp"
#include "raman_elastic_scattering.hpp"
#include <fstream>
#include <omp.h>

using namespace Eigen;
using namespace Smarties;

CalcType String2CalcType(string calc_type_string) {
  if (calc_type_string == "single")
    return CalcType::SINGLE;
  else if (calc_type_string == "double")
    return CalcType::DOUBLE;
  else if (calc_type_string == "quad")
    return CalcType::QUAD;
  else if (calc_type_string == "custom")
    return CalcType::CUSTOM;
  return CalcType::NONE;
}

string GetOption(string in_file_name, string option) {
  ifstream in_file(in_file_name, ios::in);
  if (!in_file.good())
    throw runtime_error("Error: cannot open config.txt");

  string line;
  do {
    if (in_file.peek() == EOF)
      throw runtime_error("Error: cannot find option " + option);
    getline(in_file, line);
  } while (line.find(option) == string::npos);
  in_file.close();
  return line.substr(line.find(option) + option.size());
}

bool CanWriteOutput(string in_file_name) {
  string can_write_output = GetOption(in_file_name, "Print output to file:");
  return can_write_output == "yes";
}

std::array<CalcType, 2> GetCalcType(string in_file_name) {
  string calc_types = GetOption(in_file_name, "Calculation type:");
  size_t delim_offset = calc_types.find("-");
  std::array<CalcType, 2> output;
  if (delim_offset == string::npos) {
    output[0] = String2CalcType(calc_types);
    output[1] = String2CalcType(calc_types);
  } else {
    output[0] = String2CalcType(calc_types.substr(0, delim_offset));
    output[1] = String2CalcType(calc_types.substr(delim_offset + 1));
  }
  return output;
}

int GetNumCPUs(string in_file_name) {
  string value = GetOption(in_file_name, "No. of CPUs:");
  int num_CPUs = 1;
  try {
    num_CPUs = stoi(value);
  } catch (invalid_argument&) {
    cerr << "Cannot read value of option \"No. of CPUs:\". Using default value.";
  }
  return num_CPUs;
}

void MultiPrint(string out_string, string out_file_name, bool write_output) {
  #pragma omp critical (multi_print)
  {
    cout << out_string;
    cout.flush();
    if (write_output) {
      ofstream out_file;
      out_file.open(out_file_name, ios::out | ios::app);
      if (!out_file.is_open()) {
        cerr << "Warning: cannot open output file. Some output won't be written." << endl;
        out_file.close();
        write_output = false;
      }
      out_file << out_string;
      out_file.flush();
      out_file.close();
    }
  }
}
