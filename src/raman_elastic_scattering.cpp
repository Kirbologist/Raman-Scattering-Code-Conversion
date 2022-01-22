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

std::array<CalcType, 2> GetCalcType(string in_file_name) {
  ifstream in_file;
  in_file.open(in_file_name, ios::in);
  if (!in_file.is_open())
    throw runtime_error("Error: cannot open config.txt");

  string option = "Calculation type:";
  string line;
  do {
    if (in_file.peek() == EOF)
      throw runtime_error("Error: cannot find option " + option);
    getline(in_file, line);
  } while (line.find(option) == string::npos);
  string calc_types = line.substr(line.find(option) + option.size());
  in_file.close();
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
  ifstream in_file;
  in_file.open(in_file_name, ios::in);
  if (!in_file.is_open())
    throw runtime_error("Error: cannot open config.txt");

  string option = "No. of CPUs:";
  string line;
  do {
    if (in_file.peek() == EOF)
      throw runtime_error("Error: cannot find option " + option);
    getline(in_file, line);
  } while (line.find(option) == string::npos);
  string num_CPUs = line.substr(line.find(option) + option.size());
  return (num_CPUs.size() > 0 && isdigit(num_CPUs[0])) ? stoi(num_CPUs, nullptr) : 1;
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
