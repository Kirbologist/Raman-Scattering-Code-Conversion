/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.


This code is used to calculate Raman scattering off of spheroids.
None of the code was included in the original SMARTIES package.
Many functions are used to provide I/O support to the main function.
*/

#ifndef RAMAN_ELASTIC_SCATTERING_HPP
#define RAMAN_ELASTIC_SCATTERING_HPP

#include "core.hpp"
#include "smarties.hpp"
#include <boost/lexical_cast.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <exception>

using namespace std;
using namespace Eigen;
using namespace Smarties;
using boost::lexical_cast;
using boost::bad_lexical_cast;

/* Possible floating-point types used in calculations */
enum CalcType {SINGLE, DOUBLE, QUAD, CUSTOM, NONE};

/*
Controllable parameters that are used for calculations.
_rm suffix denotes parameters used for excitation T-matrix,
absense of _rm suffix denotes parameters used for Raman T-matrix.
*/
template <class Real>
struct RamanParams {
  // Describes the minimum and maximum of the particle diameter
  Real dia_min = 1000;
  Real dia_max = 2000;

  int N_rad = 100; // Number of particle sizes to calculate for
  int N_theta_p = 19; // Number of particle orientations to calculate for
  Real phi_p = 0; // Angle of rotation of particle about the z semiaxis

  // Number of spherical coordinates used throughout the volume of the particle
  int N_r = 100;
  int N_theta = 320;
  int N_phi = 320;

  // Values of h used, where h = a/c. N_h is the number of h values
  int N_h = 1;
  Real h_min = static_cast<Real>(1.0)/3;
  Real h_max = static_cast<Real>(1.0)/3;

  // Incident light parameters. Used for stParams struct
  Real epsilon1 = 1;
  Real epsilon2 = pow(static_cast<Real>(1.35), 2);
  Real epsilon2_rm = pow(static_cast<Real>(1.344), 2);
  Real lambda = 355;
  Real lambda_rm = 403.7;

  // Used for stParams struct
  int Nb_theta = 1000;
  int Nb_theta_pst = 1;

  // Vectors of all possible radii/theta/h values
  ArrayXr<Real> rad_var;
  ArrayXr<Real> theta_p_var;
  ArrayXr<Real> h_var;
};

// All non-template functions defined in raman_elastic_scattering.cpp

/*
Get the 2 types to be used for floating-point calculations from a text file.
Inputs:
  in_file_name: path of the text file which contains the calculation type
Output:
  String array of size 2 that specifies 2 calculation types.
*/
std::array<CalcType, 2> GetCalcType(string in_file_name);

/*
Get the value given for some parameter from a text file.
Inputs:
  in_file_name: path of the text file which contains the option
  option: The exact string of the specific option/parameter to look for
Output:
  The user-defined value given to the option, as a string
*/
string GetOption(string in_file_name, string option);

/*
Determine if it's okay to write calculation info to a file from an option in a text file.
I.e. it checks the "Print output to file:" parameter.
Inputs:
  in_file_name: path of the text file which specifies if it's okay to write to output or not
Output:
  True if the file says "yes", i.e. it is okay, false otherwise
*/
bool CanWriteOutput(string in_file_name);

/*
Determine number of CPUs to use from a text file. I.e. it checks the "No. of CPUs" parameter.
Inputs:
  in_file_name: path of the text file which specifies the number of CPUs.
Output:
  The number of CPUs that are allowed
*/
int GetNumCPUs(string in_file_name);

/*
Prints a string to both standard output (the terminal) and appends it to a text ilfe (if allowed).
Inputs:
  out_string: the string to print/write
  out_file_name: path of the text file to write to
  write_output: true enables writing to output, false disables it.
*/
void MultiPrint(string out_string, string out_file_name, bool write_output = false);

/*
Converts a string that's a fraction of two floats into a single float value
Throws a bad_lexical_cast error if string can't be read or converted.
Input:
  `frac` - must either be a floating-point literal, or two floating-point literals deliminated by a '/' character.
Output:
  the calculated value as a floating-point variable.
*/
template <class Real>
Real Frac2Float(string frac) {
  Real output;
  size_t offset = frac.find("/");
  if (offset == string::npos) {
    output = static_cast<Real>(lexical_cast<Real>(frac));
    return output;
  }
  output = static_cast<Real>(lexical_cast<Real>(frac.substr(0, offset)));
  output /= static_cast<Real>(lexical_cast<Real>(frac.substr(offset + 1)));
  return output;
}

/*
Creates RamanParams struct by reading values from text file.
The file can be opened or unopened, but must not be modified during this function call.
Inputs:
  `in_file_name` - path to the file
Output:
  A unique pointer to a new RamanParams struct with member values based on the parameters given in `in_file_name`.
Dependencies:
  GetOption, Frac2Float
*/
template <class Real>
unique_ptr<RamanParams<Real>> LoadParams(string in_file_name) {
  ifstream in_file;
  in_file.open(in_file_name, ios::in);
  if (!in_file.is_open())
    throw runtime_error("Error: cannot open " + in_file_name);

  auto output = make_unique<RamanParams<Real>>();

  string line;
  // Lists all lines/options to check for. Options with floating-point values are listed first,
  // then options with integral values are checked.
  vector<string> options {
    "Minimum diameter:", "Maximum diameter:", "Particle phi:", "Minimum h:", "Maximum h:",
    "epsilon1:", "epsilon2:", "Raman epsilon2:", "lambda:", "Raman lambda:",
    "No. of particle radii:", "No. of particle thetas:", "No. of h ratios:",
    "No. of r-coordinates:", "No. of theta-coordinates:", "No. of phi-coordinates:",
    "Nb_theta:", "Nb_theta_pst:"
  };
  std::array<Real, 10> float_params = {
    output->dia_min, output->dia_max, output->phi_p, output->h_min, output->h_max,
    output->epsilon1, output->epsilon2, output->epsilon2_rm, output->lambda, output->lambda_rm
  };
  std::array<int, 8> int_params = {
    output->N_rad, output->N_theta_p, output->N_h, output->N_r, output->N_theta, output->N_phi,
    output->Nb_theta, output->Nb_theta_pst
  };
  for (size_t i = 0; i < options.size(); i++) {
    string option = options[i];
    string value = GetOption(in_file_name, option);

    try {
      if (i < float_params.size())
        float_params[i] = Frac2Float<Real>(value);
      else if (i < float_params.size() + int_params.size())
        int_params[i - float_params.size()] = lexical_cast<int>(value);
    } catch(bad_lexical_cast&) {
      cerr << "Cannot read value of option \"" << option << "\". Using default value.";
    }

    output->dia_min = float_params[0];
    output->dia_max = float_params[1];
    output->phi_p = float_params[2];
    output->h_min = float_params[3];
    output->h_max = float_params[4];
    output->epsilon1 = float_params[5];
    output->epsilon2 = float_params[6];
    output->epsilon2_rm = float_params[7];
    output->lambda = float_params[8];
    output->lambda_rm = float_params[9];
    output->N_rad = int_params[0];
    output->N_theta_p = int_params[1];
    output->N_h = int_params[2];
    output->N_r = int_params[3];
    output->N_theta = int_params[4];
    output->N_phi = int_params[5];
    output->Nb_theta = int_params[6];
    output->Nb_theta_pst = int_params[7];
  }
  ArrayXr<Real> dia_var = ArrayXr<Real>::LinSpaced(output->N_rad, output->dia_min, output->dia_max);
  output->rad_var = dia_var/2;
  output->theta_p_var = ArrayXr<Real>::LinSpaced(output->N_theta_p, 0, mp_pi<Real>()/2);
  output->h_var = ArrayXr<Real>::LinSpaced(output->N_h, output->h_min, output->h_max);
  return output;
}

/*
Takes a RamanParams struct and uses it to generate an stParams struct for use with SMARTIES functions.
stParams member values are based on RamanParams member values, and on the radius and h value
at rad_ind and h_ind of rad_var and h_var respectively.
Inputs:
  `raman_params` - unique pointer to a RamanParams struct
  `rad_ind` - the index of the entry to use in `raman_params`->rad_var
  `h_ind` - the index of the entry to use in `raman_params`->h_var
  `type` - the type of values the output struct should have. If it's "rm",
           then the values denote the parameters of the Raman T-matrix. Otherwise they denote the parameters
           of the excitation T-matrix.
Output:
  returns a unique pointer to a new stParams struct
Dependencies:
  mp_pi
*/
template <class Real>
unique_ptr<stParams<Real>> Raman2SmartiesParams(
    const unique_ptr<RamanParams<Real>>& raman_params, int rad_ind, int h_ind, string type = string()) {
  auto params = make_unique<stParams<Real>>();
  params->epsilon1 = raman_params->epsilon1;
  params->lambda = ArrayXr<Real>(1);
  params->epsilon2 = ArrayXr<Real>(1);
  params->k1 = ArrayXr<Real>(1);
  params->s = ArrayXr<Real>(1);
  if (type == "rm") {
    params->lambda(0) = raman_params->lambda_rm;
    params->epsilon2(0) = raman_params->epsilon2_rm;
    params->k1(0) = 2*mp_pi<Real>() / raman_params->lambda_rm * sqrt(raman_params->epsilon1);
    params->s(0) = sqrt(raman_params->epsilon2_rm) / sqrt(raman_params->epsilon1);
  } else {
    params->lambda(0) = raman_params->lambda;
    params->epsilon2(0) = raman_params->epsilon2;
    params->k1(0) = 2*mp_pi<Real>() / raman_params->lambda * sqrt(raman_params->epsilon1);
    params->s(0) = sqrt(raman_params->epsilon2) / sqrt(raman_params->epsilon1);
  }
  params->Nb_theta = raman_params->Nb_theta;
  params->Nb_theta_pst = raman_params->Nb_theta_pst;

  Real h = raman_params->h_var(h_ind);
  Real radius = raman_params->rad_var(rad_ind);
  if (h > 1) {
    params->a = radius;
    params->c = radius / h;
  } else {
    params->a = radius * h;
    params->c = radius;
  }
  // The greater N is, the better the T-matrix converges.
  params->N = 6 + 2*static_cast<int>(ceil(max(params->a, params->c)/40));
  return params;
}

/*
Appends a stamp with time and parameter details of Raman scattering calculations at the end of a text file.
Function opens the file itself. It does not check for any errors.
`Real1` and `Real2` are the two selectable calculation types used for Raman scattering calculations.
DEV NOTE: Can make this into a non-template function using c++20 'concepts' keyword
Inputs:
  - out_file_name: text file to write stamp to
  - raman_params: unique pointer to a RamanParams struct containing the parameters to be written
Dependencies:
  - GetTypeName
*/
template <class Real1, class Real2>
void CreateTimeStamp(string out_file_name, const unique_ptr<RamanParams<Real1>>& raman_params) {
  string main_calc_type = GetTypeName<Real1>();
  string second_calc_type = GetTypeName<Real2>();

  ofstream out_file(out_file_name, ios::out | ios::app);
  auto sys_time = chrono::system_clock::now();
  time_t sys_time_date = chrono::system_clock::to_time_t(sys_time);
  out_file << "Session began at system time: " << ctime(&sys_time_date);
  out_file << "Running RamanElasticScattering with type " << main_calc_type;
  out_file << " and " << second_calc_type << endl;
  out_file.flush();
  out_file.close();
}

/*
Converts an existing stTR struct containing entries of one type into a new stTR struct
with the same entries converted into another type. This could theoretically be
implemented as a template cast operator or template copy constructor.
`From` is the entry type of the source struct and `To` is the entry type of the new struct.
Inputs:
  - st_TR_list: the source stTR struct, with entries of type `From`
Outputs:
  - Returns a new stTR struct, with entries of type `To`.
Dependencies:
  mp_im_unit
*/
template <class From, class To>
vector<unique_ptr<stTR<To>>> ConvertStTRList(const vector<unique_ptr<stTR<From>>>& st_TR_list) {
  vector<unique_ptr<stTR<To>>> output(st_TR_list.size());
  complex<To> I = mp_im_unit<To>();
  for (size_t i = 0; i < st_TR_list.size(); i++) {
    assert(st_TR_list[i]);
    assert(st_TR_list[i]->mat_list.size() == 2);
    assert(st_TR_list[i]->mat_list[0] == "st_4M_T");
    assert(st_TR_list[i]->mat_list[1] == "st_4M_R");
    output[i] = make_unique<stTR<To>>();
    output[i]->mat_list = st_TR_list[i]->mat_list;

    // Direct casting from a complex type to a complex type containing a Boost MPFR type isn't supported.
    output[i]->st_4M_T_eo().M11 = st_TR_list[i]->st_4M_T_eo().M11.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_T_eo().M11.imag().template cast<To>();
    output[i]->st_4M_T_eo().M12 = st_TR_list[i]->st_4M_T_eo().M12.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_T_eo().M12.imag().template cast<To>();
    output[i]->st_4M_T_eo().M21 = st_TR_list[i]->st_4M_T_eo().M21.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_T_eo().M21.imag().template cast<To>();
    output[i]->st_4M_T_eo().M22 = st_TR_list[i]->st_4M_T_eo().M22.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_T_eo().M22.imag().template cast<To>();
    output[i]->st_4M_T_eo().m = st_TR_list[i]->st_4M_T_eo().m;
    output[i]->st_4M_T_eo().ind1 = st_TR_list[i]->st_4M_T_eo().ind1;
    output[i]->st_4M_T_eo().ind2 = st_TR_list[i]->st_4M_T_eo().ind2;

    output[i]->st_4M_T_oe().M11 = st_TR_list[i]->st_4M_T_oe().M11.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_T_oe().M11.imag().template cast<To>();
    output[i]->st_4M_T_oe().M12 = st_TR_list[i]->st_4M_T_oe().M12.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_T_oe().M12.imag().template cast<To>();
    output[i]->st_4M_T_oe().M21 = st_TR_list[i]->st_4M_T_oe().M21.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_T_oe().M21.imag().template cast<To>();
    output[i]->st_4M_T_oe().M22 = st_TR_list[i]->st_4M_T_oe().M22.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_T_oe().M22.imag().template cast<To>();
    output[i]->st_4M_T_oe().m = st_TR_list[i]->st_4M_T_oe().m;
    output[i]->st_4M_T_oe().ind1 = st_TR_list[i]->st_4M_T_oe().ind1;
    output[i]->st_4M_T_oe().ind2 = st_TR_list[i]->st_4M_T_oe().ind2;

    output[i]->st_4M_R_eo().M11 = st_TR_list[i]->st_4M_R_eo().M11.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_R_eo().M11.imag().template cast<To>();
    output[i]->st_4M_R_eo().M12 = st_TR_list[i]->st_4M_R_eo().M12.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_R_eo().M12.imag().template cast<To>();
    output[i]->st_4M_R_eo().M21 = st_TR_list[i]->st_4M_R_eo().M21.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_R_eo().M21.imag().template cast<To>();
    output[i]->st_4M_R_eo().M22 = st_TR_list[i]->st_4M_R_eo().M22.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_R_eo().M22.imag().template cast<To>();
    output[i]->st_4M_R_eo().m = st_TR_list[i]->st_4M_R_eo().m;
    output[i]->st_4M_R_eo().ind1 = st_TR_list[i]->st_4M_R_eo().ind1;
    output[i]->st_4M_R_eo().ind2 = st_TR_list[i]->st_4M_R_eo().ind2;

    output[i]->st_4M_R_oe().M11 = st_TR_list[i]->st_4M_R_oe().M11.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_R_oe().M11.imag().template cast<To>();
    output[i]->st_4M_R_oe().M12 = st_TR_list[i]->st_4M_R_oe().M12.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_R_oe().M12.imag().template cast<To>();
    output[i]->st_4M_R_oe().M21 = st_TR_list[i]->st_4M_R_oe().M21.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_R_oe().M21.imag().template cast<To>();
    output[i]->st_4M_R_oe().M22 = st_TR_list[i]->st_4M_R_oe().M22.real().template cast<To>() +
        I * st_TR_list[i]->st_4M_R_oe().M22.imag().template cast<To>();
    output[i]->st_4M_R_oe().m = st_TR_list[i]->st_4M_R_oe().m;
    output[i]->st_4M_R_oe().ind1 = st_TR_list[i]->st_4M_R_oe().ind1;
    output[i]->st_4M_R_oe().ind2 = st_TR_list[i]->st_4M_R_oe().ind2;
  }
  return output;
}

/*
Performs calculations of the Raman scattering of light off of spheroids with multiprocessing.
If calculations are done for multiple spheroid radii, the calculations for each radii may be
allocated to the threads using dynamic scheduling.
The parameters are based on some input file, the T-matrices are calculated using floating-point type `Real1`,
and the field calculations are calculated using floating-point type `Real2`.
It then prints a summary of results to the standard output (terminal) and in a set of files in a given directory.
Each thread prints its own calculations to its own output file, so each output file is named
based on the parameters on the number of the thread.
Inputs:
  - in_file_name: path to the file containing parameters for use in the Raman scattering calculations.
                  Parameter details are given in the README file.
  - out_dir: directory to write output files to. If empty, output files are written to.
Dependencies:
  mp_pi, ArrayMap, TensorCast, TensorConj, LoadParams, CreateTimeStamp, MultiPrint,
  Raman2SmartiesParams, ConvertStTRList, slvForT, pstScatteringMatrixOA, vshMakeIncidentParams,
  rvhGetFieldCoefficients, pstMakeStructForField, vshEgenThetaAllPhi, vshEthetaForPhi
*/
template <class Real1, class Real2>
void RamanElasticScattering(string in_file_name, string out_dir = "") {
  // Initialise all parameters and variables
  Real2 PI = mp_pi<Real2>();

  unique_ptr<RamanParams<Real1>> raman_params1 = LoadParams<Real1>(in_file_name);
  unique_ptr<RamanParams<Real2>> raman_params2 = LoadParams<Real2>(in_file_name);
  bool write_output = (out_dir != "");

  Tensor3r<Real2> sigma_yz(raman_params2->N_rad, raman_params2->N_h, raman_params2->N_theta_p);
  Tensor3r<Real2> sigma_zy(raman_params2->N_rad, raman_params2->N_h, raman_params2->N_theta_p);
  Tensor3r<Real2> sigma_zz(raman_params2->N_rad, raman_params2->N_h, raman_params2->N_theta_p);
  Tensor3r<Real2> sigma_yy(raman_params2->N_rad, raman_params2->N_h, raman_params2->N_theta_p);
  ArrayXr<Real1> C_sca = ArrayXr<Real1>::Zero(raman_params1->N_rad);
  ArrayXr<Real1> C_ext = ArrayXr<Real1>::Zero(raman_params1->N_rad);
  ArrayXr<Real1> C_abs = ArrayXr<Real1>::Zero(raman_params1->N_rad);
  vector<unique_ptr<stSM<Real2>>> stSM_list(raman_params2->N_rad);

  auto options = make_unique<stOptions>();
  options->get_R = true; // Needed for near fields and will be overridden in any case
  options->delta = 0; // Use delta=-1 to estimate it automatically.
  options->NB = 0; // NB will be estimated automatically
  options->get_symmetric_T = false;
  Real2 phi_p = raman_params2->phi_p;

  int N_r = raman_params2->N_r;
  int N_theta = raman_params2->N_theta;
  int N_phi = raman_params2->N_phi;

  // Initialises a set regularly spaced theta coordinates that are used for internal field calculations
  RowArrayXr<Real2> theta_surf = RowArrayXr<Real2>::LinSpaced(N_theta, 0, PI);
  // Initialises a non-uniform set of N_r radii inside a unit ball
  // such that spherical shells between adjacent radii are equal volume
  ArrayXr<Real2> r_surf_u = ArrayXr<Real2>::LinSpaced(
      N_r, 1/static_cast<Real2>(N_r), 1).pow(static_cast<Real2>(1.0)/3);
  ArrayXXr<Real2> r_mat_u = r_surf_u.replicate(1, N_theta);
  ArrayXXr<Real2> theta_mat = theta_surf.replicate(N_r, 1);

  // Initialise output files
  if (write_output)
    std::filesystem::create_directory(out_dir);

  for (int h_ind = 0; h_ind < raman_params1->N_h; h_ind++) {
    // Split the iterations of the next for-loop among threads for parallel processing.
    #pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < raman_params1->N_rad; k++) {
      string out_file_name = out_dir + "/a2c_" + to_string(raman_params1->h_var(0)) +
          "_to_" + to_string(raman_params1->h_var(last)) + "_dia_" + to_string(raman_params1->dia_min) +
          "_to_" + to_string(raman_params1->dia_max) + "_range_" + to_string(omp_get_thread_num()) + ".txt";
      if (write_output) {
        ofstream out_file(out_file_name, ios::out | ios::app);
        out_file << endl;
        if (!out_file.good()) {
          cerr << "Warning: cannot create or open output file. Some output won't be written." << endl;
          out_file.close();
          write_output = false;
        }
      }
      if (write_output)
        CreateTimeStamp<Real1, Real2>(out_file_name, raman_params1);
      stringstream out_stream;

      out_stream.precision(4);
      out_stream << "aspect ratio " << raman_params1->h_var(h_ind) << endl;
      MultiPrint(out_stream.str(), out_file_name, write_output);
      out_stream.str(string());

      unique_ptr<stParams<Real1>> params1 = Raman2SmartiesParams(raman_params1, k, h_ind);
      unique_ptr<stParams<Real1>> params_rm1 = Raman2SmartiesParams(raman_params1, k, h_ind, "rm");
      unique_ptr<stParams<Real2>> params2 = Raman2SmartiesParams(raman_params2, k, h_ind);
      unique_ptr<stParams<Real2>> params_rm2 = Raman2SmartiesParams(raman_params2, k, h_ind, "rm");
      Real2 a = params2->a;
      Real2 c = params2->c;

      int N = params2->N;
      // Initialises array of radii coordinates that are used for internal field calculations
      // by conforming r_mat_u to the spheroid shape. The ratios of radii along the same angle theta is kept invariant.
      ArrayXXr<Real2> r_mat = a*c*r_mat_u/sqrt(pow(c, 2)*sin(theta_mat).pow(2) + pow(a, 2)*cos(theta_mat).pow(2));
      ArrayXXr<Real2> d_theta = theta_mat(all, seq(1, last)) - theta_mat(all, seq(0, last - 1));
      ArrayXXr<Real2> dr = r_mat(seq(1, last), all) - r_mat(seq(0, last - 1), all);
      ArrayXXr<Real2> dt_dr = d_theta(seq(1, last), all) * dr(all, seq(1, last)); // Jacobian of the integral
      RowArrayXr<Real2> r_row = r_mat.reshaped().transpose();
      RowArrayXr<Real2> theta_row = theta_mat.reshaped().transpose();

      chrono::steady_clock::time_point begin, end;
      chrono::duration<double> elapsed_seconds;

      // Calculating T-matrix for excitation
      begin = chrono::steady_clock::now();
      unique_ptr<stTmatrix<Real1>> T_mat = slvForT(params1, options);
      end = chrono::steady_clock::now();
      elapsed_seconds = end - begin;
      out_stream << elapsed_seconds.count() << endl;
      out_stream << "T-matrix calculatred for excitation" << endl;
      MultiPrint(out_stream.str(), out_file_name, write_output);
      out_stream.str(string());

      // Calculating T-matrix for Raman
      begin = chrono::steady_clock::now();
      unique_ptr<stTmatrix<Real1>> T_mat_rm = slvForT(params_rm1, options);
      end = chrono::steady_clock::now();
      elapsed_seconds = end - begin;
      out_stream << elapsed_seconds.count() << endl;
      out_stream << "T-matrix calculatred for Raman" << endl;
      MultiPrint(out_stream.str(), out_file_name, write_output);
      out_stream.str(string());

      // Converting variables to secondary calculation type
      vector<unique_ptr<stTR<Real2>>> st_TR_list = ConvertStTRList<Real1, Real2>(T_mat->st_TR_list);
      vector<unique_ptr<stTR<Real2>>> st_TR_list_rm = ConvertStTRList<Real1, Real2>(T_mat_rm->st_TR_list);

      C_sca(k) = T_mat->st_C_oa->C_sca(0); // st_C_oa->sca should only contain 1 element
      C_ext(k) = T_mat->st_C_oa->C_ext(0); // st_C_oa->ext should only contain 1 element
      C_abs(k) = T_mat->st_C_oa->C_abs(0); // st_C_oa->abs should only contain 1 element
      out_stream.precision(10);
      out_stream << "C_sca " << C_sca(k) << endl;
      out_stream << "C_ext " << C_ext(k) << endl;
      out_stream << "C_abs " << C_abs(k) << endl;
      MultiPrint(out_stream.str(), out_file_name, write_output);
      out_stream.str(string());

      // Calculating scattering matrix for 2 scattering angles 0 and pi
      stSM_list[k] = pstScatteringMatrixOA(st_TR_list, params2->lambda(0), static_cast<Real2>(C_sca(k)));

      for (int t = 0; t < raman_params2->N_theta_p; t++) {
        begin = chrono::steady_clock::now();
        Real2 theta_p = raman_params2->theta_p_var(t);
        std::array<long int, 3> new_dims = {N_r, N_theta, 3};
        ArrayXr<Real2> phi_var = ArrayXr<Real2>::LinSpaced(N_phi + 1, 0, 2*PI);

        // Calculate field at excitation wavelength
        Real2 alpha_p = 0; // Defines the orientation of the electric field, in the plane orthogonal to wavevector k.
        params2->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        // Internal field calculations for the incident field defined by stIncPar
        unique_ptr<stAbcdnm<Real2>> st_abcdnm = rvhGetFieldCoefficients(N, st_TR_list, params2->inc_par);
        unique_ptr<stRes<Real2>> st_res_E = pstMakeStructForField(st_abcdnm, params2);
        ArrayXXc<Real2> c_nm = Map<RowArrayXc<Real2>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        ArrayXXc<Real2> d_nm = Map<RowArrayXc<Real2>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        unique_ptr<stEAllPhi<Real2>> st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda, st_res_E->epsilon2,
            c_nm, d_nm, r_row, theta_row, sBessel::J);
        // The field is calculated for a set of phi, from 0 and 2*pi
        Tensor4c<Real2> E_field_z(N_r, N_theta, N_phi + 1, 3);
        E_field_z.setZero();
        for (int m = 0; m <= N_phi; m++) {
          Real2 phi = phi_var(m);
          unique_ptr<stEforPhi<Real2>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->E_r.cols(), 3*st_E_for_phi->E_r.rows()};
          ArrayXXc<Real2> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->E_r.matrix().adjoint(),
              st_E_for_phi->E_t.matrix().adjoint(), st_E_for_phi->E_f.matrix().adjoint();
          E_field_z.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        alpha_p = PI/2; // Internal field calculations for different polarisation alpha_p
        params2->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        st_abcdnm = rvhGetFieldCoefficients(N, st_TR_list, params2->inc_par);
        st_res_E = pstMakeStructForField(st_abcdnm, params2);
        c_nm = Map<RowArrayXc<Real2>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        d_nm = Map<RowArrayXc<Real2>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda, st_res_E->epsilon2,
            c_nm, d_nm, r_row, theta_row, sBessel::J);
        // The field is calculated for a set of phi, from 0 and 2*pi
        Tensor4c<Real2> E_field_y(N_r, N_theta, N_phi + 1, 3);
        E_field_y.setZero();
        for (int m = 0; m <= N_phi; m++) {
          Real2 phi = phi_var(m);
          unique_ptr<stEforPhi<Real2>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->E_r.cols(), 3*st_E_for_phi->E_r.rows()};
          ArrayXXc<Real2> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->E_r.matrix().adjoint(),
              st_E_for_phi->E_t.matrix().adjoint(), st_E_for_phi->E_f.matrix().adjoint();
          E_field_y.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        // Calculate field at Raman wavelength
        alpha_p = 0;
        params_rm2->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        st_abcdnm = rvhGetFieldCoefficients(N, st_TR_list_rm, params_rm2->inc_par);
        st_res_E = pstMakeStructForField(st_abcdnm, params_rm2);
        c_nm = Map<RowArrayXc<Real2>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        d_nm = Map<RowArrayXc<Real2>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda, st_res_E->epsilon2,
            c_nm, d_nm, r_row, theta_row, sBessel::J);
        // The field is calculated for a set of phi, from 0 and 2*pi
        Tensor4c<Real2> E_field_rm_z(N_r, N_theta, N_phi + 1, 3);
        E_field_rm_z.setZero();
        for (int m = 0; m <= N_phi; m++) {
          Real2 phi = phi_var(m);
          unique_ptr<stEforPhi<Real2>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->E_r.cols(), 3*st_E_for_phi->E_r.rows()};
          ArrayXXc<Real2> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->E_r.matrix().adjoint(),
              st_E_for_phi->E_t.matrix().adjoint(), st_E_for_phi->E_f.matrix().adjoint();
          E_field_rm_z.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        alpha_p = PI/2;
        params_rm2->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        st_abcdnm = rvhGetFieldCoefficients(N, st_TR_list_rm, params_rm2->inc_par);
        st_res_E = pstMakeStructForField(st_abcdnm, params_rm2);
        c_nm = Map<RowArrayXc<Real2>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        d_nm = Map<RowArrayXc<Real2>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda, st_res_E->epsilon2,
            c_nm, d_nm, r_row, theta_row, sBessel::J);
        // The field is calculated for a set of phi, from 0 and 2*pi
        Tensor4c<Real2> E_field_rm_y(N_r, N_theta, N_phi + 1, 3);
        E_field_rm_y.setZero();
        for (int m = 0; m <= N_phi; m++) {
          Real2 phi = phi_var(m);
          unique_ptr<stEforPhi<Real2>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->E_r.cols(), 3*st_E_for_phi->E_r.rows()};
          ArrayXXc<Real2> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->E_r.matrix().adjoint(),
              st_E_for_phi->E_t.matrix().adjoint(), st_E_for_phi->E_f.matrix().adjoint();
          E_field_rm_y.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        // Calculate scattering cross-sections
        std::array<int, 1> dim3 = {3}, dim2 = {2};
        Tensor<Real2, 2> M_tens = 2*PI/N_phi*(E_field_rm_z * TensorConj(E_field_z))
            .sum(dim3).abs().pow(static_cast<Real2>(2.0)).sum(dim2);
        Map<ArrayXXr<Real2>> M_mat = ArrayMap(M_tens);
        ArrayXXr<Real2> F = M_mat * r_mat.pow(2) * sin(theta_mat);
        ArrayXXr<Real2> intergrand = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4 * dt_dr;
        sigma_zz(k, h_ind, t) = (intergrand).sum()/(PI*a*a*c*4/3);

        M_tens = 2*PI/N_phi*(E_field_rm_y * TensorConj(E_field_z))
            .sum(dim3).abs().pow(static_cast<Real2>(2.0)).sum(dim2);
        M_mat = ArrayMap(M_tens);
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        intergrand = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4 * dt_dr;
        sigma_yz(k, h_ind, t) = (intergrand).sum()/(PI*a*a*c*4/3);

        M_tens = 2*PI/N_phi*(E_field_rm_y * TensorConj(E_field_y))
            .sum(dim3).abs().pow(static_cast<Real2>(2.0)).sum(dim2);
        M_mat = ArrayMap(M_tens);
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        intergrand = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4 * dt_dr;
        sigma_yy(k, h_ind, t) = (intergrand).sum()/(PI*a*a*c*4/3);

        M_tens = 2*PI/N_phi*(E_field_rm_z * TensorConj(E_field_y))
            .sum(dim3).abs().pow(static_cast<Real2>(2.0)).sum(dim2);
        M_mat = ArrayMap(M_tens);
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        intergrand = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4 * dt_dr;
        sigma_zy(k, h_ind, t) = (intergrand).sum()/(PI*a*a*c*4/3);

        out_stream.precision(5);
        out_stream << "max diameter = " << max(a, c)*2e-3;
        out_stream.precision(10);
        out_stream << " µm, theta = " << theta_p * 180/PI << " degrees" << endl;
        out_stream.precision(11);
        out_stream << "--- relative raman zz = " << sigma_zz(k, h_ind, t) << endl;
        out_stream << "--- relative raman yz = " << sigma_yz(k, h_ind, t) << endl;
        out_stream << "--- relative raman zy = " << sigma_zy(k, h_ind, t) << endl;
        out_stream << "--- relative raman yy = " << sigma_yy(k, h_ind, t) << endl;

        end = chrono::steady_clock::now();
        elapsed_seconds = end - begin;
        out_stream << elapsed_seconds.count() << endl;

        MultiPrint(out_stream.str(), out_file_name, write_output);
        out_stream.str(string());
      }
    }
  }
}

#endif
