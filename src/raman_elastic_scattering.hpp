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

enum CalcType {SINGLE, DOUBLE, QUAD, CUSTOM, NONE};

template <class Real>
struct RamanParams {
  bool write_output = true;

  Real dia_min = 1000;
  Real dia_max = 2000;

  int N_rad = 100;
  int N_theta_p = 19;
  Real phi_p = 0;

  int N_r = 100;
  int N_theta = 320;
  int N_phi = 320;

  int N_h = 1;
  Real h_min = static_cast<Real>(1.0)/3;
  Real h_max = static_cast<Real>(1.0)/3;

  Real epsilon1 = 1;
  Real epsilon2 = pow(static_cast<Real>(1.35), 2);
  Real epsilon2_rm = pow(static_cast<Real>(1.344), 2);
  Real lambda = 355;
  Real lambda_rm = 403.7;

  int Nb_theta = 1000;
  int Nb_theta_pst = 1;

  ArrayXr<Real> rad_var;
  ArrayXr<Real> theta_p_var;
  ArrayXr<Real> h_var;
};

// defined in raman_elastic_scattering.cpp
std::array<CalcType, 2> GetCalcType(string in_file_name);
int GetNumCPUs(string in_file_name);
void MultiPrint(string out_string, string out_file_name, bool write_output = false);

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

template <class Real>
unique_ptr<RamanParams<Real>> LoadParams(string in_file_name) {
  ifstream in_file;
  in_file.open(in_file_name, ios::in);
  if (!in_file.is_open())
    throw runtime_error("Error: cannot open " + in_file_name);

  auto output = make_unique<RamanParams<Real>>();

  string line;
  vector<string> options {
    "Print output to file:",
    "Minimum diameter:", "Maximum diameter:", "Particle phi:", "Minimum h:", "Maximum h:",
    "epsilon1:", "epsilon2:", "Raman epsilon2:", "lambda:", "Raman lambda:",
    "No. of particle radii:", "No. of particle thetas:", "No. of h ratios:",
    "No. of r-coordinates:", "No. of theta-coordinates:", "No. of phi-coordinates:",
    "Nb_theta:", "Nb_theta_pst:"
  };
  std::array<bool, 1> boolParams = {output->write_output};
  std::array<Real, 10> floatParams = {
    output->dia_min, output->dia_max, output->phi_p, output->h_min, output->h_max,
    output->epsilon1, output->epsilon2, output->epsilon2_rm, output->lambda, output->lambda_rm
  };
  std::array<int, 8> intParams = {
    output->N_rad, output->N_theta_p, output->N_h, output->N_r, output->N_theta, output->N_phi,
    output->Nb_theta, output->Nb_theta_pst
  };
  for (size_t i = 0; i < options.size(); i++) {
    string option = options[i];
    do {
      if (in_file.peek() == EOF)
        throw runtime_error("Error: cannot find option " + option);
      getline(in_file, line);
    } while (line.find(option) == string::npos);
    in_file.clear();
    in_file.seekg(0, in_file.beg);

    try {
      string param = line.substr(line.find(option) + option.size());
      if (i < boolParams.size())
        boolParams[i] = (param == "yes");
      else if (i < boolParams.size() + floatParams.size())
        floatParams[i - boolParams.size()] = Frac2Float<Real>(param);
      else if (i < boolParams.size() + floatParams.size() + intParams.size())
        intParams[i - boolParams.size() - floatParams.size()] = lexical_cast<int>(param);
    } catch(bad_lexical_cast &) {
      cerr << "Cannot read value of option \"" << option << "\". Using default value.";
    }

    output->write_output = boolParams[0];
    output->dia_min = floatParams[0];
    output->dia_max = floatParams[1];
    output->phi_p = floatParams[2];
    output->h_min = floatParams[3];
    output->h_max = floatParams[4];
    output->epsilon1 = floatParams[5];
    output->epsilon2 = floatParams[6];
    output->epsilon2_rm = floatParams[7];
    output->lambda = floatParams[8];
    output->lambda_rm = floatParams[9];
    output->N_rad = intParams[0];
    output->N_theta_p = intParams[1];
    output->N_h = intParams[2];
    output->N_r = intParams[3];
    output->N_theta = intParams[4];
    output->N_phi = intParams[5];
    output->Nb_theta = intParams[6];
    output->Nb_theta_pst = intParams[7];
  }
  in_file.close();
  ArrayXr<Real> dia_var = ArrayXr<Real>::LinSpaced(output->N_rad, output->dia_min, output->dia_max);
  output->rad_var = dia_var/2;
  output->theta_p_var = ArrayXr<Real>::LinSpaced(output->N_theta_p, 0, mp_pi<Real>()/2);
  output->h_var = ArrayXr<Real>::LinSpaced(output->N_h, output->h_min, output->h_max);
  return output;
}

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
  params->N = 6 + 2*static_cast<int>(ceil(max(params->a, params->c)/40));
  return params;
}

// Can make this into a non-template function using c++20 'concepts' keyword
template <class Real1, class Real2>
void CreateTimeStamp(string out_file_name, const unique_ptr<RamanParams<Real1>>& raman_params) {
  string main_calc_type = getTypeName<Real1>();
  string second_calc_type = getTypeName<Real2>();

  ofstream out_file;
  out_file.open(out_file_name, ios::out | ios::app);
  auto sys_time = chrono::system_clock::now();
  time_t sys_time_date = chrono::system_clock::to_time_t(sys_time);
  out_file << "Session began at system time: " << ctime(&sys_time_date);
  out_file << "Running RamanElasticScattering with type " << main_calc_type;
  out_file << " and " << second_calc_type << endl;
  out_file.flush();
  out_file.close();
}

template <class From, class To>
vector<unique_ptr<stTR<To>>> ConvertstTRList(const vector<unique_ptr<stTR<From>>>& st_TR_list) {
  vector<unique_ptr<stTR<To>>> output(st_TR_list.size());
  complex<To> I = mp_im_unit<To>();
  for (size_t i = 0; i < st_TR_list.size(); i++) {
    assert(st_TR_list[i]);
    assert(st_TR_list[i]->mat_list.size() == 2);
    assert(st_TR_list[i]->mat_list[0] == "st_4M_T");
    assert(st_TR_list[i]->mat_list[1] == "st_4M_R");
    output[i] = make_unique<stTR<To>>();
    output[i]->mat_list = st_TR_list[i]->mat_list;

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

template <class Real1, class Real2>
void RamanElasticScattering(string in_file_name, string out_dir) {
  Real2 PI = mp_pi<Real2>();

  unique_ptr<RamanParams<Real1>> raman_params1 = LoadParams<Real1>(in_file_name);
  unique_ptr<RamanParams<Real2>> raman_params2 = LoadParams<Real2>(in_file_name);
  bool write_output = raman_params1->write_output;

  Tensor3r<Real2> sigma_yz(raman_params2->N_rad, raman_params2->N_h, raman_params2->N_theta_p);
  Tensor3r<Real2> sigma_zy(raman_params2->N_rad, raman_params2->N_h, raman_params2->N_theta_p);
  Tensor3r<Real2> sigma_zz(raman_params2->N_rad, raman_params2->N_h, raman_params2->N_theta_p);
  Tensor3r<Real2> sigma_yy(raman_params2->N_rad, raman_params2->N_h, raman_params2->N_theta_p);
  ArrayXr<Real1> sca = ArrayXr<Real1>::Zero(raman_params1->N_rad);
  ArrayXr<Real1> ext = ArrayXr<Real1>::Zero(raman_params1->N_rad);
  ArrayXr<Real1> abs = ArrayXr<Real1>::Zero(raman_params1->N_rad);
  vector<unique_ptr<stSM<Real2>>> stSM_list(raman_params2->N_rad);

  auto options = make_unique<stOptions>();
  options->get_R = true;
  options->delta = 0;
  options->NB = 0;
  options->get_symmetric_T = false;
  Real2 phi_p = raman_params2->phi_p;

  int N_r = raman_params2->N_r;
  int N_theta = raman_params2->N_theta;
  int N_phi = raman_params2->N_phi;

  ArrayXr<Real2> theta_surf = ArrayXr<Real2>::LinSpaced(N_theta, 0, PI);
  RowArrayXr<Real2> r_surf_u = RowArrayXr<Real2>::LinSpaced(
      N_r, 1/static_cast<Real2>(N_r), 1).pow(static_cast<Real2>(1.0)/3);
  ArrayXXr<Real2> r_mat_u = r_surf_u.replicate(N_theta, 1).transpose();
  ArrayXXr<Real2> theta_mat = theta_surf.replicate(1, N_r).transpose();

  if (write_output)
    std::filesystem::create_directory(out_dir);

  for (int h_ind = 0; h_ind < raman_params1->N_h; h_ind++) {
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

      ArrayXXr<Real2> r_mat = a*c*r_mat_u/sqrt(pow(c, 2)*sin(theta_mat).pow(2) + pow(a, 2)*cos(theta_mat).pow(2));
      ArrayXXr<Real2> d_theta = theta_mat(all, seq(1, last)) - theta_mat(all, seq(0, last - 1));
      ArrayXXr<Real2> dr = r_mat(seq(1, last), all) - r_mat(seq(0, last - 1), all);
      ArrayXXr<Real2> dt_dr = d_theta(seq(1, last), all) * dr(all, seq(1, last));
      RowArrayXr<Real2> r_row = r_mat.reshaped().transpose();
      RowArrayXr<Real2> theta_row = theta_mat.reshaped().transpose();

      chrono::steady_clock::time_point begin, end;
      chrono::duration<double> elapsed_seconds;

      begin = chrono::steady_clock::now();
      unique_ptr<stTmatrix<Real1>> T_mat = slvForT(params1, options);
      end = chrono::steady_clock::now();
      elapsed_seconds = end - begin;
      out_stream << elapsed_seconds.count() << endl;
      out_stream << "T-matrix calculatred for excitation" << endl;
      MultiPrint(out_stream.str(), out_file_name, write_output);
      out_stream.str(string());

      begin = chrono::steady_clock::now();
      unique_ptr<stTmatrix<Real1>> T_mat_rm = slvForT(params_rm1, options);
      end = chrono::steady_clock::now();
      elapsed_seconds = end - begin;
      out_stream << elapsed_seconds.count() << endl;
      out_stream << "T-matrix calculatred for Raman" << endl;
      MultiPrint(out_stream.str(), out_file_name, write_output);
      out_stream.str(string());

      vector<unique_ptr<stTR<Real2>>> st_TR_list = ConvertstTRList<Real1, Real2>(T_mat->st_TR_list);
      vector<unique_ptr<stTR<Real2>>> st_TR_list_rm = ConvertstTRList<Real1, Real2>(T_mat_rm->st_TR_list);

      sca(k) = T_mat->st_coa->sca(0); // st_coa->sca should only contain 1 element
      ext(k) = T_mat->st_coa->ext(0); // st_coa->ext should only contain 1 element
      abs(k) = T_mat->st_coa->abs(0); // st_coa->abs should only contain 1 element
      out_stream.precision(6);
      out_stream << "Csca " << sca(k) << endl;
      out_stream << "Cext " << ext(k) << endl;
      out_stream << "Cabs " << abs(k) << endl;
      MultiPrint(out_stream.str(), out_file_name, write_output);
      out_stream.str(string());

      stSM_list[k] = pstScatteringMatrixOA(st_TR_list, params2->lambda(0), static_cast<Real2>(sca(k)));

      for (int t = 0; t < raman_params2->N_theta_p; t++) {
        begin = chrono::steady_clock::now();
        Real2 theta_p = raman_params2->theta_p_var(t);
        std::array<long int, 3> new_dims = {N_r, N_theta, 3};
        ArrayXr<Real2> phi_var = ArrayXr<Real2>::LinSpaced(N_phi + 1, 0, 2*PI);

        Real2 alpha_p = 0;
        params2->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        unique_ptr<stAbcdnm<Real2>> st_abcdnm = rvhGetFieldCoefficients(N, st_TR_list, params2->inc_par);
        unique_ptr<stRes<Real2>> st_res_E = pstMakeStructForField(st_abcdnm, params2);
        ArrayXXc<Real2> c_nm = Map<RowArrayXc<Real2>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        ArrayXXc<Real2> d_nm = Map<RowArrayXc<Real2>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        unique_ptr<stEAllPhi<Real2>> st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda,
            st_res_E->epsilon2, c_nm, d_nm, r_row, theta_row, sBessel::J);
        Tensor4c<Real2> E_field_z(N_r, N_theta, N_phi + 1, 3);
        E_field_z.setZero();
        for (int m = 0; m <= N_phi; m++) {
          Real2 phi = phi_var(m);
          unique_ptr<stEforPhi<Real2>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->Er.cols(), 3*st_E_for_phi->Er.rows()};
          ArrayXXc<Real2> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->Er.matrix().adjoint(),
              st_E_for_phi->Et.matrix().adjoint(), st_E_for_phi->Ef.matrix().adjoint();
          E_field_z.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        alpha_p = PI/2;
        params2->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        st_abcdnm = rvhGetFieldCoefficients(N, st_TR_list, params2->inc_par);
        st_res_E = pstMakeStructForField(st_abcdnm, params2);
        c_nm = Map<RowArrayXc<Real2>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        d_nm = Map<RowArrayXc<Real2>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda, st_res_E->epsilon2,
            c_nm, d_nm, r_row, theta_row, sBessel::J);
        Tensor4c<Real2> E_field_y(N_r, N_theta, N_phi + 1, 3);
        E_field_y.setZero();
        for (int m = 0; m <= N_phi; m++) {
          Real2 phi = phi_var(m);
          unique_ptr<stEforPhi<Real2>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->Er.cols(), 3*st_E_for_phi->Er.rows()};
          ArrayXXc<Real2> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->Er.matrix().adjoint(),
              st_E_for_phi->Et.matrix().adjoint(), st_E_for_phi->Ef.matrix().adjoint();
          E_field_y.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        alpha_p = 0;
        params_rm2->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        st_abcdnm = rvhGetFieldCoefficients(N, st_TR_list_rm, params_rm2->inc_par);
        st_res_E = pstMakeStructForField(st_abcdnm, params_rm2);
        c_nm = Map<RowArrayXc<Real2>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        d_nm = Map<RowArrayXc<Real2>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda, st_res_E->epsilon2,
            c_nm, d_nm, r_row, theta_row, sBessel::J);
        Tensor4c<Real2> E_field_rm_z(N_r, N_theta, N_phi + 1, 3);
        E_field_rm_z.setZero();
        for (int m = 0; m <= N_phi; m++) {
          Real2 phi = phi_var(m);
          unique_ptr<stEforPhi<Real2>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->Er.cols(), 3*st_E_for_phi->Er.rows()};
          ArrayXXc<Real2> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->Er.matrix().adjoint(),
              st_E_for_phi->Et.matrix().adjoint(), st_E_for_phi->Ef.matrix().adjoint();
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
        Tensor4c<Real2> E_field_rm_y(N_r, N_theta, N_phi + 1, 3);
        E_field_rm_y.setZero();
        for (int m = 0; m <= N_phi; m++) {
          Real2 phi = phi_var(m);
          unique_ptr<stEforPhi<Real2>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->Er.cols(), 3*st_E_for_phi->Er.rows()};
          ArrayXXc<Real2> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->Er.matrix().adjoint(),
              st_E_for_phi->Et.matrix().adjoint(), st_E_for_phi->Ef.matrix().adjoint();
          E_field_rm_y.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        std::array<int, 1> dim3 = {3}, dim2 = {2};
        Tensor<Real2, 2> inter_step = 2*PI/N_phi*(E_field_rm_z * tensor_conj(E_field_z))
            .sum(dim3).abs().pow(static_cast<Real2>(2.0)).sum(dim2);
        ArrayXXr<Real2> M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        ArrayXXr<Real2> F = M_mat * r_mat.pow(2) * sin(theta_mat);
        ArrayXXr<Real2> tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_zz(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real2>(4.0)/3*PI*a*a*c);

        inter_step = 2*PI/N_phi*(E_field_rm_y * tensor_conj(E_field_z))
            .sum(dim3).abs().pow(static_cast<Real2>(2.0)).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_yz(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real2>(4.0)/3*PI*a*a*c);

        inter_step = 2*PI/N_phi*(E_field_rm_y * tensor_conj(E_field_y))
            .sum(dim3).abs().pow(static_cast<Real2>(2.0)).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_yy(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real2>(4.0)/3*PI*a*a*c);

        inter_step = 2*PI/N_phi*(E_field_rm_z * tensor_conj(E_field_y))
            .sum(dim3).abs().pow(static_cast<Real2>(2.0)).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_zy(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real2>(4.0)/3*PI*a*a*c);

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
