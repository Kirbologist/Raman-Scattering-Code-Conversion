#ifndef RAMAN_ELASTIC_SCATTERING_HPP
#define RAMAN_ELASTIC_SCATTERING_HPP

#include "core.hpp"
#include "smarties.hpp"
#include <chrono>
#include <fstream>
#include <exception>

using namespace std;
using namespace Eigen;
using namespace Smarties;

enum CalcType {SINGLE, DOUBLE, QUAD, CUSTOM, NONE};

template <class Real>
struct RamanParams {
  Real dia_min = 1000;
  Real dia_max = 2000;
  int N_rad = 100;
  int N_theta_p = 19;
  ArrayXr<Real> rad_var;
  ArrayXr<Real> theta_p_var;
  bool write_output = false;
  bool write_log = true;
};

// defined in raman_elastic_scattering.cpp
std::array<CalcType, 2> GetCalcType(string in_file_name);
int GetNumCPUs(string in_file_name);
void MultiPrint(string out_string, string out_file_name, bool write_output = false);

template <class Real>
unique_ptr<stParams<Real>> loadParam(string type = string()) {
  auto params = make_unique<stParams<Real>>();
  params->epsilon1 = 1;
  params->s = ArrayXr<Real>(1);
  if (type == "rm") {
    params->lambda = {{403.7}};
    params->epsilon2 = {{pow(1.344, 2)}};
    params->k1 = 2*mp_pi<Real>()/params->lambda;
    params->s(0) = 1.344;
  } else {
    params->lambda = {{355}};
    params->epsilon2 = {{pow(1.35, 2)}};
    params->k1 = 2*mp_pi<Real>()/params->lambda;
    params->s(0) = 1.35;
  }
  return params;
}

template <class Real>
unique_ptr<RamanParams<Real>> GetRamanParams(string in_file_name) {
  ifstream in_file;
  in_file.open(in_file_name, ios::in);
  if (!in_file.is_open())
    throw runtime_error("Error: cannot open config.txt");

  auto output = make_unique<RamanParams<Real>>();

  string line;
  vector<string> options {
    "Minimum diameter:", "Maximum diameter:", "No. of radii:", "No. of thetas:",
    "Print output to file:", "Print log to file:"
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

    string param = line.substr(line.find(option) + option.size());
    switch (i) {
      case 0 : {
        output->dia_min = (param.size() > 0 && isdigit(param[0])) ? stof(param, nullptr) : output->dia_min;
        break;
      }
      case 1 : {
        output->dia_max = (param.size() > 0 && isdigit(param[0])) ? stof(param, nullptr) : output->dia_max;
        break;
      }
      case 2 : {
        output->N_rad = (param.size() > 0 && isdigit(param[0])) ? stoi(param, nullptr) : output->N_rad;
        break;
      }
      case 3 : {
        output->N_theta_p = (param.size() > 0 && isdigit(param[0])) ? stoi(param, nullptr) : output->N_theta_p;
        break;
      }
      case 4 : {
        output->write_output = (param == "yes");
        break;
      }
      case 5 : {
        output->write_log = (param == "yes");
        break;
      }
    }
  }
  in_file.close();
  ArrayXr<Real> dia_var = ArrayXr<Real>::LinSpaced(output->N_rad, output->dia_min, output->dia_max);
  output->rad_var = dia_var/2;
  output->theta_p_var = ArrayXr<Real>::LinSpaced(output->N_theta_p, 0, mp_pi<Real>()/2);
  return output;
}

// Can make this into a non-template function using c++20 'concepts' keyword
template <class Real1, class Real2>
void CreateTimeStamp(string log_file_name, const unique_ptr<RamanParams<Real1>>& raman_params, bool is_mixed) {
  string main_calc_type = typeid(raman_params->dia_min).name(), second_calc_type = typeid(mp_pi<Real2>()).name();
  string converted_main_calc_type, converted_second_calc_type;
  if (main_calc_type.find("N5boost14multiprecision6numberINS0_8backends18mpfr_float_backend") != string::npos) {
    size_t prec_begin = main_calc_type.find("ILj") + 3;
    size_t prec_end = main_calc_type.find("ELNS");
    string precision = main_calc_type.substr(prec_begin, prec_end - prec_begin);
    converted_main_calc_type = "custom with precision of " + precision + " bits";
  } else if (main_calc_type == "f")
    converted_main_calc_type = "single";
  else if (main_calc_type == "d")
    converted_main_calc_type = "double";
  else if (main_calc_type == "e")
    converted_main_calc_type = "quad";

  if (second_calc_type.find("N5boost14multiprecision6numberINS0_8backends18mpfr_float_backend") != string::npos) {
    size_t prec_begin = second_calc_type.find("ILj") + 3;
    size_t prec_end = second_calc_type.find("ELNS");
    string precision = second_calc_type.substr(prec_begin, prec_end - prec_begin);
    converted_second_calc_type = "custom with precision of " + precision + " bits";
  } else if (second_calc_type == "f")
    converted_second_calc_type = "single";
  else if (second_calc_type == "d")
    converted_second_calc_type = "double";
  else if (second_calc_type == "e")
    converted_second_calc_type = "quad";

  ofstream log_file;
  log_file.open(log_file_name, ios::out | ios::app);
  auto sys_time = chrono::system_clock::now();
  time_t sys_time_date = chrono::system_clock::to_time_t(sys_time);
  log_file << "Session began at system time: " << ctime(&sys_time_date);
  log_file << "Running RamanElasticScattering with type " << converted_main_calc_type;
  log_file << " and " << converted_second_calc_type << endl;
  log_file << "Parameters are DIA_MIN: " << raman_params->dia_min << ", DIA_MAX: " << raman_params->dia_max <<
      ", N_RAD: " << raman_params->N_rad << ", N_THETA_P: " << raman_params->N_theta_p << endl;
  log_file.flush();
  log_file.close();
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
void RamanElasticScattering(string in_file_name, string log_file_name) {
  Real2 PI = mp_pi<Real2>();
  unique_ptr<RamanParams<Real1>> raman_params1 = GetRamanParams<Real1>(in_file_name);
  unique_ptr<RamanParams<Real2>> raman_params2 = GetRamanParams<Real2>(in_file_name);
  ArrayXr<Real1> rad_var1 = raman_params1->rad_var;
  ArrayXr<Real2> rad_var2 = raman_params2->rad_var;
  ArrayXr<Real2> theta_p_var = raman_params2->theta_p_var;
  bool write_output = raman_params1->write_output;
  bool write_log = raman_params1->write_log;
  ArrayXr<Real1> h_var1 = {{static_cast<Real1>(1.0)/3}};
  ArrayXr<Real2> h_var2 = {{static_cast<Real2>(1.0)/3}};

  ofstream log_file;
  if (write_log) {
    log_file.open(log_file_name, ios::out | ios::app);
    log_file << endl;
    if (!log_file.is_open()) {
      cerr << "Warning: cannot create or open log file. Some output won't be written." << endl;
      log_file.close();
      write_log = false;
    }
  }
  if (write_log)
    CreateTimeStamp<Real1, Real2>(log_file_name, raman_params1, false);

  Tensor3r<Real2> sigma_yz(rad_var2.size(), h_var2.size(), theta_p_var.size());
  Tensor3r<Real2> sigma_zy(rad_var2.size(), h_var2.size(), theta_p_var.size());
  Tensor3r<Real2> sigma_zz(rad_var2.size(), h_var2.size(), theta_p_var.size());
  Tensor3r<Real2> sigma_yy(rad_var2.size(), h_var2.size(), theta_p_var.size());
  ArrayXr<Real1> sca = ArrayXr<Real1>::Zero(rad_var1.size());
  ArrayXr<Real1> ext = ArrayXr<Real1>::Zero(rad_var1.size());
  ArrayXr<Real1> abs = ArrayXr<Real1>::Zero(rad_var1.size());
  vector<unique_ptr<stSM<Real2>>> stSM_list(rad_var2.size());

  auto options = make_unique<stOptions>();
  options->get_R = true;
  options->delta = 0;
  options->NB = 0;
  options->get_symmetric_T = false;
  Real2 phi_p = 0;

  int N_r = 100;
  int N_theta = 320;
  int N_phi = 320;

  ArrayXr<Real2> theta_surf = ArrayXr<Real2>::LinSpaced(N_theta, 0, PI);
  RowArrayXr<Real2> r_surf_u = RowArrayXr<Real2>::LinSpaced(
      N_r, 1/static_cast<Real2>(N_r), 1).pow(static_cast<Real2>(1.0)/3);
  ArrayXXr<Real2> r_mat_u = r_surf_u.replicate(N_theta, 1).transpose();
  ArrayXXr<Real2> theta_mat = theta_surf.replicate(1, N_r).transpose();

  int Nb_theta = 1000;
  int Nb_theta_pst = 1;

  for (int h_ind = 0; h_ind < h_var1.size(); h_ind++) {
    stringstream aspect_ratio_string;
    Real1 h1 = h_var1(h_ind);
    Real2 h2 = h_var2(h_ind);

    aspect_ratio_string.precision(4);
    aspect_ratio_string << "aspect ratio " << h1 << endl;
    MultiPrint(aspect_ratio_string.str(), log_file_name, write_log);
    aspect_ratio_string.str(string());

    #pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < rad_var1.size(); k++) {
      unique_ptr<stParams<Real1>> params1 = loadParam<Real1>();
      unique_ptr<stParams<Real1>> params_rm1 = loadParam<Real1>("rm");
      unique_ptr<stParams<Real2>> params2 = loadParam<Real2>();
      unique_ptr<stParams<Real2>> params_rm2 = loadParam<Real2>("rm");
      params1->Nb_theta = Nb_theta;
      params1->Nb_theta_pst = Nb_theta_pst;
      params_rm1->Nb_theta = Nb_theta;
      params_rm1->Nb_theta_pst = Nb_theta_pst;
      params2->Nb_theta = Nb_theta;
      params2->Nb_theta_pst = Nb_theta_pst;
      params_rm2->Nb_theta = Nb_theta;
      params_rm2->Nb_theta_pst = Nb_theta_pst;
      stringstream log_stream;

      Real1 a1, c1;
      if (h1 > 1) {
        a1 = rad_var1(k);
        c1 = rad_var1(k) / h1;
      } else {
        a1 = rad_var1(k) * h1;
        c1 = rad_var1(k);
      }
      Real2 a2, c2;
      if (h2 > 1) {
        a2 = rad_var2(k);
        c2 = rad_var2(k) / h2;
      } else {
        a2 = rad_var2(k) * h2;
        c2 = rad_var2(k);
      }

      int N = 6 + 2*static_cast<int>(ceil(max(a1, c1)/40));
      params1->N = N;
      params1->a = a1;
      params1->c = c1;
      params_rm1->N = N;
      params_rm1->a = a1;
      params_rm1->c = c1;
      params2->N = N;
      params2->a = a2;
      params2->c = c2;
      params_rm2->N = N;
      params_rm2->a = a2;
      params_rm2->c = c2;

      ArrayXXr<Real2> r_mat = a2*c2*r_mat_u/sqrt(pow(c2, 2)*sin(theta_mat).pow(2) + pow(a2, 2)*cos(theta_mat).pow(2));
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
      log_stream << elapsed_seconds.count() << endl;
      log_stream << "T-matrix calculatred for excitation" << endl;
      MultiPrint(log_stream.str(), log_file_name, write_log);
      log_stream.str(string());

      begin = chrono::steady_clock::now();
      unique_ptr<stTmatrix<Real1>> T_mat_rm = slvForT(params_rm1, options);
      end = chrono::steady_clock::now();
      elapsed_seconds = end - begin;
      log_stream << elapsed_seconds.count() << endl;
      log_stream << "T-matrix calculatred for Raman" << endl;
      MultiPrint(log_stream.str(), log_file_name, write_log);
      log_stream.str(string());

      vector<unique_ptr<stTR<Real2>>> st_TR_list = ConvertstTRList<Real1, Real2>(T_mat->st_TR_list);
      vector<unique_ptr<stTR<Real2>>> st_TR_list_rm = ConvertstTRList<Real1, Real2>(T_mat_rm->st_TR_list);

      sca(k) = T_mat->st_coa->sca(0); // st_coa->sca should only contain 1 element
      ext(k) = T_mat->st_coa->ext(0); // st_coa->ext should only contain 1 element
      abs(k) = T_mat->st_coa->abs(0); // st_coa->abs should only contain 1 element
      log_stream << "Csca " << sca(k) << endl;
      log_stream << "Cext " << ext(k) << endl;
      log_stream << "Cabs " << abs(k) << endl;
      MultiPrint(log_stream.str(), log_file_name, write_log);
      log_stream.str(string());

      stSM_list[k] = pstScatteringMatrixOA(st_TR_list, params2->lambda(0), static_cast<Real2>(sca(k)));

      for (int t = 0; t < theta_p_var.size(); t++) {
        begin = chrono::steady_clock::now();
        Real2 theta_p = theta_p_var(t);
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
        sigma_zz(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real2>(4.0)/3*PI*a2*a2*c2);

        inter_step = 2*PI/N_phi*(E_field_rm_y * tensor_conj(E_field_z))
            .sum(dim3).abs().pow(static_cast<Real2>(2.0)).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_yz(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real2>(4.0)/3*PI*a2*a2*c2);

        inter_step = 2*PI/N_phi*(E_field_rm_y * tensor_conj(E_field_y))
            .sum(dim3).abs().pow(static_cast<Real2>(2.0)).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_yy(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real2>(4.0)/3*PI*a2*a2*c2);

        inter_step = 2*PI/N_phi*(E_field_rm_z * tensor_conj(E_field_y))
            .sum(dim3).abs().pow(static_cast<Real2>(2.0)).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_zy(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real2>(4.0)/3*PI*a2*a2*c2);

        log_stream.precision(5);
        log_stream << "max diameter = " << max(a2, c2)*2e-3;
        log_stream.precision(10);
        log_stream << " µm, theta = " << theta_p * 180/PI << " degrees" << endl;
        log_stream.precision(11);
        log_stream << "--- relative raman zz = " << sigma_zz(k, h_ind, t) << endl;
        log_stream << "--- relative raman yz = " << sigma_yz(k, h_ind, t) << endl;
        log_stream << "--- relative raman zy = " << sigma_zy(k, h_ind, t) << endl;
        log_stream << "--- relative raman yy = " << sigma_yy(k, h_ind, t) << endl;

        end = chrono::steady_clock::now();
        elapsed_seconds = end - begin;
        log_stream << elapsed_seconds.count() << endl;

        MultiPrint(log_stream.str(), log_file_name, write_log);
        log_stream.str(string());
      }
    }
  }
}

#endif
