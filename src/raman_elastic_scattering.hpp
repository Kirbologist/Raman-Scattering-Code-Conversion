#ifndef RAMAN_ELASTIC_SCATTERING_HPP
#define RAMAN_ELASTIC_SCATTERING_HPP

#include "core.hpp"
#include "smarties.hpp"
#include <chrono>
#include <fstream>
#include <exception>
#include <cassert>

using namespace std;
using namespace Eigen;
using namespace Smarties;

template <class Real>
struct RamanParams {
  int cpu_n = 0;
  int cpus = 1;
  Real dia_min = 1000;
  Real dia_max = 2000;
  int N_rad = 100;
  int N_theta_p = 19;
  ArrayXr<Real> rad_var;
  ArrayXr<Real> theta_p_var;
  bool write_output = false;
};

// defined in defs_default.cpp
string GetCalcType(string in_file_name);
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
    "CPU no.:", "No. of CPUs:", "Minimum diameter:", "Maximum diameter:",
    "No. of radii:", "No. of thetas:", "Print output to file:"
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
        output->cpu_n = (param.size() > 0 && isdigit(param[0])) ? stoi(param, nullptr) : output->cpu_n;
        break;
      }
      case 1 : {
        output->cpus = (param.size() > 0 && isdigit(param[0])) ? stoi(param, nullptr) : output->cpus;
        break;
      }
      case 2 : {
        output->dia_min = (param.size() > 0 && isdigit(param[0])) ? stof(param, nullptr) : output->dia_min;
        break;
      }
      case 3 : {
        output->dia_max = (param.size() > 0 && isdigit(param[0])) ? stof(param, nullptr) : output->dia_max;
        break;
      }
      case 4 : {
        output->N_rad = (param.size() > 0 && isdigit(param[0])) ? stoi(param, nullptr) : output->N_rad;
        break;
      }
      case 5 : {
        output->N_theta_p = (param.size() > 0 && isdigit(param[0])) ? stoi(param, nullptr) : output->N_theta_p;
        break;
      }
      case 6 : {
        output->write_output = (param == "yes");
        break;
      }
    }
  }
  in_file.close();
  ArrayXr<Real> par = ArrayXr<Real>::LinSpaced(output->N_rad, output->dia_min, output->dia_max);
  int par_per_cpu = output->N_rad / output->cpus;
  ArrayXr<Real> dia_var = par(ArrayXi::LinSpaced(par_per_cpu, 0, par_per_cpu - 1) + output->cpu_n * par_per_cpu);
  output->rad_var = dia_var/2;
  output->theta_p_var = ArrayXr<Real>::LinSpaced(output->N_theta_p, 0, mp_pi<Real>()/2);
  return output;
}

// Can make this into a non-template function using c++20 'concepts' keyword
template <class Real>
void CreateTimeStamp(string out_file_name, const unique_ptr<RamanParams<Real>>& raman_params, bool is_mixed) {
  string main_calc_type = typeid(raman_params->dia_min).name();
  string converted_calc_type;
  if (main_calc_type.find("N5boost14multiprecision6numberINS0_8backends18mpfr_float_backend") != string::npos) {
    size_t prec_begin = main_calc_type.find("ILj") + 3;
    size_t prec_end = main_calc_type.find("ELNS");
    string precision = main_calc_type.substr(prec_begin, prec_end - prec_begin);
    converted_calc_type = "Boost MPFR with precision of " + precision + " bits";
  } else if (main_calc_type == "d")
    converted_calc_type = "double";
  else if (main_calc_type == "e")
    converted_calc_type = "long double";

  ofstream out_file;
  out_file.open(out_file_name, ios::out | ios::app);
  auto sys_time = chrono::system_clock::now();
  time_t sys_time_date = chrono::system_clock::to_time_t(sys_time);
  out_file << "Session began at system time: " << ctime(&sys_time_date);
  out_file << "Running RamanElasticScattering with type " << converted_calc_type;
  if (is_mixed)
    out_file << " and type double";
  out_file << endl;
  out_file << "Parameters are CPU_N: " << raman_params->cpu_n << ", CPUS: " << raman_params->cpus <<
      ", DIA_MIN: " << raman_params->dia_min << ", DIA_MAX: " << raman_params->dia_max <<
      ", N_RAD: " << raman_params->N_rad << ", N_THETA_P: " << raman_params->N_theta_p << endl;
  out_file.flush();
  out_file.close();
}

template <class Real>
void RamanElasticScattering(string in_file_name, string out_file_name) {
  Real PI = mp_pi<Real>();
  unique_ptr<RamanParams<Real>> raman_params = GetRamanParams<Real>(in_file_name);
  ArrayXr<Real> rad_var = raman_params->rad_var;
  ArrayXr<Real> theta_p_var = raman_params->theta_p_var;
  bool write_output = raman_params->write_output;
  ArrayXr<Real> h_var = {{static_cast<Real>(1.0)/3}};

  ofstream out_file;
  if (write_output) {
    out_file.open(out_file_name, ios::out | ios::app);
    out_file << endl;
    if (!out_file.is_open()) {
      cerr << "Warning: cannot create or open output file. Some output won't be written." << endl;
      out_file.close();
      write_output = false;
    }
  }
  if (write_output)
    CreateTimeStamp(out_file_name, raman_params, false);

  Tensor3r<Real> sigma_yz(rad_var.size(), h_var.size(), theta_p_var.size());
  Tensor3r<Real> sigma_zy(rad_var.size(), h_var.size(), theta_p_var.size());
  Tensor3r<Real> sigma_zz(rad_var.size(), h_var.size(), theta_p_var.size());
  Tensor3r<Real> sigma_yy(rad_var.size(), h_var.size(), theta_p_var.size());
  ArrayXr<Real> sca = ArrayXr<Real>::Zero(rad_var.size());
  ArrayXr<Real> ext = ArrayXr<Real>::Zero(rad_var.size());
  ArrayXr<Real> abs = ArrayXr<Real>::Zero(rad_var.size());
  vector<unique_ptr<stSM<Real>>> stSM_list(rad_var.size());

  auto options = make_unique<stOptions>();
  options->get_R = true;
  options->delta = 0;
  options->NB = 0;
  options->get_symmetric_T = false;
  unique_ptr<stParams<Real>> params = loadParam<Real>();
  unique_ptr<stParams<Real>> params_rm = loadParam<Real>("rm");
  Real phi_p = 0;

  int N_r = 100;
  int N_theta = 320;
  int N_phi = 320;

  ArrayXr<Real> theta_surf = ArrayXr<Real>::LinSpaced(N_theta, 0, PI);
  RowArrayXr<Real> r_surf_u = RowArrayXr<Real>::LinSpaced(
      N_r, 1/static_cast<Real>(N_r), 1).pow(static_cast<Real>(1.0)/3);
  ArrayXXr<Real> r_mat_u = r_surf_u.replicate(N_theta, 1).transpose();
  ArrayXXr<Real> theta_mat = theta_surf.replicate(1, N_r).transpose();

  int Nb_theta = 1000;
  int Nb_theta_pst = 1;
  params->Nb_theta = Nb_theta;
  params->Nb_theta_pst = Nb_theta_pst;
  params_rm->Nb_theta = Nb_theta;
  params_rm->Nb_theta_pst = Nb_theta_pst;

  for (int h_ind = 0; h_ind < h_var.size(); h_ind++) {
    stringstream out_stream;
    Real h = h_var(h_ind);

    out_stream.precision(4);
    out_stream << "aspect ratio " << h_var << endl;
    MultiPrint(out_stream.str(), out_file_name, write_output);
    out_stream.str(string());

    for (int k = 0; k < rad_var.size(); k++) {
      Real a, c;
      if (h > 1) {
        a = rad_var(k);
        c = rad_var(k) / h;
      } else {
        a = rad_var(k) * h;
        c = rad_var(k);
      }
      int N = 6 + 2*static_cast<int>(ceil(max(a, c)/40));
      params->N = N;
      params->a = a;
      params->c = c;
      params_rm->N = N;
      params_rm->a = a;
      params_rm->c = c;

      ArrayXXr<Real> r_mat = a*c*r_mat_u/sqrt(pow(c, 2)*sin(theta_mat).pow(2) + pow(a, 2)*cos(theta_mat).pow(2));
      ArrayXXr<Real> d_theta = theta_mat(all, seq(1, last)) - theta_mat(all, seq(0, last - 1));
      ArrayXXr<Real> dr = r_mat(seq(1, last), all) - r_mat(seq(0, last - 1), all);
      ArrayXXr<Real> dt_dr = d_theta(seq(1, last), all) * dr(all, seq(1, last));
      RowArrayXr<Real> r_row = r_mat.reshaped().transpose();
      RowArrayXr<Real> theta_row = theta_mat.reshaped().transpose();

      chrono::steady_clock::time_point begin, end;
      chrono::duration<double> elapsed_seconds;

      begin = chrono::steady_clock::now();
      unique_ptr<stTmatrix<Real>> T_mat = slvForT(params, options);
      end = chrono::steady_clock::now();
      elapsed_seconds = end - begin;
      out_stream << elapsed_seconds.count() << endl;
      out_stream << "T-matrix calculatred for excitation" << endl;
      MultiPrint(out_stream.str(), out_file_name, write_output);
      out_stream.str(string());

      begin = chrono::steady_clock::now();
      unique_ptr<stTmatrix<Real>> T_mat_rm = slvForT(params_rm, options);
      end = chrono::steady_clock::now();
      elapsed_seconds = end - begin;
      out_stream << elapsed_seconds.count() << endl;
      out_stream << "T-matrix calculatred for Raman" << endl;
      MultiPrint(out_stream.str(), out_file_name, write_output);
      out_stream.str(string());

      sca(k) = T_mat->st_coa->sca(0); // st_coa->sca should only contain 1 element
      ext(k) = T_mat->st_coa->ext(0); // st_coa->ext should only contain 1 element
      abs(k) = T_mat->st_coa->abs(0); // st_coa->abs should only contain 1 element
      out_stream << "Csca " << sca(k) << endl;
      out_stream << "Cext " << ext(k) << endl;
      out_stream << "Cabs " << abs(k) << endl;
      MultiPrint(out_stream.str(), out_file_name, write_output);
      out_stream.str(string());

      stSM_list[k] = pstScatteringMatrixOA(T_mat->st_TR_list, params->lambda(0), sca(k));

      for (int t = 0; t < theta_p_var.size(); t++) {
        begin = chrono::steady_clock::now();
        Real theta_p = theta_p_var(t);
        std::array<long int, 3> new_dims = {N_r, N_theta, 3};
        ArrayXr<Real> phi_var = ArrayXr<Real>::LinSpaced(N_phi + 1, 0, 2*PI);

        Real alpha_p = 0;
        params->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        unique_ptr<stAbcdnm<Real>> st_abcdnm = rvhGetFieldCoefficients(N, T_mat->st_TR_list, params->inc_par);
        unique_ptr<stRes<Real>> st_res_E = pstMakeStructForField(st_abcdnm, params);
        ArrayXXc<Real> c_nm = Map<RowArrayXc<Real>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        ArrayXXc<Real> d_nm = Map<RowArrayXc<Real>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        unique_ptr<stEAllPhi<Real>> st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda,
            st_res_E->epsilon2, c_nm, d_nm, r_row, theta_row, sBessel::J);
        Tensor4c<Real> E_field_z(N_r, N_theta, N_phi + 1, 3);
        E_field_z.setZero();
        for (int m = 0; m <= N_phi; m++) {
          Real phi = phi_var(m);
          unique_ptr<stEforPhi<Real>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->Er.cols(), 3*st_E_for_phi->Er.rows()};
          ArrayXXc<Real> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->Er.matrix().adjoint(),
              st_E_for_phi->Et.matrix().adjoint(), st_E_for_phi->Ef.matrix().adjoint();
          E_field_z.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        alpha_p = PI/2;
        params->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        st_abcdnm = rvhGetFieldCoefficients(N, T_mat->st_TR_list, params->inc_par);
        st_res_E = pstMakeStructForField(st_abcdnm, params);
        c_nm = Map<RowArrayXc<Real>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        d_nm = Map<RowArrayXc<Real>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda, st_res_E->epsilon2,
            c_nm, d_nm, r_row, theta_row, sBessel::J);
        Tensor4c<Real> E_field_y(N_r, N_theta, N_phi + 1, 3);
        E_field_y.setZero();
        for (int m = 0; m <= N_phi; m++) {
          Real phi = phi_var(m);
          unique_ptr<stEforPhi<Real>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->Er.cols(), 3*st_E_for_phi->Er.rows()};
          ArrayXXc<Real> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->Er.matrix().adjoint(),
              st_E_for_phi->Et.matrix().adjoint(), st_E_for_phi->Ef.matrix().adjoint();
          E_field_y.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        alpha_p = 0;
        params_rm->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        st_abcdnm = rvhGetFieldCoefficients(N, T_mat_rm->st_TR_list, params_rm->inc_par);
        st_res_E = pstMakeStructForField(st_abcdnm, params_rm);
        c_nm = Map<RowArrayXc<Real>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        d_nm = Map<RowArrayXc<Real>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda, st_res_E->epsilon2,
            c_nm, d_nm, r_row, theta_row, sBessel::J);
        Tensor4c<Real> E_field_rm_z(N_r, N_theta, N_phi + 1, 3);
        E_field_rm_z.setZero();
        for (int m = 0; m <= N_phi; m++) {
          Real phi = phi_var(m);
          unique_ptr<stEforPhi<Real>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->Er.cols(), 3*st_E_for_phi->Er.rows()};
          ArrayXXc<Real> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->Er.matrix().adjoint(),
              st_E_for_phi->Et.matrix().adjoint(), st_E_for_phi->Ef.matrix().adjoint();
          E_field_rm_z.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        alpha_p = PI/2;
        params_rm->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        st_abcdnm = rvhGetFieldCoefficients(N, T_mat_rm->st_TR_list, params_rm->inc_par);
        st_res_E = pstMakeStructForField(st_abcdnm, params_rm);
        c_nm = Map<RowArrayXc<Real>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        d_nm = Map<RowArrayXc<Real>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda, st_res_E->epsilon2,
            c_nm, d_nm, r_row, theta_row, sBessel::J);
        Tensor4c<Real> E_field_rm_y(N_r, N_theta, N_phi + 1, 3);
        E_field_rm_y.setZero();
        for (int m = 0; m <= N_phi; m++) {
          Real phi = phi_var(m);
          unique_ptr<stEforPhi<Real>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->Er.cols(), 3*st_E_for_phi->Er.rows()};
          ArrayXXc<Real> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->Er.matrix().adjoint(),
              st_E_for_phi->Et.matrix().adjoint(), st_E_for_phi->Ef.matrix().adjoint();
          E_field_rm_y.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        std::array<int, 1> dim3 = {3}, dim2 = {2};
        Tensor<Real, 2> inter_step = 2*PI/N_phi*(E_field_rm_z * tensor_conj(E_field_z))
            .sum(dim3).abs().pow(static_cast<Real>(2.0)).sum(dim2);
        ArrayXXr<Real> M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        ArrayXXr<Real> F = M_mat * r_mat.pow(2) * sin(theta_mat);
        ArrayXXr<Real> tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_zz(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real>(4.0)/3*PI*a*a*c);

        inter_step = 2*PI/N_phi*(E_field_rm_y * tensor_conj(E_field_z))
            .sum(dim3).abs().pow(static_cast<Real>(2.0)).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_yz(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real>(4.0)/3*PI*a*a*c);

        inter_step = 2*PI/N_phi*(E_field_rm_y * tensor_conj(E_field_y))
            .sum(dim3).abs().pow(static_cast<Real>(2.0)).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_yy(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real>(4.0)/3*PI*a*a*c);

        inter_step = 2*PI/N_phi*(E_field_rm_z * tensor_conj(E_field_y))
            .sum(dim3).abs().pow(static_cast<Real>(2.0)).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_zy(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real>(4.0)/3*PI*a*a*c);

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

// Can make this into a non-template function using c++20 'concepts' keyword
template <class Real>
vector<unique_ptr<stTR<double>>> stTRListMp2Double(const vector<unique_ptr<stTR<Real>>>& st_TR_list) {
  vector<unique_ptr<stTR<double>>> output(st_TR_list.size());
  for (size_t i = 0; i < st_TR_list.size(); i++) {
    assert(st_TR_list[i]);
    assert(st_TR_list[i]->mat_list.size() == 2);
    assert(st_TR_list[i]->mat_list[0] == "st_4M_T");
    assert(st_TR_list[i]->mat_list[1] == "st_4M_R");
    output[i] = make_unique<stTR<double>>();
    output[i]->mat_list = st_TR_list[i]->mat_list;

    output[i]->st_4M_T_eo().M11 = st_TR_list[i]->st_4M_T_eo().M11.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_T_eo().M11.imag().template cast<double>();
    output[i]->st_4M_T_eo().M12 = st_TR_list[i]->st_4M_T_eo().M12.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_T_eo().M12.imag().template cast<double>();
    output[i]->st_4M_T_eo().M21 = st_TR_list[i]->st_4M_T_eo().M21.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_T_eo().M21.imag().template cast<double>();
    output[i]->st_4M_T_eo().M22 = st_TR_list[i]->st_4M_T_eo().M22.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_T_eo().M22.imag().template cast<double>();
    output[i]->st_4M_T_eo().m = st_TR_list[i]->st_4M_T_eo().m;
    output[i]->st_4M_T_eo().ind1 = st_TR_list[i]->st_4M_T_eo().ind1;
    output[i]->st_4M_T_eo().ind2 = st_TR_list[i]->st_4M_T_eo().ind2;

    output[i]->st_4M_T_oe().M11 = st_TR_list[i]->st_4M_T_oe().M11.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_T_oe().M11.imag().template cast<double>();
    output[i]->st_4M_T_oe().M12 = st_TR_list[i]->st_4M_T_oe().M12.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_T_oe().M12.imag().template cast<double>();
    output[i]->st_4M_T_oe().M21 = st_TR_list[i]->st_4M_T_oe().M21.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_T_oe().M21.imag().template cast<double>();
    output[i]->st_4M_T_oe().M22 = st_TR_list[i]->st_4M_T_oe().M22.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_T_oe().M22.imag().template cast<double>();
    output[i]->st_4M_T_oe().m = st_TR_list[i]->st_4M_T_oe().m;
    output[i]->st_4M_T_oe().ind1 = st_TR_list[i]->st_4M_T_oe().ind1;
    output[i]->st_4M_T_oe().ind2 = st_TR_list[i]->st_4M_T_oe().ind2;

    output[i]->st_4M_R_eo().M11 = st_TR_list[i]->st_4M_R_eo().M11.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_R_eo().M11.imag().template cast<double>();
    output[i]->st_4M_R_eo().M12 = st_TR_list[i]->st_4M_R_eo().M12.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_R_eo().M12.imag().template cast<double>();
    output[i]->st_4M_R_eo().M21 = st_TR_list[i]->st_4M_R_eo().M21.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_R_eo().M21.imag().template cast<double>();
    output[i]->st_4M_R_eo().M22 = st_TR_list[i]->st_4M_R_eo().M22.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_R_eo().M22.imag().template cast<double>();
    output[i]->st_4M_R_eo().m = st_TR_list[i]->st_4M_R_eo().m;
    output[i]->st_4M_R_eo().ind1 = st_TR_list[i]->st_4M_R_eo().ind1;
    output[i]->st_4M_R_eo().ind2 = st_TR_list[i]->st_4M_R_eo().ind2;

    output[i]->st_4M_R_oe().M11 = st_TR_list[i]->st_4M_R_oe().M11.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_R_oe().M11.imag().template cast<double>();
    output[i]->st_4M_R_oe().M12 = st_TR_list[i]->st_4M_R_oe().M12.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_R_oe().M12.imag().template cast<double>();
    output[i]->st_4M_R_oe().M21 = st_TR_list[i]->st_4M_R_oe().M21.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_R_oe().M21.imag().template cast<double>();
    output[i]->st_4M_R_oe().M22 = st_TR_list[i]->st_4M_R_oe().M22.real().template cast<double>() +
        mp_im_unit<double>() * st_TR_list[i]->st_4M_R_oe().M22.imag().template cast<double>();
    output[i]->st_4M_R_oe().m = st_TR_list[i]->st_4M_R_oe().m;
    output[i]->st_4M_R_oe().ind1 = st_TR_list[i]->st_4M_R_oe().ind1;
    output[i]->st_4M_R_oe().ind2 = st_TR_list[i]->st_4M_R_oe().ind2;
  }
  return output;
}

template <class Real>
void RamanElasticScatteringMpDouble(string in_file_name, string out_file_name) {
  double PI = mp_pi<double>();
  unique_ptr<RamanParams<Real>> raman_params = GetRamanParams<Real>(in_file_name);
  ArrayXr<Real> rad_var = raman_params->rad_var;
  ArrayXd theta_p_var = raman_params->theta_p_var.template cast<double>();
  bool write_output = raman_params->write_output;
  ArrayXr<Real> h_var = {{static_cast<Real>(1.0)/3}};

  ofstream out_file;
  if (write_output) {
    out_file.open(out_file_name, ios::out | ios::app);
    out_file << endl;
    if (!out_file.is_open()) {
      cerr << "Warning: cannot create or open output file. Some output won't be written.";
      out_file.close();
      write_output = false;
    }
  }
  if (write_output)
    CreateTimeStamp(out_file_name, raman_params, true);

  Tensor3d sigma_yz(rad_var.size(), h_var.size(), theta_p_var.size());
  Tensor3d sigma_zy(rad_var.size(), h_var.size(), theta_p_var.size());
  Tensor3d sigma_zz(rad_var.size(), h_var.size(), theta_p_var.size());
  Tensor3d sigma_yy(rad_var.size(), h_var.size(), theta_p_var.size());
  ArrayXd sca = ArrayXd::Zero(rad_var.size());
  ArrayXd ext = ArrayXd::Zero(rad_var.size());
  ArrayXd abs = ArrayXd::Zero(rad_var.size());
  vector<unique_ptr<stSM<double>>> stSM_list(rad_var.size());

  auto options = make_unique<stOptions>();
  options->get_R = true;
  options->delta = 0;
  options->NB = 0;
  options->get_symmetric_T = false;
  unique_ptr<stParams<Real>> params_mp = loadParam<Real>();
  unique_ptr<stParams<Real>> params_mp_rm = loadParam<Real>("rm");
  unique_ptr<stParams<double>> params = loadParam<double>();
  unique_ptr<stParams<double>> params_rm = loadParam<double>("rm");
  double phi_p = 0;

  int N_r = 100;
  int N_theta = 320;
  int N_phi = 320;

  ArrayXd theta_surf = ArrayXd::LinSpaced(N_theta, 0, PI);
  RowArrayXd r_surf_u = RowArrayXd::LinSpaced(
      N_r, 1.0/N_r, 1.0).pow(1.0/3);
  ArrayXXd r_mat_u = r_surf_u.replicate(N_theta, 1).transpose();
  ArrayXXd theta_mat = theta_surf.replicate(1, N_r).transpose();

  int Nb_theta = 1000;
  int Nb_theta_pst = 1;
  params_mp->Nb_theta = Nb_theta;
  params_mp->Nb_theta_pst = Nb_theta_pst;
  params_mp_rm->Nb_theta = Nb_theta;
  params_mp_rm->Nb_theta_pst = Nb_theta_pst;
  params->Nb_theta = Nb_theta;
  params->Nb_theta_pst = Nb_theta_pst;
  params_rm->Nb_theta = Nb_theta;
  params_rm->Nb_theta_pst = Nb_theta_pst;

  for (int h_ind = 0; h_ind < h_var.size(); h_ind++) {
    stringstream out_stream;
    Real h = h_var(h_ind);

    out_stream.precision(4);
    out_stream << "aspect ratio " << h_var << endl;
    MultiPrint(out_stream.str(), out_file_name, write_output);
    out_stream.str(string());

    for (int k = 0; k < rad_var.size(); k++) {
      Real a_mp, c_mp;
      if (h > 1) {
        a_mp = rad_var(k);
        c_mp = rad_var(k) / h;
      } else {
        a_mp = rad_var(k) * h;
        c_mp = rad_var(k);
      }
      int N = 6 + 2*static_cast<int>(ceil(max(a_mp, c_mp)/40));
      params_mp->N = N;
      params_mp->a = a_mp;
      params_mp->c = c_mp;
      params_mp_rm->N = N;
      params_mp_rm->a = a_mp;
      params_mp_rm->c = c_mp;

      double a = static_cast<double>(a_mp);
      double c = static_cast<double>(c_mp);
      params->N = N;
      params->a = a;
      params->c = c;
      params_rm->N = N;
      params_rm->a = a;
      params_rm->c = c;

      ArrayXXd r_mat = a*c*r_mat_u/sqrt(pow(c, 2)*sin(theta_mat).pow(2) + pow(a, 2)*cos(theta_mat).pow(2));
      ArrayXXd d_theta = theta_mat(all, seq(1, last)) - theta_mat(all, seq(0, last - 1));
      ArrayXXd dr = r_mat(seq(1, last), all) - r_mat(seq(0, last - 1), all);
      ArrayXXd dt_dr = d_theta(seq(1, last), all) * dr(all, seq(1, last));
      RowArrayXd r_row = r_mat.reshaped().transpose();
      RowArrayXd theta_row = theta_mat.reshaped().transpose();

      chrono::steady_clock::time_point begin, end;
      chrono::duration<double> elapsed_seconds;

      begin = chrono::steady_clock::now();
      unique_ptr<stTmatrix<Real>> T_mat = slvForT(params_mp, options);
      end = chrono::steady_clock::now();
      elapsed_seconds = end - begin;
      out_stream << elapsed_seconds.count() << endl;
      out_stream << "T-matrix calculatred for excitation" << endl;
      MultiPrint(out_stream.str(), out_file_name, write_output);
      out_stream.str(string());

      begin = chrono::steady_clock::now();
      unique_ptr<stTmatrix<Real>> T_mat_rm = slvForT(params_mp_rm, options);
      end = chrono::steady_clock::now();
      elapsed_seconds = end - begin;
      out_stream << elapsed_seconds.count() << endl;
      out_stream << "T-matrix calculatred for Raman" << endl;
      MultiPrint(out_stream.str(), out_file_name, write_output);
      out_stream.str(string());

      vector<unique_ptr<stTR<double>>> st_TR_list = stTRListMp2Double(T_mat->st_TR_list);
      vector<unique_ptr<stTR<double>>> st_TR_list_rm = stTRListMp2Double(T_mat_rm->st_TR_list);

      sca(k) = static_cast<double>(T_mat->st_coa->sca(0)); // st_coa->sca should only contain 1 element
      ext(k) = static_cast<double>(T_mat->st_coa->ext(0)); // st_coa->ext should only contain 1 element
      abs(k) = static_cast<double>(T_mat->st_coa->abs(0)); // st_coa->abs should only contain 1 element
      out_stream << "Csca " << sca(k) << endl;
      out_stream << "Cext " << ext(k) << endl;
      out_stream << "Cabs " << abs(k) << endl;
      MultiPrint(out_stream.str(), out_file_name, write_output);
      out_stream.str(string());

      stSM_list[k] = pstScatteringMatrixOA(st_TR_list, params->lambda(0), sca(k));

      for (int t = 0; t < theta_p_var.size(); t++) {
        begin = chrono::steady_clock::now();
        double theta_p = theta_p_var(t);
        std::array<long int, 3> new_dims = {N_r, N_theta, 3};
        ArrayXd phi_var = ArrayXd::LinSpaced(N_phi + 1, 0, 2*PI);

        double alpha_p = 0;
        params->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        unique_ptr<stAbcdnm<double>> st_abcdnm = rvhGetFieldCoefficients(N, st_TR_list, params->inc_par);
        unique_ptr<stRes<double>> st_res_E = pstMakeStructForField(st_abcdnm, params);
        ArrayXXc<double> c_nm = Map<RowArrayXc<double>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        ArrayXXc<double> d_nm = Map<RowArrayXc<double>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        unique_ptr<stEAllPhi<double>> st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda,
            st_res_E->epsilon2, c_nm, d_nm, r_row, theta_row, sBessel::J);
        Tensor4c<double> E_field_z(N_r, N_theta, N_phi + 1, 3);
        E_field_z.setZero();
        for (int m = 0; m <= N_phi; m++) {
          double phi = phi_var(m);
          unique_ptr<stEforPhi<double>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->Er.cols(), 3*st_E_for_phi->Er.rows()};
          ArrayXXc<double> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->Er.matrix().adjoint(),
              st_E_for_phi->Et.matrix().adjoint(), st_E_for_phi->Ef.matrix().adjoint();
          E_field_z.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        alpha_p = PI/2;
        params->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        st_abcdnm = rvhGetFieldCoefficients(N, st_TR_list, params->inc_par);
        st_res_E = pstMakeStructForField(st_abcdnm, params);
        c_nm = Map<RowArrayXc<double>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        d_nm = Map<RowArrayXc<double>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda, st_res_E->epsilon2,
            c_nm, d_nm, r_row, theta_row, sBessel::J);
        Tensor4c<double> E_field_y(N_r, N_theta, N_phi + 1, 3);
        E_field_y.setZero();
        for (int m = 0; m <= N_phi; m++) {
          double phi = phi_var(m);
          unique_ptr<stEforPhi<double>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->Er.cols(), 3*st_E_for_phi->Er.rows()};
          ArrayXXc<double> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->Er.matrix().adjoint(),
              st_E_for_phi->Et.matrix().adjoint(), st_E_for_phi->Ef.matrix().adjoint();
          E_field_y.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        alpha_p = 0;
        params_rm->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        st_abcdnm = rvhGetFieldCoefficients(N, st_TR_list_rm, params_rm->inc_par);
        st_res_E = pstMakeStructForField(st_abcdnm, params_rm);
        c_nm = Map<RowArrayXc<double>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        d_nm = Map<RowArrayXc<double>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda, st_res_E->epsilon2,
            c_nm, d_nm, r_row, theta_row, sBessel::J);
        Tensor4c<double> E_field_rm_z(N_r, N_theta, N_phi + 1, 3);
        E_field_rm_z.setZero();
        for (int m = 0; m <= N_phi; m++) {
          double phi = phi_var(m);
          unique_ptr<stEforPhi<double>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->Er.cols(), 3*st_E_for_phi->Er.rows()};
          ArrayXXc<double> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->Er.matrix().adjoint(),
              st_E_for_phi->Et.matrix().adjoint(), st_E_for_phi->Ef.matrix().adjoint();
          E_field_rm_z.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        alpha_p = PI/2;
        params_rm->inc_par = vshMakeIncidentParams(sIncType::GENERAL, N, theta_p, phi_p, alpha_p);
        st_abcdnm = rvhGetFieldCoefficients(N, st_TR_list_rm, params_rm->inc_par);
        st_res_E = pstMakeStructForField(st_abcdnm, params_rm);
        c_nm = Map<RowArrayXc<double>>(st_res_E->c_nm.transpose().data(), st_res_E->c_nm.size());
        d_nm = Map<RowArrayXc<double>>(st_res_E->d_nm.transpose().data(), st_res_E->d_nm.size());
        st_E_surf = vshEgenThetaAllPhi(st_res_E->lambda, st_res_E->epsilon2,
            c_nm, d_nm, r_row, theta_row, sBessel::J);
        Tensor4c<double> E_field_rm_y(N_r, N_theta, N_phi + 1, 3);
        E_field_rm_y.setZero();
        for (int m = 0; m <= N_phi; m++) {
          double phi = phi_var(m);
          unique_ptr<stEforPhi<double>> st_E_for_phi = vshEthetaForPhi(st_E_surf, phi);
          long int dims[2] = {st_E_for_phi->Er.cols(), 3*st_E_for_phi->Er.rows()};
          ArrayXXc<double> E_field_phi(dims[0], dims[1]);
          E_field_phi << st_E_for_phi->Er.matrix().adjoint(),
              st_E_for_phi->Et.matrix().adjoint(), st_E_for_phi->Ef.matrix().adjoint();
          E_field_rm_y.chip(m, 2) = TensorCast(E_field_phi).reshape(new_dims);
        }

        std::array<int, 1> dim3 = {3}, dim2 = {2};
        Tensor<double, 2> inter_step = 2*PI/N_phi*(E_field_rm_z * tensor_conj(E_field_z))
            .sum(dim3).abs().pow(2.0).sum(dim2);
        ArrayXXd M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        ArrayXXd F = M_mat * r_mat.pow(2) * sin(theta_mat);
        ArrayXXd tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_zz(k, h_ind, t) = (tmp * dt_dr).sum()/(4.0/3*PI*a*a*c);

        inter_step = 2*PI/N_phi*(E_field_rm_y * tensor_conj(E_field_z))
            .sum(dim3).abs().pow(2.0).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_yz(k, h_ind, t) = (tmp * dt_dr).sum()/(4.0/3*PI*a*a*c);

        inter_step = 2*PI/N_phi*(E_field_rm_y * tensor_conj(E_field_y))
            .sum(dim3).abs().pow(2.0).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_yy(k, h_ind, t) = (tmp * dt_dr).sum()/(4.0/3*PI*a*a*c);

        inter_step = 2*PI/N_phi*(E_field_rm_z * tensor_conj(E_field_y))
            .sum(dim3).abs().pow(2.0).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_zy(k, h_ind, t) = (tmp * dt_dr).sum()/(4.0/3*PI*a*a*c);

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
