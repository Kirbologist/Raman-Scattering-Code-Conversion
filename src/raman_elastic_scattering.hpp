#ifndef RAMAN_ELASTIC_SCATTERING_HPP
#define RAMAN_ELASTIC_SCATTERING_HPP

#include "core.hpp"
#include "smarties.hpp"
#include <chrono>

using namespace std;
using namespace Eigen;
using namespace Smarties;

template <class Real>
unique_ptr<stParams<Real>> loadParam(string type = "") {
  unique_ptr<stParams<Real>> params = make_unique<stParams<Real>>();
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
void RamanElasticScattering(int argc, char** argv) {
  Real PI = mp_pi<Real>();
  int cpu_n = (argc > 1 && isdigit(argv[1][0])) ? stoi(argv[1], nullptr) : 0;
  int cpus = (argc > 2 && isdigit(argv[2][0])) ? stoi(argv[2], nullptr) : 1;
  Real dia_min = (argc > 3 && isdigit(argv[3][0])) ? stof(argv[3], nullptr) : 1000.0;
  Real dia_max = (argc > 4 && isdigit(argv[4][0])) ? stof(argv[4], nullptr) : 2000.0;
  int N_rad = (argc > 5 && isdigit(argv[5][0])) ? stoi(argv[5], nullptr) : 100;
  int N_theta_p = (argc > 6 && isdigit(argv[6][0])) ? stoi(argv[6], nullptr) : 19;
  ArrayXr<Real> par = ArrayXr<Real>::LinSpaced(N_rad, dia_min, dia_max);
  int par_per_cpu = N_rad / cpus;
  ArrayXr<Real> dia_var = par(ArrayXi::LinSpaced(par_per_cpu, 0, par_per_cpu - 1) + cpu_n*par_per_cpu);
  ArrayXr<Real> rad_var = dia_var/2;
  ArrayXr<Real> theta_p_var = ArrayXr<Real>::LinSpaced(N_theta_p, 0, PI/2);
  ArrayXr<Real> h_var = {{static_cast<Real>(1.0)/3}};

  Tensor3r<Real> sigma_yz(par_per_cpu, h_var.size(), N_theta_p);
  Tensor3r<Real> sigma_zy(par_per_cpu, h_var.size(), N_theta_p);
  Tensor3r<Real> sigma_zz(par_per_cpu, h_var.size(), N_theta_p);
  Tensor3r<Real> sigma_yy(par_per_cpu, h_var.size(), N_theta_p);
  ArrayXr<Real> sca = ArrayXr<Real>::Zero(par_per_cpu);
  ArrayXr<Real> ext = ArrayXr<Real>::Zero(par_per_cpu);
  ArrayXr<Real> abs = ArrayXr<Real>::Zero(par_per_cpu);
  vector<unique_ptr<stSM<Real>>> stSM_list(rad_var.size());

  unique_ptr<stOptions<Real>> options = make_unique<stOptions<Real>>();
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
    Real h = h_var(h_ind);
    cout.precision(4);
    cout << "aspect ratio " << h_var << endl;
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
      cout << elapsed_seconds.count() << endl;
      cout << "T-matrix calculatred for excitation" << endl;

      begin = chrono::steady_clock::now();
      unique_ptr<stTmatrix<Real>> T_mat_rm = slvForT(params_rm, options);
      end = chrono::steady_clock::now();
      elapsed_seconds = end - begin;
      cout << elapsed_seconds.count() << endl;
      cout << "T-matrix calculatred for Raman" << endl;

      sca(k) = T_mat->st_coa->sca(0); // st_coa->sca should only contain 1 element
      ext(k) = T_mat->st_coa->ext(0); // st_coa->ext should only contain 1 element
      abs(k) = T_mat->st_coa->abs(0); // st_coa->abs should only contain 1 element
      cout << "Csca " << sca(k) << endl;
      cout << "Cext " << ext(k) << endl;
      cout << "Cabs " << abs(k) << endl;

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
            .sum(dim3).abs().pow(static_cast<Real>(2)).sum(dim2);
        ArrayXXr<Real> M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        ArrayXXr<Real> F = M_mat * r_mat.pow(2) * sin(theta_mat);
        ArrayXXr<Real> tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_zz(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real>(4.0)/3*PI*a*a*c);

        inter_step = 2*PI/N_phi*(E_field_rm_y * tensor_conj(E_field_z))
            .sum(dim3).abs().pow(static_cast<Real>(2)).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_yz(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real>(4.0)/3*PI*a*a*c);

        inter_step = 2*PI/N_phi*(E_field_rm_y * tensor_conj(E_field_y))
            .sum(dim3).abs().pow(static_cast<Real>(2)).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_yy(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real>(4.0)/3*PI*a*a*c);

        inter_step = 2*PI/N_phi*(E_field_rm_z * tensor_conj(E_field_y))
            .sum(dim3).abs().pow(static_cast<Real>(2)).sum(dim2);
        M_mat = MatrixCast(inter_step, N_r, N_theta).array();
        F = M_mat * r_mat.pow(2) * sin(theta_mat);
        tmp = (F(seq(1, last), seq(1, last)) + F(seq(0, last - 1), seq(1, last)) +
            F(seq(1, last), seq(0, last - 1)) + F(seq(0, last - 1), seq(0, last - 1))) / 4;
        sigma_zy(k, h_ind, t) = (tmp * dt_dr).sum()/(static_cast<Real>(4.0)/3*PI*a*a*c);

        cout.precision(5);
        cout << "max diameter = " << max(a, c)*2e-3;
        cout.precision(10);
        cout << " µm, theta = " << theta_p * 180/PI << " degrees" << endl;
        cout.precision(11);
        cout << "--- relative raman zz = " << sigma_zz(k, h_ind, t) << endl;
        cout << "--- relative raman yz = " << sigma_yz(k, h_ind, t) << endl;
        cout << "--- relative raman zy = " << sigma_zy(k, h_ind, t) << endl;
        cout << "--- relative raman yy = " << sigma_yy(k, h_ind, t) << endl;

        end = chrono::steady_clock::now();
        elapsed_seconds = end - begin;
        cout << elapsed_seconds.count() << endl;
      }
    }
  }
}

#endif
