#include "raman_elastic_scattering.h"
#include "math.h"
#include "sph.h"
#include "rvh.h"

using namespace std;
using namespace Eigen;
using namespace Raman;

template <class Real>
unique_ptr<stParams<Real>> loadParam(string type = "") {
  unique_ptr<stParams<Real>> params = make_unique<stParams<Real>>();
  params->epsilon1 = 1;
  if (type == "rm") {
    params->lambda = 403.7;
    params->epsilon2 = pow(1.344, 2);
    params->k1 = 2*PI/params->lambda;
  } else {
    params->lambda = 355;
    params->epsilon2 = ArrayXc<Real>(1);
    params->epsilon2(0) = complex<Real>(pow(1.35, 2), 0.0);
    params->k1 = 2*PI/params->lambda;
  }
  params->s = ArrayXr<Real>(1);
  params->s(0) = 1.35;
  return params;
}

template <class Real>
void RamanElasticScattering(int argc, char** argv) {
  int cpu_n = (argc > 1 && isdigit(argv[1][0])) ? stoi(argv[1], nullptr) : 0;
  int cpus = (argc > 2 && isdigit(argv[2][0])) ? stoi(argv[2], nullptr) : 1;
  Real dia_min = (argc > 3 && isdigit(argv[3][0])) ? stof(argv[3], nullptr) : 1000.0;
  Real dia_max = (argc > 4 && isdigit(argv[4][0])) ? stof(argv[4], nullptr) : 2000.0;
  int N_rad = (argc > 5 && isdigit(argv[5][0])) ? stoi(argv[5], nullptr) : 100;
  int N_theta_p = (argc > 6 && isdigit(argv[6][0])) ? stoi(argv[6], nullptr) : 19;
  ArrayXr<Real> par = ArrayXd::LinSpaced(N_rad, dia_min, dia_max);
  int par_per_cpu = N_rad / cpus;
  ArrayXr<Real> dia_var = par(ArrayXi::LinSpaced(par_per_cpu, 0, par_per_cpu - 1) + cpu_n*par_per_cpu);
  ArrayXr<Real> rad_var = dia_var/2;
  ArrayXr<Real> theta_p_var = ArrayXd::LinSpaced(0, PI/2, N_theta_p);
  Real hvar = 1.0/3.0;

  ArrayXXr<Real> sigma_yz = ArrayXXd::Zero(par_per_cpu, N_theta_p);
  ArrayXXr<Real> sigma_zy = ArrayXXd::Zero(par_per_cpu, N_theta_p);
  ArrayXXr<Real> sigma_zz = ArrayXXd::Zero(par_per_cpu, N_theta_p);
  ArrayXXr<Real> sigma_yy = ArrayXXd::Zero(par_per_cpu, N_theta_p);
  ArrayXr<Real> sca = ArrayXd::Zero(par_per_cpu);
  ArrayXr<Real> ext = ArrayXd::Zero(par_per_cpu);
  ArrayXr<Real> abs = ArrayXd::Zero(par_per_cpu);
}

template unique_ptr<stParams<double>> loadParam(string);
template void RamanElasticScattering<double>(int, char**);

int main(int argc, char** argv) {
  RamanElasticScattering<double>(argc, argv);
  //ArrayXr<double> x = ArrayXr<double>::LinSpaced(5, 1, 5);
  //stBesselProducts<double>* test = sphGetModifiedBesselProducts<double>(10, 0.5, x, 12);
  //destructStBesselProducts<double>(test);
  unique_ptr<stRtfunc<double>> test = sphMakeGeometry<double>(5, 100, 150);
  unique_ptr<stParams<double>> params = make_unique<stParams<double>>();
  params->k1 = {{2*M_PI/355}};
  params->s = {{1.35}};
  //size_t test2 = sphEstimateNB<double>(8, test, params);
  vector<unique_ptr<stPQ<double>>> test2 = sphCalculatePQ(10, ArrayXi::LinSpaced(11, 0, 10), test, params, 12);
  cout << test2.size() << endl;
  cout << "Now trying to access test2[0]" << endl;
  cout << test2[0]->st_4M_P_oe().m << endl;
  cout << "Now trying to access test2[10]:" << endl;
  cout << test2[10]->st_4M_P_oe().m << endl;
  return 0;
}
