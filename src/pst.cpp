#include "pst.h"
#include "vsh.h"

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  unique_ptr<stRes<Real>> pstMakeStructForField(unique_ptr<stAbcdnm<Real>> st_abcdnm,
    int N_max, ArrayXr<Real> lambda, ArrayXr<Real> epsilon2, Real epsilon1,
    unique_ptr<stIncPar<Real>> inc_par, Real a, Real c) {
    if (!st_abcdnm->c_nm.size())
      cout << "Pb in  pstMakeStructForField: the structure stAbcdnm should contain c_nm and d_nm for internal fields" << endl;
    unique_ptr<stRes<Real>> output = make_unique<stRes<Real>>();
    output->N_max = N_max;
    output->lambda = lambda;
    output->epsilon1 = epsilon1;
    output->epsilon2 = epsilon2;
    output->inc_par = move(inc_par);
    if (!isnan(a))
      output->a = a;
    if (!isnan(c))
      output->a = c;
    return output;
  }

  template <class Real>
  unique_ptr<stRes<Real>> pstMakeStructForField(
      unique_ptr<stAbcdnm<Real>> st_abcdnm, unique_ptr<stParams<Real>> params) {
    if (!st_abcdnm->c_nm.size())
      cout << "Pb in  pstMakeStructForField: the structure stAbcdnm should contain c_nm and d_nm for internal fields" << endl;
    unique_ptr<stRes<Real>> output = make_unique<stRes<Real>>();
    output->N_max = params->N;
    output->lambda = params->lambda;
    output->epsilon1 = params->epsilon1;
    output->epsilon2 = params->epsilon2;
    output->a = params->a;
    output->c = params->c;
    if (params->inc_par)
      output->inc_par = move(params->inc_par);
    else {
      sIncType inc_type = params->inc_type;
      output->inc_par = vshMakeIncidentParams<Real>(inc_type, params->N);
    }
    return output;
  }

  template unique_ptr<stRes<double>> pstMakeStructForField(unique_ptr<stAbcdnm<double>>,
      int, ArrayXr<double>, ArrayXr<double>, double, unique_ptr<stIncPar<double>>, double, double);
  template unique_ptr<stRes<double>> pstMakeStructForField(
      unique_ptr<stAbcdnm<double>>, unique_ptr<stParams<double>>);
}
