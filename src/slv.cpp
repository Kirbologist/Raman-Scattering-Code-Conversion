#include "slv.h"

using namespace Eigen;
using namespace std;

namespace Smarties {
  template <class Real>
  unique_ptr<stTmatrix<Real>> slvForT(const unique_ptr<stParams<Real>>& params,
      const unique_ptr<stOptions<Real>>& options, unique_ptr<stRtfunc<Real>> stGeometry) {
    Real c = params->c;
    Real a = params->a;

    int N = params->N;
    int Nb_theta = params->Nb_theta;
    ArrayXi abs_m_vec = ArrayXi::LinSpaced(N + 1, 0, N);

    if (!stGeometry)
      stGeometry = sphMakeGeometry(Nb_theta, a, c);

    // We won't be entering this if-block for Raman scattering
    // if (options->delta < 0) {
    //  sphEstimateDelta(stGeometry, params);
    // }

    int NQ = N + options->delta;
    int NB = options->NB;

    if (NB <= 0)
      NB = sphEstimateNB(NQ, stGeometry, params);
    if (NB < NQ)
      NB = NQ;

    vector<unique_ptr<stPQ<Real>>> st_PQ_list = sphCalculatePQ(NQ, abs_m_vec, stGeometry, params, NB);
    vector<unique_ptr<stTR<Real>>> st_TR_list = rvhGetTRfromPQ(st_PQ_list, options->get_R);

    if (NQ > N)
      st_TR_list = rvhTruncateMatrices(st_TR_list, N);
    if (options->get_symmetric_T)
      st_TR_list = rvhGetSymmetricMat(st_TR_list);

    vector<vector<unique_ptr<stTR<Real>>>> st_TR_array;
    st_TR_array.push_back(move(st_TR_list));
    unique_ptr<stCrossSection<Real>> st_coa = rvhGetAverageCrossSections(params->k1, st_TR_array);
    st_TR_list = move(st_TR_array[0]);

    unique_ptr<stTmatrix<Real>> output = make_unique<stTmatrix<Real>>();
    output->st_coa = move(st_coa);
    output->st_TR_list = move(st_TR_list);

    return output;
  }

  template unique_ptr<stTmatrix<double>> slvForT(const unique_ptr<stParams<double>>&,
      const unique_ptr<stOptions<double>>&, unique_ptr<stRtfunc<double>>);
}
