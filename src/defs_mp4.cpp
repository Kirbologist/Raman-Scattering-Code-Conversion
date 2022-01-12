#include "raman_elastic_scattering.hpp"
#include "math.hpp"
#include "vsh.hpp"
#include "rvh.hpp"
#include "slv.hpp"
#include "pst.hpp"
#include "core.hpp"
#include "core_mp.hpp"

using namespace Eigen;
using namespace Smarties;
using namespace boost::multiprecision;

extern template complex<raman_float> mp_pi();

extern template Tensor4c<raman_float> tensor_conj(Tensor4c<raman_float>&);

extern template unique_ptr<stIncPar<raman_float>> vshMakeIncidentParams(sIncType, int, raman_float, raman_float, raman_float);
extern template unique_ptr<stEAllPhi<raman_float>> vshEgenThetaAllPhi(const ArrayXr<raman_float>&,
    const ArrayXr<raman_float>&, const ArrayXXc<raman_float>&, const ArrayXXc<raman_float>&,
    const RowArrayXr<raman_float>&, const RowArrayXr<raman_float>&, sBessel, unique_ptr<stPinmTaunm<raman_float>>);
extern template unique_ptr<stEforPhi<raman_float>> vshEthetaForPhi(const unique_ptr<stEAllPhi<raman_float>>&, raman_float);

extern template unique_ptr<stAbcdnm<raman_float>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<raman_float>>>&,
    const unique_ptr<stIncPar<raman_float>>&, unique_ptr<stIncEabnm<raman_float>>);

extern template unique_ptr<stTmatrix<raman_float>> slvForT(const unique_ptr<stParams<raman_float>>&,
    const unique_ptr<stOptions<raman_float>>&, unique_ptr<stRtfunc<raman_float>>);

extern template unique_ptr<stRes<raman_float>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<raman_float>>&, const unique_ptr<stParams<raman_float>>&);
extern template unique_ptr<stSM<raman_float>> pstScatteringMatrixOA(
    const vector<unique_ptr<stTR<raman_float>>>&, raman_float, raman_float, int);

template unique_ptr<stParams<raman_float>> loadParam(string);
template void RamanElasticScattering<raman_float>(int, char**);
