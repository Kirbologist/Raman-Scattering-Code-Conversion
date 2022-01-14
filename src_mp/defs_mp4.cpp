#include "../src/raman_elastic_scattering.hpp"
#include "../src/math.hpp"
#include "../src/vsh.hpp"
#include "../src/rvh.hpp"
#include "../src/slv.hpp"
#include "../src/pst.hpp"
#include "../src/core.hpp"
#include "core_mp.hpp"

using namespace Eigen;
using namespace Smarties;
using namespace boost::multiprecision;

extern template raman_float mp_pi();

extern template Tensor4c<raman_float> tensor_conj(Tensor4c<raman_float>&);

extern template unique_ptr<stIncPar<raman_float>> vshMakeIncidentParams(sIncType, int, raman_float, raman_float, raman_float);
extern template unique_ptr<stEAllPhi<raman_float>> vshEgenThetaAllPhi(const ArrayXr<raman_float>&,
    const ArrayXr<raman_float>&, const ArrayXXc<raman_float>&, const ArrayXXc<raman_float>&,
    const RowArrayXr<raman_float>&, const RowArrayXr<raman_float>&, sBessel, unique_ptr<stPinmTaunm<raman_float>>);
extern template unique_ptr<stEforPhi<raman_float>> vshEthetaForPhi(const unique_ptr<stEAllPhi<raman_float>>&, raman_float);

extern template unique_ptr<stAbcdnm<raman_float>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<raman_float>>>&,
    const unique_ptr<stIncPar<raman_float>>&, unique_ptr<stIncEabnm<raman_float>>);

extern template unique_ptr<stTmatrix<raman_float>> slvForT(const unique_ptr<stParams<raman_float>>&,
    const unique_ptr<stOptions>&, unique_ptr<stRtfunc<raman_float>>);

extern template unique_ptr<stRes<raman_float>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<raman_float>>&, const unique_ptr<stParams<raman_float>>&);
extern template unique_ptr<stSM<raman_float>> pstScatteringMatrixOA(
    const vector<unique_ptr<stTR<raman_float>>>&, raman_float, raman_float, int);


extern template double mp_pi<double>();

extern template Tensor4c<double> tensor_conj(Tensor4c<double>&);

extern template unique_ptr<stIncPar<double>> vshMakeIncidentParams(sIncType, int, double, double, double);
extern template unique_ptr<stEAllPhi<double>> vshEgenThetaAllPhi(const ArrayXd&,
    const ArrayXd&, const ArrayXXc<double>&, const ArrayXXc<double>&,
    const RowArrayXd&, const RowArrayXd&, sBessel, unique_ptr<stPinmTaunm<double>>);
extern template unique_ptr<stEforPhi<double>> vshEthetaForPhi(const unique_ptr<stEAllPhi<double>>&, double);

extern template unique_ptr<stAbcdnm<double>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<double>>>&,
    const unique_ptr<stIncPar<double>>&, unique_ptr<stIncEabnm<double>>);

extern template unique_ptr<stTmatrix<double>> slvForT(const unique_ptr<stParams<double>>&,
    const unique_ptr<stOptions>&, unique_ptr<stRtfunc<double>>);

extern template unique_ptr<stRes<double>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<double>>&, const unique_ptr<stParams<double>>&);
extern template unique_ptr<stSM<double>> pstScatteringMatrixOA(
    const vector<unique_ptr<stTR<double>>>&, double, double, int);


template unique_ptr<stParams<raman_float>> loadParam(string);
template unique_ptr<RamanParams<raman_float>> GetRamanParams(string in_file_name);
template void RamanElasticScattering<raman_float>(string, string);

extern template double mp_im_unit()
template vector<unique_ptr<stTR<double>>> stTRListMp2Double(const vector<unique_ptr<stTR<raman_float>>>&);
template void RamanElasticScatteringMpDouble<raman_float>(string, string);
