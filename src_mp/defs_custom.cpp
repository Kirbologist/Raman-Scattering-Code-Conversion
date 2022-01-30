#include "../src/smarties.hpp"
#include "../src/raman_elastic_scattering.hpp"
#include "core_mp.hpp"

using namespace Eigen;
using namespace Smarties;

template raman_float mp_pi<raman_float>();
template raman_float mp_eps<raman_float>();
template complex<raman_float> mp_im_unit<raman_float>();

template Tensor3c<raman_float> subtensor(Tensor3c<raman_float>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
template ArrayXXc<raman_float> reduceAndSlice(Tensor3c<raman_float>&, int, int);
template ArrayXXc<raman_float> invertLUcol(MatrixXc<raman_float>&);
template ArrayXi logicalSlice(ArrayXi&, ArrayXb&);
template RowArrayXr<raman_float> logicalSlice(RowArrayXr<raman_float>&, RowArrayXb&);
template ArrayXr<raman_float> logicalSlice(ArrayXr<raman_float>&, ArrayXb&);
template Tensor4c<raman_float> tensor_conj(Tensor4c<raman_float>&);
template ArrayXr<raman_float> arr_bessel_j(ArrayXr<raman_float>&, raman_float);
template ArrayXr<raman_float> arr_bessel_y(ArrayXr<raman_float>&, raman_float);

/*
extern template raman_float mp_pi<raman_float>();
extern template raman_float mp_eps<raman_float>();
*/

template unique_ptr<stGLQuad<raman_float>> auxInitLegendreQuad(int, raman_float, raman_float);
template unique_ptr<stRtfunc<raman_float>> auxPrepareIntegrals(int, sInt);

/*
extern template raman_float mp_pi<raman_float>();
extern template raman_float mp_im_unit<raman_float>();

extern template ArrayXr<raman_float> arr_bessel_j(ArrayXr<raman_float>&, raman_float);
extern template ArrayXr<raman_float> arr_bessel_y(ArrayXr<raman_float>&, raman_float);
*/

template unique_ptr<stIncPar<raman_float>> vshMakeIncidentParams(sIncType, int);
template unique_ptr<stIncPar<raman_float>> vshMakeIncidentParams(sIncType, int, raman_float, raman_float, raman_float);
template unique_ptr<stPinmTaunm<raman_float>> vshPinmTaunm(int, const ArrayXr<raman_float>&);
template unique_ptr<stIncEabnm<raman_float>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<raman_float>>&);
template unique_ptr<stZnAll<raman_float>> vshGetZnAll(int, const ArrayXr<raman_float>&, sBessel);
template unique_ptr<stEAllPhi<raman_float>> vshEgenThetaAllPhi(const ArrayXr<raman_float>&,
    const ArrayXr<raman_float>&, const ArrayXXc<raman_float>&, const ArrayXXc<raman_float>&,
    const RowArrayXr<raman_float>&, const RowArrayXr<raman_float>&, sBessel, unique_ptr<stPinmTaunm<raman_float>>);
template unique_ptr<stEforPhi<raman_float>> vshEthetaForPhi(const unique_ptr<stEAllPhi<raman_float>>&, raman_float);
template ArrayXXc<raman_float> vshRBchi(ArrayXr<raman_float>, const ArrayXr<raman_float>&);
template ArrayXXc<raman_float> vshRBpsi(ArrayXr<raman_float>, const ArrayXr<raman_float>&);

/*
extern template mp_eps<raman_float>();
extern template mp_im_unit<raman_float>();

extern template Tensor3c<raman_float> subtensor(Tensor3c<raman_float>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
extern template ArrayXXc<raman_float> reduceAndSlice(Tensor3c<raman_float>&, int, int);
extern template RowArrayXr<raman_float> logicalSlice(RowArrayXr<raman_float>&, RowArrayXb&);

extern template unique_ptr<stRtfunc<raman_float>> auxPrepareIntegrals(int, sInt);

extern template unique_ptr<stPinmTaunm<raman_float>> vshPinmTaunm(int, const ArrayXr<raman_float>&);
extern template ArrayXXc<raman_float> vshRBchi(ArrayXr<raman_float>, const ArrayXr<raman_float>&);
extern template ArrayXXc<raman_float> vshRBpsi(ArrayXr<raman_float>, const ArrayXr<raman_float>&);
*/

template unique_ptr<ArrayXXr<raman_float>> sphGetUforFp(int);
template unique_ptr<stFprow<raman_float>> sphGetFpRow(int, raman_float, const ArrayXr<raman_float>&);
template unique_ptr<stFpovx<raman_float>> sphGetFpovx(int, raman_float, const ArrayXr<raman_float>&);
template unique_ptr<stBessel<raman_float>> sphGetXiPsi(int, raman_float, const ArrayXr<raman_float>&, int);
template unique_ptr<stBesselPrimes<raman_float>> sphGetBesselProductsPrimes(const Tensor3c<raman_float>&);
template unique_ptr<stBesselProducts<raman_float>> sphGetModifiedBesselProducts(int, raman_float, const ArrayXr<raman_float>&, int);
template unique_ptr<stRtfunc<raman_float>> sphMakeGeometry(int, raman_float, raman_float, const unique_ptr<ArrayXr<raman_float>>);
template int sphCheckBesselConvergence(int, raman_float, const ArrayXr<raman_float>&, raman_float, int);
template int sphEstimateNB(int, const unique_ptr<stRtfunc<raman_float>>&, const unique_ptr<stParams<raman_float>>&, raman_float);
template vector<unique_ptr<stPQ<raman_float>>> sphCalculatePQ(int, const ArrayXi&,
    const unique_ptr<stRtfunc<raman_float>>&, const unique_ptr<stParams<raman_float>>&, int);

/*
extern template raman_float mp_pi<raman_float();

extern template ArrayXXc<raman_float> invertLUcol(MatrixXc<raman_float>&);
extern template ArrayXi logicalSlice(ArrayXi&, ArrayXb&);

extern template unique_ptr<stIncEabnm<raman_float>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<raman_float>>&);
*/

template vector<unique_ptr<stTR<raman_float>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<raman_float>>>&, bool);
template vector<unique_ptr<stTR<raman_float>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<raman_float>>>&, int);
template vector<unique_ptr<stTR<raman_float>>> rvhGetSymmetricMat(const vector<unique_ptr<stTR<raman_float>>>&, vector<string>);
template unique_ptr<stAbcdnm<raman_float>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<raman_float>>>&,
    const unique_ptr<stIncPar<raman_float>>&, unique_ptr<stIncEabnm<raman_float>>);
template unique_ptr<stCrossSection<raman_float>> rvhGetAverageCrossSections(
    const ArrayXr<raman_float>&, const vector<vector<unique_ptr<stTR<raman_float>>>>&);

/*
extern template int sphEstimateNB(int, const unique_ptr<stRtfunc<raman_float>>&,
    const unique_ptr<stParams<raman_float>>&, raman_float);
extern template vector<unique_ptr<stPQ<raman_float>>> sphCalculatePQ(int, const ArrayXi&,
    const unique_ptr<stRtfunc<raman_float>>&, const unique_ptr<stParams<raman_float>>&, int);

extern template vector<unique_ptr<stTR<raman_float>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<raman_float>>>&, bool);
extern template vector<unique_ptr<stTR<raman_float>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<raman_float>>>&, int);
extern template vector<unique_ptr<stTR<raman_float>>> rvhGetSymmetricMat(
    const vector<unique_ptr<stTR<raman_float>>>&, vector<string>);
extern template unique_ptr<stCrossSection<raman_float>> rvhGetAverageCrossSections(
    const ArrayXr<raman_float>&, const vector<vector<unique_ptr<stTR<raman_float>>>>&);
*/

template unique_ptr<stTmatrix<raman_float>> slvForT(const unique_ptr<stParams<raman_float>>&,
    const unique_ptr<stOptions>&, unique_ptr<stRtfunc<raman_float>>);

/*
extern raman_float mp_pi<raman_float>();
extern raman_float mp_eps<raman_float>();
extern raman_float mp_im_unit<raman_float>();

extern template unique_ptr<stIncPar<raman_float>> vshMakeIncidentParams(sIncType, int);
*/

template unique_ptr<stRes<raman_float>> pstMakeStructForField(const unique_ptr<stAbcdnm<raman_float>>&,
    int, ArrayXr<raman_float>, ArrayXr<raman_float>, raman_float, unique_ptr<stIncPar<raman_float>>, raman_float, raman_float);
template unique_ptr<stRes<raman_float>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<raman_float>>&, const unique_ptr<stParams<raman_float>>&);
template unique_ptr<stSM<raman_float>> pstScatteringMatrixOA(const vector<unique_ptr<stTR<raman_float>>>&,
    raman_float, raman_float, int);

/*
extern template raman_float mp_pi<raman_float>();

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
*/

extern template float mp_pi();
extern template complex<float> mp_im_unit();

extern template Tensor4c<float> tensor_conj(Tensor4c<float>&);

extern template unique_ptr<stIncPar<float>> vshMakeIncidentParams(sIncType, int, float, float, float);
extern template unique_ptr<stEAllPhi<float>> vshEgenThetaAllPhi(const ArrayXf&,
    const ArrayXf&, const ArrayXXc<float>&, const ArrayXXc<float>&,
    const RowArrayXf&, const RowArrayXf&, sBessel, unique_ptr<stPinmTaunm<float>>);
extern template unique_ptr<stEforPhi<float>> vshEthetaForPhi(const unique_ptr<stEAllPhi<float>>&, float);

extern template unique_ptr<stAbcdnm<float>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<float>>>&,
    const unique_ptr<stIncPar<float>>&, unique_ptr<stIncEabnm<float>>);

extern template unique_ptr<stTmatrix<float>> slvForT(const unique_ptr<stParams<float>>&,
    const unique_ptr<stOptions>&, unique_ptr<stRtfunc<float>>);

extern template unique_ptr<stRes<float>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<float>>&, const unique_ptr<stParams<float>>&);
extern template unique_ptr<stSM<float>> pstScatteringMatrixOA(
    const vector<unique_ptr<stTR<float>>>&, float, float, int);

extern template float Frac2Float(string);
extern template unique_ptr<RamanParams<float>> LoadParams<float>(string);
extern template unique_ptr<stParams<float>> Raman2SmartiesParams(const unique_ptr<RamanParams<float>>&, int, int, string);

extern template double mp_pi();
extern template complex<double> mp_im_unit();

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

extern template double Frac2Float(string);
extern template unique_ptr<RamanParams<double>> LoadParams<double>(string);
extern template unique_ptr<stParams<double>> Raman2SmartiesParams(const unique_ptr<RamanParams<double>>&, int, int, string);

extern template long double mp_pi<long double>();
extern template complex<long double> mp_im_unit();

extern template Tensor4c<long double> tensor_conj(Tensor4c<long double>&);

extern template unique_ptr<stIncPar<long double>> vshMakeIncidentParams(
    sIncType, int, long double, long double, long double);
extern template unique_ptr<stEAllPhi<long double>> vshEgenThetaAllPhi(const ArrayXr<long double>&,
    const ArrayXr<long double>&, const ArrayXXc<long double>&, const ArrayXXc<long double>&,
    const RowArrayXr<long double>&, const RowArrayXr<long double>&, sBessel, unique_ptr<stPinmTaunm<long double>>);
extern template unique_ptr<stEforPhi<long double>> vshEthetaForPhi(
    const unique_ptr<stEAllPhi<long double>>&, long double);

extern template unique_ptr<stAbcdnm<long double>> rvhGetFieldCoefficients(
    int, const vector<unique_ptr<stTR<long double>>>&,
    const unique_ptr<stIncPar<long double>>&, unique_ptr<stIncEabnm<long double>>);

extern template unique_ptr<stTmatrix<long double>> slvForT(const unique_ptr<stParams<long double>>&,
    const unique_ptr<stOptions>&, unique_ptr<stRtfunc<long double>>);

extern template unique_ptr<stRes<long double>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<long double>>&, const unique_ptr<stParams<long double>>&);
extern template unique_ptr<stSM<long double>> pstScatteringMatrixOA(
    const vector<unique_ptr<stTR<long double>>>&, long double, long double, int);

extern template long double Frac2Float(string);
extern template unique_ptr<RamanParams<long double>> LoadParams<long double>(string);
extern template unique_ptr<stParams<long double>> Raman2SmartiesParams(
    const unique_ptr<RamanParams<long double>>&, int, int, string);

template raman_float Frac2Float(string);
template unique_ptr<RamanParams<raman_float>> LoadParams<raman_float>(string);
template unique_ptr<stParams<raman_float>> Raman2SmartiesParams(
    const unique_ptr<RamanParams<raman_float>>&, int, int, string);
template void CreateTimeStamp<raman_float, float>(string, const unique_ptr<RamanParams<raman_float>>&);
template void CreateTimeStamp<raman_float, double>(string, const unique_ptr<RamanParams<raman_float>>&);
template void CreateTimeStamp<raman_float, long double>(string, const unique_ptr<RamanParams<raman_float>>&);
template void CreateTimeStamp<raman_float, raman_float>(string, const unique_ptr<RamanParams<raman_float>>&);
template vector<unique_ptr<stTR<float>>> ConvertstTRList(const vector<unique_ptr<stTR<raman_float>>>&);
template vector<unique_ptr<stTR<double>>> ConvertstTRList(const vector<unique_ptr<stTR<raman_float>>>&);
template vector<unique_ptr<stTR<long double>>> ConvertstTRList(const vector<unique_ptr<stTR<raman_float>>>&);
template vector<unique_ptr<stTR<raman_float>>> ConvertstTRList(const vector<unique_ptr<stTR<raman_float>>>&);
template void RamanElasticScattering<raman_float, float>(string, string, int);
template void RamanElasticScattering<raman_float, double>(string, string, int);
template void RamanElasticScattering<raman_float, long double>(string, string, int);
template void RamanElasticScattering<raman_float, raman_float>(string, string, int);

template void CreateTimeStamp<float, raman_float>(string, const unique_ptr<RamanParams<float>>&);
template void CreateTimeStamp<double, raman_float>(string, const unique_ptr<RamanParams<double>>&);
template void CreateTimeStamp<long double, raman_float>(string, const unique_ptr<RamanParams<long double>>&);
template vector<unique_ptr<stTR<raman_float>>> ConvertstTRList(const vector<unique_ptr<stTR<float>>>&);
template vector<unique_ptr<stTR<raman_float>>> ConvertstTRList(const vector<unique_ptr<stTR<double>>>&);
template vector<unique_ptr<stTR<raman_float>>> ConvertstTRList(const vector<unique_ptr<stTR<long double>>>&);
template void RamanElasticScattering<float, raman_float>(string, string, int);
template void RamanElasticScattering<double, raman_float>(string, string, int);
template void RamanElasticScattering<long double, raman_float>(string, string, int);
