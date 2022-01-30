#include "smarties.hpp"
#include "raman_elastic_scattering.hpp"

using namespace Eigen;
using namespace Smarties;

template float mp_pi<float>();
template float mp_eps<float>();
template complex<float> mp_im_unit<float>();

template Tensor3c<float> subtensor(Tensor3c<float>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
template ArrayXXc<float> reduceAndSlice(Tensor3c<float>&, int, int);
template ArrayXXc<float> invertLUcol(MatrixXc<float>&);
template ArrayXi logicalSlice(ArrayXi&, ArrayXb&);
template RowArrayXr<float> logicalSlice(RowArrayXr<float>&, RowArrayXb&);
template ArrayXr<float> logicalSlice(ArrayXr<float>&, ArrayXb&);
template Tensor4c<float> tensor_conj(Tensor4c<float>&);
template ArrayXr<float> arr_bessel_j(ArrayXr<float>&, float);
template ArrayXr<float> arr_bessel_y(ArrayXr<float>&, float);

/*
extern template float mp_pi<float>();
extern template float mp_eps<float>();
*/

template unique_ptr<stGLQuad<float>> auxInitLegendreQuad(int, float, float);
template unique_ptr<stRtfunc<float>> auxPrepareIntegrals(int, sInt);

/*
extern template float mp_pi<float>();
extern template float mp_im_unit<float>();

extern template ArrayXr<float> arr_bessel_j(ArrayXr<float>&, float);
extern template ArrayXr<float> arr_bessel_y(ArrayXr<float>&, float);
*/

template unique_ptr<stIncPar<float>> vshMakeIncidentParams(sIncType, int);
template unique_ptr<stIncPar<float>> vshMakeIncidentParams(sIncType, int, float, float, float);
template unique_ptr<stPinmTaunm<float>> vshPinmTaunm(int, const ArrayXr<float>&);
template unique_ptr<stIncEabnm<float>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<float>>&);
template unique_ptr<stZnAll<float>> vshGetZnAll(int, const ArrayXr<float>&, sBessel);
template unique_ptr<stEAllPhi<float>> vshEgenThetaAllPhi(const ArrayXr<float>&,
    const ArrayXr<float>&, const ArrayXXc<float>&, const ArrayXXc<float>&,
    const RowArrayXr<float>&, const RowArrayXr<float>&, sBessel, unique_ptr<stPinmTaunm<float>>);
template unique_ptr<stEforPhi<float>> vshEthetaForPhi(const unique_ptr<stEAllPhi<float>>&, float);
template ArrayXXc<float> vshRBchi(ArrayXr<float>, const ArrayXr<float>&);
template ArrayXXc<float> vshRBpsi(ArrayXr<float>, const ArrayXr<float>&);

/*
extern template mp_eps<float>();
extern template mp_im_unit<float>();

extern template Tensor3c<float> subtensor(Tensor3c<float>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
extern template ArrayXXc<float> reduceAndSlice(Tensor3c<float>&, int, int);
extern template RowArrayXr<float> logicalSlice(RowArrayXr<float>&, RowArrayXb&);

extern template unique_ptr<stRtfunc<float>> auxPrepareIntegrals(int, sInt);

extern template unique_ptr<stPinmTaunm<float>> vshPinmTaunm(int, const ArrayXr<float>&);
extern template ArrayXXc<float> vshRBchi(ArrayXr<float>, const ArrayXr<float>&);
extern template ArrayXXc<float> vshRBpsi(ArrayXr<float>, const ArrayXr<float>&);
*/

template unique_ptr<ArrayXXr<float>> sphGetUforFp(int);
template unique_ptr<stFprow<float>> sphGetFpRow(int, float, const ArrayXr<float>&);
template unique_ptr<stFpovx<float>> sphGetFpovx(int, float, const ArrayXr<float>&);
template unique_ptr<stBessel<float>> sphGetXiPsi(int, float, const ArrayXr<float>&, int);
template unique_ptr<stBesselPrimes<float>> sphGetBesselProductsPrimes(const Tensor3c<float>&);
template unique_ptr<stBesselProducts<float>> sphGetModifiedBesselProducts(int, float, const ArrayXr<float>&, int);
template unique_ptr<stRtfunc<float>> sphMakeGeometry(int, float, float, const unique_ptr<ArrayXr<float>>);
template int sphCheckBesselConvergence(int, float, const ArrayXr<float>&, float, int);
template int sphEstimateNB(int, const unique_ptr<stRtfunc<float>>&, const unique_ptr<stParams<float>>&, float);
template vector<unique_ptr<stPQ<float>>> sphCalculatePQ(int, const ArrayXi&,
    const unique_ptr<stRtfunc<float>>&, const unique_ptr<stParams<float>>&, int);

/*
extern template float mp_pi<float();

extern template ArrayXXc<float> invertLUcol(MatrixXc<float>&);
extern template ArrayXi logicalSlice(ArrayXi&, ArrayXb&);

extern template unique_ptr<stIncEabnm<float>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<float>>&);
*/

template vector<unique_ptr<stTR<float>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<float>>>&, bool);
template vector<unique_ptr<stTR<float>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<float>>>&, int);
template vector<unique_ptr<stTR<float>>> rvhGetSymmetricMat(const vector<unique_ptr<stTR<float>>>&, vector<string>);
template unique_ptr<stAbcdnm<float>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<float>>>&,
    const unique_ptr<stIncPar<float>>&, unique_ptr<stIncEabnm<float>>);
template unique_ptr<stCrossSection<float>> rvhGetAverageCrossSections(
    const ArrayXr<float>&, const vector<vector<unique_ptr<stTR<float>>>>&);

/*
extern template int sphEstimateNB(int, const unique_ptr<stRtfunc<float>>&,
    const unique_ptr<stParams<float>>&, float);
extern template vector<unique_ptr<stPQ<float>>> sphCalculatePQ(int, const ArrayXi&,
    const unique_ptr<stRtfunc<float>>&, const unique_ptr<stParams<float>>&, int);

extern template vector<unique_ptr<stTR<float>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<float>>>&, bool);
extern template vector<unique_ptr<stTR<float>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<float>>>&, int);
extern template vector<unique_ptr<stTR<float>>> rvhGetSymmetricMat(
    const vector<unique_ptr<stTR<float>>>&, vector<string>);
extern template unique_ptr<stCrossSection<float>> rvhGetAverageCrossSections(
    const ArrayXr<float>&, const vector<vector<unique_ptr<stTR<float>>>>&);
*/

template unique_ptr<stTmatrix<float>> slvForT(const unique_ptr<stParams<float>>&,
    const unique_ptr<stOptions>&, unique_ptr<stRtfunc<float>>);

/*
extern float mp_pi<float>();
extern float mp_eps<float>();
extern float mp_im_unit<float>();

extern template unique_ptr<stIncPar<float>> vshMakeIncidentParams(sIncType, int);
*/

template unique_ptr<stRes<float>> pstMakeStructForField(const unique_ptr<stAbcdnm<float>>&,
    int, ArrayXr<float>, ArrayXr<float>, float, unique_ptr<stIncPar<float>>, float, float);
template unique_ptr<stRes<float>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<float>>&, const unique_ptr<stParams<float>>&);
template unique_ptr<stSM<float>> pstScatteringMatrixOA(const vector<unique_ptr<stTR<float>>>&,
    float, float, int);

extern template double mp_pi<double>();
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
extern template unique_ptr<stParams<double>> Raman2SmartiesParams(
    const unique_ptr<RamanParams<double>>&, int, int, string);

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

template float Frac2Float(string);
template unique_ptr<RamanParams<float>> LoadParams<float>(string);
template unique_ptr<stParams<float>> Raman2SmartiesParams(const unique_ptr<RamanParams<float>>&, int, int, string);
template void CreateTimeStamp<float, float>(string, const unique_ptr<RamanParams<float>>&);
template void CreateTimeStamp<float, double>(string, const unique_ptr<RamanParams<float>>&);
template void CreateTimeStamp<float, long double>(string, const unique_ptr<RamanParams<float>>&);
template vector<unique_ptr<stTR<float>>> ConvertstTRList(const vector<unique_ptr<stTR<float>>>&);
template vector<unique_ptr<stTR<double>>> ConvertstTRList(const vector<unique_ptr<stTR<float>>>&);
template vector<unique_ptr<stTR<long double>>> ConvertstTRList(const vector<unique_ptr<stTR<float>>>&);
template void RamanElasticScattering<float, float>(string, string);
template void RamanElasticScattering<float, double>(string, string);
template void RamanElasticScattering<float, long double>(string, string);
