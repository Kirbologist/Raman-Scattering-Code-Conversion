#include "../src/smarties.hpp"
#include "../src/raman_elastic_scattering.hpp"
#include "core_mp.hpp"

using namespace Eigen;
using namespace Smarties;

template RamanFloat mp_pi<RamanFloat>();
template RamanFloat mp_eps<RamanFloat>();
template complex<RamanFloat> mp_im_unit<RamanFloat>();

template Tensor3c<RamanFloat> tensorSlice(Tensor3c<RamanFloat>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
template Map<ArrayXXc<RamanFloat>> subtensor2ArrMap(const Tensor3c<RamanFloat>&,
    const std::array<int, 3>&, const std::array<int, 3>&, int, int);
template ArrayXXc<RamanFloat> invertLUcol(MatrixXc<RamanFloat>&);
template ArrayXi logicalSlice(ArrayXi&, ArrayXb&);
template RowArrayXr<RamanFloat> logicalSlice(RowArrayXr<RamanFloat>&, RowArrayXb&);
template ArrayXr<RamanFloat> logicalSlice(ArrayXr<RamanFloat>&, ArrayXb&);
template Tensor4c<RamanFloat> tensor_conj(Tensor4c<RamanFloat>&);
template ArrayXr<RamanFloat> arr_bessel_j(ArrayXr<RamanFloat>&, RamanFloat);
template ArrayXr<RamanFloat> arr_bessel_y(ArrayXr<RamanFloat>&, RamanFloat);

/*
extern template RamanFloat mp_pi<RamanFloat>();
extern template RamanFloat mp_eps<RamanFloat>();
*/

template unique_ptr<stGLQuad<RamanFloat>> auxInitLegendreQuad(int, RamanFloat, RamanFloat);
template void writeBinary(const string, const ArrayXXr<RamanFloat>&);
template void readBinary(const string, ArrayXXr<RamanFloat>&);
template string getTypeName<RamanFloat>();
template void updateGLquadrature<RamanFloat>(ArrayXi, bool);
template void storeGLquadrature<RamanFloat>();
template unique_ptr<stRtfunc<RamanFloat>> auxPrepareIntegrals(int, sInt);

/*
extern template RamanFloat mp_pi<RamanFloat>();
extern template RamanFloat mp_im_unit<RamanFloat>();

extern template ArrayXr<RamanFloat> arr_bessel_j(ArrayXr<RamanFloat>&, RamanFloat);
extern template ArrayXr<RamanFloat> arr_bessel_y(ArrayXr<RamanFloat>&, RamanFloat);
*/

template unique_ptr<stIncPar<RamanFloat>> vshMakeIncidentParams(sIncType, int);
template unique_ptr<stIncPar<RamanFloat>> vshMakeIncidentParams(sIncType, int, RamanFloat, RamanFloat, RamanFloat);
template unique_ptr<stPinmTaunm<RamanFloat>> vshPinmTaunm(int, const ArrayXr<RamanFloat>&);
template unique_ptr<stIncEabnm<RamanFloat>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<RamanFloat>>&);
template unique_ptr<stZnAll<RamanFloat>> vshGetZnAll(int, const ArrayXr<RamanFloat>&, sBessel);
template unique_ptr<stEAllPhi<RamanFloat>> vshEgenThetaAllPhi(const ArrayXr<RamanFloat>&,
    const ArrayXr<RamanFloat>&, const ArrayXXc<RamanFloat>&, const ArrayXXc<RamanFloat>&,
    const RowArrayXr<RamanFloat>&, const RowArrayXr<RamanFloat>&, sBessel, unique_ptr<stPinmTaunm<RamanFloat>>);
template unique_ptr<stEforPhi<RamanFloat>> vshEthetaForPhi(const unique_ptr<stEAllPhi<RamanFloat>>&, RamanFloat);
template ArrayXXc<RamanFloat> vshRBchi(ArrayXr<RamanFloat>, const ArrayXr<RamanFloat>&);
template ArrayXXc<RamanFloat> vshRBpsi(ArrayXr<RamanFloat>, const ArrayXr<RamanFloat>&);

/*
extern template mp_eps<RamanFloat>();
extern template mp_im_unit<RamanFloat>();

extern template Tensor3c<RamanFloat> subtensor(Tensor3c<RamanFloat>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
extern template ArrayXXc<RamanFloat> reduceAndSlice(Tensor3c<RamanFloat>&, int, int);
extern template RowArrayXr<RamanFloat> logicalSlice(RowArrayXr<RamanFloat>&, RowArrayXb&);

extern template unique_ptr<stRtfunc<RamanFloat>> auxPrepareIntegrals(int, sInt);

extern template unique_ptr<stPinmTaunm<RamanFloat>> vshPinmTaunm(int, const ArrayXr<RamanFloat>&);
extern template ArrayXXc<RamanFloat> vshRBchi(ArrayXr<RamanFloat>, const ArrayXr<RamanFloat>&);
extern template ArrayXXc<RamanFloat> vshRBpsi(ArrayXr<RamanFloat>, const ArrayXr<RamanFloat>&);
*/

template unique_ptr<ArrayXXr<RamanFloat>> sphGetUforFp(int);
template unique_ptr<stFprow<RamanFloat>> sphGetFpRow(int, RamanFloat, const ArrayXr<RamanFloat>&);
template unique_ptr<stFpovx<RamanFloat>> sphGetFpovx(int, RamanFloat, const ArrayXr<RamanFloat>&);
template unique_ptr<stBessel<RamanFloat>> sphGetXiPsi(int, RamanFloat, const ArrayXr<RamanFloat>&, int);
template unique_ptr<stBesselPrimes<RamanFloat>> sphGetBesselProductsPrimes(const Tensor3c<RamanFloat>&);
template unique_ptr<stBesselProducts<RamanFloat>> sphGetModifiedBesselProducts(int, RamanFloat, const ArrayXr<RamanFloat>&, int);
template unique_ptr<stRtfunc<RamanFloat>> sphMakeGeometry(int, RamanFloat, RamanFloat, const unique_ptr<ArrayXr<RamanFloat>>);
template int sphCheckBesselConvergence(int, RamanFloat, const ArrayXr<RamanFloat>&, RamanFloat, int);
template int sphEstimateNB(int, const unique_ptr<stRtfunc<RamanFloat>>&, const unique_ptr<stParams<RamanFloat>>&, RamanFloat);
template vector<unique_ptr<stPQ<RamanFloat>>> sphCalculatePQ(int, const ArrayXi&,
    const unique_ptr<stRtfunc<RamanFloat>>&, const unique_ptr<stParams<RamanFloat>>&, int);

/*
extern template RamanFloat mp_pi<RamanFloat();

extern template ArrayXXc<RamanFloat> invertLUcol(MatrixXc<RamanFloat>&);
extern template ArrayXi logicalSlice(ArrayXi&, ArrayXb&);

extern template unique_ptr<stIncEabnm<RamanFloat>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<RamanFloat>>&);
*/

template vector<unique_ptr<stTR<RamanFloat>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<RamanFloat>>>&, bool);
template vector<unique_ptr<stTR<RamanFloat>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<RamanFloat>>>&, int);
template vector<unique_ptr<stTR<RamanFloat>>> rvhGetSymmetricMat(const vector<unique_ptr<stTR<RamanFloat>>>&, vector<string>);
template unique_ptr<stAbcdnm<RamanFloat>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<RamanFloat>>>&,
    const unique_ptr<stIncPar<RamanFloat>>&, unique_ptr<stIncEabnm<RamanFloat>>);
template unique_ptr<stCrossSection<RamanFloat>> rvhGetAverageCrossSections(
    const ArrayXr<RamanFloat>&, const vector<vector<unique_ptr<stTR<RamanFloat>>>>&);

/*
extern template int sphEstimateNB(int, const unique_ptr<stRtfunc<RamanFloat>>&,
    const unique_ptr<stParams<RamanFloat>>&, RamanFloat);
extern template vector<unique_ptr<stPQ<RamanFloat>>> sphCalculatePQ(int, const ArrayXi&,
    const unique_ptr<stRtfunc<RamanFloat>>&, const unique_ptr<stParams<RamanFloat>>&, int);

extern template vector<unique_ptr<stTR<RamanFloat>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<RamanFloat>>>&, bool);
extern template vector<unique_ptr<stTR<RamanFloat>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<RamanFloat>>>&, int);
extern template vector<unique_ptr<stTR<RamanFloat>>> rvhGetSymmetricMat(
    const vector<unique_ptr<stTR<RamanFloat>>>&, vector<string>);
extern template unique_ptr<stCrossSection<RamanFloat>> rvhGetAverageCrossSections(
    const ArrayXr<RamanFloat>&, const vector<vector<unique_ptr<stTR<RamanFloat>>>>&);
*/

template unique_ptr<stTmatrix<RamanFloat>> slvForT(const unique_ptr<stParams<RamanFloat>>&,
    const unique_ptr<stOptions>&, unique_ptr<stRtfunc<RamanFloat>>);

/*
extern RamanFloat mp_pi<RamanFloat>();
extern RamanFloat mp_eps<RamanFloat>();
extern RamanFloat mp_im_unit<RamanFloat>();

extern template unique_ptr<stIncPar<RamanFloat>> vshMakeIncidentParams(sIncType, int);
*/

template unique_ptr<stRes<RamanFloat>> pstMakeStructForField(const unique_ptr<stAbcdnm<RamanFloat>>&,
    int, ArrayXr<RamanFloat>, ArrayXr<RamanFloat>, RamanFloat, unique_ptr<stIncPar<RamanFloat>>, RamanFloat, RamanFloat);
template unique_ptr<stRes<RamanFloat>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<RamanFloat>>&, const unique_ptr<stParams<RamanFloat>>&);
template unique_ptr<stSM<RamanFloat>> pstScatteringMatrixOA(const vector<unique_ptr<stTR<RamanFloat>>>&,
    RamanFloat, RamanFloat, int);

/*
extern template RamanFloat mp_pi<RamanFloat>();

extern template Tensor4c<RamanFloat> tensor_conj(Tensor4c<RamanFloat>&);

extern template unique_ptr<stIncPar<RamanFloat>> vshMakeIncidentParams(sIncType, int, RamanFloat, RamanFloat, RamanFloat);
extern template unique_ptr<stEAllPhi<RamanFloat>> vshEgenThetaAllPhi(const ArrayXr<RamanFloat>&,
    const ArrayXr<RamanFloat>&, const ArrayXXc<RamanFloat>&, const ArrayXXc<RamanFloat>&,
    const RowArrayXr<RamanFloat>&, const RowArrayXr<RamanFloat>&, sBessel, unique_ptr<stPinmTaunm<RamanFloat>>);
extern template unique_ptr<stEforPhi<RamanFloat>> vshEthetaForPhi(const unique_ptr<stEAllPhi<RamanFloat>>&, RamanFloat);

extern template unique_ptr<stAbcdnm<RamanFloat>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<RamanFloat>>>&,
    const unique_ptr<stIncPar<RamanFloat>>&, unique_ptr<stIncEabnm<RamanFloat>>);

extern template unique_ptr<stTmatrix<RamanFloat>> slvForT(const unique_ptr<stParams<RamanFloat>>&,
    const unique_ptr<stOptions>&, unique_ptr<stRtfunc<RamanFloat>>);

extern template unique_ptr<stRes<RamanFloat>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<RamanFloat>>&, const unique_ptr<stParams<RamanFloat>>&);
extern template unique_ptr<stSM<RamanFloat>> pstScatteringMatrixOA(
    const vector<unique_ptr<stTR<RamanFloat>>>&, RamanFloat, RamanFloat, int);
*/

extern template float mp_pi();
extern template complex<float> mp_im_unit();

extern template Tensor4c<float> tensor_conj(Tensor4c<float>&);

extern template string getTypeName<float>();

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

extern template string getTypeName<double>();

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

extern template string getTypeName<long double>();

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

template RamanFloat Frac2Float(string);
template unique_ptr<RamanParams<RamanFloat>> LoadParams<RamanFloat>(string);
template unique_ptr<stParams<RamanFloat>> Raman2SmartiesParams(
    const unique_ptr<RamanParams<RamanFloat>>&, int, int, string);
template void CreateTimeStamp<RamanFloat, float>(string, const unique_ptr<RamanParams<RamanFloat>>&);
template void CreateTimeStamp<RamanFloat, double>(string, const unique_ptr<RamanParams<RamanFloat>>&);
template void CreateTimeStamp<RamanFloat, long double>(string, const unique_ptr<RamanParams<RamanFloat>>&);
template void CreateTimeStamp<RamanFloat, RamanFloat>(string, const unique_ptr<RamanParams<RamanFloat>>&);
template vector<unique_ptr<stTR<float>>> ConvertstTRList(const vector<unique_ptr<stTR<RamanFloat>>>&);
template vector<unique_ptr<stTR<double>>> ConvertstTRList(const vector<unique_ptr<stTR<RamanFloat>>>&);
template vector<unique_ptr<stTR<long double>>> ConvertstTRList(const vector<unique_ptr<stTR<RamanFloat>>>&);
template vector<unique_ptr<stTR<RamanFloat>>> ConvertstTRList(const vector<unique_ptr<stTR<RamanFloat>>>&);
template void RamanElasticScattering<RamanFloat, float>(string, string, int);
template void RamanElasticScattering<RamanFloat, double>(string, string, int);
template void RamanElasticScattering<RamanFloat, long double>(string, string, int);
template void RamanElasticScattering<RamanFloat, RamanFloat>(string, string, int);

template void CreateTimeStamp<float, RamanFloat>(string, const unique_ptr<RamanParams<float>>&);
template void CreateTimeStamp<double, RamanFloat>(string, const unique_ptr<RamanParams<double>>&);
template void CreateTimeStamp<long double, RamanFloat>(string, const unique_ptr<RamanParams<long double>>&);
template vector<unique_ptr<stTR<RamanFloat>>> ConvertstTRList(const vector<unique_ptr<stTR<float>>>&);
template vector<unique_ptr<stTR<RamanFloat>>> ConvertstTRList(const vector<unique_ptr<stTR<double>>>&);
template vector<unique_ptr<stTR<RamanFloat>>> ConvertstTRList(const vector<unique_ptr<stTR<long double>>>&);
template void RamanElasticScattering<float, RamanFloat>(string, string, int);
template void RamanElasticScattering<double, RamanFloat>(string, string, int);
template void RamanElasticScattering<long double, RamanFloat>(string, string, int);
