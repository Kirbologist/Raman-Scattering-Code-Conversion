#include "smarties.hpp"
#include "raman_elastic_scattering.hpp"

using namespace Eigen;
using namespace Smarties;

template double mp_pi<double>();
template double mp_eps<double>();
template complex<double> mp_im_unit<double>();

template Tensor3c<double> subtensor(Tensor3c<double>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
template ArrayXXc<double> reduceAndSlice(Tensor3c<double>&, int, int);
template ArrayXXc<double> invertLUcol(MatrixXc<double>&);
template ArrayXi logicalSlice(ArrayXi&, ArrayXb&);
template RowArrayXd logicalSlice(RowArrayXd&, RowArrayXb&);
template ArrayXd logicalSlice(ArrayXd&, ArrayXb&);
template Tensor4c<double> tensor_conj(Tensor4c<double>&);
template ArrayXd arr_bessel_j(ArrayXd&, double);
template ArrayXd arr_bessel_y(ArrayXd&, double);

/*
extern template double mp_pi<double>();
extern template double mp_eps<double>();
*/

template unique_ptr<stGLQuad<double>> auxInitLegendreQuad(int, double, double);
template unique_ptr<stRtfunc<double>> auxPrepareIntegrals(int, sInt);

/*
extern template double mp_pi<double>();
extern template double mp_im_unit<double>();

extern template ArrayXd arr_bessel_j(ArrayXd&, double);
extern template ArrayXd arr_bessel_y(ArrayXd&, double);
*/

template unique_ptr<stIncPar<double>> vshMakeIncidentParams(sIncType, int);
template unique_ptr<stIncPar<double>> vshMakeIncidentParams(sIncType, int, double, double, double);
template unique_ptr<stPinmTaunm<double>> vshPinmTaunm(int, const ArrayXd&);
template unique_ptr<stIncEabnm<double>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<double>>&);
template unique_ptr<stZnAll<double>> vshGetZnAll(int, const ArrayXd&, sBessel);
template unique_ptr<stEAllPhi<double>> vshEgenThetaAllPhi(const ArrayXd&,
    const ArrayXd&, const ArrayXXc<double>&, const ArrayXXc<double>&,
    const RowArrayXd&, const RowArrayXd&, sBessel, unique_ptr<stPinmTaunm<double>>);
template unique_ptr<stEforPhi<double>> vshEthetaForPhi(const unique_ptr<stEAllPhi<double>>&, double);
template ArrayXXc<double> vshRBchi(ArrayXd, const ArrayXd&);
template ArrayXXc<double> vshRBpsi(ArrayXd, const ArrayXd&);

/*
extern template mp_eps<double>();
extern template mp_im_unit<double>();

extern template Tensor3c<double> subtensor(Tensor3c<double>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
extern template ArrayXXc<double> reduceAndSlice(Tensor3c<double>&, int, int);
extern template RowArrayXd logicalSlice(RowArrayXd&, RowArrayXb&);

extern template unique_ptr<stRtfunc<double>> auxPrepareIntegrals(int, sInt);

extern template unique_ptr<stPinmTaunm<double>> vshPinmTaunm(int, const ArrayXd&);
extern template ArrayXXc<double> vshRBchi(ArrayXd, const ArrayXd&);
extern template ArrayXXc<double> vshRBpsi(ArrayXd, const ArrayXd&);
*/

template unique_ptr<ArrayXXd> sphGetUforFp(int);
template unique_ptr<stFprow<double>> sphGetFpRow(int, double, const ArrayXd&);
template unique_ptr<stFpovx<double>> sphGetFpovx(int, double, const ArrayXd&);
template unique_ptr<stBessel<double>> sphGetXiPsi(int, double, const ArrayXd&, int);
template unique_ptr<stBesselPrimes<double>> sphGetBesselProductsPrimes(const Tensor3c<double>&);
template unique_ptr<stBesselProducts<double>> sphGetModifiedBesselProducts(int, double, const ArrayXd&, int);
template unique_ptr<stRtfunc<double>> sphMakeGeometry(int, double, double, const unique_ptr<ArrayXd>);
template int sphCheckBesselConvergence(int, double, const ArrayXd&, double, int);
template int sphEstimateNB(int, const unique_ptr<stRtfunc<double>>&, const unique_ptr<stParams<double>>&, double);
template vector<unique_ptr<stPQ<double>>> sphCalculatePQ(int, const ArrayXi&,
    const unique_ptr<stRtfunc<double>>&, const unique_ptr<stParams<double>>&, int);

/*
extern template double mp_pi<double();

extern template ArrayXXc<double> invertLUcol(MatrixXc<double>&);
extern template ArrayXi logicalSlice(ArrayXi&, ArrayXb&);

extern template unique_ptr<stIncEabnm<double>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<double>>&);
*/

template vector<unique_ptr<stTR<double>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<double>>>&, bool);
template vector<unique_ptr<stTR<double>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<double>>>&, int);
template vector<unique_ptr<stTR<double>>> rvhGetSymmetricMat(const vector<unique_ptr<stTR<double>>>&, vector<string>);
template unique_ptr<stAbcdnm<double>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<double>>>&,
    const unique_ptr<stIncPar<double>>&, unique_ptr<stIncEabnm<double>>);
template unique_ptr<stCrossSection<double>> rvhGetAverageCrossSections(
    const ArrayXd&, const vector<vector<unique_ptr<stTR<double>>>>&);

/*
extern template int sphEstimateNB(int, const unique_ptr<stRtfunc<double>>&,
    const unique_ptr<stParams<double>>&, double);
extern template vector<unique_ptr<stPQ<double>>> sphCalculatePQ(int, const ArrayXi&,
    const unique_ptr<stRtfunc<double>>&, const unique_ptr<stParams<double>>&, int);

extern template vector<unique_ptr<stTR<double>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<double>>>&, bool);
extern template vector<unique_ptr<stTR<double>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<double>>>&, int);
extern template vector<unique_ptr<stTR<double>>> rvhGetSymmetricMat(
    const vector<unique_ptr<stTR<double>>>&, vector<string>);
extern template unique_ptr<stCrossSection<double>> rvhGetAverageCrossSections(
    const ArrayXd&, const vector<vector<unique_ptr<stTR<double>>>>&);
*/

template unique_ptr<stTmatrix<double>> slvForT(const unique_ptr<stParams<double>>&,
    const unique_ptr<stOptions>&, unique_ptr<stRtfunc<double>>);

/*
extern double mp_pi<double>();
extern double mp_eps<double>();
extern double mp_im_unit<double>();

extern template unique_ptr<stIncPar<double>> vshMakeIncidentParams(sIncType, int);
*/

template unique_ptr<stRes<double>> pstMakeStructForField(const unique_ptr<stAbcdnm<double>>&,
    int, ArrayXd, ArrayXd, double, unique_ptr<stIncPar<double>>, double, double);
template unique_ptr<stRes<double>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<double>>&, const unique_ptr<stParams<double>>&);
template unique_ptr<stSM<double>> pstScatteringMatrixOA(const vector<unique_ptr<stTR<double>>>&,
    double, double, int);

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

template double Frac2Float(string);
template unique_ptr<RamanParams<double>> LoadParams<double>(string);
template unique_ptr<stParams<double>> Raman2SmartiesParams(const unique_ptr<RamanParams<double>>&, int, int, string);
template void CreateTimeStamp<double, float>(string, const unique_ptr<RamanParams<double>>&);
template void CreateTimeStamp<double, double>(string, const unique_ptr<RamanParams<double>>&);
template void CreateTimeStamp<double, long double>(string, const unique_ptr<RamanParams<double>>&);
template vector<unique_ptr<stTR<float>>> ConvertstTRList(const vector<unique_ptr<stTR<double>>>&);
template vector<unique_ptr<stTR<double>>> ConvertstTRList(const vector<unique_ptr<stTR<double>>>&);
template vector<unique_ptr<stTR<long double>>> ConvertstTRList(const vector<unique_ptr<stTR<double>>>&);
template void RamanElasticScattering<double, float>(string, string);
template void RamanElasticScattering<double, double>(string, string);
template void RamanElasticScattering<double, long double>(string, string);
