#include "smarties.hpp"
#include "raman_elastic_scattering.hpp"

using namespace Eigen;
using namespace Smarties;

template long double mp_pi<long double>();
template long double mp_eps<long double>();
template complex<long double> mp_im_unit<long double>();

template Tensor3c<long double> subtensor(Tensor3c<long double>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
template ArrayXXc<long double> reduceAndSlice(Tensor3c<long double>&, int, int);
template ArrayXXc<long double> invertLUcol(MatrixXc<long double>&);
template ArrayXi logicalSlice(ArrayXi&, ArrayXb&);
template RowArrayXr<long double> logicalSlice(RowArrayXr<long double>&, RowArrayXb&);
template ArrayXr<long double> logicalSlice(ArrayXr<long double>&, ArrayXb&);
template Tensor4c<long double> tensor_conj(Tensor4c<long double>&);
template ArrayXr<long double> arr_bessel_j(ArrayXr<long double>&, long double);
template ArrayXr<long double> arr_bessel_y(ArrayXr<long double>&, long double);

/*
extern template long double mp_pi<long double>();
extern template long double mp_eps<long double>();
*/

template unique_ptr<stGLQuad<long double>> auxInitLegendreQuad(int, long double, long double);
template unique_ptr<stRtfunc<long double>> auxPrepareIntegrals(int, sInt);

/*
extern template long double mp_pi<long double>();
extern template long double mp_im_unit<long double>();

extern template ArrayXr<long double> arr_bessel_j(ArrayXr<long double>&, long double);
extern template ArrayXr<long double> arr_bessel_y(ArrayXr<long double>&, long double);
*/

template unique_ptr<stIncPar<long double>> vshMakeIncidentParams(sIncType, int);
template unique_ptr<stIncPar<long double>> vshMakeIncidentParams(sIncType, int, long double, long double, long double);
template unique_ptr<stPinmTaunm<long double>> vshPinmTaunm(int, const ArrayXr<long double>&);
template unique_ptr<stIncEabnm<long double>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<long double>>&);
template unique_ptr<stZnAll<long double>> vshGetZnAll(int, const ArrayXr<long double>&, sBessel);
template unique_ptr<stEAllPhi<long double>> vshEgenThetaAllPhi(const ArrayXr<long double>&,
    const ArrayXr<long double>&, const ArrayXXc<long double>&, const ArrayXXc<long double>&,
    const RowArrayXr<long double>&, const RowArrayXr<long double>&, sBessel, unique_ptr<stPinmTaunm<long double>>);
template unique_ptr<stEforPhi<long double>> vshEthetaForPhi(const unique_ptr<stEAllPhi<long double>>&, long double);
template ArrayXXc<long double> vshRBchi(ArrayXr<long double>, const ArrayXr<long double>&);
template ArrayXXc<long double> vshRBpsi(ArrayXr<long double>, const ArrayXr<long double>&);

/*
extern template mp_eps<long double>();
extern template mp_im_unit<long double>();

extern template Tensor3c<long double> subtensor(Tensor3c<long double>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
extern template ArrayXXc<long double> reduceAndSlice(Tensor3c<long double>&, int, int);
extern template RowArrayXr<long double> logicalSlice(RowArrayXr<long double>&, RowArrayXb&);

extern template unique_ptr<stRtfunc<long double>> auxPrepareIntegrals(int, sInt);

extern template unique_ptr<stPinmTaunm<long double>> vshPinmTaunm(int, const ArrayXr<long double>&);
extern template ArrayXXc<long double> vshRBchi(ArrayXr<long double>, const ArrayXr<long double>&);
extern template ArrayXXc<long double> vshRBpsi(ArrayXr<long double>, const ArrayXr<long double>&);
*/

template unique_ptr<ArrayXXr<long double>> sphGetUforFp(int);
template unique_ptr<stFprow<long double>> sphGetFpRow(int, long double, const ArrayXr<long double>&);
template unique_ptr<stFpovx<long double>> sphGetFpovx(int, long double, const ArrayXr<long double>&);
template unique_ptr<stBessel<long double>> sphGetXiPsi(int, long double, const ArrayXr<long double>&, int);
template unique_ptr<stBesselPrimes<long double>> sphGetBesselProductsPrimes(const Tensor3c<long double>&);
template unique_ptr<stBesselProducts<long double>> sphGetModifiedBesselProducts(int, long double, const ArrayXr<long double>&, int);
template unique_ptr<stRtfunc<long double>> sphMakeGeometry(int, long double, long double, const unique_ptr<ArrayXr<long double>>);
template int sphCheckBesselConvergence(int, long double, const ArrayXr<long double>&, long double, int);
template int sphEstimateNB(int, const unique_ptr<stRtfunc<long double>>&, const unique_ptr<stParams<long double>>&, long double);
template vector<unique_ptr<stPQ<long double>>> sphCalculatePQ(int, const ArrayXi&,
    const unique_ptr<stRtfunc<long double>>&, const unique_ptr<stParams<long double>>&, int);

/*
extern template long double mp_pi<long double();

extern template ArrayXXc<long double> invertLUcol(MatrixXc<long double>&);
extern template ArrayXi logicalSlice(ArrayXi&, ArrayXb&);

extern template unique_ptr<stIncEabnm<long double>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<long double>>&);
*/

template vector<unique_ptr<stTR<long double>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<long double>>>&, bool);
template vector<unique_ptr<stTR<long double>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<long double>>>&, int);
template vector<unique_ptr<stTR<long double>>> rvhGetSymmetricMat(const vector<unique_ptr<stTR<long double>>>&, vector<string>);
template unique_ptr<stAbcdnm<long double>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<long double>>>&,
    const unique_ptr<stIncPar<long double>>&, unique_ptr<stIncEabnm<long double>>);
template unique_ptr<stCrossSection<long double>> rvhGetAverageCrossSections(
    const ArrayXr<long double>&, const vector<vector<unique_ptr<stTR<long double>>>>&);

/*
extern template int sphEstimateNB(int, const unique_ptr<stRtfunc<long double>>&,
    const unique_ptr<stParams<long double>>&, long double);
extern template vector<unique_ptr<stPQ<long double>>> sphCalculatePQ(int, const ArrayXi&,
    const unique_ptr<stRtfunc<long double>>&, const unique_ptr<stParams<long double>>&, int);

extern template vector<unique_ptr<stTR<long double>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<long double>>>&, bool);
extern template vector<unique_ptr<stTR<long double>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<long double>>>&, int);
extern template vector<unique_ptr<stTR<long double>>> rvhGetSymmetricMat(
    const vector<unique_ptr<stTR<long double>>>&, vector<string>);
extern template unique_ptr<stCrossSection<long double>> rvhGetAverageCrossSections(
    const ArrayXr<long double>&, const vector<vector<unique_ptr<stTR<long double>>>>&);
*/

template unique_ptr<stTmatrix<long double>> slvForT(const unique_ptr<stParams<long double>>&,
    const unique_ptr<stOptions>&, unique_ptr<stRtfunc<long double>>);

/*
extern long double mp_pi<long double>();
extern long double mp_eps<long double>();
extern long double mp_im_unit<long double>();

extern template unique_ptr<stIncPar<long double>> vshMakeIncidentParams(sIncType, int);
*/

template unique_ptr<stRes<long double>> pstMakeStructForField(const unique_ptr<stAbcdnm<long double>>&,
    int, ArrayXr<long double>, ArrayXr<long double>, long double, unique_ptr<stIncPar<long double>>, long double, long double);
template unique_ptr<stRes<long double>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<long double>>&, const unique_ptr<stParams<long double>>&);
template unique_ptr<stSM<long double>> pstScatteringMatrixOA(const vector<unique_ptr<stTR<long double>>>&,
    long double, long double, int);

/*
extern template long double mp_pi<long double>();

extern template Tensor4c<long double> tensor_conj(Tensor4c<long double>&);

extern template unique_ptr<stIncPar<long double>> vshMakeIncidentParams(sIncType, int, long double, long double, long double);
extern template unique_ptr<stEAllPhi<long double>> vshEgenThetaAllPhi(const ArrayXr<long double>&,
    const ArrayXr<long double>&, const ArrayXXc<long double>&, const ArrayXXc<long double>&,
    const RowArrayXr<long double>&, const RowArrayXr<long double>&, sBessel, unique_ptr<stPinmTaunm<long double>>);
extern template unique_ptr<stEforPhi<long double>> vshEthetaForPhi(const unique_ptr<stEAllPhi<long double>>&, long double);

extern template unique_ptr<stAbcdnm<long double>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<long double>>>&,
    const unique_ptr<stIncPar<long double>>&, unique_ptr<stIncEabnm<long double>>);

extern template unique_ptr<stTmatrix<long double>> slvForT(const unique_ptr<stParams<long double>>&,
    const unique_ptr<stOptions>&, unique_ptr<stRtfunc<long double>>);

extern template unique_ptr<stRes<long double>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<long double>>&, const unique_ptr<stParams<long double>>&);
extern template unique_ptr<stSM<long double>> pstScatteringMatrixOA(
    const vector<unique_ptr<stTR<long double>>>&, long double, long double, int);
*/
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

template unique_ptr<stParams<long double>> loadParam(string);
template unique_ptr<RamanParams<long double>> GetRamanParams(string in_file_name);
template void RamanElasticScattering<long double>(string, string);

extern template complex<double> mp_im_unit();
template vector<unique_ptr<stTR<double>>> stTRListMp2Double(const vector<unique_ptr<stTR<long double>>>&);
template void RamanElasticScatteringMpDouble<long double>(string, string);