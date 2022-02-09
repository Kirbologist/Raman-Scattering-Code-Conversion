/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Based on the SMARTIES MATLAB package by W.R.C. Somerville, B. Auguié, E.C. Le Ru
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.


This code is an instantiation of all template functions using template argument 'double'.
In the case where there are two template arguments, such functions are always used in RamanElastiScattering,
where two different types are used for different parts of the calculation.
Thus, the only instantiations instantiated for those functions are for when the
template arguemnt `Real1` of RamanElasticScattering is double.

Instantiations are split into segments based on which file they are from.
Before each segment is a set of commented out external instantiations; these would've
instantiated all the dependencies of the segment. These are written so that if any of the
segments were to be instantiated on their own, these external instantiations can be used with the segment.
*/

#include "smarties.hpp"
#include "raman_elastic_scattering.hpp"

using namespace Eigen;
using namespace Smarties;

template double mp_pi<double>();
template double mp_eps<double>();
template complex<double> mp_im_unit<double>();

// TODO: instantiate ArrayMap and TensorCast functions for both Real and Complex
template string GetTypeName<double>();
template Tensor3c<double> TensorSlice(Tensor3c<double>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
template ArrayXXc<double> Subtensor2ArrMap(const Tensor3c<double>&,
    const std::array<int, 3>&, const std::array<int, 3>&, int, int);
template ArrayXXc<double> InvertLUcol(MatrixXc<double>&);
template ArrayXi LogicalSlice(ArrayXi&, ArrayXb&);
template RowArrayXd LogicalSlice(RowArrayXd&, RowArrayXb&);
template ArrayXd LogicalSlice(ArrayXd&, ArrayXb&);
template Tensor4c<double> TensorConj(Tensor4c<double>&);
template ArrayXd ArrBesselJ(ArrayXd&, double);
template ArrayXd ArrBesselY(ArrayXd&, double);

/*
extern template double mp_pi<double>();
extern template double mp_eps<double>();
*/

template unique_ptr<stGLQuad<double>> auxInitLegendreQuad(int, double, double);
extern template void WriteBinary(const string, const ArrayXi&);
template void WriteBinary(const string, const ArrayXXd&);
extern template void ReadBinary(const string, ArrayXi&);
template void ReadBinary(const string, ArrayXXd&);
template void UpdateGLquadrature<double>(ArrayXi, bool);
template void StoreGLquadrature<double>();
template unique_ptr<stRtfunc<double>> auxPrepareIntegrals(int, sInt);

/*
extern template double mp_pi<double>();
extern template double mp_im_unit<double>();

extern template ArrayXd ArrBesselJ(ArrayXd&, double);
extern template ArrayXd ArrBesselY(ArrayXd&, double);
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

extern template Tensor3c<double> TensorSlice(Tensor3c<double>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
extern template ArrayXXc<double> ReduceAndSlice(Tensor3c<double>&, int, int);
extern template RowArrayXd LogicalSlice(RowArrayXd&, RowArrayXb&);

extern template unique_ptr<stRtfunc<double>> auxPrepareIntegrals(int, sInt);

extern template unique_ptr<stPinmTaunm<double>> vshPinmTaunm(int, const ArrayXd&);
extern template ArrayXXc<double> vshRBchi(ArrayXd, const ArrayXd&);
extern template ArrayXXc<double> vshRBpsi(ArrayXd, const ArrayXd&);
*/

template ArrayXXd sphGetUforFp(int);
template unique_ptr<stFpRow<double>> sphGetFpRow(int, double, const ArrayXd&);
template unique_ptr<stFpovx<double>> sphGetFpovx(int, double, const ArrayXd&);
template unique_ptr<stBessel<double>> sphGetXiPsi(int, double, const ArrayXd&, int);
template unique_ptr<stBesselPrimes<double>> sphGetBesselProductsPrimes(const Tensor3c<double>&);
template unique_ptr<stBesselProducts<double>> sphGetModifiedBesselProducts(int, double, const ArrayXd&, int);
template unique_ptr<stRtfunc<double>> sphMakeGeometry(int, double, double, const unique_ptr<ArrayXd>);
template int sphCheckBesselConvergence(int, double, const ArrayXd&, double, int);
template int sphEstimateNB(int, const unique_ptr<stRtfunc<double>>&, const unique_ptr<stParams<double>>&, double);
template vector<unique_ptr<stPQ<double>>> sphCalculatePQ(int, const ArrayXi&,
    const unique_ptr<stRtfunc<double>>&, const unique_ptr<stParams<double>>&, int);
template stDelta<double> sphEstimateDelta(
    const unique_ptr<stRtfunc<double>>&, const unique_ptr<stParams<double>>&, int);

/*
extern template double mp_pi<double();

extern template ArrayXXc<double> InvertLUcol(MatrixXc<double>&);
extern template ArrayXi LogicalSlice(ArrayXi&, ArrayXb&);

extern template unique_ptr<stIncEabnm<double>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<double>>&);
*/

template vector<unique_ptr<stTR<double>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<double>>>&, bool);
template vector<unique_ptr<stPQ<double>>> rvhTruncateMatrices(const vector<unique_ptr<stPQ<double>>>&, int);
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
extern template stDelta<double> sphEstimateDelta(unique_ptr<stRtfunc<double>>, unique_ptr<stParams<double>>, int);

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

extern template Tensor4c<float> TensorConj(Tensor4c<float>&);

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

extern template Tensor4c<long double> TensorConj(Tensor4c<long double>&);

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
template vector<unique_ptr<stTR<float>>> ConvertStTRList(const vector<unique_ptr<stTR<double>>>&);
template vector<unique_ptr<stTR<double>>> ConvertStTRList(const vector<unique_ptr<stTR<double>>>&);
template vector<unique_ptr<stTR<long double>>> ConvertStTRList(const vector<unique_ptr<stTR<double>>>&);
template void RamanElasticScattering<double, float>(string, string);
template void RamanElasticScattering<double, double>(string, string);
template void RamanElasticScattering<double, long double>(string, string);
