/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Based on the SMARTIES MATLAB package by W.R.C. Somerville, B. Auguié, E.C. Le Ru
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.


This code is an instantiation of all template functions using
template argument 'float', i.e. single-precision floating-points.
In the case where there are two template arguments, such functions are always used in RamanElastiScattering,
where two different types are used for different parts of the calculation.
Thus, the only instantiations instantiated for those functions are for when the
template arguemnt `Real1` of RamanElasticScattering is float.

Instantiations are split into segments based on which file they are from.
Before each segment is a set of commented out external instantiations; these would've
instantiated all the dependencies of the segment. These are written so that if any of the
segments were to be instantiated on their own, these external instantiations can be used with the segment.
*/

#include "smarties.hpp"
#include "raman_elastic_scattering.hpp"

using namespace Eigen;
using namespace Smarties;

template float mp_pi<float>();
template float mp_eps<float>();
template complex<float> mp_im_unit<float>();

template Tensor3c<float> TensorSlice(Tensor3c<float>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
template ArrayXXc<float> Subtensor2ArrMap(const Tensor3c<float>&,
    const std::array<int, 3>&, const std::array<int, 3>&, int, int);
template ArrayXXc<float> InvertLUcol(MatrixXc<float>&);
template ArrayXi LogicalSlice(ArrayXi&, ArrayXb&);
template RowArrayXf LogicalSlice(RowArrayXf&, RowArrayXb&);
template ArrayXf LogicalSlice(ArrayXf&, ArrayXb&);
template Tensor4c<float> TensorConj(Tensor4c<float>&);
template ArrayXf ArrBesselJ(ArrayXf&, float);
template ArrayXf ArrBesselY(ArrayXf&, float);

/*
extern template float mp_pi<float>();
extern template float mp_eps<float>();
*/

template unique_ptr<stGLQuad<float>> auxInitLegendreQuad(int, float, float);
template unique_ptr<stRtfunc<float>> auxPrepareIntegrals(int, sInt);

/*
extern template float mp_pi<float>();
extern template float mp_im_unit<float>();

extern template ArrayXf ArrBesselJ(ArrayXf&, float);
extern template ArrayXf ArrBesselY(ArrayXf&, float);
*/

template unique_ptr<stIncPar<float>> vshMakeIncidentParams(sIncType, int);
template unique_ptr<stIncPar<float>> vshMakeIncidentParams(sIncType, int, float, float, float);
template unique_ptr<stPinmTaunm<float>> vshPinmTaunm(int, const ArrayXf&);
template unique_ptr<stIncEabnm<float>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<float>>&);
template unique_ptr<stZnAll<float>> vshGetZnAll(int, const ArrayXf&, sBessel);
template unique_ptr<stEAllPhi<float>> vshEgenThetaAllPhi(const ArrayXf&,
    const ArrayXf&, const ArrayXXc<float>&, const ArrayXXc<float>&,
    const RowArrayXf&, const RowArrayXf&, sBessel, unique_ptr<stPinmTaunm<float>>);
template unique_ptr<stEforPhi<float>> vshEthetaForPhi(const unique_ptr<stEAllPhi<float>>&, float);
template ArrayXXc<float> vshRBchi(ArrayXf, const ArrayXf&);
template ArrayXXc<float> vshRBpsi(ArrayXf, const ArrayXf&);

/*
extern template mp_eps<float>();
extern template mp_im_unit<float>();

extern template Tensor3c<float> TensorSlice(Tensor3c<float>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
extern template ArrayXXc<float> ReduceAndSlice(Tensor3c<float>&, int, int);
extern template RowArrayXf LogicalSlice(RowArrayXf&, RowArrayXb&);

extern template unique_ptr<stRtfunc<float>> auxPrepareIntegrals(int, sInt);

extern template unique_ptr<stPinmTaunm<float>> vshPinmTaunm(int, const ArrayXf&);
extern template ArrayXXc<float> vshRBchi(ArrayXf, const ArrayXf&);
extern template ArrayXXc<float> vshRBpsi(ArrayXf, const ArrayXf&);
*/

template ArrayXXf sphGetUforFp(int);
template unique_ptr<stFpRow<float>> sphGetFpRow(int, float, const ArrayXf&);
template unique_ptr<stFpovx<float>> sphGetFpovx(int, float, const ArrayXf&);
template unique_ptr<stBessel<float>> sphGetXiPsi(int, float, const ArrayXf&, int);
template unique_ptr<stBesselPrimes<float>> sphGetBesselProductsPrimes(const Tensor3c<float>&);
template unique_ptr<stBesselProducts<float>> sphGetModifiedBesselProducts(int, float, const ArrayXf&, int);
template unique_ptr<stRtfunc<float>> sphMakeGeometry(int, float, float, const unique_ptr<ArrayXf>);
template int sphCheckBesselConvergence(int, float, const ArrayXf&, float, int);
template int sphEstimateNB(int, const unique_ptr<stRtfunc<float>>&, const unique_ptr<stParams<float>>&, float);
template vector<unique_ptr<stPQ<float>>> sphCalculatePQ(int, const ArrayXi&,
    const unique_ptr<stRtfunc<float>>&, const unique_ptr<stParams<float>>&, int);

/*
extern template float mp_pi<float();

extern template ArrayXXc<float> InvertLUcol(MatrixXc<float>&);
extern template ArrayXi LogicalSlice(ArrayXi&, ArrayXb&);

extern template unique_ptr<stIncEabnm<float>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<float>>&);
*/

template vector<unique_ptr<stTR<float>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<float>>>&, bool);
template vector<unique_ptr<stTR<float>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<float>>>&, int);
template vector<unique_ptr<stTR<float>>> rvhGetSymmetricMat(const vector<unique_ptr<stTR<float>>>&, vector<string>);
template unique_ptr<stAbcdnm<float>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<float>>>&,
    const unique_ptr<stIncPar<float>>&, unique_ptr<stIncEabnm<float>>);
template unique_ptr<stCrossSection<float>> rvhGetAverageCrossSections(
    const ArrayXf&, const vector<vector<unique_ptr<stTR<float>>>>&);

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
    const ArrayXf&, const vector<vector<unique_ptr<stTR<float>>>>&);
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
    int, ArrayXf, ArrayXf, float, unique_ptr<stIncPar<float>>, float, float);
template unique_ptr<stRes<float>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<float>>&, const unique_ptr<stParams<float>>&);
template unique_ptr<stSM<float>> pstScatteringMatrixOA(const vector<unique_ptr<stTR<float>>>&,
    float, float, int);

extern template double mp_pi<double>();
extern template complex<double> mp_im_unit();

extern template Tensor4c<double> TensorConj(Tensor4c<double>&);

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

template float Frac2Float(string);
template unique_ptr<RamanParams<float>> LoadParams<float>(string);
template unique_ptr<stParams<float>> Raman2SmartiesParams(const unique_ptr<RamanParams<float>>&, int, int, string);
template void CreateTimeStamp<float, float>(string, const unique_ptr<RamanParams<float>>&);
template void CreateTimeStamp<float, double>(string, const unique_ptr<RamanParams<float>>&);
template void CreateTimeStamp<float, long double>(string, const unique_ptr<RamanParams<float>>&);
template vector<unique_ptr<stTR<float>>> ConvertStTRList(const vector<unique_ptr<stTR<float>>>&);
template vector<unique_ptr<stTR<double>>> ConvertStTRList(const vector<unique_ptr<stTR<float>>>&);
template vector<unique_ptr<stTR<long double>>> ConvertStTRList(const vector<unique_ptr<stTR<float>>>&);
template void RamanElasticScattering<float, float>(string, string);
template void RamanElasticScattering<float, double>(string, string);
template void RamanElasticScattering<float, long double>(string, string);
