#include "src/raman_elastic_scattering.hpp"
#include "src/smarties.hpp"
#include <boost/multiprecision/mpfr.hpp>

using namespace boost::multiprecision;
using namespace Eigen;
using namespace Smarties;

int precision = 34;

template mpfr_float mp_pi<mpfr_float>();
template mpfr_float mp_eps<mpfr_float>();
template complex<mpfr_float> mp_im_unit<mpfr_float>();

template Tensor3c<mpfr_float> subtensor(Tensor3c<mpfr_float>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
template ArrayXXc<mpfr_float> reduceAndSlice(Tensor3c<mpfr_float>&, int, int);
template ArrayXXc<mpfr_float> invertLUcol(MatrixXc<mpfr_float>&);
// template ArrayXi logicalSlice(ArrayXi&, ArrayXb&);
template RowArrayXr<mpfr_float> logicalSlice(RowArrayXr<mpfr_float>&, RowArrayXb&);
template ArrayXr<mpfr_float> logicalSlice(ArrayXr<mpfr_float>&, ArrayXb&);
template Tensor4c<mpfr_float> tensor_conj(Tensor4c<mpfr_float>&);
template ArrayXr<mpfr_float> arr_bessel_j(ArrayXr<mpfr_float>&, mpfr_float);
template ArrayXr<mpfr_float> arr_bessel_y(ArrayXr<mpfr_float>&, mpfr_float);

template unique_ptr<stGLQuad<mpfr_float>> auxInitLegendreQuad(int, mpfr_float, mpfr_float);
template unique_ptr<stRtfunc<mpfr_float>> auxPrepareIntegrals(int, sInt);

template unique_ptr<stIncPar<mpfr_float>> vshMakeIncidentParams(sIncType, int);
template unique_ptr<stIncPar<mpfr_float>> vshMakeIncidentParams(sIncType, int, mpfr_float, mpfr_float, mpfr_float);
template unique_ptr<stPinmTaunm<mpfr_float>> vshPinmTaunm(int, const ArrayXr<mpfr_float>&);
template unique_ptr<stIncEabnm<mpfr_float>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<mpfr_float>>&);
template unique_ptr<stZnAll<mpfr_float>> vshGetZnAll(int, const ArrayXr<mpfr_float>&, sBessel);
template unique_ptr<stEAllPhi<mpfr_float>> vshEgenThetaAllPhi(const ArrayXr<mpfr_float>&,
    const ArrayXr<mpfr_float>&, const ArrayXXc<mpfr_float>&, const ArrayXXc<mpfr_float>&,
    const RowArrayXr<mpfr_float>&, const RowArrayXr<mpfr_float>&, sBessel, unique_ptr<stPinmTaunm<mpfr_float>>);
template unique_ptr<stEforPhi<mpfr_float>> vshEthetaForPhi(const unique_ptr<stEAllPhi<mpfr_float>>&, mpfr_float);
template ArrayXXc<mpfr_float> vshRBchi(ArrayXr<mpfr_float>, const ArrayXr<mpfr_float>&);
template ArrayXXc<mpfr_float> vshRBpsi(ArrayXr<mpfr_float>, const ArrayXr<mpfr_float>&);

template unique_ptr<ArrayXXr<mpfr_float>> sphGetUforFp(int);
template unique_ptr<stFprow<mpfr_float>> sphGetFpRow(int, mpfr_float, const ArrayXr<mpfr_float>&);
template unique_ptr<stFpovx<mpfr_float>> sphGetFpovx(int, mpfr_float, const ArrayXr<mpfr_float>&);
template unique_ptr<stBessel<mpfr_float>> sphGetXiPsi(int, mpfr_float, const ArrayXr<mpfr_float>&, int);
template unique_ptr<stBesselPrimes<mpfr_float>> sphGetBesselProductsPrimes(const Tensor3c<mpfr_float>&);
template unique_ptr<stBesselProducts<mpfr_float>> sphGetModifiedBesselProducts(int, mpfr_float, const ArrayXr<mpfr_float>&, int);
template unique_ptr<stRtfunc<mpfr_float>> sphMakeGeometry(int, mpfr_float, mpfr_float, const unique_ptr<ArrayXr<mpfr_float>>);
template int sphCheckBesselConvergence(int, mpfr_float, const ArrayXr<mpfr_float>&, mpfr_float, int);
template int sphEstimateNB(int, const unique_ptr<stRtfunc<mpfr_float>>&,
    const unique_ptr<stParams<mpfr_float>>&, mpfr_float);
template vector<unique_ptr<stPQ<mpfr_float>>> sphCalculatePQ(int, const ArrayXi&,
    const unique_ptr<stRtfunc<mpfr_float>>&, const unique_ptr<stParams<mpfr_float>>&, int);

template vector<unique_ptr<stTR<mpfr_float>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<mpfr_float>>>&, bool);
template vector<unique_ptr<stTR<mpfr_float>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<mpfr_float>>>&, int);
template vector<unique_ptr<stTR<mpfr_float>>> rvhGetSymmetricMat(
    const vector<unique_ptr<stTR<mpfr_float>>>&, vector<string>);
template unique_ptr<stAbcdnm<mpfr_float>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<mpfr_float>>>&,
    const unique_ptr<stIncPar<mpfr_float>>&, unique_ptr<stIncEabnm<mpfr_float>>);
template unique_ptr<stCrossSection<mpfr_float>> rvhGetAverageCrossSections(
    const ArrayXr<mpfr_float>&, const vector<vector<unique_ptr<stTR<mpfr_float>>>>&);

template unique_ptr<stTmatrix<mpfr_float>> slvForT(const unique_ptr<stParams<mpfr_float>>&,
    const unique_ptr<stOptions<mpfr_float>>&, unique_ptr<stRtfunc<mpfr_float>>);

template unique_ptr<stRes<mpfr_float>> pstMakeStructForField(const unique_ptr<stAbcdnm<mpfr_float>>&,
    int, ArrayXr<mpfr_float>, ArrayXr<mpfr_float>, mpfr_float, unique_ptr<stIncPar<mpfr_float>>, mpfr_float, mpfr_float);
template unique_ptr<stRes<mpfr_float>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<mpfr_float>>&, const unique_ptr<stParams<mpfr_float>>&);
template unique_ptr<stSM<mpfr_float>> pstScatteringMatrixOA(const vector<unique_ptr<stTR<mpfr_float>>>&,
    mpfr_float, mpfr_float, int);

template unique_ptr<stParams<mpfr_float>> loadParam(string);
template void RamanElasticScattering<mpfr_float>(int, char**);

int main(int argc, char** argv) {
  string calc_type = (argc > 7) ? argv[7] : "double";
  RamanElasticScattering<mpfr_float>(argc, argv);
  return 0;
}
