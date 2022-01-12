#include "sph.hpp"
#include "rvh.hpp"
#include "math.hpp"
#include "smarties_aux.hpp"
#include "vsh.hpp"
#include "core.hpp"
#include "core_mp.hpp"

using namespace Eigen;
using namespace Smarties;
using namespace boost::multiprecision;

extern template raman_float mp_pi();
extern template raman_float mp_eps();
extern template complex<raman_float> mp_im_unit();

extern template Tensor3c<raman_float> subtensor(Tensor3c<raman_float>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
extern template ArrayXXc<raman_float> reduceAndSlice(Tensor3c<raman_float>&, int, int);
extern template ArrayXXc<raman_float> invertLUcol(MatrixXc<raman_float>&);
//extern template ArrayXi logicalSlice(ArrayXi&, ArrayXb&);
extern template RowArrayXr<raman_float> logicalSlice(RowArrayXr<raman_float>&, RowArrayXb&);

extern template unique_ptr<stRtfunc<raman_float>> auxPrepareIntegrals(int, sInt);

extern template unique_ptr<stIncEabnm<raman_float>> vshGetIncidentCoeffs(int, const unique_ptr<stIncPar<raman_float>>&);
extern template unique_ptr<stPinmTaunm<raman_float>> vshPinmTaunm(int, const ArrayXr<raman_float>&);
extern template ArrayXXc<raman_float> vshRBchi(ArrayXr<raman_float>, const ArrayXr<raman_float>&);
extern template ArrayXXc<raman_float> vshRBpsi(ArrayXr<raman_float>, const ArrayXr<raman_float>&);

template unique_ptr<ArrayXXr<raman_float>> sphGetUforFp(int);
template unique_ptr<stFprow<raman_float>> sphGetFpRow(int, raman_float, const ArrayXr<raman_float>&);
template unique_ptr<stFpovx<raman_float>> sphGetFpovx(int, raman_float, const ArrayXr<raman_float>&);
template unique_ptr<stBessel<raman_float>> sphGetXiPsi(int, raman_float, const ArrayXr<raman_float>&, int);
template unique_ptr<stBesselPrimes<raman_float>> sphGetBesselProductsPrimes(const Tensor3c<raman_float>&);
template unique_ptr<stBesselProducts<raman_float>> sphGetModifiedBesselProducts(int, raman_float, const ArrayXr<raman_float>&, int);
template unique_ptr<stRtfunc<raman_float>> sphMakeGeometry(int, raman_float, raman_float, const unique_ptr<ArrayXr<raman_float>>);
template int sphCheckBesselConvergence(int, raman_float, const ArrayXr<raman_float>&, raman_float, int);
template int sphEstimateNB(int, const unique_ptr<stRtfunc<raman_float>>&,
    const unique_ptr<stParams<raman_float>>&, raman_float);
template vector<unique_ptr<stPQ<raman_float>>> sphCalculatePQ(int, const ArrayXi&,
    const unique_ptr<stRtfunc<raman_float>>&, const unique_ptr<stParams<raman_float>>&, int);

template vector<unique_ptr<stTR<raman_float>>> rvhGetTRfromPQ(vector<unique_ptr<stPQ<raman_float>>>&, bool);
template vector<unique_ptr<stTR<raman_float>>> rvhTruncateMatrices(const vector<unique_ptr<stTR<raman_float>>>&, int);
template vector<unique_ptr<stTR<raman_float>>> rvhGetSymmetricMat(
    const vector<unique_ptr<stTR<raman_float>>>&, vector<string>);
template unique_ptr<stAbcdnm<raman_float>> rvhGetFieldCoefficients(int, const vector<unique_ptr<stTR<raman_float>>>&,
    const unique_ptr<stIncPar<raman_float>>&, unique_ptr<stIncEabnm<raman_float>>);
template unique_ptr<stCrossSection<raman_float>> rvhGetAverageCrossSections(
    const ArrayXr<raman_float>&, const vector<vector<unique_ptr<stTR<raman_float>>>>&);
