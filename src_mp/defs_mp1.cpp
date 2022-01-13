#include "../src/math.hpp"
#include "../src/smarties_aux.hpp"
#include "../src/vsh.hpp"
#include "../src/core.hpp"
#include "core_mp.hpp"

using namespace Eigen;
using namespace Smarties;
using namespace boost::multiprecision;

template raman_float mp_pi();
template raman_float mp_eps();
template complex<raman_float> mp_im_unit();

template Tensor3c<raman_float> subtensor(Tensor3c<raman_float>&,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>,
    ArithmeticSequence<long int, long int, long int>);
template ArrayXXc<raman_float> reduceAndSlice(Tensor3c<raman_float>&, int, int);
template ArrayXXc<raman_float> invertLUcol(MatrixXc<raman_float>&);
// template ArrayXi logicalSlice(ArrayXi&, ArrayXb&);
template RowArrayXr<raman_float> logicalSlice(RowArrayXr<raman_float>&, RowArrayXb&);
template ArrayXr<raman_float> logicalSlice(ArrayXr<raman_float>&, ArrayXb&);
template Tensor4c<raman_float> tensor_conj(Tensor4c<raman_float>&);
template ArrayXr<raman_float> arr_bessel_j(ArrayXr<raman_float>&, raman_float);
template ArrayXr<raman_float> arr_bessel_y(ArrayXr<raman_float>&, raman_float);

template unique_ptr<stGLQuad<raman_float>> auxInitLegendreQuad(int, raman_float, raman_float);
template unique_ptr<stRtfunc<raman_float>> auxPrepareIntegrals(int, sInt);

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
