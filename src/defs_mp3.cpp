#include "slv.hpp"
#include "pst.hpp"
#include "math.hpp"
#include "vsh.hpp"
#include "sph.hpp"
#include "rvh.hpp"
#include "core.hpp"
#include "core_mp.hpp"

using namespace Eigen;
using namespace Smarties;
using namespace boost::multiprecision;

extern template raman_float mp_pi();
extern template raman_float mp_eps();
extern template complex<raman_float> mp_im_unit();

extern template unique_ptr<stIncPar<raman_float>> vshMakeIncidentParams(sIncType, int);

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

template unique_ptr<stTmatrix<raman_float>> slvForT(const unique_ptr<stParams<raman_float>>&,
    const unique_ptr<stOptions<raman_float>>&, unique_ptr<stRtfunc<raman_float>>);

template unique_ptr<stRes<raman_float>> pstMakeStructForField(const unique_ptr<stAbcdnm<raman_float>>&,
    int, ArrayXr<raman_float>, ArrayXr<raman_float>, raman_float, unique_ptr<stIncPar<raman_float>>, raman_float, raman_float);
template unique_ptr<stRes<raman_float>> pstMakeStructForField(
    const unique_ptr<stAbcdnm<raman_float>>&, const unique_ptr<stParams<raman_float>>&);
template unique_ptr<stSM<raman_float>> pstScatteringMatrixOA(const vector<unique_ptr<stTR<raman_float>>>&,
    raman_float, raman_float, int);
