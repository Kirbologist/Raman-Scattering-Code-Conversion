#include "smarties.hpp"
#include "raman_elastic_scattering.hpp"
#include <fstream>

using namespace Eigen;
using namespace Smarties;

string GetCalcType(string in_file_name) {
  ifstream in_file;
  in_file.open(in_file_name, ios::in);
  if (!in_file.is_open())
    throw runtime_error("Error: cannot open config.txt");

  string option = "Calculation type: ";
  string line;
  do {
    if (in_file.peek() == EOF)
      throw runtime_error("Error: cannot find option " + option);
    getline(in_file, line);
  } while (line.find(option) == string::npos);
  string calc_type = line.substr(line.find(option) + option.size());
  in_file.close();
  return calc_type;
}

void MultiPrint(string out_string, string out_file_name, bool write_output) {
  cout << out_string;
  cout.flush();
  if (write_output) {
    ofstream out_file;
    out_file.open(out_file_name, ios::out | ios::app);
    if (!out_file.is_open()) {
      cerr << "Warning: cannot open output file. Some output won't be written." << endl;
      out_file.close();
      write_output = false;
    }
    out_file << out_string;
    out_file.flush();
    out_file.close();
  }
}

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

/*
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
*/

template unique_ptr<stParams<double>> loadParam(string);
template unique_ptr<RamanParams<double>> GetRamanParams(string);
template void CreateTimeStamp(string, const unique_ptr<RamanParams<double>>&, bool);
template void RamanElasticScattering<double>(string, string);
