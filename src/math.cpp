#include "math.h"
#include "raman_elastic_scattering.h"
#include <Eigen/LU>

using namespace Eigen;
using namespace std;

namespace Raman {
  template <class Real>
  Real mp_pi(void) {
    return 3.14159265359;
  }

  template <class Real>
  Real mp_eps(void) {
    return numeric_limits<Real>::epsilon();
  }

  template <class Real>
  complex<Real> mp_im_unit(void) {
    complex<Real> i(0.0, 1.0);
    return i;
  }

  template <class Real>
  Tensor3c<Real> subtensor(Tensor3c<Real>& tensor,
      ArithmeticSequence<long int, long int, long int> slice_dim1,
      ArithmeticSequence<long int, long int, long int> slice_dim2,
      ArithmeticSequence<long int, long int, long int> slice_dim3) {
    long int new_dim1 = slice_dim1.size(), new_dim2 = slice_dim2.size(), new_dim3 = slice_dim3.size();
    Tensor3c<Real> output(new_dim1, new_dim2, new_dim3);
    for (long int i = 0; i < new_dim1; i++) {
      for (long int j = 0; j < new_dim2; j++) {
        for (long int k = 0; k < new_dim3; k++)
          output(i, j, k) = tensor(slice_dim1[i], slice_dim2[j], slice_dim3[k]);
      }
    }
    return output;
  }

  ArrayXi seq2Array(long int first, long int last, long int stride) {
    ArithmeticSequence<long int, long int, long int> sequence = seq(first, last, stride);
    size_t rows = sequence.size();
    ArrayXi output(rows);
    for (size_t i = 0; i < rows; i++)
      output(i) = sequence[i];
    return output;
  }

  template <class Real>
  ArrayXXc<Real> reduceAndSlice(Tensor3c<Real>& tensor, int offset, int num_rows) {
    const auto dims = tensor.dimensions();
    ArrayXXc<Real> output(num_rows, dims[1]);
    for (int i = 0; i < num_rows; i++) {
      for (int j = 0; j < dims[1]; j++)
        output(i, j) = tensor(dims[0] - num_rows + i, j, offset);
    }
    return output;
  }

  template <class Real>
  ArrayXXc<Real> invertLUcol(ArrayXXc<Real>& B) {
    PartialPivLU<MatrixXc<Real>> PLU_decomp = B.matrix().lu();
    MatrixXc<Real> P = PLU_decomp.permutationP();
    MatrixXc<Real> L = PLU_decomp.matrixLU().template triangularView<UnitLower>();
    MatrixXc<Real> U = PLU_decomp.matrixLU().template triangularView<Upper>();
    MatrixXc<Real> Y = L.lu().solve(P);
    return U.lu().solve(Y).array();
  }

  template <class Real>
  ArrayXr<Real> arr_bessel_j(ArrayXr<Real>& nu, Real x) {
    ArrayXr<Real> output(nu.size());
    for (int i = 0; i < nu.size(); i++)
      output(i) = cyl_bessel_j(nu(i), x);
    return output;
  }

  template <class Real>
  ArrayXr<Real> arr_bessel_y(ArrayXr<Real>& nu, Real x) {
    ArrayXr<Real> output(nu.size());
    for (int i = 0; i < nu.size(); i++)
      output(i) = cyl_neumann(nu(i), x);
    return output;
  }

  template Tensor3c<double> subtensor(Tensor3c<double>&,
      ArithmeticSequence<long int, long int, long int>,
      ArithmeticSequence<long int, long int, long int>,
      ArithmeticSequence<long int, long int, long int>);
  template ArrayXXc<double> reduceAndSlice(Tensor3c<double>&, int, int);
  template ArrayXXc<double> invertLUcol(ArrayXXc<double>&);
  template ArrayXr<double> arr_bessel_j(ArrayXr<double>&, double);
  template ArrayXr<double> arr_bessel_y(ArrayXr<double>&, double);
}
