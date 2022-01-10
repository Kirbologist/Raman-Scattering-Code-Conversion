#ifndef SMARTIES_MATH_H
#define SMARTIES_MATH_H

#include "core.h"

using namespace std;
using namespace Eigen;

namespace Raman {
  template <class Real>
  Real mp_pi(void);

  template <class Real>
  Real mp_eps(void);

  template <class Real>
  complex<Real> mp_im_unit(void);

  template <class Real>
  Tensor3c<Real> subtensor(Tensor3c<Real>& tensor,
      ArithmeticSequence<long int, long int, long int> slice_dim1,
      ArithmeticSequence<long int, long int, long int> slice_dim2,
      ArithmeticSequence<long int, long int, long int> slice_dim3);

  ArrayXi seq2Array(long int first, long int last, long int stride);

  template <class Real>
  ArrayXXc<Real> reduceAndSlice(Tensor3c<Real>& tensor, int offset, int num_rows);

  template <class Real>
  ArrayXXc<Real> invertLUcol(MatrixXc<Real>& B);

  template <class Real>
  ArrayXr<Real> logicalSlice(ArrayXr<Real>& base_array, ArrayXb& bool_array);

  template <class Real>
  RowArrayXr<Real> logicalSlice(RowArrayXr<Real>& base_array, RowArrayXb& bool_array);

  ArrayXi logicalIndices(ArrayXb& bool_array);

  template <class Real>
  Tensor4c<Real> tensor_conj(Tensor4c<Real>& base);

  template <class Real>
  ArrayXr<Real> arr_bessel_j(ArrayXr<Real>& nu, Real x);

  template <class Real>
  ArrayXr<Real> arr_bessel_y(ArrayXr<Real>& nu, Real x);
}

#endif
