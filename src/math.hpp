#ifndef SMARTIES_MATH_HPP
#define SMARTIES_MATH_HPP

#include "core.hpp"
#include <Eigen/LU>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/constants/constants.hpp>

using namespace Eigen;
using namespace std;

namespace Smarties {
  template <class Real>
  inline Real mp_pi(void) { return boost::math::constants::pi<Real>(); }

  template <class Real>
  inline Real mp_eps(void) { return numeric_limits<Real>::epsilon(); }

  template <class Real>
  inline complex<Real> mp_im_unit(void) { return complex<Real>(0.0, 1.0); }

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
    int rows = sequence.size();
    ArrayXi output(rows);
    for (int i = 0; i < rows; i++)
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
  ArrayXXc<Real> invertLUcol(MatrixXc<Real>& B) {
    PartialPivLU<MatrixXc<Real>> PLU_decomp = B.matrix().lu();
    MatrixXc<Real> P = PLU_decomp.permutationP();
    MatrixXc<Real> L = PLU_decomp.matrixLU().template triangularView<UnitLower>();
    MatrixXc<Real> U = PLU_decomp.matrixLU().template triangularView<Upper>();
    MatrixXc<Real> Y = L.lu().solve(P);
    return U.lu().solve(Y).array();
  }

  // Expects base_array.dimensions() = bool_array.dimensions()
  template <class Real>
  ArrayXr<Real> logicalSlice(ArrayXr<Real>& base_array, ArrayXb& bool_array) {
    int output_size = bool_array.count();
    ArrayXr<Real> output(output_size);
    int base_index = 0;
    for (int i = 0; i < output_size && base_index < bool_array.size(); i++, base_index++) {
      while (!bool_array(base_index))
        base_index++;
      output(i) = base_array(base_index);
    }
    return output;
  }

  // Expects base_array.dimensions() = bool_array.dimensions()
  template <class Real>
  RowArrayXr<Real> logicalSlice(RowArrayXr<Real>& base_array, RowArrayXb& bool_array) {
    int output_size = bool_array.count();
    RowArrayXr<Real> output(output_size);
    int base_index = 0;
    for (int i = 0; i < output_size && base_index < bool_array.size(); i++, base_index++) {
      while (!bool_array(base_index))
        base_index++;
      output(i) = base_array(base_index);
    }
    return output;
  }

  ArrayXi logicalIndices(ArrayXb& bool_array) {
    int output_size = bool_array.count();
    ArrayXi output(output_size);
    int index = 0;
    for (int i = 0; i < output_size && index < bool_array.size(); i++, index++) {
      while (!bool_array(index))
        index++;
      output(i) = index;
    }
    return output;
  }

  template <class Real>
  Tensor4c<Real> tensor_conj(Tensor4c<Real>& base) {
    Tensor4c<Real> output = base;
    for (int i = 0; i < base.dimension(0); i++) {
      for (int j = 0; j < base.dimension(1); j++) {
        for (int k = 0; k < base.dimension(2); k++) {
          for (int l = 0; l < base.dimension(3); l++)
            output(i, j, k, l) = conj(base(i, j, k, l));
        }
      }
    }
    return output;
  }

  template <class Real>
  ArrayXr<Real> arr_bessel_j(ArrayXr<Real>& nu, Real x) {
    ArrayXr<Real> output(nu.size());
    for (int i = 0; i < nu.size(); i++)
      output(i) = boost::math::cyl_bessel_j<Real>(nu(i), x);
    return output;
  }

  template <class Real>
  ArrayXr<Real> arr_bessel_y(ArrayXr<Real>& nu, Real x) {
    ArrayXr<Real> output(nu.size());
    for (int i = 0; i < nu.size(); i++)
      output(i) = boost::math::cyl_neumann<Real>(nu(i), x);
    return output;
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
  template RowArrayXr<double> logicalSlice(RowArrayXr<double>&, RowArrayXb&);
  template ArrayXd logicalSlice(ArrayXd&, ArrayXb&);
  template Tensor4c<double> tensor_conj(Tensor4c<double>&);
  template ArrayXd arr_bessel_j(ArrayXd&, double);
  template ArrayXd arr_bessel_y(ArrayXd&, double);
}

#endif
