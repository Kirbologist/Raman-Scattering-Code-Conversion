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

  // View an existing Eigen::Tensor of rank 2 as an Eigen::Map<Eigen::Matrix>
  // Rows/Cols are determined from the matrix
  template<typename Scalar>
  auto ArrayMap(Eigen::Tensor<Scalar, 2> &tensor) {
      return Eigen::Map<ArrayXXr<Scalar>>(tensor.data(), tensor.dimension(0), tensor.dimension(1));
  }

  // Converts an Eigen::Matrix (or expression) to Eigen::Tensor
  // with dimensions specified in std::array
  template<typename Derived, typename T, auto rank>
  Eigen::Tensor<typename Derived::Scalar, rank>
  TensorCast(const Eigen::EigenBase<Derived> &matrix, const std::array<T, rank> &dims) {
      return Eigen::TensorMap<const Eigen::Tensor<const typename Derived::Scalar, rank>>
                  (matrix.derived().eval().data(), dims);
  }

  // Converts an Eigen::Matrix (or expression) to Eigen::Tensor
  // with dimensions as variadic arguments
  template<typename Derived, typename... Dims>
  auto TensorCast(const Eigen::EigenBase<Derived> &matrix, const Dims... dims) {
      static_assert(sizeof...(Dims) > 0, "TensorCast: sizeof... (Dims) must be larger than 0");
      return TensorCast(matrix, std::array<Eigen::Index, sizeof...(Dims)>{dims...});
  }

  // Converts an Eigen::Matrix (or expression) to Eigen::Tensor
  // with dimensions directly as arguments in a variadic template
  template<typename Derived>
  auto TensorCast(const Eigen::EigenBase<Derived> &matrix) {
    if constexpr(Derived::ColsAtCompileTime == 1 or Derived::RowsAtCompileTime == 1) {
      return TensorCast(matrix, matrix.size());
    } else {
      return TensorCast(matrix, matrix.rows(), matrix.cols());
    }
  }

  template <class Real>
  Tensor3c<Real> tensorSlice(Tensor3c<Real>& tensor,
      ArithmeticSequence<long int, long int, long int> slice_dim1,
      ArithmeticSequence<long int, long int, long int> slice_dim2,
      ArithmeticSequence<long int, long int, long int> slice_dim3) {
    long int new_dim1 = slice_dim1.size();
    long int new_dim2 = slice_dim2.size();
    long int new_dim3 = slice_dim3.size();
    Tensor3c<Real> output(new_dim1, new_dim2, new_dim3);
    for (long int i = 0; i < new_dim1; i++) {
      for (long int j = 0; j < new_dim2; j++) {
        for (long int k = 0; k < new_dim3; k++)
          output(i, j, k) = tensor(slice_dim1[i], slice_dim2[j], slice_dim3[k]);
      }
    }
    return output;
  }

  ArrayXi seq2Array(long int first, long int last, long int stride);

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

  ArrayXi logicalIndices(ArrayXb& bool_array);

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
}

#endif
