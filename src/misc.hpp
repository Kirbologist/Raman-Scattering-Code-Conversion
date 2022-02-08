/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Based on the SMARTIES MATLAB package by W.R.C. Somerville, B. Auguié, E.C. Le Ru
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.


Various small maths functions, helper functions and others used throughout SMARTIES functions that don't fit anywhere.
*/

#ifndef SMARTIES_MATH_HPP
#define SMARTIES_MATH_HPP

#include "core.hpp"
#include <Eigen/LU>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/constants/constants.hpp>

using namespace Eigen;
using namespace std;

namespace Smarties {
  /* Returns the mathematical constant pi with the amount of precision specified by template argument */
  template <class Real>
  inline Real mp_pi(void) { return boost::math::constants::pi<Real>(); }

  /* Returns the machine epsilon of the type given in the template argument */
  template <class Real>
  inline Real mp_eps(void) { return numeric_limits<Real>::epsilon(); }

  /* Returns the imaginary unit 'i' with the type given in the template argument */
  template <class Real>
  inline complex<Real> mp_im_unit(void) { return complex<Real>(0.0, 1.0); }

  /*
  View an existing Eigen::Tensor of rank 2 as an Eigen::Map<Eigen::Matrix>
  Rows/Cols are determined from the matrix
  */
  template<typename Scalar>
  auto ArrayMap(Tensor<Scalar, 2>& tensor) {
      return Eigen::Map<ArrayXXr<Scalar>>(tensor.data(), tensor.dimension(0), tensor.dimension(1));
  }

  /*
  Converts an Eigen::Matrix (or expression) to Eigen::Tensor
  with dimensions specified in std::array
  */
  template<typename Derived, typename T, auto rank>
  Tensor<typename Derived::Scalar, rank>
  TensorCast(const EigenBase<Derived>& matrix, const std::array<T, rank>& dims) {
      return Eigen::TensorMap<const Tensor<const typename Derived::Scalar, rank>>
                  (matrix.derived().eval().data(), dims);
  }

  /*
  Converts an Eigen::Matrix (or expression) to Eigen::Tensor
  with dimensions as variadic arguments
  */
  template<typename Derived, typename... Dims>
  auto TensorCast(const EigenBase<Derived>& matrix, const Dims... dims) {
      static_assert(sizeof...(Dims) > 0, "TensorCast: sizeof... (Dims) must be larger than 0");
      return TensorCast(matrix, std::array<Eigen::Index, sizeof...(Dims)>{dims...});
  }

  /*
  Converts an Eigen::Matrix (or expression) to Eigen::Tensor
  with dimensions directly as arguments in a variadic template
  */
  template<typename Derived>
  auto TensorCast(const EigenBase<Derived>& matrix) {
    if constexpr(Derived::ColsAtCompileTime == 1 or Derived::RowsAtCompileTime == 1) {
      return TensorCast(matrix, matrix.size());
    } else {
      return TensorCast(matrix, matrix.rows(), matrix.cols());
    }
  }

  /* Get a readable name/description of the type given in the template argument */
  template <class Real>
  string GetTypeName() {
    string calc_type = typeid(static_cast<Real>(0)).name();
    if (calc_type.find("N5boost14multiprecision6numberINS0_8backends18mpfr_float_backend") != string::npos) {
      size_t prec_begin = calc_type.find("ILj") + 3;
      size_t prec_end = calc_type.find("ELNS");
      string precision = calc_type.substr(prec_begin, prec_end - prec_begin);
      calc_type = "custom_" + precision;
    } else if (calc_type == "f")
      calc_type = "single";
    else if (calc_type == "d")
      calc_type = "double";
    else if (calc_type == "e")
      calc_type = "quad";
    return calc_type;
  }

  /*
  Get a slice of a rank-3 tensor, so that the indices of the slices are given by Eigen::ArithmeticSequence.
  This is just a helper function for sphCheckBesselConvergence, and is not usable for general purpose.
  */
  template <class Real>
  Tensor3c<Real> TensorSlice(Tensor3c<Real>& tensor,
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

  /*
  Initialize a linearly-spaced single-column Eigen::Array by converting from an Eigen::ArithmeticSequence.
  Parameters are the same ones used to construct an Eigen::ArithmeticSequence.
  */
  ArrayXi Seq2Array(long int first, long int last, long int stride);

  /*
  Reduce a complex rank-3 Eigen::Tensor into a complex Eigen::Array by
  reducing the third dimension of the tensor down to only coefficients whose third index is `offset`,
  and reducing the first dimenstion down to only the last `num_rows` rows.
  This was used specifically as a helper function in sphCalculatePQ.
  */
  template <class Real>
  ArrayXXc<Real> Subtensor2ArrMap(const Tensor3c<Real>& tensor,
      const std::array<int, 3>& offsets, const std::array<int, 3>& extents, int rows, int cols) {
    ArrayXXc<Real> output(rows, cols);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++)
        output(i, j) = tensor(offsets[0] + i, offsets[1], offsets[2] + j);
    }
    return output;
  }

  /* Performs matrix inversion using LU with partial pivoting of columns. This is used in rvhGetTRfromPQ. */
  template <class Real>
  ArrayXXc<Real> InvertLUcol(MatrixXc<Real>& M) {
    PartialPivLU<MatrixXc<Real>> PLU_decomp = M.matrix().lu();
    return PLU_decomp.inverse().array();
  }

  /* Performs the equivalent of MATLAB's logical indexing on a single-column Eigen::Array. */
  template <class Real>
  ArrayXr<Real> LogicalSlice(ArrayXr<Real>& base_array, ArrayXb& bool_array) {
    assert(base_array.size() == bool_array.size());
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

  /* Performs the equivalent of MATLAB's logical indexing on a single-row Eigen::Array. */
  template <class Real>
  RowArrayXr<Real> LogicalSlice(RowArrayXr<Real>& base_array, RowArrayXb& bool_array) {
    assert(base_array.size() == bool_array.size());
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

  /* Gets the indices of all 'true' components of a single-column boolean Eigen::Array. */
  ArrayXi LogicalIndices(ArrayXb& bool_array);

  /* Take a tensor and return the same tensor, but with all its coefficients conjugated. */
  template <class Real>
  Tensor4c<Real> TensorConj(Tensor4c<Real>& base) {
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

  /* Calculates the Bessel function of the first kind with many values of `nu` and a fixed value of `x`. */
  template <class Real>
  ArrayXr<Real> ArrBesselJ(ArrayXr<Real>& nu, Real x) {
    ArrayXr<Real> output(nu.size());
    for (int i = 0; i < nu.size(); i++)
      output(i) = boost::math::cyl_bessel_j<Real>(nu(i), x);
    return output;
  }

  /* Calculates the Bessel function of the second kind with many values of `nu` and a fixed value of `x`. */
  template <class Real>
  ArrayXr<Real> ArrBesselY(ArrayXr<Real>& nu, Real x) {
    ArrayXr<Real> output(nu.size());
    for (int i = 0; i < nu.size(); i++)
      output(i) = boost::math::cyl_neumann<Real>(nu(i), x);
    return output;
  }
}

#endif
