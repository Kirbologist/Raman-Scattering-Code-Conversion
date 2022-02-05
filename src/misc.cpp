/*
This file is a part of Raman-Scattering-Code-Conversion.
<https://github.com/Kirbologist/Raman-Scattering-Code-Conversion>

Written by Siwan Li for the UQ School of Maths and Physics.
Based on the SMARTIES MATLAB package by W.R.C. Somerville, B. Auguié, E.C. Le Ru
Copyright (C) 2021-2022 Siwan Li

This source code form is subject to the terms of the MIT License.
If a copy of the MIT License was not distributed with this file,
you can obtain one at <https://opensource.org/licenses/MIT>.
*/

#include "misc.hpp"

using namespace Eigen;

namespace Smarties {

  ArrayXi Seq2Array(long int first, long int last, long int stride) {
    ArithmeticSequence<long int, long int, long int> sequence = seq(first, last, stride);
    int rows = sequence.size();
    ArrayXi output(rows);
    for (int i = 0; i < rows; i++)
      output(i) = sequence[i];
    return output;
  }

  ArrayXi LogicalIndices(ArrayXb& bool_array) {
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
}
