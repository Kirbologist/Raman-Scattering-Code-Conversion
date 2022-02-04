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
