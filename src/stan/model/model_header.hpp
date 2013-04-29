#ifndef STAN__MODEL__MODEL__HEADER_HPP__
#define STAN__MODEL__MODEL__HEADER_HPP__

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/exception/all.hpp>

#include <stan/agrad/agrad.hpp>
#include <stan/agrad/matrix.hpp>
#include <stan/agrad/rev/error_handling/matrix/check_pos_definite.hpp>
#include <stan/agrad/partials_vari.hpp>
#include <stan/gm/command.hpp>
#include <stan/io/cmd_line.hpp>
#include <stan/io/dump.hpp>
#include <stan/io/reader.hpp>
#include <stan/io/writer.hpp>
#include <stan/io/csv_writer.hpp>

#include <stan/math/matrix.hpp>
#include <stan/math.hpp>
// FIXME: these should go in matrix.hpp
#include <stan/math/rep_array.hpp>
#include <stan/math/rep_vector.hpp>
#include <stan/math/rep_row_vector.hpp>
#include <stan/math/rep_matrix.hpp>
#include <stan/math/matrix/add.hpp>
#include <stan/math/matrix/block.hpp>
#include <stan/math/matrix/cholesky_decompose.hpp>
#include <stan/math/matrix/col.hpp>
#include <stan/math/matrix/columns_dot_product.hpp>
#include <stan/math/matrix/columns_dot_self.hpp>
#include <stan/math/matrix/crossprod.hpp>
#include <stan/math/matrix/cumulative_sum.hpp>
#include <stan/math/matrix/determinant.hpp>
#include <stan/math/matrix/diag_matrix.hpp>
#include <stan/math/matrix/diag_post_multiply.hpp>
#include <stan/math/matrix/diag_pre_multiply.hpp>
#include <stan/math/matrix/diagonal.hpp>
#include <stan/math/matrix/dims.hpp>
#include <stan/math/matrix/divide.hpp>
#include <stan/math/matrix/dot_product.hpp>
#include <stan/math/matrix/dot_self.hpp>
#include <stan/math/matrix/eigenvalues_sym.hpp>
#include <stan/math/matrix/eigenvectors_sym.hpp>
#include <stan/math/matrix/elt_divide.hpp>
#include <stan/math/matrix/elt_multiply.hpp>
#include <stan/math/matrix/exp.hpp>
#include <stan/math/matrix/head.hpp>
#include <stan/math/matrix/inverse.hpp>
#include <stan/math/matrix/log.hpp>
#include <stan/math/matrix/log_determinant.hpp>
#include <stan/math/matrix/max.hpp>
#include <stan/math/matrix/mdivide_left.hpp>
#include <stan/math/matrix/mdivide_left_tri.hpp>
#include <stan/math/matrix/mdivide_left_tri_low.hpp>
#include <stan/math/matrix/mdivide_right.hpp>
#include <stan/math/matrix/mdivide_right_tri.hpp>
#include <stan/math/matrix/mdivide_right_tri_low.hpp>
#include <stan/math/matrix/mean.hpp>
#include <stan/math/matrix/min.hpp>
#include <stan/math/matrix/minus.hpp>
#include <stan/math/matrix/multiply.hpp>
#include <stan/math/matrix/multiply_lower_tri_self_transpose.hpp>
#include <stan/math/matrix/prod.hpp>
#include <stan/math/matrix/row.hpp>
#include <stan/math/matrix/rows_dot_product.hpp>
#include <stan/math/matrix/rows_dot_self.hpp>
#include <stan/math/matrix/sd.hpp>
#include <stan/math/matrix/segment.hpp>
#include <stan/math/matrix/singular_values.hpp>
#include <stan/math/matrix/size.hpp>
#include <stan/math/matrix/softmax.hpp>
#include <stan/math/matrix/stan_print.hpp>
#include <stan/math/matrix/sub_col.hpp>
#include <stan/math/matrix/sub_row.hpp>
#include <stan/math/matrix/subtract.hpp>
#include <stan/math/matrix/sum.hpp>
#include <stan/math/matrix/sum_kahan.hpp>
#include <stan/math/matrix/tail.hpp>
#include <stan/math/matrix/tcrossprod.hpp>
#include <stan/math/matrix/trace.hpp>
#include <stan/math/matrix/transpose.hpp>
#include <stan/math/matrix/variance.hpp>

#include <stan/model/prob_grad_ad.hpp>
#include <stan/prob/distributions.hpp>

#endif
