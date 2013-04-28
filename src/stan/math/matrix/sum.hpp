#ifndef __STAN__MATH__MATRIX__SUM_HPP__
#define __STAN__MATH__MATRIX__SUM_HPP__

#include <vector>
#include <stan/math/matrix/Eigen.hpp>

namespace stan {
  namespace math {
   
    /**
     * Return the sum of the values in the specified
     * standard vector.
     *
     * @param xs Standard vector to sum.
     * @return Sum of elements.
     * @tparam T Type of elements summed.
     */
    template <typename T>
    inline T sum(const std::vector<T>& xs) {
      if (xs.size() == 0) return 0;
      T sum(xs[0]);
      for (size_t i = 1, size=xs.size(); i < size; ++i)
        sum += xs[i];
      return sum;
    }
 
    /**
     * Returns the sum of the coefficients of the specified
     * matrix.
     * @param v Specified matrix.
     * @return Sum of coefficients of matrix.
     */
    template <typename T, int R, int C>
    inline double sum(const Eigen::Matrix<T,R,C>& v) {
      return v.sum();
    }    
    
  }
}
#endif
