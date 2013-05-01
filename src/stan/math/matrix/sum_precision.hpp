#ifndef __STAN__MATH__MATRIX__SUM_PRECISION_HPP__
#define __STAN__MATH__MATRIX__SUM_PRECISION_HPP__

#include <vector>
#include <stan/math/matrix/Eigen.hpp>

namespace stan {
  namespace math {
   
    /**
     * Return the Kahan summation of the values in the specified
     * standard vector.
     *
     * @param xs Standard vector to sum.
     * @return Kahan summation of elements.
     * @tparam T Type of elements summed.
     */
    template <typename T>
    inline T sum_kahan(const std::vector<T>& xs)
    {
      if (xs.size() == 0) return 0;
      T sum(0), t(0), c(0), y(0);
      for (size_t i = 0, size = xs.size(); i < size; ++i)
	{
          y = xs[i] - c;
          t = sum + y;
          c = (t - sum) - y;
          sum = t;
	}
      return sum;
    }

    /**
     * Returns the Kahan summation of the coefficients of the specified
     * matrix.
     * @param xs Specified matrix.
     * @return Sum of coefficients of matrix.
    */
 template <typename T, int R, int C>
    inline T sum_kahan(const Eigen::Matrix<T,R,C>& xs)
    {
      if (xs.size() == 0) return 0;
      T sum(0), t(0), c(0), y(0);
	const T * xsBegin = xs.data();
	for (size_t i = 0, size = xs.size(); i < size; i++)
	{
          y = xsBegin[i] - c;
          t = sum + y;
          c = (t - sum) - y;
          sum = t;
	}
      return sum;
    }







    /**
     * Returns the sum of the values of the specified
     * standard vector using the algorithm described by Ogita, Rump and Oishi (2005).
     *
     * @param xs Standard vector to sum.
     * @return Kahan summation of elements.
     * @tparam T Type of elements summed.
     */
    template <typename T>
    inline T sum_ogita(const std::vector<T>& xs)
    {
      if (xs.size() == 0) return 0;
      T sum(0), x(0), z(0), inter;
      for (size_t i = 0, size = xs.size(); i < size; ++i)
	{
          inter = xs[i];
	    x = sum + inter;
          z = x - inter;
	    sum = x + ( (sum-(x-z)) + (inter-z) );
	}
      return sum;
    }


    /**
     * Returns the sum of the coefficients of the specified
     * matrix using the algorithm described by Ogita, Rump and Oishi (2005): Accurate sum and dot product.
     * @param xs Specified matrix.
     * @return Sum of coefficients of matrix.
    */
 template <typename T, int R, int C>
    inline T sum_ogita(Eigen::Matrix<T,R,C>& xs)
    {
      if (xs.size() == 0) return 0;
      T sum(0), x(0), z(0), inter;
	const T * xsBegin = xs.data();
	for (size_t i = 0, size = xs.size(); i < size; i++)
	{
          inter = xsBegin[i];
	    x = sum + inter;
          z = x - inter;
	    sum = x + ( (sum-(x-z)) + (inter-z) );
	}
      return sum;
    }
 






    /**
     * Returns the pairwise of the values of the specified
     * standard vector.
     *
     * @param xs Standard vector to sum.
     * @return pairwise summation of elements.
     * @tparam T Type of elements summed.
     */
    template <typename T>
    inline T sum_pairwise(const std::vector<T>& xs)
    {
      if (xs.size() == 0) return 0;
      T sum(0);
	size_t i = 0, sizeM1 = xs.size()-1;
      for (; i < sizeM1; ++i)
          sum += (xs[i] + xs[++i]);
	if (i == sizeM1)
          sum += xs[sizeM1]; //needed to add the last number if size is odd
      return sum;
    }


    /**
     * Returns the pairwise sum of the coefficients of the specified
     * matrix.
     * @param xs Specified matrix.
     * @return Sum of coefficients of matrix.
    */
 template <typename T, int R, int C>
    inline T sum_pairwise(Eigen::Matrix<T,R,C>& xs)
    {
      if (xs.size() == 0) return 0;
      T sum(0);
	const T * xsBegin = xs.data();
	size_t i = 0, sizeM1 = xs.size()-1;
      for (; i < sizeM1; ++i)
          sum += (xsBegin[i] + xsBegin[++i]);
	if (i == sizeM1)
          sum += xsBegin[sizeM1]; //needed to add the last number if size is odd
      return sum;
    }
    
  }
}
#endif
