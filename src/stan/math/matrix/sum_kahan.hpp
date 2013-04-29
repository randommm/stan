#ifndef __STAN__MATH__MATRIX__SUM_KAHAN_HPP__
#define __STAN__MATH__MATRIX__SUM_KAHAN_HPP__

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
    inline T sum_kahan(const std::vector<T>& xs) {
      if (xs.size() == 0) return 0;
      T sumP(0);
      T sumN(0);
      T tP(0);
      T tN(0);
      T cP(0);
      T cN(0);
      T yP(0);
      T yN(0);
      for (size_t i = 0, size = xs.size(); i < size; ++i)
	{
        if (xs[i]>0)
        {
          yP = xs[i] - cP;
	    tP = sumP + yP;
	    cP = (tP - sumP) - yP;
	    sumP = tP;
        }
	  else if (xs[i]<0)
        {
          yN = xs[i] - cN;
	    tN = sumN + yN;
	    cN = (tN - sumN) - yN;
	    sumN = tN;
        }
	}
      return sumP+sumN;
    }

    /**
     * Returns the Kahan summation of the coefficients of the specified
     * matrix.
     * @param xs Specified matrix.
     * @return Kahan summation of coefficients of vector.
    */
 template <typename T, int R, int C>
    inline T sum_kahan(const Eigen::Matrix<T,R,C>& xs) {
      if (xs.size() == 0) return 0;
      T sumP(0), sumN(0), tP(0), tN(0), cP(0), cN(0), yP(0), yN(0);
	for (size_t i = 0, size = xs.size(); i < size; i++)
	  {
	  T temporary = (*(xs.data() + i));
          if (temporary>0)
          {
            yP = temporary - cP;
	      tP = sumP + yP;
	      cP = (tP - sumP) - yP;
	      sumP = tP;
          }
	    else if (temporary<0)
          {
            yN = temporary - cN;
	      tN = sumN + yN;
	      cN = (tN - sumN) - yN;
	      sumN = tN;
          }
	  }
      return sumP+sumN;
    }
    
  }
}
#endif
