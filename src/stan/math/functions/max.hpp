#ifndef __STAN__MATH__FUNCTIONS__MAX_HPP__
#define __STAN__MATH__FUNCTIONS__MAX_HPP__

#include <boost/math/tools/promotion.hpp>

namespace stan {
  namespace math {

    template <typename T1, typename T2>
    inline typename boost::math::tools::promote_args<T1,T2>::type
    max(T1 a, T2 b) { 
      return a > b ? a : b; 
    }

    inline int max(int a, int b) { 
      return a > b ? a : b; 
    }

  }
}

#endif
