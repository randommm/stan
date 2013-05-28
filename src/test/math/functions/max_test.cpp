#include "stan/math/functions/max.hpp"
#include <gtest/gtest.h>
#include <typeinfo>

TEST(MathFunctions, max) {
  using stan::math::max;
  double a = 1.23;
  double b = 1;
  int c = 2;
  int d = 3;
  
  EXPECT_TRUE(typeid(max(a,b)) == typeid(double));
  EXPECT_TRUE(typeid(max(a,c)) == typeid(double));
  EXPECT_TRUE(typeid(max(a,2.4)) == typeid(a));
  EXPECT_TRUE(typeid(max(d,b)) == typeid(double));
  EXPECT_TRUE(typeid(max(d,a)) == typeid(b));
  EXPECT_TRUE(typeid(max(c,d)) == typeid(c));

  EXPECT_FALSE(typeid(max(a,b)) == typeid(int));
  EXPECT_FALSE(typeid(max(a,c)) == typeid(c));
  EXPECT_FALSE(typeid(max(a,2.4)) == typeid(d));
  EXPECT_FALSE(typeid(max(d,b)) == typeid(int));
  EXPECT_FALSE(typeid(max(d,a)) == typeid(d));
  EXPECT_FALSE(typeid(max(c,d)) == typeid(b));
  
  EXPECT_EQ(max(a,b), a);
  EXPECT_EQ(max(a,c), c);
  EXPECT_EQ(max(a,d), d);
  EXPECT_EQ(max(a,c), c);
  EXPECT_EQ(max(d,b), d);
  EXPECT_EQ(max(c,d), d);
}