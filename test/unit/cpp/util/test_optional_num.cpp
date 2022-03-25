#include <iostream>

#include <boost/core/ignore_unused.hpp>

#include "util/optional_num.hpp"

#include "gtest/gtest.h"

using namespace steps::util;

TEST(OptionalNum, isValue) {

  auto v3 = OptionalNum<int>(3);

  ASSERT_EQ(v3.value(), 3);
  ASSERT_EQ(*v3, 3);
  ASSERT_EQ(v3.value_or(4), 3);
  ASSERT_EQ(bool(v3), true);

  v3 = NoneType();

  ASSERT_EQ(v3.has_value(), false);
  ASSERT_EQ(bool(v3), false);
  ASSERT_EQ(v3.value_or(4), 4);
}

TEST(OptionalNum, isNone) {
  auto n = OptionalNum<int>();
  ASSERT_EQ(n.has_value(), false);
  auto n2 = OptionalNum<int>(NoneType());
  ASSERT_EQ(n2.has_value(), false);
}

#ifndef NDEBUG
TEST(OptionalNum, asserts) {
  using optint = OptionalNum<int>;
  ASSERT_DEATH({ optint var(optint::none_value()); }, "");
  auto v3 = OptionalNum<int>(3);
  ASSERT_DEATH(v3.operator=(optint::none_value()), "");
}
#endif // NDEBUG
