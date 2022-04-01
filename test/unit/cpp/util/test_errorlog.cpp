#define ENABLE_ASSERTLOG = 1

#include "easylogging++.h"
#include "util/error.hpp"
#include "util/init.hpp"
#include "gtest/gtest.h"

class ErrorLogTestEnvironment : public ::testing::Environment {
public:
  virtual void SetUp() { steps::init(); }
};

void assert_pass_test() { AssertLog((1 + 1) == 2); }
TEST(AssertLogTest, AssertPass) { EXPECT_NO_THROW(assert_pass_test()); }

void assert_fail_test() { AssertLog(1 == 2); }
TEST(AssertLogTest, AssertFail) {
  EXPECT_THROW(assert_fail_test(), steps::AssertErr);
}

void condition_errlog_test(bool error_condition) {
  if (error_condition) {
    ErrLogIf(true, "log if condition is true.");
  } else {
    ErrLogIf(false, "log if condition is true.");
  }
}
TEST(ErrLogTest, NoThrow) { EXPECT_NO_THROW(condition_errlog_test(false)); }
TEST(ErrLogTest, Throw) {
  EXPECT_THROW(condition_errlog_test(true), steps::Err);
}

void condition_not_impl_errlog_test(bool error_condition) {
  if (error_condition) {
    NotImplErrLogIf(true, "log if condition is true.");
  } else {
    NotImplErrLogIf(false, "log if condition is true.");
  }
}
TEST(NotImplErrLogTest, NoThrow) {
  EXPECT_NO_THROW(condition_not_impl_errlog_test(false));
}
TEST(NotImplErrLogTest, Throw) {
  EXPECT_THROW(condition_not_impl_errlog_test(true), steps::NotImplErr);
}

void condition_arg_errlog_test(bool error_condition) {
  if (error_condition) {
    ArgErrLogIf(true, "log if condition is true.");
  } else {
    ArgErrLogIf(false, "log if condition is true.");
  }
}
TEST(ArgErrLogTest, NoThrow) {
  EXPECT_NO_THROW(condition_arg_errlog_test(false));
}
TEST(ArgErrLogTest, Throw) {
  EXPECT_THROW(condition_arg_errlog_test(true), steps::ArgErr);
}

void condition_prog_errlog_test(bool error_condition) {
  if (error_condition) {
    ProgErrLogIf(true, "log if condition is true.");
  } else {
    ProgErrLogIf(false, "log if condition is true.");
  }
}
TEST(ProgErrLogTest, NoThrow) {
  EXPECT_NO_THROW(condition_prog_errlog_test(false));
}
TEST(ProgErrLogTest, Throw) {
  EXPECT_THROW(condition_prog_errlog_test(true), steps::ProgErr);
}

void condition_sys_errlog_test(bool error_condition) {
  if (error_condition) {
    SysErrLogIf(true, "log if condition is true.");
  } else {
    SysErrLogIf(false, "log if condition is true.");
  }
}
TEST(SysErrLogTest, NoThrow) {
  EXPECT_NO_THROW(condition_sys_errlog_test(false));
}
TEST(SysErrLogTest, Throw) {
  EXPECT_THROW(condition_sys_errlog_test(true), steps::SysErr);
}

void condition_io_errlog_test(bool error_condition) {
  if (error_condition) {
    IOErrLogIf(true, "log if condition is true.");
  } else {
    IOErrLogIf(false, "log if condition is true.");
  }
}
TEST(IOErrLogTest, NoThrow) {
  EXPECT_NO_THROW(condition_io_errlog_test(false));
}
TEST(IOErrLogTest, Throw) {
  EXPECT_THROW(condition_io_errlog_test(true), steps::IOErr);
}

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new ErrorLogTestEnvironment);
  return RUN_ALL_TESTS();
}
#undef ENABLE_ASSERTLOG