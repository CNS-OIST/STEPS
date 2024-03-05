#define ENABLE_ASSERTLOG = 1

#include "easylogging++.h"
#include "util/error.hpp"
#include "util/init.hpp"

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

static void assert_pass_test() {
    AssertLog((1 + 1) == 2);
}
TEST_CASE("AssertLogTest_AssertPass") {
    REQUIRE_NOTHROW(assert_pass_test());
}

void assert_fail_test() {
    AssertLog(1 == 2);
}
TEST_CASE("AssertLogTest_AssertFail") {
    REQUIRE_THROWS_AS(assert_fail_test(), steps::AssertErr);
}

static void condition_errlog_test(bool error_condition) {
    if (error_condition) {
        ErrLogIf(true, "log if condition is true.");
    } else {
        ErrLogIf(false, "log if condition is true.");
    }
}
TEST_CASE("ErrLogTest_NoThrow") {
    REQUIRE_NOTHROW(condition_errlog_test(false));
}
TEST_CASE("ErrLogTest_Throw") {
    REQUIRE_THROWS_AS(condition_errlog_test(true), steps::Err);
}

static void condition_not_impl_errlog_test(bool error_condition) {
    if (error_condition) {
        NotImplErrLogIf(true, "log if condition is true.");
    } else {
        NotImplErrLogIf(false, "log if condition is true.");
    }
}
TEST_CASE("NotImplErrLogTest_NoThrow") {
    REQUIRE_NOTHROW(condition_not_impl_errlog_test(false));
}
TEST_CASE("NotImplErrLogTest_Throw") {
    REQUIRE_THROWS_AS(condition_not_impl_errlog_test(true), steps::NotImplErr);
}

static void condition_arg_errlog_test(bool error_condition) {
    if (error_condition) {
        ArgErrLogIf(true, "log if condition is true.");
    } else {
        ArgErrLogIf(false, "log if condition is true.");
    }
}
TEST_CASE("ArgErrLogTest_NoThrow") {
    REQUIRE_NOTHROW(condition_arg_errlog_test(false));
}
TEST_CASE("ArgErrLogTest_Throw") {
    REQUIRE_THROWS_AS(condition_arg_errlog_test(true), steps::ArgErr);
}

static void condition_prog_errlog_test(bool error_condition) {
    if (error_condition) {
        ProgErrLogIf(true, "log if condition is true.");
    } else {
        ProgErrLogIf(false, "log if condition is true.");
    }
}
TEST_CASE("ProgErrLogTest_NoThrow") {
    REQUIRE_NOTHROW(condition_prog_errlog_test(false));
}
TEST_CASE("ProgErrLogTest_Throw") {
    REQUIRE_THROWS_AS(condition_prog_errlog_test(true), steps::ProgErr);
}

static void condition_sys_errlog_test(bool error_condition) {
    if (error_condition) {
        SysErrLogIf(true, "log if condition is true.");
    } else {
        SysErrLogIf(false, "log if condition is true.");
    }
}
TEST_CASE("SysErrLogTest_NoThrow") {
    REQUIRE_NOTHROW(condition_sys_errlog_test(false));
}
TEST_CASE("SysErrLogTest_Throw") {
    REQUIRE_THROWS_AS(condition_sys_errlog_test(true), steps::SysErr);
}

static void condition_io_errlog_test(bool error_condition) {
    if (error_condition) {
        IOErrLogIf(true, "log if condition is true.");
    } else {
        IOErrLogIf(false, "log if condition is true.");
    }
}
TEST_CASE("IOErrLogTest_NoThrow") {
    REQUIRE_NOTHROW(condition_io_errlog_test(false));
}
TEST_CASE("IOErrLogTest_Throw") {
    REQUIRE_THROWS_AS(condition_io_errlog_test(true), steps::IOErr);
}

int main(int argc, char* argv[]) {
    steps::init();
    return Catch::Session().run(argc, argv);
}
#undef ENABLE_ASSERTLOG
