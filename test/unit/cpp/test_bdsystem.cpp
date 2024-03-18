#include <cmath>
#include <limits>

#include <steps/solver/efield/bdsystem.hpp>
#include <steps/solver/efield/bdsystem_lapack.hpp>
#include <steps/solver/efield/linsystem.hpp>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "lapack_common.hpp"

using namespace steps::solver::efield;  // NOLINT

TEST_CASE("LinSystem_VVector") {
    double x[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    VVector v(8, &x[2]);

    // access through virtual methods
    AVector& a(v);

    REQUIRE_THAT(a.get(0), Catch::Matchers::WithinULP(2.0, 4));
    REQUIRE_THAT(a.get(7), Catch::Matchers::WithinULP(9.0, 4));
    REQUIRE(a.size() == 8);

    a.set(1, -1.0);
    REQUIRE_THAT(a.get(1), Catch::Matchers::WithinULP(-1.0, 4));

    // direct access through operator[]

    REQUIRE_THAT(v[0], Catch::Matchers::WithinULP(2.0, 4));
    REQUIRE_THAT(v[1], Catch::Matchers::WithinULP(-1.0, 4));
}


TEST_CASE("LinSystem_BDMatrix") {
    BDMatrix m(7, 2);  // 7x7 with half bw 2

    // direct access through operator[]

    m[0][0] = 1.0;
    m[0][2] = 2.0;
    m[2][0] = 3.0;
    m[6][5] = 4.0;

    REQUIRE_THAT(m[0][0], Catch::Matchers::WithinULP(1.0, 4));
    REQUIRE_THAT(m[0][2], Catch::Matchers::WithinULP(2.0, 4));
    REQUIRE_THAT(m[2][0], Catch::Matchers::WithinULP(3.0, 4));
    REQUIRE_THAT(m[6][5], Catch::Matchers::WithinULP(4.0, 4));

    // access through virtual methods
    AMatrix& a(m);

    REQUIRE_THAT(a.get(0, 0), Catch::Matchers::WithinULP(1.0, 4));
    REQUIRE_THAT(a.get(0, 2), Catch::Matchers::WithinULP(2.0, 4));
    REQUIRE_THAT(a.get(2, 0), Catch::Matchers::WithinULP(3.0, 4));
    REQUIRE_THAT(a.get(6, 5), Catch::Matchers::WithinULP(4.0, 4));

    REQUIRE(a.nRow() == 7);
    REQUIRE(a.nCol() == 7);

    a.set(5, 6, -1.0);
    REQUIRE_THAT(a.get(5, 6), Catch::Matchers::WithinULP(-1.0, 4));
    REQUIRE_THAT(m[5][6], Catch::Matchers::WithinULP(-1.0, 4));
}

using LinSystemImpls = std::tuple<BDSystem, BDSystemLapack>;


TEMPLATE_LIST_TEST_CASE("LinSystemImplTest_Diagonal", "", LinSystemImpls) {
    typedef TestType Impl;
    typedef typename Impl::matrix_type matrix_type;
    typedef typename Impl::vector_type vector_type;

    // diagonal-only test (halfbw=0)
    constexpr int n = 6;
    double diag[n] = {8, 2, 4, 32, 16, 64};
    double y[n] = {3, 5, 7, 9, 11, 13};

    double x_expected[n];
    for (size_t i = 0; i < n; ++i) {
        x_expected[i] = y[i] / diag[i];
    }

    Impl B(n, 0);

    matrix_type& A = B.A();
    for (int i = 0; i < n; ++i) {
        A.set(i, i, diag[i]);
    }

    vector_type& b = B.b();
    for (int i = 0; i < n; ++i) {
        b.set(i, y[i]);
    }

    B.solve();

    const vector_type& x = B.x();
    // (with powers of two in diagonal, expect result to be exact, but use
    // 4ulp test anyhow)
    for (int i = 0; i < n; ++i) {
        REQUIRE_THAT(x[i], Catch::Matchers::WithinULP(x_expected[i], 4));
    }
}

TEMPLATE_LIST_TEST_CASE("LinSystemImplTest_Bw7", "", LinSystemImpls) {
    typedef TestType Impl;
    typedef typename Impl::matrix_type matrix_type;
    typedef typename Impl::vector_type vector_type;

    constexpr int n = 8;
    constexpr int h = 3;  // half-bandwidth
    double A_full[n][n] = {{2, -1, -2, .1, 0, 0, 0, 0},
                           {-.5, 1.5, 0, -1, .2, 0, 0, 0},
                           {-9.3, -.4, -3, .2, -3, .3, 0, 0},
                           {.2, .1, 0, 4, -.2, -1, .4, 0},
                           {0, -.1, -.7, -.6, 3, 0, .2, .1},
                           {0, 0, -.3, -1.1, -2.1, 1.1, -.1, .2},
                           {0, 0, 0, -.6, 3, 0, 7, .3},
                           {0, 0, 0, 0, 2, -.3, -.5, 2}};

    double x0[n] = {1, 2, 3, 4, 5, 6, 7, 8};
    double y[n];

    // forward calculate y = A x0;
    for (int i = 0; i < n; ++i) {
        y[i] = 0;
        for (size_t j = 0; j < n; ++j) {
            y[i] += A_full[i][j] * x0[j];
        }
    }

    Impl B(n, h);

    matrix_type& A = B.A();
    for (int i = 0; i < n; ++i) {
        int jmin = std::max(0, i - h);
        int jmax = std::min(n - 1, i + h);

        for (int j = jmin; j <= jmax; ++j) {
            A.set(i, j, A_full[i][j]);
        }
    }

    vector_type& b = B.b();
    for (int i = 0; i < n; ++i) {
        b.set(i, y[i]);
    }

    B.solve();

    const vector_type& x = B.x();

    double k = estimate_condition(reinterpret_cast<const double*>(A_full), n);
    double eta = std::numeric_limits<double>::epsilon() * k;

    REQUIRE(eta < 0.25);  // our matrix should not be that ill-conditioned!

    // allow a factor of 4 for misestimation of the condition number.
    double relerr = 1 / (1 - eta * 4) - 1;

    for (int i = 0; i < n; ++i) {
        REQUIRE_THAT(x0[i], Catch::Matchers::WithinAbs(x.get(i), std::abs(x[i]) * relerr));
    }
}
