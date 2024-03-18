#include <cmath>

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <petscksp.h>

#include "steps/util/petsc.hpp"

TEST_CASE("PetscSystem_Vec") {
    Vec x, y;
    PetscInt N = 100;
    VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, N, &x);
    VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, N, &y);

    double val_x = 1., val_y = 10.;
    VecSet(x, val_x);
    VecSet(y, val_y);

    PetscScalar x_dot_y;
    VecDot(x, y, &x_dot_y);
    REQUIRE_THAT(steps::util::petsc::to_real(x_dot_y),
                 Catch::Matchers::WithinULP(N * val_x * val_y, 4));

    double vec_norm;
    VecNorm(x, NORM_2, &vec_norm);
    REQUIRE_THAT(vec_norm, Catch::Matchers::WithinULP(std::sqrt(N) * val_x, 4));
    VecNorm(y, NORM_2, &vec_norm);
    REQUIRE_THAT(vec_norm, Catch::Matchers::WithinULP(std::sqrt(N) * val_y, 4));

    VecDestroy(&x);
    VecDestroy(&y);
}


TEST_CASE("PetscSystem_KrylovSolver") {
    Mat A;
    Vec x, x_exact, b;
    KSP ksp;
    PC pc;
    PetscInt N = 100;
    double tol = 2.e-12;

    VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, N, &x);
    VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, N, &x_exact);
    VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, N, &b);
    MatCreateAIJ(
        PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, N, N, 10, PETSC_NULL, 10, PETSC_NULL, &A);
    KSPCreate(PETSC_COMM_WORLD, &ksp);

    PetscInt beg, end, loc_sz;
    VecGetOwnershipRange(x, &beg, &end);
    VecGetLocalSize(x, &loc_sz);

    PetscInt idx_x, idx_y[3];
    if (beg == 0) {
        beg = 1;
        idx_x = 0;
        idx_y[0] = 0;
        idx_y[1] = 1;
        PetscScalar vals[] = {2., -1.};
        MatSetValues(A, 1, &idx_x, 2, idx_y, vals, INSERT_VALUES);
    }
    if (end == N) {
        end = N - 1;
        idx_x = N - 1;
        idx_y[0] = N - 2;
        idx_y[1] = N - 1;
        PetscScalar vals[] = {-1., 2.};
        MatSetValues(A, 1, &idx_x, 2, idx_y, vals, INSERT_VALUES);
    }

    PetscScalar vals[] = {-1., 2., -1.};
    for (auto i = beg; i < end; ++i) {
        idx_y[0] = i - 1;
        idx_y[1] = i;
        idx_y[2] = i + 1;
        MatSetValues(A, 1, &i, 3, idx_y, vals, INSERT_VALUES);
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    VecSet(x_exact, 1.);
    MatMult(A, x_exact, b);


    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    KSPSetTolerances(ksp, tol / 10., PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    PCSetType(pc, PCGAMG);
    KSPSolve(ksp, b, x);

    VecAXPY(x, -1, x_exact);
    double err_norm;
    VecNorm(x, NORM_2, &err_norm);
    REQUIRE_THAT(err_norm, Catch::Matchers::WithinAbs(0., tol));

    VecDestroy(&x);
    VecDestroy(&x_exact);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);
}

int main(int argc, char* argv[]) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int result = Catch::Session().run(argc, argv);
    PetscFinalize();
    return result;
}
