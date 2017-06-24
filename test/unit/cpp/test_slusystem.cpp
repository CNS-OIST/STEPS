#include <iostream>
#include <limits>
#include <cmath>

#include <mpi.h>

#include "steps/solver/efield/linsystem.hpp"
#include "steps/solver/efield/slusystem.hpp"

#include "gtest/gtest.h"

#include "lapack_common.hpp"

using namespace steps::solver::efield;

constexpr int Adim=8;
double AA[Adim][Adim]={
    {   2,  -1,  -2,   .1,    0,   0,   0,  0},
    { -.5, 1.5,   0,   -1,   .2,   0,   0,  0},
    {-9.3, -.4,  -3,   .2,   -3,  .3,   0,  0},
    {  .2,  .1,   0,    4,  -.2,  -1,  .4,  0},
    {   0, -.1, -.7,  -.6,    3,   0,  .2, .1},
    {   0,   0, -.3, -1.1, -2.1, 1.1, -.1, .2},
    {   0,   0,   0,  -.6,    3,   0,   7, .3},
    {   0,   0,   0,    0,    2, -.3, -.5,  2}
};

int main(int argc, char **argv) {
    int r=0;

    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc,&argv);
    r=RUN_ALL_TESTS();
    MPI_Finalize();
    return r;
}


TEST(LinSystem,SLUMatrix) {
    sparsity_template S(Adim);
    int nnz=0;

    for (int i=0; i<Adim; ++i)
        for (int j=0; j<Adim; ++j)
            if (AA[i][j]) ++nnz,S.insert(std::make_pair(i,j));

    SLU_NCMatrix A(S);
    EXPECT_EQ(nnz,A.nNz());

    for (int i=0; i<Adim; ++i)
        for (int j=0; j<Adim; ++j)
            if (AA[i][j]) A.set(i,j,AA[i][j]);

    for (int i=0; i<Adim; ++i)
        for (int j=0; j<Adim; ++j)
            ASSERT_EQ(A.get(i,j),AA[i][j]);

}

TEST(LinSystem,SLUMatrixRedundantSparsity) {
    sparsity_template S(Adim);
    int nnz=0;

    for (int i=0; i<Adim; ++i)
        for (int j=0; j<Adim; ++j)
            if (AA[i][j]) ++nnz,S.insert(std::make_pair(i,j));

    ASSERT_EQ(nnz,S.size());
            
    bool insert_again=true;
    for (int i=0; i<Adim; ++i)
        for (int j=0; j<Adim; ++j) {
            if (AA[i][j] && insert_again) S.insert(std::make_pair(i,j));
            insert_again=!insert_again;
        }

    ASSERT_EQ(nnz+nnz/2,S.size());

    SLU_NCMatrix A(S);
    ASSERT_EQ(nnz,A.nNz());
}

TEST(LinSystem,SLUMatrixBadSparsity) {
    sparsity_template S(Adim);

    for (int i=0; i<Adim; ++i)
        for (int j=0; j<Adim; ++j)
            if (AA[i][j]) S.insert(std::make_pair(i,j));

    sparsity_template Snegidx=S;
    Snegidx.insert(std::make_pair(-1,0));

    EXPECT_THROW((void)SLU_NCMatrix(Snegidx), steps::ProgErr);

    sparsity_template Soob=S;
    Soob.insert(std::make_pair(Adim,0));

    EXPECT_THROW((void)SLU_NCMatrix(Soob), steps::ProgErr);
}

TEST(LinSystem,SLUSystem) {
    typedef SLUSystem::vector_type vector_type;
    typedef SLUSystem::matrix_type matrix_type;

    constexpr size_t n=Adim;

    double x0[n];
    double y[n];

    // fill x0;
    for (int i=0; i<n; ++i) x0[i]=i+1;

    // forward calculate y = A x0;
    for (int i=0;i<n;++i) {
        y[i]=0;
        for (size_t j=0; j<n; ++j)
            y[i]+=AA[i][j]*x0[j];
    }

    sparsity_template S(n);
    for (int i=0; i<n; ++i)
        for (int j=0; j<n; ++j)
            if (AA[i][j]) S.insert(std::make_pair(i,j));


    SLUSystem L(S,MPI_COMM_WORLD);
    matrix_type &A=L.A();

    for (int i=0; i<Adim; ++i)
        for (int j=0; j<Adim; ++j)
            if (AA[i][j]) A.set(i,j,AA[i][j]);

    vector_type &b=L.b();
    for (int i=0;i<n;++i) b.set(i,y[i]);

    L.solve();

    const vector_type &x=L.x();

    double k=estimate_condition((const double *)AA,n);
    double eta=std::numeric_limits<double>::epsilon()*k;

    ASSERT_LT(eta,0.25); // our matrix should not be that ill-conditioned!

    // allow a factor of 4 for misestimation of the condition number.
    double relerr=1/(1-eta*4)-1;

    for (int i=0;i<n;++i) {
        EXPECT_NEAR(x0[i],x.get(i),std::abs(x[i])*relerr);
    }
}

