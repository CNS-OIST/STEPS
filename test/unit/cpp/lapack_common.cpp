#include <cstddef>
#include <limits>
#include <cassert>

#include "lapack_common.hpp"

// estimate condition number of general nxn matrix A in row-major
// order, using LAPACK.

extern "C" {
double dlange_(const char *norm,int *m,int *n,double *a,int *lda,double *work);
void dgetrf_(int *m,int *n,double *a,int *lda,int *ipiv,int *info);
void dgecon_(const char *norm,int *n,double *a,int *lda,double *anorm,double *rcond,double *work,int *iwork,int *info);
}

void transpose_copy(const double *A,double *B,size_t m,size_t n) {
    // naive cache-unfriendly transpose fine for our purposes here!
    for (size_t i=0;i<m;++i)
        for (size_t j=0;j<n;++j)
            B[j*m+i]=A[i*n+j];
}

double estimate_condition(const double *A,size_t n) {
    double *T=new double[n*n];
    transpose_copy(A,T,n,n);

    auto lda= static_cast<int>(n);
    auto n_= static_cast<int>(n);
    int info=0;
    double anorm=0;
    double *work=new double[4*n];

    // Inf-norm calculation
    anorm=dlange_("I",&n_,&n_,T,&lda,work);

    // LU decomposition
    int *ipiv=new int[n];
    dgetrf_(&n_,&n_,T,&lda,ipiv,&info);
    assert(info>=0);

    // Estimate condition number
    double rcond=0;
    int *iwork=new int[n];
    dgecon_("I",&n_,T,&lda,&anorm,&rcond,work,iwork,&info);
    assert(info==0);

    delete[] ipiv;
    delete[] iwork;
    delete[] work;
    delete[] T;

    return rcond==0?std::numeric_limits<double>::infinity():1.0/rcond;
}

