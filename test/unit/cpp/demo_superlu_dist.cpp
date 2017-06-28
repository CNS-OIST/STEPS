// Still a demo at this point, not a test.
//
#include <iostream>
#include <iterator>
#include <iomanip>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <array>
#include <vector>

#include <mpi.h>

extern "C" {
#include "third_party/superlu_dist-4.1/src/superlu_ddefs.h"
}

#include "steps/solver/efield/bdsystem.hpp"
using namespace steps::solver::efield;

/*

#include "steps/solver/efield/linsystem.hpp"
#include "steps/solver/efield/bdsystem_lapack.hpp"

#include "gtest/gtest.h"
*/

// coefficients and b for solving Laplace equation Lx=b with
// mixed D and v.N. boundary conditions on 2-d grid.

struct laplace_coeffs {
    int ni,nj,N;
    struct entry_type { int r,s; double value; };

    std::vector<double> b;
    std::vector<char> clamped;

    laplace_coeffs(int ni_,int nj_): ni(ni_), nj(nj_), N(ni_*nj_), b(N,0.0), clamped(N,false) {
        if (ni<=0 || nj<=0) throw std::invalid_argument("");
    }

    size_t size() const { return N>0?(size_t)N:0; }

    void clamp(int i,int j,double v) {
        int k=i*nj+j;
        b[k]=v;
        clamped[k]=true;
    }

    void unclamp(int i,int j) {
        int k=i*nj+j;
        b[k]=0;
        clamped[k]=false;
    }
        
    void unclamp_all() {
        clamped.assign(N,false);
        b.assign(N,0);
    }

    int halfbw() const {
        return ni>1?nj:1;
    }

    struct row {
        int r;
        int nj;
        int n;
        bool id;

        std::array<entry_type,5> entries;
        size_t e_count;

        row(int r_,int nj_,int n_,bool id_): r(r_), nj(nj_), n(n_), id(id_), e_count(0) {
            if (r<0 || nj<1 || n<1) return;

            if (id) {
                entries[e_count++]={r,r,1.0};
            }
            else {
                int i=r/nj;
                int j=r-i*nj;

                if (i>=1)
                    entries[e_count++]={r,r-nj,-1.0};
                if (j>=1)
                    entries[e_count++]={r,r-1,-1.0};
                int diag_count=e_count++;
                if (j+1<nj)
                    entries[e_count++]={r,r+1,-1.0};
                if (i+nj<n)
                    entries[e_count++]={r,r+nj,-1.0};
                entries[diag_count]={r,r,(double)(e_count-1)};
            }
        }
 
        double operator[](int s) const {
            for (size_t k=0;k<e_count;++k)
                if (entries[k].s==s) return entries[k].value;
            return 0.0;
        }

        std::array<entry_type,4>::const_iterator begin() const { return entries.begin(); }
        std::array<entry_type,4>::const_iterator end() const { return entries.begin()+e_count; }
    };

    row operator[](int r) const {
        if (r<0 || r >=N) throw std::invalid_argument("");
        return row(r,nj,N,clamped[r]);
    }

    struct row_iterator: std::iterator<std::forward_iterator_tag,row> {
        const laplace_coeffs *L;
        int r;

        row_iterator(): L(nullptr), r(-1) {}

        explicit row_iterator(const laplace_coeffs &lc,int r_=0): L(&lc), r(r_) {
            if (r<0 || r>=L->N) L=nullptr;
        }

        bool operator==(const row_iterator &them) const {
            return L==them.L && (!L || r==them.r);
        }
    
        bool operator!=(const row_iterator &them) const {
            return !(*this==them);
        }

        row_iterator &operator++() {
            ++r;
            if (L && r>=L->N) L=nullptr;
            return *this;
        }

        row_iterator operator++(int) {
            row_iterator x{*this};
            ++*this;
            return x;
        }

        row operator*() const { return row(r,L->nj,L->N,L->clamped[r]); }
    };

    row_iterator begin() const { return row_iterator(*this); }
    row_iterator end() const { return row_iterator(); }
};


void emit_nonzero_coeffs(std::ostream &O,const laplace_coeffs &L) {
    for (const auto &r: L) {
        for (const auto &e: r) {
            O << e.value << " ";
        }
        O << "\n";
    }
}

void emit_matrix_elems(std::ostream &O,const laplace_coeffs &L) {
    size_t N=L.size();
    for (unsigned i=0; i<N; ++i) {
        for (unsigned j=0; j<N; ++j) {
            O << std::setw(5) << L[i][j];
        }
        O << "\n";
    }
}



int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int nproc,rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    laplace_coeffs L(9,4);
    int N=(int)L.size();

    if (rank==0) {
        std::cout << "L(9,4)\n";
        std::cout << "N=" << L.N << "\nHalf bandwidth=" << L.halfbw() << "\n";
        std::cout << "Non-zero coeffs by row:\n";
        emit_nonzero_coeffs(std::cout,L);

        std::cout << "\nAs matrix:\n";
        emit_matrix_elems(std::cout,L);
    }

    // clamp top and bottom
    int ni=L.ni;
    int nj=L.nj;
    const double top_v=5.0,bottom_v=1.0;
    for (int j=0;j<nj;++j) L.clamp(0,j,top_v),L.clamp(ni-1,j,bottom_v);

    if (rank==0) {
        std::cout << "\nL(9,4) top and bottom clamped\n";
        std::cout << "\nAs matrix:\n";
        emit_matrix_elems(std::cout,L);
        std::cout << "\nb:\n";
        std::copy(L.b.begin(),L.b.end(),std::ostream_iterator<double>(std::cout," "));
        std::cout << "\n";
    }

    // solve with bdsystem
    int h=(int)L.halfbw();
    BDSystem S(N,L.halfbw());

    typedef typename BDSystem::matrix_type matrix_type;
    typedef typename BDSystem::vector_type vector_type;

    matrix_type &A=S.A();

    for (int i=0;i<N;++i)
        for (int j=i-h;j<=i+h;++j)
            A.set(i,j,L[i][j]);

        
    vector_type &b=S.b();
    for (int i=0;i<N;++i) b.set(i,L.b[i]);
    
    S.solve();

    if (rank==0) {
        std::cout << "\nsolve x:\n";
        const auto &x=S.x();
        for (int i=0;i<N;++i)
            std::cout << x.get(i) << " ";
        std::cout << "\n";
    }

    gridinfo_t grid;
    superlu_gridinit(MPI_COMM_WORLD,nproc,1,&grid);

    superlu_options_t slu_options;
    set_default_options_dist(&slu_options);

    // set up (distributed) A matrix.
    SuperMatrix slu_A;
    slu_A.Stype=SLU_NR_loc;
    slu_A.Dtype=SLU_D;
    slu_A.Mtype=SLU_GE;
    slu_A.nrow=N;
    slu_A.ncol=N;

    std::vector<double> loc_values;
    std::vector<int> loc_col,loc_roff;

    int loc_m0=(N*grid.iam)/(grid.nprow*grid.npcol);
    int loc_m=(N*(1+grid.iam))/(grid.nprow*grid.npcol)-loc_m0;

    std::cout << "rank " << rank << ": m0=" << loc_m0 << " m=" << loc_m << "\n";

    for (int r=loc_m0;r<loc_m0+loc_m;++r) {
        loc_roff.push_back((int)loc_values.size());
        for (const auto &entry: L[r]) {
            loc_col.push_back(entry.s);
            loc_values.push_back(entry.value);
        }
    }
    loc_roff.push_back((int)loc_values.size());

    NRformat_loc data;
    data.m_loc=loc_m;
    data.fst_row=loc_m0;
    data.nnz_loc=(int)loc_values.size();
    data.nzval=(void *)&loc_values[0];
    data.rowptr=&loc_roff[0];
    data.colind=&loc_col[0];

    slu_A.Store=(void *)&data;

    // set up (distribued) B matrix.
    std::vector<double> slu_b(&L.b[loc_m0],&L.b[loc_m0]+loc_m);

    // allocate and initialise solver, permutation and LU data
    SOLVEstruct_t solve_data;
    ScalePermstruct_t perm_data;
    ScalePermstructInit(N, N, &perm_data);
    LUstruct_t lu_data;
    LUstructInit(N, &lu_data);

    // call solver, get stats
    SuperLUStat_t slu_stats;
    int slu_info;
    PStatInit(&slu_stats);
    
    double berr;
    pdgssvx(&slu_options,&slu_A,&perm_data,&slu_b[0],loc_m,1,&grid,&lu_data,&solve_data,&berr,&slu_stats,&slu_info);
    PStatPrint(&slu_options,&slu_stats,&grid);
    PStatFree(&slu_stats);

    // deallocate superlu data
    ScalePermstructFree(&perm_data);
    Destroy_LU(N,&grid,&lu_data);
    LUstructFree(&lu_data);
    if (slu_options.SolveInitialized) dSolveFinalize(&slu_options,&solve_data);
    superlu_gridexit(&grid);

    // gather data
    if (rank==0) {
        std::vector<int> dist_m(nproc);
        std::vector<int> dist_m0(nproc);
        MPI_Gather(&loc_m,1,MPI_INT,&dist_m[0],1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Gather(&loc_m0,1,MPI_INT,&dist_m0[0],1,MPI_INT,0,MPI_COMM_WORLD);

        std::vector<double> x(N);
        MPI_Gatherv(&slu_b[0],loc_m,MPI_DOUBLE,&x[0],&dist_m[0],&dist_m0[0],MPI_DOUBLE,0,MPI_COMM_WORLD);

        std::cout << "\nsuperlu solve x:\n";
        for (int i=0;i<N;++i)
            std::cout << x[i] << " ";
        std::cout << "\n";
    }
    else {
        MPI_Gather(&loc_m,1,MPI_INT,nullptr,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Gather(&loc_m0,1,MPI_INT,nullptr,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Gatherv(&slu_b[0],loc_m,MPI_DOUBLE,nullptr,nullptr,nullptr,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }

    MPI_Finalize();
}

            
        
