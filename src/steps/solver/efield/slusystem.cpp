/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */


#include <stdexcept>
#include <algorithm>

#include <mpi.h>

#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/solver/efield/slusystem.hpp"

extern "C" {
#include "third_party/superlu_dist-4.1/src/superlu_ddefs.h"
}

#include "easylogging++.h"

namespace steps {
namespace solver {
namespace efield {

struct supermatrix_nc_view {
    explicit supermatrix_nc_view(SLU_NCMatrix &S) {
        M.Stype = SLU_NC;
        M.Dtype = SLU_D;
        M.Mtype = SLU_GE;
        M.nrow = S.pN;
        M.ncol = S.pN;
        M.Store = &data;

        data.nnz = S.pNnz;
        data.nzval = &S.pValues[0];
        data.rowind = &S.pRidx[0];
        data.colptr = &S.pCoff[0];
    }

    SuperMatrix M;
    NCformat data;
};

struct SLUData {
    int N;
    bool factored;
    bool keepperm;

    superlu_options_t options;
    gridinfo_t grid;
    ScalePermstruct_t perm;
    LUstruct_t lu;
    SuperLUStat_t stat;


    SLUData(MPI_Comm mpi_comm, int N_, bool keepperm_=true):
        N(N_), factored(false), keepperm(keepperm_)
    {
        set_default_options_dist(&options);
        options.IterRefine = NOREFINE;


        int nrow = 0;
        MPI_Comm_size(mpi_comm, &nrow);
        superlu_gridinit(mpi_comm, nrow, 1, &grid);

        ScalePermstructInit(N, N, &perm);
        LUstructInit(N, &lu);
        PStatInit(&stat);
    }

    ~SLUData() {
        Destroy_LU(N, &grid, &lu);
        LUstructFree(&lu);
        PStatFree(&stat);
        ScalePermstructFree(&perm);
        superlu_gridexit(&grid);
    }
};

SLUSystem::SLUSystem(const sparsity_template &S, MPI_Comm mpi_comm):
    pN((int)S.dim()), pA(S),
    pb(pN,0.0), px(pN,0.0),
    pb_view(pN,&pb[0]), px_view(pN,&px[0]),
    pBerr(0)
{
    slu.reset(new SLUData(mpi_comm, pN, true));
}

SLUSystem::~SLUSystem() {}

void SLUSystem::solve() {
    // use copy of A...
    SLU_NCMatrix Abis(pA);

    supermatrix_nc_view slu_A(Abis);
    std::copy(pb.begin(),pb.end(),px.begin());

    if (!slu->factored)
        slu->options.Fact = DOFACT;
    else if (slu->keepperm)
        slu->options.Fact = SamePattern_SameRowPerm;
    else
        slu->options.Fact = SamePattern;

    int info;

    pdgssvx_ABglobal(&slu->options, &slu_A.M, &slu->perm, &px[0], pN, 1,
        &slu->grid, &slu->lu, &pBerr, &slu->stat, &info);

    if (info>0) {
        if (info<=pN) {
            ProgErrLog("Singular U in LU decomposition");
        }
        else {
            SysErrLog("SuperLU memory allocation failure");
        }
    }

    slu->factored = true;
}

}}} // namespace steps::solver::efield
