/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

#include "dVsolver_petsc.hpp"

#include <algorithm>
#include <numeric>
#include <fstream>
#include <mpi.h>


namespace steps {
namespace solver {
namespace efield {

/// c-tor
dVSolverPETSC::dVSolverPETSC() {
    // Initialize PETSC (also MPI if not already done)
    PetscInitialize(nullptr, nullptr, nullptr, nullptr);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create vectors for rhs and solution
    VecCreate(PETSC_COMM_WORLD, &px);

    // Create matrix for lhs
    MatCreate(PETSC_COMM_WORLD,&pA);

    // Create Krylov solver
    KSPCreate(PETSC_COMM_WORLD, &pKsp);
//    KSPSetComputeSingularValues(pKsp, PETSC_TRUE);
//    KSPSetInitialGuessNonzero(pKsp,PETSC_TRUE);
//    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
//    PetscViewerSetType(viewer, PETSCVIEWERASCII);
//    PetscLogDefaultBegin();

    int mpi_sz;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_sz);
    petsc_locsizes.resize(mpi_sz);
    petsc_displ.resize(mpi_sz);

}

/// d-tor
dVSolverPETSC::~dVSolverPETSC(){
//   PetscLogView(viewer);
   // Destroy all objects
   VecDestroy(&px);
   VecDestroy(&pb);
   MatDestroy(&pA);
   KSPDestroy(&pKsp);
   // Finalize PETSc
   PetscFinalize();
}


/// Initialize mesh and sparsity pattern
void dVSolverPETSC::initMesh(TetMesh *mesh) {

    dVSolverBase::initMesh(mesh);

//deltaV.resize(pNVerts);

    pIdxToVert.reserve(pNVerts);

    // Setup Vectors
    VecSetSizes(px,PETSC_DECIDE, pNVerts);
    VecSetType(px, VECMPI);
    VecGetOwnershipRange(px, &prbegin, &prend);
    VecDuplicate(px, &pb);
    VecGetLocalSize(px, &pNlocal);

    // Setup Matrix
    MatSetSizes(pA, pNlocal, pNlocal, pNVerts, pNVerts);
    MatSetType(pA, MATMPIAIJ);

    PetscInt idx, jdx, n_con;
    std::vector<PetscInt> d_nnz(pNlocal,0); // # nnz in rows of DIAGONAL portion of local submatrix
    std::vector<PetscInt> o_nnz(pNlocal,0); // # nnz in rows of OFF-DIAG portion of local submatrix

    for (uint i=0; i<pNVerts; ++i) {
        VertexElement* ve = mesh->getVertex(i);
        idx = ve->getIDX();     // get idx of vertex
        // NONOPTIMAL! In future 1) check if rbegin<=idx<rend 2) if yes pushback
        pIdxToVert.push_back(ve);    // now we know how is vertex number idx
        n_con = ve->getNCon();  // get how many neighbours

        //fill sparsity template
        if (idx>=prbegin && idx<prend) {
            ++d_nnz.at(idx-prbegin);
            for (PetscInt j=0; j<n_con; ++j) {
                jdx = ve->nbrIdx(j);
                if (jdx>=prbegin && jdx<prend)
                    ++d_nnz.at(idx-prbegin);
                else
                    ++o_nnz.at(idx-prbegin);
            }
        }
    }

    MatMPIAIJSetPreallocation(pA,0,d_nnz.data(),0,o_nnz.data());
    MatSetUp(pA);

    /// First, Allgather to get all the solution vector sizes
    /// FIx for the powerpc64 platform, which incorrectly converts PetscInt to int
    int pnl = pNlocal;
    MPI_Allgather(&pnl, 1, MPI_INT, &petsc_locsizes[0], 1, MPI_INT, PETSC_COMM_WORLD);

    /// Now get all the values of the solution in one global array
    petsc_displ[0] = 0;
    if (petsc_displ.size() > 1)
        std::partial_sum(petsc_locsizes.begin(), petsc_locsizes.end()-1, petsc_displ.begin()+1);

//    // Compute load imbalance
//    double avg_load = pNVerts / double(petsc_locsizes.size());
//    double max_imbalance = std::abs(petsc_locsizes[0] - avg_load);
//    if (petsc_locsizes.size()>1) {
//        for (idx = 1; idx<petsc_locsizes.size(); ++idx) {
//            double imbalance_i = std::abs(petsc_locsizes[idx] - avg_load);
//            if (imbalance_i > max_imbalance)
//                max_imbalance = imbalance_i;
//        }
//    }
//    if (!rnk) printf("PETSc Mat and Vec load imbalance: %.2f %%\n", max_imbalance/avg_load);


    // Now, let's create a vector with the idxs of the triangles tri such that at least of the vertices of tri are in the interval [prbegin, prend)
    // This is necessary to have scalable assembly time
    for (uint tri_idx=0; tri_idx<pNTris; ++tri_idx) {
        auto *triv = pMesh->getTriangle(tri_idx);
        if (std::any_of(triv, triv + 3, [this](vertex_id_t i){return i.get() >= static_cast<vertex_id_t::value_type>(prbegin) && i.get() < static_cast<vertex_id_t::value_type>(prend) ;}))
            loc_tris.push_back(tri_idx);
    }
}


/// Assemble and solve linear system to get potential at time t(n+1)
void dVSolverPETSC::advance(double dt) {

    // Get current clamp contribution
    std::copy(&pVertCurClamp[prbegin], &pVertCurClamp[prend], &pVertCur[prbegin]);
    for (auto idx : loc_tris) {
        double c = (pTriCur[idx] + pTriCurClamp[idx]) / 3.0;
        auto *triv = pMesh->getTriangle(idx);
        std::for_each(triv, triv + 3, [this, c](vertex_id_t i){ pVertCur[i.get()] +=c;});
    }

    double oodt = 1.0/dt;


    MatZeroEntries(pA);
//    VecSet(pb,0.0);
//    VecSet(px,0.0);

    std::vector<PetscInt>    idx_rhs(pNlocal);   std::iota(idx_rhs.begin(), idx_rhs.end(), prbegin);
    std::vector<double> values_rhs(pNlocal);


    // iterate over vertices in local range
    for (PetscInt i=prbegin; i<prend; ++i) {
        VertexElement *ve = pIdxToVert[i];
        // case 1: vertex is on Clamp
        if (pVertexClamp[i]) {
            values_rhs.at(i-prbegin) = 0.;
            MatSetValue(pA,i,i,1.,INSERT_VALUES);
        }
        // case 2: no clamp, get all Current Contributions
        else {
            double rhs = pVertCur[i] + pGExt[i] * (pVExt - pV[i]);
            double Aii = ve->getCapacitance()*oodt + pGExt[i];
            // indexes of columns j for row i where we want to insert
            std::vector<PetscInt> idx_columns(ve->getNCon()+1);
            // respective values we want to insert
            std::vector<double> val_columns(ve->getNCon()+1);
            for (auto inbr = 0u; inbr < ve->getNCon(); ++inbr) {
                int j = ve->nbrIdx(inbr);
                double cc = ve->getCC(inbr);
                rhs += cc * (pV[j] - pV[i]);
                Aii += cc;
                idx_columns[inbr+1] = j;
                val_columns[inbr+1] = -cc;
            }
            idx_columns[0] = i;
            val_columns[0] = Aii;
            MatSetValues(pA, 1, &i, idx_columns.size(), idx_columns.data(), val_columns.data(), INSERT_VALUES);
            values_rhs.at(i-prbegin) = rhs;
        }
    }


    // set rhs all at once to optimize
    VecSetValues(pb, pNlocal, idx_rhs.data(), values_rhs.data(), INSERT_VALUES);

// Set inital guess (???)
//VecSetValues(px, pNlocal, idx_rhs.data(), &deltaV[prbegin], INSERT_VALUES);

    // Assemble LHS, rhs, intial guess
    MatAssemblyBegin(pA,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(pA,MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(pb);
    VecAssemblyEnd(pb);
    VecAssemblyBegin(px);
    VecAssemblyEnd(px);


    // LHS used for preconditioning
    KSPSetOperators(pKsp, pA, pA);
    // Solver: Conjugate Gradient
    KSPSetType(pKsp, KSPPIPECG);
    // Preconditioner
    KSPGetPC(pKsp, &pPc);
    PCSetType(pPc, PCPBJACOBI);

    // Tolerances for iterative solver
    // KSPSetTolerances(pKsp, 1.e-1, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT);

    // Call solver and print statistics
    KSPSolve(pKsp, pb, px);


    double *larr;
    VecGetArray(px, &larr);

    std::vector<double> deltaV(pNVerts);
    MPI_Allgatherv(larr, pNlocal, MPI_DOUBLE, &deltaV[0], &petsc_locsizes[0], &petsc_displ[0], MPI_DOUBLE, PETSC_COMM_WORLD);


    // update membrane potential
    for (auto i = 0u; i < pNVerts; ++i)
        if (pVertexClamp[i] == false)
            pV[i] += deltaV[i];

    VecRestoreArray(px, &larr);

    // reset pTriCur for caller contributions
    std::fill(pTriCur.begin(), pTriCur.end(), 0.0);

}




/// Init function, here for debugging
void dVSolverPETSC::init() {
    advance(0.3);
}




}  // namespace efield
}  // namespace solver
}  // namespace steps

