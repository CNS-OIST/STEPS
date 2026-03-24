/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2026 Okinawa Institute of Science and Technology, Japan.
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
#include <mpi.h>
#include <numeric>
#include <petscsys.h>

#include "util/petsc.hpp"

namespace steps::solver::efield {

/// c-tor
dVSolverPETSC::dVSolverPETSC(MPI_Comm petsc_comm) {
    PETSC_COMM_WORLD = petsc_comm;
    // Initialize PETSC (also MPI if not already done)
    auto err = PetscInitialize(nullptr, nullptr, nullptr, nullptr);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // Create vectors for rhs and solution
    err = VecCreate(PETSC_COMM_WORLD, &px);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // Create matrix for lhs
    err = MatCreate(PETSC_COMM_WORLD, &pA);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // Create Krylov solver
    err = KSPCreate(PETSC_COMM_WORLD, &pKsp);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    //    KSPSetComputeSingularValues(pKsp, PETSC_TRUE);
    //    KSPSetInitialGuessNonzero(pKsp,PETSC_TRUE);
    //    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    //    PetscViewerSetType(viewer, PETSCVIEWERASCII);
    //    PetscLogDefaultBegin();

    int mpi_sz;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_sz);
    petsc_locsizes.resize(mpi_sz);
    petsc_displ.resize(mpi_sz);
}

/// d-tor
dVSolverPETSC::~dVSolverPETSC() {
    //   PetscLogView(viewer);
    // Destroy all objects
    auto err = VecDestroy(&px);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = VecDestroy(&pb);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = MatDestroy(&pA);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = KSPDestroy(&pKsp);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    // Finalize PETSc
    err = PetscFinalize();
    CHKERRABORT(PETSC_COMM_WORLD, err);
}

/// Initialize mesh and sparsity pattern
void dVSolverPETSC::initMesh(TetMesh* mesh) {
    dVSolverBase::initMesh(mesh);

    // deltaV.resize(pNVerts);

    pIdxToVert.reserve(pNVerts);

    // Setup Vectors
    auto err = VecSetSizes(px, PETSC_DECIDE, pNVerts);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = VecSetType(px, VECMPI);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = VecGetOwnershipRange(px, &prbegin, &prend);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = VecDuplicate(px, &pb);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = VecGetLocalSize(px, &pNlocal);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // Setup Matrix
    err = MatSetSizes(pA, pNlocal, pNlocal, pNVerts, pNVerts);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = MatSetType(pA, MATMPIAIJ);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    PetscInt idx, jdx, n_con;
    std::vector<PetscInt> d_nnz(pNlocal,
                                0);  // # nnz in rows of DIAGONAL portion of local submatrix
    std::vector<PetscInt> o_nnz(pNlocal,
                                0);  // # nnz in rows of OFF-DIAG portion of local submatrix

    for (auto vi: vertex_id_t::range(pNVerts)) {
        VertexElement* ve = mesh->getVertex(vi);
        idx = ve->getIDX();  // get idx of vertex
        // NONOPTIMAL! In future 1) check if rbegin<=idx<rend 2) if yes pushback
        pIdxToVert.push_back(ve);  // now we know how is vertex number idx
        n_con = ve->getNCon();     // get how many neighbours

        // fill sparsity template
        if (idx >= prbegin && idx < prend) {
            ++d_nnz.at(idx - prbegin);
            for (PetscInt j = 0; j < n_con; ++j) {
                jdx = ve->nbrIdx(j);
                if (jdx >= prbegin && jdx < prend) {
                    ++d_nnz.at(idx - prbegin);
                } else {
                    ++o_nnz.at(idx - prbegin);
                }
            }
        }
    }

    err = MatMPIAIJSetPreallocation(pA, 0, d_nnz.data(), 0, o_nnz.data());
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = MatSetUp(pA);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // Fill the matrix and assemble it for the first time with all non-zero values present
    // This is necessary because the first assembly can change the non-zero structure
    for (PetscInt i = prbegin; i < prend; ++i) {
        VertexElement* ve = pIdxToVert[i];
        std::vector<PetscInt> idx_columns(ve->getNCon() + 1);
        std::vector<PetscScalar> val_columns(ve->getNCon() + 1);
        for (auto inbr = 0u; inbr < ve->getNCon(); ++inbr) {
            int j = ve->nbrIdx(inbr);
            idx_columns[inbr + 1] = j;
            val_columns[inbr + 1] = -1;
        }
        idx_columns[0] = i;
        val_columns[0] = ve->getNCon();
        err = MatSetValues(
            pA, 1, &i, idx_columns.size(), idx_columns.data(), val_columns.data(), INSERT_VALUES);
        CHKERRABORT(PETSC_COMM_WORLD, err);
    }
    err = MatAssemblyBegin(pA, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = MatAssemblyEnd(pA, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    /// First, Allgather to get all the solution vector sizes
    /// FIx for the powerpc64 platform, which incorrectly converts PetscInt to int
    int pnl = pNlocal;
    MPI_Allgather(&pnl, 1, MPI_INT, &petsc_locsizes[0], 1, MPI_INT, PETSC_COMM_WORLD);

    /// Now get all the values of the solution in one global array
    petsc_displ[0] = 0;
    if (petsc_displ.size() > 1) {
        std::partial_sum(petsc_locsizes.begin(), petsc_locsizes.end() - 1, petsc_displ.begin() + 1);
    }

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


    // Now, let's create a vector with the idxs of the triangles tri such that at least of the
    // vertices of tri are in the interval [prbegin, prend) This is necessary to have scalable
    // assembly time
    for (auto tri_idx: triangle_local_id::range(pNTris)) {
        auto* triv = pMesh->getTriangle(tri_idx);
        if (std::any_of(triv, triv + 3, [this](vertex_id_t i) {
                return i.get() >= static_cast<vertex_id_t::value_type>(prbegin) &&
                       i.get() < static_cast<vertex_id_t::value_type>(prend);
            })) {
            loc_tris.push_back(tri_idx);
        }
    }
}

/// Assemble and solve linear system to get potential at time t(n+1)
void dVSolverPETSC::advance(double dt) {
    // Get current clamp contribution
    std::copy(&pVertCurClamp[prbegin], &pVertCurClamp[prend], &pVertCur[prbegin]);
    for (auto idx: loc_tris) {
        double c = (pTriCur[idx.get()] + pTriCurClamp[idx.get()]) / 3.0;
        auto* triv = pMesh->getTriangle(idx);
        std::for_each(triv, triv + 3, [this, c](vertex_id_t i) { pVertCur[i.get()] += c; });
    }

    double oodt = 1.0 / dt;

    auto err = MatZeroEntries(pA);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    //    VecSet(pb,0.0);
    //    VecSet(px,0.0);

    std::vector<PetscInt> idx_rhs(pNlocal);
    std::iota(idx_rhs.begin(), idx_rhs.end(), prbegin);
    std::vector<PetscScalar> values_rhs(pNlocal);

    // iterate over vertices in local range
    for (PetscInt i = prbegin; i < prend; ++i) {
        VertexElement* ve = pIdxToVert[i];
        // case 1: vertex is on Clamp
        if (pVertexClamp[i] != 0) {
            values_rhs.at(i - prbegin) = 0.;
            err = MatSetValue(pA, i, i, 1., INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_WORLD, err);
        }
        // case 2: no clamp, get all Current Contributions
        else {
            double rhs = pVertCur[i] + pGExt[i] * (pVExt - pV[i]);
            double Aii = ve->getCapacitance() * oodt + pGExt[i];
            // indexes of columns j for row i where we want to insert
            std::vector<PetscInt> idx_columns(ve->getNCon() + 1);
            // respective values we want to insert
            std::vector<PetscScalar> val_columns(ve->getNCon() + 1);
            for (auto inbr = 0u; inbr < ve->getNCon(); ++inbr) {
                int j = ve->nbrIdx(inbr);
                double cc = ve->getCC(inbr);
                rhs += cc * (pV[j] - pV[i]);
                Aii += cc;
                idx_columns[inbr + 1] = j;
                val_columns[inbr + 1] = -cc;
            }
            idx_columns[0] = i;
            val_columns[0] = Aii;
            err = MatSetValues(pA,
                               1,
                               &i,
                               idx_columns.size(),
                               idx_columns.data(),
                               val_columns.data(),
                               INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_WORLD, err);
            values_rhs.at(i - prbegin) = rhs;
        }
    }

    // set rhs all at once to optimize
    err = VecSetValues(pb, pNlocal, idx_rhs.data(), values_rhs.data(), INSERT_VALUES);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // Set inital guess (???)
    // VecSetValues(px, pNlocal, idx_rhs.data(), &deltaV[prbegin], INSERT_VALUES);

    // Assemble LHS, rhs, intial guess
    err = MatAssemblyBegin(pA, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = MatAssemblyEnd(pA, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = VecAssemblyBegin(pb);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = VecAssemblyEnd(pb);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = VecAssemblyBegin(px);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = VecAssemblyEnd(px);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // LHS used for preconditioning
    err = KSPSetOperators(pKsp, pA, pA);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    // Solver: Conjugate Gradient
    err = KSPSetType(pKsp, KSPPIPECG);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    // Preconditioner
    err = KSPGetPC(pKsp, &pPc);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = PCSetType(pPc, PCPBJACOBI);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // Tolerances for iterative solver
    // KSPSetTolerances(pKsp, 1.e-1, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT);

    // Call solver and print statistics
    err = KSPSolve(pKsp, pb, px);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    PetscScalar* larr;
    err = VecGetArray(px, &larr);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    std::vector<PetscScalar> deltaV(pNVerts);
    MPI_Allgatherv(larr,
                   pNlocal,
                   MPI_PETSC_SCALAR,
                   deltaV.data(),
                   petsc_locsizes.data(),
                   petsc_displ.data(),
                   MPI_PETSC_SCALAR,
                   PETSC_COMM_WORLD);

    // update membrane potential
    for (auto i = 0u; i < pNVerts; ++i) {
        if (pVertexClamp[i] == 0) {
            pV[i] += util::petsc::to_real(deltaV[i]);
        }
    }

    err = VecRestoreArray(px, &larr);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // reset pTriCur for caller contributions
    std::fill(pTriCur.begin(), pTriCur.end(), 0.0);
}

/// Init function, here for debugging
void dVSolverPETSC::init() {
    advance(0.3);
}

}  // namespace steps::solver::efield
