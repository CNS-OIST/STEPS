/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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


#ifndef STEPS_SOLVER_EFIELD_DVSOLVER_PETSC_HPP
#define STEPS_SOLVER_EFIELD_DVSOLVER_PETSC_HPP 1

#include <algorithm>
#include <memory>
#include <vector>
#include <numeric>

#include <petscksp.h>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/efield/dVsolver.hpp"

namespace steps {
namespace solver {
namespace efield {

class dVSolverPETSC: public dVSolverBase {
public:
   
    /// c-tor (*calls PetscInitialize*)
    explicit dVSolverPETSC();

    /// d-tor (*calls PetscFinalize*)
    ~dVSolverPETSC();

    /// Initialize mesh and sparsity pattern
    void initMesh(TetMesh *mesh) override final;
    
    /// Assemble and solve linear system to get potential at time t(n+1) 
    void advance(double dt) override final;
    // Delete previous implementation in dVSolverBase for safety reasons
    void _advance() = delete; 

    /// Init function, here for debugging, remove later
    void init();

private:
    PetscInt prbegin, prend;
    PetscInt pNlocal;     // number of rows handled by this processor
    std::vector<VertexElement*> pIdxToVert;  // map each idx to relative vertex
    std::vector<uint> loc_tris; // vector with all the idxs of triangles on this petsc partition
    std::vector<int> petsc_locsizes;             
    std::vector<int> petsc_displ;
    Mat pA;             // lhs
    Vec pb, px;         // rhs and approximate solution : pA * px = pb
    KSP pKsp;           // Krylov solver
    PC pPc;             // preconditioner
//std::vector<double> deltaV;
//PetscViewer viewer;
};


}}} // namespace steps::efield::solver

#endif //STEPS_SIM_EFIELD_DVSOLVER_PETSC_HPP

