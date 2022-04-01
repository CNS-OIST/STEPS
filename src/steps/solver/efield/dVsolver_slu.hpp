/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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


#ifndef STEPS_SOLVER_EFIELD_DVSOLVER_SLU_HPP
#define STEPS_SOLVER_EFIELD_DVSOLVER_SLU_HPP 1

#include <algorithm>
#include <memory>
#include <vector>

#include <mpi.h>

// STEPS headers.
#include "util/common.h"
#include "dVsolver.hpp"
#include "slusystem.hpp"

namespace steps {
namespace solver {
namespace efield {

class dVSolverSLU: public dVSolverBase {
public:
    explicit dVSolverSLU(MPI_Comm comm): pMpiComm(comm) {}

    void initMesh(TetMesh *mesh) override {
        dVSolverBase::initMesh(mesh);
        sparsity_template S(pNVerts);

        for (auto i = 0u; i < pNVerts; ++i) {
            VertexElement *ve = mesh->getVertex(i);

            int idx = ve->getIDX();
            int ncon = ve->getNCon();

            S.insert(std::make_pair(idx, idx));
            for (int j = 0; j < ncon; ++j) {
                S.insert(std::make_pair(idx, static_cast<int>(ve->nbrIdx(j))));
            }
        }

        pSLUSys.reset(new SLUSystem(S, pMpiComm));
    }

    void advance(double dt) override {
        _advance(pSLUSys.get(), dt);
    }

private:
    std::unique_ptr<SLUSystem> pSLUSys;
    MPI_Comm pMpiComm;
};


}}} // namespace steps::efield::solver

#endif // ndef STEPS_SIM_EFIELD_DVSOLVER_SLU_HPP

// END
