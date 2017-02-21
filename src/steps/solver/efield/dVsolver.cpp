/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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


// STL headers.
#include <algorithm>
#include <utility>
#include <cmath>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/efield/dVsolver.hpp"
#include "steps/solver/efield/tetmesh.hpp"

namespace steps {
namespace solver {
namespace efield {

void dVSolverBase::initMesh(TetMesh *mesh) {
    pMesh = mesh;
    pNVerts = pMesh->countVertices();
    pNTris = pMesh->getNTri();

    pV.assign(pNVerts, 0.0);
    pGExt.assign(pNVerts, 0.0);
    pVertexClamp.assign(pNVerts, false);
    pVertCur.assign(pNVerts, 0.0);
    pVertCurClamp.assign(pNVerts, 0.0);

    pTriCur.assign(pNTris, 0.0);
    pTriCurClamp.assign(pNTris, 0.0);
}

void dVSolverBase::setSurfaceConductance(double g_surface, double v_rev) {
    pVExt = v_rev;
    if (!pMesh) return;

    for (int i = 0; i < pNVerts; ++i) {
        VertexElement* ve = pMesh->getVertex(i);
        pGExt[ve->getIDX()] = g_surface * ve->getSurfaceArea();
    }
}

int dVSolverBase::meshHalfBW(TetMesh *mesh) {
    int halfbw = 0;
    int nVerts = mesh->countVertices();
    for (int i = 0; i < nVerts; ++i) {
        VertexElement *ve = mesh->getVertex(i);

        int idx = ve->getIDX();
        int ncon = ve->getNCon();
        for (int j = 0; j < ncon; ++j) {
            halfbw = std::max(halfbw, std::abs(idx - (int)ve->nbrIdx(j)));
        }
    }

    return halfbw;
}

}}} // namespace steps::solver::efield

