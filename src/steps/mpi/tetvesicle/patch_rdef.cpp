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

#include "mpi/tetvesicle/patch_rdef.hpp"

// STEPS headers.
#include "math/point.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/tri_rdef.hpp"
#include "solver/compdef.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"
#include "util/distribute.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::mpi::tetvesicle {

////////////////////////////////////////////////////////////////////////////////

PatchRDEF::PatchRDEF(solver::Patchdef* patchdef, tetmesh::Tetmesh* mesh, TetVesicleRDEF* rdef)
    : pPatchdef(patchdef)
    , pArea(0.0)
    , pMesh(mesh)
    , pRDEF(rdef) {
    AssertLog(pPatchdef != nullptr);
    pRNG = pPatchdef->statedef()->rng();
}

////////////////////////////////////////////////////////////////////////////////

PatchRDEF::~PatchRDEF() = default;

////////////////////////////////////////////////////////////////////////////////

void PatchRDEF::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pArea);
}

////////////////////////////////////////////////////////////////////////////////

void PatchRDEF::restore(std::fstream& cp_file) {
    util::compare(cp_file, pArea);
}

////////////////////////////////////////////////////////////////////////////////

void PatchRDEF::reset() const {
    def()->reset();
}

////////////////////////////////////////////////////////////////////////////////

void PatchRDEF::addTri(TriRDEF* tri) {
    AssertLog(tri->patchdef() == def());

    index_t lidx_uint = pTris.size();
    triangle_local_id lidx(lidx_uint);
    triangle_global_id gidx(tri->idx());

    pTriidcs_L_to_G.emplace(lidx, gidx);
    pTriidcs_G_to_L.emplace(gidx, lidx);

    pTris.push_back(tri);
    pArea += tri->area();

    tri->setPatchRDEF(this);
}

}  // namespace steps::mpi::tetvesicle
