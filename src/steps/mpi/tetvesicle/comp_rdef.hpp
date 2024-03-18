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

#pragma once

// STL headers.
#include <fstream>
#include <vector>

// STEPS headers.
#include "mpi/tetvesicle/tet_rdef.hpp"
#include "rng/rng.hpp"
#include "solver/compdef.hpp"

namespace steps::mpi::tetvesicle {

// Forward declarations.
class TetVesicleRDEF;

class CompRDEF {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    CompRDEF(solver::Compdef* compdef, tetmesh::Tetmesh* mesh, TetVesicleRDEF* rdef);
    ~CompRDEF();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////

    void reset() const;

    ////////////////////////////////////////////////////////////////////////

    solver::Compdef* def() const noexcept {
        return pCompdef;
    }

    // This is the baseline volume, the volume when no vesicles are present
    inline double vol() const noexcept {
        return pVol;
    }

    inline tetmesh::Tetmesh* mesh() const noexcept {
        return pMesh;
    }

    /// Return the random number generator object.
    inline const rng::RNGptr& rng() const noexcept {
        return pRNG;
    }

    ///////////////////////////// TETS /////////////////////////////////////

    /// Checks whether the Tet's compdef() corresponds to this object's
    /// CompDef. There is no check whether the Tet object has already
    /// been added to this CompRDEF object before (i.e. no duplicate checking).
    ///
    void addTet(TetRDEF* tet);

    inline TetRDEF* tet(tetrahedron_local_id tet_lidx) const noexcept {
        return pTets[tet_lidx.get()];
    }

    tetrahedron_global_id tetidx_L_to_G(tetrahedron_local_id lidx) const;
    tetrahedron_local_id tetidx_G_to_L(tetrahedron_global_id gidx) const;

    inline uint countTets() const noexcept {
        return pTets.size();
    }

    TetRDEF* pickTetByVol(double rand01) const;

    const TetRDEFPVec& tets() const noexcept {
        return pTets;
    }

    ///////////////////////////////////////////////////////////////////////

  private:
    solver::Compdef* pCompdef;
    double pVol;

    tetmesh::Tetmesh* pMesh;

    TetRDEFPVec pTets;

    std::map<tetrahedron_global_id, tetrahedron_local_id> pTetidcs_G_to_L;
    std::map<tetrahedron_local_id, tetrahedron_global_id> pTetidcs_L_to_G;

    // Store a pointer to the RNG for convenience
    const rng::RNGptr pRNG;

    // Just easier to store the solver for all the vesicle KProc stuff ??
    // TODO, still needed??
    TetVesicleRDEF* pRDEF{nullptr};
};

}  // namespace steps::mpi::tetvesicle
