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

// Standard library & STL headers.
#include <iosfwd>
#include <map>
#include <set>
#include <string>
#include <vector>

// STEPS headers.
#include "rng/rng.hpp"
#include "solver/raftdef.hpp"

namespace steps::mpi::tetvesicle {

// Forward declarations.
class Patch;
class TriRDEF;
class RaftEndocytosis;

////////////////////////////////////////////////////////////////////////////////

class RaftProxy {
  public:
    RaftProxy(solver::Raftdef* raftdef,
              TriRDEF* central_tri,
              solver::raft_individual_id unique_index);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////

    inline solver::Raftdef* def() const noexcept {
        return pDef;
    }

    inline double getDcst() const noexcept {
        return pDef->dcst();
    }

    inline double getDiam() const noexcept {
        return pDef->diameter();
    }

    inline TriRDEF* tri() const noexcept {
        return pTri;
    }

    inline solver::raft_individual_id getUniqueIndex() const noexcept {
        return pIndex;
    }

    /// Return pointer to species' counts on the raft.
    const auto& pools_global() const noexcept {
        return pPoolCount;
    }

    ////////////////////////////////////////////////////////////////////////

    // Set count of species on the raft
    // Being explicit about index type now because pools is going to hold all
    // species for possible transport, whereas reactions use local indices as
    // usual
    void setSpecCountByLidx(solver::spec_local_id slidx, uint count);
    void setSpecCountByGidx(solver::spec_global_id sgidx, uint count);

    uint getSpecCountByLidx(solver::spec_local_id slidx);
    uint getSpecCountByGidx(solver::spec_global_id sgidx);

    // Create on fly for use by MPI structs.
    std::map<index_t, uint> getSpecs();

    ////////////////////////////////////////////////////////////////////////

    // Argument is random number uniform on 0-1
    // Tri* selectDirectionTri(double unf);

    // Update mobility, which can be 0 (free-moving) or non-zero (fixed in place)
    void updImmobility(int mob_upd);

    inline int getImmobilityUpdate() const noexcept {
        return pImmobilityUpdate;
    }

    bool getRaftSReacActive(solver::raftsreac_global_id rsreacidx) const;
    void setRaftSReacInActive(solver::raftsreac_global_id rsreacidx);

    ////////////////////////////////////////////////////////////////////////

  private:
    solver::Raftdef* pDef;

    solver::raft_global_id pRaftIndex;

    solver::raft_individual_id pIndex;

    // math::point3d						pPos;
    TriRDEF* pTri;

    // Table of the populations of the species on the raft.
    std::vector<uint> pPoolCount;

    int pImmobilityUpdate;

    std::set<solver::raftsreac_global_id> pRaftSReac_inactive;
};

}  // namespace steps::mpi::tetvesicle
