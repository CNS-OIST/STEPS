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
#include <fstream>
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include "rng/rng.hpp"
#include "solver/raftdef.hpp"
#include "util/common.hpp"

namespace steps::mpi::tetvesicle {

// Forward declarations.
class PatchVesRaft;
class TriVesRaft;
class RaftEndocytosis;
class RaftDis;

////////////////////////////////////////////////////////////////////////////////

class Raft {
  public:
    Raft(solver::Raftdef* raftdef,
         PatchVesRaft* patch,
         TriVesRaft* central_tri,
         math::position_abs& pos,
         solver::raft_individual_id unique_index);

    Raft(solver::Raftdef* raftdef, PatchVesRaft* patch, std::fstream& cp_file);

    ~Raft();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////

    void checkpoint(std::fstream& cp_file);
    // NOTE: restore is done by second constructor

    ////////////////////////////////////////////////////////////////////////////////

    inline solver::Raftdef* def() const noexcept {
        return pDef;
    }

    inline solver::raft_global_id idx() const noexcept {
        return def()->gidx();
    }

    inline PatchVesRaft* patch() const noexcept {
        return pPatch;
    }

    inline double getDcst() const noexcept {
        return pDef->dcst();
    }

    double getScaledDcst() const noexcept {
        return pScaledDcst;
    }

    inline double getDiam() const noexcept {
        return pDef->diameter();
    }

    inline TriVesRaft* tri_central() const noexcept {
        return pTri_central;
    }

    std::vector<triangle_global_id> const& getOverlapVec() const noexcept {
        return pTris_overlap_vec;
    }

    inline solver::raft_individual_id getUniqueIndex() const noexcept {
        return pIndex;
    }

    /// Return pointer to species' counts on the raft.
    inline const auto& pools_global() const noexcept {
        return pPoolCount;
    }

    ////////////////////////////////////////////////////////////////////////

    // TODO check position is tested by whatever calls this function
    void setPosition(math::position_abs& new_pos, TriVesRaft* centraltri);
    inline math::position_abs const& getPosition() const noexcept {
        return pPos;
    }

    ////////////////////////////////////////////////////////////////////////

    // Set count of species on the raft
    // Being explicit about index type now because pools is going to hold all
    // species for possible transport, whereas reactions use local indices as
    // usual
    void setSpecCountByLidx(solver::spec_local_id slidx, uint count);
    void setSpecCountByGidx(solver::spec_global_id sgidx, uint count);

    // When reading from RDEF multiple Raft Proxies can have update for same
    // species, so need the inc function here
    void incSpecCountByGidx(solver::spec_global_id sgidx, uint count);

    uint getSpecCountByLidx(solver::spec_local_id slidx);
    uint getSpecCountByGidx(solver::spec_global_id sgidx);

    std::map<triangle_global_id, uint> getTriSpecCounts(solver::spec_global_id,
                                                        const rng::RNGptr rng);

    // Bit drastic, but this is the only way I can reasonably think to get this
    // working for vesMPI first version
    inline void clearSpecs() {
        std::fill(pPoolCount.begin(), pPoolCount.end(), 0);
    }

    ////////////////////////////////////////////////////////////////////////

    // This function still needed to set up special cases RaftDis and
    // RaftEndocytosis
    void setupKProcs();

    // Now endocytosis and dis stored separately, different routines for parallel
    // solver
    inline auto& raftdiss() noexcept {
        return RaftDiss;
    }

    inline auto& raftendos() noexcept {
        return RaftEndocytosiss;
    }

    // Now to be called by patch during vesicle step. Any one of these kprocs
    // may actually result in loss of this Raft, so each individual
    // applciation needs to be checked
    bool applyEndoAndDis(double raft_dt);

    // This is necessary because endos can fail, though there is no bool returned
    // from apply function
    inline void appliedEndo() noexcept {
        pAppliedEndo = true;
    }

    void setRaftSReacActive(solver::raftsreac_global_id rsridx, bool active);

    bool getRaftSReacActive(solver::raftsreac_global_id rsridx) const;

    inline const auto& raftsreacsinactive() const noexcept {
        return pRaftSReac_inactive;
    }

    ////////////////////////////////////////////////////////////////////////

    // Argument is random number uniform on 0-1
    TriVesRaft* selectDirectionTri(double unf);

    // Update mobility, which can be 0 (free-moving) or non-zero (fixed in place)
    void updImmobility(int mob_upd);
    inline uint getImmobility() const noexcept {
        return pImmobility;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    util::strongid_vector<solver::raftdis_local_id, RaftDis*> RaftDiss;
    util::strongid_vector<solver::raftendocytosis_local_id, RaftEndocytosis*> RaftEndocytosiss;

    ////////////////////////////////////////////////////////////////////////

    solver::Raftdef* pDef;
    PatchVesRaft* pPatch;
    solver::raft_individual_id pIndex;

    math::position_abs pPos;
    TriVesRaft* pTri_central;
    // Stores the overlap tris
    std::vector<triangle_global_id> pTris_overlap_vec;

    // Table of the populations of the species on the raft.
    util::strongid_vector<solver::spec_global_id, uint> pPoolCount;

    /// Properly scaled diffusivity constant.
    double pScaledDcst;

    /// Used in selecting which directory the molecule should go.
    std::array<double, 2> pCDFSelector{0.0, 0.0};

    bool pAppliedEndo;

    uint pImmobility;

    std::set<solver::raftsreac_global_id> pRaftSReac_inactive;
};

inline bool operator<(const Raft& lhs, const Raft& rhs) {
    return lhs.getUniqueIndex() < rhs.getUniqueIndex();
}

}  // namespace steps::mpi::tetvesicle
