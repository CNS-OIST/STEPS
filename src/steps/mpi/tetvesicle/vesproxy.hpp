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
#include <iostream>
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include "mpi/tetvesicle/comp_rdef.hpp"
#include "mpi/tetvesicle/tet_rdef.hpp"
#include "rng/rng.hpp"
#include "solver/vesicledef.hpp"
#include "util/common.hpp"

namespace steps::mpi::tetvesicle {

// This object stores vesicle information per tet between vesicle application.
// This information will be used by local SSA during multiple opsplit taus, and
// written back to parent vesicle object (for this tet) before next vesicle
// application.
class VesProxy {
  public:
    VesProxy(solver::Vesicledef* vesdef,
             TetRDEF* tet,
             solver::vesicle_individual_id unique_index,
             const math::position_abs& pos,
             bool contains_link);
    ~VesProxy();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // GENERAL INFORMATION
    ////////////////////////////////////////////////////////////////////////
    // No local indices for vesicles
    inline solver::vesicle_global_id idx() const noexcept {
        return def()->gidx();
    }

    inline solver::Vesicledef* def() const noexcept {
        return pDef;
    }

    inline double getDiam() const noexcept {
        return def()->diameter();
    }

    inline solver::vesicle_individual_id getUniqueIndex() const noexcept {
        return pUniqueIndex;
    }

    inline const rng::RNGptr& rng() const noexcept;

    tetmesh::Tetmesh* mesh() const noexcept;

    ////////////////////////////////////////////////////////////////////////
    // MOBILITY
    ////////////////////////////////////////////////////////////////////////

    // Update mobility, which can be 0 (free-moving) or non-zero (fixed in place)
    void updImmobility(int mob_upd);

    inline int getImmobilityUpdate() const noexcept {
        return pImmobilityUpdate;
    }

    ////////////////////////////////////////////////////////////////////////
    // TETRAHEDRONS, OVERLAP
    ////////////////////////////////////////////////////////////////////////

    inline tetrahedron_global_id getTetGidx() const noexcept {
        return pTet->idx();
    }

    inline TetRDEF* getTet() const noexcept {
        return pTet;
    }

    /////////////////// SURFACE SPECIES //////////////////////////

    // void setSpecCount_V(solver::spec_local_id spec_lidx, uint count);

    uint getSpecCount_V(solver::spec_global_id spec_gidx) const noexcept {
        auto it = pSpecs_V.find(spec_gidx);
        if (it != pSpecs_V.end()) {
            return it->second.size();
        } else {
            return 0;
        }
    }

    inline const std::map<
        solver::spec_global_id,
        std::vector<std::pair<solver::pointspec_individual_id, math::position_abs>>>&
    getSurfSpecs() const noexcept {
        return pSpecs_V;
    }

    void addSurfSpec(solver::spec_global_id spec_gidx,
                     solver::pointspec_individual_id idx,
                     const math::position_abs&);

    bool addSurfSpecs(const std::map<solver::spec_global_id, int>& specs);

    const std::vector<std::pair<solver::pointspec_individual_id, math::position_abs>>&
    getSurfaceSpecPos_absolute(solver::spec_global_id spec_gidx) const;

    // Remove one surface spec at specific array index, return position relative
    // to vesicle centre
    void removeOneSurfSpec(solver::spec_global_id spec_gidx, uint array_index);

    /////////////////// SURFACE LINKSPECIES /////////////////////


    uint getLinkSpecCount_V(solver::linkspec_global_id linkspec_gidx) const;

    void addLinkSpec(solver::linkspec_individual_id ls_unique_idx,
                     solver::linkspec_global_id ls_gidx,
                     math::position_abs ls_pos_abs);

    void changeLinkSpecGidx(solver::linkspec_global_id linkspecbefore_gidx,
                            solver::linkspec_global_id linkspecafter_gidx);

    const std::map<solver::linkspec_individual_id,
                   std::pair<solver::linkspec_global_id, math::position_abs>>&
    getLinkSpecs() const noexcept {
        return pLinkSpecs_V;
    }

    solver::linkspec_global_id getLinkSpecGidx(solver::linkspec_individual_id ls_unique_id) {
        return pLinkSpecs_V[ls_unique_id].first;
    }

    // A set of LinkSpecs that have some update this run on RDEF
    inline void reqLinkSpecUpd(solver::linkspec_individual_id linkspec_unique_idx) {
        pLinkSpecsUpd.insert(linkspec_unique_idx);
    }
    inline std::set<solver::linkspec_individual_id> const& getLinkSpecUpd() const noexcept {
        return pLinkSpecsUpd;
    }
    inline void clearLinkSpecUpd() {
        pLinkSpecsUpd.clear();
    }

    // Returns true if vesproxy has a linkspec of any flavour
    inline bool hasLinkSpec() {
        return pContainsLink || !pLinkSpecs_V.empty();
    }

    /////////////////// INNER SPECIES //////////////////////////

    void incSpecCount_I(solver::spec_global_id spec_gidx, uint count);
    // uint getSpecCount_I(solver::spec_global_id spec_gidx) const;

    const std::map<solver::spec_global_id, uint>& getSpecCounts_I() const noexcept {
        return pSpecCount_I;
    }

    ////////////////////////////////////////////////////////////////////////

    inline solver::exocytosis_global_id exoApplied() const noexcept {
        return pExoApplied;
    }

    void applyExo(solver::exocytosis_global_id exo_gidx);


    ////////////////////////////////////////////////////////////////////////

  private:
    solver::Vesicledef* pDef;

    solver::vesicle_global_id pVesicleIndex;

    solver::vesicle_individual_id pUniqueIndex;

    TetRDEF* pTet;

    // Central position of parent vesicle
    math::position_abs pPos;

    // uint                                               * pSpecCount_V;

    // Map of the populations of the species on the surface of the vesicle to
    // position, absolute
    std::map<solver::spec_global_id,
             std::vector<std::pair<solver::pointspec_individual_id, math::position_abs>>>
        pSpecs_V;
    std::vector<std::pair<solver::pointspec_individual_id, math::position_abs>> pEmptyPosVec;

    // Store all known positions to be used for adding any new species,
    // to avoid potentially long try and reject loops
    std::vector<math::position_abs> pAllPositions;

    // Table of the populations of the linkspecies on the surface of the vesicle.
    // uint                                               * pLinkSpecCount_V;

    // let's try it like this. Unique IDs shouldn't change because Comp needs to
    // track them but global ID, the type of link species, can. So makes sense to
    // store them like this:
    std::map<solver::linkspec_individual_id,
             std::pair<solver::linkspec_global_id, math::position_abs>>
        pLinkSpecs_V;

    // This is a list of LinkSpecs that changed their ID, for optimality
    std::set<solver::linkspec_individual_id> pLinkSpecsUpd;

    // This flag comes from VesRaft and will be true if a linkspec is present anywhere on the parent
    // vesicle (not necessarily the part that this vesproxy represents)
    bool pContainsLink;

    // Table of the populations of the species inside the vesicle.
    std::map<solver::spec_global_id, uint> pSpecCount_I;

    int pImmobilityUpdate;

    solver::exocytosis_global_id pExoApplied;
};

}  // namespace steps::mpi::tetvesicle
