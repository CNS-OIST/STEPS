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

#pragma once

// Standard library & STL headers.
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <vector>

// STEPS headers.
#include "math/point.hpp"
#include "mpi/tetvesicle/comp_vesraft.hpp"
#include "mpi/tetvesicle/path.hpp"
#include "mpi/tetvesicle/pointspec.hpp"
#include "mpi/tetvesicle/qtable.hpp"
#include "rng/rng.hpp"
#include "solver/fwd.hpp"
#include "solver/vesicledef.hpp"
#include "solver/vessdiffdef.hpp"

#include <fau.de/overlap.hpp>

namespace steps::mpi::tetvesicle {

// Forward declaration
class Exocytosis;
class LinkSpec;
class Path;

////////////////////////////////////////////////////////////////////////////////

class Vesicle {
  public:
    Vesicle(solver::Vesicledef* vesdef,
            CompVesRaft* comp,
            const math::position_abs& pos,
            solver::vesicle_individual_id unique_index,
            const std::map<tetrahedron_global_id, double>& overlap,
            const double diam,
            const double dcst);
    Vesicle(solver::Vesicledef* vesdef,
            CompVesRaft* comp,
            solver::vesicle_individual_id unique_index,
            std::fstream& cp_file);

    ~Vesicle();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore is done by 2nd constructor

    ////////////////////////////////////////////////////////////////////////
    // GENERAL INFORMATION
    ////////////////////////////////////////////////////////////////////////

    // No local indices for vesicles
    inline solver::vesicle_global_id idx() const noexcept {
        return def()->gidx();
    }

    // Each individual vesicle of this type's unique ID throughout the simulation
    inline solver::vesicle_individual_id getUniqueIndex() const noexcept {
        return pIndex;
    }

    inline solver::Vesicledef* def() const noexcept {
        return pDef;
    }

    inline CompVesRaft* comp() const noexcept {
        return pComp_central;
    }

    inline double getDcst() const noexcept {
        return pDcst;
    }

    inline void setDcst(double dcst) noexcept {
        pDcst = dcst;
    }

    inline double getDiam() const noexcept {
        return pDiam;
    }

    TetVesicleVesRaft* solver() const noexcept;

    const rng::RNGptr rng() const noexcept;

    tetmesh::Tetmesh* mesh() const noexcept;

    ////////////////////////////////////////////////////////////////////////
    // POSITION AND MOBILITY
    ////////////////////////////////////////////////////////////////////////

    // Will return false if it can't move because of spec overlap issue
    bool setPosition(const overlap::Vector& new_pos,
                     std::map<tetrahedron_global_id, double> const& tets_overlap_temp,
                     bool check_comps = false);

    inline const math::position_abs& getPosition() const noexcept {
        return pPos;
    }

    void updatePathBindingRates();

    // Update mobility, which can be 0 (free-moving) or non-zero (fixed in place)
    void updImmobility(int mob_upd);

    inline void setImmobility(uint immob) noexcept {
        pImmobility = immob;
    }

    inline uint getImmobility() const noexcept {
        return pImmobility;
    }

    inline void resetImmobility() noexcept {
        // used by kiss-n-run to release a vesicle that doesn't undergo full exocytosis
        pImmobility = 0;
    }

    ////////////////////////////////////////////////////////////////////////
    // TETRAHEDRONS, OVERLAP
    ////////////////////////////////////////////////////////////////////////

    inline tetrahedron_global_id getCentralTet() const noexcept {
        return pCentral_tet;
    }

    void setOverlap(const std::map<tetrahedron_global_id, double>& overlap);

    inline std::map<tetrahedron_global_id, double> const& getOverlap_gidx() const noexcept {
        return pTets_overlap_gidx;
    }

    std::vector<tetrahedron_global_id> const& getOverlapVec_gidx() noexcept {
        return pTets_overlap_vec_gidx;
    }

    inline bool overlapsTet_gidx(tetrahedron_global_id tet_gidx) const noexcept {
        return pTets_overlap_gidx.find(tet_gidx) != pTets_overlap_gidx.end();
    }

    // Very important function, widely used. Return which (of any) of
    // tets_overlap_temp overlap the point spec_x, spec_y, spec_z (intended to be
    // a species location)
    tetrahedron_global_id getTetSpecOverlap(
        const math::position_abs& pos,
        std::map<tetrahedron_global_id, double> const& tets_overlap_temp);

    ////////////////////////////////////////////////////////////////////////
    // SURFACE SPECIES AND LINK SPECIES
    ////////////////////////////////////////////////////////////////////////

    // Diffuse species and link species on the vesicle surface
    void doSurfaceDiffusion();

    void _recalcQtable_spec(solver::spec_global_id spec_gidx, double d);
    void _recalcQtable_linkspec(solver::linkspec_global_id linkspec_gidx);
    void recalcQtables_();

    double getQPhiSpec_(solver::spec_global_id spec_gidx) const noexcept;
    double getQPhiLinkspec_(solver::linkspec_global_id linkspec_gidx) const noexcept;

    /////////////////// SURFACE SPECIES //////////////////////////

    // Add (or remove if the number is lower) one type of species
    void setSurfSpecCount(solver::spec_global_id spec_gidx, uint count);

    inline uint getSurfSpecCount(solver::spec_global_id spec_gidx) const noexcept {
        const auto it = pSurfSpecs.find(spec_gidx);
        if (it != pSurfSpecs.end()) {
            return it->second.size();
        } else {
            return 0;
        }
    }

    // Return the count of species in a specific tetrahedron
    uint getSurfSpecCount(solver::spec_global_id spec_gidx, tetrahedron_global_id tet_gidx);

    // Can remove specs now for e.g. vesicle reactions
    // tet_gidx argument allows trying to do the update for a specific tet.
    // Currently only feasible if we're removing (negative numbers in specs arg)
    void addSurfSpecs(const std::map<solver::spec_global_id, int>& specs,
                      tetrahedron_global_id tet_gidx = {});

    // add a specific surface spec at a particular position, absolute
    void addOneSurfSpec(solver::spec_global_id spec_gidx,
                        solver::pointspec_individual_id spec_idx,
                        tetrahedron_global_id tet_gidx,
                        const math::position_abs& pos_abs);

    // add one surface spec at a particular position, absolute
    void addOneSurfSpec(solver::spec_global_id spec_gidx,
                        tetrahedron_global_id tet_gidx,
                        const math::position_abs& pos_abs);

    // add one surface spec at a particular position, relative
    void addOneSurfSpec(solver::spec_global_id spec_gidx, math::position_rel_to_ves pos_rel);

    void changeSurfSpecGidx(solver::spec_global_id spec_src_gidx,
                            solver::spec_global_id spec_dst_gidx);

    inline std::map<solver::spec_global_id, std::vector<PointSpec*>>& getSurfSpecs() noexcept {
        return pSurfSpecs;
    }

    std::vector<solver::pointspec_individual_id> getSurfaceSpecIndices(
        solver::spec_global_id spec_gidx) const;

    // Bit drastic, but this is the only way I can reasonably think to get this
    // working for vesMPI first version
    void clearSurfSpecs();


    // Get the positions of a particular species in spherical coordinates
    std::vector<std::vector<double>> getSurfaceSpecPosSpherical(
        solver::spec_global_id spec_gidx) const;

    // Set the position of a particular species in spherical coordinates
    void setSurfaceSpecPosSpherical(solver::spec_global_id spec_gidx,
                                    const std::vector<std::vector<double>>& pos_spherical);

    // Returns position of a surface species. Not relative to vesicle, absolute
    std::vector<std::vector<double>> getSurfaceSpecPos(solver::spec_global_id spec_gidx) const;

    // For a specific tetrahedron
    // Not relative to vesicle, absolute
    std::vector<std::vector<double>> getSurfaceSpecPos(solver::spec_global_id spec_gidx,
                                                       tetrahedron_global_id tet_gidx) const;

    /////////////////// INNER SPECIES //////////////////////////

    // Set count of species inside the vesicle. This function is needed for API
    void setInnerSpecCount(solver::spec_global_id spec_gidx, uint count);

    // Increase the inner count. Used internally, from updates from RDEF
    void incInnerSpecCount(solver::spec_global_id spec_gidx, uint count);

    // Reduce the inner count. Used internally, potentially from kiss-and-run partial release
    void redInnerSpecCounts(const std::map<solver::spec_global_id, uint>& specs);

    uint getInnerSpecCount(solver::spec_global_id spec_gidx) const;

    std::map<solver::spec_global_id, uint> const& getInnerSpecCounts() const noexcept {
        return pInnerSpecCount;
    }

    inline void clearInnerSpecCounts() noexcept {
        pInnerSpecCount.clear();
    }

    /////////////////// LINK SPECS //////////////////////////

    inline std::map<solver::linkspec_individual_id, LinkSpec*> const& getLinkSpecs()
        const noexcept {
        return pLinkSpecs;
    }

    // Add one surface link spec
    void addLinkSpec(solver::linkspec_individual_id linkspec_uniqueid, LinkSpec* link_spec);

    // Just for a bit of extra safety in the debug phase, include the expected
    // LinkSpec too
    void remLinkSpec(solver::linkspec_individual_id linkspec_uniqueid, LinkSpec* linkspec);


    LinkSpec* getLinkSpec(solver::linkspec_individual_id linkspec_uniqueid) const;

    uint getLinkSpecCount(solver::linkspec_global_id linkspec_gidx) const;

    std::vector<solver::linkspec_individual_id> getLinkSpecIndices(
        solver::linkspec_global_id linkspec_gidx) const;

    void updateLinkSpec(solver::linkspec_individual_id linkspec_uniqueid,
                        solver::LinkSpecdef* linkspec_def);


    inline bool containsLink() const noexcept {
        return !pLinkSpecs.empty();
    }

    // This represent a link spec changing identity due to a surface reaction,
    // but the overall count of link species cannot change, hence just changing
    // the index bool changeLinkSpecGidx(solver::linkspec_global_id
    // linkspecbefore_gidx, solver::linkspec_global_id linkspecafter_gidx,
    // tetrahedron_global_id tet_gidx); NOTE if this is still necessary, just
    // change the ID in the object itself, no need to create a new one. This
    // would mean storing the global_index in the LinkSpec object and removing
    // any reference to solver::LinkSpecdef (the only thing LinkSpecdef is used
    // for is to fetch the global ID)

    // Returns vector of all linkspecs with of this index
    std::vector<std::vector<double>> getLinkSpecPos(solver::linkspec_global_id linkspec_gidx) const;

    // See if the vesicle can move or not depending on linkspec length
    bool linkSpecMoveAllowed(const math::point3d& move_vector);

    ////////////////////////////////////////////////////////////////////////
    // PATHS (simple virtual actin etc)
    ////////////////////////////////////////////////////////////////////////

    void checkPathBinding(double dt);

    void checkPathUnbinding(double dt);

    inline bool onPath() const noexcept {
        return pOnPath != nullptr;
    }

    inline const std::shared_ptr<Path>& getPath() const {
        return pOnPath;
    }

    std::pair<std::vector<std::pair<double, math::position_abs>>::const_iterator,
              std::vector<std::pair<double, math::position_abs>>::const_iterator>
    getNextPositions(double dt) noexcept;

    void updatePositionOnPath(
        std::vector<std::pair<double, math::position_abs>>::const_iterator end);

    std::pair<std::shared_ptr<Path>, math::position_abs> getCurrentPathPosition() const;

    void removeFromPath();

    ////////////////////////////////////////////////////////////////////////

    inline void scheduleExocytosis(solver::exocytosis_global_id exo_gidx) {
        pAppliedExocytosis.insert(exo_gidx);
    }

    solver::exocytosis_global_id appliedExocytosis();

    inline void clearExocytosis() {
        pAppliedExocytosis.clear();
    };


  private:
    // Stores the overlap tetrahedrons by index and volume of overlap
    solver::Vesicledef* pDef;
    CompVesRaft* pComp_central;

    solver::vesicle_individual_id pIndex;
    math::position_abs pPos;

    double pDiam;
    double pDcst;

    tetrahedron_global_id pCentral_tet;

    std::map<tetrahedron_global_id, double> pTets_overlap_gidx;

    std::vector<tetrahedron_global_id> pTets_overlap_vec_gidx;

    std::map<solver::spec_global_id, std::vector<PointSpec*>> pSurfSpecs;

    std::map<solver::linkspec_individual_id, LinkSpec*> pLinkSpecs;

    // Map of the populations of the species inside the vesicle.
    std::map<solver::spec_global_id, uint> pInnerSpecCount;

    uint pImmobility;

    std::map<solver::spec_global_id, std::shared_ptr<Qtable>> pQtables_spec;
    std::map<solver::linkspec_global_id, std::shared_ptr<Qtable>> pQtables_linkspec;

    // PATHS
    util::strongid_vector<solver::path_global_id, double> pPathBindingRates;  // Binding rate to
                                                                              // each path
    double pTotalPathBindingRate;

    std::vector<std::pair<double, math::position_abs>> pPathPositions;
    std::vector<std::pair<double, math::position_abs>>::const_iterator pPath_curr_pos;
    std::vector<std::pair<double, math::position_abs>>::const_iterator
        pPath_next_pos_end;  // End iterator for the possible next positions on the path. This
                             // allows several positions to be tested, starting with the furthest.
    std::shared_ptr<Path> pOnPath;
    math::position_abs pPath_starting_shift;
    double pTime_accum;       // the accumulated time since the last changed position
    double pTime_accum_next;  // 'next' is stored because move may be invalid

    // Possible multiple exocytosis events will fire during one step. Store
    // as a list and figure out how to select one if this happens
    std::set<solver::exocytosis_global_id> pAppliedExocytosis;
};

inline bool operator<(const Vesicle& lhs, const Vesicle& rhs) {
    return lhs.getUniqueIndex() < rhs.getUniqueIndex();
}

}  // namespace steps::mpi::tetvesicle
