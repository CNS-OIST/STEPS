
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
#include <list>
#include <vector>

// STEPS headers.
#include "math/point.hpp"
#include "mpi/tetvesicle/pointspec.hpp"
#include "mpi/tetvesicle/qtable.hpp"
#include "mpi/tetvesicle/tet_vesraft.hpp"
#include "mpi/tetvesicle/vesicle.hpp"
#include "rng/rng.hpp"
#include "solver/compdef.hpp"
#include "solver/vesicledef.hpp"

// vector_t
#include <fau.de/overlap.hpp>

namespace steps::mpi::tetvesicle {

// Forward declarations.
class Vesicle;
class TetVesicleVesRaft;
class VesUnbind;

////////////////////////////////////////////////////////////////////////////////

class CompVesRaft {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    CompVesRaft(solver::Compdef* compdef, tetmesh::Tetmesh* mesh, TetVesicleVesRaft* vesraft);
    ~CompVesRaft();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////

    void reset();

    ////////////////////////////////////////////////////////////////////////

    solver::Compdef* def() const noexcept {
        return pCompdef;
    }

    // THis is the baseline volume, the volume when no vesicles are present
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
    /// been added to this CompVesRaft object before (i.e. no duplicate checking).
    ///
    void addTet(TetVesRaft* tet);

    inline TetVesRaft* tet(tetrahedron_local_id tet_lidx) const noexcept {
        return pTets[tet_lidx.get()];
    }

    tetrahedron_global_id tetidx_L_to_G(tetrahedron_local_id lidx) const;
    tetrahedron_local_id tetidx_G_to_L(tetrahedron_global_id gidx) const;

    inline uint countTets() const noexcept {
        return pTets.size();
    }

    TetVesRaft* pickTetByVol(double rand01) const;

    const TetVesRaftPVec& tets() const noexcept {
        return pTets;
    }

    ///////////////////////////// VESICLES ////////////////////////////////

    // Setup the vesicle structures
    void setupVesicles();

    // The main event: simulate vesicle diffusion and surface diffusion on
    // vesicles for time dt (s)
    void runVesicle(double dt);

    const auto& getAllVesicles() const {
        return pVesicles;
    }
    // Check position for a vesicle to move to or be placed at. There are various
    // argument options depending which functionality is invovled in calling this
    // function.
    // NOTE: tets_overlap_output and jama_vesicle_output can be modified by the
    // function so the data is available to the calling function.
    bool checkPos(vector_t* pos,
                  double diam,
                  solver::vesicle_global_id vesgidx,
                  std::map<tetrahedron_global_id, double>& tets_overlap_output,
                  Vesicle*& jama_vesicle_output,
                  tetrahedron_global_id central_tet_idx = {},
                  Vesicle* ves = nullptr,
                  math::point3d move_vector = {},
                  bool check_permcomps = true);

    // Get a random position and return overlapping tetrahedron
    tetrahedron_global_id getRandPosByTetStaticVols(math::position_abs* new_pos) const;

    // Vesicle can either be added at a given position (intended for endocytosis)
    // or at a random position within the compartment
    // -1 return means not added, otherwise its unique index
    solver::vesicle_individual_id addVesicle(solver::Vesicledef* vesdef,
                                             const math::position_abs& pos,
                                             tetrahedron_global_id tet_gidx = {});

    // Returns map of species global indices to positions in space (is that
    // necessary? Just need the numbers, surely)
    std::map<solver::spec_global_id, std::vector<PointSpec*>> removeOneVesicle(
        Vesicle* ves,
        tetrahedron_global_id tet_gidx);

    // Check if vesicle starts on a path and set up path diffusion stuff if so
    void checkVesiclePath(Vesicle* v, math::position_abs const& pos_ves) const;

    void deleteSingleVesicle(Vesicle* ves);

    // Here the vesicle is not removed from the simulation, merely it moves
    // compartment
    void loseVesicle(Vesicle* ves);

    // Similarly, not created but 'gained' from another compartment
    void gainVesicle(Vesicle* ves);

    void applyExocytosis(Vesicle* vesicle, solver::exocytosis_global_id exo_gidx);

    // Setting the count will completely wipe any existing vesicles of this
    // type already in the system and new ones will be created and positioned
    // randomly
    void setVesicleCount(solver::vesicle_global_id vidx, uint count);
    uint getVesicleCount(solver::vesicle_global_id vidx) const;

    void setVesiclePos(solver::vesicle_global_id vidx,
                       solver::vesicle_individual_id ves_unique_idx,
                       const std::vector<double>& pos,
                       bool force = false);

    // Note: needed for endocytosis, not from an API call
    void addVesicleSpecs(solver::vesicle_global_id vidx,
                         solver::vesicle_individual_id ves_unique_idx,
                         std::map<solver::spec_global_id, int> ves_specs);

    std::vector<solver::vesicle_individual_id> getVesicleIndices(
        solver::vesicle_global_id vidx) const;

    // Now returns a map/dictionary (in Python) of vesicle unique index->species
    // count for all vesicles of this type
    std::map<solver::vesicle_individual_id, uint> getVesicleSurfaceSpecCountMap(
        solver::vesicle_global_id vidx,
        solver::spec_global_id spec_gidx) const;

    // Returns the summed count
    uint getVesicleSurfaceSpecCount(solver::vesicle_global_id vidx,
                                    solver::spec_global_id spec_gidx) const;
    uint getVesicleInnerSpecCount(solver::vesicle_global_id vidx,
                                  solver::spec_global_id spec_gidx) const;

    void setVesicleTetDcst(solver::vesicle_global_id vidx, tetrahedron_global_id tidx, double dcst);

    uint getVesicleLinkSpecCount(solver::vesicle_global_id vidx,
                                 solver::linkspec_global_id linkspec_gidx) const;

    std::map<solver::vesicle_individual_id, uint> getVesicleLinkSpecCountMap(
        solver::vesicle_global_id vidx,
        solver::linkspec_global_id linkspec_gidx) const;

    // Exocytosis uses this to do a kind of boolean check, so lets return a bool
    bool getVesicleTetOverlap(solver::vesicle_global_id vidx, tetrahedron_global_id tet_gidx) const;

    void addVesiclePermittedComps(solver::vesicle_global_id vidx,
                                  const std::vector<CompVesRaft*>& comps);

    bool checkVesiclePermittedComp(solver::vesicle_global_id vidx, CompVesRaft* comp) {
        return pVesicles_permittedcomps[vidx].find(comp) != pVesicles_permittedcomps[vidx].end();
    }

    inline TetVesicleVesRaft* solverVesRaft() const noexcept {
        return pVesRaft;
    }
    ///////////////////////////////////////////////////////////////////////

  private:
    solver::Compdef* pCompdef;
    double pVol;

    tetmesh::Tetmesh* pMesh;

    TetVesRaftPVec pTets;

    std::map<tetrahedron_global_id, tetrahedron_local_id> pTetidcs_G_to_L;
    std::map<tetrahedron_local_id, tetrahedron_global_id> pTetidcs_L_to_G;

    // The vesicles are stored as a map of vesicles type by index to a
    // list of vesicles of that type. List should be more
    // efficient than a vector for insertions and removals
    std::map<solver::vesicle_global_id, std::list<Vesicle*>> pVesicles;

    // The compartments that vesicles may diffuse freely between
    std::map<solver::vesicle_global_id, std::set<CompVesRaft*>> pVesicles_permittedcomps;

    // User can optionally alter dcst per tet from default for vesicles
    std::map<solver::vesicle_global_id, std::map<tetrahedron_global_id, double>> pVes_Tetskcst;

    // Store a pointer to the RNG for convenience
    const rng::RNGptr pRNG;

    // Just easier to store the solver for all the vesicle KProc stuff
    TetVesicleVesRaft* pVesRaft{nullptr};
};

}  // namespace steps::mpi::tetvesicle
