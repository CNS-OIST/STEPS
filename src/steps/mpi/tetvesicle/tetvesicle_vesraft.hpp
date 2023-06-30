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
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <mpi.h>

// STEPS headers.
#include "geom/tetmesh.hpp"
#include "mpi/tetvesicle/comp_vesraft.hpp"
#include "mpi/tetvesicle/crstruct.hpp"
#include "mpi/tetvesicle/diffboundary.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/linkspecpair.hpp"
#include "mpi/tetvesicle/patch_vesraft.hpp"
#include "mpi/tetvesicle/path.hpp"
#include "mpi/tetvesicle/sdiffboundary.hpp"
#include "mpi/tetvesicle/tet_vesraft.hpp"
#include "mpi/tetvesicle/tetvesicle_common.hpp"
#include "mpi/tetvesicle/tri_vesraft.hpp"
#include "mpi/tetvesicle/vesicle.hpp"
#include "solver/api.hpp"
#include "solver/efield/efield.hpp"
#include "solver/statedef.hpp"
#include "util/common.hpp"
#include "util/vocabulary.hpp"
#include <fau.de/overlap.hpp>

namespace steps::mpi::tetvesicle {

class Vesicle;
class CompVesRaft;
class DiffBoundary;
class SDiffBoundary;

////////////////////////////////////////////////////////////////////////////////

class TetVesicleVesRaft: public solver::API {
  public:
    TetVesicleVesRaft(model::Model* m, wm::Geom* g, const rng::RNGptr& r, int calcMembPot);
    ~TetVesicleVesRaft();

    ////////////////////////////////////////////////////////////////////////
    // SOLVER INFORMATION
    ////////////////////////////////////////////////////////////////////////

    std::string getSolverName() const override;
    std::string getSolverDesc() const override;
    std::string getSolverAuthors() const override;
    std::string getSolverEmail() const override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS: API FUNCTIONS
    ////////////////////////////////////////////////////////////////////////

    void reset() override;
    void run(double endtime) override;

    void advance(double adv) override;
    // void advanceSteps(uint nsteps);
    void step() override;

    void setTime(double time) override;
    double getTime() const override;

    double getA0() const override;

    void setNSteps(uint nsteps) override;
    uint getNSteps() const override;

    double getVesicleDT() const override;
    void setVesicleDT(double dt) override;

    void checkpoint(std::string const& file_name) override;
    void restore(std::string const& file_name) override;

    double getTetReducedVol(tetrahedron_global_id tidx) const;

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    void setEfieldDT(double efdt) override {
        pEFDT = efdt;
    }

    inline double getEfieldDT() const noexcept override {
        return pEFDT;
    }

    void setTemp(double t) override {
        pTemp = t;
    }

    inline double getTemp() const noexcept override {
        return pTemp;
    }

    // Vesicle Paths ///////////////////////////////////////////////////////////

    void createPath(std::string const& id) override;

    void addPathPoint(std::string const& path_name,
                      uint point_id,
                      const std::vector<double>& position) override;

    void addPathBranch(std::string const& path_name,
                       uint point_id,
                       const std::map<uint, double>& dest_points) override;

    std::map<std::string, std::map<uint, std::pair<std::vector<double>, std::map<uint, double>>>>
    getAllPaths() const override;

    // MPI /////////////////////////////////////////////////////////////////////

    /// Set if the outputs are synced across all ranks.
    /// If not, set the output rank.
    void setOutputSync(bool enable_sync, int output_rank = 0);

    /// Get if the outputs are synced across all ranks.
    bool getOutputSyncStatus() const;

    /// Get the rank id for data output.
    int getOutputSyncRank() const;

    // ROI Data Access /////////////////////////////////////////////////////////

    void setROITetSpecClamped(const std::vector<tetrahedron_global_id>& triangles,
                              const std::string& s,
                              bool b);
    void setROITriSpecClamped(const std::vector<triangle_global_id>& triangles,
                              const std::string& s,
                              bool b);

    double getROITetSpecCount(const std::vector<tetrahedron_global_id>& triangles,
                              const std::string& s) const;
    double getROITriSpecCount(const std::vector<triangle_global_id>& triangles,
                              const std::string& s) const;

    void setROITetSpecCount(const std::vector<tetrahedron_global_id>& triangles,
                            const std::string& s,
                            double count);
    void setROITriSpecCount(const std::vector<triangle_global_id>& triangles,
                            const std::string& s,
                            double count);

    /// Get species counts of a list of tetrahedrons
    std::vector<double> getROITetSpecCounts(const std::string& ROI_id,
                                            std::string const& s) const override;

    /// Get species counts of a list of triangles
    std::vector<double> getROITriSpecCounts(const std::string& ROI_id,
                                            std::string const& s) const override;

    /// Get species counts of a list of tetrahedrons
    void getROITetSpecCountsNP(const std::string& ROI_id,
                               std::string const& s,
                               double* counts,
                               size_t output_size) const override;

    /// Get species counts of a list of triangles
    void getROITriSpecCountsNP(const std::string& ROI_id,
                               std::string const& s,
                               double* counts,
                               size_t output_size) const override;

    double getROIVol(const std::string& ROI_id) const override;
    double getROIArea(const std::string& ROI_id) const override;

    double getROISpecCount(const std::string& ROI_id, std::string const& s) const override;
    void setROISpecCount(const std::string& ROI_id, std::string const& s, double count) override;

    double getROISpecAmount(const std::string& ROI_id, std::string const& s) const override;
    void setROISpecAmount(const std::string& ROI_id, std::string const& s, double) override;

    double getROISpecConc(const std::string& ROI_id, std::string const& s) const override;
    void setROISpecConc(const std::string& ROI_id, std::string const& s, double conc) override;

    void setROISpecClamped(const std::string& ROI_id, std::string const& s, bool b) override;

    void setROIReacK(const std::string& ROI_id, std::string const& r, double kf) override;
    void setROISReacK(const std::string& ROI_id, std::string const& sr, double kf) override;
    void setROIDiffD(const std::string& ROI_id, std::string const& d, double dk) override;

    void setROIReacActive(const std::string& ROI_id, std::string const& r, bool a) override;
    void setROISReacActive(const std::string& ROI_id, std::string const& sr, bool a) override;
    void setROIDiffActive(const std::string& ROI_id, std::string const& d, bool act) override;
    void setROIVDepSReacActive(const std::string& ROI_id, std::string const& vsr, bool a) override;

    unsigned long long getROIReacExtent(const std::string& ROI_id,
                                        std::string const& r) const override;
    void resetROIReacExtent(const std::string& ROI_id, std::string const& r) override;

    unsigned long long getROISReacExtent(const std::string& ROI_id,
                                         std::string const& sr) const override;
    void resetROISReacExtent(const std::string& ROI_id, std::string const& sr) override;

    unsigned long long getROIDiffExtent(const std::string& ROI_id,
                                        std::string const& d) const override;
    void resetROIDiffExtent(const std::string& ROI_id, std::string const& s) override;


    // Batch Data Access //////////////////////////////////////////////////////

    std::vector<double> getBatchTetSpecCounts(const std::vector<index_t>& tets,
                                              std::string const& s) const override;

    std::vector<double> getBatchTriSpecCounts(const std::vector<index_t>& tris,
                                              std::string const& s) const override;

    void getBatchTetSpecCountsNP(const index_t* indices,
                                 size_t input_size,
                                 std::string const& s,
                                 double* counts,
                                 size_t output_size) const override;

    void getBatchTriSpecCountsNP(const index_t* indices,
                                 size_t input_size,
                                 std::string const& s,
                                 double* counts,
                                 size_t output_size) const override;


    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      CALLED FROM EXTERNAL OBJECTS
    ////////////////////////////////////////////////////////////////////////

    inline double getVesicleDT_() const {
        return pVesicledt;
    };

    // These functions are intended to be called internally, within the solver methods, and not
    // through the API
    double getTriSpecCount_(triangle_global_id tidx, solver::spec_global_id sidx) const;
    void setTriSpecCount_(triangle_global_id tidx, solver::spec_global_id sidx, double n);

    ////////////////////////////////////////////////////////////////////////
    ////////////// Internal functions, not part of API /////////////////

    // Addition of overlap from creation of a new vesicle
    void addOverlap_(std::map<tetrahedron_global_id, double>& tets_overlap, Vesicle* ves) const;

    // Remove overlap from deletion of a vesicle (probably from exocytosis)
    void removeOverlap_(std::map<tetrahedron_global_id, double>& tets_overlap, Vesicle* ves) const;

    // Get the next unique index for the vesicle
    solver::vesicle_individual_id getVesicleNextIndex_();

    void setVesicleSpecDiffD_(solver::vesicle_global_id vidx,
                              solver::spec_global_id spec_gidx,
                              double d);

    inline void recordVesicle_(solver::vesicle_individual_id ves_unique_index, Vesicle* vesicle) {
        AssertLog(pVesicles.find(ves_unique_index) == pVesicles.end());
        pVesicles[ves_unique_index] = vesicle;
    }

    inline void removeVesicle_(solver::vesicle_individual_id ves_unique_index, Vesicle* vesicle) {
        AssertLog(pVesicles[ves_unique_index] == vesicle);
        pVesicles.erase(ves_unique_index);
    }

    LinkSpec* getLinkSpec_(solver::linkspec_individual_id);

    inline void recordLinkSpec_(solver::linkspec_individual_id ls_unique_index, LinkSpec* ls) {
        AssertLog(pLinkSpecs.find(ls_unique_index) == pLinkSpecs.end());
        pLinkSpecs[ls_unique_index] = ls;
    }

    inline void removeLinkSpec_(solver::linkspec_individual_id ls_unique_index, LinkSpec* ls) {
        AssertLog(pLinkSpecs[ls_unique_index] == ls);
        pLinkSpecs.erase(ls_unique_index);
    }

    void removeLinkSpecPair_(const LinkSpecPair* lsp);

    std::vector<Path*> vesicleCrossedPaths_(const math::position_abs& ves_pos,
                                            solver::vesicle_global_id ves_gidx,
                                            double ves_rad);


    double getQPhiSpec_(solver::vesicle_global_id vidx, solver::spec_global_id spec_gidx);
    double getQPhiLinkspec_(solver::vesicle_global_id vidx,
                            solver::linkspec_global_id linkspec_gidx);

    void _recalcQtable_spec(solver::vesicle_global_id vidx,
                            solver::spec_global_id spec_gidx,
                            double d);
    void _recalcQtable_linkspec(solver::vesicle_global_id vidx,
                                solver::linkspec_global_id linkspec_gidx,
                                double d);

    void _recalcQtables();

    // Get the next unique index for a pointspec
    solver::pointspec_individual_id getPointSpecNextIndex_();

    // Get the next unique index for the raft
    solver::raft_individual_id getRaftNextIndex_();

    inline void recordRaft_(solver::raft_individual_id raft_unique_index, Raft* raft) {
        AssertLog(pRafts.find(raft_unique_index) == pRafts.end());
        pRafts[raft_unique_index] = raft;
    }

    inline void removeRaft_(solver::raft_individual_id raft_unique_index, Raft* raft) {
        AssertLog(pRafts[raft_unique_index] == raft);
        pRafts.erase(raft_unique_index);
    }

    void setSingleRaftSpecCount_(solver::raft_global_id ridx,
                                 solver::raft_individual_id raft_unique_index,
                                 solver::spec_global_id sidx,
                                 uint c);

    inline TetVesRaft* tet_(tetrahedron_global_id tidx) const {
        AssertLog(tidx.get() < pTets.size());
        return pTets[tidx];
    }

    // Tet_ext is for use by third-party overlap library
    inline Tetrahedron* tet_ext_(tetrahedron_global_id tet_gidx) const {
        AssertLog(tet_gidx.get() < pTet_ext.size());
        return pTet_ext[tet_gidx];
    }

    inline TriVesRaft* tri_(triangle_global_id tidx) const {
        AssertLog(tidx.get() < pTris.size());
        return pTris[tidx];
    }


    // register pool count sync for a tet. Needs to be public because called by comp during
    // exocytosis
    void regTetPoolSync_(tetrahedron_global_id tet_gidx,
                         solver::spec_global_id spec_gidx,
                         uint count);

    // register pool count sync for a tri
    void regTriPoolSync_(triangle_global_id tri_gidx, solver::spec_global_id spec_gidx, uint count);

    inline double getMaxWalkDistSqFact_() const {
        return maxWalkDistSqFact;
    }
    void setMaxWalkDistFact(double fact) {
        if (fact <= 1) {
            ArgErrLog(
                "The maximum walking distance factor for the walking algorithm should be higher "
                "than 1.");
        }
        maxWalkDistSqFact = fact * fact;
    }
    inline double getMinNbTetVisited_() const {
        return minNbTetVisited;
    }
    void setMinNbTetVisited(uint n) {
        if (n < 1) {
            ArgErrLog(
                "The minimum number of explored tetrahedrons for the walking algorithm should be "
                "higher or equal to 1.");
        }
        minNbTetVisited = n;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    // Private functions. Either called within API functions or otherwise
    // called only within solver
    ////////////////////////////////////////////////////////////////////////
  private:
    void _runVesicle(double dt);
    void _runRaft(double dt);

    inline tetmesh::Tetmesh* _mesh() const {
        return pMesh;
    }

  public:
    inline CompVesRaft* getComp_(solver::comp_global_id cidx) const {
        // Moved AssertLogions to this access routine to cut code duplication.
        AssertLog(cidx < statedef().countComps());
        AssertLog(statedef().countComps() == pComps.size());

        auto c = pComps[cidx];
        AssertLog(c != nullptr);
        return c;
    }

    inline PatchVesRaft* getPatch_(solver::patch_global_id pidx) const {
        // Moved AssertLogions to this access routine to cut code duplication.
        AssertLog(pidx < statedef().countPatches());
        AssertLog(statedef().countPatches() == pPatches.size());

        auto p = pPatches[pidx];
        AssertLog(p != nullptr);
        return p;
    }

  private:
    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

    double _getCompVol(solver::comp_global_id cidx) const override;

    double _getCompSpecCount(solver::comp_global_id cidx,
                             solver::spec_global_id sidx) const override;
    void _setCompSpecCount(solver::comp_global_id cidx,
                           solver::spec_global_id sidx,
                           double n) override;

    double _getCompSpecAmount(solver::comp_global_id cidx,
                              solver::spec_global_id sidx) const override;
    void _setCompSpecAmount(solver::comp_global_id cidx,
                            solver::spec_global_id sidx,
                            double a) override;

    double _getCompSpecConc(solver::comp_global_id cidx,
                            solver::spec_global_id sidx) const override;
    void _setCompSpecConc(solver::comp_global_id cidx,
                          solver::spec_global_id sidx,
                          double c) override;

    bool _getCompSpecClamped(solver::comp_global_id cidx,
                             solver::spec_global_id sidx) const override;
    void _setCompSpecClamped(solver::comp_global_id cidx,
                             solver::spec_global_id sidx,
                             bool b) override;

    double _getCompReacK(solver::comp_global_id cidx, solver::reac_global_id ridx) const override;
    void _setCompReacK(solver::comp_global_id cidx,
                       solver::reac_global_id ridx,
                       double kf) override;

    bool _getCompReacActive(solver::comp_global_id cidx,
                            solver::reac_global_id ridx) const override;
    void _setCompReacActive(solver::comp_global_id cidx,
                            solver::reac_global_id ridx,
                            bool a) override;

    double _getCompReacH(solver::comp_global_id cidx, solver::reac_global_id ridx) const override;
    double _getCompReacC(solver::comp_global_id cidx, solver::reac_global_id ridx) const override;
    long double _getCompReacA(solver::comp_global_id cidx,
                              solver::reac_global_id ridx) const override;

    unsigned long long _getCompReacExtent(solver::comp_global_id cidx,
                                          solver::reac_global_id ridx) const override;
    void _resetCompReacExtent(solver::comp_global_id cidx, solver::reac_global_id ridx) override;

    double _getCompDiffD(solver::comp_global_id cidx, solver::diff_global_id didx) const override;
    void _setCompDiffD(solver::comp_global_id cidx,
                       solver::diff_global_id didx,
                       double dk) override;

    bool _getCompDiffActive(solver::comp_global_id cidx,
                            solver::diff_global_id didx) const override;
    void _setCompDiffActive(solver::comp_global_id cidx,
                            solver::diff_global_id didx,
                            bool act) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      VESICLE-RELATED
    ////////////////////////////////////////////////////////////////////////

    uint _getCompVesicleCount(solver::comp_global_id cidx,
                              solver::vesicle_global_id vidx) const override;
    void _setCompVesicleCount(solver::comp_global_id cidx,
                              solver::vesicle_global_id vidx,
                              uint n) override;

    solver::vesicle_individual_id _addCompVesicle(solver::comp_global_id cidx,
                                                  solver::vesicle_global_id vidx) override;
    void _deleteSingleVesicle(solver::vesicle_global_id vidx,
                              solver::vesicle_individual_id ves_unique_index) override;

    std::vector<double> _getSingleSpecPosSpherical(
        solver::spec_global_id sidx,
        solver::pointspec_individual_id ps_unique_id) const override;

    uint _getSingleVesicleSurfaceLinkSpecCount(solver::vesicle_global_id vidx,
                                               solver::vesicle_individual_id ves_unique_index,
                                               solver::linkspec_global_id lsidx) const override;

    std::vector<solver::linkspec_individual_id> _getSingleVesicleSurfaceLinkSpecIndices(
        solver::vesicle_global_id vidx,
        solver::vesicle_individual_id ves_unique_index,
        solver::linkspec_global_id lsidx) const override;

    std::vector<solver::vesicle_individual_id> _getCompVesicleIndices(
        solver::comp_global_id cidx,
        solver::vesicle_global_id vidx) const override;

    solver::comp_global_id _getSingleVesicleCompartment(
        solver::vesicle_global_id vidx,
        solver::vesicle_individual_id ves_unique_index) const override;

    std::vector<double> _getSingleVesiclePos(
        solver::vesicle_global_id vidx,
        solver::vesicle_individual_id ves_unique_index) const override;

    void _setCompSingleVesiclePos(solver::comp_global_id cidx,
                                  solver::vesicle_global_id vidx,
                                  solver::vesicle_individual_id ves_unique_index,
                                  const std::vector<double>& pos,
                                  bool force) override;

    uint _getCompVesicleSurfaceSpecCount(solver::comp_global_id cidx,
                                         solver::vesicle_global_id vidx,
                                         solver::spec_global_id sidx) const override;

    uint _getCompVesicleInnerSpecCount(solver::comp_global_id cidx,
                                       solver::vesicle_global_id vidx,
                                       solver::spec_global_id sidx) const override;

    std::map<solver::vesicle_individual_id, uint> _getCompVesicleSurfaceSpecCountMap(
        solver::comp_global_id cidx,
        solver::vesicle_global_id vidx,
        solver::spec_global_id sidx) const override;

    uint _getSingleVesicleSurfaceSpecCount(solver::vesicle_global_id vidx,
                                           solver::vesicle_individual_id ves_unique_index,
                                           solver::spec_global_id sidx) const override;

    uint _getSingleVesicleInnerSpecCount(solver::vesicle_global_id vidx,
                                         solver::vesicle_individual_id ves_unique_index,
                                         solver::spec_global_id sidx) const override;

    void _setSingleVesicleSurfaceSpecCount(solver::vesicle_global_id vidx,
                                           solver::vesicle_individual_id ves_unique_index,
                                           solver::spec_global_id sidx,
                                           uint c) override;

    std::vector<solver::pointspec_individual_id> _getSingleVesicleSurfaceSpecIndices(
        solver::vesicle_global_id vidx,
        solver::vesicle_individual_id ves_unique_index,
        solver::spec_global_id sidx) const override;

    std::vector<std::vector<double>> _getSingleVesicleSurfaceSpecPos(
        solver::vesicle_global_id vidx,
        solver::vesicle_individual_id ves_unique_index,
        solver::spec_global_id sidx) override;

    std::vector<std::vector<double>> _getSingleVesicleSurfaceSpecPosSpherical(
        solver::vesicle_global_id vidx,
        solver::vesicle_individual_id ves_unique_index,
        solver::spec_global_id sidx) const override;

    void _setSingleVesicleSurfaceSpecPosSpherical(
        solver::vesicle_global_id vidx,
        solver::vesicle_individual_id ves_unique_index,
        solver::spec_global_id sidx,
        const std::vector<std::vector<double>>& pos_spherical) override;

    void _setSingleVesicleInnerSpecCount(solver::vesicle_global_id vidx,
                                         solver::vesicle_individual_id ves_unique_index,
                                         solver::spec_global_id sidx,
                                         uint c) override;

    std::map<solver::vesicle_individual_id, uint> _getCompVesicleSurfaceLinkSpecCountMap(
        solver::comp_global_id cidx,
        solver::vesicle_global_id vidx,
        solver::linkspec_global_id lsidx) const override;

    uint _getCompVesicleSurfaceLinkSpecCount(solver::comp_global_id cidx,
                                             solver::vesicle_global_id vidx,
                                             solver::linkspec_global_id lsidx) const override;

    std::vector<std::vector<double>> _getSingleVesicleSurfaceLinkSpecPos(
        solver::vesicle_global_id vidx,
        solver::vesicle_individual_id ves_unique_index,
        solver::linkspec_global_id lsidx) override;

    std::vector<double> _getSingleLinkSpecPos(
        solver::linkspec_individual_id ls_unique_id) const override;

    solver::linkspec_individual_id _getSingleLinkSpecLinkedTo(
        solver::linkspec_individual_id ls_unique_id) const override;

    solver::vesicle_individual_id _getSingleLinkSpecVes(
        solver::linkspec_individual_id ls_unique_id) const override;

    void _setVesicleSurfaceLinkSpecSDiffD(solver::vesicle_global_id vidx,
                                          solver::linkspec_global_id lsidx,
                                          double d) override;

    uint _getSingleVesicleImmobility(solver::vesicle_global_id vidx,
                                     solver::vesicle_individual_id ves_unique_index) const override;

    std::vector<tetrahedron_global_id> _getSingleVesicleOverlapTets(
        solver::vesicle_global_id vidx,
        solver::vesicle_individual_id ves_unique_index) const override;

    void _setTetVesicleDcst(tetrahedron_global_id tidx,
                            solver::vesicle_global_id vidx,
                            double dcst) override;

    void _addVesicleDiffusionGroup(
        solver::vesicle_global_id vidx,
        const std::vector<solver::comp_global_id>& comp_indices) override;


    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      VESICLE-RELATED (NOT COMPARTMENT/PATCH SPECIFIC)
    ////////////////////////////////////////////////////////////////////////

    void _setVesSReacK(solver::vessreac_global_id vsridx, double kf) override;

    uint _getVesSReacExtent(solver::vessreac_global_id vsridx) const override;

    void _setExocytosisK(solver::exocytosis_global_id exoidx, double kf) override;

    uint _getExocytosisExtent(solver::exocytosis_global_id exoidx) const override;

    std::vector<solver::ExocytosisEvent> _getExocytosisEvents(
        solver::exocytosis_global_id exoidx) override;

    uint _getRaftEndocytosisExtent(solver::raftendocytosis_global_id rendoidx) const override;

    std::vector<solver::RaftEndocytosisEvent> _getRaftEndocytosisEvents(
        solver::raftendocytosis_global_id rendoidx) override;

    void _setRaftEndocytosisK(solver::raftendocytosis_global_id rendoidx, double kcst) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS:
    //      VESICLE PATHS
    ////////////////////////////////////////////////////////////////////////

    void _addPathVesicle(std::string const& path_name,
                         solver::vesicle_global_id ves_idx,
                         double speed,
                         const std::map<solver::spec_global_id, uint>& spec_deps,
                         const std::vector<double>& stoch_stepsize) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      PATCH
    ////////////////////////////////////////////////////////////////////////

    void _setPatchRaftCount(solver::patch_global_id pidx,
                            solver::raft_global_id ridx,
                            uint n) override;
    uint _getPatchRaftCount(solver::patch_global_id pidx,
                            solver::raft_global_id ridx) const override;

    uint _getSingleRaftImmobility(solver::raft_global_id ridx,
                                  solver::raft_individual_id raft_unique_index) const override;

    std::vector<solver::raft_individual_id> _getPatchRaftIndices(
        solver::patch_global_id pidx,
        solver::raft_global_id ridx) const override;

    solver::patch_global_id _getSingleRaftPatch(
        solver::raft_global_id ridx,
        solver::raft_individual_id raft_unique_index) const override;

    std::map<solver::raft_individual_id, uint> _getPatchRaftSpecCountMap(
        solver::patch_global_id pidx,
        solver::raft_global_id ridx,
        solver::spec_global_id sidx) const override;

    uint _getPatchRaftSpecCount(solver::patch_global_id pidx,
                                solver::raft_global_id ridx,
                                solver::spec_global_id sidx) const override;

    uint _getSingleRaftSpecCount(solver::raft_global_id ridx,
                                 solver::raft_individual_id raft_unique_index,
                                 solver::spec_global_id sidx) const override;

    void _setSingleRaftSpecCount(solver::raft_global_id ridx,
                                 solver::raft_individual_id raft_unique_index,
                                 solver::spec_global_id sidx,
                                 uint c) override;

    std::vector<double> _getSingleRaftPos(
        solver::raft_global_id ridx,
        solver::raft_individual_id raft_unique_index) const override;

    double _getSingleRaftRaftEndocytosisK(
        solver::raft_global_id ridx,
        solver::raft_individual_id raft_unique_index,
        solver::raftendocytosis_global_id rendoidx) const override;

    void _setSingleRaftRaftEndocytosisK(solver::raft_global_id ridx,
                                        solver::raft_individual_id raft_unique_index,
                                        solver::raftendocytosis_global_id rendoidx,
                                        double k) override;

    void _setSingleRaftSReacActive(solver::raft_global_id ridx,
                                   solver::raft_individual_id raft_unique_index,
                                   solver::raftsreac_global_id rsreac,
                                   bool active) override;

    bool _getSingleRaftSReacActive(solver::raft_global_id ridx,
                                   solver::raft_individual_id raft_unique_index,
                                   solver::raftsreac_global_id rsreac) const override;

    inline void recordRaft(solver::raft_individual_id raft_unique_index, Raft* raft) {
        AssertLog(pRafts.find(raft_unique_index) == pRafts.end());
        pRafts[raft_unique_index] = raft;
    }

    inline void removeRaft(solver::raft_individual_id raft_unique_index, Raft* raft) {
        AssertLog(pRafts[raft_unique_index] == raft);
        pRafts.erase(raft_unique_index);
    }

    double _getPatchArea(solver::patch_global_id pidx) const override;

    double _getPatchSpecCount(solver::patch_global_id pidx,
                              solver::spec_global_id sidx) const override;
    void _setPatchSpecCount(solver::patch_global_id pidx,
                            solver::spec_global_id sidx,
                            double n) override;

    double _getPatchSpecAmount(solver::patch_global_id pidx,
                               solver::spec_global_id sidx) const override;
    void _setPatchSpecAmount(solver::patch_global_id pidx,
                             solver::spec_global_id sidx,
                             double a) override;

    bool _getPatchSpecClamped(solver::patch_global_id pidx,
                              solver::spec_global_id sidx) const override;
    void _setPatchSpecClamped(solver::patch_global_id pidx,
                              solver::spec_global_id sidx,
                              bool buf) override;

    double _getPatchSReacK(solver::patch_global_id pidx,
                           solver::sreac_global_id ridx) const override;
    void _setPatchSReacK(solver::patch_global_id pidx,
                         solver::sreac_global_id ridx,
                         double kf) override;

    bool _getPatchSReacActive(solver::patch_global_id pidx,
                              solver::sreac_global_id ridx) const override;
    void _setPatchSReacActive(solver::patch_global_id pidx,
                              solver::sreac_global_id ridx,
                              bool a) override;

    double _getPatchSReacH(solver::patch_global_id pidx,
                           solver::sreac_global_id ridx) const override;
    double _getPatchSReacC(solver::patch_global_id pidx,
                           solver::sreac_global_id ridx) const override;
    double _getPatchSReacA(solver::patch_global_id pidx,
                           solver::sreac_global_id ridx) const override;

    unsigned long long _getPatchSReacExtent(solver::patch_global_id pidx,
                                            solver::sreac_global_id ridx) const override;
    void _resetPatchSReacExtent(solver::patch_global_id pidx,
                                solver::sreac_global_id ridx) override;

    bool _getPatchVDepSReacActive(solver::patch_global_id pidx,
                                  solver::vdepsreac_global_id) const override;
    void _setPatchVDepSReacActive(solver::patch_global_id pidx,
                                  solver::vdepsreac_global_id,
                                  bool a) override;

    void _setPatchEndocyticZoneEndocytosisActive(solver::patch_global_id pidx,
                                                 std::string const& zone,
                                                 solver::endocytosis_global_id endogidx,
                                                 bool active) override;

    void _setPatchEndocyticZoneEndocytosisK(solver::patch_global_id pidx,
                                            std::string const& zone,
                                            solver::endocytosis_global_id endogidx,
                                            double k) override;

    uint _getPatchEndocyticZoneEndocytosisExtent(
        solver::patch_global_id pidx,
        std::string const& zone,
        solver::endocytosis_global_id endogidx) const override;

    std::vector<solver::EndocytosisEvent> _getPatchEndocyticZoneEndocytosisEvents(
        solver::patch_global_id pidx,
        std::string const& zone,
        solver::endocytosis_global_id endogidx) const override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      DIFFUSION BOUNDARIES
    ////////////////////////////////////////////////////////////////////////

    void _setDiffBoundarySpecDiffusionActive(solver::diffboundary_global_id dbidx,
                                             solver::spec_global_id sidx,
                                             bool act) override;
    bool _getDiffBoundarySpecDiffusionActive(solver::diffboundary_global_id dbidx,
                                             solver::spec_global_id sidx) const override;
    void _setDiffBoundarySpecDcst(solver::diffboundary_global_id dbidx,
                                  solver::spec_global_id sidx,
                                  double dcst,
                                  solver::comp_global_id direction_comp = {}) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      SURFACE DIFFUSION BOUNDARIES
    ////////////////////////////////////////////////////////////////////////

    void _setSDiffBoundarySpecDiffusionActive(solver::sdiffboundary_global_id sdbidx,
                                              solver::spec_global_id sidx,
                                              bool act) override;
    bool _getSDiffBoundarySpecDiffusionActive(solver::sdiffboundary_global_id sdbidx,
                                              solver::spec_global_id sidx) const override;
    void _setSDiffBoundarySpecDcst(solver::sdiffboundary_global_id sdbidx,
                                   solver::spec_global_id sidx,
                                   double dcst,
                                   solver::patch_global_id direction_patch = {}) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TETRAHEDRAL VOLUME ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    double _getTetVol(tetrahedron_global_id tidx) const override;
    void _setTetVol(tetrahedron_global_id tidx, double vol) override;

    bool _getTetSpecDefined(tetrahedron_global_id tidx, solver::spec_global_id sidx) const override;

    double _getTetSpecCount(tetrahedron_global_id tidx, solver::spec_global_id sidx) const override;
    void _setTetSpecCount(tetrahedron_global_id tidx,
                          solver::spec_global_id sidx,
                          double n) override;

    double _getTetSpecAmount(tetrahedron_global_id tidx,
                             solver::spec_global_id sidx) const override;
    void _setTetSpecAmount(tetrahedron_global_id tidx,
                           solver::spec_global_id sidx,
                           double m) override;

    double _getTetSpecConc(tetrahedron_global_id tidx, solver::spec_global_id sidx) const override;
    void _setTetSpecConc(tetrahedron_global_id tidx,
                         solver::spec_global_id sidx,
                         double c) override;

    bool _getTetSpecClamped(tetrahedron_global_id tidx, solver::spec_global_id sidx) const override;
    void _setTetSpecClamped(tetrahedron_global_id tidx,
                            solver::spec_global_id sidx,
                            bool buf) override;

    double _getTetReacK(tetrahedron_global_id tidx, solver::reac_global_id ridx) const override;
    void _setTetReacK(tetrahedron_global_id tidx, solver::reac_global_id ridx, double kf) override;

    bool _getTetReacActive(tetrahedron_global_id tidx, solver::reac_global_id ridx) const override;
    void _setTetReacActive(tetrahedron_global_id tidx,
                           solver::reac_global_id ridx,
                           bool act) override;

    double _getTetDiffD(tetrahedron_global_id tidx,
                        solver::diff_global_id didx,
                        tetrahedron_global_id direction_tet = {}) const override;
    void _setTetDiffD(tetrahedron_global_id tidx,
                      solver::diff_global_id didx,
                      double dk,
                      tetrahedron_global_id direction_tet = {}) override;

    bool _getTetDiffActive(tetrahedron_global_id tidx, solver::diff_global_id didx) const override;
    void _setTetDiffActive(tetrahedron_global_id tidx,
                           solver::diff_global_id didx,
                           bool act) override;

    ////////////////////////////////////////////////////////////////////////

    double _getTetReacH(tetrahedron_global_id tidx, solver::reac_global_id ridx) const override;
    double _getTetReacC(tetrahedron_global_id tidx, solver::reac_global_id ridx) const override;
    double _getTetReacA(tetrahedron_global_id tidx, solver::reac_global_id ridx) const override;

    double _getTetDiffA(tetrahedron_global_id tidx, solver::diff_global_id didx) const override;

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    double _getTetV(tetrahedron_global_id tidx) const override;
    void _setTetV(tetrahedron_global_id tidx, double v) override;
    bool _getTetVClamped(tetrahedron_global_id tidx) const override;
    void _setTetVClamped(tetrahedron_global_id tidx, bool cl) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TRIANGULAR SURFACE ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    void _setTriRaftCount(triangle_global_id tidx, solver::raft_global_id ridx, uint n) override;
    uint _getTriRaftCount(triangle_global_id tidx, solver::raft_global_id ridx) const override;

    solver::raft_individual_id _addTriRaft(triangle_global_id tidx,
                                           solver::raft_global_id ridx) override;

    double _getTriArea(triangle_global_id tidx) const override;
    void _setTriArea(triangle_global_id tidx, double area) override;

    bool _getTriSpecDefined(triangle_global_id tidx, solver::spec_global_id sidx) const override;

    double _getTriSpecCount(triangle_global_id tidx, solver::spec_global_id sidx) const override;
    void _setTriSpecCount(triangle_global_id tidx, solver::spec_global_id sidx, double n) override;

    double _getTriSpecAmount(triangle_global_id tidx, solver::spec_global_id sidx) const override;
    void _setTriSpecAmount(triangle_global_id tidx, solver::spec_global_id sidx, double m) override;

    bool _getTriSpecClamped(triangle_global_id tidx, solver::spec_global_id sidx) const override;
    void _setTriSpecClamped(triangle_global_id tidx,
                            solver::spec_global_id sidx,
                            bool buf) override;

    double _getTriSReacK(triangle_global_id tidx, solver::sreac_global_id ridx) const override;
    void _setTriSReacK(triangle_global_id tidx, solver::sreac_global_id ridx, double kf) override;

    bool _getTriSReacActive(triangle_global_id tidx, solver::sreac_global_id ridx) const override;
    void _setTriSReacActive(triangle_global_id tidx,
                            solver::sreac_global_id ridx,
                            bool act) override;

    double _getTriSDiffD(triangle_global_id tidx,
                         solver::surfdiff_global_id didx,
                         triangle_global_id direction_tri = {}) const override;
    void _setTriSDiffD(triangle_global_id tidx,
                       solver::surfdiff_global_id didx,
                       double dk,
                       triangle_global_id direction_tri = {}) override;

    // bool _getTriExocytosisActive(tetrahedron_global_id tidx, uint eidx)
    // const; void _setTriExocytosisActive(tetrahedron_global_id tidx, uint
    // eidx, bool act);

    ////////////////////////////////////////////////////////////////////////

    double _getTriSReacH(triangle_global_id tidx, solver::sreac_global_id ridx) const override;
    double _getTriSReacC(triangle_global_id tidx, solver::sreac_global_id ridx) const override;
    double _getTriSReacA(triangle_global_id tidx, solver::sreac_global_id ridx) const override;

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    double _getTriV(triangle_global_id tidx) const override;
    void _setTriV(triangle_global_id tidx, double v) override;
    bool _getTriVClamped(triangle_global_id tidx) const override;
    void _setTriVClamped(triangle_global_id tidx, bool cl) override;

    double _getTriOhmicErev(triangle_global_id tidx,
                            solver::ohmiccurr_global_id ocgidx) const override;
    void _setTriOhmicErev(triangle_global_id tidx,
                          solver::ohmiccurr_global_id ocgidx,
                          double erev) override;

    double _getTriOhmicI(triangle_global_id tidx) const override;
    double _getTriOhmicI(triangle_global_id tidx, solver::ohmiccurr_global_id ocidx) const override;

    double _getTriGHKI(triangle_global_id tidx) const override;
    double _getTriGHKI(triangle_global_id tidx, solver::ghkcurr_global_id ghkidx) const override;

    double _getTriI(triangle_global_id tidx) const override;

    double _getTriIClamp(triangle_global_id tidx) const override;
    void _setTriIClamp(triangle_global_id tidx, double cur) override;

    bool _getTriVDepSReacActive(triangle_global_id tidx,
                                solver::vdepsreac_global_id) const override;
    void _setTriVDepSReacActive(triangle_global_id tidx,
                                solver::vdepsreac_global_id,
                                bool act) override;

    void _setTriCapac(triangle_global_id tidx, double cap) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      VERTICES ELEMENTS
    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    double _getVertV(vertex_id_t vidx) const override;
    void _setVertV(vertex_id_t vidx, double v) override;

    bool _getVertVClamped(vertex_id_t vidx) const override;
    void _setVertVClamped(vertex_id_t vidx, bool cl) override;

    double _getVertIClamp(vertex_id_t tidx) const override;
    void _setVertIClamp(vertex_id_t tidx, double cur) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      MEMBRANE AND VOLUME CONDUCTOR
    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    void _setMembPotential(solver::membrane_global_id midx, double v) override;
    void _setMembCapac(solver::membrane_global_id midx, double cm) override;
    void _setMembVolRes(solver::membrane_global_id midx, double ro) override;
    void _setMembRes(solver::membrane_global_id midx, double ro, double vrev) override;

    ////////////////////////////////////////////////////////////////////////

    // Checked global to local index translations

    solver::spec_local_id _specG2L_or_throw(CompVesRaft* comp, solver::spec_global_id gidx) const;
    solver::spec_local_id _specG2L_or_throw(PatchVesRaft* patch, solver::spec_global_id gidx) const;
#if 0
    solver::spec_local_id _specG2L_or_throw(TetVesRaft *tet, solver::spec_global_id  gidx) const;
    solver::spec_local_id _specG2L_or_throw(TetVesRaft *tri, solver::spec_global_id  gidx) const;
#endif
    solver::reac_local_id _reacG2L_or_throw(CompVesRaft* comp, solver::reac_global_id gidx) const;
    solver::sreac_local_id _sreacG2L_or_throw(PatchVesRaft* patch,
                                              solver::sreac_global_id gidx) const;
    solver::diff_local_id _diffG2L_or_throw(CompVesRaft* comp, solver::diff_global_id gidx) const;
    // solver::surfdiff_local_id _sdiffG2L_or_throw(PatchVesRaft* patch, solver::surfdiff_global_id
    // gidx) const; solver::vdepsreac_local_id _vdepsreacG2L_or_throw(PatchVesRaft* patch,
    // solver::vdepsreac_global_id gidx) const;
    solver::endocytosis_local_id _endoG2L_or_throw(PatchVesRaft* patch,
                                                   solver::endocytosis_global_id gidx) const;

    ////////////////////////////////////////////////////////////////////////
    // HYBRID SOLVER METHODS
    ////////////////////////////////////////////////////////////////////////

    uint _addComp(solver::Compdef* cdef, tetmesh::Tetmesh* mesh);

    uint _addPatch(solver::Patchdef* pdef, tetmesh::Tetmesh* mesh);

    // uint _addDiffBoundary(solver::DiffBoundarydef *dbdef);

    // uint _addSDiffBoundary(solver::SDiffBoundarydef *sdbdef);

    void _addTet(tetrahedron_global_id tetidx,
                 CompVesRaft* comp,
                 double vol,
                 double a1,
                 double a2,
                 double a3,
                 double a4,
                 double d1,
                 double d2,
                 double d3,
                 double d4,
                 tetrahedron_global_id tet0,
                 tetrahedron_global_id tet1,
                 tetrahedron_global_id tet2,
                 tetrahedron_global_id tet3,
                 math::point3d baryc);

    // void _addTet(solver::comp_global_id cidx, Comp * comp,
    // double vol);

    void _addTri(triangle_global_id triidx,
                 PatchVesRaft* patch,
                 double area,
                 double l0,
                 double l1,
                 double l2,
                 double d0,
                 double d1,
                 double d2,
                 tetrahedron_global_id tinner,
                 tetrahedron_global_id touter,
                 triangle_global_id tri0,
                 triangle_global_id tri1,
                 triangle_global_id tri2,
                 math::point3d baryc,
                 math::point3d trinorm);

    // initialization internal method
    // partition the solver and the mesh
    void _partition();

    // setup stuffs for both master and client first
    // then call _setupMaster()/_setupClient()
    // based on the rank
    void _setup();

    ////////////////////////////////////////////////////////////////////////

    tetmesh::Tetmesh* pMesh;

    ////////////////////////////////////////////////////////////////////////
    // LIST OF HYBRID SOLVER OBJECTS
    ////////////////////////////////////////////////////////////////////////

    util::strongid_vector<solver::comp_global_id, CompVesRaft*> pComps;
    util::strongid_vector<solver::patch_global_id, PatchVesRaft*> pPatches;
    util::strongid_vector<triangle_global_id, TriVesRaft*> pTris;
    util::strongid_vector<tetrahedron_global_id, TetVesRaft*> pTets;
    TetVesRaft* _getTet(tetrahedron_global_id tgidx) const;
    TriVesRaft* _getTri(triangle_global_id tgidx) const;
    // For overlap
    util::strongid_vector<tetrahedron_global_id, Tetrahedron*> pTet_ext;

    double pTemp{0.0};
    double pEFDT{1.0e-5};
    ////////////////////////// ADDED FOR VESICLES ////////////////////////////

    double pVesicledt{1.0e-3};
    double pVesicleDefaultdt{1.0e-3};

    // To be used whenever anything changes in vesicles that RDEF needs to know
    // about
    bool pRequireVesicleCommunication;

    // Q values are now held at solver level and are consistent across
    // compartments

    // Table of Q values- will be set up once for every species
    // that ends up residing on the surface of these vesicles
    // dt will implicitly be constant throughout the simulation,
    // whatever is given to the vesicle at the time of construction
    // Start with a null pointer for species that are not setup,
    // and check this carefully of course
    std::map<solver::vesicle_global_id, util::strongid_vector<solver::spec_global_id, Qtable*>>
        pQtables_spec;
    uint pQtablesize_spec;

    std::map<solver::vesicle_global_id, util::strongid_vector<solver::linkspec_global_id, Qtable*>>
        pQtables_linkspec;
    uint pQtablesize_linkspec;

    // The diffusion rates for species on vesicle surfaces.
    // Start as a vector of 0s for all species, only changing those that
    // are changed by solver method setCompVesicleSpecDiffD()
    std::map<solver::vesicle_global_id, util::strongid_vector<solver::spec_global_id, double>>
        pVesSpecD;
    std::map<solver::vesicle_global_id, util::strongid_vector<solver::linkspec_global_id, double>>
        pVesLinkSpecD;

    // Also store the total number of vesicles of any type. Each vesicle now gets a unique index
    // regardless of type
    uint pVesicles_count{0};

    // Also store the total number of rafts of any type.
    uint pRafts_count{0};

    // Also give PointSpecs a unique index.
    solver::pointspec_individual_id pNextPointSpecUniqueID{0};
    solver::pointspec_individual_id pRDEFminPointSpecUniqueID;

    std::map<std::string, Path*> pPaths;

    std::map<solver::vesicle_individual_id, Vesicle*> pVesicles;
    std::map<solver::raft_individual_id, Raft*> pRafts;

    ////////////////////////////////////////////////////////////////////////

    // Vesicle diffusion algorithm parameters

    // Squared factor for maximum walking distance
    double maxWalkDistSqFact{9};
    // Minimum number of visited tetrahedrons before exiting walk
    uint minNbTetVisited{10};

    ////////////////////////////////////////////////////////////////////////

    std::set<const LinkSpecPair*, util::DerefPtrLess<LinkSpecPair>> pLinkSpecPairs;

    // Easy access to LinkSpecs from their ids, Linkspecs are still owned by vesicles
    std::map<solver::linkspec_individual_id, LinkSpec*> pLinkSpecs;

    // Quick access to vesicle unbinds that apply to a linkspec pair on specific vesicles
    // It maps a link spec pair and vesicle pair to all corresponding vesunbind reactions
    // and the sum of their rates:
    // (ves1, ves2, ls1, ls2) -> ([vub1, vub2, ...], total_rate)
    // If the rate of vesicleunbind reactions are changed, this needs to be updated
    std::map<std::tuple<solver::vesicle_global_id,
                        solver::vesicle_global_id,
                        solver::linkspec_global_id,
                        solver::linkspec_global_id>,
             std::pair<std::vector<solver::VesUnbinddef*>, double>>
        pLinkSpecPair2VesUnbinds;

    // Final assigment of unique IDs comes from this solver, since they are only unique per RDEF
    // core
    solver::linkspec_individual_id pNextLinkSpecUniqueID{uint(0)};

    ////////////////////////// ADDED FOR MPI ////////////////////////////

    std::map<tetrahedron_global_id, int> tetHosts;
    std::map<triangle_global_id, int> triHosts;

    int myRank_World{};
    int nHosts_World{};
    int vesraftRank_World{};
    int RDEFmasterRank_World{};
    bool syncOutput{false};
    int outputRank{0};

    int _getTetHost(tetrahedron_global_id tgidx) const;
    int _getTriHost(triangle_global_id tgidx) const;

    // VesRdef <-> RDEF Sync
    MPIDataTypeUtil dataTypeUtil{};

    // sync all registered counts in tetPoolCountSyncs and triPoolCountSyncs
    void _syncPools(SyncDirection direction);

    // data for sync
    std::vector<PoolCountSync> tetPoolCountSyncs_Vec;
    std::vector<PoolCountSync> triPoolCountSyncs_Vec;

    std::vector<TetV2R> tetV2R_Vec;

    // Vesicles and associated surface species and link species
    std::vector<VesProxyV2R> vesProxyV2R_Vec;
    std::vector<VesSurfSpecV2R> vesSurfSpecV2R_Vec;
    std::vector<VesLinkSpecV2R> vesLinkSpecV2R_Vec;
    void _constructVesV2R();

    std::vector<VesProxyR2V> vesProxyR2V_Vec;
    std::vector<VesSurfSpecR2V> vesSurfSpecR2V_Vec;
    std::vector<VesInnerSpecR2V> vesInnerSpecR2V_Vec;
    std::vector<VesLinkSpecPairR2V> vesLinkSpecPairR2V_Vec;
    std::vector<VesLinkSpecR2V> vesLinkSpecR2V_Vec;
    void _useVesR2V();

    // Rafts and associated species

    std::vector<RaftProxyV2R> raftProxyV2R_Vec;
    std::vector<RaftSpecV2R> raftSpecV2R_Vec;
    std::vector<RaftSReacInactiveV2R> raftSReacInactiveV2R_Vec;
    void _constructRaftV2R();

    std::vector<RaftProxyR2V> raftProxyR2V_Vec;
    std::vector<RaftSpecR2V> raftSpecR2V_Vec;
    void _useRaftR2V();

    std::vector<RaftGenCountR2V> triRaftGenR2V_Vec;
};

}  // namespace steps::mpi::tetvesicle
