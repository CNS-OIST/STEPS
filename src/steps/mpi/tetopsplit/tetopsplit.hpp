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

#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>


#include "comp.hpp"
#include "crstruct.hpp"
#include "diffboundary.hpp"
#include "kproc.hpp"
#include "patch.hpp"
#include "sdiffboundary.hpp"
#include "tet.hpp"
#include "tri.hpp"
#include "wmvol.hpp"

#include "geom/tetmesh.hpp"
#include "solver/api.hpp"
#include "solver/efield/efield.hpp"
#include "solver/statedef.hpp"

namespace steps::mpi::tetopsplit {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.

// Auxiliary declarations.
typedef solver::kproc_global_id SchedIDX;
typedef std::set<SchedIDX> SchedIDXSet;
typedef SchedIDXSet::iterator SchedIDXSetI;
typedef SchedIDXSet::const_iterator SchedIDXSetCI;
typedef std::vector<SchedIDX> SchedIDXVec;
typedef SchedIDXVec::iterator SchedIDXVecI;
typedef SchedIDXVec::const_iterator SchedIDXVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Copies the contents of a set of SchedIDX entries into a vector.
/// The contents of the vector are completely overridden.
///
extern void schedIDXSet_To_Vec(SchedIDXSet const& s, SchedIDXVec& v);

////////////////////////////////////////////////////////////////////////////////

enum SubVolType { SUB_WM, SUB_TET, SUB_TRI };

////////////////////////////////////////////////////////////////////////////////

class TetOpSplitP: public solver::API {
  public:
    TetOpSplitP(model::Model* m,
                wm::Geom* g,
                const rng::RNGptr& r,
                int calcMembPot = EF_NONE,
                std::vector<int> const& tet_hosts = {},
                const std::map<triangle_global_id, int>& tri_hosts = {},
                std::vector<int> const& wm_hosts = {});

    ~TetOpSplitP() override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER INFORMATION
    ////////////////////////////////////////////////////////////////////////

    std::string getSolverName() const override;
    std::string getSolverDesc() const override;
    std::string getSolverAuthors() const override;
    std::string getSolverEmail() const override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS
    ////////////////////////////////////////////////////////////////////////

    void reset() override;
    void run(double endtime) override;
    void advance(double adv) override;
    void step() override;

    void checkpoint(std::string const& file_name) override;
    void restore(std::string const& file_name) override;
    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    void setEfieldDT(double efdt) override;

    inline double efdt() const noexcept {
        return pEFDT;
    }

    inline double getEfieldDT() const noexcept override {
        return pEFDT;
    }

    void setTemp(double t) override;

    inline double getTemp() const noexcept override {
        return pTemp;
    }

    // save the optimal vertex indexing
    void saveMembOpt(std::string const& opt_file_name);

    ////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      GENERAL
    ////////////////////////////////////////////////////////////////////////

    double getTime() const override;

    inline double getA0() const noexcept override {
        return pA0;
    }

    uint getNSteps() const override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      ADVANCE
    //      Developer only
    ////////////////////////////////////////////////////////////////////////

    void setTime(double time) override;
    void setNSteps(uint nsteps) override;

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
    // SOLVER STATE ACCESS:
    //      PATCH
    ////////////////////////////////////////////////////////////////////////

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
                           solver::sreac_global_id sridx) const override;
    void _setPatchSReacK(solver::patch_global_id pidx,
                         solver::sreac_global_id sridx,
                         double kf) override;

    bool _getPatchSReacActive(solver::patch_global_id pidx,
                              solver::sreac_global_id sridx) const override;
    void _setPatchSReacActive(solver::patch_global_id pidx,
                              solver::sreac_global_id sridx,
                              bool a) override;

    double _getPatchSReacH(solver::patch_global_id pidx,
                           solver::sreac_global_id sridx) const override;
    double _getPatchSReacC(solver::patch_global_id pidx,
                           solver::sreac_global_id sridx) const override;
    double _getPatchSReacA(solver::patch_global_id pidx,
                           solver::sreac_global_id sridx) const override;

    unsigned long long _getPatchSReacExtent(solver::patch_global_id pidx,
                                            solver::sreac_global_id sridx) const override;
    void _resetPatchSReacExtent(solver::patch_global_id pidx,
                                solver::sreac_global_id sridx) override;

    bool _getPatchVDepSReacActive(solver::patch_global_id pidx,
                                  solver::vdepsreac_global_id vsridx) const override;
    void _setPatchVDepSReacActive(solver::patch_global_id pidx,
                                  solver::vdepsreac_global_id vsridx,
                                  bool a) override;

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

    const Tet& _getTet(tetrahedron_global_id tidx) const;
    Tet& _getTet(tetrahedron_global_id tidx);

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

    const Tri& _getTri(triangle_global_id tidx) const;
    Tri& _getTri(triangle_global_id tidx);

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

    double _getTriSReacK(triangle_global_id tidx, solver::sreac_global_id sridx) const override;
    void _setTriSReacK(triangle_global_id tidx, solver::sreac_global_id sridx, double kf) override;

    bool _getTriSReacActive(triangle_global_id tidx, solver::sreac_global_id sridx) const override;
    void _setTriSReacActive(triangle_global_id tidx,
                            solver::sreac_global_id sridx,
                            bool act) override;

    double _getTriSDiffD(triangle_global_id tidx,
                         solver::surfdiff_global_id didx,
                         triangle_global_id direction_tri) const override;

    void _setTriSDiffD(triangle_global_id tidx,
                       solver::surfdiff_global_id didx,
                       double dk,
                       triangle_global_id direction_tri) override;

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

    double _getTriVDepSReacK(triangle_global_id tidx,
                             solver::vdepsreac_global_id vsridx) const override;

    bool _getTriVDepSReacActive(triangle_global_id tidx,
                                solver::vdepsreac_global_id vsridx) const override;
    void _setTriVDepSReacActive(triangle_global_id tidx,
                                solver::vdepsreac_global_id vsridx,
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
    std::pair<double, double> _getMembRes(solver::membrane_global_id midx) const override;

    ////////////////////////////////////////////////////////////////////////

    solver::kproc_global_id addKProc(KProc* kp, bool Vdep = false);

    void addDiff(Diff* diff);
    void addSDiff(SDiff* sdiff);
    inline uint countKProcs() const noexcept {
        return pKProcs.size();
    }

    ////////////////////////////////////////////////////////////////////////

    inline tetmesh::Tetmesh* mesh() const noexcept {
        return pMesh;
    }

    inline Comp* _comp(solver::comp_global_id cidx) const {
        AssertLog(cidx < statedef().countComps());
        AssertLog(statedef().countComps() == pComps.size());
        return pComps[cidx];
    }

    inline const auto& patches() const noexcept {
        return pPatches;
    }

    inline Patch* _patch(solver::patch_global_id pidx) const {
        AssertLog(pidx < statedef().countPatches());
        AssertLog(statedef().countPatches() == pPatches.size());
        return pPatches[pidx];
    }

    inline DiffBoundary* _diffboundary(solver::diffboundary_global_id dbidx) const {
        AssertLog(dbidx < statedef().countDiffBoundaries());
        return pDiffBoundaries[dbidx.get()];
    }

    inline SDiffBoundary* _sdiffboundary(solver::sdiffboundary_global_id sdbidx) const {
        AssertLog(sdbidx < statedef().countSDiffBoundaries());
        return pSDiffBoundaries[sdbidx.get()];
    }

    inline Tet* _tet(tetrahedron_global_id tidx) const noexcept {
        return pTets[tidx];
    }

    inline Tri* _tri(triangle_global_id tidx) const noexcept {
        return pTris[tidx];
    }

    inline double a0() const noexcept {
        return pA0;
    }

    // inline bool built()
    //{ return pBuilt; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS
    ////////////////////////////////////////////////////////////////////////

    uint _addComp(solver::Compdef* cdef);

    uint _addPatch(solver::Patchdef* pdef);

    uint _addDiffBoundary(solver::DiffBoundarydef* dbdef);

    uint _addSDiffBoundary(solver::SDiffBoundarydef* sdbdef);

    void _addTet(tetrahedron_global_id tetidx,
                 Comp* comp,
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
                 tetrahedron_global_id tet3);

    void _addWmVol(solver::comp_global_id cidx, Comp* comp, double vol);

    void _addTri(triangle_global_id triidx,
                 Patch* patch,
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
                 triangle_global_id tri2);

    // called when local tet, tri, reac, sreac objects have been created
    // by constructor
    void _setup();

    void _runWithoutEField(double endtime);
    void _runWithEField(double endtime);
    // void _build();
    void _refreshEFTrisV();

    double _getRate(uint i) const {
        return pKProcs[i]->rate();
    }

    KProc* _getNext() const;

    // void _reset();

    void _executeStep(KProc* kp, double dt, double period = 0.0);
    void _updateSpec(WmVol& tet, solver::spec_global_id spec_gidx);

    /// Update the kproc's of a triangle, after a species has been changed.
    /// This does not need to update the kproc's of any neighbouring
    /// tetrahedrons.
    ///
    /// Rational: no kproc in tetrahedron depends on species changes on triangle
    ///
    /// Use depSpecTri to check if the kproc depends on the spec, therefore use
    /// spec_gidx
    void _updateSpec(Tri& tri, solver::spec_global_id spec_gidx);

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    /// Check the EField flag
    inline bool efflag() const noexcept {
        return pEFoption != EF_NONE;
    }

    void _setupEField();

    inline uint neftets() const noexcept {
        return pEFNTets;
    }

    inline uint neftris() const noexcept {
        return pEFNTris;
    }

    inline uint nefverts() const noexcept {
        return pEFNVerts;
    }

    ////////////////////////////////////////////////////////////////////////
    // Batch Data Access
    ////////////////////////////////////////////////////////////////////////

    std::vector<double> getBatchTetSpecCounts(const std::vector<index_t>& tets,
                                              std::string const& s) const override;

    std::vector<double> getBatchTriSpecCounts(const std::vector<index_t>& tris,
                                              std::string const& s) const override;

    void setBatchTetSpecConcs(const std::vector<index_t>& tets,
                              std::string const& s,
                              const std::vector<double>& concs) override;

    std::vector<double> getBatchTetSpecConcs(const std::vector<index_t>& tets,
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

    void setBatchTetSpecConcsNP(const index_t* indices,
                                size_t ntets,
                                std::string const& s,
                                const double* concs,
                                size_t output_size) override;

    void getBatchTetSpecConcsNP(const index_t* indices,
                                size_t ntets,
                                std::string const& s,
                                double* concs,
                                size_t output_size) const override;

    double sumBatchTetCountsNP(const index_t* indices, size_t input_size, std::string const& s);

    double sumBatchTriCountsNP(const index_t* indices, size_t input_size, std::string const& s);

    double sumBatchTriGHKIsNP(const index_t* indices, size_t input_size, std::string const& ghk);

    double sumBatchTriOhmicIsNP(const index_t* indices, size_t input_size, std::string const& oc);

    void getBatchTriOhmicIsNP(const index_t* indices,
                              size_t input_size,
                              std::string const& oc,
                              double* counts,
                              size_t output_size) const;

    void getBatchTriGHKIsNP(const index_t* indices,
                            size_t input_size,
                            std::string const& ghk,
                            double* counts,
                            size_t output_size) const;

    void getBatchTriVsNP(const index_t* indices,
                         size_t input_size,
                         double* counts,
                         size_t output_size) const;

    void getBatchTetVsNP(const index_t* indices,
                         size_t input_size,
                         double* counts,
                         size_t output_size) const;

    void getBatchTriBatchOhmicIsNP(const index_t* indices,
                                   size_t input_size,
                                   std::vector<std::string>& ocs,
                                   double* counts,
                                   size_t output_size) const;

    void getBatchTriBatchGHKIsNP(const index_t* indices,
                                 size_t input_size,
                                 std::vector<std::string>& ghks,
                                 double* counts,
                                 size_t output_size) const;

    ////////////////////////////////////////////////////////////////////////
    // ROI Data Access
    ////////////////////////////////////////////////////////////////////////

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
    void setROISpecAmount(const std::string& ROI_id, std::string const& s, double amount) override;

    double getROISpecConc(const std::string& ROI_id, std::string const& s) const override;
    void setROISpecConc(const std::string& ROI_id, std::string const& s, double conc) override;

    void setROISpecClamped(const std::string& ROI_id, std::string const& s, bool b) override;

    void setROIReacK(const std::string& ROI_id, std::string const& r, double kf) override;
    void setROISReacK(const std::string& ROI_id, std::string const& sr, double kf) override;
    void setROIDiffD(const std::string& ROI_id, std::string const& d, double dk) override;

    void setROIReacActive(const std::string& ROI_id, std::string const& r, bool a) override;
    void setROISReacActive(const std::string& ROI_id, std::string const& sr, bool a) override;
    void setROIDiffActive(const std::string& ROI_id, std::string const& d, bool a) override;
    void setROIVDepSReacActive(const std::string& ROI_id, std::string const& vsr, bool a) override;

    unsigned long long getROIReacExtent(const std::string& ROI_id,
                                        std::string const& r) const override;
    void resetROIReacExtent(const std::string& ROI_id, std::string const& r) override;

    unsigned long long getROISReacExtent(const std::string& ROI_id,
                                         std::string const& sr) const override;
    void resetROISReacExtent(const std::string& ROI_id, std::string const& sr) override;

    unsigned long long getROIDiffExtent(const std::string& ROI_id,
                                        std::string const& d) const override;
    void resetROIDiffExtent(const std::string& ROI_id, std::string const& d) override;

    ////////////////////////////////////////////////////////////////////////

    //////////////////////////// MPI STUFFS ////////////////////////////

    // Exposed to Python:

    void setDiffApplyThreshold(int threshold);

    unsigned long long getReacExtent(bool local = false);
    unsigned long long getDiffExtent(bool local = false);
    double getNIteration() const;

    double getUpdPeriod() {
        return updPeriod;
    }

    void repartitionAndReset(std::vector<int> const& tet_hosts,
                             std::map<uint, int> const& tri_hosts = {},
                             std::vector<int> const& wm_hosts = {});

    double getCompTime() const;
    double getSyncTime() const;
    double getIdleTime() const;
    double getEFieldTime() const;
    double getRDTime() const;
    double getDataExchangeTime() const;

    // Not currently exposed to Python:
    int getTetHostRank(tetrahedron_global_id tidx) const;
    int getTriHostRank(triangle_global_id tidx) const;
    int getWMVolHostRank(solver::comp_global_id idx) const noexcept {
        return wmHosts[idx];
    }
    // void registerSyncWmVol(WmVol * wmvol);
    // void registerSyncTet(Tet * tet);
    // void registerSyncTri(Tri * tri);
    void addNeighHost(int host);
    void registerBoundaryTet(Tet* tet);
    void registerBoundaryTri(Tri* tri);
    uint registerRemoteMoleculeChange(int svol_host,
                                      uint loc,
                                      SubVolType svol_type,
                                      unsigned long idx,
                                      solver::spec_local_id slidx,
                                      uint change);

  private:
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

    ////////////////////////////////////////////////////////////////////////

    tetmesh::Tetmesh* pMesh{nullptr};

    ////////////////////////////////////////////////////////////////////////
    // LIST OF TETEXACT SOLVER OBJECTS
    ////////////////////////////////////////////////////////////////////////

    util::strongid_vector<solver::comp_global_id, Comp*> pComps;

    util::strongid_vector<solver::patch_global_id, Patch*> pPatches;

    std::vector<DiffBoundary*> pDiffBoundaries;

    std::vector<SDiffBoundary*> pSDiffBoundaries;

    // These objects are used to describe a mesh compartment that is
    // being treated as a well-mixed volume.
    util::strongid_vector<solver::comp_global_id, WmVol*> pWmVols;

    util::strongid_vector<triangle_global_id, Tri*> pTris;

    // Now stored as base pointer
    util::strongid_vector<tetrahedron_global_id, Tet*> pTets;

    ////////////////////////////////////////////////////////////////////////
    // Diffusion Data and Methods
    ////////////////////////////////////////////////////////////////////////
    std::vector<Diff*> pDiffs;

    void _updateDiff(Diff* diff);

    // separator for non-zero and zero propensity diffusions
    uint diffSep{0};

    ////////////////////////////////////////////////////////////////////////
    // Surface Diffusion Data and Methods
    ////////////////////////////////////////////////////////////////////////
    std::vector<SDiff*> pSDiffs;

    void _updateSDiff(SDiff* sdiff);

    // separator for non-zero and zero propensity diffusions
    uint sdiffSep{0};

    ////////////////////////////////////////////////////////////////////////
    // CR SSA Kernel Data and Methods
    ////////////////////////////////////////////////////////////////////////
    uint nEntries{};
    double pSum{};
    double nSum{};
    double pA0{0.0};

    std::vector<KProc*> pKProcs;
    std::vector<KProc*> pVdepKProcs;
    std::vector<CRGroup*> nGroups;
    std::vector<CRGroup*> pGroups;

    ////////////////////////////////////////////////////////////////////////////////

    void _computeUpdPeriod();
    void _updateLocal(KProcPSet const& upd_entries);
    void _updateLocal(std::vector<KProc*> const& upd_entries);
    void _updateLocal(std::vector<uint> const& upd_entries);
    void _updateLocal(uint* upd_entries, uint buffer_size);
    void _updateLocal();
    CRGroup* _getGroup(int pow);
    void _extendPGroups(uint new_size);
    void _extendNGroups(uint new_size);
    void _extendGroup(CRGroup* group, uint size = 1024);
    void _updateSum();
    void _updateElement(KProc* kp);

    ////////////////////////////////////////////////////////////////////////
    void _partition();
    int _getHost(tetrahedron_global_id tgidx) const;
    int _getHost(triangle_global_id tgidx) const;
    ////////////////////////////////////////////////////////////////////////

    // Keeps track of whether _build() has been called
    // bool                                       pBuilt;

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    // The Efield flag. If false we don't calclulate the potential, nor include
    // any voltage-dependent transitions, ohmic or ghk currents. This means
    // the solver behaves exactly like the previous TetOpSplitP solver in this
    // case and all following members are set to null pointer or zero.
    //
    EF_solver pEFoption;

    double pTemp{0.0};

    // Pointer to the EField object
    std::unique_ptr<solver::efield::EField> pEField;

    // The Efield time-step
    double pEFDT{1.0e-5};

    // The number of vertices
    uint pEFNVerts{0};

    // The number of membrane triangles
    uint pEFNTris{0};

    std::vector<Tri*> pEFTris_vec;

    util::strongid_vector<triangle_local_id, double> EFTrisV;

    // Working space for gathering distributed computed triangle currents,
    // grouped by rank of owner.
    std::vector<double> EFTrisI_permuted;

    // Translate from permuted vector of triangle currents to local EFTri indices.
    std::vector<triangle_local_id> EFTrisI_idx;

    // Per-rank counts of EFTris.
    std::vector<int> EFTrisI_count;

    // Offsets into ETris_permuted where a rank's owned EFTri data is stored.
    std::vector<int> EFTrisI_offset;

    // The number of tetrahedrons
    uint pEFNTets{0};
    // Array of tetrahedrons

    // Table of global vertex index to EField local vertex index (0, 1, ...,
    // pEFNVerts - 1)
    util::strongid_vector<vertex_id_t, vertex_id_t> pEFVert_GtoL;

    // Table of global triangle index to EField local triangle index (0, 1, ...,
    // pEFNTris-1)
    util::strongid_vector<triangle_global_id, triangle_local_id> pEFTri_GtoL;

    // Table of global tetrahedron index to EField local tet index (0, 1, ...,
    // pEFNTets-1)
    util::strongid_vector<tetrahedron_global_id, tetrahedron_local_id> pEFTet_GtoL;

    // Table of EField local triangle index to global triangle index.
    util::strongid_vector<triangle_local_id, triangle_global_id> pEFTri_LtoG;

    ////////////////////////// MPI STUFFS ////////////////////////////

    std::map<tetrahedron_global_id, int> tetHosts;
    std::map<triangle_global_id, int> triHosts;
    util::strongid_vector<solver::comp_global_id, int> wmHosts;
    int myRank{};
    int nHosts{};
    bool recomputeUpdPeriod{true};
    unsigned long long reacExtent{0};
    unsigned long long diffExtent{0};
    double nIteration{0.0};
    double updPeriod{0.0};
    // bool                                        requireSync;
    uint diffApplyThreshold{10};

    std::set<int> neighbHosts;
    uint nNeighbHosts{};

    // mirrors of neighboring tets
    std::set<Tet*> boundaryTets;
    std::set<Tri*> boundaryTris;

    std::map<int, std::vector<uint>> remoteChanges;

    void _remoteSyncAndUpdate(void* requests,
                              std::vector<KProc*>& applied_diffs,
                              std::vector<int>& directions);

    // void _applyRemoteMoleculeChanges(std::vector<MPI_Request> & requests);
    // void _syncPoolCounts();
    // void _updateKProcRates(std::vector<KProc*> & applylist, std::vector<int> &
    // directions, std::vector<MPI_Request> & requests);
    // send, receive and process remote upadtes
    // void _updateRemoteKProcRates(std::vector<MPI_Request> & requests);

    // STL random number generator - also Mersenne twister
    std::random_device rd;
    std::mt19937 gen;

    double compTime{0.0};
    double syncTime{0.0};
    double idleTime{0.0};
    double efieldTime{0.0};
    double rdTime{0.0};
    double dataExchangeTime{0.0};
};

}  // namespace steps::mpi::tetopsplit
