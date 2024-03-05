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

#include <iostream>
#include <map>
#include <memory>
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

namespace steps::tetexact {

// Auxiliary declarations.
typedef solver::kproc_global_id SchedIDX;
typedef std::set<SchedIDX> SchedIDXSet;
typedef SchedIDXSet::iterator SchedIDXSetI;
typedef SchedIDXSet::const_iterator SchedIDXSetCI;
typedef std::vector<SchedIDX> SchedIDXVec;
typedef SchedIDXVec::iterator SchedIDXVecI;
typedef SchedIDXVec::const_iterator SchedIDXVecCI;

/// Copies the contents of a set of SchedIDX entries into a vector.
/// The contents of the vector are completely overridden.
///
extern void schedIDXSet_To_Vec(SchedIDXSet const& s, SchedIDXVec& v);

////////////////////////////////////////////////////////////////////////////////

class Tetexact: public solver::API {
  public:
    Tetexact(model::Model* m, wm::Geom* g, const rng::RNGptr& r, int calcMembPot = EF_NONE);
    ~Tetexact() override;

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
    // void advanceSteps(uint nsteps);
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
    // Batch Data Access
    ////////////////////////////////////////////////////////////////////////

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

    double _getCompReacK(solver::comp_global_id cidx, solver::reac_global_id sidx) const override;
    void _setCompReacK(solver::comp_global_id cidx,
                       solver::reac_global_id sidx,
                       double kf) override;

    bool _getCompReacActive(solver::comp_global_id cidx,
                            solver::reac_global_id sidx) const override;
    void _setCompReacActive(solver::comp_global_id cidx,
                            solver::reac_global_id sidx,
                            bool a) override;

    double _getCompReacH(solver::comp_global_id cidx, solver::reac_global_id sidx) const override;
    double _getCompReacC(solver::comp_global_id cidx, solver::reac_global_id sidx) const override;
    long double _getCompReacA(solver::comp_global_id cidx,
                              solver::reac_global_id sidx) const override;

    unsigned long long _getCompReacExtent(solver::comp_global_id cidx,
                                          solver::reac_global_id sidx) const override;
    void _resetCompReacExtent(solver::comp_global_id cidx, solver::reac_global_id sidx) override;

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
                           solver::sreac_global_id sidx) const override;
    void _setPatchSReacK(solver::patch_global_id pidx,
                         solver::sreac_global_id sidx,
                         double kf) override;

    bool _getPatchSReacActive(solver::patch_global_id pidx,
                              solver::sreac_global_id sidx) const override;
    void _setPatchSReacActive(solver::patch_global_id pidx,
                              solver::sreac_global_id sidx,
                              bool a) override;

    double _getPatchSReacH(solver::patch_global_id pidx,
                           solver::sreac_global_id sidx) const override;
    double _getPatchSReacC(solver::patch_global_id pidx,
                           solver::sreac_global_id sidx) const override;
    double _getPatchSReacA(solver::patch_global_id pidx,
                           solver::sreac_global_id sidx) const override;

    unsigned long long _getPatchSReacExtent(solver::patch_global_id pidx,
                                            solver::sreac_global_id sidx) const override;
    void _resetPatchSReacExtent(solver::patch_global_id pidx,
                                solver::sreac_global_id sidx) override;

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

    double _getTetReacK(tetrahedron_global_id tidx, solver::reac_global_id sidx) const override;
    void _setTetReacK(tetrahedron_global_id tidx, solver::reac_global_id sidx, double kf) override;

    bool _getTetReacActive(tetrahedron_global_id tidx, solver::reac_global_id sidx) const override;
    void _setTetReacActive(tetrahedron_global_id tidx,
                           solver::reac_global_id sidx,
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

    double _getTetReacH(tetrahedron_global_id tidx, solver::reac_global_id sidx) const override;
    double _getTetReacC(tetrahedron_global_id tidx, solver::reac_global_id sidx) const override;
    double _getTetReacA(tetrahedron_global_id tidx, solver::reac_global_id sidx) const override;

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

    double _getTriSReacK(triangle_global_id tidx, solver::sreac_global_id sidx) const override;
    void _setTriSReacK(triangle_global_id tidx, solver::sreac_global_id sidx, double kf) override;

    bool _getTriSReacActive(triangle_global_id tidx, solver::sreac_global_id sidx) const override;
    void _setTriSReacActive(triangle_global_id tidx,
                            solver::sreac_global_id sidx,
                            bool act) override;

    double _getTriSDiffD(triangle_global_id tidx,
                         solver::surfdiff_global_id didx,
                         triangle_global_id direction_tri) const override;

    void _setTriSDiffD(triangle_global_id tidx,
                       solver::surfdiff_global_id didx,
                       double dk,
                       triangle_global_id direction_tri) override;

    ////////////////////////////////////////////////////////////////////////

    double _getTriSReacH(triangle_global_id tidx, solver::sreac_global_id sidx) const override;
    double _getTriSReacC(triangle_global_id tidx, solver::sreac_global_id sidx) const override;
    double _getTriSReacA(triangle_global_id tidx, solver::sreac_global_id sidx) const override;

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

    // Called from local Comp or Patch objects. Add KProc to this object
    void addKProc(KProc* kp, bool Vdep = false);

    inline uint countKProcs() const {
        return pKProcs.size();
    }

    ////////////////////////////////////////////////////////////////////////

    inline const tetmesh::Tetmesh& mesh() const noexcept {
        return *pMesh;
    }

    inline tetmesh::Tetmesh& mesh() noexcept {
        return *pMesh;
    }

    inline Comp* _comp(solver::comp_global_id cidx) const {
        // Moved assertions to this access routine to cut code duplication.
        AssertLog(cidx < statedef().countComps());
        AssertLog(statedef().countComps() == pComps.size());

        auto c = pComps[cidx];
        // AssertLog(c !=0);
        return c;
    }

    inline const auto& patches() const {
        return pPatches;
    }

    inline Patch* _patch(solver::patch_global_id pidx) const {
        // Moved assertions to this access routine to cut code duplication.
        AssertLog(pidx < statedef().countPatches());
        AssertLog(statedef().countPatches() == pPatches.size());

        auto p = pPatches[pidx];
        AssertLog(p != 0);
        return p;
    }

    inline DiffBoundary* _diffboundary(solver::diffboundary_global_id dbidx) const {
        AssertLog(dbidx < statedef().countDiffBoundaries());
        return pDiffBoundaries[dbidx.get()];
    }

    inline SDiffBoundary* _sdiffboundary(solver::sdiffboundary_global_id sdbidx) const {
        AssertLog(sdbidx < statedef().countSDiffBoundaries());
        return pSDiffBoundaries[sdbidx.get()];
    }

    inline double a0() const {
        return pA0;
    }

    // Checked global to local index translations

    solver::spec_local_id specG2L_or_throw(Comp* comp, solver::spec_global_id gidx) const;
    solver::spec_local_id specG2L_or_throw(Patch* patch, solver::spec_global_id gidx) const;
#if 0
    solver::spec_local_id specG2L_or_throw(Tet *tet, solver::spec_global_id gidx) const;
    solver::spec_local_id specG2L_or_throw(Tri *tri, solver::spec_global_id gidx) const;
#endif
    solver::reac_local_id reacG2L_or_throw(Comp* comp, solver::reac_global_id gidx) const;
    solver::sreac_local_id sreacG2L_or_throw(Patch* patch, solver::sreac_global_id gidx) const;
    solver::diff_local_id diffG2L_or_throw(Comp* comp, solver::diff_global_id gidx) const;
    solver::surfdiff_local_id sdiffG2L_or_throw(Patch* patch,
                                                solver::surfdiff_global_id gidx) const;
    solver::vdepsreac_local_id vdepsreacG2L_or_throw(Patch* patch,
                                                     solver::vdepsreac_global_id gidx) const;

    ////////////////////////////////////////////////////////////////////////
    // TETEXACT SOLVER METHODS
    ////////////////////////////////////////////////////////////////////////

    solver::comp_global_id _addComp(solver::Compdef* cdef);

    solver::patch_global_id _addPatch(solver::Patchdef* pdef);

    solver::diffboundary_global_id _addDiffBoundary(solver::DiffBoundarydef* dbdef);

    solver::sdiffboundary_global_id _addSDiffBoundary(solver::SDiffBoundarydef* sdbdef);

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

    // void _build();

    double _getRate(uint i) const {
        return pKProcs[i]->rate();
    }

    KProc* _getNext() const;

    // void _reset();

    void _executeStep(KProc* kp, double dt);

    // TODO: Change the following so that only the kprocs depending on
    // the species are updated. These functions are called from interface
    // methods setting compartment or patch counts.
    /// Update the kproc's of a tet, after a species has been changed.
    /// This also updates kproc's in surrounding triangles.
    ///
    /// Currently doesn't care about the species.
    ///
    void _updateSpec(WmVol& tet);

    /// Update the kproc's of a triangle, after a species has been changed.
    /// This does not need to update the kproc's of any neighbouring
    /// tetrahedrons.
    ///
    /// Currently doesn't care about the species.
    ///
    void _updateSpec(Tri& tri);

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    /// Check the EField flag
    inline bool efflag() const {
        return pEFoption != EF_NONE;
    }

    void _setupEField();

    inline uint neftets() const {
        return pEFNTets;
    }

    inline uint neftris() const {
        return pEFNTris;
    }

    inline uint nefverts() const {
        return pEFNVerts;
    }

    ////////////////////////////////////////////////////////////////////////
  private:
    ////////////////////////////////////////////////////////////////////////

    double getROITetSpecCount(const std::vector<tetrahedron_global_id>& indices,
                              const std::string& s) const;
    double getROITriSpecCount(const std::vector<triangle_global_id>& indices,
                              const std::string& s) const;

    void setROITetSpecClamped(const std::vector<tetrahedron_global_id>& indices,
                              std::string const& s,
                              bool b);
    void setROITriSpecClamped(const std::vector<triangle_global_id>& indices,
                              std::string const& s,
                              bool b);

    void setROITriSpecCount(const std::vector<triangle_global_id>& indices,
                            std::string const& s,
                            double count);
    void setROITetSpecCount(const std::vector<tetrahedron_global_id>& indices,
                            std::string const& s,
                            double count);

    double getROIVol(const std::vector<tetrahedron_global_id>& tets) const;

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
    // CR SSA Kernel Data and Methods
    ////////////////////////////////////////////////////////////////////////
    std::size_t nEntries{};
    double pSum{};
    double nSum{};
    double pA0{0.0};

    std::vector<KProc*> pKProcs;
    std::vector<KProc*> pVdepKProcs;
    std::vector<CRGroup*> nGroups;
    std::vector<CRGroup*> pGroups;

    ////////////////////////////////////////////////////////////////////////////////

    template <typename KProcPIter>
    inline void _update(KProcPIter b, KProcPIter e) {
        while (b != e)
            _updateElement(**b++);
        _updateSum();
    }

    ////////////////////////////////////////////////////////////////////////////////

    inline void _update() {
        _update(pKProcs.begin(), pKProcs.end());
    }

    ////////////////////////////////////////////////////////////////////////////////

    inline CRGroup* _getGroup(int pow) {
#ifdef SSA_DEBUG
        CLOG(INFO, "general_log") << "SSA: get group with power " << pow << "\n";
#endif

        if (pow >= 0) {
#ifdef SSA_DEBUG
            CLOG(INFO, "general_log") << "positive group" << pow << "\n";
#endif
            return pGroups[pow];
        } else {
#ifdef SSA_DEBUG
            CLOG(INFO, "general_log") << "negative group" << -pow << "\n";
#endif
            return nGroups[-pow];
        }
#ifdef SSA_DEBUG
        CLOG(INFO, "general_log") << "--------------------------------------------------------\n";
#endif
    }

    ////////////////////////////////////////////////////////////////////////////////

    inline void _extendPGroups(uint new_size) {
        uint curr_size = pGroups.size();

#ifdef SSA_DEBUG
        CLOG(INFO, "general_log") << "SSA: extending positive group size to " << new_size;
        CLOG(INFO, "general_log") << " from " << curr_size << ".\n";
        CLOG(INFO, "general_log") << "--------------------------------------------------------\n";
#endif

        while (curr_size < new_size) {
            pGroups.push_back(new CRGroup(curr_size));
            curr_size++;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////

    inline void _extendNGroups(uint new_size) {
        uint curr_size = nGroups.size();

#ifdef SSA_DEBUG
        CLOG(INFO, "general_log") << "SSA: extending negative group size to " << new_size;
        CLOG(INFO, "general_log") << " from " << curr_size << ".\n";
        CLOG(INFO, "general_log") << "--------------------------------------------------------\n";
#endif

        while (curr_size < new_size) {
            nGroups.push_back(new CRGroup(-curr_size));
            curr_size++;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////

    inline void _extendGroup(CRGroup* group, uint size = 1024) {
#ifdef SSA_DEBUG
        CLOG(INFO, "general_log") << "SSA: extending group storage\n";
        CLOG(INFO, "general_log") << "current capacity: " << group->capacity << "\n";
#endif

        group->capacity += size;
        group->indices = static_cast<KProc**>(
            realloc(group->indices, sizeof(KProc*) * group->capacity));
        if (group->indices == nullptr) {
            SysErrLog("DirectCR: unable to allocate memory for SSA group.");
        }
#ifdef SSA_DEBUG
        CLOG(INFO, "general_log") << "capacity after extending: " << group->capacity << "\n";
        CLOG(INFO, "general_log") << "--------------------------------------------------------\n";
#endif
    }

    ////////////////////////////////////////////////////////////////////////////////

    void _updateElement(KProc& kp);

    inline void _updateSum() {
#ifdef SSA_DEBUG
        CLOG(INFO, "general_log") << "update A0 from " << pA0 << " to ";
#endif

        pA0 = 0.0;

        for (const auto& neg_grp: nGroups) {
            pA0 += neg_grp->sum;
        }
        for (const auto& pos_grp: pGroups) {
            pA0 += pos_grp->sum;
        }

#ifdef SSA_DEBUG
        CLOG(INFO, "general_log") << pA0 << "\n";
#endif
    }

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    // The Efield solve choise. If EF_NONE we don't calclulate the potential, nor
    // include any voltage-dependent transitions, ohmic or ghk currents. This
    // means the solver behaves exactly like the previous Tetexact solver in this
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
    // Array of vertices
    std::vector<double> pEFVerts;

    // The number of membrane triangles
    uint pEFNTris{0};
    // Array of membrane triangles
    std::vector<vertex_id_t> pEFTris;

    std::vector<Tri*> pEFTris_vec;

    // The number of tetrahedrons
    uint pEFNTets{0};
    // Array of tetrahedrons
    std::vector<vertex_id_t> pEFTets;

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
};

}  // namespace steps::tetexact
