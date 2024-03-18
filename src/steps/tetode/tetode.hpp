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
#include <memory>
#include <string>
#include <vector>

// STEPS headers.
#include "comp.hpp"
#include "patch.hpp"
#include "tet.hpp"
#include "tri.hpp"

#include "geom/tetmesh.hpp"
#include "solver/api.hpp"
#include "solver/efield/efield.hpp"
#include "solver/statedef.hpp"

namespace steps::tetode {


// These structures store the reaction information that can be used for every //
// species to work out it's dy/dt. This replaces a sparse matrix with a vector /
// of length of the total number of species per tetrahedron, each individual ///
// species then storing a vector of 'structA's which represent the specific ////
// represented reaction for that species in that tetrahedron. //////////////////
// The 'species index' in struct C is the 'vector index' and the total number //
// of species will be equal to sum_N: comp_n * tetsincomp_n ////////////////////

// strcutCs will store the 'vector idx' of reactant species involved in reaction
// which will dictate the rate at which the species varies by this reaction.
// Order is per spec order, so for hypothetical reaction [A+A+B-> something]
// two structCs will be built, one for A (with order 2) and one for B (with
// order 1)

struct structC {
    uint order;
    uint spec_idx;
};

// structBs are used to find the per-reaction lhs
// Simply stores a vector of structCs because there are often more than
// one reactant species involved in a reaction

struct structB {
    std::vector<structC> info;
};

// structAs used to find the per-spec reaction info
// it stores the local, microscopic, c-constant, wich will of coursr
// vary depending where in the mesh this rection occurs.
// Will also store the 'vector reaction index', which will be in terms
// of all possible reactions in the system
// The 'update value' is info about how much of a species is created
// or destroyed in this reaction, e.g. [A->B] the upd for A is -1 and for
// B is +1 for this reaction. A structA will be created for all species
// for all reactions where it changes, i.e. the upd value is non-zero.

struct structA {
    // A reaction ccst, will differ from tet to tet
    double ccst;

    // Now some matrix reaction index will be stored to allow for
    // changing the kcst (and so ccst)
    uint r_idx;

    // The update value
    int upd;

    // Structure of species 'players' information
    std::vector<structB> players;
};

// Keep CVode data structures and definitions internal
struct CVodeState;

class TetODE: public solver::API {
  public:
    TetODE(model::Model* m, wm::Geom* g, const rng::RNGptr& r, int calcMembPot = EF_NONE);
    ~TetODE() override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER INFORMATION
    ////////////////////////////////////////////////////////////////////////

    std::string getSolverName() const override;
    std::string getSolverDesc() const override;
    std::string getSolverAuthors() const override;
    std::string getSolverEmail() const override;

    void checkpoint(std::string const& file_name) override;

    void restore(std::string const& file_name) override;

    double getTime() const override;

    inline double getTemp() const override {
        return pTemp;
    }

    void setTemp(double t) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS
    ////////////////////////////////////////////////////////////////////////

    // Will of course be added- just for compilation purposes
    void reset() override;
    void run(double endtime) override;
    void advance(double adv) override;
    // void step();

    ////////////////////////////////////////////////////////////////////////

    inline tetmesh::Tetmesh* mesh() const {
        return pMesh;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

    double _getCompSpecCount(solver::comp_global_id cidx,
                             solver::spec_global_id sidx) const override;
    void _setCompSpecCount(solver::comp_global_id cidx,
                           solver::spec_global_id sidx,
                           double n) override;

    double _getCompVol(solver::comp_global_id cidx) const override;

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

    /*
   unsigned long long _getCompReacExtent(solver::comp_global_id cidx, uint ridx) const;
   void _resetCompReacExtent(solver::comp_global_id cidx, uint ridx);

   double _getCompDiffD(solver::comp_global_id cidx, uint didx) const;
   void _setCompDiffD(solver::comp_global_id cidx, uint didx, double dk);

   bool _getCompDiffActive(solver::comp_global_id cidx, uint didx) const;
   void _setCompDiffActive(solver::comp_global_id cidx, uint didx, bool act);

   ////////////////////////////////////////////////////////////////////////
   // SOLVER STATE ACCESS:
   //      PATCH
   ////////////////////////////////////////////////////////////////////////

   */
    double _getPatchSpecCount(solver::patch_global_id pidx,
                              solver::spec_global_id sidx) const override;
    void _setPatchSpecCount(solver::patch_global_id pidx,
                            solver::spec_global_id sidx,
                            double n) override;

    double _getPatchArea(solver::patch_global_id pidx) const override;

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

    /*
    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TETRAHEDRAL VOLUME ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    void _setTetVol(uint tidx, double vol);

    bool _getTetSpecDefined(uint tidx, solver::spec_global_id sidx) const;

    */
    double _getTetVol(tetrahedron_global_id tidx) const override;

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

    double _getTetReacK(tetrahedron_global_id tidx, solver::reac_global_id ridx) const override;
    void _setTetReacK(tetrahedron_global_id tidx, solver::reac_global_id ridx, double kf) override;

    /*
    bool _getTetSpecClamped(uint tidx, solver::spec_global_id sidx) const;
    void _setTetSpecClamped(uint tidx, solver::spec_global_id sidx, bool buf);

    bool _getTetReacActive(uint tidx, uint ridx) const;
    void _setTetReacActive(uint tidx, uint ridx, bool act);

    double _getTetDiffD(uint tidx, uint didx) const;
    void _setTetDiffD(uint tidx, uint didx, double dk);

    bool _getTetDiffActive(uint tidx, uint didx) const;
    void _setTetDiffActive(uint tidx, uint didx, bool act);
    */

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    double _getTetV(tetrahedron_global_id tidx) const override;
    void _setTetV(tetrahedron_global_id tidx, double v) override;
    bool _getTetVClamped(tetrahedron_global_id tidx) const override;
    void _setTetVClamped(tetrahedron_global_id tidx, bool cl) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TRIANGULAR SURFACE ELEMENTS
    ////////////////////////////////////////////////////////////////////////
    /*
    void _setTriArea(uint tidx, double area);

    bool _getTriSpecDefined(uint tidx, solver::spec_global_id sidx) const;

    */

    double _getTriArea(triangle_global_id tidx) const override;

    double _getTriSpecCount(triangle_global_id tidx, solver::spec_global_id sidx) const override;
    void _setTriSpecCount(triangle_global_id tidx, solver::spec_global_id sidx, double n) override;

    double _getTriSpecAmount(triangle_global_id tidx, solver::spec_global_id sidx) const override;
    void _setTriSpecAmount(triangle_global_id tidx, solver::spec_global_id sidx, double m) override;

    double _getTriSReacK(triangle_global_id tidx, solver::sreac_global_id ridx) const override;
    void _setTriSReacK(triangle_global_id tidx, solver::sreac_global_id ridx, double kf) override;

    /*
    bool _getTriSpecClamped(uint tidx, solver::spec_global_id sidx) const;
    void _setTriSpecClamped(uint tidx, solver::spec_global_id sidx, bool buf);

    bool _getTriSReacActive(uint tidx, uint ridx) const;
    void _setTriSReacActive(uint tidx, uint ridx, bool act);

     */

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

    /*
    double _getTriOhmicI(uint tidx) const;
    double _getTriOhmicI(uint tidx, uint ocidx) const;

    double _getTriGHKI(uint tidx) const;
    double _getTriGHKI(uint tidx, uint ghkidx) const;
    */

    double _getTriI(triangle_global_id tidx) const override;

    void _setTriIClamp(triangle_global_id tidx, double cur) override;

    ////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      VERTICES ELEMENTS
    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    double _getVertV(vertex_id_t vidx) const override;
    void _setVertV(vertex_id_t vidx, double v) override;
    bool _getVertVClamped(vertex_id_t vidx) const override;
    void _setVertVClamped(vertex_id_t vidx, bool cl) override;
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

    void _setup();

    std::size_t _addComp(solver::Compdef* cdef);

    std::size_t _addPatch(solver::Patchdef* pdef);

    void _addTet(tetrahedron_global_id tetidx,
                 steps::tetode::Comp* comp,
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

    void _addTri(triangle_global_id triidx,
                 steps::tetode::Patch* patch,
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

    /// returns properly scaled reaction constant
    ///
    double _ccst(double kcst, double vol, uint order);

    /// returns properly scaled reaction constant for surface-surface case
    ///
    double _ccst2D(double kcst, double area, uint order);

    ////////////////////////////////////////////////////////////////////////
    // CVODE FUNCTIONS
    ////////////////////////////////////////////////////////////////////////

    void setTolerances(double atol, double rtol);

    void setMaxNumSteps(uint maxn);

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

  private:
    ////////////////////////////////////////////////////////////////////////

    tetmesh::Tetmesh* pMesh{nullptr};

    util::strongid_vector<solver::comp_global_id, steps::tetode::Comp*> pComps;
    inline Comp* comps(solver::comp_global_id cidx) const noexcept {
        return pComps[cidx];
    }

    util::strongid_vector<solver::patch_global_id, steps::tetode::Patch*> pPatches;
    inline Patch* patches(solver::patch_global_id pidx) const noexcept {
        return pPatches[pidx];
    }

    std::vector<steps::tetode::Tri*> pTris;

    // Now stored as base pointer
    std::vector<steps::tetode::Tet*> pTets;

    uint pSpecs_tot{0};
    uint pReacs_tot{0};


    bool pInitialised{false};
    bool pTolsset{false};
    bool pReinit{true};

    CVodeState* pCVodeState{nullptr};

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    // The Efield flag. If false we don't calclulate the potential, nor include
    // any voltage-dependent transitions, ohmic or ghk currents. This means
    // the solver behaves exactly like the previous Tetexact solver in this case
    // and all following members are set to null pointer or zero.
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

    std::vector<steps::tetode::Tri*> pEFTris_vec;

    // The number of tetrahedrons
    uint pEFNTets{0};
    // Array of tetrahedrons
    std::vector<vertex_id_t> pEFTets;

    // Table of global vertex index to EField local vertex index (0, 1, ...,
    // pEFNVerts - 1)
    std::vector<vertex_id_t> pEFVert_GtoL;

    // Table of global triangle index to EField local triangle index (0, 1, ...,
    // pEFNTris-1)
    std::vector<triangle_local_id> pEFTri_GtoL;

    // Table of global tetrahedron index to EField local tet index (0, 1, ...,
    // pEFNTets-1)
    std::vector<tetrahedron_local_id> pEFTet_GtoL;

    // Table of EField local triangle index to global triangle index.
    std::vector<triangle_global_id> pEFTri_LtoG;
};

}  // namespace steps::tetode
