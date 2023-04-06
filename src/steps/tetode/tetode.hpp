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


#ifndef STEPS_SOLVER_TETODE_HPP
#define STEPS_SOLVER_TETODE_HPP 1


// STL headers.
#include <memory>
#include <string>
#include <vector>

// STEPS headers.
#include "comp.hpp"
#include "patch.hpp"
#include "tet.hpp"
#include "tri.hpp"

#include "util/common.h"
#include "geom/tetmesh.hpp"
#include "solver/api.hpp"
#include "solver/statedef.hpp"
#include "solver/efield/efield.hpp"


////////////////////////////////////////////////////////////////////////////////

namespace steps {
namespace tetode {


////////////////////////////////////////////////////////////////////////////////
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

struct structC
{
    uint order;
    uint spec_idx;
};


// structBs are used to find the per-reaction lhs
// Simply stores a vector of structCs because there are often more than
// one reactant species involved in a reaction

struct structB
{
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

struct structA
{
    // A reaction ccst, will differ from tet to tet
    double  ccst;

    // Now some matrix reaction index will be stored to allow for
    // changing the kcst (and so ccst)
    uint r_idx;

    // The update value
    int upd;

    // Structure of species 'players' information
    std::vector<structB> players;
};

}
}

////////////////////////////////////////////////////////////////////////////////

namespace steps {
namespace tetode {
using namespace steps::solver;

////////////////////////////////////////////////////////////////////////////////

// Keep CVode data structures and definitions internal
struct CVodeState;

class TetODE: public API
{

public:

    TetODE(steps::model::Model *m, steps::wm::Geom *g, const rng::RNGptr &r,
           int calcMembPot = EF_NONE);
    ~TetODE() override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER INFORMATION
    ////////////////////////////////////////////////////////////////////////

    std::string getSolverName() const override;
    std::string getSolverDesc() const override;
    std::string getSolverAuthors() const override;
    std::string getSolverEmail() const override;


    void checkpoint(std::string const & file_name) override;

    void restore(std::string const & file_name) override;

    double getTime() const override;

    inline double getTemp() const override
    { return pTemp; }

    void setTemp(double t) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS
    ////////////////////////////////////////////////////////////////////////

    // Will of course be added- just for compilation purposes
    void reset() override;
    void run(double endtime) override;
    void advance(double adv) override;
    //void step();

    ////////////////////////////////////////////////////////////////////////

    inline steps::tetmesh::Tetmesh * mesh() const
    { return pMesh; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

     double _getCompCount(uint cidx, uint sidx) const override;
     void _setCompCount(uint cidx, uint sidx, double n) override;

    double _getCompVol(uint cidx) const override;

    double _getCompAmount(uint cidx, uint sidx) const override;
    void _setCompAmount(uint cidx, uint sidx, double a) override;

    double _getCompConc(uint cidx, uint sidx) const override;
    void _setCompConc(uint cidx, uint sidx, double c) override;

    bool _getCompClamped(uint cidx, uint sidx) const override;
    void _setCompClamped(uint cidx, uint sidx, bool b) override;

    double _getCompReacK(uint cidx, uint ridx) const override;
    void _setCompReacK(uint cidx, uint ridx, double kf) override;

    bool _getCompReacActive(uint cidx, uint ridx) const override;
    void _setCompReacActive(uint cidx, uint ridx, bool a) override;

     /*
    unsigned long long _getCompReacExtent(uint cidx, uint ridx) const;
    void _resetCompReacExtent(uint cidx, uint ridx);

    double _getCompDiffD(uint cidx, uint didx) const;
    void _setCompDiffD(uint cidx, uint didx, double dk);

    bool _getCompDiffActive(uint cidx, uint didx) const;
    void _setCompDiffActive(uint cidx, uint didx, bool act);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      PATCH
    ////////////////////////////////////////////////////////////////////////

    */
     double _getPatchCount(uint pidx, uint sidx) const override;
    void _setPatchCount(uint pidx, uint sidx, double n) override;

    double _getPatchArea(uint pidx) const override;

    double _getPatchAmount(uint pidx, uint sidx) const override;
    void _setPatchAmount(uint pidx, uint sidx, double a) override;

    bool _getPatchClamped(uint pidx, uint sidx) const override;
    void _setPatchClamped(uint pidx, uint sidx, bool buf) override;

    double _getPatchSReacK(uint pidx, uint ridx) const override;
    void _setPatchSReacK(uint pidx, uint ridx, double kf) override;

    bool _getPatchSReacActive(uint pidx, uint ridx) const override;
    void _setPatchSReacActive(uint pidx, uint ridx, bool a) override;

    /*
    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TETRAHEDRAL VOLUME ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    void _setTetVol(uint tidx, double vol);

    bool _getTetSpecDefined(uint tidx, uint sidx) const;

    */
    double _getTetVol(tetrahedron_id_t tidx) const override;

    double _getTetCount(tetrahedron_id_t tidx, uint sidx) const override;
    void _setTetCount(tetrahedron_id_t tidx, uint sidx, double n) override;

    double _getTetAmount(tetrahedron_id_t tidx, uint sidx) const override;
    void _setTetAmount(tetrahedron_id_t tidx, uint sidx, double m) override;

    double _getTetConc(tetrahedron_id_t tidx, uint sidx) const override;
    void _setTetConc(tetrahedron_id_t tidx, uint sidx, double c) override;

    double _getTetReacK(tetrahedron_id_t tidx, uint ridx) const override;
    void _setTetReacK(tetrahedron_id_t tidx, uint ridx, double kf) override;

    /*
    bool _getTetClamped(uint tidx, uint sidx) const;
    void _setTetClamped(uint tidx, uint sidx, bool buf);

    bool _getTetReacActive(uint tidx, uint ridx) const;
    void _setTetReacActive(uint tidx, uint ridx, bool act);

    double _getTetDiffD(uint tidx, uint didx) const;
    void _setTetDiffD(uint tidx, uint didx, double dk);

    bool _getTetDiffActive(uint tidx, uint didx) const;
    void _setTetDiffActive(uint tidx, uint didx, bool act);
    */

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    double _getTetV(tetrahedron_id_t tidx) const override;
    void _setTetV(tetrahedron_id_t tidx, double v) override;
    bool _getTetVClamped(tetrahedron_id_t tidx) const override;
    void _setTetVClamped(tetrahedron_id_t tidx, bool cl) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TRIANGULAR SURFACE ELEMENTS
    ////////////////////////////////////////////////////////////////////////
    /*
    void _setTriArea(uint tidx, double area);

    bool _getTriSpecDefined(uint tidx, uint sidx) const;

    */

    double _getTriArea(triangle_id_t tidx) const override;

    double _getTriCount(triangle_id_t tidx, uint sidx) const override;
    void _setTriCount(triangle_id_t tidx, uint sidx, double n) override;

    double _getTriAmount(triangle_id_t tidx, uint sidx) const override;
    void _setTriAmount(triangle_id_t tidx, uint sidx, double m) override;

    double _getTriSReacK(triangle_id_t tidx, uint ridx) const override;
    void _setTriSReacK(triangle_id_t tidx, uint ridx, double kf) override;

    /*
    bool _getTriClamped(uint tidx, uint sidx) const;
    void _setTriClamped(uint tidx, uint sidx, bool buf);

    bool _getTriSReacActive(uint tidx, uint ridx) const;
    void _setTriSReacActive(uint tidx, uint ridx, bool act);

     */

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    double _getTriV(triangle_id_t tidx) const override;
    void _setTriV(triangle_id_t tidx, double v) override;
    bool _getTriVClamped(triangle_id_t tidx) const override;
    void _setTriVClamped(triangle_id_t tidx, bool cl) override;

    /*
    double _getTriOhmicI(uint tidx) const;
    double _getTriOhmicI(uint tidx, uint ocidx) const;

    double _getTriGHKI(uint tidx) const;
    double _getTriGHKI(uint tidx, uint ghkidx) const;
    */

    double _getTriI(triangle_id_t tidx) const override;

    void _setTriIClamp(triangle_id_t tidx, double cur) override;

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

    void _setMembPotential(uint midx, double v) override;
    void _setMembCapac(uint midx, double cm) override;
    void _setMembVolRes(uint midx, double ro) override;
    void _setMembRes(uint midx, double ro, double vrev) override;

    ////////////////////////////////////////////////////////////////////////


    void _setup();

    std::size_t _addComp(steps::solver::Compdef * cdef);

    std::size_t _addPatch(steps::solver::Patchdef * pdef);

    void _addTet(tetrahedron_id_t tetidx,
                 steps::tetode::Comp *comp,
                 double vol,
                 double a1,
                 double a2,
                 double a3,
                 double a4,
                 double d1,
                 double d2,
                 double d3,
                 double d4,
                 tetrahedron_id_t tet0,
                 tetrahedron_id_t tet1,
                 tetrahedron_id_t tet2,
                 tetrahedron_id_t tet3);

    void _addTri(triangle_id_t triidx,
                 steps::tetode::Patch *patch,
                 double area,
                 double l0,
                 double l1,
                 double l2,
                 double d0,
                 double d1,
                 double d2,
                 tetrahedron_id_t tinner,
                 tetrahedron_id_t touter,
                 triangle_id_t tri0,
                 triangle_id_t tri1,
                 triangle_id_t tri2);

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
    inline bool efflag() const noexcept
    { return pEFoption != EF_NONE; }

    void _setupEField();

    inline uint neftets() const noexcept
    { return pEFNTets; }

    inline uint neftris() const noexcept
    { return pEFNTris; }

    inline uint nefverts() const noexcept
    { return pEFNVerts; }


private:

    ////////////////////////////////////////////////////////////////////////

    steps::tetmesh::Tetmesh *                  pMesh{nullptr};

    std::vector<steps::tetode::Comp *>       pComps;
    //std::map<steps::solver::Compdef *, Comp *> pCompMap;

    std::vector<steps::tetode::Patch *>      pPatches;

    std::vector<steps::tetode::Tri *>        pTris;

    // Now stored as base pointer
    std::vector<steps::tetode::Tet *>        pTets;

    /*
    // A vector all the reaction information, hopefully ingeniously
    // removing the need for a sparse matrix at all
    std::vector<std::vector<steps::tetode::structA> >  pSpec_matrixsub;
    */

    uint                                      pSpecs_tot{0};
    uint                                      pReacs_tot{0};

    //double                                    * pCcst;

    bool                                      pInitialised{false};
    bool                                     pTolsset{false};
    bool                                      pReinit{true};

    CVodeState                              * pCVodeState{nullptr};

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    // The Efield flag. If false we don't calclulate the potential, nor include
    // any voltage-dependent transitions, ohmic or ghk currents. This means
    // the solver behaves exactly like the previous Tetexact solver in this case
    // and all following members are set to null pointer or zero.
    //
    EF_solver                                     pEFoption;

    double                                        pTemp{0.0};

    // Pointer to the EField object
    std::unique_ptr<steps::solver::efield::EField> pEField;

    // The Efield time-step
    double                                       pEFDT{1.0e-5};

    // The number of vertices
    uint                                        pEFNVerts{0};
    // Array of vertices
    double                                      * pEFVerts{nullptr};

    // The number of membrane triangles
    uint                                        pEFNTris{0};
    // Array of membrane triangles
    // \TODO use vector<tri_verts> pointer instead
    vertex_id_t                              * pEFTris{nullptr};

    std::vector<steps::tetode::Tri *>        pEFTris_vec;

    // The number of tetrahedrons
    uint                                        pEFNTets{0};
    // Array of tetrahedrons
    // \TODO use vector<tet_verts> pointer instead
    vertex_id_t                              * pEFTets{nullptr};

    // Table of global vertex index to EField local vertex index (0, 1, ..., pEFNVerts - 1)
    vertex_id_t                              * pEFVert_GtoL{0};

    // Table of global triangle index to EField local triangle index (0, 1, ..., pEFNTris-1)
    triangle_id_t                            * pEFTri_GtoL{0};

    // Table of global tetrahedron index to EField local tet index (0, 1, ..., pEFNTets-1)
    tetrahedron_id_t                         * pEFTet_GtoL{0};

    // Table of EField local triangle index to global triangle index.
    triangle_id_t                                     * pEFTri_LtoG{0};

};


////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_TETODE_TETODE_HPP

// END
