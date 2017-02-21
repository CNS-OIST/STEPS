/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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
#include "steps/common.h"
#include "steps/solver/api.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/geom/tetmesh.hpp"
#include "steps/tetode/comp.hpp"
#include "steps/tetode/patch.hpp"
#include "steps/tetode/tet.hpp"
#include "steps/tetode/tri.hpp"

#include "steps/solver/efield/efield.hpp"


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

    TetODE(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r,
            int calcMembPot = EF_NONE);
    ~TetODE(void);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER INFORMATION
    ////////////////////////////////////////////////////////////////////////

    std::string getSolverName(void) const;
    std::string getSolverDesc(void) const;
    std::string getSolverAuthors(void) const;
    std::string getSolverEmail(void) const;


    void checkpoint(std::string const & file_name);

    void restore(std::string const & file_name);

    double getTime(void) const;

    inline double getTemp(void) const
    { return pTemp; }

    void setTemp(double t);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS
    ////////////////////////////////////////////////////////////////////////

    // Will of course be added- just for compilation purposes
    void reset(void);
    void run(double endtime);
    void advance(double adv);
    //void step(void);

    ////////////////////////////////////////////////////////////////////////

    inline steps::tetmesh::Tetmesh * mesh(void) const
    { return pMesh; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

     double _getCompCount(uint cidx, uint sidx) const;
     void _setCompCount(uint cidx, uint sidx, double n);

    double _getCompVol(uint cidx) const;

    double _getCompAmount(uint cidx, uint sidx) const;
    void _setCompAmount(uint cidx, uint sidx, double a);

    double _getCompConc(uint cidx, uint sidx) const;
    void _setCompConc(uint cidx, uint sidx, double c);

    bool _getCompClamped(uint cidx, uint sidx) const;
    void _setCompClamped(uint cidx, uint sidx, bool b);

    double _getCompReacK(uint cidx, uint ridx) const;
    void _setCompReacK(uint cidx, uint ridx, double kf);

    bool _getCompReacActive(uint cidx, uint ridx) const;
    void _setCompReacActive(uint cidx, uint ridx, bool a);

     /*
    uint _getCompReacExtent(uint cidx, uint ridx) const;
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
     double _getPatchCount(uint pidx, uint sidx) const;
    void _setPatchCount(uint pidx, uint sidx, double n);

    double _getPatchArea(uint pidx) const;

    double _getPatchAmount(uint pidx, uint sidx) const;
    void _setPatchAmount(uint pidx, uint sidx, double a);

    bool _getPatchClamped(uint pidx, uint sidx) const;
    void _setPatchClamped(uint pidx, uint sidx, bool buf);

    double _getPatchSReacK(uint pidx, uint ridx) const;
    void _setPatchSReacK(uint pidx, uint ridx, double kf);

    bool _getPatchSReacActive(uint pidx, uint ridx) const;
    void _setPatchSReacActive(uint pidx, uint ridx, bool a);

    /*
    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TETRAHEDRAL VOLUME ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    void _setTetVol(uint tidx, double vol);

    bool _getTetSpecDefined(uint tidx, uint sidx) const;

    */
    double _getTetVol(uint tidx) const;

    double _getTetCount(uint tidx, uint sidx) const;
    void _setTetCount(uint tidx, uint sidx, double n);

    double _getTetAmount(uint tidx, uint sidx) const;
    void _setTetAmount(uint tidx, uint sidx, double m);

    double _getTetConc(uint tidx, uint sidx) const;
    void _setTetConc(uint tidx, uint sidx, double c);

    double _getTetReacK(uint tidx, uint ridx) const;
    void _setTetReacK(uint tidx, uint ridx, double kf);

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

    double _getTetV(uint tidx) const;
    void _setTetV(uint tidx, double v);
    bool _getTetVClamped(uint tidx) const;
    void _setTetVClamped(uint tidx, bool cl);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TRIANGULAR SURFACE ELEMENTS
    ////////////////////////////////////////////////////////////////////////
    /*
    void _setTriArea(uint tidx, double area);

    bool _getTriSpecDefined(uint tidx, uint sidx) const;

    */

    double _getTriArea(uint tidx) const;

    double _getTriCount(uint tidx, uint sidx) const;
    void _setTriCount(uint tidx, uint sidx, double n);

    double _getTriAmount(uint tidx, uint sidx) const;
    void _setTriAmount(uint tidx, uint sidx, double m);

    double _getTriSReacK(uint tidx, uint ridx) const;
    void _setTriSReacK(uint tidx, uint ridx, double kf);

    /*
    bool _getTriClamped(uint tidx, uint sidx) const;
    void _setTriClamped(uint tidx, uint sidx, bool buf);

    bool _getTriSReacActive(uint tidx, uint ridx) const;
    void _setTriSReacActive(uint tidx, uint ridx, bool act);

     */

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    double _getTriV(uint tidx) const;
    void _setTriV(uint tidx, double v);
    bool _getTriVClamped(uint tidx) const;
    void _setTriVClamped(uint tidx, bool cl);

    /*
    double _getTriOhmicI(uint tidx) const;
    double _getTriOhmicI(uint tidx, uint ocidx) const;

    double _getTriGHKI(uint tidx) const;
    double _getTriGHKI(uint tidx, uint ghkidx) const;
    */

    double _getTriI(uint tidx) const;

    void _setTriIClamp(uint tidx, double cur);

    ////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      VERTICES ELEMENTS
    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    double _getVertV(uint vidx) const;
    void _setVertV(uint vidx, double v);
    bool _getVertVClamped(uint vidx) const;
    void _setVertVClamped(uint vidx, bool cl);
    void _setVertIClamp(uint tidx, double cur);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      MEMBRANE AND VOLUME CONDUCTOR
    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    void _setMembPotential(uint midx, double v);
    void _setMembCapac(uint midx, double cm);
    void _setMembVolRes(uint midx, double ro);
    void _setMembRes(uint midx, double ro, double vrev);

    ////////////////////////////////////////////////////////////////////////


    void _setup(void);

    uint _addComp(steps::solver::Compdef * cdef);

    uint _addPatch(steps::solver::Patchdef * pdef);

    void _addTet(uint tetidx, steps::tetode::Comp * comp, double vol, double a1,
                 double a2, double a3, double a4, double d1, double d2,
                 double d3, double d4, int tet0, int tet1, int tet2, int tet3);

    void _addTri(uint triidx, steps::tetode::Patch * patch, double area,
            double l0, double l1, double l2, double d0, double d1, double d2,
            int tinner, int touter, int tri0, int tri1, int tri2);

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
    inline bool efflag() const
    { return pEFoption != EF_NONE; }

    void _setupEField(void);

    inline uint neftets(void) const
    { return pEFNTets; }

    inline uint neftris(void) const
    { return pEFNTris; }

    inline uint nefverts(void) const
    { return pEFNVerts; }


private:

    ////////////////////////////////////////////////////////////////////////

    steps::tetmesh::Tetmesh *                  pMesh;

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

    uint                                      pSpecs_tot;
    uint                                      pReacs_tot;

    //double                                    * pCcst;

    bool                                      pInitialised;
    bool                                     pTolsset;
    bool                                      pReinit;

    CVodeState                              * pCVodeState;

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    // The Efield flag. If false we don't calclulate the potential, nor include
    // any voltage-dependent transitions, ohmic or ghk currents. This means
    // the solver behaves exactly like the previous Tetexact solver in this case
    // and all following members are set to null pointer or zero.
    //
    EF_solver                                     pEFoption;

    double                                        pTemp;

    // Pointer to the EField object
    std::unique_ptr<steps::solver::efield::EField> pEField;

    // The Efield time-step
    double                                       pEFDT;

    // The number of vertices
    uint                                        pEFNVerts;
    // Array of vertices
    double                                      * pEFVerts;

    // The number of membrane triangles
    uint                                        pEFNTris;
    // Array of membrane triangles
    uint                                       * pEFTris;

    std::vector<steps::tetode::Tri *>        pEFTris_vec;

    // The number of tetrahedrons
    uint                                        pEFNTets;
    // Array of tetrahedrons
    uint                                      * pEFTets;

    // Table of global vertex index to EField local vertex index (0, 1, ..., pEFNVerts - 1)
    int                                      * pEFVert_GtoL;

    // Table of global triangle index to EField local triangle index (0, 1, ..., pEFNTris-1)
    int                                      * pEFTri_GtoL;

    // Table of global tetrahedron index to EField local tet index (0, 1, ..., pEFNTets-1)
    int                                       * pEFTet_GtoL;

    // Table of EField local triangle index to global triangle index.
    uint                                      * pEFTri_LtoG;

};


////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_TETODE_TETODE_HPP

// END
