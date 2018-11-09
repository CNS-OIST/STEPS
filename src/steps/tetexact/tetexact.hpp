/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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


#ifndef STEPS_TETEXACT_TETEXACT_HPP
#define STEPS_TETEXACT_TETEXACT_HPP 1

// STL headers.
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/api.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/geom/tetmesh.hpp"
#include "steps/tetexact/tri.hpp"
#include "steps/tetexact/tet.hpp"
#include "steps/tetexact/wmvol.hpp"
#include "steps/tetexact/kproc.hpp"
#include "steps/tetexact/comp.hpp"
#include "steps/tetexact/patch.hpp"
#include "steps/tetexact/diffboundary.hpp"
#include "steps/tetexact/sdiffboundary.hpp"
#include "steps/tetexact/crstruct.hpp"
#include "steps/solver/efield/efield.hpp"

// logging
#include "third_party/easyloggingpp/src/easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace steps {
namespace tetexact {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.


// Auxiliary declarations.
typedef uint                            SchedIDX;
typedef std::set<SchedIDX>              SchedIDXSet;
typedef SchedIDXSet::iterator           SchedIDXSetI;
typedef SchedIDXSet::const_iterator     SchedIDXSetCI;
typedef std::vector<SchedIDX>           SchedIDXVec;
typedef SchedIDXVec::iterator           SchedIDXVecI;
typedef SchedIDXVec::const_iterator     SchedIDXVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Copies the contents of a set of SchedIDX entries into a vector.
/// The contents of the vector are completely overridden.
///
extern void schedIDXSet_To_Vec(SchedIDXSet const & s, SchedIDXVec & v);

////////////////////////////////////////////////////////////////////////////////

class Tetexact: public steps::solver::API
{

public:

    Tetexact(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r,
             int calcMembPot = EF_NONE);
    ~Tetexact();


    ////////////////////////////////////////////////////////////////////////
    // SOLVER INFORMATION
    ////////////////////////////////////////////////////////////////////////

    std::string getSolverName() const;
    std::string getSolverDesc() const;
    std::string getSolverAuthors() const;
    std::string getSolverEmail() const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS
    ////////////////////////////////////////////////////////////////////////

    void reset();
    void run(double endtime);
    void advance(double adv);
    //void advanceSteps(uint nsteps);
    void step();

    void checkpoint(std::string const & file_name);
    void restore(std::string const & file_name);
    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    void setEfieldDT(double efdt);

    inline double efdt() const
    { return pEFDT; }

    inline double getEfieldDT() const
    { return pEFDT; }

    void setTemp(double t);

    inline double getTemp() const
    { return pTemp; }

    // save the optimal vertex indexing
    void saveMembOpt(std::string const & opt_file_name);

    ////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      GENERAL
    ////////////////////////////////////////////////////////////////////////

    double getTime() const;

    inline double getA0() const
    { return pA0; }

    uint getNSteps() const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      ADVANCE
    //      Developer only
    ////////////////////////////////////////////////////////////////////////

    void setTime(double time);
    void setNSteps(uint nsteps);

    ////////////////////////////////////////////////////////////////////////
     // Batch Data Access
     ////////////////////////////////////////////////////////////////////////

     std::vector<double> getBatchTetCounts(std::vector<uint> const & tets, std::string const & s) const;

     std::vector<double> getBatchTriCounts(std::vector<uint> const & tris, std::string const & s) const;

     void getBatchTetCountsNP(unsigned int* indices, int input_size, std::string const & s, double* counts, int output_size) const;

     void getBatchTriCountsNP(unsigned int* indices, int input_size, std::string const & s, double* counts, int output_size) const;

     ////////////////////////////////////////////////////////////////////////
     // ROI Data Access
     ////////////////////////////////////////////////////////////////////////

     /// Get species counts of a list of tetrahedrons
     std::vector<double> getROITetCounts(std::string ROI_id, std::string const & s) const;

     /// Get species counts of a list of triangles
     std::vector<double> getROITriCounts(std::string ROI_id, std::string const & s) const;

     /// Get species counts of a list of tetrahedrons
     void getROITetCountsNP(std::string ROI_id, std::string const & s, double* counts, int output_size) const;

     /// Get species counts of a list of triangles
     void getROITriCountsNP(std::string ROI_id, std::string const & s, double* counts, int output_size) const;

     double getROIVol(std::string ROI_id) const;
     double getROIArea(std::string ROI_id) const;

     double getROICount(std::string ROI_id, std::string const & s) const;
     void setROICount(std::string ROI_id, std::string const & s, double count);

     double getROIAmount(std::string ROI_id, std::string const & s) const;
     double getROIConc(std::string ROI_id, std::string const & s) const;
     void setROIConc(std::string ROI_id, std::string const & s, double conc);

     void setROIClamped(std::string ROI_id, std::string const & s, bool b);

     void setROIReacK(std::string ROI_id, std::string const & r, double kf);
     void setROISReacK(std::string ROI_id, std::string const & sr, double kf);
     void setROIDiffD(std::string ROI_id, std::string const & d, double dk);

     void setROIReacActive(std::string ROI_id, std::string const & r, bool a);
     void setROISReacActive(std::string ROI_id, std::string const & sr, bool a);
     void setROIDiffActive(std::string ROI_id, std::string const & d, bool a);
     void setROIVDepSReacActive(std::string ROI_id, std::string const & vsr, bool a);

     uint getROIReacExtent(std::string ROI_id, std::string const & r) const;
     void resetROIReacExtent(std::string ROI_id, std::string const & r);

     uint getROISReacExtent(std::string ROI_id, std::string const & sr) const;
     void resetROISReacExtent(std::string ROI_id, std::string const & sr);

     uint getROIDiffExtent(std::string ROI_id, std::string const & d) const;
     void resetROIDiffExtent(std::string ROI_id, std::string const & d);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

     double _getCompVol(uint cidx) const;

     double _getCompCount(uint cidx, uint sidx) const;
     void _setCompCount(uint cidx, uint sidx, double n);

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

    double _getCompReacH(uint cidx, uint ridx) const;
    double _getCompReacC(uint cidx, uint ridx) const;
    double _getCompReacA(uint cidx, uint ridx) const;

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

    double _getPatchArea(uint pidx) const;

     double _getPatchCount(uint pidx, uint sidx) const;
    void _setPatchCount(uint pidx, uint sidx, double n);

    double _getPatchAmount(uint pidx, uint sidx) const;
     void _setPatchAmount(uint pidx, uint sidx, double a);

    bool _getPatchClamped(uint pidx, uint sidx) const;
    void _setPatchClamped(uint pidx, uint sidx, bool buf);

    double _getPatchSReacK(uint pidx, uint ridx) const;
      void _setPatchSReacK(uint pidx, uint ridx, double kf);

     bool _getPatchSReacActive(uint pidx, uint ridx) const;
     void _setPatchSReacActive(uint pidx, uint ridx, bool a);

    double _getPatchSReacH(uint pidx, uint ridx) const;
    double _getPatchSReacC(uint pidx, uint ridx) const;
    double _getPatchSReacA(uint pidx, uint ridx) const;

    uint _getPatchSReacExtent(uint pidx, uint ridx) const;
    void _resetPatchSReacExtent(uint pidx, uint ridx);

    bool _getPatchVDepSReacActive(uint pidx, uint vsridx) const;
    void _setPatchVDepSReacActive(uint pidx, uint vsridx, bool a);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      DIFFUSION BOUNDARIES
    ////////////////////////////////////////////////////////////////////////

    void _setDiffBoundaryDiffusionActive(uint dbidx, uint didx, bool act);
    bool _getDiffBoundaryDiffusionActive(uint dbidx, uint didx) const;
    void _setDiffBoundaryDcst(uint dbidx, uint sidx, double dcst, uint direction_comp = std::numeric_limits<uint>::max());

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      SURFACE DIFFUSION BOUNDARIES
    ////////////////////////////////////////////////////////////////////////

    void _setSDiffBoundaryDiffusionActive(uint sdbidx, uint didx, bool act);
    bool _getSDiffBoundaryDiffusionActive(uint sdbidx, uint didx) const;
    void _setSDiffBoundaryDcst(uint sdbidx, uint sidx, double dcst, uint direction_patch = std::numeric_limits<uint>::max());

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TETRAHEDRAL VOLUME ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    double _getTetVol(uint tidx) const;
    void _setTetVol(uint tidx, double vol);

    bool _getTetSpecDefined(uint tidx, uint sidx) const;

    double _getTetCount(uint tidx, uint sidx) const;
    void _setTetCount(uint tidx, uint sidx, double n);

    double _getTetAmount(uint tidx, uint sidx) const;
    void _setTetAmount(uint tidx, uint sidx, double m);

    double _getTetConc(uint tidx, uint sidx) const;
    void _setTetConc(uint tidx, uint sidx, double c);

    bool _getTetClamped(uint tidx, uint sidx) const;
    void _setTetClamped(uint tidx, uint sidx, bool buf);

    double _getTetReacK(uint tidx, uint ridx) const;
    void _setTetReacK(uint tidx, uint ridx, double kf);

    bool _getTetReacActive(uint tidx, uint ridx) const;
    void _setTetReacActive(uint tidx, uint ridx, bool act);

    double _getTetDiffD(uint tidx, uint didx, uint direction_tet = std::numeric_limits<uint>::max()) const;
    void _setTetDiffD(uint tidx, uint didx, double dk,
                      uint direction_tet = std::numeric_limits<uint>::max());

    bool _getTetDiffActive(uint tidx, uint didx) const;
    void _setTetDiffActive(uint tidx, uint didx, bool act);

    ////////////////////////////////////////////////////////////////////////

    double _getTetReacH(uint tidx, uint ridx) const;
    double _getTetReacC(uint tidx, uint ridx) const;
    double _getTetReacA(uint tidx, uint ridx) const;

    double _getTetDiffA(uint tidx, uint didx) const;

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    double _getTetV(uint tidx) const;
    void _setTetV(uint tidx, double v);
    bool _getTetVClamped(uint tidx) const;
    void _setTetVClamped(uint tidx, bool cl);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TRIANGULAR SURFACE ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    double _getTriArea(uint tidx) const;
    void _setTriArea(uint tidx, double area);

    bool _getTriSpecDefined(uint tidx, uint sidx) const;

    double _getTriCount(uint tidx, uint sidx) const;
    void _setTriCount(uint tidx, uint sidx, double n);

    double _getTriAmount(uint tidx, uint sidx) const;
    void _setTriAmount(uint tidx, uint sidx, double m);

    bool _getTriClamped(uint tidx, uint sidx) const;
    void _setTriClamped(uint tidx, uint sidx, bool buf);

    double _getTriSReacK(uint tidx, uint ridx) const;
    void _setTriSReacK(uint tidx, uint ridx, double kf);

    bool _getTriSReacActive(uint tidx, uint ridx) const;
    void _setTriSReacActive(uint tidx, uint ridx, bool act);
    
    double _getTriSDiffD(uint tidx, uint didx, uint direction_tri = std::numeric_limits<uint>::max()) const;
    
    void _setTriSDiffD(uint tidx, uint didx, double dk,
                      uint direction_tri = std::numeric_limits<uint>::max());

    ////////////////////////////////////////////////////////////////////////

    double _getTriSReacH(uint tidx, uint ridx) const;
    double _getTriSReacC(uint tidx, uint ridx) const;
    double _getTriSReacA(uint tidx, uint ridx) const;

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    double _getTriV(uint tidx) const;
    void _setTriV(uint tidx, double v);
    bool _getTriVClamped(uint tidx) const;
    void _setTriVClamped(uint tidx, bool cl);

    double _getTriOhmicI(uint tidx);
    double _getTriOhmicI(uint tidx, uint ocidx);

    double _getTriGHKI(uint tidx);
    double _getTriGHKI(uint tidx, uint ghkidx);

    double _getTriI(uint tidx) const;

    void _setTriIClamp(uint tidx, double cur);

    bool _getTriVDepSReacActive(uint tidx, uint vsridx) const;
    void _setTriVDepSReacActive(uint tidx, uint vsridx, bool act);

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

    // Called from local Comp or Patch objects. Add KProc to this object
    void addKProc(steps::tetexact::KProc * kp);

    inline uint countKProcs() const
    { return pKProcs.size(); }

    ////////////////////////////////////////////////////////////////////////

    inline steps::tetmesh::Tetmesh * mesh() const
    { return pMesh; }

    inline steps::tetexact::Comp * _comp(uint cidx) const {
        // Moved assertions to this access routine to cut code duplication.
        AssertLog(cidx < statedef()->countComps());
        AssertLog(statedef()->countComps() == pComps.size());

        auto c = pComps[cidx];
        //AssertLog(c !=0);
        return c;
    }

    inline std::vector<steps::tetexact::Patch *>  patches() const
    { return pPatches; }

    inline steps::tetexact::Patch * _patch(uint pidx) const {
        // Moved assertions to this access routine to cut code duplication.
        AssertLog(pidx < statedef()->countPatches());
        AssertLog(statedef()->countPatches() == pPatches.size());

        auto p = pPatches[pidx];
        AssertLog(p !=0);
        return p;
    }

    inline steps::tetexact::DiffBoundary * _diffboundary(uint dbidx) const {
        AssertLog(dbidx < statedef()->countDiffBoundaries());
        return pDiffBoundaries[dbidx];
    }

    inline steps::tetexact::SDiffBoundary * _sdiffboundary(uint sdbidx) const {
        AssertLog(sdbidx < statedef()->countSDiffBoundaries());
        return pSDiffBoundaries[sdbidx];
    }

    inline steps::tetexact::Tet * _tet(uint tidx) const
    { return pTets[tidx]; }

    inline steps::tetexact::Tri * _tri(uint tidx) const
    { return pTris[tidx]; }

    inline double a0() const
    { return pA0; }

    // Checked global to local index translations
 
    uint specG2L_or_throw(Comp *comp, uint gidx) const;
    uint specG2L_or_throw(Patch *patch, uint gidx) const;
#if 0
    uint specG2L_or_throw(Tet *tet, uint gidx) const;
    uint specG2L_or_throw(Tri *tri, uint gidx) const;
#endif
    uint reacG2L_or_throw(Comp *comp, uint gidx) const;
    uint sreacG2L_or_throw(Patch *patch, uint gidx) const;
    uint diffG2L_or_throw(Comp *comp, uint gidx) const;
    uint sdiffG2L_or_throw(Patch *patch, uint gidx) const;
    uint vdepsreacG2L_or_throw(Patch *patch, uint gidx) const;

    ////////////////////////////////////////////////////////////////////////
    // TETEXACT SOLVER METHODS
    ////////////////////////////////////////////////////////////////////////

    uint _addComp(steps::solver::Compdef * cdef);

    uint _addPatch(steps::solver::Patchdef * pdef);

    uint _addDiffBoundary(steps::solver::DiffBoundarydef * dbdef);

    uint _addSDiffBoundary(steps::solver::SDiffBoundarydef * sdbdef);


    void _addTet(uint tetidx, steps::tetexact::Comp * comp, double vol, double a1,
                 double a2, double a3, double a4, double d1, double d2,
                 double d3, double d4, int tet0, int tet1, int tet2, int tet3);

    void _addWmVol(uint cidx, steps::tetexact::Comp * comp, double vol);

    void _addTri(uint triidx, steps::tetexact::Patch * patch, double area,
            double l0, double l1, double l2, double d0, double d1, double d2,
            int tinner, int touter, int tri0, int tri1, int tri2);

    // called when local tet, tri, reac, sreac objects have been created
    // by constructor
    void _setup();


    //void _build();

    double _getRate(uint i) const
    { return pKProcs[i]->rate(); }

    steps::tetexact::KProc * _getNext() const;

    //void _reset();

    void _executeStep(steps::tetexact::KProc * kp, double dt);

    // TODO: Change the following so that only the kprocs depending on
    // the species are updated. These functions are called from interface
    // methods setting compartment or patch counts.
    /// Update the kproc's of a tet, after a species has been changed.
    /// This also updates kproc's in surrounding triangles.
    ///
    /// Currently doesn't care about the species.
    ///
    void _updateSpec(steps::tetexact::WmVol * tet, uint spec_lidx);

    /// Update the kproc's of a triangle, after a species has been changed.
    /// This does not need to update the kproc's of any neighbouring
    /// tetrahedrons.
    ///
    /// Currently doesn't care about the species.
    ///
    void _updateSpec(steps::tetexact::Tri * tri, uint spec_lidx);


    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    /// Check the EField flag
    inline bool efflag() const
    { return pEFoption != EF_NONE; }

    void _setupEField();

    inline uint neftets() const
    { return pEFNTets; }

    inline uint neftris() const
    { return pEFNTris; }

    inline uint nefverts() const
    { return pEFNVerts; }

        
    ////////////////////////////////////////////////////////////////////////
private:

    ////////////////////////////////////////////////////////////////////////

    steps::tetmesh::Tetmesh *                    pMesh;

    ////////////////////////////////////////////////////////////////////////
    // LIST OF TETEXACT SOLVER OBJECTS
    ////////////////////////////////////////////////////////////////////////

    std::vector<steps::tetexact::Comp *>       pComps;
    std::map<steps::solver::Compdef *, Comp *> pCompMap;

    std::vector<steps::tetexact::Patch *>      pPatches;

    std::vector<steps::tetexact::DiffBoundary *> pDiffBoundaries;

    std::vector<steps::tetexact::SDiffBoundary *> pSDiffBoundaries;

    // These objects are used to describe a mesh compartment that is
    // being treated as a well-mixed volume.
    std::vector<steps::tetexact::WmVol *>      pWmVols;

    std::vector<steps::tetexact::Tri *>        pTris;

    // Now stored as base pointer
    std::vector<steps::tetexact::Tet *>        pTets;

    ////////////////////////////////////////////////////////////////////////
    // CR SSA Kernel Data and Methods
    ////////////////////////////////////////////////////////////////////////
    uint                                        nEntries;
    double                                      pSum;
    double                                      nSum;
    double                                      pA0;

    std::vector<KProc*>                         pKProcs;

    std::vector<CRGroup*>                       nGroups;
    std::vector<CRGroup*>                       pGroups;

    ////////////////////////////////////////////////////////////////////////////////

    template <typename KProcPIter>
    inline void _update(KProcPIter b, KProcPIter e) {
        while (b!=e) _updateElement(*b++);
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
        }
        else {
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
            curr_size ++;
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
            curr_size ++;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////

    inline void _extendGroup(CRGroup* group, uint size = 1024) {
        #ifdef SSA_DEBUG
        CLOG(INFO, "general_log") << "SSA: extending group storage\n";
        CLOG(INFO, "general_log") << "current capacity: " << group->capacity << "\n";
        #endif

        group->capacity += size;
        group->indices = (KProc**)realloc(group->indices,
                                          sizeof(KProc*) * group->capacity);
        if (group->indices == NULL) {
            SysErrLog("DirectCR: unable to allocate memory for SSA group.");
        }
        #ifdef SSA_DEBUG
        CLOG(INFO, "general_log") << "capacity after extending: " << group->capacity << "\n";
        CLOG(INFO, "general_log") << "--------------------------------------------------------\n";
        #endif
    }

    ////////////////////////////////////////////////////////////////////////////////

    void _updateElement(KProc* kp);

    inline void _updateSum() {
        #ifdef SSA_DEBUG
        CLOG(INFO, "general_log") << "update A0 from " << pA0 << " to ";
        #endif

        pA0 = 0.0;

        uint n_neg_groups = nGroups.size();
        uint n_pos_groups = pGroups.size();

        for (uint i = 0; i < n_neg_groups; i++) {
            pA0 += nGroups[i]->sum;
        }

        for (uint i = 0; i < n_pos_groups; i++) {
            pA0 += pGroups[i]->sum;
        }

        #ifdef SSA_DEBUG
        CLOG(INFO, "general_log") << pA0 << "\n";
        #endif
    }


    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    // The Efield solve choise. If EF_NONE we don't calclulate the potential, nor include
    // any voltage-dependent transitions, ohmic or ghk currents. This means
    // the solver behaves exactly like the previous Tetexact solver in this case
    // and all following members are set to null pointer or zero.
    //
    EF_solver                                    pEFoption;

    double                                       pTemp;

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

    std::vector<steps::tetexact::Tri *>        pEFTris_vec;

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

}} // namespace steps::tetexact

#endif
// STEPS_TETEXACT_TETEXACT_HPP

