/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

// logging
#include <easylogging++.h>

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

    Tetexact(steps::model::Model *m, steps::wm::Geom *g, const rng::RNGptr &r,
             int calcMembPot = EF_NONE);
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
    //void advanceSteps(uint nsteps);
    void step() override;

    void checkpoint(std::string const & file_name) override;
    void restore(std::string const & file_name) override;
    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    void setEfieldDT(double efdt) override;

    inline double efdt() const noexcept
    { return pEFDT; }

    inline double getEfieldDT() const noexcept override
    { return pEFDT; }

    void setTemp(double t) override;

    inline double getTemp() const noexcept override
    { return pTemp; }

    // save the optimal vertex indexing
    void saveMembOpt(std::string const & opt_file_name);

    ////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      GENERAL
    ////////////////////////////////////////////////////////////////////////

    double getTime() const override;

    inline double getA0() const noexcept override
    { return pA0; }

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

     std::vector<double> getBatchTetCounts(const std::vector<index_t> &tets, std::string const &s) const override;

     std::vector<double> getBatchTriCounts(const std::vector<index_t> &tris, std::string const &s) const override;

     void getBatchTetCountsNP(const index_t *indices,
                              size_t input_size,
                              std::string const &s,
                              double *counts,
                              size_t output_size) const override;

     void getBatchTriCountsNP(const index_t *indices,
                              size_t input_size,
                              std::string const &s,
                              double *counts,
                              size_t output_size) const override;

     ////////////////////////////////////////////////////////////////////////
     // ROI Data Access
     ////////////////////////////////////////////////////////////////////////

     /// Get species counts of a list of tetrahedrons
     std::vector<double> getROITetCounts(const std::string& ROI_id, std::string const & s) const override;

     /// Get species counts of a list of triangles
     std::vector<double> getROITriCounts(const std::string& ROI_id, std::string const & s) const override;

     /// Get species counts of a list of tetrahedrons
     void getROITetCountsNP(const std::string& ROI_id, std::string const & s, double* counts, size_t output_size) const override;

     /// Get species counts of a list of triangles
     void getROITriCountsNP(const std::string& ROI_id, std::string const & s, double* counts, size_t output_size) const override;

     double getROIVol(const std::string& ROI_id) const override;
     double getROIArea(const std::string& ROI_id) const override;

     double getROICount(const std::string& ROI_id, std::string const & s) const override;
     void setROICount(const std::string& ROI_id, std::string const & s, double count) override;

     double getROIAmount(const std::string& ROI_id, std::string const & s) const override;
     void setROIAmount(const std::string& ROI_id, std::string const & s, double amount) override;

     double getROIConc(const std::string& ROI_id, std::string const & s) const override;
     void setROIConc(const std::string& ROI_id, std::string const & s, double conc) override;

     void setROIClamped(const std::string& ROI_id, std::string const & s, bool b) override;

     void setROIReacK(const std::string& ROI_id, std::string const & r, double kf) override;
     void setROISReacK(const std::string& ROI_id, std::string const & sr, double kf) override;
     void setROIDiffD(const std::string& ROI_id, std::string const & d, double dk) override;

     void setROIReacActive(const std::string& ROI_id, std::string const & r, bool a) override;
     void setROISReacActive(const std::string& ROI_id, std::string const & sr, bool a) override;
     void setROIDiffActive(const std::string& ROI_id, std::string const & d, bool a) override;
     void setROIVDepSReacActive(const std::string& ROI_id, std::string const & vsr, bool a) override;

     unsigned long long getROIReacExtent(const std::string& ROI_id, std::string const & r) const override;
     void resetROIReacExtent(const std::string& ROI_id, std::string const & r) override;

     unsigned long long getROISReacExtent(const std::string& ROI_id, std::string const & sr) const override;
     void resetROISReacExtent(const std::string& ROI_id, std::string const & sr) override;

     unsigned long long getROIDiffExtent(const std::string& ROI_id, std::string const & d) const override;
     void resetROIDiffExtent(const std::string& ROI_id, std::string const & d) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

     double _getCompVol(uint cidx) const override;

     double _getCompCount(uint cidx, uint sidx) const override;
     void _setCompCount(uint cidx, uint sidx, double n) override;

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

    double _getCompReacH(uint cidx, uint ridx) const override;
    double _getCompReacC(uint cidx, uint ridx) const override;
    long double _getCompReacA(uint cidx, uint ridx) const override;

    unsigned long long _getCompReacExtent(uint cidx, uint ridx) const override;
    void _resetCompReacExtent(uint cidx, uint ridx) override;

    double _getCompDiffD(uint cidx, uint didx) const override;
    void _setCompDiffD(uint cidx, uint didx, double dk) override;

    bool _getCompDiffActive(uint cidx, uint didx) const override;
    void _setCompDiffActive(uint cidx, uint didx, bool act) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      PATCH
    ////////////////////////////////////////////////////////////////////////

    double _getPatchArea(uint pidx) const override;

     double _getPatchCount(uint pidx, uint sidx) const override;
    void _setPatchCount(uint pidx, uint sidx, double n) override;

    double _getPatchAmount(uint pidx, uint sidx) const override;
     void _setPatchAmount(uint pidx, uint sidx, double a) override;

    bool _getPatchClamped(uint pidx, uint sidx) const override;
    void _setPatchClamped(uint pidx, uint sidx, bool buf) override;

    double _getPatchSReacK(uint pidx, uint ridx) const override;
      void _setPatchSReacK(uint pidx, uint ridx, double kf) override;

     bool _getPatchSReacActive(uint pidx, uint ridx) const override;
     void _setPatchSReacActive(uint pidx, uint ridx, bool a) override;

    double _getPatchSReacH(uint pidx, uint ridx) const override;
    double _getPatchSReacC(uint pidx, uint ridx) const override;
    double _getPatchSReacA(uint pidx, uint ridx) const override;

    unsigned long long _getPatchSReacExtent(uint pidx, uint ridx) const override;
    void _resetPatchSReacExtent(uint pidx, uint ridx) override;

    bool _getPatchVDepSReacActive(uint pidx, uint vsridx) const override;
    void _setPatchVDepSReacActive(uint pidx, uint vsridx, bool a) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      DIFFUSION BOUNDARIES
    ////////////////////////////////////////////////////////////////////////

    void _setDiffBoundaryDiffusionActive(uint dbidx, uint didx, bool act) override;
    bool _getDiffBoundaryDiffusionActive(uint dbidx, uint didx) const override;
    void _setDiffBoundaryDcst(uint dbidx, uint sidx, double dcst, uint direction_comp = std::numeric_limits<uint>::max()) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      SURFACE DIFFUSION BOUNDARIES
    ////////////////////////////////////////////////////////////////////////

    void _setSDiffBoundaryDiffusionActive(uint sdbidx, uint didx, bool act) override;
    bool _getSDiffBoundaryDiffusionActive(uint sdbidx, uint didx) const override;
    void _setSDiffBoundaryDcst(uint sdbidx, uint sidx, double dcst, uint direction_patch = std::numeric_limits<uint>::max()) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TETRAHEDRAL VOLUME ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    double _getTetVol(tetrahedron_id_t tidx) const override;
    void _setTetVol(tetrahedron_id_t tidx, double vol) override;

    bool _getTetSpecDefined(tetrahedron_id_t tidx, uint sidx) const override;

    double _getTetCount(tetrahedron_id_t tidx, uint sidx) const override;
    void _setTetCount(tetrahedron_id_t tidx, uint sidx, double n) override;

    double _getTetAmount(tetrahedron_id_t tidx, uint sidx) const override;
    void _setTetAmount(tetrahedron_id_t tidx, uint sidx, double m) override;

    double _getTetConc(tetrahedron_id_t tidx, uint sidx) const override;
    void _setTetConc(tetrahedron_id_t tidx, uint sidx, double c) override;

    bool _getTetClamped(tetrahedron_id_t tidx, uint sidx) const override;
    void _setTetClamped(tetrahedron_id_t tidx, uint sidx, bool buf) override;

    double _getTetReacK(tetrahedron_id_t tidx, uint ridx) const override;
    void _setTetReacK(tetrahedron_id_t tidx, uint ridx, double kf) override;

    bool _getTetReacActive(tetrahedron_id_t tidx, uint ridx) const override;
    void _setTetReacActive(tetrahedron_id_t tidx, uint ridx, bool act) override;

    double _getTetDiffD(tetrahedron_id_t tidx, uint didx, tetrahedron_id_t direction_tet = UNKNOWN_TET) const override;
    void _setTetDiffD(tetrahedron_id_t tidx, uint didx, double dk,
                      tetrahedron_id_t direction_tet = UNKNOWN_TET) override;

    bool _getTetDiffActive(tetrahedron_id_t tidx, uint didx) const override;
    void _setTetDiffActive(tetrahedron_id_t tidx, uint didx, bool act) override;

    ////////////////////////////////////////////////////////////////////////

    double _getTetReacH(tetrahedron_id_t tidx, uint ridx) const override;
    double _getTetReacC(tetrahedron_id_t tidx, uint ridx) const override;
    double _getTetReacA(tetrahedron_id_t tidx, uint ridx) const override;

    double _getTetDiffA(tetrahedron_id_t tidx, uint didx) const override;

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    double _getTetV(tetrahedron_id_t tidx) const override;
    void _setTetV(tetrahedron_id_t tidx, double v) override;
    bool _getTetVClamped(tetrahedron_id_t tidx) const override;
    void _setTetVClamped(tetrahedron_id_t tidx, bool cl) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TRIANGULAR SURFACE ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    double _getTriArea(triangle_id_t tidx) const override;
    void _setTriArea(triangle_id_t tidx, double area) override;

    bool _getTriSpecDefined(triangle_id_t tidx, uint sidx) const override;

    double _getTriCount(triangle_id_t tidx, uint sidx) const override;
    void _setTriCount(triangle_id_t tidx, uint sidx, double n) override;

    double _getTriAmount(triangle_id_t tidx, uint sidx) const override;
    void _setTriAmount(triangle_id_t tidx, uint sidx, double m) override;

    bool _getTriClamped(triangle_id_t tidx, uint sidx) const override;
    void _setTriClamped(triangle_id_t tidx, uint sidx, bool buf) override;

    double _getTriSReacK(triangle_id_t tidx, uint ridx) const override;
    void _setTriSReacK(triangle_id_t tidx, uint ridx, double kf) override;

    bool _getTriSReacActive(triangle_id_t tidx, uint ridx) const override;
    void _setTriSReacActive(triangle_id_t tidx, uint ridx, bool act) override;
    
    double _getTriSDiffD(triangle_id_t tidx, uint didx, triangle_id_t direction_tri) const override;
    
    void _setTriSDiffD(triangle_id_t tidx, uint didx, double dk,
                       triangle_id_t direction_tri) override;

    ////////////////////////////////////////////////////////////////////////

    double _getTriSReacH(triangle_id_t tidx, uint ridx) const override;
    double _getTriSReacC(triangle_id_t tidx, uint ridx) const override;
    double _getTriSReacA(triangle_id_t tidx, uint ridx) const override;

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    double _getTriV(triangle_id_t tidx) const override;
    void _setTriV(triangle_id_t tidx, double v) override;
    bool _getTriVClamped(triangle_id_t tidx) const override;
    void _setTriVClamped(triangle_id_t tidx, bool cl) override;

    double _getTriOhmicI(triangle_id_t tidx) const override;
    double _getTriOhmicI(triangle_id_t tidx, uint ocidx) const override;

    double _getTriGHKI(triangle_id_t tidx) const override;
    double _getTriGHKI(triangle_id_t tidx, uint ghkidx) const override;

    double _getTriI(triangle_id_t tidx) const override;

    double _getTriIClamp(triangle_id_t tidx) const override;
    void _setTriIClamp(triangle_id_t tidx, double cur) override;

    bool _getTriVDepSReacActive(triangle_id_t tidx, uint vsridx) const override;
    void _setTriVDepSReacActive(triangle_id_t tidx, uint vsridx, bool act) override;

    void _setTriCapac(triangle_id_t tidx, double cap) override;

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

    void _setMembPotential(uint midx, double v) override;
    void _setMembCapac(uint midx, double cm) override;
    void _setMembVolRes(uint midx, double ro) override;
    void _setMembRes(uint midx, double ro, double vrev) override;

    ////////////////////////////////////////////////////////////////////////

    // Called from local Comp or Patch objects. Add KProc to this object
    void addKProc(steps::tetexact::KProc * kp);

    inline uint countKProcs() const
    { return pKProcs.size(); }

    ////////////////////////////////////////////////////////////////////////

    inline const steps::tetmesh::Tetmesh& mesh() const noexcept
    { return *pMesh; }

    inline steps::tetmesh::Tetmesh& mesh() noexcept
    { return *pMesh; }

    inline steps::tetexact::Comp * _comp(uint cidx) const {
        // Moved assertions to this access routine to cut code duplication.
        AssertLog(cidx < statedef().countComps());
        AssertLog(statedef().countComps() == pComps.size());

        auto c = pComps[cidx];
        //AssertLog(c !=0);
        return c;
    }

    inline const std::vector<steps::tetexact::Patch *>&  patches() const
    { return pPatches; }

    inline steps::tetexact::Patch * _patch(uint pidx) const {
        // Moved assertions to this access routine to cut code duplication.
        AssertLog(pidx < statedef().countPatches());
        AssertLog(statedef().countPatches() == pPatches.size());

        auto p = pPatches[pidx];
        AssertLog(p !=0);
        return p;
    }

    inline steps::tetexact::DiffBoundary * _diffboundary(uint dbidx) const {
        AssertLog(dbidx < statedef().countDiffBoundaries());
        return pDiffBoundaries[dbidx];
    }

    inline steps::tetexact::SDiffBoundary * _sdiffboundary(uint sdbidx) const {
        AssertLog(sdbidx < statedef().countSDiffBoundaries());
        return pSDiffBoundaries[sdbidx];
    }

    inline steps::tetexact::Tet * _tet(tetrahedron_id_t tidx) const
    { return pTets[tidx.get()]; }

    inline steps::tetexact::Tri * _tri(triangle_id_t tidx) const
    { return pTris[tidx.get()]; }

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

    std::size_t _addComp(steps::solver::Compdef * cdef);

    uint _addPatch(steps::solver::Patchdef * pdef);

    uint _addDiffBoundary(steps::solver::DiffBoundarydef * dbdef);

    uint _addSDiffBoundary(steps::solver::SDiffBoundarydef * sdbdef);


    void _addTet(tetrahedron_id_t tetidx,
                 steps::tetexact::Comp *comp,
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

    void _addWmVol(uint cidx, steps::tetexact::Comp * comp, double vol);

    void _addTri(triangle_id_t triidx,
                 steps::tetexact::Patch *patch,
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
    void _updateSpec(steps::tetexact::WmVol * tet);

    /// Update the kproc's of a triangle, after a species has been changed.
    /// This does not need to update the kproc's of any neighbouring
    /// tetrahedrons.
    ///
    /// Currently doesn't care about the species.
    ///
    void _updateSpec(steps::tetexact::Tri * tri);


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

    double getROITetCount(const std::vector<tetrahedron_id_t>& indices, const std::string& s) const;
    double getROITriCount(const std::vector<triangle_id_t>& indices, const std::string& s) const;

    void setROITetClamped(const std::vector<tetrahedron_id_t>& indices, std::string const & s, bool b);
    void setROITriClamped(const std::vector<triangle_id_t>& indices, std::string const & s, bool b);

    void setROITriCount(const std::vector<triangle_id_t>& indices, std::string const & s, double count);
    void setROITetCount(const std::vector<tetrahedron_id_t >& indices, std::string const & s, double count);

    double getROIVol(const std::vector<tetrahedron_id_t>& tets)const;

    steps::tetmesh::Tetmesh *                    pMesh{nullptr};

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
    std::size_t                                 nEntries;
    double                                      pSum;
    double                                      nSum;
    double                                      pA0{0.0};

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
        group->indices = static_cast<KProc**>(realloc(group->indices,
                                          sizeof(KProc*) * group->capacity));
        if (group->indices == nullptr) {
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

    // The Efield solve choise. If EF_NONE we don't calclulate the potential, nor include
    // any voltage-dependent transitions, ohmic or ghk currents. This means
    // the solver behaves exactly like the previous Tetexact solver in this case
    // and all following members are set to null pointer or zero.
    //
    EF_solver                                    pEFoption;

    double                                       pTemp{0.0};

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
    vertex_id_t                               * pEFTris{nullptr};

    std::vector<steps::tetexact::Tri *>        pEFTris_vec;

    // The number of tetrahedrons
    uint                                        pEFNTets{0};
    // Array of tetrahedrons
    // \TODO TCL: use vector<tet_verts> instead
    vertex_id_t                               * pEFTets{nullptr};

    // Table of global vertex index to EField local vertex index (0, 1, ..., pEFNVerts - 1)
    vertex_id_t                               * pEFVert_GtoL{nullptr};

    // Table of global triangle index to EField local triangle index (0, 1, ..., pEFNTris-1)
    triangle_id_t                             * pEFTri_GtoL{nullptr};

    // Table of global tetrahedron index to EField local tet index (0, 1, ..., pEFNTets-1)
    tetrahedron_id_t                          * pEFTet_GtoL{nullptr};

    // Table of EField local triangle index to global triangle index.
    triangle_id_t                             * pEFTri_LtoG{nullptr};


};

////////////////////////////////////////////////////////////////////////////////

}} // namespace steps::tetexact

#endif
// STEPS_TETEXACT_TETEXACT_HPP

