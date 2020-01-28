/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
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


#ifndef STEPS_MPI_TETOPSPLIT_TETOPSPLITP_HPP
#define STEPS_MPI_TETOPSPLIT_TETOPSPLITP_HPP 1


// STL headers.
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <memory>
#include <random>

// logging
#include <easylogging++.h>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/api.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/geom/tetmesh.hpp"
#include "steps/mpi/tetopsplit/tri.hpp"
#include "steps/mpi/tetopsplit/tet.hpp"
#include "steps/mpi/tetopsplit/wmvol.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
#include "steps/mpi/tetopsplit/comp.hpp"
#include "steps/mpi/tetopsplit/patch.hpp"
#include "steps/mpi/tetopsplit/diffboundary.hpp"
#include "steps/mpi/tetopsplit/sdiffboundary.hpp"
#include "steps/mpi/tetopsplit/crstruct.hpp"
#include "steps/solver/efield/efield.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace mpi{
namespace tetopsplit{

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

enum SubVolType {SUB_WM, SUB_TET, SUB_TRI};

////////////////////////////////////////////////////////////////////////////////

class TetOpSplitP: public steps::solver::API
{
public:

    TetOpSplitP(steps::model::Model *m, steps::wm::Geom *g, const rng::RNGptr &r,
                int calcMembPot = EF_NONE, std::vector<uint> const &tet_hosts = {},
                const std::map<triangle_id_t, uint> &tri_hosts = {},
                std::vector<uint> const &wm_hosts = {});
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

    double _getTriOhmicI(triangle_id_t tidx) override;
    double _getTriOhmicI(triangle_id_t tidx, uint ocidx) override;

    double _getTriGHKI(triangle_id_t tidx) override;
    double _getTriGHKI(triangle_id_t tidx, uint ghkidx) override;

    double _getTriI(triangle_id_t tidx) const override;

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

    uint addKProc(steps::mpi::tetopsplit::KProc * kp);

    void addDiff(Diff* diff);
    void addSDiff(SDiff* sdiff);
    inline uint countKProcs() const noexcept
    { return pKProcs.size(); }

    ////////////////////////////////////////////////////////////////////////

    inline steps::tetmesh::Tetmesh * mesh() const noexcept
    { return pMesh; }

    inline steps::mpi::tetopsplit::Comp * _comp(uint cidx) const
    {
        AssertLog(cidx < statedef().countComps());
        AssertLog(statedef().countComps() == pComps.size());
        return pComps[cidx];
    }

    inline const std::vector<steps::mpi::tetopsplit::Patch *>&  patches() const noexcept
    { return pPatches; }

    inline steps::mpi::tetopsplit::Patch * _patch(uint pidx) const
    {
        AssertLog(pidx < statedef().countPatches());
        AssertLog(statedef().countPatches() == pPatches.size());
        return pPatches[pidx];
    }

    inline steps::mpi::tetopsplit::DiffBoundary * _diffboundary(uint dbidx) const
    {
        AssertLog(dbidx < statedef().countDiffBoundaries());
        return pDiffBoundaries[dbidx];
    }

    inline steps::mpi::tetopsplit::SDiffBoundary * _sdiffboundary(uint sdbidx) const
    {
        AssertLog(sdbidx < statedef().countSDiffBoundaries());
        return pSDiffBoundaries[sdbidx];
    }

    inline steps::mpi::tetopsplit::Tet * _tet(tetrahedron_id_t tidx) const noexcept
    { return pTets[tidx.get()]; }

    inline steps::mpi::tetopsplit::Tri * _tri(triangle_id_t tidx) const noexcept
    { return pTris[tidx.get()]; }

    inline double a0() const noexcept
    { return pA0; }

    //inline bool built()
    //{ return pBuilt; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS
    ////////////////////////////////////////////////////////////////////////

    uint _addComp(steps::solver::Compdef * cdef);

    uint _addPatch(steps::solver::Patchdef * pdef);

    uint _addDiffBoundary(steps::solver::DiffBoundarydef * dbdef);

    uint _addSDiffBoundary(steps::solver::SDiffBoundarydef * sdbdef);

    void _addTet(tetrahedron_id_t tetidx,
                 steps::mpi::tetopsplit::Comp *comp,
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

    void _addWmVol(uint cidx, steps::mpi::tetopsplit::Comp * comp, double vol);

    void _addTri(triangle_id_t triidx,
                 steps::mpi::tetopsplit::Patch *patch,
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

    void _runWithoutEField(double endtime);
    void _runWithEField(double endtime);
    //void _build();
    void _refreshEFTrisV();

    double _getRate(uint i) const
    { return pKProcs[i]->rate(); }

    steps::mpi::tetopsplit::KProc * _getNext() const;

    //void _reset();

    void _executeStep(steps::mpi::tetopsplit::KProc * kp, double dt, double period = 0.0);
    void _updateSpec(steps::mpi::tetopsplit::WmVol * tet, uint spec_gidx);

    /// Update the kproc's of a triangle, after a species has been changed.
    /// This does not need to update the kproc's of any neighbouring
    /// tetrahedrons.
    ///
    /// Rational: no kproc in tetrahedron depends on species changes on triangle
    ///
    /// Use depSpecTri to check if the kproc depends on the spec, therefore use spec_gidx
    void _updateSpec(steps::mpi::tetopsplit::Tri * tri, uint spec_gidx);


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

    ////////////////////////////////////////////////////////////////////////
    // Batch Data Access
    ////////////////////////////////////////////////////////////////////////

    std::vector<double> getBatchTetCounts(const std::vector<index_t> &tets, std::string const &s) const override;

    std::vector<double> getBatchTriCounts(const std::vector<index_t> &tris, std::string const &s) const override;

    void setBatchTetConcs(const std::vector<index_t> &tets, std::string const &s, const std::vector<double> &concs);

    std::vector<double> getBatchTetConcs(const std::vector<index_t> &tets, std::string const &s) const;

    void getBatchTetCountsNP(const index_t *indices,
                             int input_size,
                             std::string const &s,
                             double *counts,
                             int output_size) const override;

    void getBatchTriCountsNP(const index_t *indices,
                             int input_size,
                             std::string const &s,
                             double *counts,
                             int output_size) const override;

    double sumBatchTetCountsNP(unsigned int* indices, int input_size, std::string const & s);

    double sumBatchTriCountsNP(unsigned int* indices, int input_size, std::string const & s);

    double sumBatchTriGHKIsNP(unsigned int* indices, uint input_size, std::string const & ghk);

    double sumBatchTriOhmicIsNP(unsigned int* indices, uint input_size, std::string const & oc);

    ////////////////////////////////////////////////////////////////////////
    // ROI Data Access
    ////////////////////////////////////////////////////////////////////////

    /// Get species counts of a list of tetrahedrons
    std::vector<double> getROITetCounts(const std::string& ROI_id, std::string const & s) const override;

    /// Get species counts of a list of triangles
    std::vector<double> getROITriCounts(const std::string& ROI_id, std::string const & s) const override;

    /// Get species counts of a list of tetrahedrons
    void getROITetCountsNP(const std::string& ROI_id, std::string const & s, double* counts, int output_size) const override;

    /// Get species counts of a list of triangles
    void getROITriCountsNP(const std::string& ROI_id, std::string const & s, double* counts, int output_size) const override;

    double getROIVol(const std::string& ROI_id) const override;
    double getROIArea(const std::string& ROI_id) const override;

    double getROICount(const std::string& ROI_id, std::string const & s) const override;
    void setROICount(const std::string& ROI_id, std::string const & s, double count) override;

    double getROIAmount(const std::string& ROI_id, std::string const & s) const override;
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


    //////////////////////////// MPI STUFFS ////////////////////////////

    // Exposed to Python:

    void setDiffApplyThreshold(int threshold);

    unsigned long long getReacExtent(bool local = false);
    unsigned long long getDiffExtent(bool local = false);
    double getNIteration();

    double getUpdPeriod() {return updPeriod;}

    void repartitionAndReset(std::vector<uint> const &tet_hosts,
                     std::map<uint, uint> const &tri_hosts  = {},
                     std::vector<uint> const &wm_hosts = {});

    double getCompTime();
    double getSyncTime();
    double getIdleTime();
    double getEFieldTime();
    double getRDTime();
    double getDataExchangeTime();

    // Not currently exposed to Python:
     uint getTetHostRank(uint tidx);
     uint getTriHostRank(uint tidx);
     uint getWMVolHostRank(uint idx);
 	//void registerSyncWmVol(steps::mpi::tetopsplit::WmVol * wmvol);
     //void registerSyncTet(steps::mpi::tetopsplit::Tet * tet);
     //void registerSyncTri(steps::mpi::tetopsplit::Tri * tri);
     void addNeighHost(int host);
     void registerBoundaryTet(steps::mpi::tetopsplit::Tet *tet);
     void registerBoundaryTri(steps::mpi::tetopsplit::Tri *tri);
     uint registerRemoteMoleculeChange(int svol_host,
                                       uint loc,
                                       SubVolType svol_type,
                                       unsigned long idx,
                                       uint slidx,
                                       uint change);

private:



  void setROITetClamped(const std::vector<tetrahedron_id_t>& triangles, const std::string& s, bool b);
  void setROITriClamped(const std::vector<triangle_id_t>& triangles, const std::string& s, bool b);

  double getROITetCount(const std::vector<tetrahedron_id_t>& triangles, const std::string& s) const;
  double getROITriCount(const std::vector<triangle_id_t>& triangles, const std::string& s) const;

  void setROITetCount(const std::vector<tetrahedron_id_t>& triangles, const std::string& s, double count);
  void setROITriCount(const std::vector<triangle_id_t>& triangles, const std::string& s, double count);

    ////////////////////////////////////////////////////////////////////////

    steps::tetmesh::Tetmesh *                    pMesh{nullptr};

    ////////////////////////////////////////////////////////////////////////
    // LIST OF TETEXACT SOLVER OBJECTS
    ////////////////////////////////////////////////////////////////////////

    std::vector<steps::mpi::tetopsplit::Comp *>       pComps;
    std::map<steps::solver::Compdef *, Comp *> pCompMap;

    std::vector<steps::mpi::tetopsplit::Patch *>      pPatches;

    std::vector<steps::mpi::tetopsplit::DiffBoundary *> pDiffBoundaries;

    std::vector<steps::mpi::tetopsplit::SDiffBoundary *> pSDiffBoundaries;

    // These objects are used to describe a mesh compartment that is
    // being treated as a well-mixed volume.
    std::vector<steps::mpi::tetopsplit::WmVol *>      pWmVols;

    std::vector<steps::mpi::tetopsplit::Tri *>        pTris;

    // Now stored as base pointer
    std::vector<steps::mpi::tetopsplit::Tet *>        pTets;

    ////////////////////////////////////////////////////////////////////////
    // Diffusion Data and Methods
    ////////////////////////////////////////////////////////////////////////
    std::vector<Diff*>                          pDiffs;

    void _updateDiff(Diff* diff);

    // separator for non-zero and zero propensity diffusions
    uint                                        diffSep{0};

    ////////////////////////////////////////////////////////////////////////
    // Surface Diffusion Data and Methods
    ////////////////////////////////////////////////////////////////////////
    std::vector<SDiff*>                          pSDiffs;

    void _updateSDiff(SDiff* sdiff);

    // separator for non-zero and zero propensity diffusions
    uint                                        sdiffSep{0};

    ////////////////////////////////////////////////////////////////////////
    // CR SSA Kernel Data and Methods
    ////////////////////////////////////////////////////////////////////////
    uint                                        nEntries;
    double                                      pSum;
    double                                      nSum;
    double                                      pA0{0.0};

    std::vector<KProc*>                         pKProcs;
    std::vector<CRGroup*>                       nGroups;
    std::vector<CRGroup*>                       pGroups;

    ////////////////////////////////////////////////////////////////////////////////

    void _computeUpdPeriod();
    void _updateLocal(std::set<KProc*> const & upd_entries);
    void _updateLocal(std::vector<KProc*> const & upd_entries);
    void _updateLocal(std::vector<uint> const & upd_entries);
    void _updateLocal(uint* upd_entries, uint buffer_size);
    void _updateLocal();
    CRGroup* _getGroup(int pow);
    void _extendPGroups(uint new_size);
    void _extendNGroups(uint new_size);
    void _extendGroup(CRGroup* group, uint size = 1024);
    void _updateSum();
    void _updateElement(KProc* kp);
    ////////////////////////////////////////////////////////////////////////

    // Keeps track of whether _build() has been called
    //bool                                       pBuilt;

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    // The Efield flag. If false we don't calclulate the potential, nor include
    // any voltage-dependent transitions, ohmic or ghk currents. This means
    // the solver behaves exactly like the previous TetOpSplitP solver in this case
    // and all following members are set to null pointer or zero.
    //
    EF_solver                                   pEFoption;

    double                                      pTemp{0.0};

    // Pointer to the EField object
    std::unique_ptr<steps::solver::efield::EField> pEField;

    // The Efield time-step
    double                                      pEFDT{1.0e-5};

    // The number of vertices
    uint                                        pEFNVerts{0};

    // The number of membrane triangles
    uint                                        pEFNTris{0};

    std::vector<steps::mpi::tetopsplit::Tri *>        pEFTris_vec;

    std::vector<double>                         EFTrisV;

    // Working space for gathering distributed computed triangle currents,
    // grouped by rank of owner.
    std::vector<double>                         EFTrisI_permuted;

    // Translate from permuted vector of triangle currents to local EFTri indices.
    std::vector<triangle_id_t >                 EFTrisI_idx;

    // Per-rank counts of EFTris.
    std::vector<int>                            EFTrisI_count;

    // Offsets into ETris_permuted where a rank's owned EFTri data is stored.
    std::vector<int>                            EFTrisI_offset;

    // The number of tetrahedrons
    uint                                        pEFNTets{0};
    // Array of tetrahedrons


    // Table of global vertex index to EField local vertex index (0, 1, ..., pEFNVerts - 1)
    vertex_id_t                               * pEFVert_GtoL{nullptr};

    // Table of global triangle index to EField local triangle index (0, 1, ..., pEFNTris-1)
    triangle_id_t                             * pEFTri_GtoL{nullptr};

    // Table of global tetrahedron index to EField local tet index (0, 1, ..., pEFNTets-1)
    tetrahedron_id_t                          * pEFTet_GtoL{nullptr};

    // Table of EField local triangle index to global triangle index.
    triangle_id_t                             * pEFTri_LtoG{nullptr};

    // True if our copy of tri voltages from EField solver is out of date.
    bool                                        pEFTrisVStale{true};


    ////////////////////////// MPI STUFFS ////////////////////////////

    // TODO TCL: should it be std::vector<tetrahedron_id_t>?
    std::vector<uint>                           tetHosts;
    std::map<triangle_id_t, uint>               triHosts;
    std::vector<uint>                           wmHosts;
    int                                         myRank;
    int                                         nHosts;
    bool                                        recomputeUpdPeriod{true};
    unsigned long long                          reacExtent{0};
    unsigned long long                          diffExtent{0};
    double                                      nIteration{0.0};
    double                                      updPeriod{0.0};
    //bool                                        requireSync;
    uint                                        diffApplyThreshold{10};

    std::set<int>                               neighbHosts;
    uint                                        nNeighbHosts;

    // mirrors of neighboring tets
    std::set<steps::mpi::tetopsplit::Tet *>     boundaryTets;
    std::set<steps::mpi::tetopsplit::Tri *>     boundaryTris;

    std::map<int, std::vector<uint> >           remoteChanges;

    void _remoteSyncAndUpdate(void* requests, std::vector<KProc*> & applied_diffs, std::vector<int> & directions);

    //void _applyRemoteMoleculeChanges(std::vector<MPI_Request> & requests);
    //void _syncPoolCounts();
    //void _updateKProcRates(std::vector<KProc*> & applylist, std::vector<int> & directions, std::vector<MPI_Request> & requests);
    // send, receive and process remote upadtes
    //void _updateRemoteKProcRates(std::vector<MPI_Request> & requests);

    // STL random number generator - also Mersenne twister
    std::random_device                          rd;
    std::mt19937                                gen;

    double                                      compTime{0.0};
    double                                      syncTime{0.0};
    double                                      idleTime{0.0};
    double                                      efieldTime{0.0};
    double                                      rdTime{0.0};
    double                                      dataExchangeTime{0.0};
};

////////////////////////////////////////////////////////////////////////////////

}
}
}

#endif
// STEPS_MPI_TETOPSPLIT_TETAPPROX_HPP

// END
