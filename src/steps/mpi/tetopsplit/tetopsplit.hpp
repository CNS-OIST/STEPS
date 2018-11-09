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
// logging
#include "third_party/easyloggingpp/src/easylogging++.h"
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

    TetOpSplitP(steps::model::Model *m, steps::wm::Geom *g, steps::rng::RNG *r,
            int calcMembPot = EF_NONE, std::vector<uint> const &tet_hosts = std::vector<uint>(),
            std::map<uint, uint> const &tri_hosts = std::map<uint, uint>(),
            std::vector<uint> const &wm_hosts = std::vector<uint>());
    ~TetOpSplitP();


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

    double _getTriSReacK(uint tidx, uint ridx) ;
    void _setTriSReacK(uint tidx, uint ridx, double kf);

    bool _getTriSReacActive(uint tidx, uint ridx) ;
    void _setTriSReacActive(uint tidx, uint ridx, bool act);
    
    double _getTriSDiffD(uint tidx, uint didx, uint direction_tri = std::numeric_limits<uint>::max()) ;
    
    void _setTriSDiffD(uint tidx, uint didx, double dk,
                      uint direction_tri = std::numeric_limits<uint>::max());

    ////////////////////////////////////////////////////////////////////////

    double _getTriSReacH(uint tidx, uint ridx) ;
    double _getTriSReacC(uint tidx, uint ridx) ;
    double _getTriSReacA(uint tidx, uint ridx) ;

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

    bool _getTriVDepSReacActive(uint tidx, uint vsridx) ;
    void _setTriVDepSReacActive(uint tidx, uint vsridx, bool act);

    void _setTriCapac(uint tidx, double cap);

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

    uint addKProc(steps::mpi::tetopsplit::KProc * kp);

    void addDiff(Diff* diff);
    void addSDiff(SDiff* sdiff);
    inline uint countKProcs() const
    { return pKProcs.size(); }

    ////////////////////////////////////////////////////////////////////////

    inline steps::tetmesh::Tetmesh * mesh() const
    { return pMesh; }

    inline steps::mpi::tetopsplit::Comp * _comp(uint cidx) const
    {
        AssertLog(cidx < statedef()->countComps());
        AssertLog(statedef()->countComps() == pComps.size());
        return pComps[cidx];
    }

    inline std::vector<steps::mpi::tetopsplit::Patch *>  patches() const
    { return pPatches; }

    inline steps::mpi::tetopsplit::Patch * _patch(uint pidx) const
    {
        AssertLog(pidx < statedef()->countPatches());
        AssertLog(statedef()->countPatches() == pPatches.size());
        return pPatches[pidx];
    }

    inline steps::mpi::tetopsplit::DiffBoundary * _diffboundary(uint dbidx) const
    {
        AssertLog(dbidx < statedef()->countDiffBoundaries());
        return pDiffBoundaries[dbidx];
    }

    inline steps::mpi::tetopsplit::SDiffBoundary * _sdiffboundary(uint sdbidx) const
    {
        AssertLog(sdbidx < statedef()->countSDiffBoundaries());
        return pSDiffBoundaries[sdbidx];
    }

    inline steps::mpi::tetopsplit::Tet * _tet(uint tidx) const
    { return pTets[tidx]; }

    inline steps::mpi::tetopsplit::Tri * _tri(uint tidx) const
    { return pTris[tidx]; }

    inline double a0() const
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

    void _addTet(uint tetidx, steps::mpi::tetopsplit::Comp * comp, double vol, double a1,
                 double a2, double a3, double a4, double d1, double d2,
                 double d3, double d4, int tet0, int tet1, int tet2, int tet3);

    void _addWmVol(uint cidx, steps::mpi::tetopsplit::Comp * comp, double vol);

    void _addTri(uint triidx, steps::mpi::tetopsplit::Patch * patch, double area,
            double l0, double l1, double l2, double d0, double d1, double d2,
            int tinner, int touter, int tri0, int tri1, int tri2);

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
    // Batch Data Access
    ////////////////////////////////////////////////////////////////////////
    
    std::vector<double> getBatchTetCounts(std::vector<uint> const & tets, std::string const & s) const;
    
    std::vector<double> getBatchTriCounts(std::vector<uint> const & tris, std::string const & s) const;

    void getBatchTetCountsNP(unsigned int* indices, int input_size, std::string const & s, double* counts, int output_size) const;
    
    void getBatchTriCountsNP(unsigned int* indices, int input_size, std::string const & s, double* counts, int output_size) const;
    
    double sumBatchTetCountsNP(unsigned int* indices, int input_size, std::string const & s);
    
    double sumBatchTriCountsNP(unsigned int* indices, int input_size, std::string const & s);
    
    double sumBatchTriGHKIsNP(unsigned int* indices, int input_size, std::string const & ghk);
    
    double sumBatchTriOhmicIsNP(unsigned int* indices, int input_size, std::string const & oc);
    
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

    
    //////////////////////////// MPI STUFFS ////////////////////////////
    
    // Exposed to Python:
    
    void setDiffApplyThreshold(int threshold);

    double getReacExtent(bool local = false);
    double getDiffExtent(bool local = false);
    double getNIteration();
    
    double getUpdPeriod() {return updPeriod;}
    
    void repartitionAndReset(std::vector<uint> const &tet_hosts,
                     std::map<uint, uint> const &tri_hosts  = std::map<uint, uint>(),
                     std::vector<uint> const &wm_hosts = std::vector<uint>());
    
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
     uint registerRemoteMoleculeChange(int svol_host, uint loc, SubVolType svol_type, uint idx, uint slidx, uint change);

private:

    ////////////////////////////////////////////////////////////////////////

    steps::tetmesh::Tetmesh *                    pMesh;

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
    uint                                        diffSep;
    
    ////////////////////////////////////////////////////////////////////////
    // Surface Diffusion Data and Methods
    ////////////////////////////////////////////////////////////////////////
    std::vector<SDiff*>                          pSDiffs;
    
    void _updateSDiff(SDiff* sdiff);
    
    // separator for non-zero and zero propensity diffusions
    uint                                        sdiffSep;

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

    double                                      pTemp;

    // Pointer to the EField object
    std::unique_ptr<steps::solver::efield::EField> pEField;

    // The Efield time-step
    double                                      pEFDT;

    // The number of vertices
    uint                                        pEFNVerts;

    // The number of membrane triangles
    uint                                        pEFNTris;

    std::vector<steps::mpi::tetopsplit::Tri *>        pEFTris_vec;
    
    std::vector<double>                         EFTrisV;

    // Working space for gathering distributed computed triangle currents,
    // grouped by rank of owner.
    std::vector<double>                         EFTrisI_permuted;

    // Translate from permuted vector of triangle currents to local EFTri indices.
    std::vector<int>                            EFTrisI_idx;

    // Per-rank counts of EFTris.
    std::vector<int>                            EFTrisI_count;

    // Offsets into ETris_permuted where a rank's owned EFTri data is stored.
    std::vector<int>                            EFTrisI_offset;

    // The number of tetrahedrons
    uint                                        pEFNTets;
    // Array of tetrahedrons
    

    // Table of global vertex index to EField local vertex index (0, 1, ..., pEFNVerts - 1)
    int                                      * pEFVert_GtoL;

    // Table of global triangle index to EField local triangle index (0, 1, ..., pEFNTris-1)
    int                                      * pEFTri_GtoL;

    // Table of global tetrahedron index to EField local tet index (0, 1, ..., pEFNTets-1)
    int                                       * pEFTet_GtoL;

    // Table of EField local triangle index to global triangle index.
    uint                                      * pEFTri_LtoG;
    
    // True if our copy of tri voltages from EField solver is out of date.
    bool                                        pEFTrisVStale;

    
    ////////////////////////// MPI STUFFS ////////////////////////////
    
    std::vector<uint>                           tetHosts;
    std::map<uint, uint>                        triHosts;
    std::vector<uint>                           wmHosts;
    int                                         myRank;
    int                                         nHosts;
    double                                      diffExtent;
    double                                      reacExtent;
    double                                      nIteration;
    double                                      updPeriod;
    bool                                        recomputeUpdPeriod;
    //bool                                        requireSync;
    uint                                        diffApplyThreshold;
    
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
    
    double                                      compTime;
    double                                      syncTime;
    double                                      idleTime;
    double                                      efieldTime;
    double                                      rdTime;
    double                                      dataExchangeTime;
};

////////////////////////////////////////////////////////////////////////////////

}
}
}

#endif
// STEPS_MPI_TETOPSPLIT_TETAPPROX_HPP

// END
