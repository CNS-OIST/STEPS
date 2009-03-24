////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
//
// This file is part of STEPS.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
//
//
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_TETEXACT_TETEXACT_HPP
#define STEPS_TETEXACT_TETEXACT_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>
#include <vector>
#include <set>
#include <map>

// STEPS headers.
#include <steps/common.h>
#include <steps/solver/api.hpp>
#include <steps/solver/statedef.hpp>
#include <steps/geom/tetmesh.hpp>
#include <steps/tetexact/tri.hpp>
#include <steps/tetexact/tet.hpp>
#include <steps/tetexact/kproc.hpp>
#include <steps/tetexact/comp.hpp>
#include <steps/tetexact/patch.hpp>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(tetexact)

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

    Tetexact(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r);
    ~Tetexact(void);

    steps::tetmesh::Tetmesh * mesh(void) const
    { return pMesh; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER INFORMATION
    ////////////////////////////////////////////////////////////////////////

    std::string getSolverName(void) const;
    std::string getSolverDesc(void) const;
    std::string getSolverAuthors(void) const;
    std::string getSolverEmail(void) const;

    ////////////////////////////////////////////////////////////////////////
    // SIMULATION CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////

    void saveState(std::string const & filename);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS
    ////////////////////////////////////////////////////////////////////////

    void reset(void);
    void run(double endtime);
    void step(void);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      GENERAL
    ////////////////////////////////////////////////////////////////////////

    double getTime(void) const;

    inline double getA0(void) const
    { return pA0; }

    uint getNSteps(void) const;

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
	void _setCompDiffD(uint cidx, uint didx);

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

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TETRAHEDRAL VOLUME ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    double _getTetVol(uint tidx) const;
    void _setTetVol(uint tidx, double vol);

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

    double _getTetDiffD(uint tidx, uint didx) const;
    void _setTetDiffD(uint tidx, uint didx);

    bool _getTetDiffActive(uint tidx, uint didx) const;
    void _setTetDiffActive(uint tidx, uint didx, bool act);

    ////////////////////////////////////////////////////////////////////////

    double _getTetReacH(uint tidx, uint ridx) const;
    double _getTetReacC(uint tidx, uint ridx) const;
    double _getTetReacA(uint tidx, uint ridx) const;

    double _getTetDiffA(uint tidx, uint didx) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TRIANGULAR SURFACE ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    double _getTriArea(uint tidx) const;
    void _setTriArea(uint tidx, double area);

    double _getTriCount(uint tidx, uint sidx) const;
    void _setTriCount(uint tidx, uint sidx, double n);

    bool _getTriClamped(uint tidx, uint sidx) const;
    void _setTriClamped(uint tidx, uint sidx, bool buf);

    double _getTriSReacK(uint tidx, uint ridx) const;
    void _setTriSReacK(uint tidx, uint ridx, double kf);

    bool _getTriSReacActive(uint tidx, uint ridx) const;
    void _setTriSReacActive(uint tidx, uint ridx, bool act);

    ////////////////////////////////////////////////////////////////////////

    double _getTriSReacH(uint tidx, uint ridx) const;
    double _getTriSReacC(uint tidx, uint ridx) const;
    double _getTriSReacA(uint tidx, uint ridx) const;

	////////////////////////////////////////////////////////////////////////

	// Called from local Comp or Patch objects. Ass KProc to this object
	void addKProc(steps::tetexact::KProc * kp);

	inline uint countKProcs(void) const
	{ return pKProcs.size(); }

private:

	////////////////////////////////////////////////////////////////////////

	steps::tetmesh::Tetmesh * 			pMesh;

    ////////////////////////////////////////////////////////////////////////
    // TETEXACT SOLVER METHODS
    ////////////////////////////////////////////////////////////////////////

	uint _addComp(steps::solver::Compdef * cdef);

	uint _addPatch(steps::solver::Patchdef * pdef);

	void _addTet(uint tetidx, steps::tetexact::Comp * comp, double vol, double a1,
			     double a2, double a3, double a4, double d1, double d2,
			     double d3, double d4, int tet0, int tet1, int tet2, int tet3);

	void _addTri(uint triidx, steps::tetexact::Patch * patch, double area,
				 int tinner, int touter);

	// called when local tet, tri, reac, sreac objects have been created
	// by constructor
	void _setup(void);

	void _build(void);

	double _getRate(uint i) const
	{ return pLevels[0][i]; }

	steps::tetexact::KProc * _getNext(void) const;

	void _reset(void);

	void _executeStep(steps::tetexact::KProc * kp, double dt);

	void _update(std::vector<uint> const & entries);

    /// Update the kproc's of a tet, after a species has been changed.
    /// This also updates kproc's in surrounding triangles.
    ///
    /// Currently doesn't care about the species.
    ///
    void _updateSpec(steps::tetexact::Tet * tet, uint spec_lidx);

    /// Update the kproc's of a triangle, after a species has been changed.
    /// This does not need to update the kproc's of any neighbouring
    /// tetrahedrons.
    ///
    /// Currently doesn't care about the species.
    ///
    void _updateSpec(steps::tetexact::Tri * tri, uint spec_lidx);

    ////////////////////////////////////////////////////////////////////////
    // LIST OF TETEXACT SOLVER OBJECTS
    ////////////////////////////////////////////////////////////////////////

    std::vector<steps::tetexact::KProc *>      pKProcs;

    std::vector<steps::tetexact::Comp *>       pComps;
    std::map<steps::solver::Compdef *, Comp *> pCompMap;
    inline steps::tetexact::Comp * _comp(uint cidx) const
    { return pComps[cidx]; }

    std::vector<steps::tetexact::Patch *>      pPatches;
    inline steps::tetexact::Patch * _patch(uint pidx) const
    { return pPatches[pidx]; }

    std::vector<steps::tetexact::Tet *>        pTets;

    std::vector<steps::tetexact::Tri *>        pTris;

    ////////////////////////////////////////////////////////////////////////
    // N-ARY TREE
    ////////////////////////////////////////////////////////////////////////

    double                                     pA0;

    std::vector<uint>                          pLevelSizes;

    std::vector<double *>                      pLevels;

	////////////////////////////////////////////////////////////////////////

    // Keeps track of whether _build() has been called
    bool                                       pBuilt;

	////////////////////////////////////////////////////////////////////////

    // Tables to hold update vector indices and random numbers respectively,
    // to be re-used each step.
    uint                                     * pIndices;
    uint 									   pMaxUpSize;
    double                                   * pRannum;

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetexact)
END_NAMESPACE(steps)

#endif
// STEPS_TETEXACT_TETEXACT_HPP

// END
