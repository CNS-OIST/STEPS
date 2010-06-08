////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */

#ifndef STEPS_SOLVER_WMDIRECT_HPP
#define STEPS_SOLVER_WMDIRECT_HPP 1


// STL headers.
#include <string>
#include <vector>
#include <set>
#include <map>

// STEPS headers.
#include "../common.h"
#include "../solver/api.hpp"
#include "../solver/statedef.hpp"
#include "../solver/compdef.hpp"
#include "../solver/patchdef.hpp"
#include "comp.hpp"
#include "patch.hpp"
#include "kproc.hpp"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(wmdirect)

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

class Wmdirect: public steps::solver::API
{

public:

    Wmdirect(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r);
    ~Wmdirect(void);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER INFORMATION
    ////////////////////////////////////////////////////////////////////////

    std::string getSolverName(void) const;
    std::string getSolverDesc(void) const;
    std::string getSolverAuthors(void) const;
    std::string getSolverEmail(void) const;


    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS
    ////////////////////////////////////////////////////////////////////////

    void reset(void);
    void run(double endtime);
    void advance(double adv);
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
	void _setCompVol(uint cidx, double vol);

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

	double _getCompReacC(uint cidx, uint ridx) const;

	double _getCompReacH(uint cidx, uint ridx) const;

	double _getCompReacA(uint cidx, uint ridx) const;

	uint _getCompReacExtent(uint cidx, uint ridx) const;
 	void _resetCompReacExtent(uint cidx, uint ridx);

    ////////////////////////////////////////////////////////////////////////
/*
	double _getCompDiffD(uint cidx, uint didx);
	void _setCompDiffD(uint cidx, uint didx);

	bool _getCompDiffActive(uint cidx, uint didx);
	void _setCompDiffActive(uint cidx, uint didx, bool act);
*/

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      PATCH
    ////////////////////////////////////////////////////////////////////////

	double _getPatchArea(uint pidx) const;
	void _setPatchArea(uint pidx, double area);

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

 	double _getPatchSReacC(uint pidx, uint ridx) const;

 	double _getPatchSReacH(uint pidx, uint ridx) const;

 	double _getPatchSReacA(uint pidx, uint ridx) const;

 	uint _getPatchSReacExtent(uint pidx, uint ridx) const;
	void _resetPatchSReacExtent(uint pidx, uint ridx);

	////////////////////////////////////////////////////////////////////////

	// Called from local Comp or Patch objects. Ass KProc to this object
	void addKProc(steps::wmdirect::KProc * kp);

	inline uint countKProcs(void) const
	{ return pKProcs.size(); }

	////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////
    // WMDIRECT SOLVER METHODS
    ////////////////////////////////////////////////////////////////////////

	uint _addComp(steps::solver::Compdef * cdef);

	uint _addPatch(steps::solver::Patchdef * pdef);

	// called when local comp, patch, reac, sreac objects have been created
	// by constructor
	void _setup(void);

	void _build(void);

	steps::wmdirect::KProc * _getNext(void) const;

	void _reset(void);

	void _update(SchedIDXVec const & entries);

	void _executeStep(steps::wmdirect::KProc * kp, double dt);

    ////////////////////////////////////////////////////////////////////////
    // LIST OF WMDIRECT SOLVER OBJECTS
    ////////////////////////////////////////////////////////////////////////

    std::vector<steps::wmdirect::KProc *>      pKProcs;

    std::vector<steps::wmdirect::Comp *>       pComps;
    std::map<steps::solver::Compdef *, Comp *> pCompMap;

    std::vector<steps::wmdirect::Patch *>      pPatches;

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

	////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(wmdirect)
END_NAMESPACE(steps)

#endif
// STEPS_SOLVER_WMDIRECT_HPP

// END
