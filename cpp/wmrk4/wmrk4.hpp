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

#ifndef STEPS_SOLVER_WMRK4_HPP
#define STEPS_SOLVER_WMRK4_HPP 1


// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include "../common.h"
#include "../solver/api.hpp"
#include "../solver/statedef.hpp"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(wmrk4)
USING_NAMESPACE(steps::solver);

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.

// Auxiliary declarations.
typedef std::vector<double>             dVec;
typedef dVec::iterator					dVecI;

typedef std::vector<uint>				uiVec;
typedef uiVec::iterator					uiVecI;

////////////////////////////////////////////////////////////////////////////////

class Wmrk4: public API
{

public:

    Wmrk4(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r);
    ~Wmrk4(void);

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

    void setDT(double dt);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      GENERAL
    ////////////////////////////////////////////////////////////////////////

    double getTime(void) const;

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

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////
    // WMRK4 SOLVER METHODS
    ////////////////////////////////////////////////////////////////////////

	/// initialises vectors and builds the reaction matrix.
	///
	void _setup(void);

	/// this function refills the values and flags vectors
	/// called if a flag or a count changed
	///
	void _refill(void);

	/// refill the Ccst vectors
	/// called if reaction constants are changed during simulation
	///
	void _refillCcst(void);

	/// returns properly scaled reaction constant
	///
	double _ccst(double kcst, double vol, uint order);

	/// the Runge-Kutta algorithm
	///
	void _rk4(double pdt);

	/// the simple stepper
	///
	void _rksteps(double t1, double t2);

	/// the derivatives calculator
	///
	void _setderivs(dVec& vals, dVec& dydx);

	/// update local values vector,
	/// then update state with computed counts
	///
	void _update(void);

    ////////////////////////////////////////////////////////////////////////
    // WMRK4 SOLVER MEMBERS
    ////////////////////////////////////////////////////////////////////////

	/// the reaction matrix
	uint **						        pReacMtx;

	/// update matrix of rhs - lhs
	int **						        pUpdMtx;

	/// number of species total: all species in all comps and patches
	uint						        pSpecs_tot;

	/// number of reactions total: all reactions and surface reactions
	/// in each comp and patch
	uint						        pReacs_tot;

	/// Properly scaled reaction constant vector
	dVec						        pCcst;

	/// vector holding current molecular counts (as doubles)
	dVec						        pVals;

	/// vector holding flags on species
	uiVec						        pSFlags;

	/// vector holding flags on reactions
	uiVec						        pRFlags;

	/// vector holding new, calculated counts
	dVec						        pNewVals;

	/// vector of present derivatives
	dVec						        pDyDx;

	/// matrix of derivatives from each reaction for each species
	double **					        pDyDxlhs;

	/// the time step
	double						        pDT;

	/// objects to contain temporary values important in algorithm
	dVec						        yt;
	dVec						        dyt;
	dVec						        dym;

	////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(wmrk4)
END_NAMESPACE(steps)

#endif
// STEPS_SOLVER_WMRK4_HPP

// END
