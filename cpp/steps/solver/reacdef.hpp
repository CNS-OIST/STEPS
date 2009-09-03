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

#ifndef STEPS_SOLVER_REACDEF_HPP
#define STEPS_SOLVER_REACDEF_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>

// STEPS headers.
#include <steps/common.h>
#include <steps/solver/statedef.hpp>
#include <steps/solver/api.hpp>
#include <steps/model/reac.hpp>
#include <steps/solver/types.hpp>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(solver)

// Forwards declarations
class Statedef;

////////////////////////////////////////////////////////////////////////////////

/// Defined Reaction.
class Reacdef
{

public:
    /// Constructor
    /// 
    /// \param sd State of the solver.
    /// \param idx Global index of the reaction.
    /// \param r Pointer to the associated Reac object.
	Reacdef(Statedef * sd, uint idx, steps::model::Reac * r);

    /// Destructor
	~Reacdef(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: REACTION RULE
    ////////////////////////////////////////////////////////////////////////

	/// Return the global index of this reaction rule.
	inline uint gidx(void) const
	{ return pIdx; }

    /// Return the name of the reaction.
	std::string const name(void) const;

	/// Return the order of this reaction.
	uint order(void) const;

	/// Return the MACROscopic reaction constant.
	double kcst(void) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////
    /// \todo imcompleted.
    uint lhs(uint gidx) const;
    int dep(uint gidx) const;
    uint rhs(uint gidx) const;
    int upd(uint gidx) const;
    bool reqspec(uint gidx) const;

    inline steps::solver::gidxTVecCI bgnUpdColl(void) const
    { return pSpec_UPD_Coll.begin(); }
    inline steps::solver::gidxTVecCI endUpdColl(void) const
    { return pSpec_UPD_Coll.end(); }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////
    /// Setup the object.
	void setup(void);

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

	Statedef                          * pStatedef;
	steps::model::Reac                * pReac;
	uint                                pIdx;
	bool								pSetupdone;

    ////////////////////////////////////////////////////////////////////////
    // DATA: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    int                               * pSpec_DEP;
    uint                              * pSpec_LHS;
    uint                              * pSpec_RHS;
    int                               * pSpec_UPD;
    steps::solver::gidxTVec             pSpec_UPD_Coll;

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(solver)
END_NAMESPACE(steps)

#endif
// STEPS_SOLVER_REACDEF_HPP

// END
