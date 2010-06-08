////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2010ÊOkinawa Institute of Science and Technology, Japan.
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

#ifndef STEPS_SOLVER_REACDEF_HPP
#define STEPS_SOLVER_REACDEF_HPP 1


// STL headers.
#include <string>

// STEPS headers.
#include "../common.h"
#include "statedef.hpp"
#include "api.hpp"
#include "../model/reac.hpp"
#include "../model/spec.hpp"
#include "types.hpp"

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
	uint                                pIdx;
	std::string 						pName;
	uint 								pOrder;
	double 								pKcst;

	// The stoichiometry stored as model level Spec objects.
	// To be used during setup ONLY
	steps::model::SpecPVec 				pLhs;
	steps::model::SpecPVec 				pRhs;

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
