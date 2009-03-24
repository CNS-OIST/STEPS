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

#ifndef STEPS_SOLVER_COMPDEF_HPP
#define STEPS_SOLVER_COMPDEF_HPP 1

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
#include <steps/geom/comp.hpp>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(solver)

// Forwards declarations
class Statedef;
class Reacdef;
class Patchdef;
class Compdef;

// Auxiliary declarations.
typedef Compdef *                       CompDefP;
typedef std::vector<CompDefP>           CompDefPVec;
typedef CompDefPVec::iterator           CompDefPVecI;
typedef CompDefPVec::const_iterator     CompDefPVecCI;

////////////////////////////////////////////////////////////////////////////////

class Compdef
{

public:
	Compdef(Statedef * sd, uint idx, steps::wm::Comp * c);

	~Compdef(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

	// Return the volume of this compartment
	double vol(void) const;

	// Return the global index of this compartment
	inline uint gidx(void) const
	{ return pIdx; }

	std::string const name(void) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

	void setup_references(void);

	void setup_indices(void);

	void addSpec(uint gidx);

	void addIPatchdef(steps::solver::Patchdef * p);
	void addOPatchdef(steps::solver::Patchdef * p);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

	// Set the volume of the compartment
	void setVol(double v);

	// Reset count, flags members of this compartment. Called when reset()
	// method in solver object is executed.
	void reset(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////

	// Return the total number of species that can occur in this compartment
	inline uint countSpecs(void) const
	{ return pSpecsN; }

	// Return the local species index for global index argument.
	inline uint specG2L(uint gidx) const
	{ return pSpec_G2L[gidx]; }

	// Returns pointer to flags on species for this compartment.
	inline uint * flags(void) const
	{ return pPoolFlags; }

	static const uint CLAMPED = 1;

	// Return whether a species, specified by local index argument, is
	// clamped or not
	inline bool clamped(uint slidx) const
	{ return pPoolFlags[slidx] & CLAMPED; }

	// Return pointer to species' counts in this compartment.
	inline double * pools(void) const
	{ return pPoolCount; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SPECIES
    ////////////////////////////////////////////////////////////////////////

	// Set the species count of species specified by local index argument.
	void setCount(uint slidx, double count);

	// Clamp or unclamp species specified by local index argument
	void setClamped(uint slidx, bool clamp);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: REACTIONS
    ////////////////////////////////////////////////////////////////////////

	// Return the total number of reactions that can occur in this compartment.
	inline uint countReacs(void) const
	{ return pReacsN; }

	// Return the local reaction index for global index argument.
	uint reacG2L(uint gidx) const
	{ return pReac_G2L[gidx]; }

	// Return the beginning of the lhs array of reaction specified by local
	// index argument
	uint * reac_lhs_bgn(uint rlidx) const;

	uint * reac_lhs_end(uint rlidx) const;

	// Return the beginning of the update array of reaction specified by
	// local index argument.
	int * reac_upd_bgn(uint rlidx) const;

	int * reac_upd_end(uint rlidx) const;

	int reac_dep(uint rlidx, uint slidx) const;

	// Return a pointer to reaction definition object (type Reacdef)
	// specified by local index.
	Reacdef * reacdef(uint rlidx) const;

	// Return pointer to flags on reactions for this compartment.
	inline uint * rflags(void) const
	{ return pReacFlags; }

	static const uint INACTIVATED = 1;

	// Return whether a reaction, specified by local index argument, is
	// active or not.
	inline bool active(uint rlidx) const
	{ return !(pReacFlags[rlidx] & INACTIVATED); }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: REACTIONS
    ////////////////////////////////////////////////////////////////////////

	// Activate or inactivate a reaction specified by local index argument.
	void setActive(uint rlidx, bool active);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION
    ////////////////////////////////////////////////////////////////////////

	// Return the total number of diffusion rules for this compartment.
	inline uint countDiffs(void) const
	{ return pDiffsN; }

	// Returns a pointer to Diffdef specified by local index.
	Diffdef * diffdef(uint dlidx) const;

	// Return the local diffusion index for global index argument.
	uint diffG2L(uint gidx) const
	{ return pDiff_G2L[gidx]; }

	uint diff_dep(uint dlidx, uint slidx) const;

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////
    // DATA: GENERAL
    ////////////////////////////////////////////////////////////////////////

	// A pointer to the state definition.
	Statedef                          * pStatedef;

	// A pointer to the steps:wm::Comp object this object defines.
	steps::wm::Comp                   * pComp;

	// The global index of the compartment.
	uint                                pIdx;

	// Keep track of whether setup_ has been called.
	bool                                pSetupRefsdone;
	bool                                pSetupIndsdone;

    ////////////////////////////////////////////////////////////////////////
    // DATA: SPECIES
	////////////////////////////////////////////////////////////////////////

	// Number of species that can appear in this compartment.
	uint                                pSpecsN;

	// Table to resolve species index (global -> local).
	uint                              * pSpec_G2L;

	// Table to resolve species index (local -> global).
	uint                              * pSpec_L2G;

	// Table of the populations of the species in this compartment.
	double                            * pPoolCount;

	// Table of 'clamped' flags on the species
	uint                              * pPoolFlags;

    ////////////////////////////////////////////////////////////////////////
    // DATA: REACTION RULES
    ////////////////////////////////////////////////////////////////////////

	// Number of reaction rules occuring in this compartment.
	uint                                pReacsN;

    // Table to resolve reaction rule indices (global -> local).
	uint                              * pReac_G2L;

    // Table to resolve reaction rule indices (local -> global).
	uint                              * pReac_L2G;

	// Table of 'active' flags on the reaction rules.
	uint                              * pReacFlags;

	inline uint _IDX_Reac_Spec(uint reac, uint spec) const
	{ return (pSpecsN * reac) + spec; }

	int                               * pReac_DEP_Spec;
	uint                              * pReac_LHS_Spec;
	int                               * pReac_UPD_Spec;

    ////////////////////////////////////////////////////////////////////////
    // DATA: DIFFUSION RULES
    ////////////////////////////////////////////////////////////////////////

    // Number of diffusion rules occuring in compartment.
	uint                                pDiffsN;

    // Table to resolve diffusion rule indices (global -> local).
	uint                              * pDiff_G2L;

    // Table to resolve diffusion rule indices (local -> global).
	uint                              * pDiff_L2G;

    inline uint _IDX_Diff_Spec(uint diff, uint spec) const
    { return (pSpecsN * diff) + spec; }
	uint                              * pDiff_DEP_Spec;
	uint                              * pDiff_LIG;

    ////////////////////////////////////////////////////////////////////////
    // DATA: CONNECTIVITY
    ////////////////////////////////////////////////////////////////////////

    // A list of patches surrounding
    std::vector<steps::solver::Patchdef *>   pIPatches;
    std::vector<steps::solver::Patchdef *>   pOPatches;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(solver)
END_NAMESPACE(steps)

#endif
// STEPS_SOLVER_COMPDEF_HPP

// END
