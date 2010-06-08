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

#ifndef STEPS_SOLVER_PATCHDEF_HPP
#define STEPS_SOLVER_PATCHDEF_HPP 1


// STL headers.
#include <string>

// STEPS headers.
#include "../common.h"
#include "statedef.hpp"
#include "api.hpp"
#include "../geom/patch.hpp"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(solver)

// Forwards declarations
class Statedef;
class SReacdef;
class Compdef;

// Auxiliary declarations.
typedef Patchdef *                      PatchDefP;
typedef std::vector<PatchDefP>          PatchDefPVec;
typedef PatchDefPVec::iterator          PatchDefPVecI;
typedef PatchDefPVec::const_iterator    PatchDefPVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Defined patch object.
class Patchdef
{

public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the patch.
    /// \param p Pointer to the Patch object.
	Patchdef(Statedef * sd, uint idx, steps::wm::Patch * p);

    /// Destructor
	~Patchdef(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: PATCH
    ////////////////////////////////////////////////////////////////////////

	/// Return the area of this patch.
	double area(void) const;

	/// Return the global index of this patch.
	inline uint gidx(void) const
	{ return pIdx; }

    /// Return the name of the patch.
	std::string const name(void) const;

    /// Return a pointer to the inner compartment.
	inline Compdef * icompdef(void) const
	{ return pInner; }

    /// Return a pointer to the outer compartment.
	inline Compdef * ocompdef(void) const
	{ return pOuter; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup all references.
	void setup_references(void);

    /// Setup all indices.
	void setup_indices(void);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: PATCH
    ////////////////////////////////////////////////////////////////////////

	/// Set the area of this patch.
    ///
    /// \param a Area of the patch.
	void setArea(double a);

	/// Reset count, flags members of this patch. Called when reset()
	/// method in solver object is executed.
    ///
	void reset(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the number of species defined for this surface patch.
    inline uint countSpecs(void) const
    { return pSpecsN_S; }

    /// Return the number of species defined for the inner compartment.
    /// Should not be called before Compdef::setup()
    inline uint countSpecs_I(void) const
    { return pSpecsN_I; }

    /// Return the number of species defined for the outer compartment.
    /// Should not be called before Compdef::setup()
    inline uint countSpecs_O(void) const
    { return pSpecsN_O; }

    /// Return the local species index for global index argument.
    ///
    /// \param gidx Global index of the species.
    inline uint specG2L(uint gidx) const
    { return pSpec_G2L[gidx]; }

    /*										// don't need? Just used internally, instead of using member??
    inline uint specL2G(uint lidx) const
    { return pSpec_L2G[lidx]; }
    */

    /// Auxiliary function: resolves a species gidx for the inner
    /// compartment.
    ///
    /// \param gidx Global index of the species.
    /// \return The local index or steps::sim::shared::LIDX_UNDEFINED.
    ///
    uint specG2L_I(uint gidx) const;
    /// Auxiliary function: resolves a species gidx for the outer
    /// compartment.
    ///
    /// \param gidx Global index of the species.
    /// \return The local index or steps::sim::shared::LIDX_UNDEFINED.
    ///
    uint specG2L_O(uint gidx) const;

	/// Return pointer to species' counts on this patch.
	inline double * pools(void) const
	{ return pPoolCount; }

	/// Returns pointer to flags on species for this patch.
	inline uint * flags(void) const
	{ return pPoolFlags; }

	static const uint CLAMPED = 1;

	/// Return whether a species, specified by local index argument, is
	/// clamped or not.
    ///
    /// \param slidx Local index of the species.
	inline bool clamped(uint slidx) const
	{ return pPoolFlags[slidx] & CLAMPED; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SPECIES
    ////////////////////////////////////////////////////////////////////////

	/// Set the species count of species specified by local index argument.
    ///
    /// \param slidx Local index of the species.
    /// \param count Count of species.
	void setCount(uint slidx, double count);

	/// Clamp or unclamp species specified by local index argument
    ///
    /// \param slidx Local index of the species.
    /// \param clamp Flag to clamp or unclamp species.
	void setClamped(uint slidx, bool clamp);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE REACTION RULES
    ////////////////////////////////////////////////////////////////////////

	/// Return the total number of surface reactions that can occur on this
	/// patch.
    inline uint countSReacs(void) const
    { return pSReacsN; }

	/// Return the local surface reaction index for global index argument.
    ///
    /// \param gidx Global index of the surface reaction.
    inline uint sreacG2L(uint gidx) const
    { return pSReac_G2L[gidx]; }
    /*
    inline gidxT sreacL2G(lidxT idx) const
    { return pSReac_L2G[idx]; }
	*/

	/// Return a pointer to reaction definition object (type Reacdef)
	/// specified by local index.
    SReacdef * sreacdef(uint lidx) const;

    /// Warning: these methods perform no error checking!
    /// \todo imcompleted.
    int sreac_dep_I(uint srlidx, uint splidx) const;
    int sreac_dep_S(uint srlidx, uint splidx) const;
    int sreac_dep_O(uint srlidx, uint splidx) const;

    /// Warning: these methods perform no error checking!
    ///
	// Return the beginning and end of the lhs arrays of surface reaction
    // specified by local index argument.
    uint * sreac_lhs_I_bgn(uint lidx) const;
    uint * sreac_lhs_I_end(uint lidx) const;
    uint * sreac_lhs_S_bgn(uint lidx) const;
    uint * sreac_lhs_S_end(uint lidx) const;
    uint * sreac_lhs_O_bgn(uint lidx) const;
    uint * sreac_lhs_O_end(uint lidx) const;

    /// Warning: these methods perform no error checking!
    ///
	// Return the beginning and end of the update arrays of surface reaction
    // specified by local index argument.
    int * sreac_upd_I_bgn(uint lidx) const;
    int * sreac_upd_I_end(uint lidx) const;
    int * sreac_upd_S_bgn(uint lidx) const;
    int * sreac_upd_S_end(uint lidx) const;
    int * sreac_upd_O_bgn(uint lidx) const;
    int * sreac_upd_O_end(uint lidx) const;

	/// Return pointer to flags on surface reactions for this patch.
	inline uint * srflags(void) const
	{ return pSReacFlags; }

	static const uint INACTIVATED = 1;

	/// Return whether a surface reaction, specified by local index argument,
	/// is active or not.
    ///
    /// \param rlidx Local index of the surface reaction.
	inline bool active(uint rlidx) const
	{ return !(pSReacFlags[rlidx] & INACTIVATED); }

	/// Return the kcst for a surface reaction specified by local index
    ///
    /// \param rlidx Local index of the surface reaction.
	inline double kcst(uint rlidx) const
	{ return pSReacKcst[rlidx];

	}
    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: REACTIONS
    ////////////////////////////////////////////////////////////////////////

	/// Set the kcst for surface reaction specified by local index
    ///
    /// \param srlidx Local index of the surface reaction.
    /// \param kcst Rate constant of the surface reaction.
	void setKcst(uint srlidx, double kcst);

	/// Activate or inactivate a surface reaction specified by local index.
    ///
    /// \param srlidx Local index of the surface reaction.
    /// \param active Flag to activate / inactivate the surface reaction.
	void setActive(uint srlidx, bool active);

private:

    ////////////////////////////////////////////////////////////////////////
    // DATA: GENERAL
    ////////////////////////////////////////////////////////////////////////

    /// A pointer to the state definition.
    Statedef                          * pStatedef;

	// The string identifier of the patch
	std::string                         pName;

	// The area of the patch
	double                              pArea;

    /// The index of the patch.
    uint                                pIdx;

	// The enclosed surface systems, stored as strings
	std::set<std::string> 				pPssys;

	// Pointers to geom level comps, to be used during setup and NOT later.
	steps::wm::Comp 				  * pIcomp;
	steps::wm::Comp                   * pOcomp;

    /// Pointer to inner compartment CompDef.
    Compdef                           * pInner;

    /// Pointer to outer compartment CompDef.
    Compdef                           * pOuter;

	// Keep track of whether setup methods have been called.
	bool                                pSetupRefsdone;
	bool                                pSetupIndsdone;

    ////////////////////////////////////////////////////////////////////////
    // DATA: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Number of species embedded in inner volume (_I), patch (_S)
    /// and outer volume (_O).
    uint                                pSpecsN_I;
    uint                                pSpecsN_S;
    uint                                pSpecsN_O;

    /// Table to resolve species index (global -> local).
    uint                              * pSpec_G2L;

    /// Table to resolve species index (local -> global).
    uint                              * pSpec_L2G;

	// Table of the populations of the species on this patch.
	double                            * pPoolCount;

	// Table of 'clamped' flags on the species
	uint                              * pPoolFlags;

    ////////////////////////////////////////////////////////////////////////
    // DATA: SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Number of surface reactions occuring in patch.
    uint                                pSReacsN;

    /// Table to resolve reaction rule indices (global -> local).
    uint                              * pSReac_G2L;

    /// Table to resolve reaction rule indices (local -> global).
    uint                              * pSReac_L2G;

	// Table of the K-constants of the surface reac rules in this patch
	double                            * pSReacKcst;

	// Table of 'active' flags on the surface reaction rules.
	uint                              * pSReacFlags;

    inline uint _IDX_SReac_I_Spec(uint srlidx)
    { return countSpecs_I() * srlidx; }
    inline uint _IDX_SReac_I_Spec(uint srlidx, uint splidx)
    { return (countSpecs_I() * srlidx) + splidx; }
    inline uint _IDX_SReac_S_Spec(uint srlidx)
    { return countSpecs() * srlidx; }
    inline uint _IDX_SReac_S_Spec(uint srlidx, uint splidx)
    { return (countSpecs() * srlidx) + splidx; }
    inline uint _IDX_SReac_O_Spec(uint srlidx)
    { return countSpecs_O() * srlidx; }
    inline uint _IDX_SReac_O_Spec(uint srlidx, uint splidx)
    { return (countSpecs_O() * srlidx) + splidx; }

    int                               * pSReac_DEP_I_Spec;
    int                               * pSReac_DEP_S_Spec;
    int                               * pSReac_DEP_O_Spec;
    uint                              * pSReac_LHS_I_Spec;
    uint                              * pSReac_LHS_S_Spec;
    uint                              * pSReac_LHS_O_Spec;
    int                               * pSReac_UPD_I_Spec;
    int                               * pSReac_UPD_S_Spec;
    int                               * pSReac_UPD_O_Spec;

    ////////////////////////////////////////////////////////////////////////
};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(solver)
END_NAMESPACE(steps)

#endif
// STEPS_SOLVER_PATCHDEF_HPP

// END
