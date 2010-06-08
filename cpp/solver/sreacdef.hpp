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

#ifndef STEPS_SOLVER_SREACDEF_HPP
#define STEPS_SOLVER_SREACDEF_HPP 1


// STL headers.
#include <string>

// STEPS headers.
#include "../common.h"
#include "statedef.hpp"
#include "api.hpp"
#include "../model/sreac.hpp"
#include "../model/spec.hpp"
#include "types.hpp"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(solver)

// Forwards declarations
class Statedef;

////////////////////////////////////////////////////////////////////////////////

/// Defined Surface Reaction.
/// \todo imcompleted.
class SReacdef
{

public:

    enum orientT  ///< Orientation of the reaction.
    {
        INSIDE = 0,
        OUTSIDE = 1
    };

    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the surface reaction.
    /// \param sr Pointer to the SReac object.
	SReacdef(Statedef * sd, uint idx, steps::model::SReac * sr);

    /// Destructor
	~SReacdef(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE REACTION RULE
    ////////////////////////////////////////////////////////////////////////

	/// Return the global index of this surface reaction rule.
	inline uint gidx(void) const
	{ return pIdx; }

    /// Return the name of the surface reaction.
	std::string const name(void) const;

	/// Return the order of this surface reaction.
	uint order(void) const;

	/// Return the MACROscopic reaction constant.
	double kcst(void) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup the object.
	void setup(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    /// Returns true if the left hand side of the reaction stoichiometry
    /// involves reactants on the surface and on the inside volume.
    inline bool inside(void) const
    { return (pOrient == INSIDE); }

    /// Returns true if any aspect of the surface reaction references
    /// species on the inside volume, regardless of how they are
    /// referenced. Whereas method SReacDef::inside only makes a
    /// statement about the LHS part of the reaction stoichiometry, this
    /// method checks everything, including the right hand side.
    ///
    /// As such, this method will always return true if
    /// SReacDef::inside is true. The converse, however, is not the
    /// the case: SReacDef::inside does not necessarily return true
    /// if this routine returns true.
    ///
    /// It basically polls SReacDef::req_I for each possible species.

    ///
    bool reqInside(void) const;

    /// Returns true if the left hand side of the reaction stoichiometry
    /// involves reactants on the surface and on the outside volume. This
    /// method is mutually exclusive with SReacDef::inside, but not
    /// with SReacDef::insideRef.
    ///
    inline bool outside(void) const
    { return (pOrient == OUTSIDE); }

    /// Returns true if any aspect of the surface reaction references
    /// species on the outside volume, regardless of how they are
    /// referenced. Whereas method SReacDef::outside only makes a
    /// statement about the LHS part of the reaction stoichiometry, this
    /// method checks everything, including the right hand side.
    ///
    /// As such, this method will always return true if
    /// SReacDef::outside is true. The converse, however, is not the
    /// the case: SReacDef::outside does not necessarily return true
    /// if this routine returns true.
    ///
    /// It basically polls SReacDef::req_O for each possible species.
    ///
    bool reqOutside(void) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    /// Returns the number of molecules of species idx required in
    /// the inner volume (_I), outer volume (_O) or surface patch (_S)
    /// to have one occurence of this surface reaction.
    ///
    uint lhs_I(uint gidx) const;
    uint lhs_S(uint gidx) const;
    uint lhs_O(uint gidx) const;

    /// Returns a description of how an occurence of this surface reaction
    /// depends on some species, defined by its global index idx, to occur.
    /// See steps/sim/shared/types.hpp for more information on the return
    /// type. This method is distinct from the SReacDef::req_I,
    /// SReacDef::req_S and SReacDef::req_O methods.
    ///
    depT dep_I(uint gidx) const;
    depT dep_S(uint gidx) const;
    depT dep_O(uint gidx) const;

    /// Returns how many molecules of some species, specified by its
    /// global index, are produced after a single occurence of this
    /// surface reaction. '_I' returns this number for the inner volume,
    /// '_S' for the surface patch and '_O' for the outer volume.
    ///
    uint rhs_I(uint gidx) const;
    uint rhs_S(uint gidx) const;
    uint rhs_O(uint gidx) const;

    /// Returns how the amount of a species, specified by its global index,
    /// changes as the result of a single occurence of this surface
    /// reaction on the inside volume (_I), outer volume (_O) or
    /// surface patch (_S).
    ///
    int upd_I(uint gidx) const;
    int upd_S(uint gidx) const;
    int upd_O(uint gidx) const;

    /// Returns whether the surface reaction rule references a species,
    /// specified by its global index, on the inner volume side (_I),
    /// outer volume (_O) or surface patch (_S).
    ///
    bool reqspec_I(uint gidx) const;
    bool reqspec_S(uint gidx) const;
    bool reqspec_O(uint gidx) const;

    inline gidxTVecCI beginUpdColl_I(void) const
    { return pSpec_I_UPD_Coll.begin(); }
    inline gidxTVecCI endUpdColl_I(void) const
    { return pSpec_I_UPD_Coll.end(); }
    inline gidxTVecCI beginUpdColl_S(void) const
    { return pSpec_S_UPD_Coll.begin(); }
    inline gidxTVecCI endUpdColl_S(void) const
    { return pSpec_S_UPD_Coll.end(); }
    inline gidxTVecCI beginUpdColl_O(void) const
    { return pSpec_O_UPD_Coll.begin(); }
    inline gidxTVecCI endUpdColl_O(void) const
    { return pSpec_O_UPD_Coll.end(); }

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
	steps::model::SpecPVec				pIlhs;
	steps::model::SpecPVec				pOlhs;
	steps::model::SpecPVec				pSlhs;

	steps::model::SpecPVec 				pIrhs;
	steps::model::SpecPVec 				pOrhs;
	steps::model::SpecPVec				pSrhs;

	bool								pSetupdone;


    /// Does the left-hand side of the stoichiometry involve molecules
    /// on the inside or on the outside?
    orientT                             pOrient;

    ////////////////////////////////////////////////////////////////////////
    // DATA: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    /// Array vector describing dependencies for inner volume (_I_),
    /// outer volume (_O_) and surface (_S_) species. Dependencies can
    /// be stoichiometric or in the rate function (this is not implemented
    /// yet) -- see 'steps/sim/shared/types.hpp'. The vector must be
    /// indexed through global indices, i.e. it runs over all species in
    /// the entire model.
    ///
    depT                              * pSpec_I_DEP;
    depT                              * pSpec_S_DEP;
    depT                              * pSpec_O_DEP;

    /// Vector describing the left hand (reactant) side of the reaction
    /// stoichiometry, for species in the inner volume (_I_), outer
    /// volume (_O_) and surface (_S_) species. The vector must be indexed
    /// through global indices, i.e. it runs over all species in the entire
    /// model.
    ///
    uint                              * pSpec_I_LHS;
    uint                              * pSpec_S_LHS;
    uint                              * pSpec_O_LHS;

    /// An array vector describing the right hand (reaction product) side
    /// of the surface reaction stoichiometry, for species in the inner
    /// volume (_I_), outer volume (_O_) and surface (_S_) species. The
    /// vector must be indexed through global indices, i.e. it runs over
    /// all species in the entire model.
    ///
    uint                              * pSpec_I_RHS;
    uint                              * pSpec_S_RHS;
    uint                              * pSpec_O_RHS;

    /// An array describing the update vector (i.e. RHS[] - LHS[]) of
    /// the surface reaction, for species in the inner volume (_I),
    /// outer volume (_O_) and patch surface (_S_). The vector must be
    /// indexed through global indices, i.e. it runs over all species in
    /// the entire model.
    ///
    int                               * pSpec_I_UPD;
    int                               * pSpec_S_UPD;
    int                               * pSpec_O_UPD;

    /// A vector collecting the global indices of all species that are
    /// updated when this surface reaction rule occurs.
    gidxTVec                            pSpec_I_UPD_Coll;
    gidxTVec                            pSpec_S_UPD_Coll;
    gidxTVec                            pSpec_O_UPD_Coll;

    ////////////////////////////////////////////////////////////////////////
};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(solver)
END_NAMESPACE(steps)

#endif
// STEPS_SOLVER_SREACDEF_HPP

// END
