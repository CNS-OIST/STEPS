////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
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
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SIM_SHARED_SREACDEF_HPP
#define STEPS_SIM_SHARED_SREACDEF_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/types.hpp>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(sim)

// Forward declarations.
class SReacDef;
class StateDef;

// Auxiliary declarations.
typedef SReacDef *                      SReacDefP;
typedef std::vector<SReacDefP>          SReacDefPVec;
typedef SReacDefPVec::iterator          SReacDefPVecI;
typedef SReacDefPVec::const_iterator    SReacDefPVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Stores the information related to a given surface reaction rule.
///
/// Objects of this class should only be created and added to a StateDef 
/// object, after all possible species in the system have already been 
/// defined and added to the StateDef object.

class SReacDef
{

public:

    enum orientT
    {
        INSIDE = 0,
        OUTSIDE = 1
    };
    
    SReacDef(StateDef * sdef, gidxT idx, std::string const & name, orientT o);
    ~SReacDef(void);
    
    ////////////////////////////////////////////////////////////////////////
    // SREACDEF SETUP
    ////////////////////////////////////////////////////////////////////////
    
    uint incLHS_I(gidxT idx);
    uint incLHS_S(gidxT idx);
    uint incLHS_O(gidxT idx);
    void setLHS_I(gidxT idx, uint n);
    void setLHS_S(gidxT idx, uint n);
    void setLHS_O(gidxT idx, uint n);
    
    uint incRHS_I(gidxT idx);
    uint incRHS_S(gidxT idx);
    uint incRHS_O(gidxT idx);
    void setRHS_I(gidxT idx, uint n);
    void setRHS_S(gidxT idx, uint n);
    void setRHS_O(gidxT idx, uint n);
    
    /// Gets called when the definition of all components in the entire
    /// state has finished.
    ///
    /// Currently, this method does the following things:
    /// <UL>
    /// <LI> Computes update (UPD) and dependency (DEP) vectors.
    /// <LI> Computes the "collected" updated species.
    /// <LI> Re-computes order.
    /// </UL>
    ///
    void setupFinal(void);
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: GENERAL
    ////////////////////////////////////////////////////////////////////////
    
    inline StateDef * statedef(void) const
    { return pStateDef; }
    
    inline gidxT gidx(void) const
    { return pGIDX; }
    
    inline std::string const & name(void) const
    { return pName; }

    ////////////////////////////////////////////////////////////////////////

    /// Returns true if the left hand side of the reaction stoichiometry
    /// involves reactants on the surface and on the inside volume. This
    /// method is mutually exclusive with SReacDef::outside, but not
    /// with SReacDef::reqOutside.
    ///
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
    /// This method should only be called after SReacDef::setupFinal has
    /// been called.
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
    /// This method should only be called after SReacDef::setupFinal has
    /// been called.
    ///
    bool reqOutside(void) const;
    
    inline uint order(void) const
    { return pOrder; }

    /// Returns the kinetic constant of this surface reaction. Its exact
    /// units (the way in which it's interpreted by the rest of STEPS) 
    /// depend on the order:
    /// <UL>
    /// <LI>Order 0 -> 1 / s
    /// <LI>Order 1 -> 1 / s
    /// <LI>Order 2 -> 1 / M*s
    /// <LI>Order 3 -> 1 / M*s^2
    /// </UL>
    ///
    double kcst(void) const;
    
    /// Sets the kinetic constant to a new value. Should only be called
    /// before SReacDef::setupFinal.
    ///
    void setKcst(double const & k);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////
    
    //@{
    /// Returns the number of molecules of species idx required in 
    /// the inner volume (_I), outer volume (_O) or surface patch (_S)
    /// to have one occurence of this surface reaction.
    ///
    uint lhs_I(gidxT idx) const;
    uint lhs_S(gidxT idx) const;
    uint lhs_O(gidxT idx) const;
    //@}
    
    //@{
    /// Returns a description of how an occurence of this surface reaction 
    /// depends on some species, defined by its global index idx, to occur. 
    /// See steps/sim/shared/types.hpp for more information on the return 
    /// type. This method is distinct from the SReacDef::req_I, 
    /// SReacDef::req_S and SReacDef::req_O methods.
    ///
    /// This method should only be called after SReacDef::setupFinal has
    /// been called.
    ///
    depT dep_I(gidxT idx) const;
    depT dep_S(gidxT idx) const;
    depT dep_O(gidxT idx) const;
    //@}
    
    //@{
    /// Returns how many molecules of some species, specified by its
    /// global index, are produced after a single occurence of this
    /// surface reaction. '_I' returns this number for the inner volume,
    /// '_S' for the surface patch and '_O' for the outer volume.
    ///
    uint rhs_I(gidxT idx) const;
    uint rhs_S(gidxT idx) const;
    uint rhs_O(gidxT idx) const;
    //@}
    
    //@{
    /// Returns how the amount of a species, specified by its global index,
    /// changes as the result of a single occurence of this surface
    /// reaction on the inside volume (_I), outer volume (_O) or 
    /// surface patch (_S).
    /// 
    /// This method should only be called after SReacDef::setupFinal has 
    /// been called.
    ///
    int upd_I(gidxT idx) const;
    int upd_S(gidxT idx) const;
    int upd_O(gidxT idx) const;
    //@}
    
    //@{
    /// Returns whether the surface reaction rule references a species,
    /// specified by its global index, on the inner volume side (_I),
    /// outer volume (_O) or surface patch (_S). 
    /// 
    /// This method should only be called after SReacDef::setupFinal has
    /// been called.
    ///
    bool req_I(gidxT idx) const;
    bool req_S(gidxT idx) const;
    bool req_O(gidxT idx) const;
    //@}
    
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
    // DATA: GENERAL
    ////////////////////////////////////////////////////////////////////////
    
    StateDef *                  pStateDef;
    
    /// The index of the reaction rule.
    gidxT                       pGIDX;
    /// The name of the reaction rule.
    std::string                 pName;
    
    /// Does the left-hand side of the stoichiometry involve molecules
    /// on the inside or on the outside?
    orientT                     pOrient;
    /// Auxiliary method to compute the order of a reaction.
    void computeOrder(void);
    /// The order of the reaction.
    uint                        pOrder;
    /// Default (MACROscopic) reaction constant.
    double                      pKcst;
    
    /// Tracks whether SReacDef::setupFinal has been called.
    bool                        pFinalSetupDone;
    
    ////////////////////////////////////////////////////////////////////////
    // DATA: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////
    
    //@{
    /// Array vector describing dependencies for inner volume (_I_), 
    /// outer volume (_O_) and surface (_S_) species. Dependencies can 
    /// be stoichiometric or in the rate function (this is not implemented 
    /// yet) -- see 'steps/sim/shared/types.hpp'. The vector must be 
    /// indexed through global indices, i.e. it runs over all species in 
    /// the entire model.
    ///
    /// \todo{Change the name of pSpec_I_DEP into pDEP_I_Spec (similar for
    /// the other two) to promote consistency with reacdef.hpp.}
    depT *                      pSpec_I_DEP;
    depT *                      pSpec_S_DEP;
    depT *                      pSpec_O_DEP;
    //@}
    
    //@{
    /// Vector describing the left hand (reactant) side of the reaction 
    /// stoichiometry, for species in the inner volume (_I_), outer 
    /// volume (_O_) and surface (_S_) species. The vector must be indexed 
    /// through global indices, i.e. it runs over all species in the entire 
    /// model.
    ///
    /// \todo{Rename pSpec_I_LHS to pLHS_I_Spec (similar for the other 
    /// two cases), for consistency with ReacDef.}
    uint *                      pSpec_I_LHS;
    uint *                      pSpec_S_LHS;
    uint *                      pSpec_O_LHS;
    //@}
    
    //@{
    /// An array vector describing the right hand (reaction product) side
    /// of the surface reaction stoichiometry, for species in the inner 
    /// volume (_I_), outer volume (_O_) and surface (_S_) species. The 
    /// vector must be indexed through global indices, i.e. it runs over 
    /// all species in the entire model.
    ///
    /// \todo{Rename pSpec_I_RHS to pRHS_I_Spec, for consistency with
    /// ReacDef (same for the other two in this group).}
    ///
    uint *                      pSpec_I_RHS;
    uint *                      pSpec_S_RHS;
    uint *                      pSpec_O_RHS;
    //@}
    
    //@{
    /// An array describing the update vector (i.e. RHS[] - LHS[]) of 
    /// the surface reaction, for species in the inner volume (_I), 
    /// outer volume (_O_) and patch surface (_S_). The vector must be
    /// indexed through global indices, i.e. it runs over all species in
    /// the entire model.
    /// 
    /// \todo{Rename pSpec_I_UPD to pUPD_I_Spec, for consistency with
    /// ReacDef (same for the other two in this documentation group).}
    int *                       pSpec_I_UPD;
    int *                       pSpec_S_UPD;
    int *                       pSpec_O_UPD;
    //@}
    
    //@{
    /// A vector collecting the global indices of all species that are
    /// updated when this surface reaction rule occurs. 
    gidxTVec                    pSpec_I_UPD_Coll;
    gidxTVec                    pSpec_S_UPD_Coll;
    gidxTVec                    pSpec_O_UPD_Coll;
    //@}

    ////////////////////////////////////////////////////////////////////////
    
};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(sim)
END_NAMESPACE(steps)

#endif
// STEPS_SIM_SHARED_SREACDEF_HPP

// END
