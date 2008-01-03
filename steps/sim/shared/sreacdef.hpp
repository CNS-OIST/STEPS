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
    /// Currently, this method only pre-computes the order of the reaction.
    ///
    void setupFinal(void);
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: GENERAL
    ////////////////////////////////////////////////////////////////////////
    
    StateDef * statedef(void) const
    { return pStateDef; }
    
    gidxT gidx(void) const
    { return pGIDX; }
    
    std::string const & name(void) const
    { return pName; }

    ////////////////////////////////////////////////////////////////////////

    /// Returns true if the left hand side of the reaction stoichiometry
    /// involves reactants on the surface and on the inside volume.
    ///
    bool inside(void) const
    { return (pOrient == INSIDE); }
    /// Returns true if the left hand side of the reaction stoichiometry
    /// involves reactants on the surface and on the outside volume.
    ///
    bool outside(void) const
    { return (pOrient == OUTSIDE); }
    
    uint order(void) const
    { return pOrder; }

    double kcst(void) const
    { return pKcst; }
    
    void setKcst(double const & k)
    { pKcst = k; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////
    
    uint lhs_I(gidxT idx) const;
    uint lhs_S(gidxT idx) const;
    uint lhs_O(gidxT idx) const;
    
    depT dependsOnSpec_I(gidxT idx) const;
    depT dependsOnSpec_S(gidxT idx) const;
    depT dependsOnSpec_O(gidxT idx) const;
    
    uint rhs_I(gidxT idx) const;
    uint rhs_S(gidxT idx) const;
    uint rhs_O(gidxT idx) const;
    
    int upd_I(gidxT idx) const;
    int upd_S(gidxT idx) const;
    int upd_O(gidxT idx) const;
    
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
    
    ////////////////////////////////////////////////////////////////////////
    // DATA: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////
    
    depT *                      pSpec_I_DEP;
    depT *                      pSpec_S_DEP;
    depT *                      pSpec_O_DEP;
    uint *                      pSpec_I_LHS;
    uint *                      pSpec_S_LHS;
    uint *                      pSpec_O_LHS;
    uint *                      pSpec_I_RHS;
    uint *                      pSpec_S_RHS;
    uint *                      pSpec_O_RHS;
    int *                       pSpec_I_UPD;
    int *                       pSpec_S_UPD;
    int *                       pSpec_O_UPD;
    gidxTVec                    pSpec_I_UPD_Coll;
    gidxTVec                    pSpec_S_UPD_Coll;
    gidxTVec                    pSpec_O_UPD_Coll;

    ////////////////////////////////////////////////////////////////////////
    
};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(sim)
END_NAMESPACE(steps)

#endif
// STEPS_SIM_SHARED_SREACDEF_HPP

// END
