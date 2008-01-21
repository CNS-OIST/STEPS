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

#ifndef STEPS_SIM_SHARED_REACDEF_HPP
#define STEPS_SIM_SHARED_REACDEF_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/types.hpp>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(sim)

// Forward declarations.
class ReacDef;
class StateDef;

// Auxiliary declarations.
typedef ReacDef *                       ReacDefP;
typedef std::vector<ReacDefP>           ReacDefPVec;
typedef ReacDefPVec::iterator           ReacDefPVecI;
typedef ReacDefPVec::const_iterator     ReacDefPVecCI;

////////////////////////////////////////////////////////////////////////////////

class ReacDef
{

public:

    ReacDef(StateDef * sdef, gidxT idx, std::string const & name);

    ~ReacDef(void);
    
    ////////////////////////////////////////////////////////////////////////
    // REACDEF SETUP
    ////////////////////////////////////////////////////////////////////////
    
    uint incLHS(gidxT idx);
    void setLHS(gidxT idx, uint n);
    
    uint incRHS(gidxT idx);
    void setRHS(gidxT idx, uint n);
    
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

    uint order(void) const
    { return pOrder; }

    double kcst(void) const
    { return pKcst; }
    
    void setKcst(double const & k)
    { pKcst = k; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////
    
    uint lhs(gidxT idx) const;
    depT dep(gidxT idx) const;
    uint rhs(gidxT idx) const;
    int upd(gidxT idx) const;
    bool req(gidxT idx) const;
    
    inline gidxTVecCI bgnUpdColl(void) const
    { return pSpec_UPD_Coll.begin(); }
    inline gidxTVecCI endUpdColl(void) const
    { return pSpec_UPD_Coll.end(); }
    
    ////////////////////////////////////////////////////////////////////////
    
    /// Check whether occurence of a reaction rule <EM>depends</EM> on 
    /// changes to the concentration of some species, specified by its
    /// global index (gidx).
    ///
    /// Currently, this only has to check whether the species gidx is
    /// present in the LHS vector. As we might have concentration-
    /// dependent reaction coefficients in the future (or other fancy
    /// extensions to the basic reaction rule idea), this function will
    /// have to check for more things.
    ///
    //bool dependsOnSpec(gidxT idx) const;
    
    /// Check whether the occurence of the reaction rule <EM>affects</EM>
    /// the concentration of some species, specified by its global index
    /// (gidx).
    ///
    /// Currently, this only has to check whether the stoichiometry vectors,
    /// to see whether rhs[gidx] - lhs[gidx] is zero or not. (If it's zero,
    /// it's not affected either because it's not used, or because the
    /// <EM>net effect</EM> of the reaction cancels out changes. If it's
    /// not zero, this means the species <EM>is</EM> in fact affected.)
    ///
    //bool affectsSpec(gidxT idx) const;
    
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
    
    depT *                      pSpec_DEP;
    uint *                      pSpec_LHS;
    uint *                      pSpec_RHS;
    int *                       pSpec_UPD;
    gidxTVec                    pSpec_UPD_Coll;

    ////////////////////////////////////////////////////////////////////////
    
};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(sim)
END_NAMESPACE(steps)

#endif
// STEPS_SIM_SHARED_REACDEF_HPP

// END
