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

////////////////////////////////////////////////////////////////////////////////

class StateDef;

class ReacDef
{

public:

    ReacDef(StateDef * sdef, uint gidx, std::string const & name);

    StateDef * statedef(void) const
    { return pStateDef; }
    
    uint gidx(void) const
    { return pGIDX; }
    
    std::string const & name(void) const
    { return pName; }

    ////////////////////////////////////////////////////////////////////////
    
    /// Gets called when the definition of all components in the entire state
    /// has finished.
    ///
    /// Currently, this method only pre-computes the order of the reaction.
    ///
    void setupFinal(void);

    ////////////////////////////////////////////////////////////////////////

    uint order(void) const
    { return pOrder; }

    double kcst(void) const
    { return pKcst; }
    
    void setKcst(double const & k)
    { pKcst = k; }

    ////////////////////////////////////////////////////////////////////////
    
    std::vector<uint>::const_iterator beginLHS(void) const
    { return pLHS.begin(); }
    
    std::vector<uint>::const_iterator endLHS(void) const
    { return pLHS.end(); }
    
    uint lhs(uint gidx) const;
    
    uint incLHS(uint gidx);
    
    void setLHS(uint gidx, uint n);
    
    ////////////////////////////////////////////////////////////////////////
    
    std::vector<uint>::const_iterator beginRHS(void) const
    { return pRHS.begin(); }
    
    std::vector<uint>::const_iterator endRHS(void) const
    { return pRHS.end(); }
    
    uint rhs(uint gidx) const;
    
    uint incRHS(uint gidx);
    
    void setRHS(uint gidx, uint n);

private:

    /// Auxiliary method to compute the order of a reaction.
    ///
    void computeOrder(void);
    
    StateDef *                  pStateDef;
    
    ///
    
    uint                        pGIDX;
    
    std::string                 pName;
    
    /// The order of the reaction.
    ///
    uint                        pOrder;
    
    /// Default (MACROscopic) reaction constant.
    ///
    double                      pKcst;
    
    std::vector<uint>           pLHS;

    std::vector<uint>           pRHS;

};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_SHARED_REACDEF_HPP

// END
