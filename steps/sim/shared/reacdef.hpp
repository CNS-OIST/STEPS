////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SIM_SHARED_REACDEF_HPP
#define STEPS_SIM_SHARED_REACDEF_HPP 1

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
    void finalSetup(void);

    ////////////////////////////////////////////////////////////////////////

    uint order(void) const
    { return pOrder; }

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

    void computeOrder(void);
    
    StateDef *                  pStateDef;
    
    uint                        pGIDX;
    
    std::string                 pName;
    
    uint                        pOrder;
    
    double                      pKf;
    
    std::vector<uint>           pLHS;

    std::vector<uint>           pRHS;

};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_SHARED_REACDEF_HPP

// END
