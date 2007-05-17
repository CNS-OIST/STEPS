////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SIM_SHARED_SPECDEF_HPP
#define STEPS_SIM_SHARED_SPECDEF_HPP 1

// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>

////////////////////////////////////////////////////////////////////////////////

class StateDef;

class SpecDef
{
    
public:

    /// Constructor.
    SpecDef(StateDef * sdef, uint gidx, std::string const & name);
    /// Destructor.
    ~SpecDef(void);

    StateDef * statedef(void) const 
    { return pStateDef; }

    uint gidx(void) const
    { return pGIDX; }

    std::string const & name(void) const
    { return pName; }

    /// Gets called when the definition of all components in the entire state
    /// has finished.
    void finalSetup(void);

private:

    ///
    StateDef *                  pStateDef;

    /// The global (not compartment/patch-specific) index of the species.
    uint                        pGIDX;
    
    /// The name of the species.
    std::string                 pName;

};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_SHARED_SPECDEF_HPP

// END
