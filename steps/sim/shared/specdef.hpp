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

#ifndef STEPS_SIM_SHARED_SPECDEF_HPP
#define STEPS_SIM_SHARED_SPECDEF_HPP 1

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
    ///
    /// Currently, this method doesn't do anything.
    ///
    void setupFinal(void);

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
