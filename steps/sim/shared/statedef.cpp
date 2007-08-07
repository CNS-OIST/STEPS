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

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/diffdef.hpp>
#include <steps/sim/shared/reacdef.hpp>
#include <steps/sim/shared/specdef.hpp>
#include <steps/sim/shared/statedef.hpp>

USING(std, string);

////////////////////////////////////////////////////////////////////////////////

StateDef::StateDef(void)
: pMode(StateDef::NAIVE_MODE)
, pFinalSetupFinished(false)
, pSpecs()
, pReacs()
, pComps()
{
}

////////////////////////////////////////////////////////////////////////////////

StateDef::~StateDef(void)
{
    for (std::vector<SpecDef*>::iterator i = pSpecs.begin(); 
         i != pSpecs.end(); ++i)
    {
        delete *i;
    }
    for (std::vector<ReacDef*>::iterator i = pReacs.begin(); 
         i != pReacs.end(); ++i)
    {
        delete *i;
    }
    for (std::vector<CompDef*>::iterator i = pComps.begin(); 
         i != pComps.end(); ++i)
    {
        delete *i;
    }
}

////////////////////////////////////////////////////////////////////////////////

void StateDef::setupFinal(void)
{
    // This can only be done once! Check the flag.
    assert(pFinalSetupFinished == false);
    
    // Loop over all simulation variables, kinetics processes and rules
    // that have been declared, and allow them to do some final (internal) 
    // setting up...
    for (std::vector<SpecDef*>::iterator i = pSpecs.begin(); 
         i != pSpecs.end(); ++i)
    {
        (*i)->setupFinal();
    }
    for (std::vector<ReacDef*>::iterator i = pReacs.begin(); 
         i != pReacs.end(); ++i)
    {
        (*i)->setupFinal();
    }
    for (std::vector<DiffDef*>::iterator i = pDiffs.begin();
         i != pDiffs.end(); ++i)
    {
        (*i)->setupFinal();
    }
    
    // Make local indices for species, reactions, diffusion rules, ... in 
    // compartments and then in patches.
    for (std::vector<CompDef*>::iterator i = pComps.begin(); 
         i != pComps.end(); ++i)
    {
        (*i)->setupLocalIndices();
    }
    //for (std::vector<PatchDef*>::iterator i = pPatches.begin();
    //     i != pPatches.end(); ++i)
    //{
    //    (*i)->setupLocalIndices();
    //}
    
    // Resolve update dependencies for various kinetic processes and rules,
    // again first in compartments and then in patches.
    for (std::vector<CompDef*>::iterator i = pComps.begin();
         i != pComps.end(); ++i)
    {
        (*i)->setupDependencies();
    }
    //for (std::vector<PatchDef*>::iterator i = pPatches.begin();
    //     i != pPatches.end(); ++i)
    //{
    //    (*i)->setupDependencies();
    //}
    
    // Set the 'finished' flag.
    pFinalSetupFinished = true;
}

////////////////////////////////////////////////////////////////////////////////

SpecDef * StateDef::createSpecDef(string const & name)
{
    SpecDef * spec = new SpecDef(this, countSpecs(), name);
    assert(spec != 0);
    pSpecs.push_back(spec);
    assert((spec->gidx() + 1) == countSpecs());
    return spec;
}

////////////////////////////////////////////////////////////////////////////////

ReacDef * StateDef::createReacDef(string const & name)
{
    ReacDef * reac = new ReacDef(this, countReacs(), name);
    assert(reac != 0);
    pReacs.push_back(reac);
    assert((reac->gidx() + 1) == countReacs());
    return reac;
}

////////////////////////////////////////////////////////////////////////////////

DiffDef * StateDef::createDiffDef(string const & name)
{
    DiffDef * diff = new DiffDef(this, countDiffs(), name);
    assert(diff != 0);
    pDiffs.push_back(diff);
    assert((diff->gidx() + 1) == countDiffs());
    return diff;
}

////////////////////////////////////////////////////////////////////////////////

CompDef * StateDef::createCompDef(string const & name)
{
    CompDef * comp = new CompDef(this, countComps(), name);
    assert(comp != 0);
    pComps.push_back(comp);
    assert((comp->gidx() + 1) == countComps());
    return comp;
}

////////////////////////////////////////////////////////////////////////////////

// END
