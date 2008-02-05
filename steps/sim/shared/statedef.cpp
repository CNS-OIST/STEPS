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
#include <cassert>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/diffdef.hpp>
#include <steps/sim/shared/patchdef.hpp>
#include <steps/sim/shared/reacdef.hpp>
#include <steps/sim/shared/specdef.hpp>
#include <steps/sim/shared/sreacdef.hpp>
#include <steps/sim/shared/statedef.hpp>
#include <steps/sim/shared/types.hpp>

USING(std, string);
USING(std, vector);
USING_NAMESPACE(steps::sim);

////////////////////////////////////////////////////////////////////////////////

StateDef::StateDef(void)
: pMode(StateDef::NAIVE_MODE)
, pFinalSetupFinished(false)
, pSpecs()
, pComps()
, pReacs()
, pDiffs()
, pPatches()
, pSReacs()
{
}

////////////////////////////////////////////////////////////////////////////////

StateDef::~StateDef(void)
{
    for (SpecDefPVecI i = pSpecs.begin(); i != pSpecs.end(); ++i) 
        delete *i;
    for (CompDefPVecI i = pComps.begin(); i != pComps.end(); ++i) 
        delete *i;
    for (ReacDefPVecI i = pReacs.begin(); i != pReacs.end(); ++i) 
        delete *i;
    for (DiffDefPVecI i = pDiffs.begin(); i != pDiffs.end(); ++i) 
        delete *i;
    for (PatchDefPVecI i = pPatches.begin(); i != pPatches.end(); ++i) 
        delete *i;
    for (SReacDefPVecI i = pSReacs.begin(); i != pSReacs.end(); ++i) 
        delete *i;
}

////////////////////////////////////////////////////////////////////////////////

void StateDef::setupFinal(void)
{
    // This is an important routine to understand all the *Def code in 
    // directory steps/sim/shared.
    
    // This can only be done once! Check the flag.
    assert(pFinalSetupFinished == false);
    
    // Loop over all simulation variables, kinetics processes and rules
    // that have been declared, and allow them to do some final (internal) 
    // setting up...
    for (SpecDefPVecI i = pSpecs.begin(); i != pSpecs.end(); ++i)
        (*i)->setupFinal();
    for (ReacDefPVecI i = pReacs.begin(); i != pReacs.end(); ++i)
        (*i)->setupFinal();
    for (DiffDefPVecI i = pDiffs.begin(); i != pDiffs.end(); ++i)
        (*i)->setupFinal();
    for (SReacDefPVecI i = pSReacs.begin(); i != pSReacs.end(); ++i)
        (*i)->setupFinal();
    
    // Make sure all variables required by the kinetic processes and rules
    // are present in each compartment and patch.
    for (CompDefPVecI i = pComps.begin(); i != pComps.end(); ++i)
        (*i)->addReferences();
    for (PatchDefPVecI i = pPatches.begin(); i != pPatches.end(); ++i)
        (*i)->addReferences();
    
    // Make local indices for species, reactions, diffusion rules, ... in 
    // compartments and then in patches.
    for (CompDefPVecI i = pComps.begin(); i != pComps.end(); ++i)
        (*i)->setupLocalIndices();
    for (PatchDefPVecI i = pPatches.begin(); i != pPatches.end(); ++i)
        (*i)->setupLocalIndices();
    
    // Resolve update dependencies for various kinetic processes and rules,
    // again first in compartments and then in patches.
    for (CompDefPVecI i = pComps.begin(); i != pComps.end(); ++i)
        (*i)->setupDependencies();
    for (PatchDefPVecI i = pPatches.begin(); i != pPatches.end(); ++i)
        (*i)->setupDependencies();
    
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

CompDef * StateDef::createCompDef(string const & name)
{
    CompDef * comp = new CompDef(this, countComps(), name);
    assert(comp != 0);
    pComps.push_back(comp);
    assert((comp->gidx() + 1) == countComps());
    return comp;
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

PatchDef * StateDef::createPatchDef
(
    string const & name, 
    CompDef * inner, 
    CompDef * outer
)
{
    PatchDef * patch = new PatchDef(this, countPatches(), name, inner, outer);
    assert(patch != 0);
    pPatches.push_back(patch);
    assert((patch->gidx() + 1) == countPatches());
    return patch;
}

////////////////////////////////////////////////////////////////////////////////

SReacDef * StateDef::createSReacDef(string const & name, bool inside)
{
    SReacDef::orientT o = 
        (inside == true ? SReacDef::INSIDE : SReacDef::OUTSIDE);
    SReacDef * sreac = new SReacDef(this, countSReacs(), name, o);
    assert(sreac != 0);
    pSReacs.push_back(sreac);
    assert((sreac->gidx() + 1) == countSReacs());
    return sreac;
}

////////////////////////////////////////////////////////////////////////////////

// END
