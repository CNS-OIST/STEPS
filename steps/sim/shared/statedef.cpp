////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
////////////////////////////////////////////////////////////////////////////////

// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/compdef.hpp>
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

void StateDef::finalSetup(void)
{
    assert(pFinalSetupFinished == false);
    for (std::vector<SpecDef*>::iterator i = pSpecs.begin(); 
         i != pSpecs.end(); ++i)
    {
        (*i)->finalSetup();
    }
    for (std::vector<ReacDef*>::iterator i = pReacs.begin(); 
         i != pReacs.end(); ++i)
    {
        (*i)->finalSetup();
    }
    for (std::vector<CompDef*>::iterator i = pComps.begin(); 
         i != pComps.end(); ++i)
    {
        (*i)->finalSetup();
    }
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
