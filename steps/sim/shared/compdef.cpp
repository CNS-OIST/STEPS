////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

// STL headers.
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/reacdef.hpp>
#include <steps/sim/shared/specdef.hpp>
#include <steps/sim/shared/statedef.hpp>

USING(std, string);
USING(std, vector);

////////////////////////////////////////////////////////////////////////////////

CompUpd::CompUpd(void)
: pLReacs()
{
}

////////////////////////////////////////////////////////////////////////////////

CompUpd::~CompUpd(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void CompUpd::merge(CompUpd const & upd)
{
    // Merge local reactions.
    std::set<uint> reacset;
    std::copy(pLReacs.begin(), pLReacs.end(), 
        std::inserter(reacset, reacset.begin()));
    std::copy(upd.pLReacs.begin(), upd.pLReacs.end(), 
        std::inserter(reacset, reacset.begin()));
    pLReacs.clear();
    std::copy(reacset.begin(), reacset.end(), 
        std::inserter(pLReacs, pLReacs.begin()));
}

////////////////////////////////////////////////////////////////////////////////

CompDef::CompDef(StateDef * sdef, uint gidx, string const & name)
: pStateDef(sdef)
, pGIDX(gidx)
, pName(name)
, pG2LSpec(sdef->countSpecs())
, pL2GSpec(0)
, pSpecUpd(0)
, pG2LReac(sdef->countReacs())
, pL2GReac(0)
, pReacSpecDeps(0)
, pReacSpecUpds(0)
, pReacUpd(0)
{
    std::fill(pG2LSpec.begin(), pG2LSpec.end(), 0xFFFF);
    std::fill(pG2LReac.begin(), pG2LReac.end(), 0xFFFF);
}

////////////////////////////////////////////////////////////////////////////////

CompDef::~CompDef(void)
{
    delete[] pSpecUpd;
    delete[] pReacSpecDeps;
    delete[] pReacSpecUpds;
    delete[] pReacUpd;
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::finalSetup(void)
{
    // Make indices for species.
    uint nspecies = std::count(pG2LSpec.begin(), pG2LSpec.end(), 0);
    pL2GSpec.resize(nspecies);
    uint cur = 0;
    for (uint i = 0; i < pG2LSpec.size(); ++i)
    {
        if (pG2LSpec[i] == 0xFFFF) continue;
        assert(pG2LSpec[i] == 0);
        pG2LSpec[i] = cur;
        pL2GSpec[cur++] = i;
    }
    assert(cur == pL2GSpec.size());
    
    // Make indices for reactions.
    uint nreactions = std::count(pG2LReac.begin(), pG2LReac.end(), 0);
    pL2GReac.resize(nreactions);
    cur = 0;
    for (uint i = 0; i < pG2LReac.size(); ++i)
    {
        if (pG2LReac[i] == 0xFFFF) continue;
        assert(pG2LReac[i] == 0);
        pG2LReac[i] = cur;
        pL2GReac[cur++] = i;
    }
    assert(cur == pL2GReac.size());
    
    // For each reaction, build 
    pSpecUpd = new CompUpd[nspecies];
    pReacSpecDeps = new uint[nspecies * nreactions];
    pReacSpecUpds = new int[nspecies * nreactions];
    for (uint l_ridx = 0; l_ridx < nreactions; ++l_ridx)
    {
        ReacDef * rdef = reac(l_ridx);
        uint idx = l_ridx * nspecies;
        for (uint l_sidx = 0; l_sidx < nspecies; ++l_sidx)
        {
            uint g_sidx = specL2G(l_sidx);
            uint lhs = rdef->lhs(g_sidx);
            uint rhs = rdef->rhs(g_sidx);
            pReacSpecDeps[idx] = lhs;
            pReacSpecUpds[idx] = static_cast<int>(rhs) - static_cast<int>(lhs);
            if (lhs > 0) pSpecUpd[l_sidx].pLReacs.push_back(l_ridx);
            idx++;
        }
    }
    
    //
    pReacUpd = new CompUpd[nreactions];
    for (uint l_ridx = 0; l_ridx < nreactions; ++l_ridx)
    {
        uint idx = l_ridx * nspecies;
        for (uint l_sidx = 0; l_sidx < nspecies; ++l_sidx)
        {
            if (pReacSpecUpds[idx] != 0)
                pReacUpd[l_ridx].merge(pSpecUpd[l_sidx]);
            idx++;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

SpecDef * CompDef::spec(uint lidx) const
{
    return statedef()->spec(specL2G(lidx));
}

////////////////////////////////////////////////////////////////////////////////

ReacDef * CompDef::reac(uint lidx) const
{
    return statedef()->reac(reacL2G(lidx));
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::addSpec(uint gidx)
{
    assert(gidx < pG2LSpec.size());
    pG2LSpec[gidx] = 0;
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::addReac(uint gidx)
{
    assert(gidx < pG2LReac.size());
    pG2LReac[gidx] = 0;
    
    // Also add all the referenced species in the reaction.
    ReacDef * reac = pStateDef->reac(gidx);
    assert(reac != 0);
    for (std::vector<uint>::const_iterator i = reac->beginLHS(); 
         i != reac->endLHS(); ++i)
    {
        if (*i > 0) addSpec(*i);
    }
    for (std::vector<uint>::const_iterator i = reac->beginRHS();
         i != reac->endRHS(); ++i)
    {
        if (*i > 0) addSpec(*i);
    }
}

////////////////////////////////////////////////////////////////////////////////

// END
