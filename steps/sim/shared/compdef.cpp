////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

// STL headers.
#include <algorithm>
#include <iostream>
#include <iterator>
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
USING(std, vector);

////////////////////////////////////////////////////////////////////////////////

CompUpd::CompUpd(void)
: pLReacs()
, pLDiffs()
{
}

////////////////////////////////////////////////////////////////////////////////

CompUpd::CompUpd(CompUpd const & c)
: pLReacs(c.pLReacs)
, pLDiffs(c.pLDiffs)
{
}

////////////////////////////////////////////////////////////////////////////////

CompUpd::~CompUpd(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void CompUpd::compact(void)
{
    // Compact local reactions.
    std::set<uint> reacset;
    std::copy(pLReacs.begin(), pLReacs.end(), 
        std::inserter(reacset, reacset.begin()));
    pLReacs.clear();
    std::copy(reacset.begin(), reacset.end(), 
        std::inserter(pLReacs, pLReacs.begin()));
    
    // Compact local diffusion rules.
    std::set<uint> diffset;
    std::copy(pLDiffs.begin(), pLDiffs.end(),
        std::inserter(diffset, diffset.begin()));
    pLDiffs.clear();
    std::copy(diffset.begin(), diffset.end(),
        std::inserter(pLDiffs, pLDiffs.begin()));
}

////////////////////////////////////////////////////////////////////////////////

void CompUpd::merge(CompUpd const & upd)
{
    // Merge local reactions -- eliminating duplicates.
    std::set<uint> reacset;
    std::copy(pLReacs.begin(), pLReacs.end(), 
        std::inserter(reacset, reacset.begin()));
    std::copy(upd.pLReacs.begin(), upd.pLReacs.end(), 
        std::inserter(reacset, reacset.begin()));
    pLReacs.clear();
    std::copy(reacset.begin(), reacset.end(), 
        std::inserter(pLReacs, pLReacs.begin()));
    
    // Merge local diffusion rules -- eliminating duplicates.
    std::set<uint> diffset;
    std::copy(pLDiffs.begin(), pLDiffs.end(),
        std::inserter(diffset, diffset.begin()));
    std::copy(upd.pLDiffs.begin(), upd.pLDiffs.end(),
        std::inserter(diffset, diffset.begin()));
    pLDiffs.clear();
    std::copy(diffset.begin(), diffset.end(),
        std::inserter(pLDiffs, pLDiffs.begin()));
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
, pG2LDiff(sdef->countDiffs())
, pL2GDiff(0)
, pDiffUpd(0)
{
    std::fill(pG2LSpec.begin(), pG2LSpec.end(), 0xFFFF);
    std::fill(pG2LReac.begin(), pG2LReac.end(), 0xFFFF);
    std::fill(pG2LDiff.begin(), pG2LDiff.end(), 0xFFFF);
}

////////////////////////////////////////////////////////////////////////////////

CompDef::~CompDef(void)
{
    delete[] pSpecUpd;
    delete[] pReacSpecDeps;
    delete[] pReacSpecUpds;
    delete[] pReacUpd;
    delete[] pDiffUpd;
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::setupLocalIndices(void)
{
    // Make local indices for species that have to be present in voxels
    // belonging to this compartment. These species have been marked earlier
    // with '0' in array pG2LSpec (e.g. by methods such as addReac() or 
    // addDiff()), as opposed in 0xFFFF.
    uint nspecs = std::count(pG2LSpec.begin(), pG2LSpec.end(), 0);
    pL2GSpec.resize(nspecs);
    uint cur = 0;
    for (uint i = 0; i < pG2LSpec.size(); ++i)
    {
        if (pG2LSpec[i] == 0xFFFF) continue;
        assert(pG2LSpec[i] == 0);
        pG2LSpec[i] = cur;
        pL2GSpec[cur++] = i;
    }
    assert(cur == pL2GSpec.size());
    
    // Make local indices for reactions occuring in this compartment. Again,
    // these reactions have been marked with '0' earlier (i.e. in method
    // addReac()).
    uint nreacs = std::count(pG2LReac.begin(), pG2LReac.end(), 0);
    pL2GReac.resize(nreacs);
    cur = 0;
    for (uint i = 0; i < pG2LReac.size(); ++i)
    {
        if (pG2LReac[i] == 0xFFFF) continue;
        assert(pG2LReac[i] == 0);
        pG2LReac[i] = cur;
        pL2GReac[cur++] = i;
    }
    assert(cur == pL2GReac.size());
    
    // Make local indices for diffusion rules occuring in this compartment.
    // Same system as with species and reactions.
    uint ndiffs = std::count(pG2LDiff.begin(), pG2LDiff.end(), 0);
    pL2GDiff.resize(ndiffs);
    cur = 0;
    for (uint i = 0; i < pG2LDiff.size(); ++i)
    {
        if (pG2LDiff[i] == 0xFFFF) continue;
        assert(pG2LDiff[i] == 0);
        pG2LDiff[i] = cur;
        pL2GDiff[cur++] = i;
    }
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::setupDependencies(void)
{
    // Prefetch a couple of useful variables.
    uint nspecs = countSpecs();
    uint nreacs = countReacs();
    uint ndiffs = countDiffs();
    
    // Create a list of CompUpd objects, one for each species variable,
    // that summarize which (local) processes and rules need to be 
    // updated when the species is changed.
    pSpecUpd = new CompUpd[nspecs];
    
    // In the following loop, a number of things are set up:
    //   - for each reaction, its stoichiometric dependencies and update 
    //     vectors. These are stored as C arrays of signed and unsigned
    //     integers, respectively.
    //   - add the reaction to the CompUpd objects of all variables on
    //     which the reaction is dependent.
    pReacSpecDeps = new uint[nspecs * nreacs];
    pReacSpecUpds = new int[nspecs * nreacs];
    for (uint l_ridx = 0; l_ridx < nreacs; ++l_ridx)
    {
        ReacDef * rdef = reac(l_ridx);
        uint idx = l_ridx * nspecs;
        for (uint l_sidx = 0; l_sidx < nspecs; ++l_sidx)
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
    
    // For each diffusion rule, add it to the updates for the ligand
    // species.
    for (uint l_didx = 0; l_didx < ndiffs; ++l_didx)
    {
        DiffDef * ddef = diff(l_didx);
        uint l_sidx = diffG2L(ddef->lig());
        pSpecUpd[l_sidx].pLDiffs.push_back(l_didx);
    }
    
    // For each reaction, build the set of required updates by merging
    // the CompUpd objects of each species that is non-zero in the
    // reaction's stoichiometric update vector.
    pReacUpd = new CompUpd[nreacs];
    for (uint l_ridx = 0; l_ridx < nreacs; ++l_ridx)
    {
        uint idx = l_ridx * nspecs;
        for (uint l_sidx = 0; l_sidx < nspecs; ++l_sidx)
        {
            if (pReacSpecUpds[idx] != 0)
                pReacUpd[l_ridx].merge(pSpecUpd[l_sidx]);
            idx++;
        }
    }
    
    // For each diffusion, build the set of required updates by copying
    // the CompUpd object of the ligand.
    pDiffUpd = new CompUpd[ndiffs];
    for (uint l_didx = 0; l_didx < ndiffs; ++l_didx)
    {
        DiffDef * ddef = diff(l_didx);
        uint l_sidx = diffG2L(ddef->lig());
        pDiffUpd[l_didx].merge(pSpecUpd[l_sidx]);
    }
    
    // Perform compacting step on all CompUpd objects (just to be sure).
    for (uint l_sidx = 0; l_sidx < nspecs; ++l_sidx)
        pSpecUpd[l_sidx].compact();
    for (uint l_ridx = 0; l_ridx < nreacs; ++l_ridx)
        pSpecUpd[l_ridx].compact();
    for (uint l_didx = 0; l_didx < ndiffs; ++l_didx)
        pSpecUpd[l_didx].compact();
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::addSpec(uint gidx)
{
    assert(gidx < pG2LSpec.size());
    pG2LSpec[gidx] = 0;
}

////////////////////////////////////////////////////////////////////////////////

SpecDef * CompDef::spec(uint lidx) const
{
    return statedef()->spec(specL2G(lidx));
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::addReac(uint gidx)
{
    assert(gidx < pG2LReac.size());
    pG2LReac[gidx] = 0;
    
    // Also add all the referenced species in the reaction.
    ReacDef * reac = statedef()->reac(gidx);
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

ReacDef * CompDef::reac(uint lidx) const
{
    return statedef()->reac(reacL2G(lidx));
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::addDiff(uint gidx)
{
    assert(gidx < pG2LDiff.size());
    pG2LDiff[gidx] = 0;
    
    // Also add the ligand to the compartment.
    DiffDef * diff = statedef()->diff(gidx);
    assert(diff != 0);
    addSpec(diff->lig());
}

////////////////////////////////////////////////////////////////////////////////

DiffDef * CompDef::diff(uint lidx) const
{
    return statedef()->diff(diffL2G(lidx));
}

////////////////////////////////////////////////////////////////////////////////

// END
