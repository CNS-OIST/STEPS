////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

// STL headers.
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/constants.hpp>
#include <steps/rng/rng.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/reacdef.hpp>
#include <steps/sim/shared/specdef.hpp>
#include <steps/sim/shared/statedef.hpp>
#include <steps/sim/wmdirect/state.hpp>

NAMESPACE_ALIAS(steps::math, smath);
NAMESPACE_ALIAS(steps::rng, srng);

////////////////////////////////////////////////////////////////////////////////

State::State(void)
: fTime(0.0)
, fNSteps(0)
, fPoolCount(0)
, fPoolFlags(0)
, fReacKcsts(0)
, fReacCcsts(0)
, fReacProps(0)
, fReacFlags(0)
, fCompVols(0)
, pStateDef(0)
{
    pStateDef = new StateDef();
}

////////////////////////////////////////////////////////////////////////////////

State::~State(void)
{
    uint ncomps = def()->countComps();
    for (uint i = 0; i < ncomps; ++i)
    {
        delete[] fPoolCount[i];
        delete[] fPoolFlags[i];
        delete[] fReacKcsts[i];
        delete[] fReacCcsts[i];
        delete[] fReacProps[i];
        delete[] fReacFlags[i];
    }
    delete[] fPoolCount;
    delete[] fPoolFlags;
    delete[] fReacKcsts;
    delete[] fReacCcsts;
    delete[] fReacProps;
    delete[] fReacFlags;
    delete[] fCompVols;
       
    // Finally, delete the state definition.
    delete pStateDef;
}

////////////////////////////////////////////////////////////////////////////////

void State::setup(void)
{
    // First create pools & reactions.
    uint ncomps = pStateDef->countComps();
    fPoolCount = new uint*[ncomps];
    fPoolFlags = new uint*[ncomps];
    fReacKcsts = new double*[ncomps];
    fReacCcsts = new double*[ncomps];
    fReacProps = new double*[ncomps];
    fReacFlags = new uint*[ncomps];
    for (uint i = 0; i < ncomps; ++i)
    {
        CompDef * cdef = pStateDef->comp(i);
        
        // Create pools.
        uint nspecs = cdef->countSpecs();
        fPoolCount[i] = new uint[nspecs];
        fPoolFlags[i] = new uint[nspecs];
        
        // Create reaction stuff.
        uint nreacs = cdef->countReacs();
        fReacKcsts[i] = new double[nreacs];
        fReacCcsts[i] = new double[nreacs];
        fReacProps[i] = new double[nreacs];
        fReacFlags[i] = new uint[nreacs];
    }
    
    // Setup the compartment volumes.
    fCompVols = new double[ncomps];
    for (uint i = 0; i < ncomps; ++i)
    {
        fCompVols[i] = 0.0;
    }
}

////////////////////////////////////////////////////////////////////////////////

void State::reset(void)
{
    fTime = 0.0;
    fNSteps = 0;

    uint ncomps = pStateDef->countComps();
    for (uint i = 0; i < ncomps; ++i)
    {
        CompDef * cdef = pStateDef->comp(i);
        
        // Set pools to 0 molecules.
        uint nspecs = cdef->countSpecs();
        for (uint j = 0; j < nspecs; ++j) fPoolCount[i][j] = 0;
        // Reset flags.
        for (uint j = 0; j < nspecs; ++j) fPoolFlags[i][j] = PoolFlagDefault;
        
        // Set propensities to 0.0.
        uint nreacs = cdef->countReacs();
        for (uint j = 0; j < nreacs; ++j) fReacKcsts[i][j] = 0.0;
        for (uint j = 0; j < nreacs; ++j) fReacCcsts[i][j] = 0.0;
        for (uint j = 0; j < nreacs; ++j) fReacProps[i][j] = 0.0;
        // Reset flags.
        for (uint j = 0; j < nreacs; ++j) fReacFlags[i][j] = ReacFlagDefault;
    }
}

////////////////////////////////////////////////////////////////////////////////

void State::computeCcst(uint cidx, uint l_ridx)
{
//std::cout << "COMPUTE CCST\n";
    // Scaling base.
    double vscale = 1.0e3 * fCompVols[cidx] * smath::AVOGADRO;
//std::cout << "vscale = " << vscale << std::endl;
    // Fetch order of reaction.
    uint g_ridx = def()->comp(cidx)->reacL2G(l_ridx);
    int o1 = static_cast<int>(def()->reac(g_ridx)->order()) - 1;
//std::cout << "o1 = " << o1 << std::endl;
    if (o1 < 0) o1 = 0;
    vscale = pow(vscale, static_cast<double>(-o1));
//std::cout << "vscale' = " << vscale << std::endl;
    // Set the order-dependent scaled reaction constant.
//std::cout << "Kf = " << fReacKcsts[cidx][l_ridx] << std::endl;
    fReacCcsts[cidx][l_ridx] = fReacKcsts[cidx][l_ridx] * vscale;
//std::cout << "Cf = " << fReacCcsts[cidx][l_ridx] << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

void State::computeReacProp(uint cidx, uint l_ridx)
{
    // Check whether the reaction is active.
    if ((fReacFlags[cidx][l_ridx] & State::ACTIVE_REACFLAG) == 0)
    {
        fReacProps[cidx][l_ridx] = 0.0;
        return;
    }
    
    // Prefetch some variables.
    CompDef * cdef = def()->comp(cidx);
    uint nspecs = cdef->countSpecs();
    uint * deps = cdef->reacSpecDeps(l_ridx);
    uint * cnts = fPoolCount[cidx];
    
    // Compute combinatorial part.
    double h_mu = 1.0;
    for (uint pool = 0; pool < nspecs; ++pool)
    {
        uint dep = deps[pool];
        if (dep == 0) continue;
        uint cnt = cnts[pool];
        if (dep > cnt) 
        {
            h_mu = 0.0;
            break;
        }
        switch (dep)
        {
            case 4:
            {
                h_mu *= static_cast<double>(cnt - 3);
            }
            case 3:
            {
                h_mu *= static_cast<double>(cnt - 2);
            }
            case 2:
            {
                h_mu *= static_cast<double>(cnt - 1);
            }
            case 1:
            {
                h_mu *= static_cast<double>(cnt);
                break;
            }
            default:
            {
                assert(0);
                return;
            }
        }
    }

//std::cout << "h_mu = " << h_mu << std::endl;
//std::cout << "c_mu = " << fReacCcsts[cidx][l_ridx] << std::endl;

    // Multiply with scaled reaction constant.
    fReacProps[cidx][l_ridx] = h_mu * fReacCcsts[cidx][l_ridx];
//std::cout << "a_mu = " << fReacProps[cidx][l_ridx] << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

double State::computeZeroProp(void) const
{
    double a0 = 0.0;
    
    // Loop over all compartments.
    uint ncomps = def()->countComps();
    for (uint c = 0; c < ncomps; ++c)
    {
        uint nreacs = def()->comp(c)->countReacs();
        for (uint r = 0; r < nreacs; ++r)
        {
            a0 += fReacProps[c][r];
        }
    }
    return a0;
}

////////////////////////////////////////////////////////////////////////////////

// END
