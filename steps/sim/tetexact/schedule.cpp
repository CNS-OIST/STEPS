////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/constants.hpp>
#include <steps/rng/rng.hpp>
#include <steps/sim/kproc.hpp>
#include <steps/sim/schedule.hpp>
#include <steps/sim/simenv.hpp>
#include <steps/tools/cplusplus.hpp>

////////////////////////////////////////////////////////////////////////////////

#define SCHEDULEWIDTH 16
#define MAXLEVELS 10

////////////////////////////////////////////////////////////////////////////////

// Standard library.
USING(std, map);
USING(std, set);
USING(std, vector);

// STEPS library.
NAMESPACE_ALIAS(steps::math, smath);
NAMESPACE_ALIAS(steps::rng, srng);
USING(srng, RNG);
USING(srng, RNGPtr);
NAMESPACE_ALIAS(steps::sim, ssim);
USING(ssim, KProc);
USING(ssim, KProcPtr);
USING(ssim, KProcPtrVector);
USING(ssim, KProcPtrVectorIt);
USING(ssim, KProcPtrVectorCtIt);
USING(ssim, Schedule);
USING(ssim, ScheduleIDX);
USING(ssim, ScheduleIDXVector);
USING(ssim, ScheduleIDXVectorIt);
USING(ssim, SimEnv);
NAMESPACE_ALIAS(steps::tools, stools);

////////////////////////////////////////////////////////////////////////////////

Schedule::Schedule(void)
: pKProcs()
, pA0(0.0)
, pLevelSizes()
, pLevels()
{
}

////////////////////////////////////////////////////////////////////////////////

Schedule::~Schedule(void)
{
    std::for_each(pLevels.begin(), pLevels.end(), stools::DeleteArray());
}

////////////////////////////////////////////////////////////////////////////////

void Schedule::addKProc(KProc * kproc)
{
    // Check preconditions.
    assert(kproc != 0);
    
    // Add the new entry.
    ScheduleIDX nidx = pKProcs.size();
    pKProcs.push_back(kproc);
    kproc->setScheduleIDX(nidx);
}

////////////////////////////////////////////////////////////////////////////////

void Schedule::build(void)
{
    // Setup level.
    uint clevel = 0;
    uint clsize = pKProcs.size();
    
    // Work up.
    while (clsize > 1)
    {
        // Make sure the new size is a multiple of SCHEDULEWIDTH.
        uint extra = clsize % SCHEDULEWIDTH;
        if (extra != 0) clsize += SCHEDULEWIDTH - extra;
        
        // Create the level and add it.
        double * level = new double[clsize];
        std::fill_n(level, clsize, 0.0);
        pLevelSizes.push_back(clsize);
        pLevels.push_back(level);
    
        // Prepare for next level.
        clevel++;
        clsize = clsize / SCHEDULEWIDTH;
    }
    
    // Set top level.
    pA0 = 0.0;
}

////////////////////////////////////////////////////////////////////////////////

//
// THIS VERSION IS NUMERICALLY MORE ROBUST!
//

KProc * Schedule::getNext(SimEnv & simenv) const
{
    // Start at top level.
    uint clevel = pLevels.size();
    // And start at the first node of that level.
    uint cur_node = 0;
    
    // Prepare random numbers.
    double rannum[16];
    for (uint i = 0; i < clevel; ++i)
    {
        rannum[i] = simenv.fRNG->getUnfIE();
    }
    
    // Run until top level.
    double a0 = pA0;
    double * level = 0;
    while (clevel != 0)
    {
        // Decrease the current level.
        clevel--;
        // and start looking in the right place.
        cur_node *= SCHEDULEWIDTH;
        uint max_node = cur_node + SCHEDULEWIDTH;
        
        // Fetch the level.
        level = pLevels[clevel];
        
        // Compute local selector.
        double selector = rannum[clevel] * a0;
        
        // Compare.
        double accum = 0.0;
        double old = 0.0;
        double curval = 0.0;
        for (uint i = 0; i < SCHEDULEWIDTH; ++i)
        {
            curval = level[cur_node];
            if (selector < curval + accum) break;
            accum += curval;
            old = accum;
            cur_node++;
        }
        
        // Checks.
        assert(cur_node < max_node);
        assert(curval > 0.0);
        a0 = curval;
    }
    
    // Check.
    assert(cur_node < pKProcs.size());
    return pKProcs[cur_node];
}

////////////////////////////////////////////////////////////////////////////////

void Schedule::reset(void)
{
    // Reset the basic level: compute rates.
    double * oldlevel = pLevels[0];
    uint cur_node = 0;
    vector<KProc*>::iterator kp_end = pKProcs.end();
    for (vector<KProc*>::iterator kp = pKProcs.begin(); kp != kp_end; ++kp)
    {
        oldlevel[cur_node++] = (*kp)->computeRate();
    }
    
    // Work up.
    for (uint cur_level = 1; cur_level < pLevels.size(); ++cur_level)
    {
        // Compute the number of nodes to reset on this level.
        uint numnodes = pLevelSizes[cur_level - 1] / SCHEDULEWIDTH;
        
        // Fetch a pointer to this level.
        double * level = pLevels[cur_level];
        
        // Recompute them.
        uint child_node = 0;
        for (cur_node = 0; cur_node < numnodes; ++cur_node)
        {
            double val = 0.0;
            for (uint i = 0; i < SCHEDULEWIDTH; ++i)
            {
                val += oldlevel[child_node++];
            }
            level[cur_node] = val;
        }
        
        // Copy the level.
        oldlevel = level;
    }
    
    // Compute zero propensity.
    pA0 = 0.0;
    for (uint i = 0; i < SCHEDULEWIDTH; ++i)
    {
        pA0 += oldlevel[i];
    }
}

////////////////////////////////////////////////////////////////////////////////

void Schedule::recomp(void)
{
    // Setup.
    double * oldlevel = pLevels[0];
    uint cur_node = 0;
    
    // Work up.
    for (uint cur_level = 1; cur_level < pLevels.size(); ++cur_level)
    {
        // Compute the number of nodes to reset on this level.
        uint numnodes = pLevelSizes[cur_level - 1] / SCHEDULEWIDTH;
        
        // Fetch a pointer to this level.
        double * level = pLevels[cur_level];
        
        // Recompute them.
        uint child_node = 0;
        for (cur_node = 0; cur_node < numnodes; ++cur_node)
        {
            double val = 0.0;
            for (uint i = 0; i < SCHEDULEWIDTH; ++i)
            {
                val += oldlevel[child_node++];
            }
            level[cur_node] = val;
        }
        
        // Copy the level.
        oldlevel = level;
    }
    
    // Compute zero propensity.
    pA0 = 0.0;
    for (uint i = 0; i < SCHEDULEWIDTH; ++i)
    {
        pA0 += oldlevel[i];
    }
}

////////////////////////////////////////////////////////////////////////////////

void Schedule::update(ScheduleIDXVector const & entries)
{
    // Prefetch zero level.
    double * level0 = pLevels[0];
    // Number of entries.
    assert(entries.size() <= 16);

    // Recompute rates.
    uint indices[16];
    ScheduleIDXVectorCtIt sidx_end = entries.end();
    uint prev_e = 0xFFFFFFFF;
    uint cur_e = 0;
    for (ScheduleIDXVectorCtIt sidx = entries.begin(); sidx != sidx_end; ++sidx)
    {
        // Fetch index.
        uint idx = *sidx;
        // Recompute rate, get difference, and store.
        double newrate = pKProcs[idx]->computeRate();
        level0[idx] = newrate;
        
        // Store and collapse if possible.
        idx /= SCHEDULEWIDTH;
        if (prev_e == 0xFFFFFFFF)
        {
            prev_e = 0;
            indices[cur_e++] = idx;
        }
        else if (indices[prev_e] != idx)
        {
            prev_e = cur_e;
            indices[cur_e++] = idx;
        }
    }
    uint nentries = cur_e;

    // Update upper levels.
    uint nlevels = pLevels.size();
    double * prevlevel = pLevels[0];
    for (uint l = 1; l < nlevels; ++l)
    {
        // Update the first entry.
        cur_e = 0;
        prev_e = 0xFFFFFFFF;
        
        // Fetch a pointer to the current level.
        double * currlevel = pLevels[l];
           
        // Recompute the entries.
        for (uint e = 0; e < nentries; ++e)
        {
            // Fetch index.
            uint idx = indices[e];
            
            // Recompute.
            double val = 0.0;
            uint idx2 = idx * SCHEDULEWIDTH;
            for (uint i = 0; i < SCHEDULEWIDTH; ++i)
            {
                val += prevlevel[idx2++];
            }
            currlevel[idx] = val;
            
            // Store and collapse if possible.
            idx /= SCHEDULEWIDTH;
            if (prev_e == 0xFFFFFFFF)
            {
                prev_e = 0;
                indices[cur_e++] = idx;
            }
            else if (indices[prev_e] != idx)
            {
                prev_e = cur_e;
                indices[cur_e++] = idx;
            }
        }

        // Update the pointer to the previous level.
        prevlevel = currlevel;

        // cur_e now is the new number of entries to handle.
        nentries = cur_e;
    }

    // Update zero propensity.
    double * toplevel = pLevels[pLevels.size() - 1];
    pA0 = 0.0;
    for (uint i = 0; i < SCHEDULEWIDTH; ++i)
    {
        pA0 += toplevel[i];
    }
}

////////////////////////////////////////////////////////////////////////////////

// END
