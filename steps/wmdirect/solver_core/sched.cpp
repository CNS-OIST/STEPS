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
#include <steps/wmdirect/solver_core/kproc.hpp>
#include <steps/wmdirect/solver_core/sched.hpp>
#include <steps/wmdirect/solver_core/state.hpp>

////////////////////////////////////////////////////////////////////////////////

#define SCHEDULEWIDTH 32
#define MAXLEVELS 10

////////////////////////////////////////////////////////////////////////////////

// Standard library.
USING(std, map);
USING(std, set);
USING(std, vector);

////////////////////////////////////////////////////////////////////////////////

void schedIDXSet_To_Vec(SchedIDXSet const & s, SchedIDXVec & v)
{
    v.resize(s.size());
    std::copy(s.begin(), s.end(), v.begin());
}

////////////////////////////////////////////////////////////////////////////////

/// Unary function that calls the array delete[] operator on pointers. Easy
/// to use with STL/Boost (see steps::tools::DeletePointer).
///
struct DeleteArray
{
    template <typename Type> void operator() (Type * pointer) const
    {
        delete[] pointer;
    }
};

////////////////////////////////////////////////////////////////////////////////

Sched::Sched(void)
: pKProcs()
, pA0(0.0)
, pLevelSizes()
, pLevels()
{
}

////////////////////////////////////////////////////////////////////////////////

Sched::~Sched(void)
{
    std::for_each(pLevels.begin(), pLevels.end(), DeleteArray());
}

////////////////////////////////////////////////////////////////////////////////

void Sched::addKProc(KProc * kproc)
{
    // Check preconditions.
    assert(kproc != 0);
    
    // Add the new entry.
    SchedIDX nidx = pKProcs.size();
    pKProcs.push_back(kproc);
    kproc->setSchedIDX(nidx);
}

////////////////////////////////////////////////////////////////////////////////

void Sched::build(void)
{
    // Setup level.
    uint clevel = 0;
    uint clsize = pKProcs.size();
    
    // Work up.
    while (clsize > 1)
    {
        // Make sure the new size is a multiple of SchedWIDTH.
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

KProc * Sched::getNext(State * state) const
{
    assert(pA0 >= 0.0);
    // Quick check to see whether nothing is there.
    // DEBUG 04-Sep-2007
    if (pA0 == 0.0) return 0;
    
    // Start at top level.
    uint clevel = pLevels.size();
    // And start at the first node of that level.
    uint cur_node = 0;
    
    // Prepare random numbers.
    double rannum[16];
    for (uint i = 0; i < clevel; ++i)
    {
        rannum[i] = state->rng()->getUnfIE();
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

void Sched::reset(void)
{
    // Reset the basic level: compute rates.
    double * oldlevel = pLevels[0];
    uint cur_node = 0;
    vector<KProc*>::iterator kp_end = pKProcs.end();
    for (vector<KProc*>::iterator kp = pKProcs.begin(); kp != kp_end; ++kp)
    {
        oldlevel[cur_node++] = (*kp)->rate();
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

void Sched::recomp(void)
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

void Sched::update(SchedIDXVec const & entries)
{
    // Prefetch zero level.
    double * level0 = pLevels[0];
    // Number of entries.
    assert(entries.size() <= 16);

    // Recompute rates.
    uint indices[16];
    SchedIDXVecCI sidx_end = entries.end();
    uint prev_e = 0xFFFFFFFF;
    uint cur_e = 0;
    for (SchedIDXVecCI sidx = entries.begin(); sidx != sidx_end; ++sidx)
    {
        // Fetch index.
        uint idx = *sidx;
        // Recompute rate, get difference, and store.
        double newrate = pKProcs[idx]->rate();
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
