/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */


// Standard library & STL headers.
// #include <vector>
#include <algorithm>

// STEPS headers.
#include "comp.hpp"
#include "reac.hpp"
#include "wmrssa.hpp"

// logging
#include <easylogging++.h>
#include "util/error.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace swmrssa = steps::wmrssa;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

swmrssa::Comp::Comp(steps::solver::Compdef * compdef)
: pCompdef(compdef)
, pKProcs()
, pIPatches()
, pOPatches()
{
    assert (pCompdef != nullptr);
    uint nspecs = compdef->countSpecs();
    pPoolLB = new double[nspecs]();
    pPoolUB = new double[nspecs]();
}

////////////////////////////////////////////////////////////////////////////////

swmrssa::Comp::~Comp()
{
    for (auto const& k :pKProcs) {
        delete k;
    }
    delete[] pPoolLB;
    delete[] pPoolUB;
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Comp::checkpoint(std::fstream & cp_file)
{
    for (auto const& k : pKProcs) {
        k->checkpoint(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Comp::restore(std::fstream & cp_file)
{
    for (auto const& k : pKProcs) {
        k->restore(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Comp::reset()
{
    for (auto const& k: pKProcs) {
        k->reset();
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Comp::setupKProcs(swmrssa::Wmrssa * wmd)
{
    // Create reaction kproc's.
    uint nreacs = def()->countReacs();
    pKProcs.resize(nreacs);
    for (uint i = 0; i < nreacs; ++i)
    {
        ssolver::Reacdef * rdef = def()->reacdef(i);
        auto * r = new swmrssa::Reac(rdef, this);
        pKProcs[i] = r;
        wmd->addKProc(r);
    }
}

////////////////////////////////////////////////////////////////////////////////

steps::wmrssa::KProc * swmrssa::Comp::reac(uint lridx) const
{
    assert (lridx < pKProcs.size());
    return pKProcs[lridx];
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Comp::setupDeps()
{
    for (auto const& k: pKProcs) {
        k->setupDeps();
    }
}

////////////////////////////////////////////////////////////////////////////////


void swmrssa::Comp::addIPatch(swmrssa::Patch * p)
{
    AssertLog(std::find(pIPatches.begin(), pIPatches.end(), p) == pIPatches.end());
    pIPatches.push_back(p);
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Comp::addOPatch(swmrssa::Patch * p)
{
    AssertLog(std::find(pOPatches.begin(), pOPatches.end(), p) == pOPatches.end());
    pOPatches.push_back(p);
}
////////////////////////////////////////////////////////////////////////////////

void swmrssa::Comp::setBounds(uint i, int nc)
{
    const double delta = .05;
    if (nc > 3/delta)
    {
        pPoolLB[i] = nc*(1 - delta);
        pPoolUB[i] = nc*(1 + delta);
    }
    else if (nc > 3)
    {
        pPoolLB[i] = nc - 3;
        pPoolUB[i] = nc + 3;
    }
    else if (nc > 0)
    {
        pPoolLB[i] = 1;
        pPoolUB[i] = 2*nc;
    }
    else{
        pPoolLB[i] = 0;
        pPoolUB[i] = 0;
    }
    pPoolLB[i] -= delta;
    pPoolUB[i] += delta;
    //pPoolLB[i] = std::max(std::min(nc*(1 - delta), nc - 3.), 0.); // nc/(7 - delta - 2*(3 - delta)/(1 + 2./(nc+1))); //*(1 - delta);//
    //pPoolUB[i] = nc > 0 ? nc*(7 - delta - 2*(3 - delta)/(1 + 2./(nc+1))) : 3; //std::max(nc*(1 + delta), nc + 3.); // (1 + delta);//
}
////////////////////////////////////////////////////////////////////////////////

bool swmrssa::Comp::isOutOfBound(uint i, int nc)
{
    AssertLog(i < def()->countSpecs());
    if (nc > pPoolLB[i] && nc < pPoolUB[i]) {
        return false;
}
    setBounds(i, nc);
    return true;
}

double* swmrssa::Comp::pools(steps::wmrssa::PropensityRSSA prssa) const
{
    switch(prssa)
    {
        case steps::wmrssa::CURRENT:
            return def()->pools();
        case steps::wmrssa::LOWERBOUND:
            return pPoolLB;
        case steps::wmrssa::BOUNDS:
            return pPoolUB;
        default:
            AssertLog(false);
    }
}

void swmrssa::Comp::setupSpecDeps()
{
    uint nspecs = def()->countSpecs();
    localSpecUpdKProcs.resize(nspecs);
    for (uint slidx = 0; slidx < nspecs; slidx++) {
        uint sgidx = def()->specL2G(slidx);
        for (auto const& k : pKProcs) {
            if (k->depSpecComp(sgidx, this)) {
                localSpecUpdKProcs[slidx].push_back(k);
            }
        }
        for (auto const& ip : pIPatches) {
            for (auto const& k : ip->kprocs()) {
                if (k->depSpecComp(sgidx, this))
                    localSpecUpdKProcs[slidx].push_back(k);
            }
        }
        for (auto const& ip: pOPatches) {
            for (auto const& k : ip->kprocs()) {
                if (k->depSpecComp(sgidx, this))
                    localSpecUpdKProcs[slidx].push_back(k);
            }
        }
    }
}

std::vector<swmrssa::KProc*> const & swmrssa::Comp::getSpecUpdKProcs(uint slidx)
{
    return localSpecUpdKProcs[slidx];
}

// END
