/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/math/constants.hpp"
#include "steps/solver/compdef.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/solver/reacdef.hpp"
#include "steps/tetexact/diff.hpp"
#include "steps/tetexact/kproc.hpp"
#include "steps/tetexact/reac.hpp"
#include "steps/tetexact/tet.hpp"
#include "steps/tetexact/tetexact.hpp"
#include "steps/tetexact/tri.hpp"
#include "steps/tetexact/wmvol.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

stex::WmVol::WmVol
  (
    tetrahedron_id_t idx, solver::Compdef *cdef, double vol
  )
: pIdx(idx)
, pCompdef(cdef)
, pVol(vol)
{
    AssertLog(pCompdef != nullptr);
    AssertLog(pVol > 0.0);

    // Based on compartment definition, build other structures.
    auto nspecs = compdef()->countSpecs();
    pPoolCount.resize(nspecs, 0);
    pPoolFlags.resize(nspecs, 0);
    pKProcs.resize(compdef()->countReacs());

}

////////////////////////////////////////////////////////////////////////////////

stex::WmVol::~WmVol()
{
    // Delete reaction rules.
    for (auto const& i : pKProcs) {
      delete i;
    }
}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::checkpoint(std::fstream & cp_file)
{
    steps::checkpoint(cp_file, pPoolCount, false /* with_size */);
    steps::checkpoint(cp_file, pPoolFlags, false /* with_size */);
}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::restore(std::fstream & cp_file)
{
    const auto nspecs = compdef()->countSpecs();
    steps::restore(cp_file, nspecs, pPoolCount);
    steps::restore(cp_file, nspecs, pPoolFlags);
}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::setNextTri(stex::Tri * t)
{
    pNextTris.push_back(t);
}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::setupKProcs(stex::Tetexact * tex)
{

    uint j = 0;

    // Note: ignoring diffusion KProcs

    // Create reaction kproc's.
    uint nreacs = compdef()->countReacs();
    for (uint i = 0; i < nreacs; ++i)
    {
        ssolver::Reacdef * rdef = compdef()->reacdef(i);
        auto * r = new stex::Reac(rdef, this);
        pKProcs[j++] = r;
        tex->addKProc(r);
    }

}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::reset()
{
    std::fill(pPoolCount.begin(), pPoolCount.end(), 0);
    std::fill(pPoolFlags.begin(), pPoolFlags.end(), 0);
    for (auto kproc: pKProcs) {
        kproc->reset();
    }
}

////////////////////////////////////////////////////////////////////////////////

double stex::WmVol::conc(uint gidx) const
{
    uint lspidx = compdef()->specG2L(gidx);
    double n = pPoolCount[lspidx];
    return (n/(1.0e3*pVol*steps::math::AVOGADRO));
}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::setCount(uint lidx, uint count)
{
    AssertLog(lidx < compdef()->countSpecs());
    pPoolCount[lidx] = count;
}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::incCount(uint lidx, int inc)
{
    AssertLog(lidx < compdef()->countSpecs());
#ifndef NDEBUG
    uint old_count = pPoolCount[lidx];
#endif

    pPoolCount[lidx] += inc;

#ifndef NDEBUG
    uint new_count = pPoolCount[lidx];
    AssertLog((inc >= 0 && new_count >= old_count) || (inc < 0 && new_count < old_count));
#endif
}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::setClamped(uint lidx, bool clamp)
{
    if (clamp) { pPoolFlags[lidx] |= CLAMPED;
    } else { pPoolFlags[lidx] &= ~CLAMPED;
}
}

////////////////////////////////////////////////////////////////////////////////

stex::Reac * stex::WmVol::reac(uint lidx) const
{
    AssertLog(lidx < compdef()->countReacs());
    return dynamic_cast<stex::Reac*>(pKProcs[lidx]);
}

////////////////////////////////////////////////////////////////////////////////

// END
