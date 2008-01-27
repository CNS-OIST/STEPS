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

#ifndef STEPS_WMDIRECT_SOLVER_CORE_PATCH_HPP
#define STEPS_WMDIRECT_SOLVER_CORE_PATCH_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/patchdef.hpp>
#include <steps/wmdirect/solver_core/kproc.hpp>

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Comp;
class Patch;
class Sched;
class SReac;

// Auxiliary declarations.
typedef Patch *                         PatchP;
typedef std::vector<PatchP>             PatchPVec;
typedef PatchPVec::iterator             PatchPVecI;
typedef PatchPVec::const_iterator       PatchPVecCI;

////////////////////////////////////////////////////////////////////////////////

class Patch
{
    
public:
    
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    
    Patch(steps::sim::PatchDef * patchdef, Comp * icomp, Comp * ocomp);
    ~Patch(void);
    
    void setupKProcs(Sched * s);
    void setupDeps(void);
    
    void reset(void);
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////
    
    inline steps::sim::PatchDef * def(void) const
    { return pPatchDef; }
    
    inline double area(void) const
    { return pPatchDef->area(); }
    
    ////////////////////////////////////////////////////////////////////////
    
    inline uint * pools(void) const
    { return pPoolCount; }
    
    static const uint CLAMPED = 1;
    
    inline bool clamped(uint lidx) const
    { return pPoolFlags[lidx] & CLAMPED; }
    void setClamped(uint lidx, bool clamp);

    ////////////////////////////////////////////////////////////////////////

    inline KProcPVecCI kprocBegin(void) const
    { return pKProcs.begin(); }
    inline KProcPVecCI kprocEnd(void) const
    { return pKProcs.end(); }
    inline uint countKProcs(void) const
    { return pKProcs.size(); }
    
    SReac * sreac(uint lidx) const;
    
    ////////////////////////////////////////////////////////////////////////
    
    inline Comp * iComp(void) const
    { return pIComp; }
    
    inline Comp * oComp(void) const
    { return pOComp; }
    
    ////////////////////////////////////////////////////////////////////////
    
private:
    
    ////////////////////////////////////////////////////////////////////////
    
    steps::sim::PatchDef *      pPatchDef;
    
    /// Numbers of molecules -- stored as machine word integers.
    uint *                      pPoolCount;
    /// Flags on these pools -- stored as machine word flags.
    uint *                      pPoolFlags;
    
    /// The kinetic processes.
    KProcPVec                   pKProcs;
    
    Comp *                      pIComp;
    Comp *                      pOComp;
    
    ////////////////////////////////////////////////////////////////////////
    
    // Disable constructors.
    Patch(void);
    Patch(Patch const &);
    
    ////////////////////////////////////////////////////////////////////////
    
};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_WMDIRECT_SOLVER_CORE_PATCH_HPP

// END
