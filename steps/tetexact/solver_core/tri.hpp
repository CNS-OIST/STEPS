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

#ifndef STEPS_TETEXACT_SOLVER_CORE_TRI_HPP
#define STEPS_TETEXACT_SOLVER_CORE_TRI_HPP 1

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
#include <steps/tetexact/solver_core/kproc.hpp>

// Forward declarations.
class Sched;
class SReac;
class Tet;
class Tri;

// Auxiliary declarations.
typedef Tri *                           TriP;
typedef std::vector<TriP>               TriPVec;
typedef TriPVec::iterator               TriPVecI;
typedef TriPVec::const_iterator         TriPVecCI;

////////////////////////////////////////////////////////////////////////////////
 
class Tri
{

public:

    ////////////////////////////////////////////////////////////////////////
    
    /// Constructor.
    Tri(steps::sim::PatchDef * pdef, double area);
    
    /// Destructor.
    ~Tri(void);

    ////////////////////////////////////////////////////////////////////////
    // TRIANGLE SETUP
    ////////////////////////////////////////////////////////////////////////
    
    /// Set pointer to the 'inside' neighbouring tetrahedron.
    ///
    void setInnerTet(Tet * t);
    
    /// Set pointer to the 'outside' neighbouring tetrahedron.
    ///
    void setOuterTet(Tet * t);
    
    /// Create the kinetic processes -- to be called when all tetrahedrons 
    /// and triangles have been fully declared and connected.
    ///
    void setupKProcs(Sched * s);
    
    /// Once all the kinetic processes have been created in all tetrahedrons
    /// and triangles, their interdepencies are pre-computed so that they
    /// are immediately available during actual simulation. This is the
    /// final step of the initialization. (Except of course calling 
    /// reset()).
    ///
    void setupDeps(void);
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: GENERAL
    ////////////////////////////////////////////////////////////////////////
    
    inline steps::sim::PatchDef * patchdef(void) const
    { return pPatchDef; }
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SHAPE & CONNECTIVITY
    ////////////////////////////////////////////////////////////////////////
    
    inline double area(void) const
    { return pArea; }
    
    inline Tet * iTet(void) const
    { return pInnerTet; }
    
    inline Tet * oTet(void) const
    { return pOuterTet; }
    
    ////////////////////////////////////////////////////////////////////////
    // MAIN FUNCTIONALITY
    ////////////////////////////////////////////////////////////////////////
    
    /// Set all pool flags and molecular populations to zero.
    void reset(void);
    
    inline uint * pools(void) const
    { return pPoolCount; }
    
    ///inline uint poolCount(uint lidx) const
    ///{ return pPoolCount[lidx]; }
    ///inline void setPoolCount(uint lidx, uint num)
    ///{ pPoolCount[lidx] = num; }
    ///inline void incPoolCount(uint lidx, int count) const
    ///{ pPoolCount[lidx] += count; }
    ///inline uint poolFlags(uint lidx) const
    ///{ return pPoolFlags[lidx]; }

    ////////////////////////////////////////////////////////////////////////

    inline KProcPVecCI kprocBegin(void) const
    { return pKProcs.begin(); }
    inline KProcPVecCI kprocEnd(void) const
    { return pKProcs.end(); }
    inline uint countKProcs(void) const
    { return pKProcs.size(); }
    
    SReac * sreac(uint lidx) const;
    
    ////////////////////////////////////////////////////////////////////////
    
private:

    ////////////////////////////////////////////////////////////////////////
    
    /// Type of tetrahedron.
    steps::sim::PatchDef *      pPatchDef;
    
    /// Pointers to neighbouring tetrahedra.
    Tet *                       pInnerTet;
    Tet *                       pOuterTet;
    
    double                      pArea;
    
    /// Numbers of molecules -- stored as machine word integers.
    uint *                      pPoolCount;
    /// Flags on these pools -- stored as machine word flags.
    uint *                      pPoolFlags;
    
    /// The kinetic processes.
    KProcPVec                   pKProcs;
    
    ////////////////////////////////////////////////////////////////////////
    
};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_TETEXACT_SOLVER_CORE_TRI_HPP

// END
