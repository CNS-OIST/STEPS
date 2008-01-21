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

#ifndef STEPS_TETEXACT_SOLVER_CORE_PATCH_HPP
#define STEPS_TETEXACT_SOLVER_CORE_PATCH_HPP 1

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
#include <steps/tetexact/solver_core/tri.hpp>

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Patch;

// Auxiliary declarations.
typedef Patch *                         PatchP;
typedef std::vector<PatchP>             PatchPVec;
typedef PatchPVec::iterator             PatchPVecI;
typedef PatchPVec::const_iterator       PatchPVecCI;


////////////////////////////////////////////////////////////////////////////////

/// Currently, this class is only a lightweight container for grouping 
/// a number of triangles. There is no internal ordering scheme. Tri
/// objects themselves do not at this point refer back to objects of this 
/// class (rather they refer directly to the PatchDef object for their 
/// color).
///
/// It also stores the total patch area as computed by summing the
/// areas of constituent triangles.
///
/// Later, when we transform the implementation of the Python-to-C-solver
/// interface into a Boost::Python/template based system, this class might
/// become the place in which certain patch-related concepts are
/// implemented (e.g. getting/setting the concentration in a whole
/// patch etc). When that happens, this class will probably become
/// related to PatchDef through some inheritance mechanism, removing the
/// double storage (and potentially double meaning) of the patch area.

class Patch
{
    
public:
    
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    
    Patch(steps::sim::PatchDef * patchdef);
    ~Patch(void);
    
    /// Checks whether Tri::patchdef() corresponds to this object's
    /// PatchDef. There is no check whether the Tri object has already
    /// been added to this Patch object before (i.e. no duplicate 
    /// checking).
    ///
    void addTri(Tri * tri);
    
    /// (Re-)computes the total area of the patch by adding the 
    /// triangular areas.
    ///
    void computeArea(void);
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////
    
    inline steps::sim::PatchDef * def(void) const
    { return pPatchDef; }
    
    inline double area(void) const
    { return pArea; }
    
    inline uint countTris(void) const
    { return pTris.size(); }
    
    inline TriPVecCI bgnTri(void) const
    { return pTris.begin(); }
    inline TriPVecCI endTri(void) const
    { return pTris.end(); }
    
    ////////////////////////////////////////////////////////////////////////
    
private:
    
    ////////////////////////////////////////////////////////////////////////
    
    steps::sim::PatchDef *      pPatchDef;
    double                      pArea;
    TriPVec                     pTris;

    ////////////////////////////////////////////////////////////////////////
    
    // Disable constructors.
    Patch(void);
    Patch(Patch const &);
    
    ////////////////////////////////////////////////////////////////////////
    
};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_TETEXACT_SOLVER_CORE_PATCH_HPP

// END
