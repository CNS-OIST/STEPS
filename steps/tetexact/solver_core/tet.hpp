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
// $Id:tet.hpp 64 2007-08-20 06:25:41Z stefan $
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_TETEXACT_SOLVER_CORE_TET_HPP
#define STEPS_TETEXACT_SOLVER_CORE_TET_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/compdef.hpp>
#include <steps/tetexact/solver_core/kproc.hpp>

// Forward declarations.
class Diff;
class Reac;
class Sched;
class Tet;
class Tri;

// Auxiliary declarations.
typedef Tet *                           TetP;
typedef std::vector<TetP>               TetPVec;
typedef TetPVec::iterator               TetPVecI;
typedef TetPVec::const_iterator         TetPVecCI;

////////////////////////////////////////////////////////////////////////////////
 
class Tet
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    
    Tet
    (
        steps::sim::CompDef * cdef, double vol, 
        double a0, double a1, double a2, double a3, 
        double d0, double d1, double d2, double d3
    );
    ~Tet(void);

    /// Set pointer to the next neighbouring tetrahedron.
    ///
    void setNextTet(uint i, Tet * t);
    
    /// Set pointer to the next neighbouring triangle.
    ///
    void setNextTri(uint i, Tri * t);
    
    /// Create the kinetic processes -- to be called when all tetrahedrons 
    /// and triangles have been fully declared and connected.
    ///
    void setupKProcs(Sched * s);
    
    /// Once all the kinetic processes have been created in all tetraedrons
    /// and triangles, their interdepencies are pre-computed so that they
    /// are immediately available during actual simulation. This is the
    /// final step of the initialization. (Except of course calling 
    /// reset()).
    ///
    void setupDeps(void);
    
    ////////////////////////////////////////////////////////////////////////
    
    void reset(void);
    
    ////////////////////////////////////////////////////////////////////////
    // GENERAL INFORMATION
    ////////////////////////////////////////////////////////////////////////
    
    inline steps::sim::CompDef * compdef(void) const
    { return pCompDef; }
    
    ////////////////////////////////////////////////////////////////////////
    // SHAPE & CONNECTIVITY INFORMATION.
    ////////////////////////////////////////////////////////////////////////
    
    /// Get pointer to the next neighbouring tetrahedron.
    ///
    inline Tet * nextTet(uint i) const
    { return pNextTet[i]; }
    
    /// Get pointer to the next neighbouring triangle.
    ///
    inline Tri * nextTri(uint i) const
    { return pNextTri[i]; }
    
    /// Get the volume.
    ///
    inline double vol(void) const
    { return pVol; }
    
    /// Get the area of a boundary triangle.
    ///
    inline double area(uint i) const
    { return pAreas[i]; }
    
    /// Get the distance to the centroid of the next neighbouring 
    /// tetrahedron. 
    ///
    inline double dist(uint i) const
    { return pDist[i]; }
    
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
    
    Diff * diff(uint lidx) const;
    Reac * reac(uint lidx) const;
    
    ////////////////////////////////////////////////////////////////////////
    
private:

    ////////////////////////////////////////////////////////////////////////
    
    /// Type of tetrahedron.
    steps::sim::CompDef *       pCompDef;
    
    /// Pointers to neighbouring tetrahedra.
    Tet *                       pNextTet[4];
    /// Pointers to neighbouring triangles.
    Tri *                       pNextTri[4];
    
    double                      pVol;
    double                      pAreas[4];
    double                      pDist[4];
    
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
// STEPS_TETEXACT_SOLVER_CORE_TET_HPP

// END
