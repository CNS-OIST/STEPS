////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
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
//
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_TETEXACT_TET_HPP
#define STEPS_TETEXACT_TET_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/solver/compdef.hpp>
#include <steps/tetexact/kproc.hpp>
#include <steps/solver/types.hpp>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(tetexact)

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::tetexact, stex);

////////////////////////////////////////////////////////////////////////////////

// Forward declarations
class Tet;
class Diff;
class Tri;
class Reac;
class Tetexact;

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
        steps::solver::Compdef * cdef, double vol,
        double a0, double a1, double a2, double a3,
        double d0, double d1, double d2, double d3,
        int tet0, int tet1, int tet2, int tet3
    );
    ~Tet(void);

    ////////////////////////////////////////////////////////////////////////
    // SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Set pointer to the next neighbouring tetrahedron.
    ///
    void setNextTet(uint i, stex::Tet * t);

    /// Set pointer to the next neighbouring triangle.
    ///
    void setNextTri(uint i, stex::Tri * t);

    /// Create the kinetic processes -- to be called when all tetrahedrons
    /// and triangles have been fully declared and connected.
    ///
    void setupKProcs(stex::Tetexact * tex);

    ////////////////////////////////////////////////////////////////////////

    void reset(void);

    ////////////////////////////////////////////////////////////////////////
    // GENERAL INFORMATION
    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::Compdef * compdef(void) const
    { return pCompdef; }

    ////////////////////////////////////////////////////////////////////////
    // SHAPE & CONNECTIVITY INFORMATION.
    ////////////////////////////////////////////////////////////////////////

    /// Get pointer to the next neighbouring tetrahedron.
    ///
    inline stex::Tet * nextTet(uint i) const
    { return pNextTet[i]; }

    /// Get pointer to the next neighbouring triangle.
    ///
    inline stex::Tri * nextTri(uint i) const
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
    void setCount(uint lidx, uint count);

    static const uint CLAMPED = 1;

    inline bool clamped(uint lidx) const
    { return pPoolFlags[lidx] & CLAMPED; }
    void setClamped(uint lidx, bool clamp);

    ////////////////////////////////////////////////////////////////////////

    inline std::vector<stex::KProc *>::const_iterator kprocBegin(void) const
    { return pKProcs.begin(); }
    inline std::vector<stex::KProc *>::const_iterator kprocEnd(void) const
    { return pKProcs.end(); }
    inline uint countKProcs(void) const
    { return pKProcs.size(); }

    stex::Diff * diff(uint lidx) const;
    stex::Reac * reac(uint lidx) const;

    inline int tet(uint t)
    { return pTets[t]; }

    /*
    inline uint tri(uint t)
    { return pTris[t]; }
	*/

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    /// Type of tetrahedron.
    steps::solver::Compdef            * pCompdef;

    // Indices of neighbouring tetrahedra.
    int                                 pTets[4];
    /*
    // Indices of neighbouring triangles.
    uint                                pTris[4];
	*/

    /// Pointers to neighbouring tetrahedra.
    stex::Tet                         * pNextTet[4];
    /// Pointers to neighbouring triangles.
    stex::Tri                         * pNextTri[4];

    double                              pVol;
    double                              pAreas[4];
    double                              pDist[4];

    /// Numbers of molecules -- stored as uint.
    uint                              * pPoolCount;
    /// Flags on these pools -- stored as machine word flags.
    uint                              * pPoolFlags;

    /// The kinetic processes.
    std::vector<stex::KProc *>          pKProcs;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetexact)
END_NAMESPACE(steps)

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_TETEXACT_TET_HPP

// END
