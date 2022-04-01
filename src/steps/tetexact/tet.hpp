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


#ifndef STEPS_TETEXACT_TET_HPP
#define STEPS_TETEXACT_TET_HPP 1

// STL headers.
#include <cassert>
#include <vector>

// logging
#include <easylogging++.h>

// STEPS headers.
#include "util/common.h"
#include "kproc.hpp"
#include "wmvol.hpp"
#include "solver/compdef.hpp"
#include "solver/types.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace tetexact{

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

class Tet: public WmVol
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Tet
      (
        tetrahedron_id_t idx, steps::solver::Compdef *cdef, double vol,
        double a0, double a1, double a2, double a3,
        double d0, double d1, double d2, double d3,
        tetrahedron_id_t tet0, tetrahedron_id_t tet1, tetrahedron_id_t tet2, tetrahedron_id_t tet3
      );
    ~Tet() override;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file) override;

    /// restore data
    void restore(std::fstream & cp_file) override;

    ////////////////////////////////////////////////////////////////////////
    // SETUP
    ////////////////////////////////////////////////////////////////////////


    /// Set pointer to the next neighbouring tetrahedron.
    ///
    void setNextTet(uint i, Tet * t);


    /// Set pointer to the next neighbouring triangle.
    void setNextTri(uint i, Tri *t);

    // This method only asserts this method is not called on derived object
    void setNextTri(Tri *t) override;


    /// Create the kinetic processes -- to be called when all tetrahedrons
    /// and triangles have been fully declared and connected.
    ///
    void setupKProcs(Tetexact * tex) override;


    ////////////////////////////////////////////////////////////////////////
    // SHAPE & CONNECTIVITY INFORMATION.
    ////////////////////////////////////////////////////////////////////////

    /// Get pointer to the next neighbouring triangle.
    ///
    inline Tri * nextTri(uint i) const
    {
        AssertLog(i < 4);
        return pNextTris[i];
    }

    /// Get pointer to the next neighbouring tetrahedron.
    ///
    inline Tet * nextTet(uint i) const noexcept
    { return pNextTet[i]; }

    /// Get the area of a boundary triangle.
    ///
    inline double area(uint i) const noexcept
    { return pAreas[i]; }

    /// Get the distance to the centroid of the next neighbouring
    /// tetrahedron.
    ///
    inline double dist(uint i) const noexcept
    { return pDist[i]; }

    /// Find the direction index towards a neighbor tetrahedron.
    ///
    int getTetDirection(tetrahedron_id_t tidx) const;

    ////////////////////////////////////////////////////////////////////////

    // Set whether a direction is a diffusion boundary
    void setDiffBndDirection(uint i);

    inline bool getDiffBndDirection(uint idx) const noexcept
    { return pDiffBndDirection[idx]; }


    Diff * diff(uint lidx) const;

    inline tetrahedron_id_t tet(uint t) const noexcept
    { return pTets[t]; }

    /*
    inline uint tri(uint t)
    { return pTris[t]; }
    */

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    // Indices of neighbouring tetrahedra.
    tetrahedron_id_t                    pTets[4];

    /// Pointers to neighbouring tetrahedra.
    Tet                         * pNextTet[4];

    double                              pAreas[4];
    double                              pDist[4];

    bool                                pDiffBndDirection[4];


    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_TETEXACT_TET_HPP

// END
