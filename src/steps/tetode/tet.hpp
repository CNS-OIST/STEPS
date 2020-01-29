/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
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



#ifndef STEPS_TETODE_TET_HPP
#define STEPS_TETODE_TET_HPP 1

// STL headers.
#include <array>
#include <cassert>
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include <steps/geom/fwd.hpp>
#include "steps/error.hpp"
#include "steps/solver/types.hpp"

// logging
#include <easylogging++.h>
////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace tetode {

////////////////////////////////////////////////////////////////////////////////

namespace stode = steps::tetode;

////////////////////////////////////////////////////////////////////////////////

// Forward declarations
class Tet;
class TetODE;
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
        tetrahedron_id_t idx, steps::solver::Compdef *cdef, double vol,
        double a0, double a1, double a2, double a3,
        double d0, double d1, double d2, double d3,
        tetrahedron_id_t tet0, tetrahedron_id_t tet1, tetrahedron_id_t tet2, tetrahedron_id_t tet3
      );
    ~Tet() = default;

    inline steps::solver::Compdef * compdef() const
    { return pCompdef; }

    ////////////////////////////////////////////////////////////////////////
    // SETUP
    ////////////////////////////////////////////////////////////////////////


    /// Set pointer to the next neighbouring tetrahedron.
    ///
    void setNextTet(uint i, stode::Tet * t);

    void setNextTri(uint i, stode::Tri *t);

    /* Not needed?
    // This method only asserts this method is not called on derived object
    void setNextTri(stode::Tri *t);
    */

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // SHAPE & CONNECTIVITY INFORMATION.
    ////////////////////////////////////////////////////////////////////////

    /// Get the volume.
    ///
    inline double vol() const noexcept
    { return pVol; }

    /// get the index
    inline tetrahedron_id_t idx() const  noexcept
    { return pIdx; }

    /// Get pointer to the next neighbouring triangle.
    ///
    inline stode::Tri * nextTri(uint i) const
    {
        AssertLog(i < 4);
        return pNextTri[i];
    }

    /// Get pointer to the next neighbouring tetrahedron.
    ///
    inline stode::Tet * nextTet(uint i) const noexcept
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

    ////////////////////////////////////////////////////////////////////////

    /*
    // Set whether a direction is a diffusion boundary
    void setDiffBndDirection(uint i);

    inline bool getDiffBndDirection(uint idx) const
    { return pDiffBndDirection[idx]; }
    */


    inline tetrahedron_id_t tet(uint t) const noexcept
    { return pTets[t]; }

    /*
    inline uint tri(uint t)
    { return pTris[t]; }
    */

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::Compdef            * pCompdef;

    tetrahedron_id_t                    pIdx;

    double                              pVol;

    // Indices of neighbouring tetrahedra.
    std::array<tetrahedron_id_t, 4>     pTets;

    // Indices of neighbouring triangles.
    //uint                                pTris[4];



    /// Pointers to neighbouring triangles
    stode::Tri                         * pNextTri[4];


    /// Pointers to neighbouring tetrahedra.
    stode::Tet                         * pNextTet[4];

    double                              pAreas[4];
    double                              pDist[4];

    // bool                              pDiffBndDirection[4];

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_TETODE_TET_HPP

// END
