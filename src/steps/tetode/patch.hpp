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


#ifndef STEPS_TETODE_PATCH_HPP
#define STEPS_TETODE_PATCH_HPP 1


// STL headers.
#include <cassert>
#include <vector>
#include <fstream>
#include <map>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/patchdef.hpp"
#include "steps/tetode/tri.hpp"
#include "steps/solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace tetode {

////////////////////////////////////////////////////////////////////////////////

namespace stode = steps::tetode;

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Patch;

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

    Patch(steps::solver::Patchdef * patchdef);
    ~Patch();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    /// Checks whether Tri::patchdef() corresponds to this object's
    /// PatchDef. There is no check whether the Tri object has already
    /// been added to this Patch object before (i.e. no duplicate
    /// checking).
    ///
    void addTri(stode::Tri * tri);

    ////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    inline double area() const noexcept
    { return pArea; }

    // Return the local index of a tri given by global index
    triangle_id_t getTri_GtoL(triangle_id_t gidx);

    // Return the tri of a given local index
    Tri * getTri(uint lidx);

    inline steps::solver::Patchdef * def() const noexcept
    { return pPatchdef; }

    inline uint countTris() const noexcept
    { return pTris.size(); }


    inline TriPVecCI bgnTri() const noexcept
    { return pTris.begin(); }
    inline TriPVecCI endTri() const noexcept
    { return pTris.end(); }
    inline const TriPVec& tris() const noexcept
    { return pTris; }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::Patchdef           * pPatchdef;

    TriPVec                             pTris;

    double                                 pArea;

    // A map storing global index to local
    std::map<triangle_id_t, triangle_id_t>                    pTris_GtoL;

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_TETODE_PATCH_HPP

// END
