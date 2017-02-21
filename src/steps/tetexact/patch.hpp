/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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


#ifndef STEPS_TETEXACT_PATCH_HPP
#define STEPS_TETEXACT_PATCH_HPP 1


// STL headers.
#include <cassert>
#include <vector>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/patchdef.hpp"
#include "steps/tetexact/tri.hpp"
#include "steps/solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace tetexact {

////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;

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
    ~Patch(void);

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
    void addTri(stex::Tri * tri);

    ////////////////////////////////////////////////////////////////////////

    inline void reset(void)
    { def()->reset(); }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::Patchdef * def(void) const
    { return pPatchdef; }

    inline double area(void) const
    { return pArea; }

    //void setArea(double a);

    inline double * pools(void) const
    { return def()->pools(); }

    void modCount(uint slidx, double count);


    inline uint countTris(void) const
    { return pTris.size(); }

    stex::Tri * pickTriByArea(double rand01) const;

    inline TriPVecCI bgnTri(void) const
    { return pTris.begin(); }
    inline TriPVecCI endTri(void) const
    { return pTris.end(); }

    inline const TriPVec &tris() const { return pTris; }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::Patchdef           * pPatchdef;
    double                              pArea;

    TriPVec                             pTris;

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_TETEXACT_PATCH_HPP

// END
