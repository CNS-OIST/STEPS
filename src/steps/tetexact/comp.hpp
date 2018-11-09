/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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


#ifndef STEPS_TETEXACT_COMP_HPP
#define STEPS_TETEXACT_COMP_HPP 1

// Standard headers.
#include <cassert>
#include <vector>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/tetexact/tet.hpp"
#include "steps/solver/compdef.hpp"
#include "steps/solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace tetexact {

////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Comp;

// Auxiliary declarations.
typedef Comp *                          CompP;
typedef std::vector<CompP>              CompPVec;
typedef CompPVec::iterator              CompPVecI;
typedef CompPVec::const_iterator        CompPVecCI;

////////////////////////////////////////////////////////////////////////////////

class Comp
{
public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Comp(steps::solver::Compdef * compdef);
    ~Comp();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    /// Checks whether the Tet's compdef() corresponds to this object's
    /// CompDef. There is no check whether the Tet object has already
    /// been added to this Comp object before (i.e. no duplicate checking).
    ///
    void addTet(stex::WmVol * tet);

    ////////////////////////////////////////////////////////////////////////

    void reset() { def()->reset(); }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    steps::solver::Compdef *def() const { return pCompdef; }

    double vol() const { return pVol; }

    double * pools() const { return def()->pools(); }

    void modCount(uint slidx, double count);

    uint countTets() const { return pTets.size(); }

    stex::WmVol * pickTetByVol(double rand01) const;

    WmVolPVecCI bgnTet() const { return pTets.begin(); }
    WmVolPVecCI endTet() const { return pTets.end(); }
    const WmVolPVec &tets() const { return pTets; }

    ////////////////////////////////////////////////////////////////////////

private:

    steps::solver::Compdef                    * pCompdef;
    double                                     pVol;

    WmVolPVec                                pTets;
};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif

// STEPS_TETEXACT_COMP_HPP

// END
