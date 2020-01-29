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

#ifndef STEPS_MPI_TETOPSPLIT_COMP_HPP
#define STEPS_MPI_TETOPSPLIT_COMP_HPP 1


// STL headers.
#include <cassert>
#include <vector>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/mpi/tetopsplit/tet.hpp"
#include "steps/solver/compdef.hpp"
#include "steps/solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace mpi{
namespace tetopsplit{

////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;

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

    explicit Comp(steps::solver::Compdef * compdef);

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
    void addTet(smtos::WmVol * tet);

    ////////////////////////////////////////////////////////////////////////

    inline void reset()
    { def()->reset(); }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::Compdef * def() const noexcept
    { return pCompdef; }

    inline double vol() const noexcept
    { return pVol; }

    inline double * pools() const noexcept
    { return def()->pools(); }

    void modCount(uint slidx, double count);

    inline uint countTets() const noexcept
    { return pTets.size(); }

    smtos::WmVol * pickTetByVol(double rand01) const;

    inline WmVolPVecCI bgnTet() const noexcept
    { return pTets.begin(); }
    inline WmVolPVecCI endTet() const noexcept
    { return pTets.end(); }
    inline const WmVolPVec& tets() const noexcept
    { return pTets; }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::Compdef                    * pCompdef;
    double                                     pVol{0.0};

    WmVolPVec                                pTets;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

#endif

// STEPS_MPI_TETOPSPLIT_COMP_HPP

// END
