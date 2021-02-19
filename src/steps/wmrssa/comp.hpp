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


#ifndef STEPS_WMRSSA_COMP_HPP
#define STEPS_WMRSSA_COMP_HPP 1


// STL headers.
#include <cassert>
#include <vector>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/wmrssa/kproc.hpp"
#include "steps/wmrssa/patch.hpp"
#include "steps/solver/compdef.hpp"
#include "steps/solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace wmrssa{

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Comp;
class Wmrssa;

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

    ////////////////////////////////////////////////////////////////////////

    void setupKProcs(Wmrssa * wmd);
    void setupDeps();

    void reset();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::Compdef * def() const
    { return pCompdef; }

    ////////////////////////////////////////////////////////////////////////

    inline std::vector<steps::wmrssa::KProc*>::const_iterator begin() const noexcept {
        return pKProcs.begin();
    }
    inline std::vector<steps::wmrssa::KProc*>::const_iterator end() const noexcept {
        return pKProcs.end();
    }
    inline uint countKProcs() const noexcept
    { return pKProcs.size(); }
    inline const std::vector<steps::wmrssa::KProc *> kprocs() const noexcept
    { return pKProcs; }

    steps::wmrssa::KProc * reac(uint lridx) const;

    ////////////////////////////////////////////////////////////////////////

    void addIPatch(Patch * p);

    inline uint countIPatches() const noexcept
    { return pIPatches.size(); }

    inline std::vector<steps::wmrssa::Patch *>::const_iterator beginIPatches() const noexcept
    { return pIPatches.begin(); }
    inline std::vector<steps::wmrssa::Patch *>::const_iterator  endIPatches() const noexcept
    { return pIPatches.end(); }
    inline const std::vector<steps::wmrssa::Patch *>& ipatches() const noexcept
    { return pIPatches; }

    void addOPatch(Patch * p);

    inline uint countOPatches() const noexcept
    { return pOPatches.size(); }

    inline std::vector<steps::wmrssa::Patch *>::const_iterator  beginOPatches() const noexcept
    { return pOPatches.begin(); }
    inline std::vector<steps::wmrssa::Patch *>::const_iterator  endOPatches() const noexcept
    { return pOPatches.end(); }
    inline const std::vector<steps::wmrssa::Patch *>& opatches() const noexcept
    { return pOPatches; }

    ////////////////////////////////////////////////////////////////////////

    void setBounds(uint i, int nc);
    bool isOutOfBound(uint i, int nc);
    double* pools(steps::wmrssa::PropensityRSSA prssa) const;
    void setupSpecDeps();
    std::vector<steps::wmrssa::KProc*> const & getSpecUpdKProcs(uint slidx);

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::Compdef *              pCompdef;

    /// The kinetic processes.
    std::vector<steps::wmrssa::KProc *> pKProcs;

    std::vector<steps::wmrssa::Patch *> pIPatches;
    std::vector<steps::wmrssa::Patch *> pOPatches;
    // Table of the lower bounds on populations of the species in this compartment.
    double                            * pPoolLB;
    // Table of the upper bounds on populations of the species in this compartment.
    double                            * pPoolUB;
    std::vector<std::vector<steps::wmrssa::KProc *>> localSpecUpdKProcs;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_WMRSSA_COMP_HPP

// END

