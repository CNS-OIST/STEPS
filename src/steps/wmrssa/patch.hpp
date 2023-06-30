/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

#pragma once

// STL headers.
#include <cassert>
#include <fstream>
#include <vector>

// STEPS headers.
#include "comp.hpp"
#include "kproc.hpp"
#include "solver/patchdef.hpp"
#include "solver/types.hpp"
#include "sreac.hpp"
#include "util/common.hpp"

namespace steps::wmrssa {

// Forward declarations.
class Patch;
class Wmrssa;

// Auxiliary declarations.
typedef Patch* PatchP;
typedef std::vector<PatchP> PatchPVec;
typedef PatchPVec::iterator PatchPVecI;
typedef PatchPVec::const_iterator PatchPVecCI;

////////////////////////////////////////////////////////////////////////////////

class Patch {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Patch(solver::Patchdef* patchdef, Comp* icomp, Comp* ocomp);
    ~Patch();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////

    void setupKProcs(Wmrssa* wmd);
    void setupDeps();

    void reset();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    inline solver::Patchdef* def() const {
        return pPatchdef;
    }

    ////////////////////////////////////////////////////////////////////////

    inline std::vector<wmrssa::KProc*>::const_iterator begin() const {
        return pKProcs.begin();
    }
    inline std::vector<wmrssa::KProc*>::const_iterator end() const {
        return pKProcs.end();
    }
    inline uint countKProcs() const {
        return pKProcs.size();
    }
    inline const std::vector<wmrssa::KProc*>& kprocs() const noexcept {
        return pKProcs;
    }
    wmrssa::KProc* sreac(solver::sreac_local_id lsridx) const;

    ////////////////////////////////////////////////////////////////////////

    inline Comp* iComp() const {
        return pIComp;
    }

    inline Comp* oComp() const {
        return pOComp;
    }

    ////////////////////////////////////////////////////////////////////////

    void setBounds(solver::spec_local_id i, int nc);
    bool isOutOfBound(solver::spec_local_id i, int nc);
    const util::strongid_vector<solver::spec_local_id, double>& pools(
        wmrssa::PropensityRSSA prssa) const;
    void setupSpecDeps();
    inline std::vector<wmrssa::KProc*> const& getSpecUpdKProcs(
        solver::spec_local_id slidx) noexcept {
        return localSpecUpdKProcs[slidx.get()];
    }

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::Patchdef* pPatchdef;

    /// The kinetic processes.
    std::vector<wmrssa::KProc*> pKProcs;

    Comp* pIComp;
    Comp* pOComp;
    // Table of the lower bounds on populations of the species in this compartment.
    util::strongid_vector<solver::spec_local_id, double> pPoolLB;
    // Table of the upper bounds on populations of the species in this compartment.
    util::strongid_vector<solver::spec_local_id, double> pPoolUB;
    std::vector<std::vector<wmrssa::KProc*>> localSpecUpdKProcs;
};

}  // namespace steps::wmrssa
