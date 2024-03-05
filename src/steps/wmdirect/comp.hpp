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
#include "kproc.hpp"
#include "patch.hpp"
#include "solver/compdef.hpp"

namespace steps::wmdirect {

// Forward declarations.
class Comp;
class Wmdirect;

// Auxiliary declarations.
typedef Comp* CompP;
typedef std::vector<CompP> CompPVec;
typedef CompPVec::iterator CompPVecI;
typedef CompPVec::const_iterator CompPVecCI;

class Comp {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Comp(solver::Compdef* compdef, Wmdirect* solver);
    ~Comp();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////

    void setupKProcs(Wmdirect* wmd);
    void setupDeps();

    void reset();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    inline solver::Compdef* def() const noexcept {
        return pCompdef;
    }

    inline Wmdirect* solver() const noexcept {
        return pSolver;
    }

    ////////////////////////////////////////////////////////////////////////

    inline std::vector<steps::wmdirect::KProc*>::const_iterator begin() const noexcept {
        return pKProcs.begin();
    }
    inline std::vector<steps::wmdirect::KProc*>::const_iterator end() const noexcept {
        return pKProcs.end();
    }
    inline uint countKProcs() const noexcept {
        return static_cast<uint>(pKProcs.size());
    }
    inline const std::vector<steps::wmdirect::KProc*>& kprocs() const noexcept {
        return pKProcs;
    }
    inline std::vector<steps::wmdirect::KProc*>& kprocs() noexcept {
        return pKProcs;
    }

    steps::wmdirect::KProc* reac(solver::reac_local_id lridx) const;
    steps::wmdirect::KProc* reac(solver::complexreac_local_id lridx) const;

    ////////////////////////////////////////////////////////////////////////

    void addIPatch(Patch* p);

    inline uint countIPatches() const noexcept {
        return static_cast<uint>(pIPatches.size());
    }

    inline std::vector<steps::wmdirect::Patch*>::const_iterator beginIPatches() const noexcept {
        return pIPatches.begin();
    }
    inline std::vector<steps::wmdirect::Patch*>::const_iterator endIPatches() const noexcept {
        return pIPatches.end();
    }
    inline const std::vector<steps::wmdirect::Patch*>& ipatches() const noexcept {
        return pIPatches;
    }
    inline std::vector<steps::wmdirect::Patch*>& ipatches() noexcept {
        return pIPatches;
    }

    void addOPatch(Patch* p);

    inline uint countOPatches() const noexcept {
        return static_cast<uint>(pOPatches.size());
    }

    inline std::vector<steps::wmdirect::Patch*>::const_iterator beginOPatches() const noexcept {
        return pOPatches.begin();
    }
    inline std::vector<steps::wmdirect::Patch*>::const_iterator endOPatches() const noexcept {
        return pOPatches.end();
    }
    inline const std::vector<steps::wmdirect::Patch*>& opatches() const noexcept {
        return pOPatches;
    }
    inline std::vector<steps::wmdirect::Patch*>& opatches() noexcept {
        return pOPatches;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    Wmdirect* pSolver;

    solver::Compdef* pCompdef;

    /// The kinetic processes.
    std::vector<steps::wmdirect::KProc*> pKProcs;

    std::vector<steps::wmdirect::Patch*> pIPatches;
    std::vector<steps::wmdirect::Patch*> pOPatches;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::wmdirect
