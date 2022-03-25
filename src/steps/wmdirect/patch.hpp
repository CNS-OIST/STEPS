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


#ifndef STEPS_WMDIRECT_PATCH_HPP
#define STEPS_WMDIRECT_PATCH_HPP 1


// STL headers.
#include <cassert>
#include <vector>
#include <fstream>

// STEPS headers.
#include "kproc.hpp"
#include "comp.hpp"
#include "sreac.hpp"
#include "util/common.h"
#include "solver/patchdef.hpp"
#include "solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps {
namespace wmdirect {

////////////////////////////////////////////////////////////////////////////////

namespace swmd = steps::wmdirect;

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Patch;
class Wmdirect;

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

    Patch(steps::solver::Patchdef * patchdef, swmd::Comp * icomp, swmd::Comp * ocomp);
    ~Patch();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////

    void setupKProcs(Wmdirect * wmd);
    void setupDeps();

    void reset();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::Patchdef * def() const noexcept
    { return pPatchdef; }

    ////////////////////////////////////////////////////////////////////////

    inline std::vector<steps::wmdirect::KProc*>::const_iterator begin() const noexcept {
        return pKProcs.begin();
    }
    inline std::vector<steps::wmdirect::KProc*>::const_iterator end() const noexcept {
        return pKProcs.end();
    }
    inline uint countKProcs() const noexcept
    { return static_cast<uint>(pKProcs.size()); }
    inline const std::vector<steps::wmdirect::KProc *>& kprocs() const noexcept {
        return pKProcs;
    }
    inline std::vector<steps::wmdirect::KProc *>& kprocs() noexcept {
      return pKProcs;
    }

    steps::wmdirect::KProc * sreac(uint lsridx) const;

    ////////////////////////////////////////////////////////////////////////

    inline swmd::Comp * iComp() const noexcept
    { return pIComp; }

    inline swmd::Comp * oComp() const noexcept
    { return pOComp; }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::Patchdef             * pPatchdef;

    /// The kinetic processes.
    std::vector<steps::wmdirect::KProc *> pKProcs;

    swmd::Comp                          * pIComp;
    swmd::Comp                          * pOComp;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_WMDIRECT_PATCH_HPP

// END

