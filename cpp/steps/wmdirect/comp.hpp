////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_WMDIRECT_COMP_HPP
#define STEPS_WMDIRECT_COMP_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/wmdirect/kproc.hpp>
#include <steps/wmdirect/patch.hpp>
#include <steps/solver/compdef.hpp>
#include <steps/solver/types.hpp>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(wmdirect)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Comp;
class Wmdirect;

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
    ~Comp(void);

    ////////////////////////////////////////////////////////////////////////

    void setupKProcs(Wmdirect * wmd);
    void setupDeps(void);

    void reset(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::Compdef * def(void) const
    { return pCompdef; }

    ////////////////////////////////////////////////////////////////////////

    inline std::vector<steps::wmdirect::KProc *>::const_iterator kprocBegin(void) const
    { return pKProcs.begin(); }
    inline std::vector<steps::wmdirect::KProc *>::const_iterator kprocEnd(void) const
    { return pKProcs.end(); }
    inline uint countKProcs(void) const
    { return pKProcs.size(); }

    steps::wmdirect::KProc * reac(uint lridx) const;

    ////////////////////////////////////////////////////////////////////////

    void addIPatch(Patch * p);

    inline uint countIPatches(void) const
    { return pIPatches.size(); }

    inline std::vector<steps::wmdirect::Patch *>::const_iterator beginIPatches(void) const
    { return pIPatches.begin(); }
    inline std::vector<steps::wmdirect::Patch *>::const_iterator  endIPatches(void) const
    { return pIPatches.end(); }

    void addOPatch(Patch * p);

    inline uint countOPatches(void) const
    { return pOPatches.size(); }

    inline std::vector<steps::wmdirect::Patch *>::const_iterator  beginOPatches(void) const
    { return pOPatches.begin(); }
    inline std::vector<steps::wmdirect::Patch *>::const_iterator  endOPatches(void) const
    { return pOPatches.end(); }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::Compdef *              pCompdef;

    /// The kinetic processes.
    std::vector<steps::wmdirect::KProc *> pKProcs;

    std::vector<steps::wmdirect::Patch *> pIPatches;
    std::vector<steps::wmdirect::Patch *> pOPatches;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(wmdirect)
END_NAMESPACE(steps)

#endif
// STEPS_WMDIRECT_COMP_HPP

// END

