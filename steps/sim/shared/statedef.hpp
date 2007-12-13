////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// This file is part of STEPS.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SIM_SHARED_STATEDEF_HPP
#define STEPS_SIM_SHARED_STATEDEF_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/types.hpp>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(sim)

// Forward declarations.
class CompDef;
class DiffDef;
class PatchDef;
class ReacDef;
class SpecDef;
class SReacDef;

////////////////////////////////////////////////////////////////////////////////

/// This class was designed for fast lookup during simulation
/// rather than accomodating change; the latter is the task of the Python
/// objects in package steps.model.
///
class StateDef
{

public:
    
    enum Mode
    {
        NAIVE_MODE,
        SETUP_MODE,
        SETUP_VAR_MODE,
        SETUP_REAC_MODE,
        SETUP_DIFF_MODE,
        SETUP_COMP_MODE,
        READY_MODE
    };
    
    ////////////////////////////////////////////////////////////////////////
    
    /// Default constructor.
    StateDef(void);
    /// Destructor.
    ~StateDef(void);

    ////////////////////////////////////////////////////////////////////////

    Mode mode(void) const 
    {
        return pMode;
    }
    
    /// We're not checking whether the specified mode transition makes
    /// sense; this should be left to the discretion of the caller.
    /// The reason for this, is that the semantics of this might change
    /// in the future, when solvers might allow changes in the state
    /// (new compartments, ...) during a simulation.
    ///
    void setMode(Mode m)
    {
        pMode = m;
    }
    
    ////////////////////////////////////////////////////////////////////////
    
    /// Gets called when all components in the entire state description
    /// have been fully defined. This allows creation of local state
    /// descriptions for compartments and patches, and resolution of
    /// interdependencies between various compartments and patches.
    ///
    /// Can be called only once!
    ///
    void setupFinal(void);
    
    ////////////////////////////////////////////////////////////////////////

    ///
    SpecDef * createSpecDef(std::string const & name);
    
    /// Returns the number of species in the overall simulation state.
    ///
    uint countSpecs(void) const
    { return pSpecs.size(); }
    
    /// Check whether the specified global index refers to a valid species.
    ///
    bool isValidSpec(gidxT idx) const
    { return (idx < countSpecs()); }
    
    /// Fetches a species by its idx.
    ///
    SpecDef * spec(gidxT idx) const
    { return (isValidSpec(idx) ? pSpecs[idx] : 0 ); }
    
    ///
    std::vector<SpecDef*>::const_iterator beginSpec(void) const
    { return pSpecs.begin(); }
    
    ///
    std::vector<SpecDef*>::const_iterator endSpec(void) const
    { return pSpecs.end(); }

    ////////////////////////////////////////////////////////////////////////

    /// 
    CompDef * createCompDef(std::string const & name);
    
    /// Returns the number of compartments in the simulation state.
    uint countComps(void) const
    { return pComps.size(); }
    
    /// Check whether the specified global index refers to a valid compartment.
    bool isValidComp(gidxT idx) const
    { return (idx < countComps()); }
    
    /// Fetches a compartment by its idx.
    CompDef * comp(gidxT idx) const
    { return (isValidComp(idx) ? pComps[idx] : 0); }

    ///
    std::vector<CompDef*>::const_iterator beginComp(void) const
    { return pComps.begin(); }
    
    ///
    std::vector<CompDef*>::const_iterator endComp(void) const
    { return pComps.end(); }
    
    ////////////////////////////////////////////////////////////////////////

    ///
    ReacDef * createReacDef(std::string const & name);
    
    /// Returns the number of reactions in the overall simulation state.
    uint countReacs(void) const
    { return pReacs.size(); }
    
    /// Check whether the specified global index refers to a valid reaction.
    bool isValidReac(gidxT idx) const
    { return (idx < countReacs()); }
    
    /// Fetches a reaction by its idx.
    ReacDef * reac(gidxT idx) const
    { return (isValidReac(idx) ? pReacs[idx] : 0 ); }
    
    ///
    std::vector<ReacDef*>::const_iterator beginReac(void) const
    { return pReacs.begin(); }
    
    ///
    std::vector<ReacDef*>::const_iterator endReac(void) const
    { return pReacs.end(); }
    
    ////////////////////////////////////////////////////////////////////////
    
    ///
    DiffDef * createDiffDef(std::string const & name);
    
    /// Returns the number of diffusion rules in the overall 
    /// simulation state.
    uint countDiffs(void) const
    { return pDiffs.size(); }
    
    /// Check whether the specified global index refers to a valid 
    /// diffusion rule.
    bool isValidDiff(gidxT idx) const
    { return (idx < countDiffs()); }
    
    /// Fetches a diffusion rule by its idx.
    DiffDef * diff(gidxT idx) const
    { return (isValidDiff(idx) ? pDiffs[idx] : 0); }
    
    ///
    std::vector<DiffDef*>::const_iterator beginDiff(void) const
    { return pDiffs.begin(); }
    
    ///
    std::vector<DiffDef*>::const_iterator endDiff(void) const
    { return pDiffs.end(); }
    
    ////////////////////////////////////////////////////////////////////////
    
    PatchDef * createPatchDef(std::string const & name, 
        CompDef * inner, CompDef * outer);
    
    uint countPatches(void) const
    { return pPatches.size(); }
    
    bool isValidPatch(gidxT idx) const
    { return (idx < countPatches()); }
    
    PatchDef * patch(gidxT idx) const
    { return (isValidPatch(idx) ? pPatches[idx] : 0); }
    
    std::vector<PatchDef*>::const_iterator beginPatch(void) const
    { return pPatches.begin(); }
    
    std::vector<PatchDef*>::const_iterator endPatch(void) const
    { return pPatches.end(); }
    
    ////////////////////////////////////////////////////////////////////////
    
    SReacDef * createSReacDef(std::string const & name, bool inside = true);
    
    uint countSReacs(void) const
    { return pSReacs.size(); }
    
    bool isValidSReac(gidxT idx) const
    { return (idx < countSReacs()); }
    
    SReacDef * sreac(gidxT idx) const
    { return (isValidSReac(idx) ? pSReacs[idx] : 0); }
    
    std::vector<SReacDef*>::const_iterator beginSReac(void) const
    { return pSReacs.begin(); }
    
    std::vector<SReacDef*>::const_iterator endSReac(void) const
    { return pSReacs.end(); }
    
    ////////////////////////////////////////////////////////////////////////
    
private:
    
    ////////////////////////////////////////////////////////////////////////
    
    Mode                        pMode;
    bool                        pFinalSetupFinished;
    
    ////////////////////////////////////////////////////////////////////////
    
    std::vector<SpecDef*>       pSpecs;
    std::vector<CompDef*>       pComps;
    std::vector<ReacDef*>       pReacs;
    std::vector<DiffDef*>       pDiffs;
    std::vector<PatchDef*>      pPatches;
    std::vector<SReacDef*>      pSReacs;
    
    ////////////////////////////////////////////////////////////////////////
    
};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(sim)
END_NAMESPACE(steps)

#endif
// STEPS_SIM_SHARED_STATEDEF_HPP

// END
