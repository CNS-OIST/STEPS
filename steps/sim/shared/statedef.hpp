////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SIM_SHARED_STATEDEF_HPP
#define STEPS_SIM_SHARED_STATEDEF_HPP 1

// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class CompDef;
class ReacDef;
class SpecDef;

////////////////////////////////////////////////////////////////////////////////

/// This class was designed for fast lookup during simulation
/// rather than accomodating change; the latter is the task of the Python
/// objects in package steps.model.

class StateDef
{

public:
    
    enum Mode
    {
        NAIVE_MODE,
        SETUP_MODE,
        SETUP_VAR_MODE,
        SETUP_REAC_MODE,
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
    
    /// Gets called when the definition of all components in the entire state
    /// has finished.
    void finalSetup(void);
    
    ////////////////////////////////////////////////////////////////////////

    ///
    SpecDef * createSpecDef(std::string const & name);
    
    /// Returns the number of species in the overall simulation state.
    uint countSpecs(void) const
    { return pSpecs.size(); }
    
    /// Check whether the specified global index refers to a valid species.
    bool isValidSpec(uint gidx) const
    { return (gidx < countSpecs()); }
    
    /// Fetches a species by its gidx.
    SpecDef * spec(uint gidx) const
    { return (isValidSpec(gidx) ? pSpecs[gidx] : 0 ); }
    
    ///
    std::vector<SpecDef*>::const_iterator beginSpec(void) const
    { return pSpecs.begin(); }
    
    ///
    std::vector<SpecDef*>::const_iterator endSpec(void) const
    { return pSpecs.end(); }
    
    ////////////////////////////////////////////////////////////////////////

    ///
    ReacDef * createReacDef(std::string const & name);
    
    /// Returns the number of reactions in the overall simulation state.
    uint countReacs(void) const
    { return pReacs.size(); }
    
    /// Check whether the specified global index refers to a valid reaction.
    bool isValidReac(uint gidx) const
    { return (gidx < countReacs()); }
    
    /// Fetches a reaction by its gidx.
    ReacDef * reac(uint gidx) const
    { return (isValidReac(gidx) ? pReacs[gidx] : 0 ); }
    
    ///
    std::vector<ReacDef*>::const_iterator beginReac(void) const
    { return pReacs.begin(); }
    
    ///
    std::vector<ReacDef*>::const_iterator endReac(void) const
    { return pReacs.end(); }
    
    ////////////////////////////////////////////////////////////////////////

    /// 
    CompDef * createCompDef(std::string const & name);
    
    /// Returns the number of compartments in the simulation state.
    uint countComps(void) const
    { return pComps.size(); }
    
    /// Check whether the specified global index refers to a valid compartment.
    bool isValidComp(uint gidx) const
    { return (gidx < countComps()); }
    
    /// Fetches a compartment by its gidx.
    CompDef * comp(uint gidx) const
    { return (isValidComp(gidx) ? pComps[gidx] : 0); }

    ///
    std::vector<CompDef*>::const_iterator beginComp(void) const
    { return pComps.begin(); }
    
    ///
    std::vector<CompDef*>::const_iterator endComp(void) const
    { return pComps.end(); }
    
private:
    
    Mode                        pMode;
    bool                        pFinalSetupFinished;
    
    /// A list of species.
    std::vector<SpecDef*>       pSpecs;
    
    /// A list of reactions.
    std::vector<ReacDef*>       pReacs;
           
    /// A list of compartments.
    std::vector<CompDef*>       pComps;
    
};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_SHARED_STATEDEF_HPP

// END
