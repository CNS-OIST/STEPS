////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SIM_SHARED_COMPDEF_HPP
#define STEPS_SIM_SHARED_COMPDEF_HPP 1

// STL headers.
#include <set>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class ReacDef;
class SpecDef;
class StateDef;

////////////////////////////////////////////////////////////////////////////////

///
class CompUpd
{

public:

    CompUpd(void);
    ~CompUpd(void);
    
    void merge(CompUpd const & upd);
    
    ////////////////////////////////////////////////////////////////////////
    
    std::vector<uint>::const_iterator beginLReacs(void) const
    { return pLReacs.begin(); }
    
    std::vector<uint>::const_iterator endLReacs(void) const
    { return pLReacs.end(); }

private:

    /// Mind, these are local reaction indices!
    std::vector<uint>           pLReacs;

    ////////////////////////////////////////////////////////////////////////

    friend class CompDef;

};

////////////////////////////////////////////////////////////////////////////////

///
class CompDef
{

public:

    /// Constructor.
    CompDef(StateDef * sdef, uint gidx, std::string const & name);
    /// Destructor.
    ~CompDef(void);
    
    StateDef * statedef(void) const
    { return pStateDef; }
    
    uint gidx(void) const
    { return pGIDX; }
    
    std::string const & name(void) const
    { return pName; }
    
    ////////////////////////////////////////////////////////////////////////
    
    /// Gets called when the definition of all components in the entire state
    /// has finished.
    void finalSetup(void);

    ////////////////////////////////////////////////////////////////////////

    ///
    void addSpec(uint gidx);
    
    uint countSpecs(void) const
    { return pL2GSpec.size(); }
    
    uint specG2L(uint gidx) const
    { return pG2LSpec[gidx]; }
    
    uint specL2G(uint lidx) const
    { return pL2GSpec[lidx]; }
    
    SpecDef * spec(uint lidx) const;
    
    /// Return a list of things that need to be updated when the
    /// local concentration of a species has changed.
    CompUpd * updateSpec(uint lidx) const
    { return & pSpecUpd[lidx]; }
    
    ////////////////////////////////////////////////////////////////////////
    
    /// 
    void addReac(uint gidx);
    
    uint countReacs(void) const
    { return pL2GReac.size(); }
    
    uint reacG2L(uint gidx) const
    { return pG2LReac[gidx]; }
    
    uint reacL2G(uint lidx) const
    { return pL2GReac[lidx]; }
    
    ReacDef * reac(uint lidx) const;
    
    uint * reacSpecDeps(uint lidx) const
    { return pReacSpecDeps + (lidx * countSpecs()); }
    
    int * reacSpecUpds(uint lidx) const
    { return pReacSpecUpds + (lidx * countSpecs()); }
    
    /// Return a list of things that need to be updated when a compartment
    /// reaction has occured.
    ///
    CompUpd * updateReac(uint lidx) const
    { return & pReacUpd[lidx]; }
    
    ////////////////////////////////////////////////////////////////////////
    
private:
    
    /// A pointer to the state definition.
    StateDef *                  pStateDef;
    
    /// The index of the compartment.
    uint                        pGIDX;
    
    /// The name of the compartment.
    std::string                 pName;
        
    /// The volume of the compartment.
    double                      pVolume;
        
    ///
    std::vector<uint>           pG2LSpec;
    ///
    std::vector<uint>           pL2GSpec;
    ///
    CompUpd *                   pSpecUpd;
    
    ///
    std::vector<uint>           pG2LReac;
    ///
    std::vector<uint>           pL2GReac;
    ///
    uint *                      pReacSpecDeps;
    ///
    int *                       pReacSpecUpds;
    ///
    CompUpd *                   pReacUpd;
    
};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_SHARED_COMPDEF_HPP

// END
