////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SIM_SHARED_DIFFDEF_HPP
#define STEPS_SIM_SHARED_DIFFDEF_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>

////////////////////////////////////////////////////////////////////////////////

class StateDef;

class DiffDef
{

public:

    /// Constructor. 
    ///
    DiffDef(StateDef * sdef, uint gidx, std::string const & name);

    StateDef * statedef(void) const
    { return pStateDef; }
    
    uint gidx(void) const
    { return pGIDX; }
    
    std::string const & name(void) const
    { return pName; }

    ////////////////////////////////////////////////////////////////////////
    
    /// Gets called when the definition of all components in the entire state
    /// has finished.
    ///
    /// Currently, this method doesn't really have to do anything.
    ///
    void setupFinal(void);

    ////////////////////////////////////////////////////////////////////////

    /// Returns the default diffusion constant of this diffusion rule.
    ///
    double dcst(void) const
    { return pDcst; }
    
    /// Set the default diffusion constant for this diffusion rule.
    ///
    void setDcst(double const & d);

    ////////////////////////////////////////////////////////////////////////
    
    /// Return the global index of the ligand species.
    ///
    uint lig(void) const
    { return pLig; }

    /// Set the ligand of a diffusion rule by its global index.
    ///
    void setLig(uint gidx);

private:
    
    StateDef *                  pStateDef;
    
    /// Keeps track of whether method setupFinal() has already been called.
    ///
    bool                        pFinalSetupFinished;
    
    uint                        pGIDX;
    
    std::string                 pName;

    /// Default (MACROscopic) diffusion constant.
    ///
    double                      pDcst;
    
    uint                        pLig;

};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_SHARED_DIFFDEF_HPP

// END
