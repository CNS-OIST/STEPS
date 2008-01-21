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
#include <steps/sim/shared/types.hpp>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(sim)

// Forward declarations.
class DiffDef;
class StateDef;

// Auxiliary declarations.
typedef DiffDef *                       DiffDefP;
typedef std::vector<DiffDefP>           DiffDefPVec;
typedef DiffDefPVec::iterator           DiffDefPVecI;
typedef DiffDefPVec::const_iterator     DiffDefPVecCI;

////////////////////////////////////////////////////////////////////////////////

class DiffDef
{

public:

    ////////////////////////////////////////////////////////////////////////
    
    /// Constructor. 
    ///
    DiffDef(StateDef * sdef, gidxT idx, std::string const & name);
    
    /// Destructor
    ///
    ~DiffDef(void);
    
    ////////////////////////////////////////////////////////////////////////
    // DIFFDEF SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Set the ligand of a diffusion rule by its global index.
    ///
    void setLig(gidxT idx);
    
    /// Gets called when the definition of all components in the entire state
    /// has finished.
    ///
    /// Currently, this method doesn't really have to do anything.
    ///
    void setupFinal(void);
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: GENERAL
    ////////////////////////////////////////////////////////////////////////
    
    StateDef * statedef(void) const
    { return pStateDef; }
    
    gidxT gidx(void) const
    { return pGIDX; }
    
    std::string const & name(void) const
    { return pName; }

    ////////////////////////////////////////////////////////////////////////

    /// Returns the default diffusion constant of this diffusion rule.
    ///
    double dcst(void) const
    { return pDcst; }
    
    /// Set the default diffusion constant for this diffusion rule.
    ///
    void setDcst(double const & d);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: LIGAND
    ////////////////////////////////////////////////////////////////////////
    
    /// Return the global index of the ligand species.
    ///
    inline gidxT lig(void) const
    { return pSpec_LIG; }
    depT dep(gidxT idx) const;
    bool req(gidxT idx) const;

    ////////////////////////////////////////////////////////////////////////
    
    /*
    /// Check whether occurence of a diffusion rule <EM>depends</EM> on 
    /// local changes to the concentration of some species, specified by 
    /// its global index (gidx).
    ///
    /// Currently, this only has to check whether the specified gidx equals
    /// the ligand for which the diffusion rule is defined. As we might 
    /// have concentration-dependent diffusion coefficients in the future 
    /// (or other fancy extensions to the basic diffusion rule idea), this 
    /// function will have to check for more things.
    ///
    bool dependsOnSpec(uint gidx) const;
    
    /// Check whether the occurence of the diffusion rule <EM>affects</EM>
    /// affects the concentration of some species, specified by its global 
    /// index (gidx).
    ///
    /// Currently, this only has to check whether the specified gidx equals
    /// the ligand for which the diffusion rule is defined.
    ///
    bool affectsSpec(uint gidx) const;
    */
    
    ////////////////////////////////////////////////////////////////////////
    
private:
    
    ////////////////////////////////////////////////////////////////////////
    // DATA: GENERAL
    ////////////////////////////////////////////////////////////////////////
    
    StateDef *                  pStateDef;
    
    /// Keeps track of whether method setupFinal() has already been called.
    ///
    bool                        pFinalSetupDone;
    
    gidxT                       pGIDX;
    
    std::string                 pName;

    /// Default (MACROscopic) diffusion constant.
    ///
    double                      pDcst;
    
    ////////////////////////////////////////////////////////////////////////
    // DATA: LIGAND
    ////////////////////////////////////////////////////////////////////////
    
    depT *                      pSpec_DEP;
    gidxT                       pSpec_LIG;

    ////////////////////////////////////////////////////////////////////////
    
};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(sim)
END_NAMESPACE(steps)

#endif
// STEPS_SIM_SHARED_DIFFDEF_HPP

// END
