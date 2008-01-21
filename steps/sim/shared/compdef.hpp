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

#ifndef STEPS_SIM_SHARED_COMPDEF_HPP
#define STEPS_SIM_SHARED_COMPDEF_HPP 1

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
class ReacDef;
class SpecDef;
class StateDef;

// Auxiliary declarations.
typedef CompDef *                       CompDefP;
typedef std::vector<CompDefP>           CompDefPVec;
typedef CompDefPVec::iterator           CompDefPVecI;
typedef CompDefPVec::const_iterator     CompDefPVecCI;

////////////////////////////////////////////////////////////////////////////////

class CompDef
{

public:

    /// Constructor.
    ///
    CompDef(StateDef * sdef, gidxT idx, std::string const & name);
    
    /// Destructor.
    ///
    ~CompDef(void);
    
    ////////////////////////////////////////////////////////////////////////
    // COMPDEF: SETUP
    ////////////////////////////////////////////////////////////////////////
    
    /// Declare that a species, specified by its global index, can occur
    /// in this compartment.
    ///
    void addSpec(gidxT idx);
    
    /// Declare that a certain reaction rule, specified by its global
    /// index, can occur in this compartment.
    ///
    /// Currently, this method automatically adds all species involved
    /// in the reaction rule to the list of species that can occur in
    /// this compartment (the alternative would be to report an error).
    ///
    void addReac(gidxT idx);
    
    /// Declare that a diffusion rule, specified by its global index, 
    /// can occur in this compartment.
    ///
    /// Currently, this method adds the ligand to the list of species
    /// that can occur in this compartment (the alternative would be
    /// to report an error saying that the ligand was not explicitly
    /// added to the compartment).
    ///
    void addDiff(gidxT idx);
    
    /// Build a set of local indices for species, reaction and 
    /// diffusion rules (called during setup, by StateDef::setupFinal()).
    ///
    void setupLocalIndices(void);
    
    /// Create CompDef objects for all species, reaction and diffusion
    /// rules (called during setup, by StateDef::setupFinal(); requires
    /// that CompDef::setupLocalIndices() has been called prior to this).
    ///
    void setupDependencies(void);
    
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
    
    double vol(void) const
    { return pVolume; }
    
    void setVol(double const & vol)
    { pVolume = vol; }
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////
    
    uint countSpecs(void) const
    { return pSpec_N; }
    
    lidxT specG2L(gidxT idx) const
    { return pSpec_G2L[idx]; }
    
    gidxT specL2G(lidxT idx) const
    { return pSpec_L2G[idx]; }
    
    /// Return a species definition (a pointer to an object of type 
    /// SpecDef) by its local index.
    ///
    SpecDef * spec(lidxT idx) const;
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: REACTION RULES
    ////////////////////////////////////////////////////////////////////////
    
    uint countReacs(void) const
    { return pReac_N; }
    
    lidxT reacG2L(gidxT idx) const
    { return pReac_G2L[idx]; }
    
    gidxT reacL2G(lidxT idx) const
    { return pReac_L2G[idx]; }
    
    depT reac_dep(lidxT reac, lidxT spec) const;
    uint * reac_lhs_bgn(lidxT reac) const;
    uint * reac_lhs_end(lidxT reac) const;
    int * reac_upd_bgn(lidxT reac) const;
    int * reac_upd_end(lidxT reac) const;
    
    /// Return a reaction definition (a pointer to an object of type 
    /// ReacDef) by its local index.
    ///
    ReacDef * reac(lidxT idx) const;
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION RULES
    ////////////////////////////////////////////////////////////////////////
    
    uint countDiffs(void) const
    { return pDiff_N; }
    
    lidxT diffG2L(gidxT idx) const
    { return pDiff_G2L[idx]; }
    
    gidxT diffL2G(lidxT idx) const
    { return pDiff_L2G[idx]; }
    
    depT diff_dep(lidxT diff, lidxT spec) const;
    lidxT diff_lig(lidxT diff) const;
    
    /// Return a diffusion definition (a pointer to an object of type 
    /// DiffDef) by its local index.
    ///
    DiffDef * diff(lidxT idx) const;
    
    ////////////////////////////////////////////////////////////////////////
    
private:
    
    ////////////////////////////////////////////////////////////////////////
    // DATA: GENERAL
    ////////////////////////////////////////////////////////////////////////
    
    /// A pointer to the state definition.
    StateDef *                  pStateDef;
    
    /// The index of the compartment.
    gidxT                       pGIDX;
    /// The name of the compartment.
    std::string                 pName;
        
    /// The volume of the compartment.
    double                      pVolume;
    
    ////////////////////////////////////////////////////////////////////////
    // DATA: SPECIES
    ////////////////////////////////////////////////////////////////////////
    
    /// Number of species embedded in compartment.
    uint                        pSpec_N;
    /// Table to resolve species index (global -> local).
    lidxT *                     pSpec_G2L;
    /// Table to resolve species index (local -> global).
    gidxT *                     pSpec_L2G;
    
    ////////////////////////////////////////////////////////////////////////
    // DATA: REACTION RULES
    ////////////////////////////////////////////////////////////////////////
    
    /// Number of reaction rules occuring in compartment.
    uint                        pReac_N;
    /// Table to resolve reaction rule indices (global -> local).
    lidxT *                     pReac_G2L;
    /// Table to resolve reaction rule indices (local -> global).
    gidxT *                     pReac_L2G;

    inline uint _IDX_Reac_Spec(lidxT reac) const
    { return countSpecs() * reac; }
    inline uint _IDX_Reac_Spec(lidxT reac, lidxT spec) const
    { return (countSpecs() * reac) + spec; }
        
    depT *                      pReac_DEP_Spec;
    uint *                      pReac_LHS_Spec;
    int *                       pReac_UPD_Spec;
    
    ////////////////////////////////////////////////////////////////////////
    // DATA: DIFFUSION RULES
    ////////////////////////////////////////////////////////////////////////
    
    /// Number of diffusion rules occuring in compartment.
    uint                        pDiff_N;
    /// Table to resolve diffusion rule indices (global -> local).
    lidxT *                     pDiff_G2L;
    /// Table to resolve diffusion rule indices (local -> global).
    gidxT *                     pDiff_L2G;
    
    inline uint _IDX_Diff_Spec(lidxT diff) const
    { return countSpecs() * diff; }
    inline uint _IDX_Diff_Spec(lidxT diff, lidxT spec) const
    { return (countSpecs() * diff) + spec; }
    
    depT *                      pDiff_DEP_Spec;
    lidxT *                     pDiff_LIG;
    
    ////////////////////////////////////////////////////////////////////////
    
};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(sim)
END_NAMESPACE(steps)

#endif
// STEPS_SIM_SHARED_COMPDEF_HPP

// END
