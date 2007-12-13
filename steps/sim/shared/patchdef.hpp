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

#ifndef STEPS_SIM_SHARED_PATCHDEF_HPP
#define STEPS_SIM_SHARED_PATCHDEF_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>

// STEPS headers.
#include <steps/common.h>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(sim)

// Forward declarations.
class CompDef;
class PatchDef;
class SpecDef;
class SReacDef;
class StateDef;

// Auxiliary declarations.
typedef PatchDef *                      PatchDefP;
typedef std::vector<PatchDefP>          PatchDefPVec;
typedef PatchDefPVec::iterator          PatchDefPVecI;
typedef PatchDefPVec::const_iterator    PatchDefPVecCI;

////////////////////////////////////////////////////////////////////////////////

class PatchDef
{

public:

    /// Constructor.
    ///
    PatchDef
    (
        StateDef * sdef, gidxT idx, 
        std::string const & name, 
        CompDef * inner, CompDef * outer
    );
    
    /// Destructor.
    ///
    ~PatchDef(void);
    
    ////////////////////////////////////////////////////////////////////////
    // PATCHDEF SETUP
    ////////////////////////////////////////////////////////////////////////
    
    /// Declare that a chemical species, specified by its global index, 
    /// can be embedded on this patch. This is not required for species
    /// embedded on the in- or outside compartment.
    ///
    void addSpec(gidxT idx);
    
    /// Declare that a certain sreaction rule, specified by its global
    /// index, can occur in this compartment.
    ///
    void addSReac(gidxT idx);
    
    /// Build a set of local indices for species and sreaction rules 
    /// (called during setup, by StateDef::setupFinal()).
    ///
    void setupLocalIndices(void);
    
    /// Create update objects for all species and sreaction rules.
    ///
    /// NOTE: Called during setup, by StateDef::setupFinal(). Requires 
    /// that CompDef::setupLocalIndices() and PatchDef::setupLocalIndices
    /// has been called for all compartments and patches prior to this.
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
    
    double area(void) const
    { return pArea; }
    
    void setArea(double const & area)
    { pArea = area; }
    
    ////////////////////////////////////////////////////////////////////////

    CompDef * icompdef(void) const
    { return pInner; }
    
    CompDef * ocompdef(void) const
    { return pOuter; }
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////
    
    uint countSpecs(void) const
    { return pSpec_N; }
    
    lidxT specG2L(gidxT idx) const
    { return pSpec_G2L[idx]; }
    
    gidxT specL2G(lidxT idx) const
    { return pSpec_L2G[idx]; }
    
    /// Auxiliary function: resolves a species gidx for the outer
    /// compartment.
    ///
    /// \return The local index or steps::sim::shared::LIDX_UNDEFINED.
    ///
    lidxT specG2L_O(gidxT idx) const
    {
        if (ocompdef() == 0) return LIDX_UNDEFINED;
        return ocompdef()->specG2L(idx);
    }
    
    /// Auxiliary function: resolves a species lidx for the outer
    /// compartment.
    ///
    /// \return The global index or steps::sim::shared::GIDX_UNDEFINED.
    ///
    gidxT specL2G_O(lidxT idx) const
    {
        if (ocompdef() == 0) return GIDX_UNDEFINED;
        return ocompdef()->specL2G(idx);
    }
    
    /// Auxiliary function: resolves a species gidx for the inner
    /// compartment.
    ///
    /// \return The local index or steps::sim::shared::LIDX_UNDEFINED.
    ///
    lidxT specG2L_I(gidxT idx) const
    {
        if (icompdef() == 0) return LIDX_UNDEFINED;
        return icompdef()->specG2L(idx);
    }
    
    /// Auxiliary function: resolves a species lidx for the inner
    /// compartment.
    ///
    /// \return The local index or steps::sim::shared::GIDX_UNDEFINED.
    ///
    gidxT specL2G_I(lidxT idx) const
    {
        if (icompdef() == 0) return GIDX_UNDEFINED;
        return icompdef()->specL2G(idx);
    }
    
    /// Return a species definition (a pointer to an object of type 
    /// SpecDef) by its local index.
    ///
    SpecDef * spec(lidxT idx) const;
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE REACTION RULES
    ////////////////////////////////////////////////////////////////////////
    
    uint countSReacs(void) const
    { return pSReac_N; }
    
    lidxT sreacG2L(gidxT idx) const
    { return pSReac_G2L[idx]; }
    
    gidxT sreacL2G(lidxT idx) const
    { return pSReac_L2G[idx]; }
    
    /// Return a surface reaction definition (a pointer to an object of 
    /// type SReacDef) by its local index.
    ///
    SReacDef * sreac(lidxT idx) const;
    
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
    double                      pArea;
    
    /// Pointer to inner compartment CompDef.
    CompDef *                   pInner;
    /// Pointer to outer compartment CompDef.
    CompDef *                   pOuter;
    
    ////////////////////////////////////////////////////////////////////////
    // DATA: SPECIES
    ////////////////////////////////////////////////////////////////////////
    
    /// Number of species embedded in patch.
    uint                        pSpec_N;
    /// Table to resolve species index (global -> local).
    lidxT *                     pSpec_G2L;
    /// Table to resolve species index (local -> global).
    gidxT *                     pSpec_L2G;
    
    ////////////////////////////////////////////////////////////////////////
    // DATA: SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////
    
    /// Number of surface reactions occuring in patch.
    uint                        pSReac_N;
    /// Table to resolve reaction rule indices (global -> local).
    lidxT *                     pSReac_G2L;
    /// Table to resolve reaction rule indices (local -> global).
    gidxT *                     pSReac_L2G;

    depT *                      pSReac_DEP_I_Spec;
    depT *                      pSReac_DEP_S_Spec;
    depT *                      pSReac_DEP_O_Spec;
    uint *                      pSReac_LHS_I_Spec;
    uint *                      pSReac_LHS_S_Spec;
    uint *                      pSReac_LHS_O_Spec;
    int *                       pSReac_UPD_I_Spec;
    int *                       pSReac_UPD_S_Spec;
    int *                       pSReac_UPD_O_Spec;
    
    ////////////////////////////////////////////////////////////////////////
    
};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(sim)
END_NAMESPACE(steps)

#endif
// STEPS_SIM_SHARED_PATCHDEF_HPP

// END
