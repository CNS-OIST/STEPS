////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
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
#include <set>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class DiffDef;
class ReacDef;
class SpecDef;
class StateDef;

////////////////////////////////////////////////////////////////////////////////

/// Describes the compartment-associated processes and rules that need to
/// be updated because of some event. All 'updatees' are therefore described
/// using local indices.
///
class CompUpd
{

public:

    /// Default constructor.
    ///
    CompUpd(void);
    
    /// Copy constructor.
    ///
    CompUpd(CompUpd const & c);
    
    /// Destructor.
    ///
    ~CompUpd(void);

    /// Remove duplicates and order all indices.
    ///
    void compact(void);

    /// This important operation allows a CompUpd object to be build
    /// incrementally.
    ///
    void merge(CompUpd const & upd);
    
    ////////////////////////////////////////////////////////////////////////
    
    std::vector<uint>::const_iterator beginLReacs(void) const
    { return pLReacs.begin(); }
    
    std::vector<uint>::const_iterator endLReacs(void) const
    { return pLReacs.end(); }

    std::vector<uint>::const_iterator beginLDiffs(void) const
    { return pLDiffs.begin(); }
    
    std::vector<uint>::const_iterator endLDiffs(void) const
    { return pLDiffs.end(); }

    ////////////////////////////////////////////////////////////////////////

private:

    /// A list of local reactions that need to be updated.
    /// These are local reaction indices!
    ///
    std::vector<uint>           pLReacs;
    
    /// A list of local diffusion rules that need to be updated.
    /// These are local diffusion indices!
    ///
    std::vector<uint>           pLDiffs;

    ////////////////////////////////////////////////////////////////////////

    /// Used intimately by class CompDef.
    ///
    friend class CompDef;

};

////////////////////////////////////////////////////////////////////////////////

///
class CompDef
{

public:

    /// Constructor.
    ///
    CompDef(StateDef * sdef, uint gidx, std::string const & name);
    
    // Copy constructor (object manages pointers).
    //
    // CompDef(CompDef const & c);
    
    /// Destructor.
    ///
    ~CompDef(void);
    
    ////////////////////////////////////////////////////////////////////////
    
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
    
    StateDef * statedef(void) const
    { return pStateDef; }
    
    uint gidx(void) const
    { return pGIDX; }
    
    std::string const & name(void) const
    { return pName; }
    
    ////////////////////////////////////////////////////////////////////////
    
    double vol(void) const
    { return pVolume; }
    
    void setVol(double const & vol)
    { pVolume = vol; }
    
    ////////////////////////////////////////////////////////////////////////

    /// Declare that a species, specified by its global index, is being 
    /// used by processes in this compartment.
    ///
    void addSpec(uint gidx);
    
    /// It doesn't make sense to call this method, until after method
    /// CompDef::setupLocalIndices() has been called.
    ///
    uint countSpecs(void) const
    { return pL2GSpec.size(); }
    
    /// It doesn't make sense to call this method, until after method
    /// CompDef::setupLocalIndices() has been called.
    ///
    uint specG2L(uint gidx) const
    { return pG2LSpec[gidx]; }
    
    /// It doesn't make sense to call this method, until after method
    /// CompDef::setupLocalIndices() has been called.
    ///
    uint specL2G(uint lidx) const
    { return pL2GSpec[lidx]; }
    
    /// Return a species definition (a pointer to an object of type 
    /// SpecDef) by its local index.
    ///
    SpecDef * spec(uint lidx) const;
    
    /// Return a list of things that need to be updated when the
    /// local concentration of a species has changed.
    ///
    CompUpd * updateSpec(uint lidx) const
    { return & pSpecUpd[lidx]; }
    
    ////////////////////////////////////////////////////////////////////////
    
    /// Declare that a certain reaction rule, specified by its global
    /// index, can occur in this compartment.
    ///
    /// Currently, this method automatically adds all species involved
    /// in the reaction rule to the list of species that can occur in
    /// this compartment (the alternative would be to report an error).
    ///
    void addReac(uint gidx);
    
    /// It doesn't make sense to call this method, until after method
    /// CompDef::setupLocalIndices() has been called.
    ///
    uint countReacs(void) const
    { return pL2GReac.size(); }
    
    /// It doesn't make sense to call this method, until after method
    /// CompDef::setupLocalIndices() has been called.
    ///
    uint reacG2L(uint gidx) const
    { return pG2LReac[gidx]; }
    
    /// It doesn't make sense to call this method, until after method
    /// CompDef::setupLocalIndices() has been called.
    ///
    uint reacL2G(uint lidx) const
    { return pL2GReac[lidx]; }
    
    /// Return a reaction definition (a pointer to an object of type 
    /// ReacDef) by its local index.
    ///
    ReacDef * reac(uint lidx) const;
    
    /// For a given reaction (specified by its compartment-specific index)
    /// return its stoichiometric dependency vector (= left hand side of
    /// the reaction rule). 
    ///
    /// This is returned as an array of unsigned integers; the length of 
    /// this array is determined by the number of species that can occur 
    /// in this compartment.
    ///
    uint * reacSpecDeps(uint lidx) const
    { return pReacSpecDeps + (lidx * countSpecs()); }
    
    /// For a given reaction (specified by its compartment-specific index)
    /// return its stoichiometric update vector (= {right hand side} minus
    /// {left hand side} of the reaction rule). 
    /// 
    /// This is returned as an array of signed integers; the length of 
    /// this array is determined by the number of species that can occur in 
    /// this compartment.
    ///
    int * reacSpecUpds(uint lidx) const
    { return pReacSpecUpds + (lidx * countSpecs()); }
    
    /// Return a list of things that need to be updated when a compartment
    /// reaction has occured.
    ///
    CompUpd * updateReac(uint lidx) const
    { return & pReacUpd[lidx]; }
    
    ////////////////////////////////////////////////////////////////////////
    
    /// Declare that a diffusion rule, specified by its global index, 
    /// can occur in this compartment.
    ///
    /// Currently, this method adds the ligand to the list of species
    /// that can occur in this compartment (the alternative would be
    /// to report an error saying that the ligand was not explicitly
    /// added to the compartment).
    ///
    void addDiff(uint gidx);
    
    /// It doesn't make sense to call this method, until after method
    /// CompDef::setupLocalIndices() has been called.
    ///
    uint countDiffs(void) const
    { return pL2GDiff.size(); }
    
    /// It doesn't make sense to call this method, until after method
    /// CompDef::setupLocalIndices() has been called.
    ///
    uint diffG2L(uint gidx) const
    { return pG2LDiff[gidx]; }
    
    /// It doesn't make sense to call this method, until after method
    /// CompDef::setupLocalIndices() has been called.
    ///
    uint diffL2G(uint lidx) const
    { return pL2GDiff[lidx]; }
    
    /// Return a diffusion definition (a pointer to an object of type 
    /// DiffDef) by its local index.
    ///
    DiffDef * diff(uint lidx) const;
    
    /// Return a list of things that need to be updated when a "diffusion
    /// event" (e.g. in the form of a jump) has occured. The diffusion 
    /// rule is specified by its local index.
    ///
    CompUpd * updateDiff(uint lidx) const
    { return & pDiffUpd[lidx]; }
    
    ////////////////////////////////////////////////////////////////////////
    
private:
    
    /// A pointer to the state definition.
    ///
    StateDef *                  pStateDef;
    
    /// The index of the compartment.
    ///
    uint                        pGIDX;
    
    /// The name of the compartment.
    ///
    std::string                 pName;
        
    /// The volume of the compartment.
    ///
    double                      pVolume;
    
    ////////////////////////////////////////////////////////////////////////
    
    /// Table to resolve species index (global -> local).
    ///
    std::vector<uint>           pG2LSpec;
    
    /// Table to resolve species index (local -> global).
    ///
    std::vector<uint>           pL2GSpec;
    
    /// The local processes and rules that need to be updated when a 
    /// given species is changed.
    ///
    CompUpd *                   pSpecUpd;
    
    ////////////////////////////////////////////////////////////////////////
    
    /// Table to resolve reaction rule indices (global -> local).
    ///
    std::vector<uint>           pG2LReac;
    
    /// Table to resolve reaction rule indices (local -> global).
    ///
    std::vector<uint>           pL2GReac;
    
    /// For each reaction in the compartment, its stoichiometric LHS
    /// vector (based on local indices).
    ///
    uint *                      pReacSpecDeps;
    
    /// For each reaction in the compartment, its stoichiometric 
    /// 'update vector' (based on local indices).
    ///
    int *                       pReacSpecUpds;
    
    /// The local processes and rules that need to be updated when a
    /// given reaction rule is executed.
    ///
    CompUpd *                   pReacUpd;
    
    ////////////////////////////////////////////////////////////////////////
    
    /// Table to resolve diffusion rule indices (global -> local).
    ///
    std::vector<uint>           pG2LDiff;
    
    /// Table to resolve diffusion rule indices (local -> global).
    ///
    std::vector<uint>           pL2GDiff;
    
    /// The local processes and rules that need to be updated when a 
    /// given diffusion event is executed.
    ///
    CompUpd *                   pDiffUpd;
    
};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_SHARED_COMPDEF_HPP

// END
