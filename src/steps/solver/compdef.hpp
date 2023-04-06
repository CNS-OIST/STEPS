/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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


#ifndef STEPS_SOLVER_COMPDEF_HPP
#define STEPS_SOLVER_COMPDEF_HPP 1


// STL headers.
#include <string>
#include <fstream>
#include <set>

// STEPS headers.
#include "util/common.h"
#include "api.hpp"
#include "diffdef.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace solver{

// Forwards declarations
class Statedef;
class Reacdef;
class Patchdef;
class Compdef;

// Auxiliary declarations.
typedef Compdef *                       CompDefP;
typedef std::vector<CompDefP>           CompDefPVec;
typedef CompDefPVec::iterator           CompDefPVecI;
typedef CompDefPVec::const_iterator     CompDefPVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Compdef object defines a compartment object.
class Compdef
{

public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the compartment.
    /// \param c Associated Comp object.
    Compdef(Statedef * sd, uint idx, steps::wm::Comp * c);

    /// Destructor
    ~Compdef();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);


    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

    /// Return the volume of this compartment.
    inline double vol() const noexcept
    { return pVol; }

    /// Return the global index of this compartment
    inline uint gidx() const noexcept
    { return pIdx; }

    /// Return the name of this compartment.
    inline const std::string& name() const noexcept
    { return pName; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup all reference.
    void setup_references();

    /// Setup all indices.
    void setup_indices();

    /// Add a species with global index gidx.
    ///
    /// \param gidx Index of the species.
    void addSpec(uint gidx);

    /// Add an inner patch.
    ///
    /// \param p Pointer to the inner patch.
    void addIPatchdef(steps::solver::Patchdef * p);

    /// Add an outer patch.
    ///
    /// \param p Pointer to the outer patch.
    void addOPatchdef(steps::solver::Patchdef * p);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

    /// Set the volume of the compartment.
    ///
    /// \param v Volume of the compartment.
    void setVol(double v);

    /// Reset count, flags members of this compartment. Called when reset()
    /// method in solver object is executed.
    void reset();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of species that can occur in this compartment
    inline uint countSpecs() const noexcept
    { return pSpecsN; }

    /// Return the local species index for global index argument.
    ///
    /// \param gidx Global index of the species.
    inline uint specG2L(uint gidx) const noexcept
    { return pSpec_G2L[gidx]; }

    /// Return the global species index for local index argument.
    ///
    /// \param lidx local index of the species.
    inline uint specL2G(uint lidx) const noexcept
    { return pSpec_L2G[lidx]; }

    /// Returns pointer to flags on species for this compartment.
    inline uint * flags() const noexcept
    { return pPoolFlags; }

    static const uint CLAMPED = 1;

    /// Return whether a species, specified by local index argument, is
    /// clamped or not.
    ///
    /// \param sildx Local index of the species.
    inline bool clamped(uint slidx) const noexcept
    { return pPoolFlags[slidx] & CLAMPED; }

    /// Return pointer to species' counts in this compartment.
    inline double * pools() const noexcept
    { return pPoolCount; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Set the species count of species specified by local index argument.
    ///
    /// \param slidx Local index of the species.
    /// \param count Count of molecules of a species.
    void setCount(uint slidx, double count);

    /// Clamp or unclamp species specified by local index argument
    ///
    /// \param slidx Local index of the species.
    /// \param clamp Flag to set if clamping or not.
    void setClamped(uint slidx, bool clamp);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of reactions that can occur in this compartment.
    inline uint countReacs() const noexcept
    { return pReacsN; }

    /// Return the local reaction index for global index argument.
    ///
    /// \param gidx Global index of the reaction.
    uint reacG2L(uint gidx) const noexcept
    { return pReac_G2L[gidx]; }

    /// Return the beginning of the lhs array of reaction specified by local
    /// index argument.
    ///
    /// \param rlidx Local index of the reaction.
    uint * reac_lhs_bgn(uint rlidx) const;


    /// Return the end of the lhs array of reaction specified by local
    /// index argument.
    ///
    /// \param rlidx Local index of the reaction.
    uint * reac_lhs_end(uint rlidx) const;

    /// Return the beginning of the update array of reaction specified by
    /// local index argument.
    ///
    /// \param rlidx Local index of the reaction.
    int * reac_upd_bgn(uint rlidx) const;

    /// Return the end of the update array of reaction specified by
    /// local index argument.
    ///
    /// \param rlidx Local index of the reaction.
    int * reac_upd_end(uint rlidx) const;

    /// Return the local index of species of reaction specified by
    /// local index argument.
    ///
    /// \param rlidx Local index of the reaction.
    /// \param slidx Local index of the species.
    /// \todo make sure this is correct.
    int reac_dep(uint rlidx, uint slidx) const;

    /// Return a pointer to reaction definition object (type Reacdef)
    /// specified by local index.
    ///
    /// \param rlidx Local index of the reaction.
    Reacdef * reacdef(uint rlidx) const;

    /// Return pointer to flags on reactions for this compartment.
    inline uint * rflags() const noexcept
    { return pReacFlags; }

    static const uint INACTIVATED = 1;

    /// Return whether a reaction, specified by local index argument, is
    /// active or not.
    ///
    /// \param rlidx Local index of the reaction.
    inline bool active(uint rlidx) const noexcept
    { return !(pReacFlags[rlidx] & INACTIVATED); }

    /// Return the kcst for a reaction specified by local index
    ///
    /// \param rlidx Local index of the reaction.
    inline double kcst(uint rlidx) const noexcept
    { return pReacKcst[rlidx]; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Set the reaction kcst for reaction specified by local index
    ///
    /// \param rlidx Local index of the reaction.
    /// \param kcst Rate constant of the reaction.
    void setKcst(uint rlidx, double kcst);

    /// Activate or inactivate a reaction specified by local index argument.
    ///
    /// \param rlidx Local index of the reaction.
    /// \param Flag to activate or inactivate a reaction.
    void setActive(uint rlidx, bool active);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of diffusion rules for this compartment.
    inline uint countDiffs() const noexcept
    { return pDiffsN; }

    /// Returns a pointer to Diffdef specified by local index.
    ///
    /// \param dlidx Local index of the difusion.
    Diffdef * diffdef(uint dlidx) const;

    /// Return the local diffusion index for global index argument.
    ///
    /// \param gidx Global index of the diffusion.
    uint diffG2L(uint gidx) const noexcept
    { return pDiff_G2L[gidx]; }

    /// Return the local index of species of diffusion specified by
    /// local index argument.
    ///
    /// \param dlidx Local index of the diffusion rule.
    /// \param slidx Local index of the species.
    /// \todo make sure this is correct.
    uint diff_dep(uint dlidx, uint slidx) const;

    /// Return the rate constant of diffusion by local index argument.
    ///
    /// \param dlidx Local index of the diffusion.
    inline double dcst(uint dlidx) const noexcept
    { return pDiffDcst[dlidx]; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Set the diffusion dcst for diffusion rule specified by local index
    ///
    /// \param dlidx Local index of the diffusion.
    /// \param dcst Rate constant of the diffusion.
    void setDcst(uint dlidx, double dcst);

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////
    // DATA: GENERAL
    ////////////////////////////////////////////////////////////////////////

    // A pointer to the state definition.
    Statedef                          * pStatedef;

    // The string identifier of the compartment
    std::string                         pName;

    // The volume of the compartment
    double                              pVol;

    // The global index of the compartment.
    uint                                pIdx;

    // The enclosed volume systems, stored as strings
    std::set<std::string>                 pCvsys;

    // Keep track of whether setup_ has been called.
    bool                                pSetupRefsdone;
    bool                                pSetupIndsdone;

    ////////////////////////////////////////////////////////////////////////
    // DATA: SPECIES
    ////////////////////////////////////////////////////////////////////////

    // Number of species that can appear in this compartment.
    uint                                pSpecsN;

    // Table to resolve species index (global -> local).
    uint                              * pSpec_G2L;

    // Table to resolve species index (local -> global).
    uint                              * pSpec_L2G;

    // Table of the populations of the species in this compartment.
    double                            * pPoolCount;

    // Table of 'clamped' flags on the species
    uint                              * pPoolFlags;

    ////////////////////////////////////////////////////////////////////////
    // DATA: REACTION RULES
    ////////////////////////////////////////////////////////////////////////

    // Number of reaction rules occuring in this compartment.
    uint                                pReacsN;

    // Table to resolve reaction rule indices (global -> local).
    uint                              * pReac_G2L;

    // Table to resolve reaction rule indices (local -> global).
    uint                              * pReac_L2G;

    // Table of the K-constants of the reaction rules in this compartment
    double                            * pReacKcst;

    // Table of 'active' flags on the reaction rules.
    uint                              * pReacFlags;

    inline uint _IDX_Reac_Spec(uint reac, uint spec) const
    { return (pSpecsN * reac) + spec; }

    int                               * pReac_DEP_Spec;
    uint                              * pReac_LHS_Spec;
    int                               * pReac_UPD_Spec;

    ////////////////////////////////////////////////////////////////////////
    // DATA: DIFFUSION RULES
    ////////////////////////////////////////////////////////////////////////

    // Number of diffusion rules occuring in compartment.
    uint                                pDiffsN;

    // Table to resolve diffusion rule indices (global -> local).
    uint                              * pDiff_G2L;

    // Table to resolve diffusion rule indices (local -> global).
    uint                              * pDiff_L2G;

    // Table of the D-constants of the diffusion rules in this compartment
    double                            * pDiffDcst;

    inline uint _IDX_Diff_Spec(uint diff, uint spec) const noexcept
    { return (pSpecsN * diff) + spec; }
    uint                              * pDiff_DEP_Spec;
    uint                              * pDiff_LIG;

    ////////////////////////////////////////////////////////////////////////
    // DATA: CONNECTIVITY
    ////////////////////////////////////////////////////////////////////////

    // A list of patches surrounding
    std::vector<Patchdef *>   pIPatches;
    std::vector<Patchdef *>   pOPatches;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_SOLVER_COMPDEF_HPP

// END
