/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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


/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

#ifndef STEPS_SOLVER_PATCHDEF_HPP
#define STEPS_SOLVER_PATCHDEF_HPP 1


// STL headers.
#include <string>
#include <fstream>

// STEPS headers.
#include "util/common.h"
#include "statedef.hpp"
#include "api.hpp"
#include "geom/patch.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace solver{

////////////////////////////////////////////////////////////////////////////////

// Forwards declarations
class Statedef;
class SReacdef;
class Diffdef;
class Compdef;
class VDepTransdef;
class VDepSReacdef;
class OhmicCurrdef;
class GHKcurrdef;

// Auxiliary declarations.
typedef Patchdef *                      PatchDefP;
typedef std::vector<PatchDefP>          PatchDefPVec;
typedef PatchDefPVec::iterator          PatchDefPVecI;
typedef PatchDefPVec::const_iterator    PatchDefPVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Defined patch object.
class Patchdef
{

public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the patch.
    /// \param p Pointer to the Patch object.
    Patchdef(Statedef * sd, uint idx, steps::wm::Patch * p);

    /// Destructor
    ~Patchdef();
    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: PATCH
    ////////////////////////////////////////////////////////////////////////

    /// Return the area of this patch.
    double area() const;

    /// Return the global index of this patch.
    inline uint gidx() const noexcept
    { return pIdx; }

    /// Return the name of the patch.
    const std::string& name() const noexcept
    { return pName; }

    /// Return a pointer to the inner compartment.
    inline Compdef * icompdef() const noexcept
    { return pInner; }

    /// Return a pointer to the outer compartment.
    inline Compdef * ocompdef() const noexcept
    { return pOuter; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup all references.
    void setup_references();

    /// Setup all indices.
    void setup_indices();

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: PATCH
    ////////////////////////////////////////////////////////////////////////

    /// Set the area of this patch.
    ///
    /// \param a Area of the patch.
    void setArea(double a);

    /// Get the area of the patch
    ///
    /// \return Area of the patch
    inline double getArea() const noexcept
    { return pArea; }

    /// Reset count, flags members of this patch. Called when reset()
    /// method in solver object is executed.
    ///
    void reset();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the number of species defined for this surface patch.
    inline uint countSpecs() const noexcept
    { return pSpecsN_S; }

    /// Return the number of species defined for the inner compartment.
    /// Should not be called before Compdef::setup()
    inline uint countSpecs_I() const noexcept
    { return pSpecsN_I; }

    /// Return the number of species defined for the outer compartment.
    /// Should not be called before Compdef::setup()
    inline uint countSpecs_O() const noexcept
    { return pSpecsN_O; }

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


    /// Auxiliary function: resolves a species gidx for the inner
    /// compartment.
    ///
    /// \param gidx Global index of the species.
    /// \return The local index or steps::sim::shared::LIDX_UNDEFINED.
    ///
    uint specG2L_I(uint gidx) const;
    /// Auxiliary function: resolves a species gidx for the outer
    /// compartment.
    ///
    /// \param gidx Global index of the species.
    /// \return The local index or steps::sim::shared::LIDX_UNDEFINED.
    ///
    uint specG2L_O(uint gidx) const;

    /// Return pointer to species' counts on this patch.
    inline double * pools() const noexcept
    { return pPoolCount; }

    /// Returns pointer to flags on species for this patch.
    inline uint * flags() const noexcept
    { return pPoolFlags; }

    static const uint CLAMPED = 1;

    /// Return whether a species, specified by local index argument, is
    /// clamped or not.
    ///
    /// \param slidx Local index of the species.
    inline bool clamped(uint slidx) const noexcept
    { return pPoolFlags[slidx] & CLAMPED; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Set the species count of species specified by local index argument.
    ///
    /// \param slidx Local index of the species.
    /// \param count Count of species.
    void setCount(uint slidx, double count);

    /// Clamp or unclamp species specified by local index argument
    ///
    /// \param slidx Local index of the species.
    /// \param clamp Flag to clamp or unclamp species.
    void setClamped(uint slidx, bool clamp);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE REACTION RULES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of surface reactions that can occur on this
    /// patch.
    inline uint countSReacs() const noexcept
    { return pSReacsN; }

    /// Return the local surface reaction index for global index argument.
    ///
    /// \param gidx Global index of the surface reaction.
    inline uint sreacG2L(uint gidx) const noexcept
    { return pSReac_G2L[gidx]; }
    /*
    inline gidxT sreacL2G(lidxT idx) const
    { return pSReac_L2G[idx]; }
    */

    /// Return a pointer to reaction definition object (type SReacdef)
    /// specified by local index.
    SReacdef * sreacdef(uint lidx) const;

    /// Warning: these methods perform no error checking!
    /// \todo imcompleted.
    int sreac_dep_I(uint srlidx, uint splidx) const;
    int sreac_dep_S(uint srlidx, uint splidx) const;
    int sreac_dep_O(uint srlidx, uint splidx) const;

    /// Warning: these methods perform no error checking!
    ///
    // Return the beginning and end of the lhs arrays of surface reaction
    // specified by local index argument.
    uint * sreac_lhs_I_bgn(uint lidx) const;
    uint * sreac_lhs_I_end(uint lidx) const;
    uint * sreac_lhs_S_bgn(uint lidx) const;
    uint * sreac_lhs_S_end(uint lidx) const;
    uint * sreac_lhs_O_bgn(uint lidx) const;
    uint * sreac_lhs_O_end(uint lidx) const;

    /// Warning: these methods perform no error checking!
    ///
    // Return the beginning and end of the update arrays of surface reaction
    // specified by local index argument.
    int * sreac_upd_I_bgn(uint lidx) const;
    int * sreac_upd_I_end(uint lidx) const;
    int * sreac_upd_S_bgn(uint lidx) const;
    int * sreac_upd_S_end(uint lidx) const;
    int * sreac_upd_O_bgn(uint lidx) const;
    int * sreac_upd_O_end(uint lidx) const;

    /// Return pointer to flags on surface reactions for this patch.
    inline uint * srflags() const noexcept
    { return pSReacFlags; }

    static const uint INACTIVATED = 1;

    /// Return whether a surface reaction, specified by local index argument,
    /// is active or not.
    ///
    /// \param rlidx Local index of the surface reaction.
    inline bool active(uint rlidx) const noexcept
    { return !(pSReacFlags[rlidx] & INACTIVATED); }

    /// Return the kcst for a surface reaction specified by local index
    ///
    /// \param rlidx Local index of the surface reaction.
    inline double kcst(uint rlidx) const noexcept
    { return pSReacKcst[rlidx]; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of surface diffusion rules for this patch.
    inline uint countSurfDiffs() const noexcept
    { return pSurfDiffsN; }

    /// Returns a pointer to Diffdef specified by local index.
    ///
    /// \param dlidx Local index of the surface diffusion.
    Diffdef * surfdiffdef(uint dlidx) const;

    /// Return the local surface diffusion index for global index argument.
    ///
    /// \param gidx Global index of the surface diffusion.
    uint surfdiffG2L(uint gidx) const noexcept
    { return pSurfDiff_G2L[gidx]; }

    /// Return the local index of species of surface diffusion specified by
    /// local index argument.
    ///
    /// \param dlidx Local index of the surface diffusion rule.
    /// \param slidx Local index of the species.
    /// \todo make sure this is correct.
    uint surfdiff_dep(uint dlidx, uint slidx) const;

    /// Return the rate constant of surface diffusion by local index argument.
    ///
    /// \param dlidx Local index of the surface diffusion.
    inline double dcst(uint dlidx) const noexcept
    { return pSurfDiffDcst[dlidx]; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: OHMIC CURRENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return the number of ohmic currents defined in this patch.
    inline uint countOhmicCurrs() const noexcept
    { return pOhmicCurrsN; }

    /// Return the local index of ohmic current with global index gidx.
    inline uint ohmiccurrG2L(uint gidx) const noexcept
    { return pOhmicCurr_G2L[gidx]; }

    /// Return the global index of ohmic current with local index lidx.
    inline uint ohmiccurrL2G(uint lidx) const noexcept
    { return pOhmicCurr_L2G[lidx]; }

    /// Return a pointer to ohmic current definition object (type OhmicCurrdef)
    /// specified by local index oclidx.
    OhmicCurrdef * ohmiccurrdef(uint oclidx) const;

    /// Return dependency information of ohmic current with local index oclidx
    /// on species with local index splidx.
    int ohmiccurr_dep_S(uint oclidx, uint splidx) const;

    /// Return the local index of the channel state associated with the
    /// ohmic current with local index oclidx
    uint ohmiccurr_chanstate(uint oclidx) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: GHK CURRENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return the number of GHK currents defined in this patch.
    inline uint countGHKcurrs() const noexcept
    { return pGHKcurrsN; }

    /// Return the local index of GHK current with global index gidx.
    inline uint ghkcurrG2L(uint gidx) const noexcept
    { return pGHKcurr_G2L[gidx]; }

    /// Return the global index of GHK current with local index lidx.
    inline uint GHKcurrL2G(uint lidx) const noexcept
    { return pGHKcurr_L2G[lidx]; }

    /// Return a pointer to GHK current definition object (type GHKcurrdef)
    /// specified by local index lidx.
    GHKcurrdef * ghkcurrdef(uint ghklidx) const;

    /// Return dependency information of ghk current with local index ghklidx
    /// on species with local index splidx.
    int ghkcurr_dep_S(uint ghklidx, uint splidx) const;

    /// Return the local index of the channel state associated with the
    /// ghk current with local index ghklidx
    uint ghkcurr_chanstate(uint ghklidx) const;

    /// Return the local index of the ION associated with the ghk current
    // with local index ghklidx
    uint ghkcurr_ion(uint ghklidx) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: V-DEPENDENT TRANSITIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the number of voltage-dependent transitions defined in this patch
    inline uint countVDepTrans() const noexcept
    { return pVDepTransN; }

    /// Return the local index of v-dep transition with global index gidx.
    inline uint vdeptransG2L(uint gidx) const noexcept
    { return pVDepTrans_G2L[gidx]; }

    /// Return the global index of v-dep transition with local index lidx
    inline uint vdeptransL2G(uint lidx) const noexcept
    { return pVDepTrans_L2G[lidx]; }

    /// Return a pointer to vdeptrans definition object (type VDepTransdef)
    /// specified by local index lidx.
    VDepTransdef * vdeptransdef(uint vdtlidx) const;

    /// Return dependency information of v-dep transition with local index
    /// vdtlidx on species with local index splidx
    int vdeptrans_dep_S(uint vdtlidx, uint splidx) const;

    /// Return the local index of the 'source' channel state of v-dep trans
    /// with local index vdtlidx
    uint vdeptrans_srcchanstate(uint vdtlidx) const;
    uint vdeptrans_dstchanstate(uint vdtlidx) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VOLTAGE-DEPENDENT REACTION RULES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of voltage-dependent reactions that can occur on this
    /// patch.
    inline uint countVDepSReacs() const noexcept
    { return pVDepSReacsN; }

    /// Return the local voltage-dependent reaction index for global index argument.
    ///
    /// \param gidx Global index of the voltage-dependent reaction.
    inline uint vdepsreacG2L(uint gidx) const noexcept
    { return pVDepSReac_G2L[gidx]; }
    /*
    inline gidxT sreacL2G(lidxT idx) const
    { return pSReac_L2G[idx]; }
    */

    /// Return a pointer to reaction definition object (type SReacdef)
    /// specified by local index.
    VDepSReacdef * vdepsreacdef(uint lidx) const;

    /// Warning: these methods perform no error checking!
    /// \todo imcompleted.
    int vdepsreac_dep_I(uint vdsrlidx, uint splidx) const;
    int vdepsreac_dep_S(uint vdsrlidx, uint splidx) const;
    int vdepsreac_dep_O(uint vdsrlidx, uint splidx) const;

    /// Warning: these methods perform no error checking!
    ///
    // Return the beginning and end of the lhs arrays of surface reaction
    // specified by local index argument.
    uint * vdepsreac_lhs_I_bgn(uint lidx) const;
    uint * vdepsreac_lhs_I_end(uint lidx) const;
    uint * vdepsreac_lhs_S_bgn(uint lidx) const;
    uint * vdepsreac_lhs_S_end(uint lidx) const;
    uint * vdepsreac_lhs_O_bgn(uint lidx) const;
    uint * vdepsreac_lhs_O_end(uint lidx) const;

    /// Warning: these methods perform no error checking!
    ///
    // Return the beginning and end of the update arrays of surface reaction
    // specified by local index argument.
    int * vdepsreac_upd_I_bgn(uint lidx) const;
    int * vdepsreac_upd_I_end(uint lidx) const;
    int * vdepsreac_upd_S_bgn(uint lidx) const;
    int * vdepsreac_upd_S_end(uint lidx) const;
    int * vdepsreac_upd_O_bgn(uint lidx) const;
    int * vdepsreac_upd_O_end(uint lidx) const;


    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Set the kcst for surface reaction specified by local index
    ///
    /// \param srlidx Local index of the surface reaction.
    /// \param kcst Rate constant of the surface reaction.
    void setKcst(uint srlidx, double kcst);

    /// Activate or inactivate a surface reaction specified by local index.
    ///
    /// \param srlidx Local index of the surface reaction.
    /// \param active Flag to activate / inactivate the surface reaction.
    void setActive(uint srlidx, bool active);

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////
    // DATA: GENERAL
    ////////////////////////////////////////////////////////////////////////

    /// A pointer to the state definition.
    Statedef                          * pStatedef;

    // A pointer to the steps::wm::Patch object this object defines.
    // steps::wm::Patch                  * pPatch;

    // The string identifier of the patch
    std::string                         pName;

    // The area of the patch
    double                              pArea;

    /// The index of the patch.
    uint                                pIdx;

    // The enclosed surface systems, stored as strings
    std::set<std::string>                 pPssys;

    // Pointers to geom level comps, to be used during setup and NOT later.
    steps::wm::Comp                   * pIcomp;
    steps::wm::Comp                   * pOcomp;

    /// Pointer to inner compartment CompDef.
    Compdef                           * pInner;

    /// Pointer to outer compartment CompDef.
    Compdef                           * pOuter;

    // Keep track of whether setup methods have been called.
    bool                                pSetupRefsdone;
    bool                                pSetupIndsdone;

    ////////////////////////////////////////////////////////////////////////
    // DATA: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Number of species embedded in inner volume (_I), patch (_S)
    /// and outer volume (_O).
    uint                                pSpecsN_I;
    uint                                pSpecsN_S;
    uint                                pSpecsN_O;

    /// Table to resolve species index (global -> local).
    uint                              * pSpec_G2L;

    /// Table to resolve species index (local -> global).
    uint                              * pSpec_L2G;

    // Table of the populations of the species on this patch.
    double                            * pPoolCount;

    // Table of 'clamped' flags on the species
    uint                              * pPoolFlags;

    ////////////////////////////////////////////////////////////////////////
    // DATA: SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Number of surface reactions occurring in patch.
    uint                                pSReacsN;

    /// Table to resolve reaction rule indices (global -> local).
    uint                              * pSReac_G2L;

    /// Table to resolve reaction rule indices (local -> global).
    uint                              * pSReac_L2G;

    // Table of the K-constants of the surface reac rules in this patch
    double                            * pSReacKcst;

    // Table of 'active' flags on the surface reaction rules.
    uint                              * pSReacFlags;

    inline uint _IDX_SReac_I_Spec(uint srlidx)
    { return countSpecs_I() * srlidx; }
    inline uint _IDX_SReac_I_Spec(uint srlidx, uint splidx)
    { return (countSpecs_I() * srlidx) + splidx; }
    inline uint _IDX_SReac_S_Spec(uint srlidx)
    { return countSpecs() * srlidx; }
    inline uint _IDX_SReac_S_Spec(uint srlidx, uint splidx)
    { return (countSpecs() * srlidx) + splidx; }
    inline uint _IDX_SReac_O_Spec(uint srlidx)
    { return countSpecs_O() * srlidx; }
    inline uint _IDX_SReac_O_Spec(uint srlidx, uint splidx)
    { return (countSpecs_O() * srlidx) + splidx; }

    int                               * pSReac_DEP_I_Spec;
    int                               * pSReac_DEP_S_Spec;
    int                               * pSReac_DEP_O_Spec;
    uint                              * pSReac_LHS_I_Spec;
    uint                              * pSReac_LHS_S_Spec;
    uint                              * pSReac_LHS_O_Spec;
    int                               * pSReac_UPD_I_Spec;
    int                               * pSReac_UPD_S_Spec;
    int                               * pSReac_UPD_O_Spec;

    ////////////////////////////////////////////////////////////////////////
    // DATA: SURFACE DIFFUSION RULES
    ////////////////////////////////////////////////////////////////////////

    // Number of surface diffusion rules occurring in patch.
    uint                                pSurfDiffsN;

    // Table to resolve diffusion rule indices (global -> local).
    uint                              * pSurfDiff_G2L;

    // Table to resolve diffusion rule indices (local -> global).
    uint                              * pSurfDiff_L2G;

    // Table of the D-constants of the diffusion rules in this compartment
    double                            * pSurfDiffDcst;

    inline uint _IDX_SurfDiff_Spec(uint sdiff, uint spec) const
    { return (pSpecsN_S * sdiff) + spec; }
    uint                              * pSurfDiff_DEP_Spec;
    uint                              * pSurfDiff_LIG;

    ////////////////////////////////////////////////////////////////////////
    // DATA: VOLTAGE-DEPENDENT REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Number of voltage-dependent reactions occurring in patch.
    uint                                pVDepSReacsN;

    /// Table to resolve reaction rule indices (global -> local).
    uint                              * pVDepSReac_G2L;

    /// Table to resolve reaction rule indices (local -> global).
    uint                              * pVDepSReac_L2G;


    inline uint _IDX_VDepSReac_I_Spec(uint vdsrlidx)
    { return countSpecs_I() * vdsrlidx; }
    inline uint _IDX_VDepSReac_I_Spec(uint vdsrlidx, uint splidx)
    { return (countSpecs_I() * vdsrlidx) + splidx; }
    inline uint _IDX_VDepSReac_S_Spec(uint vdsrlidx)
    { return countSpecs() * vdsrlidx; }
    inline uint _IDX_VDepSReac_S_Spec(uint vdsrlidx, uint splidx)
    { return (countSpecs() * vdsrlidx) + splidx; }
    inline uint _IDX_VDepSReac_O_Spec(uint vdsrlidx)
    { return countSpecs_O() * vdsrlidx; }
    inline uint _IDX_VDepSReac_O_Spec(uint vdsrlidx, uint splidx)
    { return (countSpecs_O() * vdsrlidx) + splidx; }

    int                               * pVDepSReac_DEP_I_Spec;
    int                               * pVDepSReac_DEP_S_Spec;
    int                               * pVDepSReac_DEP_O_Spec;
    uint                              * pVDepSReac_LHS_I_Spec;
    uint                              * pVDepSReac_LHS_S_Spec;
    uint                              * pVDepSReac_LHS_O_Spec;
    int                               * pVDepSReac_UPD_I_Spec;
    int                               * pVDepSReac_UPD_S_Spec;
    int                               * pVDepSReac_UPD_O_Spec;

    ////////////////////////////////////////////////////////////////////////
    // DATA: OHMIC CURRENTS
    ////////////////////////////////////////////////////////////////////////

    uint                                 pOhmicCurrsN;
    uint                               * pOhmicCurr_G2L;
    uint                               * pOhmicCurr_L2G;

    inline uint _IDX_OhmicCurr_Spec(uint ohmicclidx)
    { return countOhmicCurrs() * ohmicclidx; }
    inline uint _IDX_OhmicCurr_Spec(uint ohmicclidx, uint speclidx)
    { return (countOhmicCurrs() * ohmicclidx) + speclidx; }

    int                               * pOhmicCurr_DEP_Spec;
    uint                               * pOhmicCurr_CHANSTATE;

    ////////////////////////////////////////////////////////////////////////
    // DATA: GHK CURRENTS
    ////////////////////////////////////////////////////////////////////////

    uint                                 pGHKcurrsN;
    uint                               * pGHKcurr_G2L;
    uint                               * pGHKcurr_L2G;

    inline uint _IDX_GHKcurr_Spec(uint ghkclidx)
    { return countGHKcurrs() * ghkclidx; }
    inline uint _IDX_GHKcurr_Spec(uint ghkclidx, uint speclidx)
    { return (countGHKcurrs() * ghkclidx) + speclidx; }

    int                               * pGHKcurr_DEP_Spec;
    uint                               * pGHKcurr_CHANSTATE;
    uint                               * pGHKcurr_ION;

    ////////////////////////////////////////////////////////////////////////
    // DATA: VOLTAGE-DEPENDENT TRANSITIONS
    ////////////////////////////////////////////////////////////////////////

    uint                                  pVDepTransN;
    uint                               * pVDepTrans_G2L;
    uint                               * pVDepTrans_L2G;

    inline uint _IDX_VDepTrans_Spec(uint vdtlidx)
    { return countVDepTrans() * vdtlidx; }
    inline uint _IDX_VDepTrans_Spec(uint vdtlidx, uint speclidx)
    { return (countVDepTrans() * vdtlidx) + speclidx; }

    int                                 * pVDepTrans_DEP_Spec;
    uint                              * pVDepTrans_SRCCHANSTATE;
    uint                              * pVDepTrans_DSTCHANSTATE;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_SOLVER_PATCHDEF_HPP

// END
