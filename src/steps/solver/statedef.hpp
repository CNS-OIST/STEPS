/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

#ifndef STEPS_SOLVER_STATEDEF_HPP
#define STEPS_SOLVER_STATEDEF_HPP 1


// STL headers.
#include <string>
#include <vector>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/geom/geom.hpp"
#include "steps/model/model.hpp"
#include "steps/rng/rng.hpp"
#include "steps/solver/api.hpp"
#include "steps/geom/patch.hpp"
#include "steps/geom/diffboundary.hpp"
#include "steps/geom/sdiffboundary.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace solver{

////////////////////////////////////////////////////////////////////////////////

// Forwards declarations
class Compdef;
class Patchdef;
class Specdef;
class Reacdef;
class SReacdef;
class Diffdef;
class Chandef;
class VDepTransdef;
class VDepSReacdef;
class OhmicCurrdef;
class GHKcurrdef;
class DiffBoundarydef;
class SDiffBoundarydef;

// Auxiliary declarations.

typedef Specdef *                        SpecdefP;
typedef std::vector<SpecdefP>            SpecdefPVec;
typedef SpecdefPVec::iterator            SpecdefPVecI;
typedef SpecdefPVec::const_iterator      SpecdefPVecCI;

typedef Compdef *                        CompdefP;
typedef std::vector<CompdefP>            CompdefPVec;
typedef CompdefPVec::iterator            CompdefPVecI;
typedef CompdefPVec::const_iterator      CompdefPVecCI;

typedef Patchdef *                       PatchdefP;
typedef std::vector<PatchdefP>           PatchdefPVec;
typedef PatchdefPVec::iterator           PatchdefPVecI;
typedef PatchdefPVec::const_iterator     PatchdefPVecCI;

typedef Reacdef *                        ReacdefP;
typedef std::vector<ReacdefP>            ReacdefPVec;
typedef ReacdefPVec::iterator            ReacdefPVecI;
typedef ReacdefPVec::const_iterator      ReacdefPVecCI;

typedef SReacdef *                       SReacdefP;
typedef std::vector<SReacdefP>           SReacdefPVec;
typedef SReacdefPVec::iterator           SReacdefPVecI;
typedef SReacdefPVec::const_iterator     SReacdefPVecCI;

typedef Diffdef *                        DiffdefP;
typedef std::vector<DiffdefP>            DiffdefPVec;
typedef DiffdefPVec::iterator            DiffdefPVecI;
typedef DiffdefPVec::const_iterator      DiffdefPVecCI;

typedef Diffdef *                        SurfDiffdefP;
typedef std::vector<DiffdefP>            SurfDiffdefPVec;
typedef DiffdefPVec::iterator            SurfDiffdefPVecI;
typedef DiffdefPVec::const_iterator      SurfDiffdefPVecCI;

typedef Chandef *                        ChandefP;
typedef std::vector<ChandefP>            ChandefPVec;
typedef ChandefPVec::iterator            ChandefPVecI;
typedef ChandefPVec::const_iterator      ChandefPVecCI;

typedef VDepTransdef *                   VDepTransdefP;
typedef std::vector<VDepTransdefP>       VDepTransdefPVec;
typedef VDepTransdefPVec::iterator       VDepTransdefPVecI;
typedef VDepTransdefPVec::const_iterator VDepTransdefPVecCI;

typedef VDepSReacdef *                   VDepSReacdefP;
typedef std::vector<VDepSReacdefP>       VDepSReacdefPVec;
typedef VDepSReacdefPVec::iterator       VDepSReacdefPVecI;
typedef VDepSReacdefPVec::const_iterator VDepSReacdefPVecCI;

typedef OhmicCurrdef *                   OhmicCurrdefP;
typedef std::vector<OhmicCurrdefP>       OhmicCurrdefPVec;
typedef OhmicCurrdefPVec::iterator       OhmicCurrdefPVecI;
typedef OhmicCurrdefPVec::const_iterator OhmicCurrdefPVecCI;

typedef GHKcurrdef *                     GHKcurrdefP;
typedef std::vector<GHKcurrdefP>         GHKcurrdefPVec;
typedef GHKcurrdefPVec::iterator         GHKcurrdefPVecI;
typedef GHKcurrdefPVec::const_iterator   GHKcurrdefPVecCI;

typedef DiffBoundarydef *                    DiffBoundarydefP;
typedef std::vector<DiffBoundarydefP>        DiffBoundarydefPVec;
typedef DiffBoundarydefPVec::iterator        DiffBoundarydefPVecI;
typedef DiffBoundarydefPVec::const_iterator  DiffBoundarydefPVecCI;

typedef SDiffBoundarydef *                    SDiffBoundarydefP;
typedef std::vector<SDiffBoundarydefP>        SDiffBoundarydefPVec;
typedef SDiffBoundarydefPVec::iterator        SDiffBoundarydefPVecI;
typedef SDiffBoundarydefPVec::const_iterator  SDiffBoundarydefPVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Defined State
class Statedef
{

public:

    enum PoolFlags
    {
        CLAMPED_POOLFLAG       = 1
    };
    static const uint PoolFlagDefault = 0;

    enum ReacFlags
    {
        INACTIVE_REACFLAG         = 1,
        KCONST_REACFLAG         = 2
    };
    static const uint ReacFlagDefault = INACTIVE_REACFLAG | KCONST_REACFLAG;

    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param m Pointer to the model object.
    /// \param g Pointer to the geometry container.
    /// \param r Pointer to the random number generator.
    Statedef(steps::model::Model *m, steps::wm::Geom *g, const rng::RNGptr &r);

    /// Destructor
    ~Statedef();

    uint getMembIdx(std::string const & m) const;

    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPARTMENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return pointer to Compdef object specified by global index argument.
    ///
    /// \param gidx Global index of the Compdef object.
    Compdef * compdef(uint gidx) const;

    /// Return the total number of compartments in the simulation state.
    inline uint countComps() const noexcept
    { return pCompdefs.size(); }

    /// Return the global index of compartment identified by string argument.
    ///
    /// \param c Name of the compartment.
    /// \exception Throw exception if geometry does not contain comp with this identifier.
    uint getCompIdx(std::string const & c) const;

    /// Return the global index of compartment identified by  object argument.
    ///
    /// \param comp Pointer to the Comp object..
    /// \exception Throw exception if geometry does not contain comp with this identifier.
    uint getCompIdx(steps::wm::Comp * comp) const;

    /// Return the beginning iterator of the Compdefs objects.
    inline std::vector<Compdef*>::const_iterator bgnComp() const noexcept
    { return pCompdefs.begin(); }

    /// Return the end iterator of the Compdefs objects.
    inline std::vector<Compdef*>::const_iterator endComp() const noexcept
    { return pCompdefs.end(); }
    inline const std::vector<Compdef*>& comps() const noexcept
    { return pCompdefs; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: PATCHES
    ////////////////////////////////////////////////////////////////////////

    /// Return pointer to Patchdef object specified by global index argument.
    ///
    /// \param gidx Global index of the patch.
    Patchdef * patchdef(uint gidx) const;

    /// Return the total number of patches in the simulation state.
    inline uint countPatches() const noexcept
    { return pPatchdefs.size(); }

    /// Return the global index of patch identified by string argument.
    ///
    /// \param p Name of the patch.
    /// \exception Throw exception if geometry does not contain patch with this identifier.
    uint getPatchIdx(std::string const & p) const;

    /// Return the global index of patch identified by string argument.
    ///
    /// \param patch Pointer to the patch.
    /// \exception Throw exception if geometry does not contain patch with this identifier.
    uint getPatchIdx(steps::wm::Patch * patch) const;

    /// Return the beginning iterator of the Patchdefs objects.
    inline std::vector<Patchdef*>::const_iterator bgnPatch() const noexcept
    { return pPatchdefs.begin(); }

    inline const std::vector<Patchdef*>& patches() const noexcept
    { return pPatchdefs; }

    /// Return the end iterator of the Patchdefs objects.
    inline std::vector<Patchdef*>::const_iterator endPatch() const noexcept
    { return pPatchdefs.end(); }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of species in the simulation state.
    inline uint countSpecs() const noexcept
    { return pSpecdefs.size(); }

    /// Return pointer to Specdef object specified by global index argument.
    ///
    /// \param gidx Global index of the species.
    Specdef * specdef(uint gidx) const;

    /// Return the global index of species identified by string argument.
    ///
    /// \param s Name of the species.
    /// \exception Throw exception if model does not contain species with this identifier.
    uint getSpecIdx(std::string const & s) const;
    uint getSpecIdx(steps::model::Spec * spec) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of reactions in the simulation state.
    inline uint countReacs() const noexcept
    { return pReacdefs.size(); }

    /// Return pointer to Reacdef object specified by global index argument.
    ///
    /// \param gidx Global index of the reaction.
    Reacdef * reacdef(uint gidx) const;

    /// Return the global index of reac identified by string argument.
    ///
    /// \param r Name of the reaction.
    /// \exception Throw exception if model does not contain reac with this identifier.
    uint getReacIdx(std::string const & r) const;

    /// Return the global index of reac identified by object argument.
    ///
    /// \param reac Pointer to the reaction object.
    /// \exception Throw exception if model does not contain reac with this identifier.
    uint getReacIdx(steps::model::Reac * reac) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of surface reactions in the simulation state.
    inline uint countSReacs() const noexcept
    { return pSReacdefs.size(); }

    /// Return pointer to SReacdef object specified by global index argument.
    ///
    /// \param gidx Global index of the surface reaction.
    SReacdef * sreacdef(uint gidx) const;

    /// Return the global index of surface reaction identified by string argument.
    ///
    /// \param sr Name of the surface reaction.
    /// \exception Throw exception if model does not contain sreac with this identifier.
    uint getSReacIdx(std::string const & sr) const;
    /// Return the global index of surface reaction identified by object argument.
    ///
    /// \param sreac Pointer to the surface reaction object..
    /// \exception Throw exception if model does not contain sreac with this identifier.
    uint getSReacIdx(steps::model::SReac * sreac) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of diffusion rules in the simulation state.
    inline uint countDiffs() const noexcept
    { return pDiffdefs.size(); }

    /// Return pointer to Diffdef object specified by global index argument.
    ///
    /// \param gidx Global index of the diffusion.
    Diffdef * diffdef(uint gidx) const;

    /// Return the global index of diffusion identified by string argument.
    ///
    /// \param d Name of the diffusion.
    /// \exception Throw exception if model does not contain diff with this identifier.
    uint getDiffIdx(std::string const & d) const;

    /// Return the global index of diffusion identified by object argument.
    ///
    /// \param diff Pointer to the diffusion object.
    /// \exception Throw exception if model does not contain this diff.
    uint getDiffIdx(steps::model::Diff * diff) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of surface diffusion rules in the simulation state.
    inline uint countSurfDiffs() const noexcept
    { return pSurfDiffdefs.size(); }

    /// Return pointer to SurfDiffdef object specified by global index argument.
    ///
    /// \param gidx Global index of the surface diffusion.
    Diffdef * surfdiffdef(uint gidx) const;

    /// Return the global index of surface diffusion identified by string argument.
    ///
    /// \param d Name of the surface diffusion.
    /// \exception Throw exception if model does not contain diff with this identifier.
    uint getSurfDiffIdx(std::string const & d) const;

    /// Return the global index of surface diffusion identified by object argument.
    ///
    /// \param diff Pointer to the surface diffusion object.
    /// \exception Throw exception if model does not contain this diff.
    uint getSurfDiffIdx(steps::model::Diff * diff) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VOLTAGE-DEPENDENT TRANSITIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of voltage-dependent transitions in the simulation state.
    inline uint countVDepTrans() const noexcept
    { return pVDepTransdefs.size(); }

    /// Return pointer to VDepTransdef object specified by global index argument.
    ///
    /// \param gidx Global index of the voltage-dependent transition.
    VDepTransdef * vdeptransdef(uint gidx) const;

    /// Return the global index of voltage-dependent transition identified by string argument.
    ///
    /// \param vdt Name of the voltage-dependent transition.
    /// \exception Throw exception if model does not contain vdeptrans with this identifier.
    uint getVDepTransIdx(std::string const & vdt) const;

    /// Return the global index of voltage-dependent transition identified by object argument.
    ///
    /// \param sreac Pointer to the voltage-dependent transition object..
    /// \exception Throw exception if model does not contain this vdeptrans.
    uint getVDepTransIdx(steps::model::VDepTrans * vdeptrans) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VOLTAGE-DEPENDENT REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of voltage-dependent transitions in the simulation state.
    inline uint countVDepSReacs() const noexcept
    { return pVDepSReacdefs.size(); }

    /// Return pointer to VDepSReacdef object specified by global index argument.
    ///
    /// \param gidx Global index of the voltage-dependent reaction.
    VDepSReacdef * vdepsreacdef(uint gidx) const;

    /// Return the global index of voltage-dependent reaction identified by string argument.
    ///
    /// \param vdt Name of the voltage-dependent reaction.
    /// \exception Throw exception if model does not contain vdepsreac with this identifier.
    uint getVDepSReacIdx(std::string const & vdt) const;

    /// Return the global index of voltage-dependent reaction identified by object argument.
    ///
    /// \param sreac Pointer to the voltage-dependent reaction object..
    /// \exception Throw exception if model does not contain this vdepsreac.
    uint getVDepSReacIdx(steps::model::VDepSReac * vdepsreac) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: OHMIC CURRENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of ohmic currents in the simulation state.
    inline uint countOhmicCurrs() const noexcept
    { return pOhmicCurrdefs.size(); }

    /// Return pointer to OhmicCurr object specified by global index argument.
    ///
    /// \param gidx Global index of the ohmic current.
    OhmicCurrdef * ohmiccurrdef(uint gidx) const;

    /// Return the global index of ohmic current identified by string argument.
    ///
    /// \param oc Name of the ohmic current.
    /// \exception Throw exception if model does not contain ohmic current with this identifier.
    uint getOhmicCurrIdx(std::string const & oc) const;

    /// Return the global index of ohmic current identified by object argument.
    ///
    /// \param ocurr Pointer to the ohmic current object..
    /// \exception Throw exception if model does not contain this ohmiccurr.
    uint getOhmicCurrIdx(steps::model::OhmicCurr * ohmiccurr) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: GHK CURRENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of ghk currents in the simulation state.
    inline uint countGHKcurrs() const noexcept
    { return pGHKcurrdefs.size(); }

    /// Return pointer to GHKcurr object specified by global index argument.
    ///
    /// \param gidx Global index of the ghk current.
    GHKcurrdef * ghkcurrdef(uint gidx) const;

    /// Return the global index of ghk current identified by string argument.
    ///
    /// \param ghk Name of the ghk current.
    /// \exception Throw exception if model does not contain ghk current with this identifier.
    uint getGHKcurrIdx(std::string const & ghk) const;

    /// Return the global index of ghk current identified by object argument.
    ///
    /// \param ocurr Pointer to the ghk current object..
    /// \exception Throw exception if model does not contain this ghkcurr.
    uint getGHKcurrIdx(steps::model::GHKcurr * ghkcurr) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION BOUNDARY
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of diffusion boundaries in the simulation state.
    inline uint countDiffBoundaries() const noexcept
    { return pDiffBoundarydefs.size(); }

    /// Return pointer to DiffBoundarydef object specified by global index argument.
    ///
    /// \param gidx Global index of the diffusion boundary.
    DiffBoundarydef * diffboundarydef(uint gidx) const;

    /// Return the global index of diffusion boundary identified by string argument.
    ///
    /// \param d Name of the diffusion boundary.
    /// \exception Throw exception if geometry does not contain diff boundary with this identifier.
    uint getDiffBoundaryIdx(std::string const & d) const;

    /// Return the global index of diffusion boundary identified by object argument.
    ///
    /// \param diff Pointer to the diffusion boundary object.
    /// \exception Throw exception if geoemtry does not contain diff boundary with this identifier.
    uint getDiffBoundaryIdx(steps::tetmesh::DiffBoundary * diffb) const;

    /// Return the beginning iterator of the Diffusion Boundary objects.
    inline std::vector<DiffBoundarydef *>::const_iterator bgnDiffBoundary() const noexcept
    { return pDiffBoundarydefs.begin(); }

    /// Return the end iterator of the Diffusion Boundary objects.
    inline std::vector<DiffBoundarydef *>::const_iterator endDiffBoundary() const noexcept
    { return pDiffBoundarydefs.end(); }
    inline const std::vector<DiffBoundarydef*>& diffBoundaries() const noexcept
    { return pDiffBoundarydefs; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE DIFFUSION BOUNDARY
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of surface diffusion boundaries in the simulation state.
    inline uint countSDiffBoundaries() const noexcept
    { return pSDiffBoundarydefs.size(); }

    /// Return pointer to SDiffBoundarydef object specified by global index argument.
    ///
    /// \param gidx Global index of the surface diffusion boundary.
    SDiffBoundarydef * sdiffboundarydef(uint gidx) const;

    /// Return the global index of surface diffusion boundary identified by string argument.
    ///
    /// \param d Name of the surface diffusion boundary.
    /// \exception Throw exception if geometry does not contain sdiff boundary with this identifier.
    uint getSDiffBoundaryIdx(std::string const & d) const;

    /// Return the global index of surface diffusion boundary identified by object argument.
    ///
    /// \param sdiff Pointer to the surface diffusion boundary object.
    /// \exception Throw exception if geometry does not contain sdiff boundary with this identifier.
    uint getSDiffBoundaryIdx(steps::tetmesh::SDiffBoundary * sdiffb) const;

    /// Return the beginning iterator of the Surface Diffusion Boundary objects.
    inline std::vector<SDiffBoundarydef *>::const_iterator bgnSDiffBoundary() const noexcept
    { return pSDiffBoundarydefs.begin(); }

    /// Return the end iterator of the Surface Diffusion Boundary objects.
    inline std::vector<SDiffBoundarydef *>::const_iterator endSDiffBoundary() const noexcept
    { return pSDiffBoundarydefs.end(); }

    inline const std::vector<SDiffBoundarydef*>& sdiffBoundaries() const noexcept
    { return pSDiffBoundarydefs; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STATE
    ////////////////////////////////////////////////////////////////////////

    /// Return the current simulation time.
    inline double time() const noexcept
    { return pTime; }

    /// Return the model object.
    inline steps::model::Model * model() const noexcept
    { return pModel; }

    /// Return the random number generator object.
    inline const steps::rng::RNGptr& rng() const noexcept
    { return pRNG; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: STATE
    ////////////////////////////////////////////////////////////////////////

    /// Set the current simulation time.
    ///
    /// \param t Simulation time.
    void setTime(double t);

    /// Increase the current simulation time.
    ///
    /// \param dt Discrete time.
    void incTime(double dt);

    /// Reset the simulation time to 0s.
    inline void resetTime() noexcept
    { pTime = 0.0; }

    /// Increase the time step.
    ///
    /// \param Time step to be increased.
    void incNSteps(uint i = 1);

    /// Reset the time step to 0.
    inline void resetNSteps() noexcept
    { pNSteps = 0; }

    /// Return current simulation time step.
    inline uint nsteps() const noexcept
    { return pNSteps; }

    inline void setNSteps(uint nsteps) noexcept
    { pNSteps = nsteps; }

    ////////////////////////////////////////////////////////////////////////

private:

    steps::model::Model               * pModel;
    steps::wm::Geom                   * pGeom;
    const steps::rng::RNGptr            pRNG;

    double                              pTime;

    uint                                pNSteps;

    std::vector<Specdef *>              pSpecdefs;
    std::vector<Chandef *>              pChandefs;
    std::vector<Compdef *>              pCompdefs;
    std::vector<Patchdef *>             pPatchdefs;
    std::vector<Reacdef *>              pReacdefs;
    std::vector<SReacdef *>             pSReacdefs;

    std::vector<Diffdef *>              pDiffdefs;
    std::vector<Diffdef *>              pSurfDiffdefs;

    std::vector<DiffBoundarydef *>      pDiffBoundarydefs;
    std::vector<SDiffBoundarydef *>      pSDiffBoundarydefs;
    std::vector<VDepTransdef *>         pVDepTransdefs;
    std::vector<VDepSReacdef *>         pVDepSReacdefs;
    std::vector<OhmicCurrdef *>         pOhmicCurrdefs;
    std::vector<GHKcurrdef *>           pGHKcurrdefs;

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_SOLVER_STATEDEF_HPP

// END
