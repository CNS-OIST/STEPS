////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2010ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */

#ifndef STEPS_SOLVER_STATEDEF_HPP
#define STEPS_SOLVER_STATEDEF_HPP 1


// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include "../common.h"
#include "../geom/geom.hpp"
#include "../model/model.hpp"
#include "../rng/rng.hpp"
#include "api.hpp"
#include "../geom/patch.hpp"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(solver)

// Forwards declarations
class Compdef;
class Patchdef;
class Specdef;
class Reacdef;
class SReacdef;
class Diffdef;

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
    Statedef(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r);

    /// Destructor
    ~Statedef(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPARTMENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return pointer to Compdef object specified by global index argument.
    ///
    /// \param gidx Global index of the Compdef object.
    Compdef * compdef(uint gidx) const;

    /// Return the total number of compartments in the simulation state.
    inline uint countComps(void) const
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
    std::vector<Compdef*>::const_iterator bgnComp(void) const
    { return pCompdefs.begin(); }

    /// Return the end iterator of the Compdefs objects.
    std::vector<Compdef*>::const_iterator endComp(void) const
    { return pCompdefs.end(); }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: PATCHES
    ////////////////////////////////////////////////////////////////////////

    /// Return pointer to Patchdef object specified by global index argument.
    ///
    /// \param gidx Global index of the patch.
    Patchdef * patchdef(uint gidx) const;

    /// Return the total number of patches in the simulation state.
    uint countPatches(void) const
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
    std::vector<Patchdef*>::const_iterator bgnPatch(void) const
    { return pPatchdefs.begin(); }

    /// Return the end iterator of the Patchdefs objects.
    std::vector<Patchdef*>::const_iterator endPatch(void) const
    { return pPatchdefs.end(); }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of species in the simulation state.
    uint countSpecs(void) const
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

    /// Return the global index of spec identified by object argument.
    ///
    /// \param spec Pointer to the species object.
    /// \exception Throw exception if model does not contain spec with this identifier.
    uint getSpecIdx(steps::model::Spec * spec) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of reactions in the simulation state.
    uint countReacs(void) const
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
    uint countSReacs(void) const
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
    uint countDiffs(void) const
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
    /// \exception Throw exception if model does not contain diff with this identifier.
    uint getDiffIdx(steps::model::Diff * diff) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STATE
    ////////////////////////////////////////////////////////////////////////

    /// Return the current simulation time.
    inline double time(void) const
    { return pTime; }

    /// Return the model object.
    inline steps::model::Model * model(void) const
    { return pModel; }

    /// Return the random number generator object.
    inline steps::rng::RNG * rng(void) const
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
    inline void resetTime(void)
    { pTime = 0.0; }

    /// Increase the time step.
    ///
    /// \param Time step to be increased.
    void incNSteps(uint i = 1);

    /// Reset the time step to 0.
    inline void resetNSteps(void)
    { pNSteps = 0; }

    /// Return current simulation time step.
    inline uint nsteps(void) const
    { return pNSteps; }
    
    inline void setNSteps(uint nsteps)
    { pNSteps = nsteps; }

    ////////////////////////////////////////////////////////////////////////

private:

	steps::model::Model               * pModel;
	steps::wm::Geom                   * pGeom;
	steps::rng::RNG                   * pRNG;

	double                              pTime;

	uint                                pNSteps;

	std::vector<Specdef *>              pSpecdefs;
	std::vector<Compdef *>              pCompdefs;
	std::vector<Patchdef *>             pPatchdefs;
	std::vector<Reacdef *>              pReacdefs;
	std::vector<SReacdef *>             pSReacdefs;
	std::vector<Diffdef *>              pDiffdefs;

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(solver)
END_NAMESPACE(steps)

#endif
// STEPS_SOLVER_STATEDEF_HPP

// END
