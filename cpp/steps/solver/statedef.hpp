////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
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
//
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SOLVER_STATEDEF_HPP
#define STEPS_SOLVER_STATEDEF_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/geom/geom.hpp>
#include <steps/model/model.hpp>
#include <steps/rng/rng.hpp>
#include <steps/solver/api.hpp>
#include <steps/geom/patch.hpp>

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

    Statedef(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r);
    ~Statedef(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPARTMENTS
    ////////////////////////////////////////////////////////////////////////

    // Return pointer to Compdef object specified by global index argument.
    Compdef * compdef(uint gidx) const;

    // Return the total number of compartments in the simulation state.
    inline uint countComps(void) const
    { return pCompdefs.size(); }

    // Return the global index of compartment identified by string argument,
    // or object argument.
    // Throw exception if geometry does not contain comp with this identifier.
    uint getCompIdx(std::string const & c) const;
    uint getCompIdx(steps::wm::Comp * comp) const;


    std::vector<Compdef*>::const_iterator bgnComp(void) const
    { return pCompdefs.begin(); }

    std::vector<Compdef*>::const_iterator endComp(void) const
    { return pCompdefs.end(); }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: PATCHES
    ////////////////////////////////////////////////////////////////////////

    // Return pointer to Patchdef object specified by global index argument.
    Patchdef * patchdef(uint gidx) const;

    // Return the total number of patches in the simulation state.
    uint countPatches(void) const
    { return pPatchdefs.size(); }

    // Return the global index of patch identified by string argument,
    // or object argument.
    // Throw exception if geometry does not contain patch with this identifier.
    uint getPatchIdx(std::string const & p) const;
    uint getPatchIdx(steps::wm::Patch * patch) const;

    std::vector<Patchdef*>::const_iterator bgnPatch(void) const
    { return pPatchdefs.begin(); }

    std::vector<Patchdef*>::const_iterator endPatch(void) const
    { return pPatchdefs.end(); }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    // Return the total number of species in the simulation state.
    uint countSpecs(void) const
    { return pSpecdefs.size(); }

    // Return pointer to Specdef object specified by global index argument.
    Specdef * specdef(uint gidx) const;

    // Return the global index of species identified by string argument,
    // or object argument.
    // Throw exception if model does not contain species with this identifier.
    uint getSpecIdx(std::string const & s) const;
    uint getSpecIdx(steps::model::Spec * spec) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: REACTIONS
    ////////////////////////////////////////////////////////////////////////

    // Return the total number of reactions in the simulation state.
    uint countReacs(void) const
    { return pReacdefs.size(); }

    // Return pointer to Reacdef object specified by global index argument.
    Reacdef * reacdef(uint gidx) const;

    // Return the global index of reac identified by string argument,
    // or object argument.
    // Throw exception if model does not contain reac with this identifier.
    uint getReacIdx(std::string const & r) const;
    uint getReacIdx(steps::model::Reac * reac) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    // Return the total number of surface reactions in the simulation state.
    uint countSReacs(void) const
    { return pSReacdefs.size(); }

    // Return pointer to SReacdef object specified by global index argument.
    SReacdef * sreacdef(uint gidx) const;

    // Return the global index of surface reaction identified by string
    // argument, or object argument.
    // Throw exception if model does not contain sreac with this identifier.
    uint getSReacIdx(std::string const & sr) const;
    uint getSReacIdx(steps::model::SReac * sreac) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    // Return the total number of diffusion rules in the simulation state.
    uint countDiffs(void) const
    { return pDiffdefs.size(); }

    // Return pointer to Diffdef object specified by global index argument.
    Diffdef * diffdef(uint gidx) const;

    // Return the global index of diffusion identified by string argument,
    // or object argument.
    // Throw exception if model does not contain diff with this identifier.
    uint getDiffIdx(std::string const & d) const;
    uint getDiffIdx(steps::model::Diff * diff) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STATE
    ////////////////////////////////////////////////////////////////////////

    // Return the current simulation time.
    inline double time(void) const
    { return pTime; }

    inline steps::model::Model * model(void) const
    { return pModel; }

    inline steps::rng::RNG * rng(void) const
    { return pRNG; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: STATE
    ////////////////////////////////////////////////////////////////////////

    // Set the current simulation time.
    void setTime(double t);

    // Increase the current simulation time
    void incTime(double dt);

    // Reset the simulation time to 0s.
    inline void resetTime(void)
    { pTime = 0.0; }

    void incNSteps(uint i = 1);

    inline void resetNSteps(void)
    { pNSteps = 0; }

    inline uint nsteps(void) const
    { return pNSteps; }

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
