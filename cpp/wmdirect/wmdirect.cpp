////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
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

// Standard library & STL headers.
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <sstream>
#include <algorithm>

// STEPS headers.
#include "../common.h"
#include "wmdirect.hpp"
#include "kproc.hpp"
#include "comp.hpp"
#include "patch.hpp"
#include "reac.hpp"
#include "sreac.hpp"
#include "../math/constants.hpp"
#include "../error.hpp"
#include "../solver/statedef.hpp"
#include "../solver/compdef.hpp"
#include "../solver/patchdef.hpp"
#include "../solver/reacdef.hpp"
#include "../solver/sreacdef.hpp"
#include "../solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

#define SCHEDULEWIDTH 32
#define MAXLEVELS 10

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::wmdirect, swmd);
NAMESPACE_ALIAS(steps::solver, ssolver);
NAMESPACE_ALIAS(steps::math, smath);

////////////////////////////////////////////////////////////////////////////////

/// Unary function that calls the array delete[] operator on pointers. Easy
/// to use with STL/Boost (see steps::tools::DeletePointer).
///
struct DeleteArray
{
    template <typename Type> void operator() (Type * pointer) const
    {
        delete[] pointer;
    }
};

////////////////////////////////////////////////////////////////////////////////

void swmd::schedIDXSet_To_Vec(swmd::SchedIDXSet const & s, swmd::SchedIDXVec & v)
{
    v.resize(s.size());
    std::copy(s.begin(), s.end(), v.begin());
}

////////////////////////////////////////////////////////////////////////////////

swmd::Wmdirect::Wmdirect(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r)
: API(m, g, r)
, pKProcs()
, pComps()
, pCompMap()
, pPatches()
, pA0(0.0)
, pLevelSizes()
, pLevels()
, pBuilt(false)
, pIndices(0)
, pMaxUpSize(0)
, pRannum(0)
{
	assert (model() != 0);
	assert (geom() != 0);
	assert (rng() != 0);
	ssolver::CompDefPVecCI c_end = statedef()->endComp();
    for (ssolver::CompDefPVecCI c = statedef()->bgnComp(); c != c_end; ++c)
    {
        uint compdef_gidx = (*c)->gidx();
        uint comp_idx = _addComp(*c);
        assert(compdef_gidx == comp_idx);
    }

    // Create the actual patches.
    ssolver::PatchDefPVecCI p_end = statedef()->endPatch();
    for (ssolver::PatchDefPVecCI p = statedef()->bgnPatch(); p != p_end; ++p)
    {
        uint patchdef_gidx = (*p)->gidx();
        uint patch_idx = _addPatch(*p);
        assert(patchdef_gidx == patch_idx);
    }

    _setup();
}

////////////////////////////////////////////////////////////////////////////////

swmd::Wmdirect::~Wmdirect(void)
{

    CompPVecCI comp_e = pComps.end();
    for (CompPVecCI c = pComps.begin(); c != comp_e; ++c) delete *c;
    PatchPVecCI patch_e = pPatches.end();
    for (PatchPVecCI p = pPatches.begin(); p != patch_e; ++p) delete *p;

    std::for_each(pLevels.begin(), pLevels.end(), DeleteArray());
	delete[] pIndices;
    delete[] pRannum;
}

////////////////////////////////////////////////////////////////////////////////

uint swmd::Wmdirect::_addComp(steps::solver::Compdef * cdef)
{
    swmd::Comp * comp = new Comp(cdef);
    assert(comp != 0);
    uint compidx = pComps.size();
    pComps.push_back(comp);
    pCompMap[cdef] = comp;
    return compidx;
}

////////////////////////////////////////////////////////////////////////////////

uint swmd::Wmdirect::_addPatch(steps::solver::Patchdef * pdef)
{
    Comp * icomp = 0;
    Comp * ocomp = 0;
    if (pdef->icompdef()) icomp = pCompMap[pdef->icompdef()];
    if (pdef->ocompdef()) ocomp = pCompMap[pdef->ocompdef()];
    swmd::Patch * patch = new Patch(pdef, icomp, ocomp);
    assert(patch != 0);
    uint patchidx = pPatches.size();
    pPatches.push_back(patch);
    return patchidx;
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_setup(void)
{
	CompPVecCI c_end = pComps.end();
	for (CompPVecCI c = pComps.begin(); c != c_end; ++c)
	{
		(*c)->setupKProcs(this);
	}

	PatchPVecCI p_end = pPatches.end();
	for (PatchPVecCI p = pPatches.begin(); p != p_end; ++p)
	{

		(*p)->setupKProcs(this);
	}

	// Resolve all dependencies
	for (CompPVecCI c = pComps.begin(); c != c_end; ++c)
	{
		KProcPVecCI kprocend = (*c)->kprocEnd();
		for (KProcPVecCI k = (*c)->kprocBegin(); k != kprocend; ++k)
		{
		    (*k)->setupDeps();
		}
	}
	for (PatchPVecCI p = pPatches.begin(); p != p_end; ++p)
	{

	    KProcPVecCI kprocend = (*p)->kprocEnd();
	    for (KProcPVecCI k = (*p)->kprocBegin(); k != kprocend; ++k)
	    {
	        (*k)->setupDeps();
	    }
	}

	_build();
}

////////////////////////////////////////////////////////////////////////////////

std::string swmd::Wmdirect::getSolverName(void) const
{
	return "wmdirect";
}

////////////////////////////////////////////////////////////////////////////////

std::string swmd::Wmdirect::getSolverDesc(void) const
{
	return "SSA Direct Method in well-mixed conditions";
}

////////////////////////////////////////////////////////////////////////////////

std::string swmd::Wmdirect::getSolverAuthors(void) const
{
	return "Stefan Wils and Iain Hepburn";
}

////////////////////////////////////////////////////////////////////////////////

std::string swmd::Wmdirect::getSolverEmail(void) const
{
	return "stefan@tnb.ua.ac.be, ihepburn@oist.jp";
}


////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::reset(void)
{
	uint comps = statedef()->countComps();
	for (uint i=0; i < comps; ++i) statedef()->compdef(i)->reset();
	uint patches = statedef()->countPatches();
	for (uint i=0; i < patches; ++i) statedef()->patchdef(i)->reset();

    std::for_each(pComps.begin(), pComps.end(), std::mem_fun(&Comp::reset));
    std::for_each(pPatches.begin(), pPatches.end(), std::mem_fun(&Patch::reset));

	statedef()->resetTime();
	statedef()->resetNSteps();

	_reset();
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::run(double endtime)
{
	if (endtime < statedef()->time())
	{
		std::ostringstream os;
		os << "Endtime is before current simulation time";
	    throw steps::ArgErr(os.str());
	}
	while (statedef()->time() < endtime)
	{
		swmd::KProc * kp = _getNext();
		if (kp == 0) break;
		double a0 = getA0();
		if (a0 == 0.0) break;
		double dt = rng()->getExp(a0);
		if ((statedef()->time() + dt) > endtime) break;
		_executeStep(kp, dt);
	}
	statedef()->setTime(endtime);
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::advance(double adv)
{
	if (adv < 0.0)
	{
		std::ostringstream os;
		os << "Time to advance cannot be negative";
	    throw steps::ArgErr(os.str());
	}

	double endtime = statedef()->time() + adv;
	run(endtime);
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::step(void)
{
	swmd::KProc * kp = _getNext();
	if (kp == 0) return;
	double a0 = getA0();
	if (a0 == 0.0) return;
	double dt = rng()->getExp(a0);
	_executeStep(kp, dt);
}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::getTime(void) const
{
	return statedef()->time();
}

////////////////////////////////////////////////////////////////////////

uint swmd::Wmdirect::getNSteps(void) const
{
    return statedef()->nsteps();
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::setTime(double time)
{
	statedef()->setTime(time);
}


////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::setNSteps(uint nsteps)
{
    statedef()->setNSteps(nsteps);
}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getCompVol(uint cidx) const
{
	assert(cidx < statedef()->countComps());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	return comp->vol();
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_setCompVol(uint cidx, double vol)
{
	assert(cidx < statedef()->countComps());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	comp->setVol(vol);

	// Reset the reaction C constants
	swmd::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);
    std::for_each(lcomp->kprocBegin(), lcomp->kprocEnd(), std::mem_fun(&KProc::resetCcst));
    // Rates have changed
    _reset();
}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getCompCount(uint cidx, uint sidx) const
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint slidx = comp->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED) return 0.0;
	return comp->pools()[slidx];
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_setCompCount(uint cidx, uint sidx, double n)
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint slidx = comp->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED) return;

	assert (n >= 0.0);
	double n_int = std::floor(n);
	double n_frc = n-n_int;
	uint c = static_cast<uint>(n_int);
	if (n_frc > 0.0)
	{
		double rand01 = rng()->getUnfIE();
		if (rand01 < n_frc)c++;
	}
	double n_double = static_cast<double>(c);
	comp->setCount(slidx, n_double);
	// Rates have changed
	_reset();
}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getCompAmount(uint cidx, uint sidx) const
{
	// the following method does all the necessary argument checking
	double count = _getCompCount(cidx, sidx);
	return (count / smath::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_setCompAmount(uint cidx, uint sidx, double a)
{
	// convert amount in mols to number of molecules
	double a2 = a * steps::math::AVOGADRO;
	// the following method does all the necessary argument checking
	_setCompCount(cidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getCompConc(uint cidx, uint sidx) const
{
	// the following method does all the necessary argument checking
	double count = _getCompCount(cidx, sidx);
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	double vol = comp->vol();
	return count/ (1.0e3 * vol * steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_setCompConc(uint cidx, uint sidx, double c)
{
	assert(c >= 0.0);
	assert (cidx < statedef()->countComps());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	double count = c * (1.0e3 * comp->vol() * steps::math::AVOGADRO);
	// the following method does all the necessary argument checking
	_setCompCount(cidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////

bool swmd::Wmdirect::_getCompClamped(uint cidx, uint sidx) const
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lsidx = comp->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return false;

    return comp->clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_setCompClamped(uint cidx, uint sidx, bool b)
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lsidx = comp->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return;

    comp->setClamped(lsidx, b);

}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getCompReacK(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	return comp->kcst(lridx);
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_setCompReacK(uint cidx, uint ridx, double kf)
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	assert (kf >= 0.0);
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	comp->setKcst(lridx, kf);

	// Reset the reaction C constants
	swmd::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);

	steps::wmdirect::KProc * lreac = lcomp->reac(lridx);
	assert (lreac->defr() == comp->reacdef(lridx));
	lreac->resetCcst();

	// Rates have changed
	_reset();
}

////////////////////////////////////////////////////////////////////////

bool swmd::Wmdirect::_getCompReacActive(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return false;

	return (comp->active(lridx));
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_setCompReacActive(uint cidx, uint ridx, bool a)
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	comp->setActive(lridx, a);

    // It's cheaper to just recompute everything.
    _reset();
}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getCompReacC(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// The 'local' Comp object has same index as solver::Compdef object
	swmd::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);
	// Reacdef local indices in Compdef object also have same index
	// as Reacs in Comp object
	swmd::KProc * lreac = lcomp->reac(lridx);
	assert (lreac->defr() == comp->reacdef(lridx));

	return lreac->c();
}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getCompReacH(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// The 'local' Comp object has same index as solver::Compdef object
	swmd::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);
	// Reacdef local indices in Compdef object also have same index
	// as Reacs in Comp object
	swmd::KProc * lreac = lcomp->reac(lridx);
	assert (lreac->defr() == comp->reacdef(lridx));

	return lreac->h();
}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getCompReacA(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// The 'local' Comp object has same index as solver::Compdef object
	swmd::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);
	// Reacdef local indices in Compdef object also have same index
	// as Reacs in Comp object
	swmd::KProc * lreac = lcomp->reac(lridx);
	assert (lreac->defr() == comp->reacdef(lridx));

	return lreac->rate();
}

////////////////////////////////////////////////////////////////////////

uint swmd::Wmdirect::_getCompReacExtent(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0;

	// The 'local' Comp object has same index as solver::Compdef object
	swmd::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);
	// Reacdef local indices in Compdef object also have same index
	// as Reacs in Comp object
	swmd::KProc * lreac = lcomp->reac(lridx);
	assert (lreac->defr() == comp->reacdef(lridx));

	return lreac->getExtent();
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_resetCompReacExtent(uint cidx, uint ridx)
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	// The 'local' Comp object has same index as solver::Compdef object
	swmd::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);
	// Reacdef local indices in Compdef object also have same index
	// as Reacs in Comp object
	swmd::KProc * lreac = lcomp->reac(lridx);
	assert (lreac->defr() == comp->reacdef(lridx));

	lreac->resetExtent();
}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getPatchArea(uint pidx) const
{
	assert(pidx < statedef()->countPatches());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	return patch->area();
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_setPatchArea(uint pidx, double area)
{
	assert(pidx < statedef()->countPatches());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	patch->setArea(area);
}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getPatchCount(uint pidx, uint sidx) const
{
	assert(pidx < statedef()->countPatches());
	assert(sidx < statedef()->countSpecs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint slidx = patch->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED) return 0.0;

	return patch->pools()[slidx];
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_setPatchCount(uint pidx, uint sidx, double n)
{
	assert(pidx < statedef()->countPatches());
	assert(sidx< statedef()->countSpecs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint slidx = patch->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED) return;

	double n_int = std::floor(n);
	double n_frc = n-n_int;
	uint c = static_cast<uint>(n_int);
	if (n_frc > 0.0)
	{
		double rand01 = statedef()->rng()->getUnfIE();
		if (rand01 < n_frc)c++;
	}
	n_int = static_cast<double>(c);
	patch->setCount(slidx, n_int);
	_reset();
}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getPatchAmount(uint pidx, uint sidx) const
{
	// the following method does all the necessary argument checking
	double count = _getPatchCount(pidx, sidx);
	return (count / steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_setPatchAmount(uint pidx, uint sidx, double a)
{
	assert(a > 0.0);
	// convert amount in mols to number of molecules
	double a2 = a * steps::math::AVOGADRO;
	// the following method does all the necessary argument checking
	_setPatchCount(pidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////

bool swmd::Wmdirect::_getPatchClamped(uint pidx, uint sidx) const
{
	assert(pidx < statedef()->countPatches());
	assert(sidx < statedef()->countSpecs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lsidx = patch->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return false;

    return patch->clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_setPatchClamped(uint pidx, uint sidx, bool buf)
{
	assert(pidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lsidx = patch->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return;

    patch->setClamped(lsidx, buf);

}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getPatchSReacK(uint pidx, uint ridx) const
{
	assert(pidx < statedef()->countPatches());
	assert(ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lridx = patch->sreacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	return patch->kcst(lridx);
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_setPatchSReacK(uint pidx, uint ridx, double kf)
{
	assert(pidx < statedef()->countPatches());
	assert(ridx < statedef()->countSReacs());
	assert (kf >= 0.0);
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lridx = patch->sreacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	patch->setKcst(lridx, kf);

	// The 'local' Patch object has same index as solver::Patchdef object
	swmd::Patch * lpatch = pPatches[pidx];
	assert(lpatch->def() == patch);
	// SReacdef local indices in Patchdef object also have same index
	//  as SReacs in Patch object
	swmd::KProc * lsreac = lpatch->sreac(lridx);
	assert (lsreac->defsr() == patch->sreacdef(lridx));
	lsreac->resetCcst();

	// Rates have changed
	_reset();
}

////////////////////////////////////////////////////////////////////////

bool swmd::Wmdirect::_getPatchSReacActive(uint pidx, uint ridx) const
{
	assert(pidx < statedef()->countPatches());
	assert(ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lridx = patch->sreacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return false;

	return (patch->active(lridx));
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_setPatchSReacActive(uint pidx, uint ridx, bool a)
{
	assert(pidx < statedef()->countPatches());
	assert(ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lridx = patch->sreacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	patch->setActive(lridx, a);

    // It's cheaper to just recompute everything
	_reset();
}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getPatchSReacC(uint pidx, uint ridx) const
{
	assert(pidx < statedef()->countPatches());
	assert(ridx < statedef()->countSReacs());

	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lridx = patch->sreacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// The 'local' Patch object has same index as solver::Patchdef object
	swmd::Patch * lpatch = pPatches[pidx];
	assert(lpatch->def() == patch);
	// SReacdef local indices in Patchdef object also have same index
	//  as SReacs in Patch object
	swmd::KProc * lsreac = lpatch->sreac(lridx);
	assert (lsreac->defsr() == patch->sreacdef(lridx));

	return lsreac->c();
}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getPatchSReacH(uint pidx, uint ridx) const
{
	assert(pidx < statedef()->countPatches());
	assert(ridx < statedef()->countSReacs());

	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lridx = patch->sreacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// The 'local' Patch object has same index as solver::Patchdef object
	swmd::Patch * lpatch = pPatches[pidx];
	assert(lpatch->def() == patch);
	// SReacdef local indices in Patchdef object also have same index
	//  as SReacs in Patch object
	swmd::KProc * lsreac = lpatch->sreac(lridx);
	assert (lsreac->defsr() == patch->sreacdef(lridx));

	return lsreac->h();
}

////////////////////////////////////////////////////////////////////////

double swmd::Wmdirect::_getPatchSReacA(uint pidx, uint ridx) const
{
	assert(pidx < statedef()->countPatches());
	assert(ridx < statedef()->countSReacs());

	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lridx = patch->sreacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// The 'local' Patch object has same index as solver::Patchdef object
	swmd::Patch * lpatch = pPatches[pidx];
	assert(lpatch->def() == patch);
	// SReacdef local indices in Patchdef object also have same index
	//  as SReacs in Patch object
	swmd::KProc * lsreac = lpatch->sreac(lridx);
	assert (lsreac->defsr() == patch->sreacdef(lridx));

	return lsreac->rate();

}

////////////////////////////////////////////////////////////////////////

uint swmd::Wmdirect::_getPatchSReacExtent(uint pidx, uint ridx) const
{
	assert(pidx < statedef()->countPatches());
	assert(ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lridx = patch->sreacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0;

	// The 'local' Patch object has same index as solver::Patchdef object
	swmd::Patch * lpatch = pPatches[pidx];
	assert(lpatch->def() == patch);
	// SReacdef local indices in Patchdef object also have same index
	//  as SReacs in Patch object
	swmd::KProc * lsreac = lpatch->sreac(lridx);
	assert (lsreac->defsr() == patch->sreacdef(lridx));

	return lsreac->getExtent();

}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_resetPatchSReacExtent(uint pidx, uint ridx)
{
	assert(pidx < statedef()->countPatches());
	assert(ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lridx = patch->sreacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	// The 'local' Patch object has same index as solver::Patchdef object
	swmd::Patch * lpatch = pPatches[pidx];
	assert(lpatch->def() == patch);
	// SReacdef local indices in Patchdef object also have same index
	//  as SReacs in Patch object
	swmd::KProc * lsreac = lpatch->sreac(lridx);
	assert (lsreac->defsr() == patch->sreacdef(lridx));

	lsreac->resetExtent();

}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::addKProc(KProc * kp)
{
	assert (kp != 0);

	SchedIDX nidx = pKProcs.size();
	pKProcs.push_back(kp);
	kp->setSchedIDX(nidx);
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_build(void)
{
	assert (pBuilt == false);

    // Setup level.
    uint clsize = pKProcs.size();
    if (clsize == 0) return;

    // Work up.
    uint clevel = 0;
    do
    {
        // Make sure the new size is a multiple of SCHEDULEWIDTH.
        uint extra = clsize % SCHEDULEWIDTH;
        if (extra != 0) clsize += SCHEDULEWIDTH - extra;

        // Create the level and add it.
        double * level = new double[clsize];
        std::fill_n(level, clsize, 0.0);
        pLevelSizes.push_back(clsize);
        pLevels.push_back(level);

        // Prepare for next level.
        clevel++;
        clsize = clsize / SCHEDULEWIDTH;
    }
    while (clsize > 1);

    // Set top level.
    pA0 = 0.0;

	// Time to create ONE indices table to hold the run's present reaction of
    // choice's update vector. This will be re-used and replace old-version's
    // hard-coded table in _update. Size is the maximum possible, found by looping
    // over all KProcs. This little bit of computational time is well worth all
    // that needless memory allocation
	uint maxupvecsize = 0;
	KProcPVecCI kproc_end = pKProcs.end();
	for (KProcPVecCI kproc = pKProcs.begin(); kproc != kproc_end; ++kproc)
	{
		if ((*kproc)->updVecSize() > maxupvecsize) maxupvecsize = (*kproc)->updVecSize();
	}

	pMaxUpSize = maxupvecsize;
	pIndices = new uint[pMaxUpSize];

    // Also let's create a random number holder-table,
    // size of number of KProcs % SCHEDULEWIDTH or pLevels.size()
    // This will be re-used in _getNext as opposed to hard-coded (again maximum
    // limit).
    uint lsize = pLevels.size();
    pRannum = new double[lsize];


    pBuilt = true;
}

////////////////////////////////////////////////////////////////////////

swmd::KProc * swmd::Wmdirect::_getNext(void) const
{
    assert(pA0 >= 0.0);
    // Quick check to see whether nothing is there.
    if (pA0 == 0.0) return 0;

    // Start at top level.
    uint clevel = pLevels.size();
    // And start at the first node of that level.
    uint cur_node = 0;

    // Prepare random numbers.
    for (uint i = 0; i < clevel; ++i)
    {
        pRannum[i] = rng()->getUnfIE();
    }

    // Run until top level.
    double a0 = pA0;
    double * level = 0;
    while (clevel != 0)
    {
        // Decrease the current level.
        clevel--;
        // and start looking in the right place.
        cur_node *= SCHEDULEWIDTH;
        uint max_node = cur_node + SCHEDULEWIDTH;

        // Fetch the level.
        level = pLevels[clevel];

        // Compute local selector.
        double selector = pRannum[clevel] * a0;

        // Compare.
        double accum = 0.0;
        // 27/10/09 I.H. 'old' removed from for loop because not used.
        // double old = 0.0;
        double curval = 0.0;
        for (uint i = 0; i < SCHEDULEWIDTH; ++i)
        {
            curval = level[cur_node];
            if (selector < curval + accum) break;
            accum += curval;
            // old = accum;
            cur_node++;
        }

        // Checks.
        assert(cur_node < max_node);
        assert(curval > 0.0);
        a0 = curval;
    }

    // Check.
    assert(cur_node < pKProcs.size());
    return pKProcs[cur_node];
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_reset(void)
{
    if (pKProcs.size() == 0) return;

    // Reset the basic level: compute rates.
    double * oldlevel = pLevels[0];
    uint cur_node = 0;
    std::vector<KProc*>::iterator kp_end = pKProcs.end();
    for (std::vector<KProc*>::iterator kp = pKProcs.begin(); kp != kp_end; ++kp)
    {
        oldlevel[cur_node++] = (*kp)->rate();
    }

    // Work up.
    for (uint cur_level = 1; cur_level < pLevels.size(); ++cur_level)
    {
        // Compute the number of nodes to reset on this level.
        uint numnodes = pLevelSizes[cur_level - 1] / SCHEDULEWIDTH;

        // Fetch a pointer to this level.
        double * level = pLevels[cur_level];

        // Recompute them.
        uint child_node = 0;
        for (cur_node = 0; cur_node < numnodes; ++cur_node)
        {
            double val = 0.0;
            for (uint i = 0; i < SCHEDULEWIDTH; ++i)
            {
                val += oldlevel[child_node++];
            }
            level[cur_node] = val;
        }

        // Copy the level.
        oldlevel = level;
    }

    // Compute zero propensity.
    pA0 = 0.0;
    for (uint i = 0; i < SCHEDULEWIDTH; ++i)
    {
        pA0 += oldlevel[i];
    }
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_update(SchedIDXVec const & entries)
{
    if (countKProcs() == 0) return;

    // Prefetch zero level.
    double * level0 = pLevels[0];
    // Number of entries.
    assert(entries.size() <= pMaxUpSize);											/////////

    // Recompute rates.
    SchedIDXVecCI sidx_end = entries.end();
    uint prev_e = 0xFFFFFFFF;
    uint cur_e = 0;
    for (SchedIDXVecCI sidx = entries.begin(); sidx != sidx_end; ++sidx)
    {
        // Fetch index.
        uint idx = *sidx;
        // Recompute rate, get difference, and store.
        double newrate = pKProcs[idx]->rate();
        level0[idx] = newrate;

        // Store and collapse if possible.
        idx /= SCHEDULEWIDTH;
        if (prev_e == 0xFFFFFFFF)
        {
            prev_e = 0;
            pIndices[cur_e++] = idx;
        }
        else if (pIndices[prev_e] != idx)
        {
            prev_e = cur_e;
            pIndices[cur_e++] = idx;
        }
    }
    uint nentries = cur_e;

    // Update upper levels.
    uint nlevels = pLevels.size();
    double * prevlevel = pLevels[0];
    for (uint l = 1; l < nlevels; ++l)
    {
        // Update the first entry.
        cur_e = 0;
        prev_e = 0xFFFFFFFF;

        // Fetch a pointer to the current level.
        double * currlevel = pLevels[l];

        // Recompute the entries.
        for (uint e = 0; e < nentries; ++e)
        {
            // Fetch index.
            uint idx = pIndices[e];

            // Recompute.
            double val = 0.0;
            uint idx2 = idx * SCHEDULEWIDTH;
            for (uint i = 0; i < SCHEDULEWIDTH; ++i)
            {
                val += prevlevel[idx2++];
            }
            currlevel[idx] = val;

            // Store and collapse if possible.
            idx /= SCHEDULEWIDTH;
            if (prev_e == 0xFFFFFFFF)
            {
                prev_e = 0;
                pIndices[cur_e++] = idx;
            }
            else if (pIndices[prev_e] != idx)
            {
                prev_e = cur_e;
                pIndices[cur_e++] = idx;
            }
        }

        // Update the pointer to the previous level.
        prevlevel = currlevel;

        // cur_e now is the new number of entries to handle.
        nentries = cur_e;
    }

    // Update zero propensity.
    double * toplevel = pLevels[pLevels.size() - 1];
    pA0 = 0.0;
    for (uint i = 0; i < SCHEDULEWIDTH; ++i)
    {
        pA0 += toplevel[i];
    }
}

////////////////////////////////////////////////////////////////////////

void swmd::Wmdirect::_executeStep(swmd::KProc * kp, double dt)
{
	SchedIDXVec const & upd = kp->apply();
	_update(upd);
	statedef()->incTime(dt);
	statedef()->incNSteps(1);
}

////////////////////////////////////////////////////////////////////////

// END


