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

// STEPS headers.
#include "../common.h"
#include "wmrk4.hpp"
#include "../math/constants.hpp"
#include "../error.hpp"
#include "../solver/statedef.hpp"
#include "../solver/compdef.hpp"
#include "../solver/patchdef.hpp"
#include "../solver/reacdef.hpp"
#include "../solver/sreacdef.hpp"
#include "../solver/types.hpp"

NAMESPACE_ALIAS(steps::wmrk4, swmrk4);
NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

swmrk4::Wmrk4::Wmrk4(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r)
: API(m, g, r)
, pReacMtx(0)
, pUpdMtx(0)
, pSpecs_tot(0)
, pReacs_tot(0)
, pCcst()
, pVals()
, pSFlags()
, pRFlags()
, pNewVals()
, pDyDx()
, pDyDxlhs(0)
, pDT(0.0)
, yt()
, dyt()
, dym()
{
	assert (statedef() != 0);
	assert (model() != 0);
	assert (geom() != 0);
	assert (rng() != 0);

	uint nspecstot = 0;
	uint nreacstot = 0;
	for(uint i=0; i< statedef()->countComps(); ++i)
	{
		nspecstot += statedef()->compdef(i)->countSpecs();
		nreacstot += statedef()->compdef(i)->countReacs();
	}
	for(uint i=0; i< statedef()->countPatches(); ++i)
	{
		nspecstot += statedef()->patchdef(i)->countSpecs();
		nreacstot += statedef()->patchdef(i)->countSReacs();
	}

	pReacMtx = new uint * [nreacstot];
	for (uint i=0; i< nreacstot; ++i)
	{
		pReacMtx[i]= new uint[nspecstot];
		std::fill_n(pReacMtx[i], nspecstot, 0);
	}

	pUpdMtx = new int * [nreacstot];
	for (uint i=0; i< nreacstot; ++i)
	{
		pUpdMtx[i] = new int[nspecstot];
		std::fill_n(pUpdMtx[i], nspecstot, 0);
	}

	pDyDxlhs = new double * [nspecstot];
	for (uint i=0; i< nspecstot; ++i)
	{
		pDyDxlhs[i] = new double[nreacstot];
		std::fill_n(pDyDxlhs[i], nreacstot, 0.0);
	}

	_setup();
	_refill();
	_refillCcst();
}

///////////////////////////////////////////////////////////////////////////////

swmrk4::Wmrk4::~Wmrk4(void)
{
    for(uint i=0; i< pReacs_tot; ++i) delete[] pReacMtx[i];
	delete[] pReacMtx;
	for(uint i=0; i< pReacs_tot; ++i) delete[] pUpdMtx[i];
	delete[] pUpdMtx;
	for (uint i=0; i< pSpecs_tot; ++i) delete[] pDyDxlhs[i];
	delete[] pDyDxlhs;
}

///////////////////////////////////////////////////////////////////////////////

std::string swmrk4::Wmrk4::getSolverName(void) const
{
	return "wmrk4";
}

///////////////////////////////////////////////////////////////////////////////

std::string swmrk4::Wmrk4::getSolverDesc(void) const
{
    return "Runge-Kutta Method in well-mixed conditions";
}

///////////////////////////////////////////////////////////////////////////////

std::string swmrk4::Wmrk4::getSolverAuthors(void) const
{
    return "Iain Hepburn and Stefan Wils";
}

///////////////////////////////////////////////////////////////////////////////

std::string swmrk4::Wmrk4::getSolverEmail(void) const
{
    return "ihepburn@oist.jp";
}


///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::reset(void)
{
	uint comps = statedef()->countComps();
	for (uint i=0; i < comps; ++i) statedef()->compdef(i)->reset();
	uint patches = statedef()->countPatches();
	for (uint i=0; i < patches; ++i) statedef()->compdef(i)->reset();
	statedef()->resetTime();
	// recompute flags and counts vectors in Wmrk4 object
	_refill();

}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::run(double endtime)
{
	if (endtime < statedef()->time())
	{
		std::ostringstream os;
		os << "Endtime is before current simulation time";
	    throw steps::ArgErr(os.str());
	}
	_rksteps(statedef()->time(), endtime);
	statedef()->setTime(endtime);
}

////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::advance(double adv)
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

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::step(void)
{
	assert(pDT > 0.0);
	_rksteps(statedef()->time(), statedef()->time() + pDT);
	statedef()->setTime(statedef()->time() + pDT);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::setDT(double dt)
{
	if (dt < 0.0)
	{
		std::ostringstream os;
		os << "Time step cannot be negative or zero.";
		throw steps::ArgErr(os.str());
	}
	pDT = dt;
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::getTime(void) const
{
	return statedef()->time();
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getCompVol(uint cidx) const
{
	assert(cidx < statedef()->countComps());
	Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	return comp->vol();
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setCompVol(uint cidx, double vol)
{
	assert(cidx < statedef()->countComps());
	Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	comp->setVol(vol);

	// recompute the scaled reaction constants
	_refillCcst();
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getCompCount(uint cidx, uint sidx) const
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint slidx = comp->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED) return 0.0;
	return comp->pools()[slidx];
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setCompCount(uint cidx, uint sidx, double n)
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	assert (n >= 0.0);
	Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint slidx = comp->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED) return;
	comp->setCount(slidx, n);
	// easier to recompute all counts with _refill method
	_refill();					/// may be a better way of doing this
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getCompAmount(uint cidx, uint sidx) const
{
	// the following method does all the necessary argument checking
	double count = _getCompCount(cidx, sidx);
	return count / steps::math::AVOGADRO;
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setCompAmount(uint cidx, uint sidx, double a)
{
	assert(a > 0.0);
	// convert amount in mols to number of molecules
	double a2 = a * steps::math::AVOGADRO;
	// the following method does all the necessary argument checking
	_setCompCount(cidx, sidx, a2);
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getCompConc(uint cidx, uint sidx) const
{
	// the following method does all the necessary argument checking
	double count = _getCompCount(cidx, sidx);
	Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	double vol = comp->vol();
	return count/ (1.0e3 * vol * steps::math::AVOGADRO);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setCompConc(uint cidx, uint sidx, double c)
{
	assert(c >= 0.0);
	assert (cidx < statedef()->countComps());
	Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	double count = c * (1.0e3 * comp->vol() * steps::math::AVOGADRO);
	// the following method does all the necessary argument checking
	_setCompCount(cidx, sidx, count);
}

///////////////////////////////////////////////////////////////////////////////

bool swmrk4::Wmrk4::_getCompClamped(uint cidx, uint sidx) const
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lsidx = comp->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return false;

    return comp->clamped(lsidx);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setCompClamped(uint cidx, uint sidx, bool b)
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lsidx = comp->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return;

    comp->setClamped(lsidx, b);

    // copy flags to this solver
    _refill();
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getCompReacK(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	return comp->kcst(lridx);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setCompReacK(uint cidx, uint ridx, double kf)
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	assert (kf >= 0.0);
	Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	comp->setKcst(lridx, kf);

	// recompute the reaction constants
	_refillCcst();
}

///////////////////////////////////////////////////////////////////////////////

bool swmrk4::Wmrk4::_getCompReacActive(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return false;

	return (comp->active(lridx));
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setCompReacActive(uint cidx, uint ridx, bool a)
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	comp->setActive(lridx, a);

	// copy flags to this solver
	_refill();
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getPatchArea(uint pidx) const
{
	assert(pidx < statedef()->countPatches());
	Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	return patch->area();
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setPatchArea(uint pidx, double area)
{
	assert(pidx < statedef()->countPatches());
	Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	patch->setArea(area);
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getPatchCount(uint pidx, uint sidx) const
{
	assert(pidx < statedef()->countPatches());
	assert(sidx < statedef()->countSpecs());
	Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint slidx = patch->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED) return 0.0;

	return patch->pools()[slidx];
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setPatchCount(uint pidx, uint sidx, double n)
{
	assert(pidx < statedef()->countPatches());
	assert(sidx< statedef()->countSpecs());
	Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint slidx = patch->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED) return;
	patch->setCount(slidx, n);
	// easier to recompute all counts with _refill method
	_refill();
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getPatchAmount(uint pidx, uint sidx) const
{
	// the following method does all the necessary argument checking
	double count = _getPatchCount(pidx, sidx);
	return (count / steps::math::AVOGADRO);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setPatchAmount(uint pidx, uint sidx, double a)
{
	assert(a > 0.0);
	// convert amount in mols to number of molecules
	double a2 = a * steps::math::AVOGADRO;
	// the following method does all the necessary argument checking
	_setPatchCount(pidx, sidx, a2);
}

///////////////////////////////////////////////////////////////////////////////

bool swmrk4::Wmrk4::_getPatchClamped(uint pidx, uint sidx) const
{
	assert(pidx < statedef()->countPatches());
	assert(sidx < statedef()->countSpecs());
	Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lsidx = patch->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return false;

    return patch->clamped(lsidx);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setPatchClamped(uint pidx, uint sidx, bool buf)
{
	assert(pidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lsidx = patch->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return;

    patch->setClamped(lsidx, buf);

    // copy flags to this solver
    _refill();
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getPatchSReacK(uint pidx, uint ridx) const
{
	assert(pidx < statedef()->countPatches());
	assert(ridx < statedef()->countSReacs());
	Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lridx = patch->sreacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	return patch->kcst(lridx);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setPatchSReacK(uint pidx, uint ridx, double kf)
{
	assert(pidx < statedef()->countPatches());
	assert(ridx < statedef()->countSReacs());
	assert (kf >= 0.0);
	Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lridx = patch->sreacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	patch->setKcst(lridx, kf);

	// recompute the reaction constants
	_refillCcst();
}

///////////////////////////////////////////////////////////////////////////////

bool swmrk4::Wmrk4::_getPatchSReacActive(uint pidx, uint ridx) const
{
	assert(pidx < statedef()->countPatches());
	assert(ridx < statedef()->countSReacs());
	Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lridx = patch->sreacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return false;

	return (patch->active(lridx));
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setPatchSReacActive(uint pidx, uint ridx, bool a)
{
	assert(pidx < statedef()->countPatches());
	assert(ridx < statedef()->countSReacs());
	Patchdef * patch = statedef()->patchdef(pidx);
	assert(patch != 0);
	uint lridx = patch->sreacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	patch->setActive(lridx, a);

	// copy flags to this solver
	_refill();
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_ccst(double kcst, double vol, uint order)
{
    double vscale = 1.0e3 * vol * steps::math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;
    if (o1 < 0) o1 = 0;
    return kcst * pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setup(void)

{
	/// cumulative number of species from each compartment and patch
    pSpecs_tot = 0;
	/// cumulative number of reacs from each compartment and patch
	pReacs_tot = 0;

	uint Comps_N = statedef()->countComps();
	uint Patches_N = statedef()->countPatches();

	for(uint i=0; i< Comps_N; ++i)
	{
		pSpecs_tot += statedef()->compdef(i)->countSpecs();
		pReacs_tot += statedef()->compdef(i)->countReacs();
	}
	for(uint i=0; i< Patches_N; ++i)
	{
		pSpecs_tot += statedef()->patchdef(i)->countSpecs();
		pReacs_tot += statedef()->patchdef(i)->countSReacs();
	}
	assert(pSpecs_tot > 0);
	assert(pReacs_tot > 0);

	/// initialise vectors
	for(uint i=0; i< pReacs_tot; ++i)
	{
		pRFlags.push_back(0);
	}
	for(uint i=0; i< pSpecs_tot; ++i)
	{
		pVals.push_back(0.0);
		pSFlags.push_back(0);
		pNewVals.push_back(0.0);
		pDyDx.push_back(0.0);
		yt.push_back(0.0);
		dyt.push_back(0.0);
		dym.push_back(0.0);
	}

	/// fill the reaction matrix
	/// loop over compartments,
	/// then comp reacs and copy compdef LHS values to correct index

	/// set row marker to beginning of matrix for first compartment
	uint rowp = 0;
	/// set column marker to beginning of matrix for first compartment
	uint colp = 0;

	for (uint i=0; i< Comps_N; ++i)
	{
		uint compReacs_N = statedef()->compdef(i)->countReacs();
		uint compSpecs_N = statedef()->compdef(i)->countSpecs();

        for(uint j=0; j< compReacs_N; ++j)
		{
			for(uint k=0; k< compSpecs_N; ++k)
			{
				uint lhs = statedef()->compdef(i)->reac_lhs_bgn(j)[k];
				int upd = statedef()->compdef(i)->reac_upd_bgn(j)[k];
				pReacMtx[rowp + j][colp + k] = lhs;
				pUpdMtx[rowp + j][colp + k] = upd;
			}
			/// set scaled reaction constant
			double reac_kcst = statedef()->compdef(i)->kcst(j);
			double comp_vol = statedef()->compdef(i)->vol();
			uint reac_order = statedef()->compdef(i)->reacdef(j)->order();
			pCcst.push_back(_ccst(reac_kcst, comp_vol, reac_order));
		}
		/// step up markers for next compartment
		rowp += compReacs_N;
		colp += compSpecs_N;
	}

	/// now loop over patches,
	/// then sreacs, filling for correct inner and outer compartments
	for(uint i=0; i< Patches_N; ++i)
	{
		uint patchReacs_N = statedef()->patchdef(i)->countSReacs();
		uint patchSpecs_N_S = statedef()->patchdef(i)->countSpecs();
	    uint patchSpecs_N_I = statedef()->patchdef(i)->countSpecs_I();
		uint patchSpecs_N_O = statedef()->patchdef(i)->countSpecs_O();

		for (uint j=0; j< patchReacs_N; ++j)
		{
			for(uint k=0; k< patchSpecs_N_S; ++k)
			{   uint slhs = statedef()->patchdef(i)->sreac_lhs_S_bgn(j)[k];
				int supd = statedef()->patchdef(i)->sreac_upd_S_bgn(j)[k];
				pReacMtx[rowp + j][colp + k] = slhs;
				pUpdMtx[rowp + j][colp + k] = supd;
			}

			/// fill for inner and outer compartments involved in sreac j
			/// only perform if inner comp exists, similarly for outer comp
			if (statedef()->patchdef(i)->icompdef() != 0)
			{
				/// fetch global index of inner compartment
				uint icompidx = statedef()->patchdef(i)->icompdef()->gidx();
				// marker for correct position of inner compartment in matrix
				uint mtx_icompidx = 0;
				/// step up marker to correct comp
				for (uint l=0; l< icompidx; ++l)
				{
					mtx_icompidx += statedef()->compdef(l)->countSpecs();
				}
				for(uint k=0; k< patchSpecs_N_I; ++k)
				{
					uint ilhs = statedef()->patchdef(i)->sreac_lhs_I_bgn(j)[k];
					int iupd = statedef()->patchdef(i)->sreac_upd_I_bgn(j)[k];
					pReacMtx[rowp + j][mtx_icompidx + k] = ilhs;
					pUpdMtx[rowp + j][mtx_icompidx + k] = iupd;
				}
			}
			if (statedef()->patchdef(i)->ocompdef() != 0)
			{
				uint ocompidx = statedef()->patchdef(i)->ocompdef()->gidx();
				uint mtx_ocompidx =0;
				for (uint l=0; l< ocompidx; ++l)
				{
					mtx_ocompidx += statedef()->compdef(l)->countSpecs();
				}
				for(uint k=0; k< patchSpecs_N_O; ++k)
				{
					uint olhs = statedef()->patchdef(i)->sreac_lhs_O_bgn(j)[k];
					int oupd = statedef()->patchdef(i)->sreac_upd_O_bgn(j)[k];
					pReacMtx[rowp + j][mtx_ocompidx + k] = olhs;
					pUpdMtx[rowp + j][mtx_ocompidx + k] = oupd;
				}
			}
			/// set scaled reaction constant
			/// depends on volume of lhs reaction compartment
			double vol;
			if (statedef()->patchdef(i)->sreacdef(j)->inside() == true)
			{
				assert(statedef()->patchdef(i)->icompdef() != 0);
				vol = statedef()->patchdef(i)->icompdef()->vol();
			}
			else
			{
				assert(statedef()->patchdef(i)->ocompdef() != 0);
				vol = statedef()->patchdef(i)->ocompdef()->vol();
			}
			double sreac_kcst = statedef()->patchdef(i)->kcst(j);
			uint sreac_order = statedef()->patchdef(i)->sreacdef(j)->order();
			pCcst.push_back(_ccst(sreac_kcst, vol, sreac_order));
		}
		/// move markers to next point in matrix
		rowp += patchReacs_N;
		colp += patchSpecs_N_S;
	}

	assert (pCcst.size() == pReacs_tot);
}

////////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_refill(void)
{
	uint Comps_N = statedef()->countComps();
	uint Patches_N = statedef()->countPatches();
	assert (Comps_N > 0);

	uint c_marker = 0;
	uint r_marker = 0;

    for(uint i=0; i< Comps_N; ++i)
	{
		uint comp_Specs_N = statedef()->compdef(i)->countSpecs();
		uint comp_Reacs_N = statedef()->compdef(i)->countReacs();
		Compdef * comp = statedef()->compdef(i);
		assert(comp != 0);
		for (uint j=0; j < comp_Specs_N; ++j)
		{
			pVals[c_marker + j] = comp->pools()[j];
			pSFlags[c_marker + j] = comp->flags()[j];
		}
		for(uint k=0; k< comp_Reacs_N; ++k)
		{
			pRFlags[r_marker + k] = comp->rflags()[k];
		}
		c_marker += comp_Specs_N;
		r_marker += comp_Reacs_N;
	}

	for(uint i=0; i< Patches_N; ++i)
	{
		uint patch_Specs_N = statedef()->patchdef(i)->countSpecs();
		uint patch_Reacs_N = statedef()->patchdef(i)->countSReacs();
		Patchdef * patch = statedef()->patchdef(i);
		assert(patch != 0);
		for (uint j=0; j< patch_Specs_N; ++j)
		{
			pVals[c_marker +j] = patch->pools()[j];
			pSFlags[c_marker + j] = patch->flags()[j];
		}
		for (uint k=0; k< patch_Reacs_N; ++k)
		{
			pRFlags[r_marker + k] = patch->srflags()[k];
		}
		c_marker += patch_Specs_N;
		r_marker += patch_Reacs_N;
	}

	assert(c_marker == pVals.size());
	assert(pVals.size() == pSFlags.size());
	assert(r_marker == pRFlags.size());
	assert(pReacs_tot == pRFlags.size());
	assert(pSFlags.size() == pSpecs_tot);
}

////////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_refillCcst(void)
{
	uint Comps_N = statedef()->countComps();
	uint Patches_N = statedef()->countPatches();
	assert (Comps_N > 0);

	// simplest to reset Ccst vector
	pCcst = std::vector<double>();

	for (uint i=0; i< Comps_N; ++i)
	{
		uint compReacs_N = statedef()->compdef(i)->countReacs();
		uint compSpecs_N = statedef()->compdef(i)->countSpecs();

        for(uint j=0; j< compReacs_N; ++j)
		{
			/// set scaled reaction constant
        	// DEBUG 8/4/09: reaction constants were found from model level objects
        	// so didn't take into account sim-level changes
			double reac_kcst = statedef()->compdef(i)->kcst(j);
			double comp_vol = statedef()->compdef(i)->vol();
			uint reac_order = statedef()->compdef(i)->reacdef(j)->order();
			pCcst.push_back(_ccst(reac_kcst, comp_vol, reac_order));
		}
	}

	/// now loop over patches,
	/// then sreacs, filling for correct inner and outer compartments
	for(uint i=0; i< Patches_N; ++i)
	{
		uint patchReacs_N = statedef()->patchdef(i)->countSReacs();
		uint patchSpecs_N_S = statedef()->patchdef(i)->countSpecs();
	    uint patchSpecs_N_I = statedef()->patchdef(i)->countSpecs_I();
		uint patchSpecs_N_O = statedef()->patchdef(i)->countSpecs_O();

		for (uint j=0; j< patchReacs_N; ++j)
		{
			/// set scaled reaction constant
			/// depends on volume of lhs reaction compartment
			double vol;
			if (statedef()->patchdef(i)->sreacdef(j)->inside() == true)
			{
				assert(statedef()->patchdef(i)->icompdef() != 0);
				vol = statedef()->patchdef(i)->icompdef()->vol();
			}
			else
			{
				assert(statedef()->patchdef(i)->ocompdef() != 0);
				vol = statedef()->patchdef(i)->ocompdef()->vol();
			}
        	// DEBUG 8/4/09: reaction constants were found from model level objects
        	// so didn't take into account sim-level changes
			double sreac_kcst = statedef()->patchdef(i)->kcst(j);
			uint sreac_order = statedef()->patchdef(i)->sreacdef(j)->order();
			pCcst.push_back(_ccst(sreac_kcst, vol, sreac_order));
		}
	}

	assert (pCcst.size() == pReacs_tot);
}

////////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setderivs(dVec & vals, dVec & dydx)
{
	for (uint n=0; n< pSpecs_tot; ++n)
	{
		for(uint r=0; r< pReacs_tot; ++r)
		{
			/// check reaction flags
			if (pRFlags[r] & Statedef::INACTIVE_REACFLAG)						/////////
			{
				pDyDxlhs[n][r] = 0.0;
				continue;
			}

			if (pUpdMtx[r][n] == 0)
			{
				pDyDxlhs[n][r] = 0.0;
				continue;
			}

			pDyDxlhs[n][r] = pUpdMtx[r][n] * pCcst[r];
			for (uint i=0; i < pSpecs_tot; ++i)
			{
				double dydx_lhs_temp = 1.0;
				double val = vals[i];
				switch (pReacMtx[r][i])
				{
					case 3: dydx_lhs_temp *= val;
					case 2: dydx_lhs_temp *= val;
					case 1: dydx_lhs_temp *= val;
					case 0: break;
					/// allow maximum 3 molecules of one species in reaction
					default: assert(0);
				}
				pDyDxlhs[n][r] *= dydx_lhs_temp;
			}
		}
	}

	for (uint n=0; n< pSpecs_tot; ++n)
	{
		dydx[n] = 0.0;
		for (uint r=0; r< pReacs_tot; ++r)
		{
			dydx[n] += pDyDxlhs[n][r];
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_rk4(double pdt)
{
	double dt_2 = pdt/2.0;
	double dt_6 = pdt/6.0;

	for(uint i=0; i< pSpecs_tot; ++i) yt[i] = pVals[i] + (dt_2 * pDyDx[i]);
	_setderivs(yt, dyt);
	for(uint i =0; i< pSpecs_tot; ++i) yt[i]= pVals[i] + (dt_2 * dyt[i]);
	_setderivs(yt, dym);
	for (uint i=0; i< pSpecs_tot; ++i)
	{
		yt[i] = pVals[i] + (pdt * dym[i]);
		dym[i] += dyt[i];
	}
	_setderivs(yt, dyt);
	for (uint i=0; i< pSpecs_tot; ++i)
	{
		pNewVals[i]=pVals[i]+dt_6*(pDyDx[i]+dyt[i]+(2.0*dym[i]));
	}
}

////////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_rksteps(double t1, double t2)
{
	if (t1 == t2) return;
	assert(t1 < t2);
	double t = t1;
	if (pDT <= 0.0)
	{
		std::ostringstream os;
		os << "dt is zero or negative. Call setDT() method.";
		throw steps::ArgErr(os.str());
	}
	if (pDT >= (t2-t1))
	{
		std::ostringstream os;
		os << "dt is larger than simulation step.";
		throw steps::ArgErr(os.str());
	}

	/// step up until over maximum time
	while(t < t2)
	{
		if ((t+pDT) > t2) break;

		_setderivs(pVals, pDyDx);
		_rk4(pDT);
		_update();
		t += pDT;
	}

	////////////////////////////////////////////////////////////////////////////
	// DEBUG: 25/03/09. This is general fix for this solver, inspired by a
	// problem with evaluating two supposedly equal double values as non-equal.
	// Now any discrepancy between simulation time and the end step time is dealt
	// with by solving to the simulation time by passing in the fraction (see below)
	// Changed _rk4() to take the dt as a double argument
	double tfrac = t2-t;
	assert (tfrac >= 0.0);

	// Lets only concern ourselves with fractions greater than 1 percent
	if (tfrac != 0.0 && tfrac/pDT >= 0.01)
	{
		assert (tfrac < pDT);
		_setderivs(pVals, pDyDx);
		_rk4(tfrac);
		_update();
	}
	////////////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_update(void)
{
	/// update local values vector with computed counts
	for (uint i=0; i< pSpecs_tot; ++i)
	{
		/// check clamped flag and only update if not clamped
		if (pSFlags[i] & Statedef::CLAMPED_POOLFLAG)
		{
			continue;
		}
		else
		{
			pVals[i] = pNewVals[i];
		}
	}

	/// update pools with computed values
	uint Comps_N = statedef()->countComps();
	uint Patches_N = statedef()->countPatches();
	uint c_marker = 0;

	for(uint i=0; i< Comps_N; ++i)
	{
		uint comp_Specs_N = statedef()->compdef(i)->countSpecs();
		for (uint j=0; j< comp_Specs_N; ++j)
		{
			double c = pVals[c_marker + j];
		    statedef()->compdef(i)->setCount(j, c);
		}
		c_marker += comp_Specs_N;
	}

	for(uint i=0; i< Patches_N; ++i)
	{
		uint patch_Specs_N = statedef()->patchdef(i)->countSpecs();
		for (uint j=0; j< patch_Specs_N; ++j)
		{
			double c = pVals[c_marker + j];
			statedef()->patchdef(i)->setCount(j, c);
		}
		c_marker += patch_Specs_N;
	}
}

////////////////////////////////////////////////////////////////////////////////

// END



