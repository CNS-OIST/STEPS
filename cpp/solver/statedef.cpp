
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

// STL headers.
#include <string>
#include <sstream>
#include <cassert>

// STEPS headers.
#include "../common.h"
#include "../error.hpp"
#include "../geom/geom.hpp"
#include "../model/model.hpp"
#include "../rng/rng.hpp"
#include "statedef.hpp"
#include "specdef.hpp"
#include "compdef.hpp"
#include "patchdef.hpp"
#include "reacdef.hpp"
#include "sreacdef.hpp"
#include "diffdef.hpp"

NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

ssolver::Statedef::Statedef(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r)
: pModel(m)
, pGeom(g)
, pRNG(r)
, pTime(0.0)
, pNSteps(0)
, pSpecdefs()
, pCompdefs()
, pPatchdefs()
, pReacdefs()
, pSReacdefs()
, pDiffdefs()

{
    assert(pModel != 0);
    assert(pGeom != 0);
    assert(pRNG != 0);

    uint nspecs = pModel->_countSpecs();
    assert (nspecs > 0);
    for (uint sidx = 0; sidx < nspecs; ++sidx)
    {
    	ssolver::Specdef * specdef = new Specdef(this, sidx,  pModel->_getSpec(sidx));
    	assert (specdef != 0);
    	pSpecdefs.push_back(specdef);
    }

    uint nreacs = pModel->_countReacs();
    for (uint ridx = 0; ridx < nreacs; ++ridx)
    {
    	ssolver::Reacdef * reacdef = new Reacdef(this, ridx, pModel->_getReac(ridx));
    	assert (reacdef != 0);
    	pReacdefs.push_back(reacdef);
    }

    uint nsreacs = pModel->_countSReacs();
    for (uint sridx = 0; sridx < nsreacs; ++sridx)
    {
      	ssolver::SReacdef * sreacdef = new SReacdef(this, sridx, pModel->_getSReac(sridx));
       	assert (sreacdef != 0);
       	pSReacdefs.push_back(sreacdef);
    }

    uint ndiffs = pModel->_countDiffs();
    for (uint didx = 0; didx < ndiffs; ++didx)
    {
       	ssolver::Diffdef * diffdef = new Diffdef(this, didx, pModel->_getDiff(didx));
       	assert (diffdef != 0);
       	pDiffdefs.push_back(diffdef);
    }

    uint ncomps = pGeom->_countComps();
    assert(ncomps >0);
    for (uint cidx = 0; cidx < ncomps; ++cidx)
    {
    	ssolver::Compdef * compdef = new Compdef(this, cidx, pGeom->_getComp(cidx));
    	assert (compdef != 0);
    	pCompdefs.push_back(compdef);
    }

    uint npatches = pGeom->_countPatches();
    for (uint pidx = 0; pidx < npatches; ++pidx)
    {
    	ssolver::Patchdef * patchdef = new Patchdef(this, pidx, pGeom->_getPatch(pidx));
    	assert (patchdef != 0);
    	pPatchdefs.push_back(patchdef);
    }

    // Now setup all the def objects. This can't be achieved purely with		/// check the order
    // the constructors, e.g.  a patch may need to add species from its
    // surface reactions to inner, outer comp
    for (SpecdefPVecI s = pSpecdefs.begin(); s != pSpecdefs.end(); ++s)
    	(*s)->setup();
    for (ReacdefPVecI r = pReacdefs.begin(); r != pReacdefs.end(); ++r)
    	(*r)->setup();
    for (SReacdefPVecI sr = pSReacdefs.begin(); sr != pSReacdefs.end(); ++sr)
    	(*sr)->setup();
    for (DiffdefPVecI d = pDiffdefs.begin(); d != pDiffdefs.end(); ++d)
    	(*d)->setup();
    for (CompdefPVecI c = pCompdefs.begin(); c != pCompdefs.end(); ++c)
    	(*c)->setup_references();
    for (PatchdefPVecI p = pPatchdefs.begin(); p != pPatchdefs.end(); ++p)
    	(*p)->setup_references();
    // Make local indices for species, (surface) reactions, diffusion rules
    // in compartments then patches. Separate method necessary since e.g.
    // Patchdef::setup_references can add species to Compdef
    for (CompdefPVecI c = pCompdefs.begin(); c != pCompdefs.end(); ++c)
    	(*c)->setup_indices();
    for (PatchdefPVecI p = pPatchdefs.begin(); p != pPatchdefs.end(); ++p)
    	(*p)->setup_indices();
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Statedef::~Statedef()
{
	SpecdefPVecCI s_end = pSpecdefs.end();
	for (SpecdefPVecCI s = pSpecdefs.begin(); s != s_end; ++s) delete *s;

	CompdefPVecCI c_end = pCompdefs.end();
	for (CompdefPVecCI c = pCompdefs.begin(); c != c_end; ++c) delete *c;

	PatchdefPVecCI p_end = pPatchdefs.end();
	for (PatchdefPVecCI p = pPatchdefs.begin(); p != p_end; ++p) delete *p;

	ReacdefPVecCI r_end = pReacdefs.end();
	for (ReacdefPVecCI r = pReacdefs.begin(); r != r_end; ++r) delete *r;

	SReacdefPVecCI sr_end = pSReacdefs.end();
	for (SReacdefPVecCI sr = pSReacdefs.begin(); sr != sr_end; ++sr) delete *sr;

	DiffdefPVecCI d_end = pDiffdefs.end();
	for (DiffdefPVecCI d = pDiffdefs.begin(); d != d_end; ++d) delete *d;
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Compdef * ssolver::Statedef::compdef(uint gidx) const
{
    assert(gidx < pCompdefs.size());
    return pCompdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getCompIdx(std::string const & c) const
{
	uint maxcidx = pCompdefs.size();
	assert (maxcidx > 0);
	assert (maxcidx == pGeom->_countComps());
	uint cidx = 0;
    while(cidx < maxcidx)
    {
    	if (c == pGeom->_getComp(cidx)->getID()) return cidx;
    	++cidx;
    }
    std::ostringstream os;
    os << "Geometry does not contain comp with string identifier '" << c << "'.";
    throw steps::ArgErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getCompIdx(steps::wm::Comp * comp) const
{
	uint maxcidx = pCompdefs.size();
	assert (maxcidx > 0);
	assert (maxcidx == pGeom->_countComps());
	uint cidx = 0;
    while(cidx < maxcidx)
    {
    	if (comp == pGeom->_getComp(cidx)) return cidx;
    	++cidx;
    }
    // Argument should be valid so we should not get here
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Patchdef * ssolver::Statedef::patchdef(uint gidx) const
{
	assert(gidx < pPatchdefs.size());
	return pPatchdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getPatchIdx(std::string const & p) const
{
	uint maxpidx = pPatchdefs.size();
	assert (maxpidx == pGeom->_countPatches());
	uint pidx = 0;
	while(pidx < maxpidx)
	{
		if (p == pGeom->_getPatch(pidx)->getID()) return pidx;
		++pidx;
	}
	std::ostringstream os;
	os << "Geometry does not contain patch with string identifier '" << p << "'.";
	throw steps::ArgErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getPatchIdx(steps::wm::Patch * patch) const
{
	uint maxpidx = pPatchdefs.size();
	assert (maxpidx == pGeom->_countPatches());
	uint pidx = 0;
	while(pidx < maxpidx)
	{
		if (patch == pGeom->_getPatch(pidx)) return pidx;
		++pidx;
	}
    // Argument should be valid so we should not get here
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Specdef * ssolver::Statedef::specdef(uint gidx) const
{
	assert(gidx < pSpecdefs.size());
	return pSpecdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getSpecIdx(std::string const & s) const
{
	uint maxsidx = pSpecdefs.size();
	assert (maxsidx > 0);
	assert(maxsidx == pModel->_countSpecs());
	uint sidx = 0;
	while(sidx < maxsidx)
	{
		if (s == pModel->_getSpec(sidx)->getID()) return sidx;
		++sidx;
	}
	std::ostringstream os;
	os << "Model does not contain species with string identifier '" << s << "'.";
	throw steps::ArgErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getSpecIdx(steps::model::Spec * spec) const
{
	uint maxsidx = pSpecdefs.size();
	assert (maxsidx > 0);
	assert(maxsidx == pModel->_countSpecs());
	uint sidx = 0;
	while(sidx < maxsidx)
	{
		if (spec == pModel->_getSpec(sidx)) return sidx;
		++sidx;
	}
    // Argument should be valid so we should not get here
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Reacdef * ssolver::Statedef::reacdef(uint gidx) const
{
	assert(gidx < pReacdefs.size());
	return pReacdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getReacIdx(std::string const & r) const
{
	uint maxridx = pReacdefs.size();
	assert (maxridx == pModel->_countReacs());
	uint ridx = 0;
	while(ridx < maxridx)
	{
		if (r == pModel->_getReac(ridx)->getID()) return ridx;
		++ridx;
	}
	std::ostringstream os;
	os << "Model does not contain reac with string identifier '" << r <<"'.";
	throw steps::ArgErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getReacIdx(steps::model::Reac * reac) const
{
	uint maxridx = pReacdefs.size();
	assert (maxridx == pModel->_countReacs());
	uint ridx = 0;
	while(ridx < maxridx)
	{
		if (reac == pModel->_getReac(ridx)) return ridx;
		++ridx;
	}
    // Argument should be valid so we should not get here
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::SReacdef * ssolver::Statedef::sreacdef(uint gidx) const
{
	assert(gidx < pSReacdefs.size());
	return pSReacdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getSReacIdx(std::string const & sr) const
{
	uint maxsridx = pSReacdefs.size();
	assert (maxsridx == pModel->_countSReacs());
	uint sridx = 0;
	while(sridx < maxsridx)
	{
		if (sr == pModel->_getSReac(sridx)->getID()) return sridx;
		++sridx;
	}
	std::ostringstream os;
	os << "Model does not contain sreac with string identifier '" << sr <<"'.";
	throw steps::ArgErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getSReacIdx(steps::model::SReac * sreac) const
{
	uint maxsridx = pSReacdefs.size();
	assert (maxsridx == pModel->_countSReacs());
	uint sridx = 0;
	while(sridx < maxsridx)
	{
		if (sreac == pModel->_getSReac(sridx)) return sridx;
		++sridx;
	}
    // Argument should be valid so we should not get here
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Diffdef * ssolver::Statedef::diffdef(uint gidx) const
{
	assert(gidx < pDiffdefs.size());
	return pDiffdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getDiffIdx(std::string const & d) const
{
	uint maxdidx = pDiffdefs.size();
	assert (maxdidx == pModel->_countDiffs());
	uint didx = 0;
	while(didx < maxdidx)
	{
		if (d == pModel->_getDiff(didx)->getID()) return didx;
		++didx;
	}
	std::ostringstream os;
	os << "Model does not contain diff with string identifier '" << d <<"'.";
	throw steps::ArgErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getDiffIdx(steps::model::Diff * diff) const
{
	uint maxdidx = pDiffdefs.size();
	assert (maxdidx == pModel->_countDiffs());
	uint didx = 0;
	while(didx < maxdidx)
	{
		if (diff == pModel->_getDiff(didx)) return didx;
		++didx;
	}
    // Argument should be valid so we should not get here
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Statedef::setTime(double t)
{
	assert (t >= 0.0);
	pTime = t;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Statedef::incTime(double dt)
{
	assert (dt >= 0.0);
	pTime += dt;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Statedef::incNSteps(uint i)
{
	assert (i != 0);
	pNSteps += i;
}

////////////////////////////////////////////////////////////////////////////////

// END
