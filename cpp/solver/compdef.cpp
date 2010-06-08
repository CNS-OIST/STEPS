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
#include <cassert>
#include <algorithm>

// STEPS headers.
#include "../common.h"
#include "types.hpp"
#include "../error.hpp"
#include "statedef.hpp"
#include "compdef.hpp"
#include "patchdef.hpp"
#include "reacdef.hpp"
#include "diffdef.hpp"
#include "../geom/comp.hpp"
#include "../model/reac.hpp"
#include "../model/diff.hpp"

NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

ssolver::Compdef::Compdef(Statedef * sd, uint idx, steps::wm::Comp * c)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pVol()
, pCvsys()
, pSetupRefsdone(false)
, pSetupIndsdone(false)
, pPoolFlags(0)
, pPoolCount(0)
, pReacKcst(0)
, pDiffDcst(0)
, pReacFlags(0)
, pSpecsN(0)
, pSpec_G2L(0)
, pSpec_L2G(0)
, pReacsN(0)
, pReac_G2L(0)
, pReac_L2G(0)
, pReac_DEP_Spec(0)
, pReac_LHS_Spec(0)
, pReac_UPD_Spec(0)
, pDiffsN(0)
, pDiff_G2L(0)
, pDiff_L2G(0)
, pDiff_DEP_Spec(0)
, pDiff_LIG(0)
{
    assert(pStatedef != 0);
    assert(c != 0);

    pName = c->getID();
    pVol = c->getVol();
    pCvsys = c->getVolsys();

    uint nspecs = pStatedef->countSpecs();
    pSpec_G2L = new uint[nspecs];
    std::fill_n(pSpec_G2L, nspecs, LIDX_UNDEFINED);
    uint nreacs = pStatedef->countReacs();
    pReac_G2L = new uint[nreacs];
    std::fill_n(pReac_G2L, nreacs, LIDX_UNDEFINED);
    uint ndiffs = pStatedef->countDiffs();
    pDiff_G2L = new uint[ndiffs];
    std::fill_n(pDiff_G2L, ndiffs, LIDX_UNDEFINED);

}

////////////////////////////////////////////////////////////////////////////////

ssolver::Compdef::~Compdef(void)
{
	delete[] pSpec_G2L;
	delete[] pReac_G2L;
	delete[] pDiff_G2L;
	delete[] pSpec_L2G;
	delete[] pReac_L2G;
	delete[] pDiff_L2G;
	delete[] pReac_DEP_Spec;
	delete[] pReac_LHS_Spec;
	delete[] pReac_UPD_Spec;
	delete[] pDiff_DEP_Spec;
	delete[] pDiff_LIG;
	delete[] pPoolFlags;
	delete[] pPoolCount;
	delete[] pReacKcst;
	delete[] pReacFlags;
	delete[] pDiffDcst;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Compdef::setup_references()
{
	assert (pSetupRefsdone == false);
	assert (pSetupIndsdone == false);

	uint ngspecs = pStatedef->countSpecs();
    uint ngreacs = pStatedef->countReacs();
    uint ngdiffs = pStatedef->countDiffs();

    ////////############ this currently assumes that all species from patch sreacs have
    /////// been added first. Statedef calls setup on patches, which add Specs to their inner
    ///// and outer compartments.

	// set up local reac indices ////vsys also has _countReacs and Reac * _getReac(lidx)
    // The local reacs index count (pReacsN) starts at 0.
    std::set<std::string>::const_iterator v_end = pCvsys.end();
    for(std::set<std::string>::const_iterator v = pCvsys.begin();
		v != v_end; ++v)
    {
    	std::map<std::string, steps::model::Reac *> vreacs = pStatedef->model()->getVolsys(*v)->_getAllReacs();
		std::map<std::string, steps::model::Reac*>::const_iterator r_end = vreacs.end();
       	for (std::map<std::string, steps::model::Reac*>::const_iterator r = vreacs.begin(); r != r_end; ++r)
       	{
       		uint gidx = pStatedef->getReacIdx((r->second));
       		assert(gidx < ngreacs);
       		if (pReac_G2L[gidx] != LIDX_UNDEFINED) continue;
       		pReac_G2L[gidx] = pReacsN++;
      	}
       	std::map<std::string, steps::model::Diff *> vdiffs = pStatedef->model()->getVolsys(*v)->_getAllDiffs();
       	std::map<std::string, steps::model::Diff*>::const_iterator d_end = vdiffs.end();
       	for (std::map<std::string, steps::model::Diff*>::const_iterator d = vdiffs.begin(); d != d_end; ++d)
       	{
       		uint gidx = pStatedef->getDiffIdx((d->second));
       		assert(gidx < ngdiffs);
       		if (pDiff_G2L[gidx] != LIDX_UNDEFINED) continue;
       		pDiff_G2L[gidx] = pDiffsN++;
       	}
    }

    // now add all species that appear in all reactions, diffusions that
    // can occur in this compartment
    // NOTE: Patchdef setups have called addSpec() to already add some
    // species (from surface reactions)
    for (uint r = 0; r < ngreacs; ++r)
    {
    	if (pReac_G2L[r] == LIDX_UNDEFINED) continue;
    	Reacdef * rdef = pStatedef->reacdef(r);
    	assert (rdef != 0);
    	for (uint s = 0; s < ngspecs; ++s)
    	{
    		if (rdef->reqspec(s) == true) addSpec(s);
    	}
    }
    for (uint d = 0; d < ngdiffs; ++d)
    {
    	if (pDiff_G2L[d] == LIDX_UNDEFINED) continue;
    	Diffdef * ddef = pStatedef->diffdef(d);
    	assert (ddef != 0);
    	for (uint s = 0; s < ngspecs; ++s)
    	{
    		if (ddef->reqspec(s) == true) addSpec(s);
    	}
    }

    pSetupRefsdone = true;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Compdef::setup_indices(void)
{
	assert (pSetupRefsdone == true);
	assert (pSetupIndsdone == false);

	uint ngspecs = pStatedef->countSpecs();
    uint ngreacs = pStatedef->countReacs();
    uint ngdiffs = pStatedef->countDiffs();

    // Set up local indices
    if (pSpecsN != 0)
    {
    	pSpec_L2G = new uint[pSpecsN];
    	for (uint i=0; i < ngspecs; ++i)
    	{
    		uint lidx = pSpec_G2L[i];
    		if (lidx ==  LIDX_UNDEFINED) continue;
    		pSpec_L2G[lidx] = i;
    	}
    }

    if (pReacsN != 0)
    {
    	pReac_L2G = new uint[pReacsN];
    	for (uint i=0; i < ngreacs; ++i)
    	{
    		uint lidx = pReac_G2L[i];
    		if (lidx == LIDX_UNDEFINED) continue;
    		pReac_L2G[lidx] = i;
    	}
    	uint arrsize = pSpecsN * pReacsN;
    	pReac_DEP_Spec = new int[arrsize];
    	pReac_LHS_Spec = new uint[arrsize];
    	pReac_UPD_Spec = new int[arrsize];
        std::fill_n(pReac_DEP_Spec, arrsize, 0);
        std::fill_n(pReac_LHS_Spec, arrsize, 0);
        std::fill_n(pReac_UPD_Spec, arrsize, 0);
        for(uint ri = 0; ri < pReacsN; ++ri)
        {
        	Reacdef * rdef = reacdef(ri);
        	for (uint si = 0; si < ngspecs; ++si)
        	{
        		if (rdef->reqspec(si) == false) continue;
        		uint sil = pSpec_G2L[si];
        		assert(sil != LIDX_UNDEFINED);

        		uint aridx = _IDX_Reac_Spec(ri, sil);
        		pReac_DEP_Spec[aridx] = rdef->dep(si);
        		pReac_LHS_Spec[aridx] = rdef->lhs(si);
        		pReac_UPD_Spec[aridx] = rdef->upd(si);
        	}
        }
    }

    if (pDiffsN != 0)
    {
    	pDiff_L2G = new uint[pDiffsN];
    	for (uint i = 0; i < ngdiffs; ++i)
    	{
    		uint lidx = pDiff_G2L[i];
    		if (lidx == LIDX_UNDEFINED) continue;
    		pDiff_L2G[lidx] = i;
    	}

    	uint arrsize = pSpecsN * pDiffsN;
    	pDiff_DEP_Spec = new uint[arrsize];
    	std::fill_n(pDiff_DEP_Spec, arrsize, 0);
    	pDiff_LIG = new uint[pDiffsN];
    	for (uint di = 0; di < pDiffsN; ++di)
    	{
    		Diffdef * ddef = diffdef(di);
    		pDiff_LIG[di] = pSpec_G2L[ddef->lig()];
    		for (uint si = 0; si < ngspecs; ++si)
    		{
    			if (ddef->reqspec(si) == false) continue;
    			uint sil = pSpec_G2L[si];
    			assert(sil != LIDX_UNDEFINED);
    			uint aridx = _IDX_Diff_Spec(di, sil);
    			pDiff_DEP_Spec[aridx] = ddef->dep(si);
    		}
    	}
    }

    // Initialise the pools and flags members to zeros.
    pPoolCount = new double[pSpecsN];
    pPoolFlags = new uint[pSpecsN];
    pReacFlags = new uint[pReacsN];
    std::fill_n(pPoolCount, pSpecsN, 0.0);
    std::fill_n(pPoolFlags, pSpecsN, 0);
	std::fill_n(pReacFlags, pReacsN, 0);

	// Finally initialise constants to user-supplied values
    pReacKcst = new double[pReacsN];

	for (uint i = 0; i < pReacsN; ++i)
	{
		// reacdef() returns global Reacdef by local index
		ssolver::Reacdef * reac = reacdef(i);
		pReacKcst[i] = reac->kcst();
	}

	pDiffDcst = new double[pDiffsN];

	for (uint i = 0; i <pDiffsN; ++i)
	{
		// diffdef() returns global Diffdef by local index
		ssolver::Diffdef * diff = diffdef(i);
		pDiffDcst[i] = diff->dcst();
	}

    pSetupIndsdone = true;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Compdef::addIPatchdef(ssolver::Patchdef * p)
{
    // Make some checks.
    assert(p != 0);
    assert(p->ocompdef() == this);
    // Check whether it's already included.
    ssolver::PatchDefPVecI ip_end = pIPatches.end();
    if (std::find(pIPatches.begin(), ip_end, p) != ip_end) return;
    ssolver::PatchDefPVecI op_end = pOPatches.end();
    assert(std::find(pOPatches.begin(), op_end, p) == op_end);
    // Include.
    pIPatches.push_back(p);
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Compdef::addOPatchdef(ssolver::Patchdef * p)
{
    // Make some checks.
    assert(p != 0);
    assert(p->icompdef() == this);
    // Check whether it's already included.
    ssolver::PatchDefPVecI op_end = pOPatches.end();
    if (std::find(pOPatches.begin(), op_end, p) != op_end) return;
    ssolver::PatchDefPVecI ip_end = pIPatches.end();
    assert(std::find(pIPatches.begin(), ip_end, p) == ip_end);
    // Include.
    pOPatches.push_back(p);
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Compdef::addSpec(uint gidx)
{
	assert (pSetupIndsdone == false);
	assert (pStatedef->specdef(gidx) != 0);
	if (pSpec_G2L[gidx] != LIDX_UNDEFINED) return;
	pSpec_G2L[gidx] = pSpecsN++;
}

////////////////////////////////////////////////////////////////////////////////

double ssolver::Compdef::vol(void) const
{
	return pVol;
}

////////////////////////////////////////////////////////////////////////////////

std::string const ssolver::Compdef::name(void) const
{
	return pName;

}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Compdef::setVol(double v)
{
    assert (v > 0.0);
    pVol = v;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Compdef::reset(void)
{
	assert(pSetupRefsdone == true);
	assert(pSetupIndsdone == true);
    std::fill_n(pPoolCount, pSpecsN, 0.0);
    std::fill_n(pPoolFlags, pSpecsN, 0);
	std::fill_n(pReacFlags, pReacsN, 0);
	for (uint i = 0; i < pReacsN; ++i)
	{
		ssolver::Reacdef * reac = reacdef(i);
		pReacKcst[i] = reac->kcst();
	}

	for (uint i = 0; i <pDiffsN; ++i)
	{
		ssolver::Diffdef * diff = diffdef(i);
		pDiffDcst[i] = diff->dcst();
	}
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Compdef::setCount(uint slidx, double count)
{
	assert(pSetupRefsdone == true);
	assert(pSetupIndsdone == true);
	assert(slidx < pSpecsN);
	assert (count >= 0.0);
	pPoolCount[slidx] = count;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Compdef::setClamped(uint slidx, bool clamp)
{
	assert(pSetupRefsdone == true);
	assert(pSetupIndsdone == true);
	assert(slidx < pSpecsN);
    if (clamp == true) pPoolFlags[slidx] |= CLAMPED;
    else pPoolFlags[slidx] &= ~CLAMPED;
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Compdef::reac_lhs_bgn(uint rlidx) const
{
	assert (rlidx < pReacsN);
	return pReac_LHS_Spec + (rlidx * pSpecsN);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Compdef::reac_lhs_end(uint rlidx) const
{
	assert (rlidx < pReacsN);
	return pReac_LHS_Spec + ((rlidx+1) * pSpecsN);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Compdef::reac_upd_bgn(uint rlidx) const
{
	assert (rlidx < pReacsN);
	return pReac_UPD_Spec + ((rlidx) * pSpecsN);
}

////////////////////////////////////////////////////////////////////////////////


int * ssolver::Compdef::reac_upd_end(uint rlidx) const
{
	assert (rlidx < pReacsN);
	return pReac_UPD_Spec + ((rlidx+1) * pSpecsN);
}
////////////////////////////////////////////////////////////////////////////////

int ssolver::Compdef::reac_dep(uint rlidx, uint slidx) const
{
	return pReac_DEP_Spec[slidx + ((rlidx) * pSpecsN)];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Compdef::diff_dep(uint dlidx, uint slidx) const
{
	return pDiff_DEP_Spec[slidx + ((dlidx) * pSpecsN)];
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Reacdef * ssolver::Compdef::reacdef(uint rlidx) const
{
	assert(pSetupRefsdone == true);
	assert (rlidx < pReacsN);
	return pStatedef->reacdef(pReac_L2G[rlidx]);
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Compdef::setKcst(uint rlidx, double kcst)
{
	assert(pSetupRefsdone == true);
	assert(pSetupIndsdone == true);
	assert(rlidx < pReacsN);
	assert(kcst >= 0.0);
	pReacKcst[rlidx] = kcst;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Compdef::setDcst(uint dlidx, double dcst)
{
	assert(pSetupRefsdone == true);
	assert(pSetupIndsdone == true);
	assert(dlidx < pDiffsN);
	assert(dcst >= 0.0);
	pDiffDcst[dlidx] = dcst;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Compdef::setActive(uint rlidx, bool active)
{
	assert(pSetupRefsdone == true);
	assert(pSetupIndsdone == true);
	assert(rlidx < pReacsN);
	if (active == true) pReacFlags[rlidx] &= ~INACTIVATED;
	else pReacFlags[rlidx] |= INACTIVATED;
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Diffdef * ssolver::Compdef::diffdef(uint dlidx) const
{
	assert(pSetupRefsdone == true);
	assert(dlidx < pDiffsN);
	return pStatedef->diffdef(pDiff_L2G[dlidx]);
}

////////////////////////////////////////////////////////////////////////////////

// END
