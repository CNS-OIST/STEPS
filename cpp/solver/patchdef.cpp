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

// STEPS headers.
#include "../common.h"
#include "types.hpp"
#include "../error.hpp"
#include "statedef.hpp"
#include "patchdef.hpp"
#include "sreacdef.hpp"
#include "compdef.hpp"
#include "../geom/patch.hpp"
#include "../model/sreac.hpp"

NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

ssolver::Patchdef::Patchdef(Statedef * sd, uint idx, steps::wm::Patch * p)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pArea()
, pPssys()
, pIcomp(0)
, pOcomp(0)
, pInner(0)
, pOuter(0)
, pSetupRefsdone(false)
, pSetupIndsdone(false)
, pSpecsN_I(0)
, pSpecsN_S(0)
, pSpecsN_O(0)
, pSpec_G2L(0)
, pSpec_L2G(0)
, pPoolCount(0)
, pPoolFlags(0)
, pSReacsN(0)
, pSReac_G2L(0)
, pSReac_L2G(0)
, pSReacKcst(0)
, pSReacFlags(0)
, pSReac_DEP_I_Spec(0)
, pSReac_DEP_S_Spec(0)
, pSReac_DEP_O_Spec(0)
, pSReac_LHS_I_Spec(0)
, pSReac_LHS_S_Spec(0)
, pSReac_LHS_O_Spec(0)
, pSReac_UPD_I_Spec(0)
, pSReac_UPD_S_Spec(0)
, pSReac_UPD_O_Spec(0)
{
	assert(pStatedef != 0);
	assert(p != 0);

	pName = p->getID();
	pArea = p->getArea();
	pPssys = p->getSurfsys();
	pIcomp = p->getIComp();
	pOcomp = p->getOComp();

	uint nspecs = pStatedef->countSpecs();
	pSpec_G2L = new uint[nspecs];
	std::fill_n(pSpec_G2L, nspecs, LIDX_UNDEFINED);
    uint nsreacs = pStatedef->countSReacs();
    pSReac_G2L = new lidxT[nsreacs];
    std::fill_n(pSReac_G2L, nsreacs, LIDX_UNDEFINED);

}

////////////////////////////////////////////////////////////////////////////////

ssolver::Patchdef::~Patchdef(void)
{
	delete[] pSpec_G2L;
	delete[] pSReac_G2L;
	delete[] pSpec_L2G;
	delete[] pSReac_L2G;
	delete[] pSReac_DEP_S_Spec;
	delete[] pSReac_LHS_S_Spec;
	delete[] pSReac_UPD_S_Spec;
	delete[] pSReac_DEP_I_Spec;
	delete[] pSReac_LHS_I_Spec;
	delete[] pSReac_UPD_I_Spec;
	delete[] pSReac_DEP_O_Spec;
	delete[] pSReac_LHS_O_Spec;
	delete[] pSReac_UPD_O_Spec;
	delete[] pPoolCount;
	delete[] pPoolFlags;
	delete[] pSReacFlags;
	delete[] pSReacKcst;

}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::setup_references(void)
{
	assert(pSetupRefsdone == false);
	assert(pSetupIndsdone == false);

	// first find the inner and outer comps of this patch
	assert (pIcomp != 0);
	uint icompidx = pStatedef->getCompIdx(pIcomp);				///// bit long-winded, add new method to statedef??
	pInner = pStatedef->compdef(icompidx);
	if (pOcomp != 0)
	{
		uint ocompidx = pStatedef->getCompIdx(pOcomp);
		pOuter = pStatedef->compdef(ocompidx);
	}

	uint ngspecs = pStatedef->countSpecs();
	uint ngsreacs = pStatedef->countSReacs();

	// set up local sreac indices
	std::set<std::string>::const_iterator s_end = pPssys.end();
	for(std::set<std::string>::const_iterator s = pPssys.begin();
		s != s_end; ++s)
	{
		std::map<std::string, steps::model::SReac *> ssreacs = pStatedef->model()->getSurfsys(*s)->_getAllSReacs();
		std::map<std::string, steps::model::SReac*>::const_iterator sr_end = ssreacs.end();
		for(std::map<std::string, steps::model::SReac *>::const_iterator sr = ssreacs.begin(); sr != sr_end; ++sr)
		{
			uint gidx = pStatedef->getSReacIdx((sr->second));
			assert(gidx < ngsreacs);
			if(pSReac_G2L[gidx] != LIDX_UNDEFINED) continue;
			pSReac_G2L[gidx] = pSReacsN++;
		}
	}

	// Now add all species that appear in all surface reactions that can occur
	// on this patch: to the patch, inner or outer compartment.
	for(uint sr = 0; sr < ngsreacs; ++sr)
	{
		if(pSReac_G2L[sr] == LIDX_UNDEFINED) continue;
		SReacdef * srdef = pStatedef->sreacdef(sr);
		assert(srdef != 0);
		for (uint s = 0; s < ngspecs; ++s)
		{
			if (srdef->reqspec_S(s) == true)
			{
				assert (pStatedef->specdef(s) != 0);
				if (pSpec_G2L[s] == LIDX_UNDEFINED) pSpec_G2L[s] = pSpecsN_S++;
			}
			if (srdef->reqspec_I(s) == true)
			{
				assert(pInner != 0);
				pInner->addSpec(s);
			}
			if(srdef->reqspec_O(s) == true)
			{
				assert(pOuter != 0);
				pOuter->addSpec(s);
			}
		}
	}

	pSetupRefsdone = true;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::setup_indices(void)
{
	assert(pSetupRefsdone == true);
	assert(pSetupIndsdone == false);

    // 1 -- DEAL WITH PATCH SPECIES
    // (Only if any species have been added to the patch)
    //   -> Setup local indices for all species.
    //
    // 2 -- COPY #SPECIES FOR INNER AND OUTER COMPS
    //   (These are required a lot during simulation, so it's sound
    //   to have them ready here to avoid an extra level of pointer
    //   lookup.)
    //
    // 3 -- DEAL WITH PATCH SREAC'S
    // (Only if any surface reactions have been added to the patch)
    //   -> Setup local indices for all surface reactions.
    //   -> The SReac objects have LHS, DEP and UPD vectors expressed in
    //      global species indices. Transform this to local indices:
    //      -> Pre-create pSReac_DEP, _LHS and _UPD vectors
    //          -> Always for surface (_S)
    //          -> For inner comp (_I) if inner comp has been defined
    //          -> For outer comp (_O) if outer comp has been defined
    //      -> Loop over the SReacDef objects added to this patch:
    //          -> Fill out the newly created vectors by appropriate
    //             copying of the vectors defined the SReacDef object.
    //          -> While doing this, check whether everything can be
    //             resolved.

    // 1 -- DEAL WITH PATCH SPECIES
    uint ngspecs = pStatedef->countSpecs();
    if (pSpecsN_S != 0)
    {
        pSpec_L2G = new uint[pSpecsN_S];
        for (uint i = 0; i < ngspecs; ++i)
        {
            uint lidx = pSpec_G2L[i];
            if (lidx == LIDX_UNDEFINED) continue;
            pSpec_L2G[lidx] = i;
        }
    }

    // 2 -- COPY #SPECS FOR INNER AND OUTER COMPS
    if (pInner != 0) pSpecsN_I = pInner->countSpecs();
    if (pOuter != 0) pSpecsN_O = pOuter->countSpecs();

    // 3 -- DEAL WITH PATCH SREAC'S
    if (pSReacsN != 0)
    {
        // Set up local indices.
        pSReac_L2G = new uint[pSReacsN];
        uint ngsreacs = pStatedef->countSReacs();
        for (uint i = 0; i < ngsreacs; ++i)
        {
            uint lidx = pSReac_G2L[i];
            if (lidx == LIDX_UNDEFINED) continue;
            pSReac_L2G[lidx] = i;
        }

        // Create _DEP, _LHS and _UPD vectors.
        uint arrsize_i = 0;
        uint arrsize_s = pSpecsN_S * pSReacsN;
        uint arrsize_o = 0;
        pSReac_DEP_S_Spec = new int[arrsize_s];
        pSReac_LHS_S_Spec = new uint[arrsize_s];
        pSReac_UPD_S_Spec = new int[arrsize_s];
        std::fill_n(pSReac_DEP_S_Spec, arrsize_s, 0);
        std::fill_n(pSReac_LHS_S_Spec, arrsize_s, 0);
        std::fill_n(pSReac_UPD_S_Spec, arrsize_s, 0);
        if (pInner != 0) // Only create if inner comp exists.
        {
            arrsize_i = pSpecsN_I * pSReacsN;
            pSReac_DEP_I_Spec = new int[arrsize_i];
            pSReac_LHS_I_Spec = new uint[arrsize_i];
            pSReac_UPD_I_Spec = new int[arrsize_i];
            std::fill_n(pSReac_DEP_I_Spec, arrsize_i, 0);
            std::fill_n(pSReac_LHS_I_Spec, arrsize_i, 0);
            std::fill_n(pSReac_UPD_I_Spec, arrsize_i, 0);
        }
        if (pOuter != 0) // Only create if outer comp exists.
        {
            arrsize_o = pSpecsN_O * pSReacsN;
            pSReac_DEP_O_Spec = new int[arrsize_o];
            pSReac_LHS_O_Spec = new uint[arrsize_o];
            pSReac_UPD_O_Spec = new int[arrsize_o];
            std::fill_n(pSReac_DEP_O_Spec, arrsize_o, 0);
            std::fill_n(pSReac_LHS_O_Spec, arrsize_o, 0);
            std::fill_n(pSReac_UPD_O_Spec, arrsize_o, 0);
        }

        // Fill the vectors with all kinds of useful information.
        for (uint ri = 0; ri < pSReacsN; ++ri)
        {
            SReacdef * srdef = sreacdef(ri);

            // Handle surface stuff.
            for (uint si = 0; si < ngspecs; ++si)
            {
                if (srdef->reqspec_S(si) == false) continue;

                // TODO: turn into error check?
                uint sil = pSpec_G2L[si];
                assert(sil != LIDX_UNDEFINED);

                uint aridx = _IDX_SReac_S_Spec(ri, sil);
                pSReac_DEP_S_Spec[aridx] = srdef->dep_S(si);
                pSReac_LHS_S_Spec[aridx] = srdef->lhs_S(si);
                pSReac_UPD_S_Spec[aridx] = srdef->upd_S(si);
            }

            // Handle the inside comp stuff.
            if (srdef->reqInside() == true)
            {
                // TODO: turn into real error check?
                assert(pInner != 0);

                for (uint si = 0; si < ngspecs; ++si)
                {
                    if (srdef->reqspec_I(si) == false) continue;

                    // TODO: turn into error check?
                    uint sil = specG2L_I(si);
                    assert(sil != LIDX_UNDEFINED);

                    uint aridx = _IDX_SReac_I_Spec(ri, sil);
                    pSReac_DEP_I_Spec[aridx] = srdef->dep_I(si);
                    pSReac_LHS_I_Spec[aridx] = srdef->lhs_I(si);
                    pSReac_UPD_I_Spec[aridx] = srdef->upd_I(si);
                }
            }

            // Handle the outside comp stuff.
            if (srdef->reqOutside() == true)
            {
                // TODO: turn into real error check?
                assert(pOuter != 0);

                for (uint si = 0; si < ngspecs; ++si)
                {
                    if (srdef->reqspec_O(si) == false) continue;

                    // TODO: turn into error check?
                    uint sil = specG2L_O(si);
                    assert(sil != LIDX_UNDEFINED);

                    uint aridx = _IDX_SReac_O_Spec(ri, sil);
                    pSReac_DEP_O_Spec[aridx] = srdef->dep_O(si);
                    pSReac_LHS_O_Spec[aridx] = srdef->lhs_O(si);
                    pSReac_UPD_O_Spec[aridx] = srdef->upd_O(si);
                }
            }
        }
    }

    // Initialise the pools and flags members to zeros.
    pPoolCount = new double[pSpecsN_S];
    pPoolFlags = new uint[pSpecsN_S];
    pSReacFlags = new uint[pSReacsN];
    std::fill_n(pPoolCount, pSpecsN_S, 0.0);
    std::fill_n(pPoolFlags, pSpecsN_S, 0);
	std::fill_n(pSReacFlags, pSReacsN, 0);

	// Finally initialise Kcsts to user-supplied values
    pSReacKcst = new double[pSReacsN];

	for (uint i = 0; i < pSReacsN; ++i)
	{
		// sreacdef() returns global Reacdef by local index
		ssolver::SReacdef * sreac = sreacdef(i);
		pSReacKcst[i] = sreac->kcst();
	}

	pSetupIndsdone = true;
}

////////////////////////////////////////////////////////////////////////////////

double ssolver::Patchdef::area(void) const
{
	return pArea;
}

////////////////////////////////////////////////////////////////////////////////

std::string const ssolver::Patchdef::name(void) const
{
	return pName;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::setArea(double a)
{
	assert (a > 0.0);
	pArea = a;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::reset(void)
{
	assert(pSetupRefsdone == true);
	assert(pSetupIndsdone == true);
    std::fill_n(pPoolCount, pSpecsN_S, 0.0);
    std::fill_n(pPoolFlags, pSpecsN_S, 0);
	std::fill_n(pSReacFlags, pSReacsN, 0);

	for (uint i = 0; i < pSReacsN; ++i)
	{
		// sreacdef() returns global Reacdef by local index
		ssolver::SReacdef * sreac = sreacdef(i);
		pSReacKcst[i] = sreac->kcst();
	}
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::setCount(uint slidx, double count)
{
	assert(pSetupRefsdone == true);
	assert(pSetupIndsdone == true);
	assert(slidx < pSpecsN_S);
	assert (count >= 0.0);
	pPoolCount[slidx] = count;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::setClamped(uint slidx, bool clamp)
{
	assert(pSetupRefsdone == true);
	assert(pSetupIndsdone == true);
	assert(slidx < pSpecsN_S);
    if (clamp == true) pPoolFlags[slidx] |= CLAMPED;
    else pPoolFlags[slidx] &= ~CLAMPED;
}

////////////////////////////////////////////////////////////////////////////////

ssolver::SReacdef * ssolver::Patchdef::sreacdef(uint lidx) const
{
	assert(pSetupRefsdone == true);
	assert (lidx < pSReacsN);
	return pStatedef->sreacdef(pSReac_L2G[lidx]);
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Patchdef::sreac_dep_I(uint srlidx, uint splidx) const
{
    return pSReac_DEP_I_Spec[splidx + (srlidx * pSpecsN_I)];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Patchdef::sreac_dep_S(uint srlidx, uint splidx) const
{
    return pSReac_DEP_S_Spec[splidx + (srlidx * pSpecsN_S)];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Patchdef::sreac_dep_O(uint srlidx, uint splidx) const
{
    return pSReac_DEP_O_Spec[splidx + (srlidx * pSpecsN_O)];
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::sreac_lhs_I_bgn(uint lidx) const
{
    return pSReac_LHS_I_Spec + (lidx * pSpecsN_I);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::sreac_lhs_I_end(uint lidx) const
{
    return pSReac_LHS_I_Spec + ((lidx + 1) * pSpecsN_I);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::sreac_lhs_S_bgn(uint lidx) const
{
    return pSReac_LHS_S_Spec + (lidx * pSpecsN_S);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::sreac_lhs_S_end(uint lidx) const
{
    return pSReac_LHS_S_Spec + ((lidx + 1) * pSpecsN_S);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::sreac_lhs_O_bgn(uint lidx) const
{
    return pSReac_LHS_O_Spec + (lidx * pSpecsN_O);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::sreac_lhs_O_end(uint lidx) const
{
    return pSReac_LHS_O_Spec + ((lidx + 1) * pSpecsN_O);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::sreac_upd_I_bgn(uint lidx) const
{
    return pSReac_UPD_I_Spec + (lidx * pSpecsN_I);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::sreac_upd_I_end(uint lidx) const
{
    return pSReac_UPD_I_Spec + ((lidx + 1) * pSpecsN_I);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::sreac_upd_S_bgn(uint lidx) const
{
    return pSReac_UPD_S_Spec + (lidx * pSpecsN_S);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::sreac_upd_S_end(uint lidx) const
{
    return pSReac_UPD_S_Spec + ((lidx + 1) * pSpecsN_S);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::sreac_upd_O_bgn(uint lidx) const
{
    return pSReac_UPD_O_Spec + (lidx * pSpecsN_O);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::sreac_upd_O_end(uint lidx) const
{
    return pSReac_UPD_O_Spec + ((lidx + 1) * pSpecsN_O);
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::setKcst(uint srlidx, double kcst)
{
	assert(pSetupRefsdone == true);
	assert(pSetupIndsdone == true);
	assert(srlidx < pSReacsN);
	assert (kcst >= 0.0);
	pSReacKcst[srlidx] = kcst;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::setActive(uint srlidx, bool active)
{
	assert(pSetupRefsdone == true);
	assert(pSetupIndsdone == true);
	assert(srlidx < pSReacsN);
	if (active == true) pSReacFlags[srlidx] &= ~INACTIVATED;
	else pSReacFlags[srlidx] |= INACTIVATED;
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Patchdef::specG2L_I(uint gidx) const
{
    if (pInner == 0) return LIDX_UNDEFINED;
    return pInner->specG2L(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Patchdef::specG2L_O(uint gidx) const
{
    if (pOuter == 0) return LIDX_UNDEFINED;
    return pOuter->specG2L(gidx);
}

////////////////////////////////////////////////////////////////////////////////

// END
