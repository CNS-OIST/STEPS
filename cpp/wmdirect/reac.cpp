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
// #include <vector>

// STEPS headers.
#include "../common.h"
#include "../math/constants.hpp"
#include "reac.hpp"
#include "comp.hpp"
#include "kproc.hpp"
#include "wmdirect.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::wmdirect, swmd);
NAMESPACE_ALIAS(steps::solver, ssolver);
NAMESPACE_ALIAS(steps::math, smath);

////////////////////////////////////////////////////////////////////////////////

static inline double comp_ccst(double kcst, double vol, uint order)
{
    double vscale = 1.0e3 * vol * smath::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;
    if (o1 < 0) o1 = 0;
    return kcst * pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

swmd::Reac::Reac(ssolver::Reacdef * rdef, swmd::Comp * comp)
: KProc()
, pReacdef(rdef)
, pComp(comp)
, pUpdVec()
, pCcst(0.0)
{
	assert (pReacdef != 0);
	assert (pComp != 0);
	uint lridx = pComp->def()->reacG2L(pReacdef->gidx());
	double kcst = pComp->def()->kcst(lridx);
	pCcst = comp_ccst(kcst, pComp->def()->vol(), pReacdef->order());
	assert (pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

swmd::Reac::~Reac(void)
{
}

////////////////////////////////////////////////////////////////////////////////

bool swmd::Reac::active(void) const
{
	uint lridx = pComp->def()->reacG2L(defr()->gidx());
	return pComp->def()->active(lridx);
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Reac::reset(void)
{
    resetExtent();
	uint lridx = pComp->def()->reacG2L(defr()->gidx());
    pComp->def()->setActive(lridx, true);
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Reac::resetCcst(void)
{
	uint lridx = pComp->def()->reacG2L(pReacdef->gidx());
	double kcst = pComp->def()->kcst(lridx);
	pCcst = comp_ccst(kcst, pComp->def()->vol(), pReacdef->order());
	assert (pCcst >= 0);
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Reac::setupDeps(void)
{
    SchedIDXSet updset;
    ssolver::gidxTVecCI sbgn = defr()->bgnUpdColl();
    ssolver::gidxTVecCI send = defr()->endUpdColl();

    // Search in local compartment.
    KProcPVecCI kprocend = pComp->kprocEnd();
    for (KProcPVecCI k = pComp->kprocBegin(); k != kprocend; ++k)
    {
        for (ssolver::gidxTVecCI s = sbgn; s != send; ++s)
        {
            if ((*k)->depSpecComp(*s, pComp) == true)
                updset.insert((*k)->schedIDX());
        }
    }

    // Search in neighbouring patches.
    for (PatchPVecCI p = pComp->beginIPatches(); p != pComp->endIPatches(); ++p)
    {
        kprocend = (*p)->kprocEnd();
        for (KProcPVecCI k = (*p)->kprocBegin(); k != kprocend; ++k)
        {
            for (ssolver::gidxTVecCI s = sbgn; s != send; ++s)
            {
                if ((*k)->depSpecComp(*s, pComp) == true)
                    updset.insert((*k)->schedIDX());
            }
        }
    }
    for (PatchPVecCI p = pComp->beginOPatches(); p != pComp->endOPatches(); ++p)
    {
        kprocend = (*p)->kprocEnd();
        for (KProcPVecCI k = (*p)->kprocBegin(); k != kprocend; ++k)
        {
            for (ssolver::gidxTVecCI s = sbgn; s != send; ++s)
            {
                if ((*k)->depSpecComp(*s, pComp) == true)
                    updset.insert((*k)->schedIDX());
            }
        }
    }

    swmd::schedIDXSet_To_Vec(updset, pUpdVec);
}

////////////////////////////////////////////////////////////////////////////////

bool swmd::Reac::depSpecComp(uint gidx, swmd::Comp * comp)
{
	if (pComp != comp) return false;
	return defr()->dep(gidx);
}

////////////////////////////////////////////////////////////////////////////////

bool swmd::Reac::depSpecPatch(uint gidx, swmd::Patch * patch)
{
	return false;
}

////////////////////////////////////////////////////////////////////////////////

double swmd::Reac::rate(void) const
{
	if (inactive()) return 0.0;

    // Prefetch some variables.
    ssolver::Compdef * cdef = pComp->def();
    uint nspecs = cdef->countSpecs();
    uint * lhs_vec = cdef->reac_lhs_bgn(cdef->reacG2L(defr()->gidx()));
    double * cnt_vec = cdef->pools();

    // Compute combinatorial part.
        double h_mu = 1.0;
        for (uint pool = 0; pool < nspecs; ++pool)
        {
            uint lhs = lhs_vec[pool];
            if (lhs == 0) continue;
            uint cnt = static_cast<uint>(cnt_vec[pool]);
            if (lhs > cnt)
            {
                h_mu = 0.0;
                break;
            }
            switch (lhs)
            {
                case 4:
                {
                    h_mu *= static_cast<double>(cnt - 3);
                }
                case 3:
                {
                    h_mu *= static_cast<double>(cnt - 2);
                }
                case 2:
                {
                    h_mu *= static_cast<double>(cnt - 1);
                }
                case 1:
                {
                    h_mu *= static_cast<double>(cnt);
                    break;
                }
                default:
                {
                    assert(0);
                    return 0.0;
                }
            }
        }

        // Multiply with scaled reaction constant.
        return h_mu * pCcst;

}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> const & swmd::Reac::apply(void)
{
    ssolver::Compdef * cdef = pComp->def();
    double * local = cdef->pools();
    uint l_ridx = cdef->reacG2L(defr()->gidx());
    int * upd_vec = cdef->reac_upd_bgn(l_ridx);
    uint nspecs = cdef->countSpecs();
    for (uint i=0; i < nspecs; ++i)
    {
    	if (cdef->clamped(i) == true) continue;
    	int j = upd_vec[i];
    	if (j == 0) continue;
    	int nc = static_cast<int>(local[i]) + j;
    	cdef->setCount(i, static_cast<double>(nc));
    }
    rExtent++;
    return pUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

// END

