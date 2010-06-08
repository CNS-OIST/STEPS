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
#include <vector>

// STEPS headers.
#include "../common.h"
#include "../math/constants.hpp"
#include "../solver/reacdef.hpp"
#include "reac.hpp"
#include "tet.hpp"
#include "kproc.hpp"
#include "tetexact.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::tetexact, stex);
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

stex::Reac::Reac(ssolver::Reacdef * rdef, stex::Tet * tet)
: KProc()
, pReacdef(rdef)
, pTet(tet)
, pUpdVec()
, pCcst(0.0)
, pKcst(0.0)
{
	assert (pReacdef != 0);
	assert (pTet != 0);

	uint lridx = pTet->compdef()->reacG2L(pReacdef->gidx());
	double kcst = pTet->compdef()->kcst(lridx);
	pKcst = kcst;
	pCcst = comp_ccst(kcst, pTet->vol(), pReacdef->order());
	assert (pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

stex::Reac::~Reac(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void stex::Reac::reset(void)
{
    resetExtent();
    resetCcst();
	setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Reac::resetCcst(void)
{
	uint lridx = pTet->compdef()->reacG2L(pReacdef->gidx());
	double kcst = pTet->compdef()->kcst(lridx);
	// Also reset kcst
	pKcst = kcst;
	pCcst = comp_ccst(kcst, pTet->vol(), pReacdef->order());
	assert (pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Reac::setKcst(double k)
{
	assert (k >= 0.0);
	pKcst = k;
	pCcst = comp_ccst(k, pTet->vol(), pReacdef->order());
	assert (pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Reac::setupDeps(void)
{
    SchedIDXSet updset;
    ssolver::gidxTVecCI sbgn = pReacdef->bgnUpdColl();
    ssolver::gidxTVecCI send = pReacdef->endUpdColl();

    // Search in local tetrahedron.
    KProcPVecCI kprocend = pTet->kprocEnd();
    for (KProcPVecCI k = pTet->kprocBegin(); k != kprocend; ++k)
    {
        for (ssolver::gidxTVecCI s = sbgn; s != send; ++s)
        {
            if ((*k)->depSpecTet(*s, pTet) == true)
                updset.insert((*k)->schedIDX());
        }
    }

    // Search in neighbouring triangles.
    for (uint i = 0; i < 4; ++i)
    {
        // Fetch next triangle, if it exists.
        stex::Tri * next = pTet->nextTri(i);
        if (next == 0) continue;

        kprocend = next->kprocEnd();
        for (KProcPVecCI k = next->kprocBegin(); k != kprocend; ++k)
        {
            for (ssolver::gidxTVecCI s = sbgn; s != send; ++s)
            {
                if ((*k)->depSpecTet(*s, pTet) == true)
                    updset.insert((*k)->schedIDX());
            }
        }
    }

    schedIDXSet_To_Vec(updset, pUpdVec);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Reac::depSpecTet(uint gidx, stex::Tet * tet)
{
    if (pTet != tet) return false;
    return pReacdef->dep(gidx);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Reac::depSpecTri(uint gidx, stex::Tri * tri)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Reac::rate(void) const
{
    if (inactive()) return 0.0;

    // Prefetch some variables.
    ssolver::Compdef * cdef = pTet->compdef();
    uint nspecs = cdef->countSpecs();
    uint * lhs_vec = cdef->reac_lhs_bgn(cdef->reacG2L(pReacdef->gidx()));
    uint * cnt_vec = pTet->pools();

    // Compute combinatorial part.
    double h_mu = 1.0;
    for (uint pool = 0; pool < nspecs; ++pool)
    {
        uint lhs = lhs_vec[pool];
        if (lhs == 0) continue;
        uint cnt = cnt_vec[pool];
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

std::vector<uint> const & stex::Reac::apply(steps::rng::RNG * rng)
{
    uint * local = pTet->pools();
    ssolver::Compdef * cdef = pTet->compdef();
    uint l_ridx = cdef->reacG2L(pReacdef->gidx());
    int * upd_vec = cdef->reac_upd_bgn(l_ridx);
    uint nspecs = cdef->countSpecs();
    for (uint i = 0; i < nspecs; ++i)
    {
        if (pTet->clamped(i) == true) continue;
        int j = upd_vec[i];
        if (j == 0) continue;
        int nc = static_cast<int>(local[i]) + j;
        pTet->setCount(i, static_cast<uint>(nc));
    }
    rExtent++;
    return pUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

// END
