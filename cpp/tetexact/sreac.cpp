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
#include "sreac.hpp"
#include "tri.hpp"
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

stex::SReac::SReac(ssolver::SReacdef * srdef, stex::Tri * tri)
: KProc()
, pSReacdef(srdef)
, pTri(tri)
, pUpdVec()
, pCcst(0.0)
, pKcst(0.0)
{
	assert (pSReacdef != 0);
	assert (pTri != 0);
	double vol;
	if (pSReacdef->inside() == true)
	{
		assert(pTri->iTet() != 0);
		vol = pTri->iTet()->vol();
	}
	else
	{
		assert (pTri->oTet() != 0);
		vol = pTri->oTet()->vol();
	}

	uint lsridx = pTri->patchdef()->sreacG2L(pSReacdef->gidx());
	double kcst = pTri->patchdef()->kcst(lsridx);
	pKcst = kcst;
	pCcst = comp_ccst(kcst, vol, pSReacdef->order());
	assert (pCcst >= 0);
}

////////////////////////////////////////////////////////////////////////////////

stex::SReac::~SReac(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void stex::SReac::reset(void)
{
    resetExtent();
    resetCcst();
	setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void stex::SReac::resetCcst(void)
{
	double vol;
	if (pSReacdef->inside() == true)
	{
		assert(pTri->iTet() != 0);
		vol = pTri->iTet()->vol();
	}
	else
	{
		assert (pTri->oTet() != 0);
		vol = pTri->oTet()->vol();
	}
	uint lsridx = pTri->patchdef()->sreacG2L(pSReacdef->gidx());
	double kcst = pTri->patchdef()->kcst(lsridx);
	pKcst = kcst;
	pCcst = comp_ccst(kcst, vol, pSReacdef->order());
	assert (pCcst >= 0);
}

////////////////////////////////////////////////////////////////////////////////

void stex::SReac::setKcst(double k)
{
	assert (k >= 0.0);
	pKcst = k;

	double vol;
	if (pSReacdef->inside() == true)
	{
		assert(pTri->iTet() != 0);
		vol = pTri->iTet()->vol();
	}
	else
	{
		assert (pTri->oTet() != 0);
		vol = pTri->oTet()->vol();
	}

	pCcst = comp_ccst(k, vol, pSReacdef->order());
	assert(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void stex::SReac::setupDeps(void)
{
    // For all non-zero entries gidx in SReacDef's UPD_S:
    //   Perform depSpecTri(gidx,tri()) for:
    //     All kproc's of tri()
    //     NOTE: we currently don't test kproc's of inner and external tet,
    //     because that's not strictly necessary.
    //
    // If inner tetrahedron exists:
    //   For all non-zero entries gidx in SReacDef's UPD_I:
    //     Perform depSpecTet(gidx,itet) for:
    //       All kproc's of itet
    //       All kproc's of triangles next to itet
    //
    // If outer tetrahedron exists:
    //   Similar to inner tet.
    //
    // All dependencies are first collected into a std::set, to sort them
    // and to eliminate duplicates. At the end of the routine, they are
    // copied into the vector that will be returned during execution.

    Tet * itet = pTri->iTet();
    Tet * otet = pTri->oTet();

    ssolver::gidxTVecCI s_beg = pSReacdef->beginUpdColl_S();
    ssolver::gidxTVecCI s_end = pSReacdef->endUpdColl_S();
    ssolver::gidxTVecCI i_beg = pSReacdef->beginUpdColl_I();
    ssolver::gidxTVecCI i_end = pSReacdef->endUpdColl_I();
    ssolver::gidxTVecCI o_beg = pSReacdef->beginUpdColl_O();
    ssolver::gidxTVecCI o_end = pSReacdef->endUpdColl_O();

    SchedIDXSet updset;

    KProcPVecCI kprocend = pTri->kprocEnd();
    for (KProcPVecCI k = pTri->kprocBegin(); k != kprocend; ++k)
    {
        for (ssolver::gidxTVecCI spec = s_beg; spec != s_end; ++spec)
        {
            if ((*k)->depSpecTri(*spec, pTri) == true)
                updset.insert((*k)->schedIDX());
        }
    }

    if (itet != 0)
    {
        kprocend = itet->kprocEnd();
        for (KProcPVecCI k = itet->kprocBegin(); k != kprocend; ++k)
        {
            for (ssolver::gidxTVecCI spec = i_beg; spec != i_end; ++spec)
            {
                if ((*k)->depSpecTet(*spec, itet) == true)
                    updset.insert((*k)->schedIDX());
            }
        }

        for (uint i = 0; i < 4; ++i)
        {
            stex::Tri * tri = itet->nextTri(i);
            if (tri == 0) continue;
            kprocend = tri->kprocEnd();
            for (KProcPVecCI k = tri->kprocBegin(); k != kprocend; ++k)
            {
                for (ssolver::gidxTVecCI spec = i_beg; spec != i_end; ++spec)
                {
                    if ((*k)->depSpecTet(*spec, itet) == true)
                        updset.insert((*k)->schedIDX());
                }
            }
        }
    }

    if (otet != 0)
    {
        kprocend = otet->kprocEnd();
        for (KProcPVecCI k = otet->kprocBegin(); k != kprocend; ++k)
        {
            for (ssolver::gidxTVecCI spec = o_beg; spec != o_end; ++spec)
            {
                if ((*k)->depSpecTet(*spec, otet) == true)
                    updset.insert((*k)->schedIDX());
            }
        }

        for (uint i = 0; i < 4; ++i)
        {
            stex::Tri * tri = otet->nextTri(i);
            if (tri == 0) continue;
            kprocend = tri->kprocEnd();
            for (KProcPVecCI k = tri->kprocBegin(); k != kprocend; ++k)
            {
                for (ssolver::gidxTVecCI spec = o_beg; spec != o_end; ++spec)
                {
                    if ((*k)->depSpecTet(*spec, otet) == true)
                        updset.insert((*k)->schedIDX());
                }
            }
        }
    }

    stex::schedIDXSet_To_Vec(updset, pUpdVec);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::SReac::depSpecTet(uint gidx, stex::Tet * tet)
{
    // We need to check whether the tet is inside or outside.
    //   -> If inside: check dependency using SReacDef's I_DEP
    //   -> If outside: check dependency using SReacDef's O_DEP
    //   -> If neither, return.
    //
    // NOTE: DEP_NONE is defined in steps/sim/shared/types.hpp

    if (tet == pTri->iTet())
    {
        return (pSReacdef->dep_I(gidx) != ssolver::DEP_NONE);
    }
    else if (tet == pTri->oTet())
    {
        return (pSReacdef->dep_O(gidx) != ssolver::DEP_NONE);
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool stex::SReac::depSpecTri(uint gidx, stex::Tri * triangle)
{
    if (triangle != pTri) return false;
    return (pSReacdef->dep_S(gidx) != ssolver::DEP_NONE);
}

////////////////////////////////////////////////////////////////////////////////

double stex::SReac::rate(void) const
{
	   if (inactive()) return 0.0;

	    // First we compute the combinatorial part.
	    //   1/ for the surface part of the stoichiometry
	    //   2/ for the inner or outer volume part of the stoichiometry,
	    //      depending on whether the sreac is inner() or outer()
	    // Then we multiply with mesoscopic constant.

	    ssolver::Patchdef * pdef = pTri->patchdef();
	    uint lidx = pdef->sreacG2L(pSReacdef->gidx());

	    double h_mu = 1.0;

	    uint * lhs_s_vec = pdef->sreac_lhs_S_bgn(lidx);
	    uint * cnt_s_vec = pTri->pools();
	    uint nspecs_s = pdef->countSpecs();
	    for (uint s = 0; s < nspecs_s; ++s)
	    {
	        uint lhs = lhs_s_vec[s];
	        if (lhs == 0) continue;
	        uint cnt = cnt_s_vec[s];
	        if (lhs > cnt)
	        {
	            return 0.0;
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

	    if (pSReacdef->inside())
	    {
	        uint * lhs_i_vec = pdef->sreac_lhs_I_bgn(lidx);
	        uint * cnt_i_vec = pTri->iTet()->pools();
	        uint nspecs_i = pdef->countSpecs_I();
	        for (uint s = 0; s < nspecs_i; ++s)
	        {
	            uint lhs = lhs_i_vec[s];
	            if (lhs == 0) continue;
	            uint cnt = cnt_i_vec[s];
	            if (lhs > cnt)
	            {
	                return 0.0;
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
	    }
	    else if (pSReacdef->outside())
	    {
	        uint * lhs_o_vec = pdef->sreac_lhs_O_bgn(lidx);
	        uint * cnt_o_vec = pTri->oTet()->pools();
	        uint nspecs_o = pdef->countSpecs_O();
	        for (uint s = 0; s < nspecs_o; ++s)
	        {
	            uint lhs = lhs_o_vec[s];
	            if (lhs == 0) continue;
	            uint cnt = cnt_o_vec[s];
	            if (lhs > cnt)
	            {
	                return 0.0;
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
	    }

	    return h_mu * pCcst;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> const & stex::SReac::apply(steps::rng::RNG * rng)
{
    ssolver::Patchdef * pdef = pTri->patchdef();
    uint lidx = pdef->sreacG2L(pSReacdef->gidx());

    // Update triangle pools.
    int * upd_s_vec = pdef->sreac_upd_S_bgn(lidx);
    uint * cnt_s_vec = pTri->pools();
    uint nspecs_s = pdef->countSpecs();
    for (uint s = 0; s < nspecs_s; ++s)
    {
        if (pTri->clamped(s) == true) continue;
        int upd = upd_s_vec[s];
        if (upd == 0) continue;
        int nc = static_cast<int>(cnt_s_vec[s]) + upd;
        assert(nc >= 0);
        pTri->setCount(s, static_cast<uint>(nc));
    }

    // Update inner tet pools.
    stex::Tet * itet = pTri->iTet();
    if (itet != 0)
    {
        int * upd_i_vec = pdef->sreac_upd_I_bgn(lidx);
        uint * cnt_i_vec = itet->pools();
        uint nspecs_i = pdef->countSpecs_I();
        for (uint s = 0; s < nspecs_i; ++s)
        {
            if (itet->clamped(s) == true) continue;
            int upd = upd_i_vec[s];
            if (upd == 0) continue;
            int nc = static_cast<int>(cnt_i_vec[s]) + upd;
            assert(nc >= 0);
            itet->setCount(s, static_cast<uint>(nc));
        }
    }

    // Update outer tet pools.
    stex::Tet * otet = pTri->oTet();
    if (otet != 0)
    {
        int * upd_o_vec = pdef->sreac_upd_O_bgn(lidx);
        uint * cnt_o_vec = otet->pools();
        uint nspecs_o = pdef->countSpecs_O();
        for (uint s = 0; s < nspecs_o; ++s)
        {
            if (otet->clamped(s) == true) continue;
            int upd = upd_o_vec[s];
            if (upd == 0) continue;
            int nc = static_cast<int>(cnt_o_vec[s]) + upd;
            assert(nc >= 0);
            otet->setCount(s, static_cast<uint>(nc));
        }
    }

    rExtent++;

    // Update outer tet pools.
    return pUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

// END
