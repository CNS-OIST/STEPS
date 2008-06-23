////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
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
// $Id:func_tetmesh.cpp 64 2007-08-20 06:25:41Z stefan $
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/constants.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/patchdef.hpp>
#include <steps/sim/swiginf/func_tetmesh.hpp>
#include <steps/tetexact/solver_core/diff.hpp>
#include <steps/tetexact/solver_core/reac.hpp>
#include <steps/tetexact/solver_core/sched.hpp>
#include <steps/tetexact/solver_core/sreac.hpp>
#include <steps/tetexact/solver_core/state.hpp>
#include <steps/tetexact/solver_core/tet.hpp>
#include <steps/tetexact/solver_core/tri.hpp>

NAMESPACE_ALIAS(steps::math, smath);
NAMESPACE_ALIAS(steps::sim, ssim);

////////////////////////////////////////////////////////////////////////////////

void siBeginTetmeshDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siEndTetmeshDef(State *s)
{
    s->setupTetmesh();
}

////////////////////////////////////////////////////////////////////////////////

void siBeginTetDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siEndTetDef(State * s)
{
    assert(s != 0);
    CompPVecCI c_end = s->endComp();
    for (CompPVecCI c = s->bgnComp(); c != c_end; ++c)
    {
        (*c)->computeVol();
    }
}

////////////////////////////////////////////////////////////////////////////////

uint siNewTet(State * s, uint cidx, double vol, 
    double a1, double a2, double a3, double a4,
    double d1, double d2, double d3, double d4)
{
    assert(s != 0);
    Comp * comp = s->comp(cidx);
    assert(comp != 0);
    return s->addTet(comp, vol, a1, a2, a3, a4, d1, d2, d3, d4);
}

////////////////////////////////////////////////////////////////////////////////

void siBeginTriDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siEndTriDef(State * s)
{
    assert(s != 0);
    PatchPVecCI p_end = s->endPatch();
    for (PatchPVecCI p = s->bgnPatch(); p != p_end; ++p)
    {
        (*p)->computeArea();
    }
}

////////////////////////////////////////////////////////////////////////////////

uint siNewTri(State * s, uint pidx, double area)
{
    assert(s != 0);
    Patch * patch = s->patch(pidx);
    assert(patch != 0);
    assert(area >= 0.0);
    return s->addTri(patch, area);
}

////////////////////////////////////////////////////////////////////////////////

void siBeginConnectDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siEndConnectDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siConnectTetTet(State * s, uint side, uint tidx1, uint tidx2)
{
    Tet * t1 = s->tet(tidx1);
    Tet * t2 = s->tet(tidx2);
    t1->setNextTet(side, t2);
}

////////////////////////////////////////////////////////////////////////////////

void siConnectTetTri(State * s, uint side, uint tetidx, uint triidx)
{
    TetP tet = s->tet(tetidx);
    assert(tet != 0);
    TriP tri = s->tri(triidx);
    assert(tri != 0);
    tet->setNextTri(side, tri);
}

////////////////////////////////////////////////////////////////////////////////

void siConnectTriTetInner(State * s, uint triidx, uint tetidx)
{
    TetP tet = s->tet(tetidx);
    assert(tet != 0);
    TriP tri = s->tri(triidx);
    assert(tri != 0);
    tri->setInnerTet(tet);
}

////////////////////////////////////////////////////////////////////////////////

void siConnectTriTetOuter(State * s, uint triidx, uint tetidx)
{
    TetP tet = s->tet(tetidx);
    assert(tet != 0);
    TriP tri = s->tri(triidx);
    assert(tri != 0);
    tri->setOuterTet(tet);
}

////////////////////////////////////////////////////////////////////////////////

double siGetTetVol(State * s, uint tidx)
{
    TetP tet = s->tet(tidx);
    assert(tet != 0);
    return tet->vol();
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetVol(State * s, uint tidx, double vol)
{
    // Will probably never be implemented. On the other hand, might be
    // a (very) cheap way of implementing dynamic morphology (as in 
    // structural plasticity).
}

////////////////////////////////////////////////////////////////////////////////

uint siGetTetCount(State * s, uint tidx, uint sidx)
{
    assert(s != 0);
    Tet * tet = s->tet(tidx);
    assert(tet != 0);
    
    uint l_sidx = tet->compdef()->specG2L(sidx);
    if (l_sidx == ssim::LIDX_UNDEFINED) return 0;
    return tet->pools()[l_sidx];
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetCount(State * s, uint tidx, uint sidx, uint n)
{
    assert(s != 0);
    Tet * tet = s->tet(tidx);
    assert(tet != 0);
    
    // Apply the change.
    uint l_sidx = tet->compdef()->specG2L(sidx);
    if (l_sidx == ssim::LIDX_UNDEFINED) return;
    tet->pools()[l_sidx] = n;
    
    // Send the list of kprocs that need to be updated to the schedule.
    s->sched()->updateSpec(tet, l_sidx);
}

////////////////////////////////////////////////////////////////////////////////

double siGetTetAmount(State * s, uint tidx, uint sidx)
{
    double count = static_cast<double>(siGetTetCount(s, tidx, sidx));
    return count / smath::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetAmount(State * s, uint tidx, uint sidx, double m)
{
    assert(s != 0);
    assert(m >= 0.0);
    
    double m2 = m * smath::AVOGADRO;
    double m_int = std::floor(m2);
    double m_frc = m2 - m_int;
    uint c = static_cast<uint>(m_int);
    if (m_frc > 0.0)
    {
        double rand01 = s->rng()->getUnfIE();
        if (rand01 < m_frc) c++;
    }
    
    siSetTetCount(s, tidx, sidx, c);
}

////////////////////////////////////////////////////////////////////////////////

double siGetTetConc(State * s, uint tidx, uint sidx)
{
    assert(s != 0);
    Tet * tet = s->tet(tidx);
    assert(tet != 0);
    
    double count = static_cast<double>(siGetTetCount(s, tidx, sidx));
    double vol = tet->vol();
    return count / (1.0e3 * vol * smath::AVOGADRO); 
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetConc(State * s, uint tidx, uint sidx, double c)
{
    assert(s != 0);
    assert(c >= 0.0);

    Tet * tet = s->tet(tidx);
    assert(tet != 0);    
    
    double c2 = c * (1.0e3 * tet->vol() * smath::AVOGADRO);
    double c_int = std::floor(c2);
    double c_frc = c2 - c_int;
    uint count = static_cast<uint>(c_int);
    if (c_frc > 0.0)
    {
        double rand01 = s->rng()->getUnfIE();
        if (rand01 < c_frc) count++;
    }
    
    siSetTetCount(s, tidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool siGetTetClamped(State * s, uint tidx, uint sidx)
{
    assert(s != 0);
    Tet * tet = s->tet(tidx);
    assert(tet != 0);
    
    uint lsidx = tet->compdef()->specG2L(sidx);
    if (lsidx == ssim::LIDX_UNDEFINED) return false;
    
    return tet->clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetClamped(State * s, uint tidx, uint sidx, bool buf)
{
    assert(s != 0);
    Tet * tet = s->tet(tidx);
    assert(tet != 0);
    
    uint lsidx = tet->compdef()->specG2L(sidx);
    if (lsidx == ssim::LIDX_UNDEFINED) return;
    
    tet->setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double siGetTetReacK(State * s, uint tidx, uint ridx)
{
    // Currently not implemented.
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetReacK(State * s, uint tidx, uint ridx, double kf)
{
    // Currently not implemented.
}

////////////////////////////////////////////////////////////////////////////////

double siGetTetReacA(State * s, uint tidx, uint ridx)
{
    assert(s != 0);
    Tet * tet = s->tet(tidx);
    assert(tet != 0);
    
    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssim::LIDX_UNDEFINED) return 0.0;
    
    return tet->reac(lridx)->rate();
}

////////////////////////////////////////////////////////////////////////////////

bool siGetTetReacActive(State * s, uint tidx, uint ridx)
{
    assert(s != 0);
    Tet * tet = s->tet(tidx);
    assert(tet != 0);
    
    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssim::LIDX_UNDEFINED) return false;
    
    if (tet->reac(lridx)->inactive() == true) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetReacActive(State * s, uint tidx, uint ridx, bool act)
{
    assert(s != 0);
    Tet * tet = s->tet(tidx);
    assert(tet != 0);
    
    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssim::LIDX_UNDEFINED) return;
    
    tet->reac(lridx)->setActive(act);
    
    SchedIDXVec updvec;
    updvec.push_back(tet->reac(lridx)->schedIDX());
    s->sched()->update(updvec);
}

////////////////////////////////////////////////////////////////////////////////

double siGetTetDiffD(State * s, uint tidx, uint didx)
{
    // Currently not implemented.
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetDiffD(State * s, uint tidx, uint didx)
{
    // Currently not implemented.
}

////////////////////////////////////////////////////////////////////////////////

double siGetTetDiffA(State * s, uint tidx, uint didx)
{
    assert(s != 0);
    Tet * tet = s->tet(tidx);
    assert(tet != 0);
    
    uint ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx == ssim::LIDX_UNDEFINED) return 0.0;
    
    return tet->diff(ldidx)->rate();
}

////////////////////////////////////////////////////////////////////////////////

bool siGetTetDiffActive(State * s, uint tidx, uint didx)
{
    assert(s != 0);
    Tet * tet = s->tet(tidx);
    assert(tet != 0);
    
    uint ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx == ssim::LIDX_UNDEFINED) return false;
    
    if (tet->diff(ldidx)->inactive() == true) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetDiffActive(State * s, uint tidx, uint didx, bool act)
{
    assert(s != 0);
    Tet * tet = s->tet(tidx);
    assert(tet != 0);
    
    uint ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx == ssim::LIDX_UNDEFINED) return;
    
    tet->diff(ldidx)->setActive(act);
    
    SchedIDXVec updvec;
    updvec.push_back(tet->diff(ldidx)->schedIDX());
    s->sched()->update(updvec);
}

////////////////////////////////////////////////////////////////////////////////

double siGetTriArea(State * s, uint tidx)
{
    assert(s != 0);
    Tri * tri = s->tri(tidx);
    assert(tri != 0);
    return tri->area();
}

////////////////////////////////////////////////////////////////////////////////

void siSetTriArea(State * s, uint tidx, double area)
{
    // Will probably never be implemented. On the other hand, might be
    // a (very) cheap way of implementing dynamic morphology (as in 
    // structural plasticity).
}

////////////////////////////////////////////////////////////////////////////////

uint siGetTriCount(State * s, uint tidx, uint sidx)
{
    assert(s != 0);
    Tri * tri = s->tri(tidx);
    assert(tri != 0);
    
    uint l_sidx = tri->patchdef()->specG2L(sidx);
    if (l_sidx == ssim::LIDX_UNDEFINED) return 0;
    return tri->pools()[l_sidx];
}

////////////////////////////////////////////////////////////////////////////////

void siSetTriCount(State * s, uint tidx, uint sidx, uint n)
{
    assert(s != 0);
    Tri * tri = s->tri(tidx);
    assert(tri != 0);
    
    // Apply the change.
    uint lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx == ssim::LIDX_UNDEFINED) return;
    tri->pools()[lsidx] = n;
    
    // Send the list of kprocs that need to be updated to the schedule.
    s->sched()->updateSpec(tri, lsidx);
}

////////////////////////////////////////////////////////////////////////////////

bool siGetTriClamped(State * s, uint tidx, uint sidx)
{
    assert(s != 0);
    Tri * tri = s->tri(tidx);
    assert(tri != 0);
    
    uint lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx == ssim::LIDX_UNDEFINED) return false;
    
    return tri->clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////////////

void siSetTriClamped(State * s, uint tidx, uint sidx, bool buf)
{
    assert(s != 0);
    Tri * tri = s->tri(tidx);
    assert(tri != 0);
    
    uint lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx == ssim::LIDX_UNDEFINED) return;
    
    tri->setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double siGetTriSReacK(State * s, uint tidx, uint ridx)
{
    // Currently not implemented.
}

////////////////////////////////////////////////////////////////////////////////

void siSetTriSReacK(State * s, uint tidx, uint ridx, double kf)
{
    // Currently not implemented.
}

////////////////////////////////////////////////////////////////////////////////

bool siGetTriSReacActive(State * s, uint tidx, uint ridx)
{
    assert(s != 0);
    Tri * tri = s->tri(tidx);
    assert(tri != 0);
    
    uint lridx = tri->patchdef()->sreacG2L(ridx);
    if (lridx == ssim::LIDX_UNDEFINED) return false;
    
    if (tri->sreac(lridx)->inactive() == true) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void siSetTriSReacActive(State * s, uint tidx, uint ridx, bool act)
{
    assert(s != 0);
    Tri * tri = s->tri(tidx);
    assert(tri != 0);
    
    uint lridx = tri->patchdef()->sreacG2L(ridx);
    if (lridx == ssim::LIDX_UNDEFINED) return;
    
    tri->sreac(lridx)->setActive(act);
    
    SchedIDXVec updvec;
    updvec.push_back(tri->sreac(lridx)->schedIDX());
    s->sched()->update(updvec);
}

////////////////////////////////////////////////////////////////////////////////

// END
