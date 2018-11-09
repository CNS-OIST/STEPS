/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */



// Standard library & STL headers.
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/math/constants.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
#include "steps/mpi/tetopsplit/reac.hpp"
#include "steps/mpi/tetopsplit/tet.hpp"
#include "steps/mpi/tetopsplit/tetopsplit.hpp"
#include "steps/mpi/tetopsplit/wmvol.hpp"
#include "steps/solver/reacdef.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;
namespace ssolver = steps::solver;
namespace smath = steps::math;

////////////////////////////////////////////////////////////////////////////////

static inline double comp_ccst(double kcst, double vol, uint order, double compvol)
{
    double vscale = 1.0e3 * vol * smath::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;

    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;

    double ccst = kcst * pow(vscale, static_cast<double>(-o1));

    return ccst;
}

////////////////////////////////////////////////////////////////////////////////

smtos::Reac::Reac(ssolver::Reacdef * rdef, smtos::WmVol * tet)
: 
 pReacdef(rdef)
, pTet(tet)
, localUpdVec()
, remoteUpdVec()
, pCcst(0.0)
, pKcst(0.0)
{
    AssertLog(pReacdef != 0);
    AssertLog(pTet != 0);
    type = KP_REAC;
    uint lridx = pTet->compdef()->reacG2L(pReacdef->gidx());
    double kcst = pTet->compdef()->kcst(lridx);
    pKcst = kcst;
    pCcst = comp_ccst(kcst, pTet->vol(), pReacdef->order(), pTet->compdef()->vol());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

smtos::Reac::~Reac()
= default;

////////////////////////////////////////////////////////////////////////////////

void smtos::Reac::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&rExtent, sizeof(uint));
    cp_file.write((char*)&pFlags, sizeof(uint));

    cp_file.write((char*)&pCcst, sizeof(double));
    cp_file.write((char*)&pKcst, sizeof(double));

    cp_file.write((char*)&(crData.recorded), sizeof(bool));
    cp_file.write((char*)&(crData.pow), sizeof(int));
    cp_file.write((char*)&(crData.pos), sizeof(unsigned));
    cp_file.write((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Reac::restore(std::fstream & cp_file)
{
    cp_file.read((char*)&rExtent, sizeof(uint));
    cp_file.read((char*)&pFlags, sizeof(uint));

    cp_file.read((char*)&pCcst, sizeof(double));
    cp_file.read((char*)&pKcst, sizeof(double));

    cp_file.read((char*)&(crData.recorded), sizeof(bool));
    cp_file.read((char*)&(crData.pow), sizeof(int));
    cp_file.read((char*)&(crData.pos), sizeof(unsigned));
    cp_file.read((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Reac::reset()
{

    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    resetExtent();
    resetCcst();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Reac::resetCcst()
{
    uint lridx = pTet->compdef()->reacG2L(pReacdef->gidx());
    double kcst = pTet->compdef()->kcst(lridx);
    // Also reset kcst
    pKcst = kcst;
    pCcst = comp_ccst(kcst, pTet->vol(), pReacdef->order(), pTet->compdef()->vol());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Reac::setKcst(double k)
{
    AssertLog(k >= 0.0);
    pKcst = k;
    pCcst = comp_ccst(k, pTet->vol(), pReacdef->order(), pTet->compdef()->vol());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Reac::setupDeps()
{
    AssertLog(pTet->getInHost());
    std::set<smtos::KProc*> updset;
    ssolver::gidxTVecCI sbgn = pReacdef->bgnUpdColl();
    ssolver::gidxTVecCI send = pReacdef->endUpdColl();

    // Search in local tetrahedron.
    uint nkprocs = pTet->countKProcs();
    uint startKProcIdx = pTet->getStartKProcIdx();
    
    for (uint k = 0; k < nkprocs; k++)
    {
        for (ssolver::gidxTVecCI s = sbgn; s != send; ++s)
        {
            if (pTet->KProcDepSpecTet(k, pTet, *s) == true) {
                updset.insert(pTet->getKProc(k));
            }

        }
    }

    std::vector<smtos::Tri *>::const_iterator tri_end = pTet->nexttriEnd();
    for (std::vector<smtos::Tri *>::const_iterator tri = pTet->nexttriBegin();
             tri != tri_end; ++tri)
    {
        if ((*tri) == 0) continue;
        if (pTet->getHost() != (*tri)->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << (*tri)->idx() << " and its compartment tetrahedron " << pTet->idx()  << " belong to different hosts.\n";
            NotImplErrLog(os.str());
        }
        nkprocs = (*tri)->countKProcs();
        startKProcIdx = (*tri)->getStartKProcIdx();
        for (uint sk = 0; sk < nkprocs; sk++)
        {
            for (ssolver::gidxTVecCI s = sbgn; s != send; ++s)
            {
                if ((*tri)->KProcDepSpecTet(sk, pTet, *s) == true) {
                    updset.insert((*tri)->getKProc(sk));
                }

            }
        }
    }

    localUpdVec.assign(updset.begin(), updset.end());
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::Reac::depSpecTet(uint gidx, smtos::WmVol * tet)
{
    if (pTet != tet) { return false;
}
    return pReacdef->dep(gidx) != 0;
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::Reac::depSpecTri(uint gidx, smtos::Tri * tri)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::Reac::rate(smtos::TetOpSplitP * solver)
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
        if (lhs == 0) { continue;
}
        uint cnt = cnt_vec[pool];
        if (lhs > cnt) {
            h_mu = 0.0;
            break;
        }
        AssertLog((lhs <= 4) && (lhs >0));
        h_mu *= (cnt - lhs + 1);
    }

    // Multiply with scaled reaction constant.
    return h_mu * pCcst;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Reac::apply(steps::rng::RNG * rng, double dt, double simtime, double period)
{
    uint * local = pTet->pools();
    ssolver::Compdef * cdef = pTet->compdef();
    uint l_ridx = cdef->reacG2L(pReacdef->gidx());
    int * upd_vec = cdef->reac_upd_bgn(l_ridx);
    uint nspecs = cdef->countSpecs();
    for (uint i = 0; i < nspecs; ++i)
    {
        if (pTet->clamped(i) == true) { continue;
}
        int j = upd_vec[i];
        int nc = static_cast<int>(local[i]) + j;
        AssertLog(nc >= 0);
        pTet->setCount(i, static_cast<uint>(nc), period);
    }

    rExtent++;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> const & smtos::Reac::getRemoteUpdVec(int direction)
{
    return remoteUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<smtos::KProc*> const & smtos::Reac::getLocalUpdVec(int direction)
{
    return localUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Reac::resetOccupancies()
{
    pTet->resetPoolOccupancy();
}
////////////////////////////////////////////////////////////////////////////////

// END
