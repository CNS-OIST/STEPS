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
#include <iostream>
#include <vector>
// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/math/constants.hpp"
#include "steps/solver/reacdef.hpp"
#include "steps/tetexact/kproc.hpp"
#include "steps/tetexact/reac.hpp"
#include "steps/tetexact/tet.hpp"
#include "steps/tetexact/tetexact.hpp"
#include "steps/tetexact/wmvol.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;
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

stex::Reac::Reac(ssolver::Reacdef * rdef, stex::WmVol * tet)
: 
 pReacdef(rdef)
, pTet(tet)
, pUpdVec()
, pCcst(0.0)
, pKcst(0.0)
{
    AssertLog(pReacdef != 0);
    AssertLog(pTet != 0);

    uint lridx = pTet->compdef()->reacG2L(pReacdef->gidx());
    double kcst = pTet->compdef()->kcst(lridx);
    pKcst = kcst;
    pCcst = comp_ccst(kcst, pTet->vol(), pReacdef->order(), pTet->compdef()->vol());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

stex::Reac::~Reac()
= default;

////////////////////////////////////////////////////////////////////////////////

void stex::Reac::checkpoint(std::fstream & cp_file)
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

void stex::Reac::restore(std::fstream & cp_file)
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

void stex::Reac::reset()
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

void stex::Reac::resetCcst()
{
    uint lridx = pTet->compdef()->reacG2L(pReacdef->gidx());
    double kcst = pTet->compdef()->kcst(lridx);
    // Also reset kcst
    pKcst = kcst;
    pCcst = comp_ccst(kcst, pTet->vol(), pReacdef->order(), pTet->compdef()->vol());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Reac::setKcst(double k)
{
    AssertLog(k >= 0.0);
    pKcst = k;
    pCcst = comp_ccst(k, pTet->vol(), pReacdef->order(), pTet->compdef()->vol());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Reac::setupDeps()
{
    std::set<stex::KProc*> updset;
    ssolver::gidxTVecCI sbgn = pReacdef->bgnUpdColl();
    ssolver::gidxTVecCI send = pReacdef->endUpdColl();

    // Search in local tetrahedron.
    KProcPVecCI kprocend = pTet->kprocEnd();
    for (KProcPVecCI k = pTet->kprocBegin(); k != kprocend; ++k)
    {
        for (ssolver::gidxTVecCI s = sbgn; s != send; ++s)
        {
            if ((*k)->depSpecTet(*s, pTet) == true) {
                //updset.insert((*k)->getSSARef());
                updset.insert(*k);
            }
        }
    }

    std::vector<stex::Tri *>::const_iterator tri_end = pTet->nexttriEnd();
    for (std::vector<stex::Tri *>::const_iterator tri = pTet->nexttriBegin();
             tri != tri_end; ++tri)
    {
        if ((*tri) == 0) continue;

        kprocend = (*tri)->kprocEnd();
        for (KProcPVecCI k = (*tri)->kprocBegin(); k != kprocend; ++k)
        {
            for (ssolver::gidxTVecCI s = sbgn; s != send; ++s)
            {
                if ((*k)->depSpecTet(*s, pTet) == true) {
                    //updset.insert((*k)->getSSARef());
                    updset.insert(*k);
                }
            }
        }
    }

    pUpdVec.assign(updset.begin(), updset.end());
    //pUpdObjVec.assign(updset_obj.begin(), updset_obj.end());
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Reac::depSpecTet(uint gidx, stex::WmVol * tet)
{
    if (pTet != tet) { return false;
}
    return pReacdef->dep(gidx) != 0;
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Reac::depSpecTri(uint gidx, stex::Tri * tri)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Reac::rate(steps::tetexact::Tetexact * solver)
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
                AssertLog(0);
                return 0.0;
            }
        }
    }

    // Multiply with scaled reaction constant.
    return h_mu * pCcst;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<stex::KProc*> const & stex::Reac::apply(steps::rng::RNG * rng, double dt, double simtime)
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
