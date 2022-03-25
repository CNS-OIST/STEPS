/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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
#include "reac.hpp"
#include "tet.hpp"
#include "tetexact.hpp"
#include "wmvol.hpp"
#include "math/constants.hpp"

// logging
#include <easylogging++.h>
#include "util/error.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;
namespace ssolver = steps::solver;
namespace smath = steps::math;

////////////////////////////////////////////////////////////////////////////////

static inline double comp_ccst(double kcst, double vol, uint order, double /*compvol*/)
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
    AssertLog(pReacdef != nullptr);
    AssertLog(pTet != nullptr);

    uint lridx = pTet->compdef()->reacG2L(pReacdef->gidx());
    double kcst = pTet->compdef()->kcst(lridx);
    pKcst = kcst;
    pCcst = comp_ccst(kcst, pTet->vol(), pReacdef->order(), pTet->compdef()->vol());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

stex::Reac::~Reac() = default;

////////////////////////////////////////////////////////////////////////////////

void stex::Reac::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(&rExtent), sizeof(unsigned long long));
    cp_file.write(reinterpret_cast<char*>(&pFlags), sizeof(uint));

    cp_file.write(reinterpret_cast<char*>(&pCcst), sizeof(double));
    cp_file.write(reinterpret_cast<char*>(&pKcst), sizeof(double));

    cp_file.write(reinterpret_cast<char*>(&crData.recorded), sizeof(bool));
    cp_file.write(reinterpret_cast<char*>(&crData.pow), sizeof(int));
    cp_file.write(reinterpret_cast<char*>(&crData.pos), sizeof(unsigned));
    cp_file.write(reinterpret_cast<char*>(&crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Reac::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(&rExtent), sizeof(unsigned long long));
    cp_file.read(reinterpret_cast<char*>(&pFlags), sizeof(uint));

    cp_file.read(reinterpret_cast<char*>(&pCcst), sizeof(double));
    cp_file.read(reinterpret_cast<char*>(&pKcst), sizeof(double));

    cp_file.read(reinterpret_cast<char*>(&crData.recorded), sizeof(bool));
    cp_file.read(reinterpret_cast<char*>(&crData.pow), sizeof(int));
    cp_file.read(reinterpret_cast<char*>(&crData.pos), sizeof(unsigned));
    cp_file.read(reinterpret_cast<char*>(&crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Reac::reset()
{

    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    resetExtent();
    _resetCcst();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Reac::_resetCcst()
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

    // Search in local tetrahedron.
    for (auto const& k : pTet->kprocs()) {
        for (auto const& s : pReacdef->UPD_Coll()) {
            if (k->depSpecTet(s, pTet)) {
                //updset.insert((*k)->getSSARef());
                updset.insert(k);
            }
        }
    }

    for (auto const& tri : pTet->nexttris()) {
        if (tri == nullptr) continue;

        for (auto const& k : tri->kprocs()) {
            for (auto const& s : pReacdef->UPD_Coll()) {
                if (k->depSpecTet(s, pTet) == true) {
                    //updset.insert((*k)->getSSARef());
                    updset.insert(k);
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

bool stex::Reac::depSpecTri(uint /*gidx*/, stex::Tri * /*tri*/)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Reac::rate(steps::tetexact::Tetexact * /*solver*/)
{
    if (inactive()) return 0.0;

    // Prefetch some variables.
    ssolver::Compdef * cdef = pTet->compdef();
    uint nspecs = cdef->countSpecs();
    uint * lhs_vec = cdef->reac_lhs_bgn(cdef->reacG2L(pReacdef->gidx()));
    auto const& cnt_vec = pTet->pools();

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
            STEPS_FALLTHROUGH;
            case 3:
            {
                h_mu *= static_cast<double>(cnt - 2);
            }
            STEPS_FALLTHROUGH;
            case 2:
            {
                h_mu *= static_cast<double>(cnt - 1);
            }
            STEPS_FALLTHROUGH;
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

std::vector<stex::KProc*> const & stex::Reac::apply(const rng::RNGptr &/*rng*/, double /*dt*/, double /*simtime*/)
{
    auto const& local = pTet->pools();
    ssolver::Compdef * cdef = pTet->compdef();
    uint l_ridx = cdef->reacG2L(pReacdef->gidx());
    int * upd_vec = cdef->reac_upd_bgn(l_ridx);
    uint nspecs = cdef->countSpecs();
    for (uint i = 0; i < nspecs; ++i)
    {
        if (pTet->clamped(i)) continue;
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
