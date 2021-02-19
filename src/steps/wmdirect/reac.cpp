/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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
// #include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/math/constants.hpp"
#include "steps/wmdirect/comp.hpp"
#include "steps/wmdirect/kproc.hpp"
#include "steps/wmdirect/reac.hpp"
#include "steps/wmdirect/wmdirect.hpp"
// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace swmd = steps::wmdirect;
namespace ssolver = steps::solver;
namespace smath = steps::math;

////////////////////////////////////////////////////////////////////////////////

static inline double comp_ccst(double kcst, double vol, uint order)
{
    double vscale = 1.0e3 * vol * smath::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;

    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;

    return kcst * pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

swmd::Reac::Reac(ssolver::Reacdef * rdef, swmd::Comp * comp)
: 
 pReacdef(rdef)
, pComp(comp)
{
    AssertLog(pReacdef != nullptr);
    AssertLog(pComp != nullptr);
    uint lridx = pComp->def()->reacG2L(pReacdef->gidx());
    double kcst = pComp->def()->kcst(lridx);
    pCcst = comp_ccst(kcst, pComp->def()->vol(), pReacdef->order());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

swmd::Reac::~Reac() = default;

////////////////////////////////////////////////////////////////////////////////

void swmd::Reac::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(&pCcst), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Reac::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(&pCcst), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

bool swmd::Reac::active() const
{
    uint lridx = pComp->def()->reacG2L(defr()->gidx());
    return pComp->def()->active(lridx);
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Reac::reset()
{
    resetExtent();
    uint lridx = pComp->def()->reacG2L(defr()->gidx());
    pComp->def()->setActive(lridx, true);
    resetCcst();
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Reac::resetCcst()
{
    uint lridx = pComp->def()->reacG2L(pReacdef->gidx());
    double kcst = pComp->def()->kcst(lridx);
    pCcst = comp_ccst(kcst, pComp->def()->vol(), pReacdef->order());
    AssertLog(pCcst >= 0);

}

////////////////////////////////////////////////////////////////////////////////

void swmd::Reac::setupDeps()
{
    SchedIDXSet updset;

    // Search in local compartment.
    for (auto const& k : pComp->kprocs()) {
        for (auto const& s : defr()->updColl()) {
            if (k->depSpecComp(s, pComp))
                updset.insert(k->schedIDX());
        }
    }

    // Search in neighbouring patches.
    for (auto const& p : pComp->ipatches()) {
        for (auto const& k : p->kprocs()) {
            for (auto const& s : defr()->updColl()) {
                if (k->depSpecComp(s, pComp))
                    updset.insert(k->schedIDX());
            }
        }
    }
    for (auto const& p : pComp->opatches()) {
        for (auto const& k : p->kprocs()) {
            for (auto const& s : defr()->updColl()) {
                if (k->depSpecComp(s, pComp))
                    updset.insert(k->schedIDX());
            }
        }
    }

    swmd::schedIDXSet_To_Vec(updset, pUpdVec);
}

////////////////////////////////////////////////////////////////////////////////

bool swmd::Reac::depSpecComp(uint gidx, swmd::Comp * comp)
{
    if (pComp != comp) { return false;
}
    return defr()->dep(gidx) != 0;
}

////////////////////////////////////////////////////////////////////////////////

bool swmd::Reac::depSpecPatch(uint /*gidx*/, swmd::Patch * /*patch*/)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

double swmd::Reac::rate() const
{
    if (inactive()) { return 0.0;
}

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
            if (lhs == 0) { continue;
}
            auto cnt = static_cast<uint>(cnt_vec[pool]);
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
                    AssertLog(false);
                }
            }
        }

        // Multiply with scaled reaction constant.
        return h_mu * pCcst;

}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> const & swmd::Reac::apply()
{
    ssolver::Compdef * cdef = pComp->def();
    double * local = cdef->pools();
    uint l_ridx = cdef->reacG2L(defr()->gidx());
    int * upd_vec = cdef->reac_upd_bgn(l_ridx);
    uint nspecs = cdef->countSpecs();
    for (uint i=0; i < nspecs; ++i)
    {
        if (cdef->clamped(i)) continue;
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

