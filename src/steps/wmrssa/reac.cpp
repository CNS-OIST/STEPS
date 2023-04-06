/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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
#include "reac.hpp"
#include "comp.hpp"
#include "wmrssa.hpp"
#include "math/constants.hpp"

// logging
#include <easylogging++.h>
#include "util/error.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace swmrssa = steps::wmrssa;
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

swmrssa::Reac::Reac(ssolver::Reacdef * rdef, swmrssa::Comp * comp)
:
 pReacdef(rdef)
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

swmrssa::Reac::~Reac() = default;

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Reac::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(&pCcst), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Reac::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(&pCcst), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

bool swmrssa::Reac::active() const
{
    uint lridx = pComp->def()->reacG2L(defr()->gidx());
    return pComp->def()->active(lridx);
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Reac::reset()
{
    resetExtent();
    uint lridx = pComp->def()->reacG2L(defr()->gidx());
    pComp->def()->setActive(lridx, true);
    resetCcst();
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Reac::resetCcst()
{
    uint lridx = pComp->def()->reacG2L(pReacdef->gidx());
    double kcst = pComp->def()->kcst(lridx);
    pCcst = comp_ccst(kcst, pComp->def()->vol(), pReacdef->order());
    assert (pCcst >= 0);

}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Reac::setupDeps()
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

    swmrssa::schedIDXSet_To_Vec(updset, pUpdVec);
}

////////////////////////////////////////////////////////////////////////////////

bool swmrssa::Reac::depSpecComp(uint gidx, swmrssa::Comp * comp)
{
    if (pComp != comp) { return false;
}
    return defr()->dep(gidx) != 0;
}

////////////////////////////////////////////////////////////////////////////////

bool swmrssa::Reac::depSpecPatch(uint /*gidx*/, swmrssa::Patch * /*patch*/)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

double swmrssa::Reac::rate(steps::wmrssa::PropensityRSSA prssa)
{
    if (inactive()) { return 0.0;
}

    if (prssa == steps::wmrssa::BOUNDS) {
        pPropensityLB = rate(steps::wmrssa::LOWERBOUND);
}

    // Prefetch some variables.
    ssolver::Compdef * cdef = pComp->def();
    uint nspecs = cdef->countSpecs();
    uint * lhs_vec = cdef->reac_lhs_bgn(cdef->reacG2L(defr()->gidx()));
    double * cnt_vec = pComp->pools(prssa);

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
                    AssertLog(false);
                }
            }
        }

        // Multiply with scaled reaction constant.
        return h_mu * pCcst;

}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> const & swmrssa::Reac::apply()
{
    SchedIDXSet updset;
    //bool returnUpdVec = false;
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
        if (pComp->isOutOfBound(i, nc))
        {
            std::vector<steps::wmrssa::KProc *> dependentReacs = pComp->getSpecUpdKProcs(i);
            for (auto &dependentReac : dependentReacs)
                updset.insert(dependentReac->schedIDX());
            //returnUpdVec = true;
        }
    }
    rExtent++;
    swmrssa::schedIDXSet_To_Vec(updset, pUpdVec);
    return pUpdVec; //returnUpdVec ? pUpdVec : emptyVec;
}

////////////////////////////////////////////////////////////////////////////////

// END

