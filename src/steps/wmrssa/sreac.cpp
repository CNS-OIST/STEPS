/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
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
#include "steps/math/constants.hpp"
#include "steps/wmrssa/comp.hpp"
#include "steps/wmrssa/patch.hpp"
#include "steps/wmrssa/sreac.hpp"
#include "steps/wmrssa/wmrssa.hpp"

// logging
#include "easylogging++.h"
#include "steps/error.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace swmrssa = steps::wmrssa;
namespace ssolver = steps::solver;
namespace smath = steps::math;

////////////////////////////////////////////////////////////////////////////////

static inline double comp_ccst_vol(double kcst, double vol, uint order)
{
    double vscale = 1.0e3 * vol * smath::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;
    // I.H 5/1/2011 Removed this strange special behaviour for zero-order
    // if (o1 < 0) o1 = 0;
    return kcst * pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

static inline double comp_ccst_area(double kcst, double area, uint order)
{
    double ascale = area * smath::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;
    // I.H 5/1/2011 Removed this strange special behaviour for zero-order
    // if (o1 < 0) o1 = 0;
    return kcst * pow(ascale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

swmrssa::SReac::SReac(ssolver::SReacdef * srdef, swmrssa::Patch * patch)
: 
 pSReacdef(srdef)
, pPatch(patch)
{
    assert (pSReacdef != nullptr);
    assert (pPatch != nullptr);

    uint lsridx = pPatch->def()->sreacG2L(defsr()->gidx());
    double kcst = pPatch->def()->kcst(lsridx);

    if (defsr()->surf_surf() == false)
    {
        double vol;
        if (defsr()->inside())
        {
            AssertLog(pPatch->iComp() != nullptr);
            vol = pPatch->iComp()->def()->vol();
        }
        else
        {
            assert (pPatch->oComp() != 0);
            vol = pPatch->oComp()->def()->vol();
        }

        pCcst = comp_ccst_vol(kcst, vol, defsr()->order());
    }
    else
    {
        double area;
        area = pPatch->def()->area();
        pCcst = comp_ccst_area(kcst, area, defsr()->order());
    }

    assert (pCcst >= 0);
}

////////////////////////////////////////////////////////////////////////////////

swmrssa::SReac::~SReac() = default;

////////////////////////////////////////////////////////////////////////////////

void swmrssa::SReac::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(&pCcst), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::SReac::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(&pCcst), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

bool swmrssa::SReac::active() const
{
    uint lsridx = pPatch->def()->sreacG2L(defsr()->gidx());
    return pPatch->def()->active(lsridx);
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::SReac::setupDeps()
{
    Comp * icomp = pPatch->iComp();
    Comp * ocomp = pPatch->oComp();

    SchedIDXSet updset;

    for (auto const& k : pPatch->kprocs()) {
        for (auto const& spec : defsr()->updColl_S()) {
            if (k->depSpecPatch(spec, pPatch))
                updset.insert(k->schedIDX());
        }
    }

    if (icomp != nullptr)
    {
        for (auto const& k : icomp->kprocs()) {
            for (auto const& spec : defsr()->updColl_I()) {
                if (k->depSpecComp(spec, icomp))
                    updset.insert(k->schedIDX());
            }
        }

        for (auto const& ip : icomp->ipatches()) {
            for (auto const& k : ip->kprocs()) {
                for (auto const& spec : defsr()->updColl_I()) {
                    if (k->depSpecComp(spec, icomp))
                        updset.insert(k->schedIDX());
                }
            }
        }

        for (auto const& op : icomp->opatches()) {
            for (auto const& k : op->kprocs()) {
                for (auto const& spec : defsr()->updColl_I()) {
                    if (k->depSpecComp(spec, icomp))
                        updset.insert(k->schedIDX());
                }
            }
        }
    }

    if (ocomp != nullptr)
    {
        for (auto const& k : ocomp->kprocs()) {
            for (auto const& spec : defsr()->updColl_O()) {
                if (k->depSpecComp(spec, ocomp))
                    updset.insert(k->schedIDX());
            }
        }

        for (auto const& ip : ocomp->ipatches()) {
            for (auto const& k : ip->kprocs()) {
                for (auto const& spec : defsr()->updColl_O()) {
                    if (k->depSpecComp(spec, ocomp))
                        updset.insert(k->schedIDX());
                }
            }
        }

        for (auto const& op : ocomp->opatches())
        {
            for (auto const& k : op->kprocs()) {
                for (auto const& spec : defsr()->updColl_O()) {
                    if (k->depSpecComp(spec, ocomp))
                        updset.insert(k->schedIDX());
                }
            }
        }
    }

    swmrssa::schedIDXSet_To_Vec(updset, pUpdVec);
}

////////////////////////////////////////////////////////////////////////////////

bool swmrssa::SReac::depSpecComp(uint gidx, swmrssa::Comp * comp)
{
    if (comp == pPatch->iComp())
    {
        return (defsr()->dep_I(gidx) != ssolver::DEP_NONE);
    }
    else if (comp == pPatch->oComp())
    {
        return (defsr()->dep_O(gidx) != ssolver::DEP_NONE);
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool swmrssa::SReac::depSpecPatch(uint gidx, swmrssa::Patch * patch)
{
    if (patch != pPatch) { return false;
}
    return (defsr()->dep_S(gidx) != ssolver::DEP_NONE);
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::SReac::reset()
{
    resetExtent();
    uint lsridx = pPatch->def()->sreacG2L(defsr()->gidx());
    pPatch->def()->setActive(lsridx, true);
    resetCcst();
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::SReac::resetCcst()
{
    uint lsridx = pPatch->def()->sreacG2L(defsr()->gidx());
    double kcst = pPatch->def()->kcst(lsridx);

    if (defsr()->surf_surf() == false)
    {
        double vol;
        if (defsr()->inside())
        {
            AssertLog(pPatch->iComp() != nullptr);
            vol = pPatch->iComp()->def()->vol();
        }
        else
        {
            assert (pPatch->oComp() != nullptr);
            vol = pPatch->oComp()->def()->vol();
        }

        pCcst = comp_ccst_vol(kcst, vol, defsr()->order());
    }
    else
    {
        double area;
        area = pPatch->def()->area();
        pCcst = comp_ccst_area(kcst, area, defsr()->order());
    }

    assert (pCcst >= 0);
}

////////////////////////////////////////////////////////////////////////////////

double swmrssa::SReac::rate(steps::wmrssa::PropensityRSSA prssa)
{
    if (inactive()) { return 0.0;
}

    // First we compute the combinatorial part.
    //   1/ for the surface part of the stoichiometry
    //   2/ for the inner or outer volume part of the stoichiometry, pool
    //      depending on whether the sreac is inner() or outer()
    // Then we multiply with mesoscopic constant.

    if (prssa == steps::wmrssa::BOUNDS) {
        pPropensityLB = rate(steps::wmrssa::LOWERBOUND);
}

    ssolver::Patchdef * pdef = pPatch->def();
    uint lidx = pdef->sreacG2L(defsr()->gidx());

    double h_mu = 1.0;

    uint * lhs_s_vec = pdef->sreac_lhs_S_bgn(lidx);
    double * cnt_s_vec = pPatch->pools(prssa);
    uint nspecs_s = pdef->countSpecs();
    for (uint s = 0; s < nspecs_s; ++s)
    {
        uint lhs = lhs_s_vec[s];
        if (lhs == 0) { continue;
}
        auto cnt = static_cast<uint>(cnt_s_vec[s]);
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
                AssertLog(false);
            }
        }
    }

    if (defsr()->inside())
    {
        uint * lhs_i_vec = pdef->sreac_lhs_I_bgn(lidx);
        double * cnt_i_vec = pPatch->iComp()->pools(prssa);
        uint nspecs_i = pdef->countSpecs_I();
        for (uint s = 0; s < nspecs_i; ++s)
        {
            uint lhs = lhs_i_vec[s];
            if (lhs == 0) { continue;
}
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
                    AssertLog(false);
                }
            }
        }
    }
    else if (defsr()->outside())
    {
        uint * lhs_o_vec = pdef->sreac_lhs_O_bgn(lidx);
        double * cnt_o_vec = pPatch->oComp()->pools(prssa);
        uint nspecs_o = pdef->countSpecs_O();
        for (uint s = 0; s < nspecs_o; ++s)
        {
            uint lhs = lhs_o_vec[s];
            if (lhs == 0) { continue;
}
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
                    AssertLog(false);
                }
            }
        }
    }

    return h_mu * pCcst;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> const & swmrssa::SReac::apply()
{
    SchedIDXSet updset;
    //bool returnUpdVec = false;
    ssolver::Patchdef * pdef = pPatch->def();
    uint lidx = pdef->sreacG2L(defsr()->gidx());

    // Update patch pools.
    int * upd_s_vec = pdef->sreac_upd_S_bgn(lidx);
    double * cnt_s_vec = pdef->pools();
    uint nspecs_s = pdef->countSpecs();

    for (uint s = 0; s < nspecs_s; ++s)
    {
        if (pdef->clamped(s)) continue;
        int upd = upd_s_vec[s];
        if (upd == 0) continue;
        int nc = static_cast<int>(cnt_s_vec[s]) + upd;
        AssertLog(nc >= 0);
        pdef->setCount(s, static_cast<double>(nc));
        if (pPatch->isOutOfBound(s, nc))
        {
            std::vector<steps::wmrssa::KProc *> dependentReacs = pPatch->getSpecUpdKProcs(s);
            for (auto &dependentReac : dependentReacs)
              updset.insert(dependentReac->schedIDX());
            //returnUpdVec = true;
        }
    }

    // Update inner comp pools.
    Comp * icomp = pPatch->iComp();
    if (icomp != nullptr)
    {
        int * upd_i_vec = pdef->sreac_upd_I_bgn(lidx);
        double * cnt_i_vec = icomp->def()->pools();
        uint nspecs_i = pdef->countSpecs_I();
        for (uint s = 0; s < nspecs_i; ++s)
        {
            if (icomp->def()->clamped(s)) continue;
            int upd = upd_i_vec[s];
            if (upd == 0) continue;
            int nc = static_cast<int>(cnt_i_vec[s]) + upd;
            AssertLog(nc >= 0);
            icomp->def()->setCount(s, static_cast<double>(nc));
            if (icomp->isOutOfBound(s, nc))
            {
                std::vector<steps::wmrssa::KProc *> dependentReacs = icomp->getSpecUpdKProcs(s);
                for (auto &dependentReac : dependentReacs)
                  updset.insert(dependentReac->schedIDX());
                //returnUpdVec = true;
            }
        }
    }

    // Update outer comp pools.
    Comp * ocomp = pPatch->oComp();
    if (ocomp != nullptr)
    {
        int * upd_o_vec = pdef->sreac_upd_O_bgn(lidx);
        double * cnt_o_vec = ocomp->def()->pools();
        uint nspecs_o = pdef->countSpecs_O();
        for (uint s = 0; s < nspecs_o; ++s)
        {
            if (ocomp->def()->clamped(s)) continue;
            int upd = upd_o_vec[s];
            if (upd == 0) continue;
            int nc = static_cast<int>(cnt_o_vec[s]) + upd;
            AssertLog(nc >= 0);
            ocomp->def()->setCount(s, static_cast<double>(nc));
            if (ocomp->isOutOfBound(s, nc))
            {
                std::vector<steps::wmrssa::KProc *> dependentReacs = ocomp->getSpecUpdKProcs(s);
                for (auto &dependentReac : dependentReacs)
                  updset.insert(dependentReac->schedIDX());
                //returnUpdVec = true;
            }
        }
    }

    rExtent++;
    swmrssa::schedIDXSet_To_Vec(updset, pUpdVec);
    return pUpdVec; //returnUpdVec ? pUpdVec : emptyVec;
}

////////////////////////////////////////////////////////////////////////////////

// END

