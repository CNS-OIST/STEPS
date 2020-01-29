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
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>


// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/math/constants.hpp"
#include "steps/solver/compdef.hpp"
#include "steps/solver/patchdef.hpp"
#include "steps/solver/reacdef.hpp"
#include "steps/solver/sreacdef.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/types.hpp"
#include "steps/wmrk4/wmrk4.hpp"

// logging
#include "easylogging++.h"

namespace swmrk4 = steps::wmrk4;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

swmrk4::Wmrk4::Wmrk4(steps::model::Model *m, steps::wm::Geom *g, const rng::RNGptr &r)
: API(m, g, r)
, pSpecs_tot(0)
, pReacs_tot(0)
, pVals()
, pSFlags()
, pNewVals()
, pDyDx()
, pDT(0.0)
, yt()
, dyt()
, dym()
{
    uint nspecstot = 0;
    uint nreacstot = 0;
    for(uint i=0; i< statedef().countComps(); ++i)
    {
        nspecstot += statedef().compdef(i)->countSpecs();
        nreacstot += statedef().compdef(i)->countReacs();
    }
    for(uint i=0; i< statedef().countPatches(); ++i)
    {
        nspecstot += statedef().patchdef(i)->countSpecs();
        nreacstot += statedef().patchdef(i)->countSReacs();
    }

    _setup();
}

///////////////////////////////////////////////////////////////////////////////

swmrk4::Wmrk4::~Wmrk4()
= default;

///////////////////////////////////////////////////////////////////////////////

std::string swmrk4::Wmrk4::getSolverName() const
{
    return "wmrk4";
}

///////////////////////////////////////////////////////////////////////////////

std::string swmrk4::Wmrk4::getSolverDesc() const
{
    return "Runge-Kutta Method in well-mixed conditions";
}

///////////////////////////////////////////////////////////////////////////////

std::string swmrk4::Wmrk4::getSolverAuthors() const
{
    return "Sam Melchior, Iain Hepburn and Stefan Wils";
}

///////////////////////////////////////////////////////////////////////////////

std::string swmrk4::Wmrk4::getSolverEmail() const
{
    return "ihepburn@oist.jp";
}


///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::reset()
{
    uint comps = statedef().countComps();
    for (uint i=0; i < comps; ++i) statedef().compdef(i)->reset();
    uint patches = statedef().countPatches();
    for (uint i=0; i < patches; ++i) statedef().patchdef(i)->reset();
    statedef().resetTime();
    // recompute flags and counts vectors in Wmrk4 object
    _refill();

}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::run(double endtime)
{
    if (endtime < statedef().time())
    {
        std::ostringstream os;
        os << "Endtime is before current simulation time";
        ArgErrLog(os.str());
    }
    _rksteps(statedef().time(), endtime);
    statedef().setTime(endtime);
}

////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::advance(double adv)
{
    if (adv < 0.0)
    {
        std::ostringstream os;
        os << "Time to advance cannot be negative";
        ArgErrLog(os.str());
    }

    double endtime = statedef().time() + adv;
    run(endtime);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::step()
{
    AssertLog(pDT > 0.0);
    _rksteps(statedef().time(), statedef().time() + pDT);
    statedef().setTime(statedef().time() + pDT);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::setRk4DT(double dt)
{
    if (dt < 0.0)
    {
        std::ostringstream os;
        os << "Time step cannot be negative or zero.";
        ArgErrLog(os.str());
    }
    pDT = dt;
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::getTime() const
{
    return statedef().time();
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::checkpoint(std::string const & file_name)
{
    std::fstream cp_file;

    cp_file.open(file_name.c_str(),
                std::fstream::out | std::fstream::binary | std::fstream::trunc);
    double state_buffer[3];
    state_buffer[0] = static_cast<double>(pSpecs_tot);
    state_buffer[1] = static_cast<double>(pReacs_tot);
    state_buffer[2] = pDT;

    cp_file.write(reinterpret_cast<char*>(&state_buffer), sizeof(double) * 3);

    cp_file.write(reinterpret_cast<char*>(&pVals.front()), sizeof(double) * pVals.size());

    cp_file.write(reinterpret_cast<char*>(&pSFlags.front()), sizeof(uint) * pSFlags.size());

    cp_file.write(reinterpret_cast<char*>(&pNewVals.front()), sizeof(double) * pNewVals.size());

    cp_file.write(reinterpret_cast<char*>(&pDyDx.front()), sizeof(double) * pDyDx.size());

    cp_file.write(reinterpret_cast<char*>(&yt.front()), sizeof(double) * yt.size());

    cp_file.write(reinterpret_cast<char*>(&dyt.front()), sizeof(double) * dyt.size());

    cp_file.write(reinterpret_cast<char*>(&dym.front()), sizeof(double) * dym.size());

    statedef().checkpoint(cp_file);

    cp_file.close();
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::restore(std::string const & file_name)
{
    std::fstream cp_file;

    cp_file.open(file_name.c_str(),
                std::fstream::in | std::fstream::binary);

    cp_file.seekg(0);
    double state_buffer[3];
    cp_file.read(reinterpret_cast<char*>(&state_buffer), sizeof(double) * 3);

    if (static_cast<uint>(state_buffer[0]) != pSpecs_tot) {
        std::ostringstream os;
        os << "checkpoint data mismatch with simulator parameters: pSpecs_tot.";
        ArgErrLog(os.str());
    }

    if (static_cast<uint>(state_buffer[1]) != pReacs_tot) {
        std::ostringstream os;
        os << "checkpoint data mismatch with simulator parameters: pReacs_tot.";
        ArgErrLog(os.str());
    }

    pDT = state_buffer[2];

    cp_file.read(reinterpret_cast<char*>(&pVals.front()), sizeof(double) * pVals.size());

    cp_file.read(reinterpret_cast<char*>(&pSFlags.front()), sizeof(uint) * pSFlags.size());

    cp_file.read(reinterpret_cast<char*>(&pNewVals.front()), sizeof(double) * pNewVals.size());

    cp_file.read(reinterpret_cast<char*>(&pDyDx.front()), sizeof(double) * pDyDx.size());

    cp_file.read(reinterpret_cast<char*>(&yt.front()), sizeof(double) * yt.size());

    cp_file.read(reinterpret_cast<char*>(&dyt.front()), sizeof(double) * dyt.size());

    cp_file.read(reinterpret_cast<char*>(&dym.front()), sizeof(double) * dym.size());

    statedef().restore(cp_file);

    cp_file.close();
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getCompVol(uint cidx) const
{
    AssertLog(cidx < statedef().countComps());
    Compdef * comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    return comp->vol();
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setCompVol(uint cidx, double vol)
{
    AssertLog(cidx < statedef().countComps());
    Compdef * comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    comp->setVol(vol);

    // recompute the scaled reaction constants
    _refillCcst();
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getCompCount(uint cidx, uint sidx) const
{
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    Compdef * comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    uint slidx = comp->specG2L(sidx);
    if (slidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }
    return comp->pools()[slidx];
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setCompCount(uint cidx, uint sidx, double n)
{
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(n >= 0.0);
    Compdef * comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    uint slidx = comp->specG2L(sidx);
    if (slidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }
    comp->setCount(slidx, n);
    // easier to recompute all counts with _refill method
    _refill();                    /// may be a better way of doing this
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getCompAmount(uint cidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getCompCount(cidx, sidx);
    return count / steps::math::AVOGADRO;
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setCompAmount(uint cidx, uint sidx, double a)
{
    AssertLog(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * steps::math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setCompCount(cidx, sidx, a2);
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getCompConc(uint cidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getCompCount(cidx, sidx);
    Compdef * comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    double vol = comp->vol();
    return count/ (1.0e3 * vol * steps::math::AVOGADRO);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setCompConc(uint cidx, uint sidx, double c)
{
    AssertLog(c >= 0.0);
    AssertLog(cidx < statedef().countComps());
    Compdef * comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    double count = c * (1.0e3 * comp->vol() * steps::math::AVOGADRO);
    // the following method does all the necessary argument checking
    _setCompCount(cidx, sidx, count);
}

///////////////////////////////////////////////////////////////////////////////

bool swmrk4::Wmrk4::_getCompClamped(uint cidx, uint sidx) const
{
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    Compdef * comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    uint lsidx = comp->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    return comp->clamped(lsidx);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setCompClamped(uint cidx, uint sidx, bool b)
{
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    Compdef * comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    uint lsidx = comp->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    comp->setClamped(lsidx, b);

    // copy flags to this solver
    _refill();
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getCompReacK(uint cidx, uint ridx) const
{
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    Compdef * comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    uint lridx = comp->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    return comp->kcst(lridx);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setCompReacK(uint cidx, uint ridx, double kf)
{
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    AssertLog(kf >= 0.0);
    Compdef * comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    uint lridx = comp->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    comp->setKcst(lridx, kf);

    // recompute the reaction constants
    _refillCcst();
}

///////////////////////////////////////////////////////////////////////////////

bool swmrk4::Wmrk4::_getCompReacActive(uint cidx, uint ridx) const
{
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    Compdef * comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    uint lridx = comp->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    return (comp->active(lridx));
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setCompReacActive(uint cidx, uint ridx, bool a)
{
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    Compdef * comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    uint lridx = comp->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    comp->setActive(lridx, a);

    // copy flags to this solver
    _refill();
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getPatchArea(uint pidx) const
{
    AssertLog(pidx < statedef().countPatches());
    Patchdef * patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    return patch->area();
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setPatchArea(uint pidx, double area)
{
    AssertLog(pidx < statedef().countPatches());
    Patchdef * patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    patch->setArea(area);
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getPatchCount(uint pidx, uint sidx) const
{
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countSpecs());
    Patchdef * patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    uint slidx = patch->specG2L(sidx);
    if (slidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    return patch->pools()[slidx];
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setPatchCount(uint pidx, uint sidx, double n)
{
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx< statedef().countSpecs());
    Patchdef * patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    uint slidx = patch->specG2L(sidx);
    if (slidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }
    patch->setCount(slidx, n);
    // easier to recompute all counts with _refill method
    _refill();
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getPatchAmount(uint pidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getPatchCount(pidx, sidx);
    return (count / steps::math::AVOGADRO);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setPatchAmount(uint pidx, uint sidx, double a)
{
    AssertLog(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * steps::math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setPatchCount(pidx, sidx, a2);
}

///////////////////////////////////////////////////////////////////////////////

bool swmrk4::Wmrk4::_getPatchClamped(uint pidx, uint sidx) const
{
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countSpecs());
    Patchdef * patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    uint lsidx = patch->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    return patch->clamped(lsidx);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setPatchClamped(uint pidx, uint sidx, bool buf)
{
    AssertLog(pidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    Patchdef * patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    uint lsidx = patch->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    patch->setClamped(lsidx, buf);

    // copy flags to this solver
    _refill();
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_getPatchSReacK(uint pidx, uint ridx) const
{
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    Patchdef * patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    uint lridx = patch->sreacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    return patch->kcst(lridx);
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setPatchSReacK(uint pidx, uint ridx, double kf)
{
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    AssertLog(kf >= 0.0);
    Patchdef * patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    uint lridx = patch->sreacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    patch->setKcst(lridx, kf);

    // recompute the reaction constants
    _refillCcst();
}

///////////////////////////////////////////////////////////////////////////////

bool swmrk4::Wmrk4::_getPatchSReacActive(uint pidx, uint ridx) const
{
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    Patchdef * patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    uint lridx = patch->sreacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    return (patch->active(lridx));
}

///////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setPatchSReacActive(uint pidx, uint ridx, bool a)
{
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    Patchdef * patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    uint lridx = patch->sreacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    patch->setActive(lridx, a);

    // copy flags to this solver
    _refill();
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_ccst(double kcst, double vol, uint order)
{
    double vscale = 1.0e3 * vol * steps::math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;

    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;
    return kcst * pow(vscale, static_cast<double>(-o1));
}

///////////////////////////////////////////////////////////////////////////////

double swmrk4::Wmrk4::_ccst2D(double kcst, double area, uint order)
{
    double vscale = area * steps::math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;
    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;
    return kcst * pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setup()

{
    /// cumulative number of species from each compartment and patch
    pSpecs_tot = 0;
    /// cumulative number of reacs from each compartment and patch
    pReacs_tot = 0;

    uint Comps_N = statedef().countComps();
    uint Patches_N = statedef().countPatches();

    for(uint i=0; i< Comps_N; ++i)
    {
        pSpecs_tot += statedef().compdef(i)->countSpecs();
        pReacs_tot += statedef().compdef(i)->countReacs();
    }
    for(uint i=0; i< Patches_N; ++i)
    {
        pSpecs_tot += statedef().patchdef(i)->countSpecs();
        pReacs_tot += statedef().patchdef(i)->countSReacs();
    }
    AssertLog(pSpecs_tot > 0);
    AssertLog(pReacs_tot > 0);

    for(uint i=0; i< pSpecs_tot; ++i)
    {
        pVals.push_back(0.0);
        pSFlags.push_back(0);
        pNewVals.push_back(0.0);
        pDyDx.push_back(0.0);
        yt.push_back(0.0);
        dyt.push_back(0.0);
        dym.push_back(0.0);
    }

    /// fill the reaction matrix
    /// loop over compartments,
    /// then comp reacs and copy compdef LHS values to correct index

    /// set row marker to beginning of matrix for first compartment
    uint rowp = 0;
    /// set column marker to beginning of matrix for first compartment
    uint colp = 0;

    for (uint i=0; i< Comps_N; ++i)
    {
        uint compReacs_N = statedef().compdef(i)->countReacs();
        uint compSpecs_N = statedef().compdef(i)->countSpecs();

        for(uint j=0; j< compReacs_N; ++j)
        {
            Reaction reaction;
            for(uint k=0; k< compSpecs_N; ++k)
            {
                uint lhs = statedef().compdef(i)->reac_lhs_bgn(j)[k];
                int upd = statedef().compdef(i)->reac_upd_bgn(j)[k];
                reaction.addSpecies(colp + k, lhs, upd);
            }
            /// set scaled reaction constant
            double reac_kcst = statedef().compdef(i)->kcst(j);
            double comp_vol = statedef().compdef(i)->vol();
            uint reac_order = statedef().compdef(i)->reacdef(j)->order();
            reaction.c = _ccst(reac_kcst, comp_vol, reac_order);
            reactions.push_back(reaction);
        }
        /// step up markers for next compartment
        rowp += compReacs_N;
        colp += compSpecs_N;
    }

    /// now loop over patches,
    /// then sreacs, filling for correct inner and outer compartments
    for(uint i=0; i< Patches_N; ++i)
    {
        const auto patch = statedef().patchdef(i);
        uint patchReacs_N = patch->countSReacs();
        uint patchSpecs_N_S = patch->countSpecs();
        uint patchSpecs_N_I = patch->countSpecs_I();
        uint patchSpecs_N_O = patch->countSpecs_O();

        for (uint j=0; j< patchReacs_N; ++j)
        {
            Reaction reaction;
            for(uint k=0; k< patchSpecs_N_S; ++k)
            {   uint slhs = patch->sreac_lhs_S_bgn(j)[k];
                int supd = patch->sreac_upd_S_bgn(j)[k];
                reaction.addSpecies(colp + k, slhs, supd);
            }

            /// fill for inner and outer compartments involved in sreac j
            /// only perform if inner comp exists, similarly for outer comp
            if (patch->icompdef() != nullptr)
            {
                /// fetch global index of inner compartment
                uint icompidx = patch->icompdef()->gidx();
                // marker for correct position of inner compartment in matrix
                uint mtx_icompidx = 0;
                /// step up marker to correct comp
                for (uint l=0; l< icompidx; ++l)
                {
                    mtx_icompidx += statedef().compdef(l)->countSpecs();
                }
                for(uint k=0; k< patchSpecs_N_I; ++k)
                {
                    uint ilhs = patch->sreac_lhs_I_bgn(j)[k];
                    int iupd = patch->sreac_upd_I_bgn(j)[k];
                    reaction.addSpecies(mtx_icompidx + k, ilhs, iupd);
                }
            }
            if (patch->ocompdef() != nullptr)
            {
                uint ocompidx = patch->ocompdef()->gidx();
                uint mtx_ocompidx =0;
                for (uint l=0; l< ocompidx; ++l)
                {
                    mtx_ocompidx += statedef().compdef(l)->countSpecs();
                }
                for(uint k=0; k< patchSpecs_N_O; ++k)
                {
                    uint olhs = patch->sreac_lhs_O_bgn(j)[k];
                    int oupd = patch->sreac_upd_O_bgn(j)[k];
                    reaction.addSpecies(mtx_ocompidx + k, olhs, oupd);
                }
            }
            if (patch->sreacdef(j)->surf_surf() == false)
            {    /// set scaled reaction constant
                /// depends on volume of lhs reaction compartment
                double vol;
                if (patch->sreacdef(j)->inside())
                {
                  const auto comp = patch->icompdef();
                  AssertLog(comp != nullptr);
                  vol = comp->vol();
                }
                else
                {
                    const auto comp = patch->ocompdef();
                    AssertLog(comp != nullptr);
                    vol = comp->vol();
                }
                double sreac_kcst = patch->kcst(j);
                uint sreac_order = patch->sreacdef(j)->order();
                reaction.c = _ccst(sreac_kcst, vol, sreac_order);
            }
            else
            {
                /// 2D reaction
                double area = patch->area();
                double sreac_kcst = patch->kcst(j);
                uint sreac_order = patch->sreacdef(j)->order();
                reaction.c = _ccst2D(sreac_kcst, area, sreac_order);
            }
            reactions.push_back(reaction);

        }
        /// move markers to next point in matrix
        rowp += patchReacs_N;
        colp += patchSpecs_N_S;
    }

    AssertLog(reactions.size() == pReacs_tot);
    _refill();
    _refillCcst();
}

////////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_refill()
{
    uint Comps_N = statedef().countComps();
    uint Patches_N = statedef().countPatches();
    AssertLog(Comps_N > 0);

    uint c_marker = 0;
    uint r_marker = 0;

    for(uint i=0; i< Comps_N; ++i)
    {
        uint comp_Specs_N = statedef().compdef(i)->countSpecs();
        uint comp_Reacs_N = statedef().compdef(i)->countReacs();
        Compdef * comp = statedef().compdef(i);
        AssertLog(comp != nullptr);
        for (uint j=0; j < comp_Specs_N; ++j)
        {
            pVals[c_marker + j] = comp->pools()[j];
            pSFlags[c_marker + j] = comp->flags()[j];
        }
        for(uint k=0; k< comp_Reacs_N; ++k)
        {
            reactions[r_marker + k].isActivated = comp->active(k);
        }
        c_marker += comp_Specs_N;
        r_marker += comp_Reacs_N;
    }

    for(uint i=0; i< Patches_N; ++i)
    {
        uint patch_Specs_N = statedef().patchdef(i)->countSpecs();
        uint patch_Reacs_N = statedef().patchdef(i)->countSReacs();
        Patchdef * patch = statedef().patchdef(i);
        AssertLog(patch != nullptr);
        for (uint j=0; j< patch_Specs_N; ++j)
        {
            pVals[c_marker +j] = patch->pools()[j];
            pSFlags[c_marker + j] = patch->flags()[j];
        }
        for (uint k=0; k< patch_Reacs_N; ++k)
        {
            reactions[r_marker + k].isActivated = patch->active(k);
        }
        c_marker += patch_Specs_N;
        r_marker += patch_Reacs_N;
    }

    AssertLog(c_marker == pVals.size());
    AssertLog(pVals.size() == pSFlags.size());
    AssertLog(pSFlags.size() == pSpecs_tot);
}

////////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_refillCcst()
{
    uint Comps_N = statedef().countComps();
    uint Patches_N = statedef().countPatches();
    AssertLog(Comps_N > 0);

    uint r_marker = 0;

    for (uint i=0; i< Comps_N; ++i)
    {
        uint compReacs_N = statedef().compdef(i)->countReacs();
        // uint compSpecs_N = statedef().compdef(i)->countSpecs();

        for(uint j=0; j< compReacs_N; ++j)
        {
            /// set scaled reaction constant
            // DEBUG 8/4/09: reaction constants were found from model level objects
            // so didn't take into account sim-level changes
            double reac_kcst = statedef().compdef(i)->kcst(j);
            double comp_vol = statedef().compdef(i)->vol();
            uint reac_order = statedef().compdef(i)->reacdef(j)->order();
            reactions[r_marker + j].c = _ccst(reac_kcst, comp_vol, reac_order);
        }
        r_marker += compReacs_N;
    }

    /// now loop over patches,
    /// then sreacs, filling for correct inner and outer compartments
    for(uint i=0; i< Patches_N; ++i)
    {
        uint patchReacs_N = statedef().patchdef(i)->countSReacs();

        for (uint j=0; j< patchReacs_N; ++j)
        {
            if (statedef().patchdef(i)->sreacdef(j)->surf_surf() == false)
            {
                /// set scaled reaction constant
                /// depends on volume of lhs reaction compartment
                double vol;
                if (statedef().patchdef(i)->sreacdef(j)->inside())
                {
                    AssertLog(statedef().patchdef(i)->icompdef() != nullptr);
                    vol = statedef().patchdef(i)->icompdef()->vol();
                }
                else
                {
                    AssertLog(statedef().patchdef(i)->ocompdef() != nullptr);
                    vol = statedef().patchdef(i)->ocompdef()->vol();
                }
                // DEBUG 8/4/09: reaction constants were found from model level objects
                // so didn't take into account sim-level changes
                double sreac_kcst = statedef().patchdef(i)->kcst(j);
                uint sreac_order = statedef().patchdef(i)->sreacdef(j)->order();
                reactions[r_marker + j].c = _ccst(sreac_kcst, vol, sreac_order);
                }
            else
            {
                /// 2D reaction
                double area = statedef().patchdef(i)->area();
                double sreac_kcst = statedef().patchdef(i)->kcst(j);
                uint sreac_order = statedef().patchdef(i)->sreacdef(j)->order();
                reactions[r_marker + j].c = _ccst2D(sreac_kcst, area, sreac_order);
            }
        }
        r_marker += patchReacs_N;
    }

}

////////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_setderivs(dVec & vals, dVec & dydx)
{
    std::fill(dydx.begin(), dydx.end(), 0);
    for(const auto &reaction: reactions)
    {
        if (reaction.isActivated)
        {
            double numberOfReactionFirings = reaction.c;
            for (const auto &reactant: reaction.reactants)
            {
                const double population = vals[reactant.globalIndex];
                switch (reactant.order)
                {
                    case 4: numberOfReactionFirings *= population;
                    case 3: numberOfReactionFirings *= population;
                    case 2: numberOfReactionFirings *= population;
                    case 1: numberOfReactionFirings *= population;
                    case 0: break;
                    /// allow maximum 4 molecules of one species in reaction
                    default: AssertLog(0);
                }
            }
            for (const auto &specie: reaction.affectedSpecies)
            {
                if (pSFlags[specie.globalIndex] & Statedef::CLAMPED_POOLFLAG)
                    continue;
                dydx[specie.globalIndex] += specie.populationChange*numberOfReactionFirings;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_rk4(double pdt)
{
    double dt_2 = pdt/2.0;
    double dt_6 = pdt/6.0;

    for(uint i=0; i< pSpecs_tot; ++i) yt[i] = pVals[i] + (dt_2 * pDyDx[i]);
    _setderivs(yt, dyt);
    for(uint i =0; i< pSpecs_tot; ++i) yt[i]= pVals[i] + (dt_2 * dyt[i]);
    _setderivs(yt, dym);
    for (uint i=0; i< pSpecs_tot; ++i)
    {
        yt[i] = pVals[i] + (pdt * dym[i]);
        dym[i] += dyt[i];
    }
    _setderivs(yt, dyt);
    for (uint i=0; i< pSpecs_tot; ++i)
    {
        pNewVals[i]=pVals[i]+dt_6*(pDyDx[i]+dyt[i]+(2.0*dym[i]));
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_rksteps(double t1, double t2)
{
    if (t1 == t2) return;
    AssertLog(t1 < t2);
    double t = t1;
    if (pDT <= 0.0)
    {
        std::ostringstream os;
        os << "dt is zero or negative. Call setDT() method.";
        ArgErrLog(os.str());
    }
    // Bugfix - step() calls this function with one dt step, so these two values can be equal
    if (pDT > (t2-t1))
    {
        std::ostringstream os;
        os << "dt is larger than simulation step.";
        ArgErrLog(os.str());
    }

    /// step up until over maximum time
    while(t < t2)
    {
        if ((t+pDT) > t2) break;

        _setderivs(pVals, pDyDx);
        _rk4(pDT);
        _update();
        t += pDT;
    }

    ////////////////////////////////////////////////////////////////////////////
    // DEBUG: 25/03/09. This is general fix for this solver, inspired by a
    // problem with evaluating two supposedly equal double values as non-equal.
    // Now any discrepancy between simulation time and the end step time is dealt
    // with by solving to the simulation time by passing in the fraction (see below)
    // Changed _rk4() to take the dt as a double argument
    double tfrac = t2-t;
    AssertLog(tfrac >= 0.0);

    // Lets only concern ourselves with fractions greater than 1 percent
    if (tfrac != 0.0) // && tfrac/pDT >= 0.01)
    {
        AssertLog(tfrac < pDT);
        _setderivs(pVals, pDyDx);
        _rk4(tfrac);
        _update();
    }
    ////////////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////

void swmrk4::Wmrk4::_update()
{
    /// update local values vector with computed counts
    for (uint i=0; i< pSpecs_tot; ++i)
    {
        /// check clamped flag and only update if not clamped
        if (pSFlags[i] & Statedef::CLAMPED_POOLFLAG)
        {
            continue;
        }
        else
        {
            double newval = pNewVals[i];
            if (newval < 0.0) newval = 0.0;
            pVals[i] = newval;
        }
    }

    /// update pools with computed values
    uint Comps_N = statedef().countComps();
    uint Patches_N = statedef().countPatches();
    uint c_marker = 0;

    for(uint i=0; i< Comps_N; ++i)
    {
        uint comp_Specs_N = statedef().compdef(i)->countSpecs();
        for (uint j=0; j< comp_Specs_N; ++j)
        {
            statedef().compdef(i)->setCount(j, pVals[c_marker + j]);
        }
        c_marker += comp_Specs_N;
    }

    for(uint i=0; i< Patches_N; ++i)
    {
        uint patch_Specs_N = statedef().patchdef(i)->countSpecs();
        for (uint j=0; j< patch_Specs_N; ++j)
        {
            statedef().patchdef(i)->setCount(j, pVals[c_marker + j]);
        }
        c_marker += patch_Specs_N;
    }
}

////////////////////////////////////////////////////////////////////////////////

// END

