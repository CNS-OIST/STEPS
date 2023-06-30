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

// STEPS headers.
#include "wmrk4.hpp"

#include "math/constants.hpp"
#include "solver/compdef.hpp"
#include "solver/patchdef.hpp"
#include "solver/reacdef.hpp"
#include "solver/sreacdef.hpp"
#include "solver/statedef.hpp"
#include "solver/types.hpp"

// logging
#include "util/error.hpp"
#include <easylogging++.h>

#include "util/checkpointing.hpp"

namespace steps::wmrk4 {

Wmrk4::Wmrk4(model::Model* m, wm::Geom* g, const rng::RNGptr& r)
    : API(m, g, r)
    , pSpecs_tot(0)
    , pReacs_tot(0)
    , pDT(0.0) {
    _setup();
}

///////////////////////////////////////////////////////////////////////////////

Wmrk4::~Wmrk4() = default;

///////////////////////////////////////////////////////////////////////////////

std::string Wmrk4::getSolverName() const {
    return "wmrk4";
}

///////////////////////////////////////////////////////////////////////////////

std::string Wmrk4::getSolverDesc() const {
    return "Runge-Kutta Method in well-mixed conditions";
}

///////////////////////////////////////////////////////////////////////////////

std::string Wmrk4::getSolverAuthors() const {
    return "Sam Melchior, Iain Hepburn and Stefan Wils";
}

///////////////////////////////////////////////////////////////////////////////

std::string Wmrk4::getSolverEmail() const {
    return "ihepburn@oist.jp";
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::reset() {
    uint comps = statedef().countComps();
    for (auto i: solver::comp_global_id::range(comps)) {
        statedef().compdef(i)->reset();
    }
    uint patches = statedef().countPatches();
    for (auto i: solver::patch_global_id::range(patches)) {
        statedef().patchdef(i)->reset();
    }
    statedef().resetTime();
    // recompute flags and counts vectors in Wmrk4 object
    _refill();
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::run(double endtime) {
    if (endtime < statedef().time()) {
        std::ostringstream os;
        os << "Endtime is before current simulation time";
        ArgErrLog(os.str());
    }
    _rksteps(statedef().time(), endtime);
    statedef().setTime(endtime);
}

////////////////////////////////////////////////////////////////////////

void Wmrk4::advance(double adv) {
    if (adv < 0.0) {
        std::ostringstream os;
        os << "Time to advance cannot be negative";
        ArgErrLog(os.str());
    }

    double endtime = statedef().time() + adv;
    run(endtime);
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::step() {
    AssertLog(pDT > 0.0);
    _rksteps(statedef().time(), statedef().time() + pDT);
    statedef().setTime(statedef().time() + pDT);
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::setRk4DT(double dt) {
    if (dt < 0.0) {
        std::ostringstream os;
        os << "Time step cannot be negative or zero.";
        ArgErrLog(os.str());
    }
    pDT = dt;
}

///////////////////////////////////////////////////////////////////////////////

double Wmrk4::getTime() const {
    return statedef().time();
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::checkpoint(std::string const& file_name) {
    std::fstream cp_file;

    cp_file.open(file_name.c_str(), std::fstream::out | std::fstream::binary | std::fstream::trunc);

    API::checkpoint(cp_file);

    util::checkpoint(cp_file, pSpecs_tot);
    util::checkpoint(cp_file, pReacs_tot);

    util::checkpoint(cp_file, pDT);
    util::checkpoint(cp_file, pVals);
    util::checkpoint(cp_file, pSFlags);
    util::checkpoint(cp_file, pNewVals);
    util::checkpoint(cp_file, pDyDx);
    util::checkpoint(cp_file, yt);
    util::checkpoint(cp_file, dyt);
    util::checkpoint(cp_file, dym);

    statedef().checkpoint(cp_file);

    for (auto& reaction: reactions) {
        reaction.checkpoint(cp_file);
    }

    cp_file.close();
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::restore(std::string const& file_name) {
    std::fstream cp_file;

    cp_file.open(file_name.c_str(), std::fstream::in | std::fstream::binary);

    cp_file.seekg(0);

    API::restore(cp_file);

    util::compare(cp_file, pSpecs_tot, "Mismatched pSpecs_tot restore value.");
    util::compare(cp_file, pReacs_tot, "Mismatched pReacs_tot restore value.");

    util::restore(cp_file, pDT);
    util::restore(cp_file, pVals);
    util::restore(cp_file, pSFlags);
    util::restore(cp_file, pNewVals);
    util::restore(cp_file, pDyDx);
    util::restore(cp_file, yt);
    util::restore(cp_file, dyt);
    util::restore(cp_file, dym);
    statedef().restore(cp_file);

    for (auto& reaction: reactions) {
        reaction.restore(cp_file);
    }

    cp_file.close();
}

///////////////////////////////////////////////////////////////////////////////

double Wmrk4::_getCompVol(solver::comp_global_id cidx) const {
    AssertLog(cidx < statedef().countComps());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    return comp->vol();
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setCompVol(solver::comp_global_id cidx, double vol) {
    AssertLog(cidx < statedef().countComps());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    comp->setVol(vol);

    // recompute the scaled reaction constants
    _refillCcst();
}

///////////////////////////////////////////////////////////////////////////////

double Wmrk4::_getCompSpecCount(solver::comp_global_id cidx, solver::spec_global_id sidx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    solver::spec_local_id slidx = comp->specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }
    return comp->pools()[slidx];
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setCompSpecCount(solver::comp_global_id cidx, solver::spec_global_id sidx, double n) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(n >= 0.0);
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    solver::spec_local_id slidx = comp->specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }
    comp->setCount(slidx, n);
    // easier to recompute all counts with _refill method
    _refill();  /// may be a better way of doing this
}

///////////////////////////////////////////////////////////////////////////////

double Wmrk4::_getCompSpecAmount(solver::comp_global_id cidx, solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getCompSpecCount(cidx, sidx);
    return count / math::AVOGADRO;
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setCompSpecAmount(solver::comp_global_id cidx, solver::spec_global_id sidx, double a) {
    AssertLog(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setCompSpecCount(cidx, sidx, a2);
}

///////////////////////////////////////////////////////////////////////////////

double Wmrk4::_getCompSpecConc(solver::comp_global_id cidx, solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getCompSpecCount(cidx, sidx);
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    double vol = comp->vol();
    return count / (1.0e3 * vol * math::AVOGADRO);
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setCompSpecConc(solver::comp_global_id cidx, solver::spec_global_id sidx, double c) {
    AssertLog(c >= 0.0);
    AssertLog(cidx < statedef().countComps());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    double count = c * (1.0e3 * comp->vol() * math::AVOGADRO);
    // the following method does all the necessary argument checking
    _setCompSpecCount(cidx, sidx, count);
}

///////////////////////////////////////////////////////////////////////////////

bool Wmrk4::_getCompSpecClamped(solver::comp_global_id cidx, solver::spec_global_id sidx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    solver::spec_local_id lsidx = comp->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    return comp->clamped(lsidx);
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setCompSpecClamped(solver::comp_global_id cidx, solver::spec_global_id sidx, bool b) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    solver::spec_local_id lsidx = comp->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    comp->setClamped(lsidx, b);

    // copy flags to this solver
    _refill();
}

///////////////////////////////////////////////////////////////////////////////

double Wmrk4::_getCompReacK(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = comp->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    return comp->kcst(lridx);
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setCompReacK(solver::comp_global_id cidx, solver::reac_global_id ridx, double kf) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    AssertLog(kf >= 0.0);
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = comp->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    comp->setKcst(lridx, kf);

    // recompute the reaction constants
    _refillCcst();
}

///////////////////////////////////////////////////////////////////////////////

bool Wmrk4::_getCompReacActive(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = comp->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    return (comp->active(lridx));
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setCompReacActive(solver::comp_global_id cidx, solver::reac_global_id ridx, bool a) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = comp->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    comp->setActive(lridx, a);

    // copy flags to this solver
    _refill();
}

///////////////////////////////////////////////////////////////////////////////

double Wmrk4::_getPatchArea(solver::patch_global_id pidx) const {
    AssertLog(pidx < statedef().countPatches());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    return patch->area();
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setPatchArea(solver::patch_global_id pidx, double area) {
    AssertLog(pidx < statedef().countPatches());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    patch->setArea(area);
}

///////////////////////////////////////////////////////////////////////////////

double Wmrk4::_getPatchSpecCount(solver::patch_global_id pidx, solver::spec_global_id sidx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countSpecs());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    solver::spec_local_id slidx = patch->specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    return patch->pools()[slidx];
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setPatchSpecCount(solver::patch_global_id pidx,
                               solver::spec_global_id sidx,
                               double n) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countSpecs());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    solver::spec_local_id slidx = patch->specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }
    patch->setCount(slidx, n);
    // easier to recompute all counts with _refill method
    _refill();
}

///////////////////////////////////////////////////////////////////////////////

double Wmrk4::_getPatchSpecAmount(solver::patch_global_id pidx, solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getPatchSpecCount(pidx, sidx);
    return (count / math::AVOGADRO);
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setPatchSpecAmount(solver::patch_global_id pidx,
                                solver::spec_global_id sidx,
                                double a) {
    AssertLog(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setPatchSpecCount(pidx, sidx, a2);
}

///////////////////////////////////////////////////////////////////////////////

bool Wmrk4::_getPatchSpecClamped(solver::patch_global_id pidx, solver::spec_global_id sidx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countSpecs());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    solver::spec_local_id lsidx = patch->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    return patch->clamped(lsidx);
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setPatchSpecClamped(solver::patch_global_id pidx,
                                 solver::spec_global_id sidx,
                                 bool buf) {
    AssertLog(pidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    solver::spec_local_id lsidx = patch->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    patch->setClamped(lsidx, buf);

    // copy flags to this solver
    _refill();
}

///////////////////////////////////////////////////////////////////////////////

double Wmrk4::_getPatchSReacK(solver::patch_global_id pidx, solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lridx = patch->sreacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    return patch->kcst(lridx);
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setPatchSReacK(solver::patch_global_id pidx, solver::sreac_global_id ridx, double kf) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    AssertLog(kf >= 0.0);
    solver::Patchdef* patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lridx = patch->sreacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    patch->setKcst(lridx, kf);

    // recompute the reaction constants
    _refillCcst();
}

///////////////////////////////////////////////////////////////////////////////

bool Wmrk4::_getPatchSReacActive(solver::patch_global_id pidx, solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lridx = patch->sreacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    return (patch->active(lridx));
}

///////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setPatchSReacActive(solver::patch_global_id pidx,
                                 solver::sreac_global_id ridx,
                                 bool a) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lridx = patch->sreacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    patch->setActive(lridx, a);

    // copy flags to this solver
    _refill();
}

///////////////////////////////////////////////////////////////////////////////

double Wmrk4::_ccst(double kcst, double vol, uint order) {
    double vscale = 1.0e3 * vol * math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;

    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;
    return kcst * pow(vscale, static_cast<double>(-o1));
}

///////////////////////////////////////////////////////////////////////////////

double Wmrk4::_ccst2D(double kcst, double area, uint order) {
    double vscale = area * math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;
    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;
    return kcst * pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setup()

{
    /// cumulative number of species from each compartment and patch
    pSpecs_tot = 0;
    /// cumulative number of reacs from each compartment and patch
    pReacs_tot = 0;

    uint Comps_N = statedef().countComps();
    uint Patches_N = statedef().countPatches();

    for (auto i: solver::comp_global_id::range(Comps_N)) {
        auto* compdef = statedef().compdef(i);
        pSpecs_tot += compdef->countSpecs();
        pReacs_tot += compdef->countReacs();
    }
    for (auto i: solver::patch_global_id::range(Patches_N)) {
        auto* patchdef = statedef().patchdef(i);
        pSpecs_tot += patchdef->countSpecs();
        pReacs_tot += patchdef->countSReacs();
    }
    AssertLog(pSpecs_tot > 0);
    AssertLog(pReacs_tot > 0);

    for (uint i = 0; i < pSpecs_tot; ++i) {
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

    /// set column marker to beginning of matrix for first compartment
    uint colp = 0;

    for (auto i: solver::comp_global_id::range(Comps_N)) {
        const auto compdef = statedef().compdef(i);
        const uint compReacs_N = compdef->countReacs();
        const uint compSpecs_N = compdef->countSpecs();

        for (auto j: solver::reac_local_id::range(compReacs_N)) {
            const auto& reac_lhs = compdef->reac_lhs(j);
            const auto& reac_upd = compdef->reac_upd(j);
            Reaction reaction;
            for (auto k: reac_lhs.range()) {
                uint lhs = reac_lhs[k];
                int upd = reac_upd[k];
                reaction.addSpecies(colp + k.get(), lhs, upd);
            }
            /// set scaled reaction constant
            double reac_kcst = compdef->kcst(j);
            double comp_vol = compdef->vol();
            uint reac_order = compdef->reacdef(j)->order();
            reaction.c = _ccst(reac_kcst, comp_vol, reac_order);
            reactions.push_back(reaction);
        }
        /// step up markers for next compartment
        colp += compSpecs_N;
    }

    /// now loop over patches,
    /// then sreacs, filling for correct inner and outer compartments
    for (auto i: solver::patch_global_id::range(Patches_N)) {
        const auto patch = statedef().patchdef(i);

        for (auto j: solver::sreac_local_id::range(patch->countSReacs())) {
            Reaction reaction;
            const auto& sreac_lhs_S = patch->sreac_lhs_S(j);
            const auto& sreac_upd_S = patch->sreac_upd_S(j);
            for (auto k: sreac_lhs_S.range()) {
                uint slhs = sreac_lhs_S[k];
                int supd = sreac_upd_S[k];
                reaction.addSpecies(colp + k.get(), slhs, supd);
            }

            /// fill for inner and outer compartments involved in sreac j
            /// only perform if inner comp exists, similarly for outer comp
            if (patch->icompdef() != nullptr) {
                const auto& sreac_lhs_I = patch->sreac_lhs_I(j);
                const auto& sreac_upd_I = patch->sreac_upd_I(j);

                /// fetch global index of inner compartment
                solver::comp_global_id icompidx = patch->icompdef()->gidx();
                // marker for correct position of inner compartment in matrix
                uint mtx_icompidx = 0;
                /// step up marker to correct comp
                for (auto l: icompidx.range()) {
                    mtx_icompidx += statedef().compdef(l)->countSpecs();
                }
                for (auto k: sreac_lhs_I.range()) {
                    uint ilhs = sreac_lhs_I[k];
                    int iupd = sreac_upd_I[k];
                    reaction.addSpecies(mtx_icompidx + k.get(), ilhs, iupd);
                }
            }
            if (patch->ocompdef() != nullptr) {
                const auto& sreac_lhs_O = patch->sreac_lhs_O(j);
                const auto& sreac_upd_O = patch->sreac_upd_O(j);

                solver::comp_global_id ocompidx = patch->ocompdef()->gidx();
                uint mtx_ocompidx = 0;
                for (auto l: ocompidx.range()) {
                    mtx_ocompidx += statedef().compdef(l)->countSpecs();
                }
                for (auto k: sreac_lhs_O.range()) {
                    uint olhs = sreac_lhs_O[k];
                    int oupd = sreac_upd_O[k];
                    reaction.addSpecies(mtx_ocompidx + k.get(), olhs, oupd);
                }
            }
            if (patch->sreacdef(j)->surf_surf() == false) {  /// set scaled reaction constant
                /// depends on volume of lhs reaction compartment
                double vol;
                if (patch->sreacdef(j)->inside()) {
                    const auto comp = patch->icompdef();
                    AssertLog(comp != nullptr);
                    vol = comp->vol();
                } else {
                    const auto comp = patch->ocompdef();
                    AssertLog(comp != nullptr);
                    vol = comp->vol();
                }
                double sreac_kcst = patch->kcst(j);
                uint sreac_order = patch->sreacdef(j)->order();
                reaction.c = _ccst(sreac_kcst, vol, sreac_order);
            } else {
                /// 2D reaction
                double area = patch->area();
                double sreac_kcst = patch->kcst(j);
                uint sreac_order = patch->sreacdef(j)->order();
                reaction.c = _ccst2D(sreac_kcst, area, sreac_order);
            }
            reactions.push_back(reaction);
        }
        /// move markers to next point in matrix
        colp += patch->countSpecs();
    }

    AssertLog(reactions.size() == pReacs_tot);
    _refill();
    _refillCcst();
}

////////////////////////////////////////////////////////////////////////////////

void Wmrk4::_refill() {
    uint Comps_N = statedef().countComps();
    uint Patches_N = statedef().countPatches();
    AssertLog(Comps_N > 0);

    uint c_marker = 0;
    uint r_marker = 0;

    for (auto i: solver::comp_global_id::range(Comps_N)) {
        const auto compdef = statedef().compdef(i);
        const uint comp_Specs_N = compdef->countSpecs();
        const uint comp_Reacs_N = compdef->countReacs();
        for (auto j: compdef->pools().range()) {
            pVals[c_marker + j.get()] = compdef->pools()[j];
            pSFlags[c_marker + j.get()] = compdef->flags()[j];
        }
        for (auto k: solver::reac_local_id::range(comp_Reacs_N)) {
            reactions[r_marker + k.get()].isActivated = compdef->active(k);
        }
        c_marker += comp_Specs_N;
        r_marker += comp_Reacs_N;
    }

    for (auto i: solver::patch_global_id::range(Patches_N)) {
        solver::Patchdef* patch = statedef().patchdef(i);
        uint patch_Specs_N = patch->countSpecs();
        uint patch_Reacs_N = patch->countSReacs();

        AssertLog(patch != nullptr);
        for (auto j: solver::spec_local_id::range(patch_Specs_N)) {
            pVals[c_marker + j.get()] = patch->pools()[j];
            pSFlags[c_marker + j.get()] = patch->flags()[j];
        }
        for (auto k: solver::sreac_local_id::range(patch_Reacs_N)) {
            reactions[r_marker + k.get()].isActivated = patch->active(k);
        }
        c_marker += patch_Specs_N;
        r_marker += patch_Reacs_N;
    }

    AssertLog(c_marker == pVals.size());
    AssertLog(pVals.size() == pSFlags.size());
    AssertLog(pSFlags.size() == pSpecs_tot);
}

////////////////////////////////////////////////////////////////////////////////

void Wmrk4::_refillCcst() {
    uint Comps_N = statedef().countComps();
    uint Patches_N = statedef().countPatches();
    AssertLog(Comps_N > 0);

    uint r_marker = 0;

    for (auto i: solver::comp_global_id::range(Comps_N)) {
        const auto compdef = statedef().compdef(i);
        const uint compReacs_N = compdef->countReacs();
        // uint compSpecs_N = compdef->countSpecs();

        for (auto j: solver::reac_local_id::range(compReacs_N)) {
            /// set scaled reaction constant
            // DEBUG 8/4/09: reaction constants were found from model level objects
            // so didn't take into account sim-level changes
            double reac_kcst = compdef->kcst(j);
            double comp_vol = compdef->vol();
            uint reac_order = compdef->reacdef(j)->order();
            reactions[r_marker + j.get()].c = _ccst(reac_kcst, comp_vol, reac_order);
        }
        r_marker += compReacs_N;
    }

    /// now loop over patches,
    /// then sreacs, filling for correct inner and outer compartments
    for (auto i: solver::patch_global_id::range(Patches_N)) {
        const auto& patchdef = statedef().patchdef(i);
        uint patchReacs_N = patchdef->countSReacs();

        for (auto j: solver::sreac_local_id::range(patchReacs_N)) {
            if (patchdef->sreacdef(j)->surf_surf() == false) {
                /// set scaled reaction constant
                /// depends on volume of lhs reaction compartment
                double vol;
                if (patchdef->sreacdef(j)->inside()) {
                    AssertLog(patchdef->icompdef() != nullptr);
                    vol = patchdef->icompdef()->vol();
                } else {
                    AssertLog(patchdef->ocompdef() != nullptr);
                    vol = patchdef->ocompdef()->vol();
                }
                // DEBUG 8/4/09: reaction constants were found from model level objects
                // so didn't take into account sim-level changes
                double sreac_kcst = patchdef->kcst(j);
                uint sreac_order = patchdef->sreacdef(j)->order();
                reactions[r_marker + j.get()].c = _ccst(sreac_kcst, vol, sreac_order);
            } else {
                /// 2D reaction
                double area = patchdef->area();
                double sreac_kcst = patchdef->kcst(j);
                uint sreac_order = patchdef->sreacdef(j)->order();
                reactions[r_marker + j.get()].c = _ccst2D(sreac_kcst, area, sreac_order);
            }
        }
        r_marker += patchReacs_N;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Wmrk4::_setderivs(dVec& vals, dVec& dydx) {
    std::fill(dydx.begin(), dydx.end(), 0);
    for (const auto& reaction: reactions) {
        if (reaction.isActivated) {
            double numberOfReactionFirings = reaction.c;
            for (const auto& reactant: reaction.reactants) {
                const double population = vals[reactant.globalIndex];
                switch (reactant.order) {
                case 4:
                    numberOfReactionFirings *= population;
                    STEPS_FALLTHROUGH;
                case 3:
                    numberOfReactionFirings *= population;
                    STEPS_FALLTHROUGH;
                case 2:
                    numberOfReactionFirings *= population;
                    STEPS_FALLTHROUGH;
                case 1:
                    numberOfReactionFirings *= population;
                    STEPS_FALLTHROUGH;
                case 0:
                    break;
                /// allow maximum 4 molecules of one species in reaction
                default:
                    AssertLog(0);
                }
            }
            for (const auto& specie: reaction.affectedSpecies) {
                if ((pSFlags[specie.globalIndex] & solver::Statedef::CLAMPED_POOLFLAG) != 0) {
                    continue;
                }
                dydx[specie.globalIndex] += specie.populationChange * numberOfReactionFirings;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void Wmrk4::_rk4(double pdt) {
    double dt_2 = pdt / 2.0;
    double dt_6 = pdt / 6.0;

    for (uint i = 0; i < pSpecs_tot; ++i) {
        yt[i] = pVals[i] + (dt_2 * pDyDx[i]);
    }
    _setderivs(yt, dyt);
    for (uint i = 0; i < pSpecs_tot; ++i) {
        yt[i] = pVals[i] + (dt_2 * dyt[i]);
    }
    _setderivs(yt, dym);
    for (uint i = 0; i < pSpecs_tot; ++i) {
        yt[i] = pVals[i] + (pdt * dym[i]);
        dym[i] += dyt[i];
    }
    _setderivs(yt, dyt);
    for (uint i = 0; i < pSpecs_tot; ++i) {
        pNewVals[i] = pVals[i] + dt_6 * (pDyDx[i] + dyt[i] + (2.0 * dym[i]));
    }
}

////////////////////////////////////////////////////////////////////////////////

void Wmrk4::_rksteps(double t1, double t2) {
    if (t1 == t2) {
        return;
    }
    AssertLog(t1 < t2);
    double t = t1;
    if (pDT <= 0.0) {
        std::ostringstream os;
        os << "dt is zero or negative. Call setDT() method.";
        ArgErrLog(os.str());
    }
    // Bugfix - step() calls this function with one dt step, so these two values
    // can be equal
    if (pDT > (t2 - t1)) {
        std::ostringstream os;
        os << "dt is larger than simulation step.";
        ArgErrLog(os.str());
    }

    /// step up until over maximum time
    while (t < t2) {
        if ((t + pDT) > t2) {
            break;
        }

        _setderivs(pVals, pDyDx);
        _rk4(pDT);
        _update();
        t += pDT;
    }

    ////////////////////////////////////////////////////////////////////////////
    // DEBUG: 25/03/09. This is general fix for this solver, inspired by a
    // problem with evaluating two supposedly equal double values as non-equal.
    // Now any discrepancy between simulation time and the end step time is dealt
    // with by solving to the simulation time by passing in the fraction (see
    // below) Changed _rk4() to take the dt as a double argument
    double tfrac = t2 - t;
    AssertLog(tfrac >= 0.0);

    // Lets only concern ourselves with fractions greater than 1 percent
    if (tfrac != 0.0)  // && tfrac/pDT >= 0.01)
    {
        AssertLog(tfrac < pDT);
        _setderivs(pVals, pDyDx);
        _rk4(tfrac);
        _update();
    }
    ////////////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////

void Wmrk4::_update() {
    /// update local values vector with computed counts
    for (uint i = 0; i < pSpecs_tot; ++i) {
        /// check clamped flag and only update if not clamped
        if ((pSFlags[i] & solver::Statedef::CLAMPED_POOLFLAG) != 0) {
            continue;
        } else {
            double newval = pNewVals[i];
            if (newval < 0.0) {
                newval = 0.0;
            }
            pVals[i] = newval;
        }
    }

    /// update pools with computed values
    uint Comps_N = statedef().countComps();
    uint Patches_N = statedef().countPatches();
    uint c_marker = 0;

    for (auto i: solver::comp_global_id::range(Comps_N)) {
        const auto compdef = statedef().compdef(i);
        const uint comp_Specs_N = compdef->countSpecs();
        for (auto j: solver::spec_local_id::range(comp_Specs_N)) {
            compdef->setCount(j, pVals[c_marker + j.get()]);
        }
        c_marker += comp_Specs_N;
    }

    for (auto i: solver::patch_global_id::range(Patches_N)) {
        const auto patchdef = statedef().patchdef(i);
        for (auto j: solver::spec_local_id::range(patchdef->countSpecs())) {
            patchdef->setCount(j, pVals[c_marker + j.get()]);
        }
        c_marker += patchdef->countSpecs();
    }
}

}  // namespace steps::wmrk4
