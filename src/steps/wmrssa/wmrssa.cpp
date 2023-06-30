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
#include "wmrssa.hpp"
#include "math/constants.hpp"
#include "model/reac.hpp"
#include "solver/reacdef.hpp"
#include "solver/types.hpp"
#include "sreac.hpp"

// logging
#include "util/error.hpp"
#include <easylogging++.h>

#include "util/checkpointing.hpp"

#define SCHEDULEWIDTH 32
#define MAXLEVELS     10

namespace steps::wmrssa {

void schedIDXSet_To_Vec(SchedIDXSet const& s, SchedIDXVec& v) {
    v.resize(s.size());
    std::copy(s.begin(), s.end(), v.begin());
}

////////////////////////////////////////////////////////////////////////////////

Wmrssa::Wmrssa(model::Model* m, wm::Geom* g, const rng::RNGptr& r)
    : API(m, g, r) {
    if (rng() == nullptr) {
        std::ostringstream os;
        os << "No RNG provided to solver initializer function";
        ArgErrLog(os.str());
    }

    for (auto const& c: statedef().comps()) {
        solver::comp_global_id compdef_gidx = c->gidx();
        uint comp_idx = _addComp(c);
        AssertLog(compdef_gidx.get() == comp_idx);
    }

    // Create the actual patches.
    for (auto const& p: statedef().patches()) {
        solver::patch_global_id patchdef_gidx = p->gidx();
        uint patch_idx = _addPatch(p);
        AssertLog(patchdef_gidx.get() == patch_idx);
    }

    _setup();
    // force update for zero order reactions
    _reset();
}

////////////////////////////////////////////////////////////////////////////////

Wmrssa::~Wmrssa() {
    for (auto& comp: pComps) {
        delete comp;
    }
    for (auto& patch: pPatches) {
        delete patch;
    }
}

///////////////////////////////////////////////////////////////////////////////

void Wmrssa::checkpoint(std::string const& file_name) {
    std::fstream cp_file;

    cp_file.open(file_name.c_str(), std::fstream::out | std::fstream::binary | std::fstream::trunc);

    API::checkpoint(cp_file);

    for (auto& comp: pComps) {
        comp->checkpoint(cp_file);
    }
    for (auto& patch: pPatches) {
        patch->checkpoint(cp_file);
    }

    statedef().checkpoint(cp_file);

    cp_file.close();
}

///////////////////////////////////////////////////////////////////////////////

void Wmrssa::restore(std::string const& file_name) {
    std::fstream cp_file;

    cp_file.open(file_name.c_str(), std::fstream::in | std::fstream::binary);

    cp_file.seekg(0);

    API::restore(cp_file);

    for (auto& comp: pComps) {
        comp->restore(cp_file);
    }

    for (auto& patch: pPatches) {
        patch->restore(cp_file);
    }

    statedef().restore(cp_file);

    cp_file.close();

    _reset();
}

////////////////////////////////////////////////////////////////////////////////

uint Wmrssa::_addComp(solver::Compdef* cdef) {
    auto* comp = new Comp(cdef);
    AssertLog(comp != nullptr);
    uint compidx = pComps.size();
    pComps.push_back(comp);
    pCompMap[cdef] = comp;
    return compidx;
}

////////////////////////////////////////////////////////////////////////////////

uint Wmrssa::_addPatch(solver::Patchdef* pdef) {
    Comp* icomp = nullptr;
    Comp* ocomp = nullptr;
    if (pdef->icompdef() != nullptr) {
        icomp = pCompMap[pdef->icompdef()];
    }
    if (pdef->ocompdef() != nullptr) {
        ocomp = pCompMap[pdef->ocompdef()];
    }
    auto* patch = new Patch(pdef, icomp, ocomp);
    AssertLog(patch != nullptr);
    uint patchidx = pPatches.size();
    pPatches.push_back(patch);
    return patchidx;
}

////////////////////////////////////////////////////////////////////////////////

void Wmrssa::_setup() {
    for (auto const& c: pComps) {
        c->setupKProcs(this);
    }

    for (auto const& p: pPatches) {
        p->setupKProcs(this);
    }

    // Resolve all dependencies
    for (auto const& c: pComps) {
        for (auto const& k: c->kprocs()) {
            k->setupDeps();  // This is no longer useful, right?
        }
        c->setupSpecDeps();
    }
    for (auto const& p: pPatches) {
        for (auto const& k: p->kprocs()) {
            k->setupDeps();
        }
        p->setupSpecDeps();
    }

    _build();
}

////////////////////////////////////////////////////////////////////////////////

std::string Wmrssa::getSolverName() const {
    return "wmrssa";
}

////////////////////////////////////////////////////////////////////////////////

std::string Wmrssa::getSolverDesc() const {
    return "Rejection-based SSA Method in well-mixed conditions, based on Thanh "
           "V, Zunino R, Priami C (n.d.) On the rejection-based algorithm for "
           "simulation and analysis of large-scale reaction networks. The "
           "Journal of Chemical Physics 142:244106";
}

////////////////////////////////////////////////////////////////////////////////

std::string Wmrssa::getSolverAuthors() const {
    return "Samuel Melchior";
}

////////////////////////////////////////////////////////////////////////////////

std::string Wmrssa::getSolverEmail() const {
    return "Please visit our website for more information "
           "(https://steps.sourceforge.net)";
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::reset() {
    uint comps = statedef().countComps();
    for (auto i: solver::comp_global_id::range(comps)) {
        statedef().compdef(i)->reset();
    }
    uint patches = statedef().countPatches();
    for (auto i: solver::patch_global_id::range(patches)) {
        statedef().patchdef(i)->reset();
    }

    for (auto const& comp: pComps) {
        comp->reset();
    }

    for (auto const& patch: pPatches) {
        patch->reset();
    }

    statedef().resetTime();
    statedef().resetNSteps();

    _reset();
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::run(double endtime) {
    if (endtime < statedef().time()) {
        std::ostringstream os;
        os << "Endtime is before current simulation time";
        ArgErrLog(os.str());
    }
    while (statedef().time() < endtime) {
        if (pA0 == 0.0) {
            break;
        }
        bool isRejected = true;
        double erlangFactor = 1;
        KProc* kp;
        while (isRejected) {
            uint cur_node = _getNext();
            kp = pKProcs[cur_node];
            if (kp == nullptr) {
                break;
            }
            double randnum = rng()->getUnfIE() * pLevels[0][cur_node];
            if (randnum <= kp->propensityLB() || randnum <= kp->rate()) {
                isRejected = false;
            }
            erlangFactor *= rng()->getUnfIE();
        }
        double dt = -1 / pA0 * log(erlangFactor);
        if ((statedef().time() + dt) > endtime) {
            break;
        }
        _executeStep(kp, dt);
    }
    statedef().setTime(endtime);
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::advance(double adv) {
    if (adv < 0.0) {
        std::ostringstream os;
        os << "Time to advance cannot be negative";
        ArgErrLog(os.str());
    }

    double endtime = statedef().time() + adv;
    run(endtime);
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::step() {
    bool isRejected = true;
    double erlangFactor = 1;
    KProc* kp = nullptr;
    while (isRejected) {
        if (pA0 == 0.0) {
            break;
        }
        uint cur_node = _getNext();
        kp = pKProcs[cur_node];
        if (kp == nullptr) {
            break;
        }
        double randnum = rng()->getUnfIE() * pLevels[0][cur_node];
        if (randnum <= kp->propensityLB() || randnum <= kp->rate()) {
            isRejected = false;
        }
        erlangFactor *= rng()->getUnfIE();
    }
    AssertLog(kp != nullptr);
    AssertLog(pA0 != 0.0);
    double dt = -1 / pA0 * std::log(erlangFactor);
    _executeStep(kp, dt);
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::getTime() const {
    return statedef().time();
}

////////////////////////////////////////////////////////////////////////

uint Wmrssa::getNSteps() const {
    return statedef().nsteps();
}

////////////////////////////////////////////////////////////////////////////////

void Wmrssa::setTime(double time) {
    statedef().setTime(time);
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::setNSteps(uint nsteps) {
    statedef().setNSteps(nsteps);
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::_getCompVol(solver::comp_global_id cidx) const {
    AssertLog(cidx < statedef().countComps());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    return comp->vol();
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::_setCompVol(solver::comp_global_id cidx, double vol) {
    AssertLog(cidx < statedef().countComps());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    comp->setVol(vol);

    // Reset the reaction C constants
    Comp* lcomp = pComps[cidx.get()];
    assert(lcomp->def() == comp);

    for (auto& ck: lcomp->kprocs()) {
        ck->resetCcst();
    }

    for (auto& p: lcomp->ipatches()) {
        for (auto& pk: p->kprocs()) {
            pk->resetCcst();
        }
    }
    for (auto& p: lcomp->opatches()) {
        for (auto& pk: p->kprocs()) {
            pk->resetCcst();
        }
    }
    // Rates have changed
    _reset();
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::_getCompSpecCount(solver::comp_global_id cidx, solver::spec_global_id sidx) const {
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

////////////////////////////////////////////////////////////////////////

void Wmrssa::_setCompSpecCount(solver::comp_global_id cidx, solver::spec_global_id sidx, double n) {
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
    if (n > std::numeric_limits<unsigned int>::max()) {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << std::numeric_limits<unsigned int>::max() << ").\n";
        ArgErrLog(os.str());
    }

    assert(n >= 0.0);
    double n_int = std::floor(n);
    double n_frc = n - n_int;
    uint c = static_cast<uint>(n_int);
    if (n_frc > 0.0) {
        double rand01 = rng()->getUnfIE();
        if (rand01 < n_frc) {
            c++;
        }
    }
    auto n_double = static_cast<double>(c);
    comp->setCount(slidx, n_double);
    // Rates have changed
    _reset();
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::_getCompSpecAmount(solver::comp_global_id cidx, solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getCompSpecCount(cidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::_setCompSpecAmount(solver::comp_global_id cidx,
                                solver::spec_global_id sidx,
                                double a) {
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setCompSpecCount(cidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::_getCompSpecConc(solver::comp_global_id cidx, solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getCompSpecCount(cidx, sidx);
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    double vol = comp->vol();
    return count / (1.0e3 * vol * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::_setCompSpecConc(solver::comp_global_id cidx, solver::spec_global_id sidx, double c) {
    AssertLog(c >= 0.0);
    assert(cidx < statedef().countComps());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    double count = c * (1.0e3 * comp->vol() * math::AVOGADRO);
    // the following method does all the necessary argument checking
    _setCompSpecCount(cidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////

bool Wmrssa::_getCompSpecClamped(solver::comp_global_id cidx, solver::spec_global_id sidx) const {
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

////////////////////////////////////////////////////////////////////////

void Wmrssa::_setCompSpecClamped(solver::comp_global_id cidx, solver::spec_global_id sidx, bool b) {
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
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::_getCompReacK(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
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

////////////////////////////////////////////////////////////////////////

void Wmrssa::_setCompReacK(solver::comp_global_id cidx, solver::reac_global_id ridx, double kf) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    assert(kf >= 0.0);
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = comp->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    comp->setKcst(lridx, kf);

    // Reset the reaction C constants
    Comp* lcomp = pComps[cidx.get()];
    assert(lcomp->def() == comp);

    wmrssa::KProc* lreac = lcomp->reac(lridx);
    assert(lreac->defr() == comp->reacdef(lridx));
    lreac->resetCcst();

    // Rates have changed
    _reset();
}

////////////////////////////////////////////////////////////////////////

bool Wmrssa::_getCompReacActive(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
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

////////////////////////////////////////////////////////////////////////

void Wmrssa::_setCompReacActive(solver::comp_global_id cidx, solver::reac_global_id ridx, bool a) {
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

    // It's cheaper to just recompute everything.
    _reset();
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::_getCompReacC(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
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

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx.get()];
    assert(lcomp->def() == comp);
    // Reacdef local indices in Compdef object also have same index
    // as Reacs in Comp object
    KProc* lreac = lcomp->reac(lridx);
    assert(lreac->defr() == comp->reacdef(lridx));

    return lreac->c();
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::_getCompReacH(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
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

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx.get()];
    assert(lcomp->def() == comp);
    // Reacdef local indices in Compdef object also have same index
    // as Reacs in Comp object
    KProc* lreac = lcomp->reac(lridx);
    assert(lreac->defr() == comp->reacdef(lridx));

    return lreac->h();
}

////////////////////////////////////////////////////////////////////////

long double Wmrssa::_getCompReacA(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
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

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx.get()];
    assert(lcomp->def() == comp);
    // Reacdef local indices in Compdef object also have same index
    // as Reacs in Comp object
    KProc* lreac = lcomp->reac(lridx);
    assert(lreac->defr() == comp->reacdef(lridx));

    return static_cast<long double>(lreac->rate());
}

////////////////////////////////////////////////////////////////////////

unsigned long long Wmrssa::_getCompReacExtent(solver::comp_global_id cidx,
                                              solver::reac_global_id ridx) const {
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

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx.get()];
    assert(lcomp->def() == comp);
    // Reacdef local indices in Compdef object also have same index
    // as Reacs in Comp object
    KProc* lreac = lcomp->reac(lridx);
    assert(lreac->defr() == comp->reacdef(lridx));

    return lreac->getExtent();
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::_resetCompReacExtent(solver::comp_global_id cidx, solver::reac_global_id ridx) {
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

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx.get()];
    assert(lcomp->def() == comp);
    // Reacdef local indices in Compdef object also have same index
    // as Reacs in Comp object
    KProc* lreac = lcomp->reac(lridx);
    assert(lreac->defr() == comp->reacdef(lridx));

    lreac->resetExtent();
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::_getPatchArea(solver::patch_global_id pidx) const {
    AssertLog(pidx < statedef().countPatches());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    return patch->area();
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::_setPatchArea(solver::patch_global_id pidx, double area) {
    AssertLog(pidx < statedef().countPatches());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    patch->setArea(area);
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::_getPatchSpecCount(solver::patch_global_id pidx, solver::spec_global_id sidx) const {
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

////////////////////////////////////////////////////////////////////////

void Wmrssa::_setPatchSpecCount(solver::patch_global_id pidx,
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
    if (n > std::numeric_limits<unsigned int>::max()) {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << std::numeric_limits<unsigned int>::max() << ").\n";
        ArgErrLog(os.str());
    }

    double n_int = std::floor(n);
    double n_frc = n - n_int;
    uint c = static_cast<uint>(n_int);
    if (n_frc > 0.0) {
        double rand01 = statedef().rng()->getUnfIE();
        if (rand01 < n_frc) {
            c++;
        }
    }
    n_int = static_cast<double>(c);
    patch->setCount(slidx, n_int);
    _reset();
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::_getPatchSpecAmount(solver::patch_global_id pidx,
                                   solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getPatchSpecCount(pidx, sidx);
    return (count / math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::_setPatchSpecAmount(solver::patch_global_id pidx,
                                 solver::spec_global_id sidx,
                                 double a) {
    AssertLog(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setPatchSpecCount(pidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////

bool Wmrssa::_getPatchSpecClamped(solver::patch_global_id pidx, solver::spec_global_id sidx) const {
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

////////////////////////////////////////////////////////////////////////

void Wmrssa::_setPatchSpecClamped(solver::patch_global_id pidx,
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
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::_getPatchSReacK(solver::patch_global_id pidx, solver::sreac_global_id ridx) const {
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

////////////////////////////////////////////////////////////////////////

void Wmrssa::_setPatchSReacK(solver::patch_global_id pidx,
                             solver::sreac_global_id ridx,
                             double kf) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    assert(kf >= 0.0);
    solver::Patchdef* patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lridx = patch->sreacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    patch->setKcst(lridx, kf);

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx.get()];
    AssertLog(lpatch->def() == patch);
    // SReacdef local indices in Patchdef object also have same index
    //  as SReacs in Patch object
    KProc* lsreac = lpatch->sreac(lridx);
    assert(lsreac->defsr() == patch->sreacdef(lridx));
    lsreac->resetCcst();

    // Rates have changed
    _reset();
}

////////////////////////////////////////////////////////////////////////

bool Wmrssa::_getPatchSReacActive(solver::patch_global_id pidx,
                                  solver::sreac_global_id ridx) const {
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

////////////////////////////////////////////////////////////////////////

void Wmrssa::_setPatchSReacActive(solver::patch_global_id pidx,
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

    // It's cheaper to just recompute everything
    _reset();
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::_getPatchSReacC(solver::patch_global_id pidx, solver::sreac_global_id ridx) const {
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

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx.get()];
    AssertLog(lpatch->def() == patch);
    // SReacdef local indices in Patchdef object also have same index
    //  as SReacs in Patch object
    KProc* lsreac = lpatch->sreac(lridx);
    assert(lsreac->defsr() == patch->sreacdef(lridx));

    return lsreac->c();
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::_getPatchSReacH(solver::patch_global_id pidx, solver::sreac_global_id ridx) const {
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

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx.get()];
    AssertLog(lpatch->def() == patch);
    // SReacdef local indices in Patchdef object also have same index
    //  as SReacs in Patch object
    KProc* lsreac = lpatch->sreac(lridx);
    assert(lsreac->defsr() == patch->sreacdef(lridx));

    return lsreac->h();
}

////////////////////////////////////////////////////////////////////////

double Wmrssa::_getPatchSReacA(solver::patch_global_id pidx, solver::sreac_global_id ridx) const {
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

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx.get()];
    AssertLog(lpatch->def() == patch);
    // SReacdef local indices in Patchdef object also have same index
    //  as SReacs in Patch object
    KProc* lsreac = lpatch->sreac(lridx);
    assert(lsreac->defsr() == patch->sreacdef(lridx));

    return lsreac->rate();
}

////////////////////////////////////////////////////////////////////////

unsigned long long Wmrssa::_getPatchSReacExtent(solver::patch_global_id pidx,
                                                solver::sreac_global_id ridx) const {
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

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx.get()];
    AssertLog(lpatch->def() == patch);
    // SReacdef local indices in Patchdef object also have same index
    //  as SReacs in Patch object
    KProc* lsreac = lpatch->sreac(lridx);
    assert(lsreac->defsr() == patch->sreacdef(lridx));

    return lsreac->getExtent();
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::_resetPatchSReacExtent(solver::patch_global_id pidx, solver::sreac_global_id ridx) {
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

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx.get()];
    AssertLog(lpatch->def() == patch);
    // SReacdef local indices in Patchdef object also have same index
    //  as SReacs in Patch object
    KProc* lsreac = lpatch->sreac(lridx);
    assert(lsreac->defsr() == patch->sreacdef(lridx));

    lsreac->resetExtent();
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::addKProc(KProc* kp) {
    AssertLog(kp != nullptr);

    SchedIDX nidx(static_cast<uint>(pKProcs.size()));  // because pKProcs.size() is ulong
    pKProcs.push_back(kp);
    kp->setSchedIDX(nidx);
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::_build() {
    AssertLog(pBuilt == false);

    // Setup level.
    uint clsize = pKProcs.size();
    if (clsize == 0) {
        return;
    }

    // Work up.
    do {
        // Make sure the new size is a multiple of SCHEDULEWIDTH.
        uint extra = clsize % SCHEDULEWIDTH;
        if (extra != 0) {
            clsize += SCHEDULEWIDTH - extra;
        }

        // Create the level and add it.
        pLevels.emplace_back(clsize);

        // Prepare for next level.
        clsize = clsize / SCHEDULEWIDTH;
    } while (clsize > 1);

    // Set top level.
    pA0 = 0.0;

    // Time to create ONE indices table to hold the run's present reaction of
    // choice's update vector. This will be re-used and replace old-version's
    // hard-coded table in _update. Size is the maximum possible, found by looping
    // over all KProcs. This little bit of computational time is well worth all
    // that needless memory allocation
    uint maxupvecsize = 0;
    auto kproc_end = pKProcs.end();
    for (auto kproc = pKProcs.begin(); kproc != kproc_end; ++kproc) {
        if ((*kproc)->updVecSize() > maxupvecsize) {
            maxupvecsize = (*kproc)->updVecSize();
        }
    }

    pMaxUpSize = maxupvecsize;
    pIndices.resize(pMaxUpSize);

    // Also let's create a random number holder-table,
    // size of number of KProcs % SCHEDULEWIDTH or pLevels.size()
    // This will be re-used in _getNext as opposed to hard-coded (again maximum
    // limit).
    uint lsize = pLevels.size();
    pRannum.resize(lsize);

    pBuilt = true;
}

////////////////////////////////////////////////////////////////////////

uint Wmrssa::_getNext() {
    AssertLog(pA0 >= 0.0);
    // Quick check to see whether nothing is there.
    if (pA0 == 0.0) {
        return 0;
    }

    // Start at top level.
    uint clevel = pLevels.size();
    // And start at the first node of that level.
    uint cur_node = 0;

    // Prepare random numbers.
    for (uint i = 0; i < clevel; ++i) {
        pRannum[i] = rng()->getUnfIE();
    }

    // Run until top level.
    double a0 = pA0;
    while (clevel != 0) {
        // Decrease the current level.
        clevel--;
        // and start looking in the right place.
        cur_node *= SCHEDULEWIDTH;
        uint max_node = cur_node + SCHEDULEWIDTH;

        // Fetch the level.
        const auto& level = pLevels[clevel];

        // Compute local selector.
        double selector = pRannum[clevel] * a0;

        // Compare.
        double accum = 0.0;
        // 27/10/09 I.H. 'old' removed from for loop because not used.
        // double old = 0.0;
        double curval = 0.0;
        for (uint i = 0; i < SCHEDULEWIDTH; ++i) {
            curval = level[cur_node];
            if (selector < curval + accum) {
                break;
            }
            accum += curval;
            // old = accum;
            cur_node++;
        }

        // Checks.
        AssertLog(cur_node < max_node);
        AssertLog(curval > 0.0);
        a0 = curval;
    }

    // Check.
    AssertLog(cur_node < pKProcs.size());
    return cur_node;
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::_reset() {
    if (pKProcs.empty()) {
        return;
    }

    for (auto c: pComps) {
        const auto& cnt_vec = c->def()->pools();
        for (auto pool: cnt_vec.range()) {
            c->setBounds(pool, static_cast<int>(cnt_vec[pool]));
        }
    }
    for (const auto& patch: pPatches) {
        const auto& cnt_vec = patch->def()->pools();
        for (auto pool: cnt_vec.range()) {
            patch->setBounds(pool, static_cast<int>(cnt_vec[pool]));
        }
    }

    // Reset the basic level: compute rates.
    auto& oldlevel = pLevels[0];
    uint cur_node = 0;
    for (auto const& kp: pKProcs) {
        oldlevel[cur_node++] = kp->rate(BOUNDS);
    }

    // Work up.
    for (uint cur_level = 1; cur_level < pLevels.size(); ++cur_level) {
        // Compute the number of nodes to reset on this level.
        uint numnodes = pLevels[cur_level - 1].size() / SCHEDULEWIDTH;

        // Fetch a pointer to this level.
        auto& level = pLevels[cur_level];

        // Recompute them.
        uint child_node = 0;
        for (cur_node = 0; cur_node < numnodes; ++cur_node) {
            double val = 0.0;
            for (uint i = 0; i < SCHEDULEWIDTH; ++i) {
                val += oldlevel[child_node++];
            }
            level[cur_node] = val;
        }

        // Copy the level.
        oldlevel = level;
    }

    // Compute zero propensity.
    pA0 = 0.0;
    for (uint i = 0; i < SCHEDULEWIDTH; ++i) {
        pA0 += oldlevel[i];
    }
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::_update(SchedIDXVec const& entries) {
    if (countKProcs() == 0) {
        return;
    }

    // Prefetch zero level.
    auto& level0 = pLevels[0];
    // Number of entries.
    AssertLog(entries.size() <= pMaxUpSize);  /////////

    // Recompute rates.
    uint prev_e = 0xFFFFFFFF;
    uint cur_e = 0;
    for (auto const& sidx: entries) {
        // Fetch index.
        uint idx = sidx.get();
        // Recompute rate, get difference, and store.
        double newrate = pKProcs[idx]->rate(BOUNDS);
        level0[idx] = newrate;

        // Store and collapse if possible.
        idx /= SCHEDULEWIDTH;
        if (prev_e == 0xFFFFFFFF) {
            prev_e = 0;
            pIndices[cur_e++] = idx;
        } else if (pIndices[prev_e] != idx) {
            prev_e = cur_e;
            pIndices[cur_e++] = idx;
        }
    }
    uint nentries = cur_e;

    // Update upper levels.
    uint nlevels = pLevels.size();
    auto& prevlevel = pLevels[0];
    for (uint l = 1; l < nlevels; ++l) {
        // Update the first entry.
        cur_e = 0;
        prev_e = 0xFFFFFFFF;

        // Fetch a pointer to the current level.
        auto& currlevel = pLevels[l];

        // Recompute the entries.
        for (uint e = 0; e < nentries; ++e) {
            // Fetch index.
            uint idx = pIndices[e];

            // Recompute.
            double val = 0.0;
            uint idx2 = idx * SCHEDULEWIDTH;
            for (uint i = 0; i < SCHEDULEWIDTH; ++i) {
                val += prevlevel[idx2++];
            }
            currlevel[idx] = val;

            // Store and collapse if possible.
            idx /= SCHEDULEWIDTH;
            if (prev_e == 0xFFFFFFFF) {
                prev_e = 0;
                pIndices[cur_e++] = idx;
            } else if (pIndices[prev_e] != idx) {
                prev_e = cur_e;
                pIndices[cur_e++] = idx;
            }
        }

        // Update the pointer to the previous level.
        prevlevel = currlevel;

        // cur_e now is the new number of entries to handle.
        nentries = cur_e;
    }

    // Update zero propensity.
    auto& toplevel = pLevels[pLevels.size() - 1];
    pA0 = 0.0;
    for (uint i = 0; i < SCHEDULEWIDTH; ++i) {
        pA0 += toplevel[i];
    }
}

////////////////////////////////////////////////////////////////////////

void Wmrssa::_executeStep(KProc* kp, double dt) {
    SchedIDXVec const& upd = kp->apply();
    if (upd.size() > 0) {
        _update(upd);
    }
    statedef().incTime(dt);
    statedef().incNSteps(1);
}

}  // namespace steps::wmrssa
