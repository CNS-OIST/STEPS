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
#include "wmdirect.hpp"

#include "math/constants.hpp"
#include "model/reac.hpp"
#include "rng/rng.hpp"
#include "solver/reacdef.hpp"
#include "sreac.hpp"
// logging
#include "util/error.hpp"

#include "util/checkpointing.hpp"

#define SCHEDULEWIDTH 32
#define MAXLEVELS     10

////////////////////////////////////////////////////////////////////////////////

namespace steps::wmdirect {


////////////////////////////////////////////////////////////////////////////////

/// Unary function that calls the array delete[] operator on pointers. Easy
/// to use with STL/Boost (see steps::tools::DeletePointer).
///
struct DeleteArray {
    template <typename Type>
    void operator()(Type* pointer) const {
        delete[] pointer;
    }
};

////////////////////////////////////////////////////////////////////////////////

void schedIDXSet_To_Vec(SchedIDXSet const& s, SchedIDXVec& v) {
    v.resize(s.size());
    std::copy(s.begin(), s.end(), v.begin());
}

////////////////////////////////////////////////////////////////////////////////

Wmdirect::Wmdirect(model::Model* m, wm::Geom* g, const rng::RNGptr& r)
    : API(*m, *g, r) {
    if (rng() == nullptr) {
        std::ostringstream os;
        os << "No RNG provided to solver initializer function";
        ArgErrLog(os.str());
    }

    for (auto const& c: statedef().comps()) {
        solver::comp_global_id compdef_gidx = c->gidx();
        uint comp_idx = _addComp(c.get());
        AssertLog(compdef_gidx == comp_idx);
    }

    // Create the actual patches.
    for (auto const& p: statedef().patches()) {
        solver::patch_global_id patchdef_gidx = p->gidx();
        uint patch_idx = _addPatch(p.get());
        AssertLog(patchdef_gidx == patch_idx);
    }

    _setup();

    // force update for zero order reactions
    _reset();
}

////////////////////////////////////////////////////////////////////////////////

Wmdirect::~Wmdirect() {
    for (auto const& c: pComps) {
        delete c;
    }
    for (auto const& p: pPatches) {
        delete p;
    }
}

///////////////////////////////////////////////////////////////////////////////

void Wmdirect::checkpoint(std::string const& file_name) {
    std::fstream cp_file;

    cp_file.open(file_name.c_str(), std::fstream::out | std::fstream::binary | std::fstream::trunc);

    API::checkpoint(cp_file);

    for (auto const& c: pComps) {
        c->checkpoint(cp_file);
    }
    for (auto const& p: pPatches) {
        p->checkpoint(cp_file);
    }

    statedef().checkpoint(cp_file);

    cp_file.close();
}

///////////////////////////////////////////////////////////////////////////////

void Wmdirect::restore(std::string const& file_name) {
    std::fstream cp_file;

    cp_file.open(file_name.c_str(), std::fstream::in | std::fstream::binary);

    cp_file.seekg(0);

    API::restore(cp_file);

    for (auto const& c: pComps) {
        c->restore(cp_file);
    }
    for (auto const& p: pPatches) {
        p->restore(cp_file);
    }

    statedef().restore(cp_file);

    if (cp_file.fail()) {
        ArgErrLog("Checkpoint restoration failed.");
    }

    cp_file.close();

    _reset();
}

////////////////////////////////////////////////////////////////////////////////

uint Wmdirect::_addComp(solver::Compdef* cdef) {
    auto* comp = new Comp(cdef, this);
    AssertLog(comp != nullptr);
    auto compidx = pComps.size();
    pComps.push_back(comp);
    pCompMap[cdef] = comp;
    return compidx;
}

////////////////////////////////////////////////////////////////////////////////

uint Wmdirect::_addPatch(solver::Patchdef* pdef) {
    Comp* icomp = nullptr;
    Comp* ocomp = nullptr;
    if (pdef->icompdef() != nullptr) {
        icomp = pCompMap[pdef->icompdef()];
    }
    if (pdef->ocompdef() != nullptr) {
        ocomp = pCompMap[pdef->ocompdef()];
    }
    auto* patch = new Patch(pdef, icomp, ocomp, this);
    auto patchidx = pPatches.size();
    pPatches.push_back(patch);
    return patchidx;
}

////////////////////////////////////////////////////////////////////////////////

void Wmdirect::_setup() {
    for (auto const& c: pComps) {
        c->setupKProcs(this);
    }

    for (auto const& p: pPatches) {
        p->setupKProcs(this);
    }

    // Resolve all dependencies
    for (auto const& c: pComps) {
        for (auto const& k: c->kprocs()) {
            k->setupDeps();
        }
    }
    for (auto const& p: pPatches) {
        for (auto const& k: p->kprocs()) {
            k->setupDeps();
        }
    }

    _build();
}

////////////////////////////////////////////////////////////////////////////////

std::string Wmdirect::getSolverName() const {
    return "wmdirect";
}

////////////////////////////////////////////////////////////////////////////////

std::string Wmdirect::getSolverDesc() const {
    return "SSA Direct Method in well-mixed conditions";
}

////////////////////////////////////////////////////////////////////////////////

std::string Wmdirect::getSolverAuthors() const {
    return "Stefan Wils and Iain Hepburn";
}

////////////////////////////////////////////////////////////////////////////////

std::string Wmdirect::getSolverEmail() const {
    return "stefan@tnb.ua.ac.be, ihepburn@oist.jp";
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::reset() {
    for (auto& comp: statedef().comps()) {
        comp->reset();
    }
    for (auto& patch: statedef().patches()) {
        patch->reset();
    }

    for (auto const& comp: pComps) {
        comp->reset();
    }
    for (auto const& patch: pPatches) {
        patch->reset();
    }

    statedef().resetTime();
    statedef().resetNSteps();

    for (auto const& kp: pKProcs) {
        kp->reset();
    }

    _reset();
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::run(double endtime) {
    if (endtime < statedef().time()) {
        std::ostringstream os;
        os << "Endtime is before current simulation time";
        ArgErrLog(os.str());
    }
    while (statedef().time() < endtime) {
        KProc* kp = _getNext();
        if (kp == nullptr) {
            break;
        }
        double a0 = getA0();
        if (a0 == 0.0) {
            break;
        }
        double dt = rng()->getExp(a0);
        if ((statedef().time() + dt) > endtime) {
            break;
        }
        _executeStep(kp, dt);
    }
    statedef().setTime(endtime);
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::advance(double adv) {
    if (adv < 0.0) {
        std::ostringstream os;
        os << "Time to advance cannot be negative";
        ArgErrLog(os.str());
    }

    double endtime = statedef().time() + adv;
    run(endtime);
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::step() {
    KProc* kp = _getNext();
    if (kp == nullptr) {
        return;
    }
    double a0 = getA0();
    if (a0 == 0.0) {
        return;
    }
    double dt = rng()->getExp(a0);
    _executeStep(kp, dt);
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::getTime() const {
    return statedef().time();
}

////////////////////////////////////////////////////////////////////////

uint Wmdirect::getNSteps() const {
    return statedef().nsteps();
}

////////////////////////////////////////////////////////////////////////////////

void Wmdirect::setTime(double time) {
    statedef().setTime(time);
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::setNSteps(uint nsteps) {
    statedef().setNSteps(nsteps);
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getCompVol(solver::comp_global_id cidx) const {
    AssertLog(cidx < statedef().countComps());
    return statedef().compdef(cidx).vol();
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setCompVol(solver::comp_global_id cidx, double vol) {
    AssertLog(cidx < statedef().countComps());
    auto& comp = statedef().compdef(cidx);
    comp.setVol(vol);

    // Reset the reaction C constants
    Comp* lcomp = pComps[cidx.get()];
    AssertLog(lcomp->def() == &comp);

    for (auto const& ck: lcomp->kprocs()) {
        ck->resetCcst();
    }

    for (auto const& patch: pPatches) {
        patch->reset();
    }
    for (auto& p: lcomp->ipatches()) {
        for (auto const& pk: p->kprocs()) {
            pk->resetCcst();
        }
    }
    for (auto& p: lcomp->opatches()) {
        for (auto const& pk: p->kprocs()) {
            pk->resetCcst();
        }
    }
    // Rates have changed
    _reset();
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getCompSpecCount(solver::comp_global_id cidx, solver::spec_global_id sidx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    const auto& comp = statedef().compdef(cidx);
    solver::spec_local_id slidx = comp.specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }
    return comp.pools()[slidx];
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setCompSpecCount(solver::comp_global_id cidx,
                                 solver::spec_global_id sidx,
                                 double n) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    auto& comp = statedef().compdef(cidx);
    solver::spec_local_id slidx = comp.specG2L(sidx);
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

    AssertLog(n >= 0.0);
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
    comp.setCount(slidx, n_double);
    // Rates have changed
    _reset();
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getCompSpecAmount(solver::comp_global_id cidx,
                                    solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getCompSpecCount(cidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setCompSpecAmount(solver::comp_global_id cidx,
                                  solver::spec_global_id sidx,
                                  double a) {
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setCompSpecCount(cidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getCompSpecConc(solver::comp_global_id cidx, solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getCompSpecCount(cidx, sidx);
    const auto& comp = statedef().compdef(cidx);
    double vol = comp.vol();
    return count / (1.0e3 * vol * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setCompSpecConc(solver::comp_global_id cidx,
                                solver::spec_global_id sidx,
                                double c) {
    AssertLog(c >= 0.0);
    AssertLog(cidx < statedef().countComps());
    auto& comp = statedef().compdef(cidx);
    double count = c * (1.0e3 * comp.vol() * math::AVOGADRO);
    // the following method does all the necessary argument checking
    _setCompSpecCount(cidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////

bool Wmdirect::_getCompSpecClamped(solver::comp_global_id cidx, solver::spec_global_id sidx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    const auto& comp = statedef().compdef(cidx);
    solver::spec_local_id lsidx = comp.specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    return comp.clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setCompSpecClamped(solver::comp_global_id cidx,
                                   solver::spec_global_id sidx,
                                   bool b) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    auto& comp = statedef().compdef(cidx);
    solver::spec_local_id lsidx = comp.specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    comp.setClamped(lsidx, b);
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getCompComplexCount(
    solver::comp_global_id cidx,
    solver::complex_global_id cmplIdx,
    const std::vector<
        util::strongid_vector<solver::complex_substate_id, model::SubunitStateFilter>>& f) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(cmplIdx < statedef().countComplexes());
    solver::Compdef& comp = statedef().compdef(cidx);

    auto filt = comp.GetFilter(cmplIdx, f);
    const auto& states = comp.complexStates(cmplIdx);
    filt->processUpdates(states);
    return filt->nbMatches(states);
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setCompComplexCount(
    solver::comp_global_id cidx,
    solver::complex_global_id cmplIdx,
    const util::strongid_vector<solver::complex_substate_id, uint>& i,
    double n) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(cmplIdx < statedef().countComplexes());
    solver::Compdef& comp = statedef().compdef(cidx);

    if (n > std::numeric_limits<unsigned int>::max()) {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << std::numeric_limits<unsigned int>::max() << ").\n";
        ArgErrLog(os.str());
    }

    AssertLog(n >= 0.0);
    double n_int = std::floor(n);
    double n_frc = n - n_int;
    uint c = static_cast<uint>(n_int);
    if (n_frc > 0.0) {
        double rand01 = rng()->getUnfIE();
        if (rand01 < n_frc) {
            c++;
        }
    }

    std::vector<solver::complex_individual_id> existingInds;
    for (auto state: comp.complexStates(cmplIdx)) {
        if (state.second == i) {
            existingInds.push_back(state.first);
        }
    }
    // If we need to remove some of the complexes
    if (existingInds.size() > c) {
        uint nbToDel = existingInds.size() - c;
        for (auto it = existingInds.begin(); it != existingInds.begin() + nbToDel; ++it) {
            comp.removeComplex(cmplIdx, *it);
        }
    } else {
        for (uint nb = 0; nb < c - existingInds.size(); ++nb) {
            comp.addComplex(cmplIdx, getNextComplexInd(), i);
        }
    }

    // Rates have changed
    _reset();
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getCompComplexAmount(
    solver::comp_global_id cidx,
    solver::complex_global_id sidx,
    const std::vector<
        util::strongid_vector<solver::complex_substate_id, model::SubunitStateFilter>>& f) const {
    // the following method does all the necessary argument checking
    double count = _getCompComplexCount(cidx, sidx, f);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setCompComplexAmount(
    solver::comp_global_id cidx,
    solver::complex_global_id sidx,
    const util::strongid_vector<solver::complex_substate_id, uint>& i,
    double a) {
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setCompComplexCount(cidx, sidx, i, a2);
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getCompComplexConc(
    solver::comp_global_id cidx,
    solver::complex_global_id sidx,
    const std::vector<
        util::strongid_vector<solver::complex_substate_id, model::SubunitStateFilter>>& f) const {
    // the following method does all the necessary argument checking
    double count = _getCompComplexCount(cidx, sidx, f);
    solver::Compdef& comp = statedef().compdef(cidx);
    double vol = comp.vol();
    return count / (1.0e3 * vol * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setCompComplexConc(
    solver::comp_global_id cidx,
    solver::complex_global_id sidx,
    const util::strongid_vector<solver::complex_substate_id, uint>& i,
    double c) {
    AssertLog(c >= 0.0);
    AssertLog(cidx < statedef().countComps());
    solver::Compdef& comp = statedef().compdef(cidx);
    double count = c * (1.0e3 * comp.vol() * math::AVOGADRO);
    // the following method does all the necessary argument checking
    _setCompComplexCount(cidx, sidx, i, count);
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getCompComplexSUSCount(
    solver::comp_global_id cidx,
    solver::complex_global_id sidx,
    const std::vector<
        util::strongid_vector<solver::complex_substate_id, model::SubunitStateFilter>>& f,
    solver::complex_substate_id m) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countComplexes());
    solver::Compdef& comp = statedef().compdef(cidx);

    auto filt = comp.GetFilter(sidx, f);
    filt->processUpdates(comp.complexStates(sidx));
    return filt->nbSus(m, comp.complexStates(sidx));
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getCompComplexSUSConc(
    solver::comp_global_id cidx,
    solver::complex_global_id sidx,
    const std::vector<
        util::strongid_vector<solver::complex_substate_id, model::SubunitStateFilter>>& f,
    solver::complex_substate_id m) const {
    double count = _getCompComplexSUSCount(cidx, sidx, f, m);
    solver::Compdef& comp = statedef().compdef(cidx);
    double vol = comp.vol();
    return count / (1.0e3 * vol * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getCompComplexSUSAmount(
    solver::comp_global_id cidx,
    solver::complex_global_id sidx,
    const std::vector<
        util::strongid_vector<solver::complex_substate_id, model::SubunitStateFilter>>& f,
    solver::complex_substate_id m) const {
    double count = _getCompComplexSUSCount(cidx, sidx, f, m);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getCompReacK(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    const auto& comp = statedef().compdef(cidx);
    solver::reac_local_id lridx = comp.reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    return comp.kcst(lridx);
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setCompReacK(solver::comp_global_id cidx, solver::reac_global_id ridx, double kf) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    AssertLog(kf >= 0.0);
    auto& comp = statedef().compdef(cidx);
    solver::reac_local_id lridx = comp.reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    comp.setKcst(lridx, kf);

    // Reset the reaction C constants
    Comp* lcomp = pComps[cidx.get()];
    AssertLog(lcomp->def() == &comp);

    steps::wmdirect::KProc* lreac = lcomp->reac(lridx);
    AssertLog(lreac->defr() == &comp.reacdef(lridx));
    lreac->resetCcst();

    // Rates have changed
    _reset();
}

////////////////////////////////////////////////////////////////////////

bool Wmdirect::_getCompReacActive(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    const auto& comp = statedef().compdef(cidx);
    solver::reac_local_id lridx = comp.reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    return comp.active(lridx);
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setCompReacActive(solver::comp_global_id cidx,
                                  solver::reac_global_id ridx,
                                  bool a) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    auto& comp = statedef().compdef(cidx);
    solver::reac_local_id lridx = comp.reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    comp.setActive(lridx, a);

    // It's cheaper to just recompute everything.
    _reset();
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getCompReacC(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    const auto& comp = statedef().compdef(cidx);
    solver::reac_local_id lridx = comp.reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx.get()];
    AssertLog(lcomp->def() == &comp);
    // Reacdef local indices in Compdef object also have same index
    // as Reacs in Comp object
    KProc* lreac = lcomp->reac(lridx);
    AssertLog(lreac->defr() == &comp.reacdef(lridx));

    return lreac->c();
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getCompReacH(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    const auto& comp = statedef().compdef(cidx);
    solver::reac_local_id lridx = comp.reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx.get()];
    AssertLog(lcomp->def() == &comp);
    // Reacdef local indices in Compdef object also have same index
    // as Reacs in Comp object
    KProc* lreac = lcomp->reac(lridx);
    AssertLog(lreac->defr() == &comp.reacdef(lridx));

    return lreac->h();
}

////////////////////////////////////////////////////////////////////////

long double Wmdirect::_getCompReacA(solver::comp_global_id cidx,
                                    solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    const auto& comp = statedef().compdef(cidx);
    solver::reac_local_id lridx = comp.reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx.get()];
    AssertLog(lcomp->def() == &comp);
    // Reacdef local indices in Compdef object also have same index
    // as Reacs in Comp object
    KProc* lreac = lcomp->reac(lridx);
    AssertLog(lreac->defr() == &comp.reacdef(lridx));

    return static_cast<long double>(lreac->rate());
}

////////////////////////////////////////////////////////////////////////

unsigned long long Wmdirect::_getCompReacExtent(solver::comp_global_id cidx,
                                                solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    const auto& comp = statedef().compdef(cidx);
    solver::reac_local_id lridx = comp.reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx.get()];
    AssertLog(lcomp->def() == &comp);
    // Reacdef local indices in Compdef object also have same index
    // as Reacs in Comp object
    KProc* lreac = lcomp->reac(lridx);
    AssertLog(lreac->defr() == &comp.reacdef(lridx));

    return lreac->getExtent();
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_resetCompReacExtent(solver::comp_global_id cidx, solver::reac_global_id ridx) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    auto& comp = statedef().compdef(cidx);
    solver::reac_local_id lridx = comp.reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx.get()];
    AssertLog(lcomp->def() == &comp);
    // Reacdef local indices in Compdef object also have same index
    // as Reacs in Comp object
    KProc* lreac = lcomp->reac(lridx);
    AssertLog(lreac->defr() == &comp.reacdef(lridx));

    lreac->resetExtent();
}

////////////////////////////////////////////////////////////////////////

unsigned long long Wmdirect::_getCompComplexReacExtent(solver::comp_global_id cidx,
                                                       solver::complexreac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countComplexReacs());
    solver::Compdef& comp = statedef().compdef(cidx);
    solver::complexreac_local_id lridx = comp.complexreacG2L(ridx);
    if (lridx.unknown()) {
        ArgErrLog("Complex reaction undefined in compartment.\n");
    }

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx.get()];
    AssertLog(lcomp->def() == &comp);
    // Reacdef local indices in Compdef object also have same index
    // as Reacs in Comp object
    KProc* lreac = lcomp->reac(lridx);
    AssertLog(&lreac->defcr() == &comp.complexreacdef(lridx));

    return lreac->getExtent();
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getPatchArea(solver::patch_global_id pidx) const {
    AssertLog(pidx < statedef().countPatches());
    return statedef().patchdef(pidx).area();
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setPatchArea(solver::patch_global_id pidx, double area) {
    AssertLog(pidx < statedef().countPatches());
    statedef().patchdef(pidx).setArea(area);
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getPatchSpecCount(solver::patch_global_id pidx,
                                    solver::spec_global_id sidx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countSpecs());
    const auto& patch = statedef().patchdef(pidx);
    solver::spec_local_id slidx = patch.specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    return patch.pools()[slidx];
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setPatchSpecCount(solver::patch_global_id pidx,
                                  solver::spec_global_id sidx,
                                  double n) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countSpecs());
    auto& patch = statedef().patchdef(pidx);
    solver::spec_local_id slidx = patch.specG2L(sidx);
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
    patch.setCount(slidx, n_int);
    _reset();
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getPatchSpecAmount(solver::patch_global_id pidx,
                                     solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getPatchSpecCount(pidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setPatchSpecAmount(solver::patch_global_id pidx,
                                   solver::spec_global_id sidx,
                                   double a) {
    AssertLog(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setPatchSpecCount(pidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////

bool Wmdirect::_getPatchSpecClamped(solver::patch_global_id pidx,
                                    solver::spec_global_id sidx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countSpecs());
    const auto& patch = statedef().patchdef(pidx);
    solver::spec_local_id lsidx = patch.specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    return patch.clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setPatchSpecClamped(solver::patch_global_id pidx,
                                    solver::spec_global_id sidx,
                                    bool buf) {
    AssertLog(pidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    auto& patch = statedef().patchdef(pidx);
    solver::spec_local_id lsidx = patch.specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    patch.setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getPatchComplexCount(
    solver::patch_global_id pidx,
    solver::complex_global_id cmplIdx,
    const std::vector<
        util::strongid_vector<solver::complex_substate_id, model::SubunitStateFilter>>& f) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(cmplIdx < statedef().countComplexes());
    solver::Patchdef& patch = statedef().patchdef(pidx);

    auto filt = patch.GetFilter(cmplIdx, f);
    const auto& states = patch.complexStates(cmplIdx);
    filt->processUpdates(states);
    return filt->nbMatches(states);
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setPatchComplexCount(
    solver::patch_global_id pidx,
    solver::complex_global_id cmplIdx,
    const util::strongid_vector<solver::complex_substate_id, uint>& i,
    double n) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(cmplIdx < statedef().countComplexes());
    solver::Patchdef& patch = statedef().patchdef(pidx);

    if (n > std::numeric_limits<unsigned int>::max()) {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << std::numeric_limits<unsigned int>::max() << ").\n";
        ArgErrLog(os.str());
    }

    AssertLog(n >= 0.0);
    double n_int = std::floor(n);
    double n_frc = n - n_int;
    uint c = static_cast<uint>(n_int);
    if (n_frc > 0.0) {
        double rand01 = rng()->getUnfIE();
        if (rand01 < n_frc) {
            c++;
        }
    }

    std::vector<solver::complex_individual_id> existingInds;
    for (auto state: patch.complexStates(cmplIdx)) {
        if (state.second == i) {
            existingInds.push_back(state.first);
        }
    }
    // If we need to remove some of the complexes
    if (existingInds.size() > c) {
        uint nbToDel = existingInds.size() - c;
        for (auto it = existingInds.begin(); it != existingInds.begin() + nbToDel; ++it) {
            patch.removeComplex(cmplIdx, *it);
        }
    } else {
        for (uint nb = 0; nb < c - existingInds.size(); ++nb) {
            patch.addComplex(cmplIdx, getNextComplexInd(), i);
        }
    }

    // Rates have changed
    _reset();
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getPatchComplexAmount(
    solver::patch_global_id pidx,
    solver::complex_global_id sidx,
    const std::vector<
        util::strongid_vector<solver::complex_substate_id, model::SubunitStateFilter>>& f) const {
    // the following method does all the necessary argument checking
    double count = _getPatchComplexCount(pidx, sidx, f);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setPatchComplexAmount(
    solver::patch_global_id pidx,
    solver::complex_global_id sidx,
    const util::strongid_vector<solver::complex_substate_id, uint>& i,
    double a) {
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setPatchComplexCount(pidx, sidx, i, a2);
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getPatchComplexSUSCount(
    solver::patch_global_id pidx,
    solver::complex_global_id sidx,
    const std::vector<
        util::strongid_vector<solver::complex_substate_id, model::SubunitStateFilter>>& f,
    solver::complex_substate_id m) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countComplexes());
    solver::Patchdef& patch = statedef().patchdef(pidx);

    auto filt = patch.GetFilter(sidx, f);
    filt->processUpdates(patch.complexStates(sidx));
    return filt->nbSus(m, patch.complexStates(sidx));
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getPatchSReacK(solver::patch_global_id pidx, solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    const auto& patch = statedef().patchdef(pidx);
    solver::sreac_local_id lridx = patch.sreacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    return patch.kcst(lridx);
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setPatchSReacK(solver::patch_global_id pidx,
                               solver::sreac_global_id ridx,
                               double kf) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    AssertLog(kf >= 0.0);
    auto& patch = statedef().patchdef(pidx);
    solver::sreac_local_id lridx = patch.sreacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    patch.setKcst(lridx, kf);

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx.get()];
    AssertLog(lpatch->def() == &patch);
    // SReacdef local indices in Patchdef object also have same index
    //  as SReacs in Patch object
    KProc* lsreac = lpatch->sreac(lridx);
    AssertLog(lsreac->defsr() == &patch.sreacdef(lridx));
    lsreac->resetCcst();

    // Rates have changed
    _reset();
}

////////////////////////////////////////////////////////////////////////

bool Wmdirect::_getPatchSReacActive(solver::patch_global_id pidx,
                                    solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    const auto& patch = statedef().patchdef(pidx);
    solver::sreac_local_id lridx = patch.sreacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    return patch.active(lridx);
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_setPatchSReacActive(solver::patch_global_id pidx,
                                    solver::sreac_global_id ridx,
                                    bool a) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    auto& patch = statedef().patchdef(pidx);
    solver::sreac_local_id lridx = patch.sreacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    patch.setActive(lridx, a);

    // It's cheaper to just recompute everything
    _reset();
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getPatchSReacC(solver::patch_global_id pidx, solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());

    const auto& patch = statedef().patchdef(pidx);
    solver::sreac_local_id lridx = patch.sreacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx.get()];
    AssertLog(lpatch->def() == &patch);
    // SReacdef local indices in Patchdef object also have same index
    //  as SReacs in Patch object
    KProc* lsreac = lpatch->sreac(lridx);
    AssertLog(lsreac->defsr() == &patch.sreacdef(lridx));

    return lsreac->c();
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getPatchSReacH(solver::patch_global_id pidx, solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());

    const auto& patch = statedef().patchdef(pidx);
    solver::sreac_local_id lridx = patch.sreacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx.get()];
    AssertLog(lpatch->def() == &patch);
    // SReacdef local indices in Patchdef object also have same index
    //  as SReacs in Patch object
    KProc* lsreac = lpatch->sreac(lridx);
    AssertLog(lsreac->defsr() == &patch.sreacdef(lridx));

    return lsreac->h();
}

////////////////////////////////////////////////////////////////////////

double Wmdirect::_getPatchSReacA(solver::patch_global_id pidx, solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());

    const auto& patch = statedef().patchdef(pidx);
    solver::sreac_local_id lridx = patch.sreacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx.get()];
    AssertLog(lpatch->def() == &patch);
    // SReacdef local indices in Patchdef object also have same index
    //  as SReacs in Patch object
    KProc* lsreac = lpatch->sreac(lridx);
    AssertLog(lsreac->defsr() == &patch.sreacdef(lridx));

    return lsreac->rate();
}

////////////////////////////////////////////////////////////////////////

unsigned long long Wmdirect::_getPatchSReacExtent(solver::patch_global_id pidx,
                                                  solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    const auto& patch = statedef().patchdef(pidx);
    solver::sreac_local_id lridx = patch.sreacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx.get()];
    AssertLog(lpatch->def() == &patch);
    // SReacdef local indices in Patchdef object also have same index
    //  as SReacs in Patch object
    KProc* lsreac = lpatch->sreac(lridx);
    AssertLog(lsreac->defsr() == &patch.sreacdef(lridx));

    return lsreac->getExtent();
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_resetPatchSReacExtent(solver::patch_global_id pidx, solver::sreac_global_id ridx) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    auto& patch = statedef().patchdef(pidx);
    solver::sreac_local_id lridx = patch.sreacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx.get()];
    AssertLog(lpatch->def() == &patch);
    // SReacdef local indices in Patchdef object also have same index
    //  as SReacs in Patch object
    KProc* lsreac = lpatch->sreac(lridx);
    AssertLog(lsreac->defsr() == &patch.sreacdef(lridx));

    lsreac->resetExtent();
}

////////////////////////////////////////////////////////////////////////

unsigned long long Wmdirect::_getPatchComplexSReacExtent(
    solver::patch_global_id pidx,
    solver::complexsreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countComplexSReacs());
    solver::Patchdef& patch = statedef().patchdef(pidx);
    solver::complexsreac_local_id lridx = patch.complexsreacG2L(ridx);
    if (lridx.unknown()) {
        ArgErrLog("Complex surface reaction undefined in patch.\n");
    }

    Patch* lpatch = pPatches[pidx.get()];
    AssertLog(lpatch->def() == &patch);
    KProc* lreac = lpatch->sreac(lridx);
    AssertLog(&lreac->defcsr() == &patch.complexsreacdef(lridx));

    return lreac->getExtent();
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::addKProc(KProc* kp) {
    AssertLog(kp != nullptr);

    SchedIDX nidx(static_cast<uint>(pKProcs.size()));  // because pKProcs.size() is ulong
    pKProcs.push_back(kp);
    kp->setSchedIDX(nidx);
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_build() {
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
        auto level = new double[clsize];
        std::fill_n(level, clsize, 0.0);
        pLevelSizes.push_back(clsize);
        pLevels.push_back(level);

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

KProc* Wmdirect::_getNext() {
    AssertLog(pA0 >= 0.0);
    // Stefan's remark on using one random number for each branch instead of just
    // one for the whole search
    /*
    //The actual implementation is slightly more complicated because we need take
    care of issues that arise
    //from the limited precision of floating point arithmetic. {These issues,
    while of a technical nature and
    //not really crucial, are interesting enough to expand but I’ll do this later.
    They arise, in different
    //ways, in every implementation of SSA and are therefore of general interest.
    */

    // Quick check to see whether nothing is there.
    if (pA0 == 0.0) {
        return nullptr;
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
    double* level = nullptr;
    while (clevel != 0) {
        // Decrease the current level.
        clevel--;
        // and start looking in the right place.
        cur_node *= SCHEDULEWIDTH;
        uint max_node = cur_node + SCHEDULEWIDTH;

        // Fetch the level.
        level = pLevels[clevel];

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
    return pKProcs[cur_node];
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_reset() {
    if (pKProcs.size() == 0) {
        return;
    }

    // Reset the basic level: compute rates.
    double* oldlevel = pLevels[0];
    uint cur_node = 0;
    for (auto const& kp: pKProcs) {
        oldlevel[cur_node++] = kp->rate();
    }
    _clearComplexFilterUpdates();

    // Work up.
    for (uint cur_level = 1; cur_level < pLevels.size(); ++cur_level) {
        // Compute the number of nodes to reset on this level.
        uint numnodes = pLevelSizes[cur_level - 1] / SCHEDULEWIDTH;

        // Fetch a pointer to this level.
        double* level = pLevels[cur_level];

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

void Wmdirect::_update(SchedIDXVec const& entries) {
    if (countKProcs() == 0) {
        return;
    }

    // Prefetch zero level.
    double* level0 = pLevels[0];
    // Number of entries.
    AssertLog(entries.size() <= pMaxUpSize);  /////////

    // Recompute rates.
    uint prev_e = 0xFFFFFFFF;
    uint cur_e = 0;
    for (auto const& sidx: entries) {
        // Fetch index.
        uint idx = sidx.get();
        // Recompute rate, get difference, and store.
        double newrate = pKProcs[idx]->rate();
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
    _clearComplexFilterUpdates();

    // Update upper levels.
    uint nlevels = pLevels.size();
    double* prevlevel = pLevels[0];
    for (uint l = 1; l < nlevels; ++l) {
        // Update the first entry.
        cur_e = 0;
        prev_e = 0xFFFFFFFF;

        // Fetch a pointer to the current level.
        double* currlevel = pLevels[l];

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
    double* toplevel = pLevels[pLevels.size() - 1];
    pA0 = 0.0;
    for (uint i = 0; i < SCHEDULEWIDTH; ++i) {
        pA0 += toplevel[i];
    }
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_clearComplexFilterUpdates() {
    uint comps = statedef().countComps();
    uint patches = statedef().countPatches();
    for (auto c: solver::comp_global_id::range(comps)) {
        statedef().compdef(c).clearComplexFilterUpdates();
    }
    for (auto p: solver::patch_global_id::range(patches)) {
        statedef().patchdef(p).clearComplexFilterUpdates();
    }
}

////////////////////////////////////////////////////////////////////////

void Wmdirect::_executeStep(KProc* kp, double dt) {
    SchedIDXVec const& upd = kp->apply();
    _update(upd);
    statedef().incTime(dt);
    statedef().incNSteps(1);
}

}  // namespace steps::wmdirect
