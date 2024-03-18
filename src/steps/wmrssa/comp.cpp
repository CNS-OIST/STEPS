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
#include "comp.hpp"
#include "reac.hpp"
#include "wmrssa.hpp"

// logging
#include "util/error.hpp"

#include "util/checkpointing.hpp"

namespace steps::wmrssa {

////////////////////////////////////////////////////////////////////////////////

Comp::Comp(solver::Compdef* compdef)
    : pCompdef(compdef) {
    assert(pCompdef != nullptr);
    uint nspecs = compdef->countSpecs();
    pPoolLB.container().resize(nspecs);
    pPoolUB.container().resize(nspecs);
}

////////////////////////////////////////////////////////////////////////////////

Comp::~Comp() {
    for (auto const& k: pKProcs) {
        delete k;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Comp::checkpoint(std::fstream& cp_file) {
    for (auto const& k: pKProcs) {
        k->checkpoint(cp_file);
    }
    util::checkpoint(cp_file, pPoolLB);
    util::checkpoint(cp_file, pPoolUB);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::restore(std::fstream& cp_file) {
    for (auto const& k: pKProcs) {
        k->restore(cp_file);
    }
    util::restore(cp_file, pPoolLB);
    util::restore(cp_file, pPoolUB);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::reset() {
    for (auto const& k: pKProcs) {
        k->reset();
    }
}

////////////////////////////////////////////////////////////////////////////////

void Comp::setupKProcs(Wmrssa* wmd) {
    // Create reaction kproc's.
    uint nreacs = def()->countReacs();
    pKProcs.resize(nreacs);
    for (auto i: solver::reac_local_id::range(nreacs)) {
        auto& rdef = def()->reacdef(i);
        auto* r = new Reac(&rdef, this);
        pKProcs[i.get()] = r;
        wmd->addKProc(r);
    }
}

////////////////////////////////////////////////////////////////////////////////

wmrssa::KProc* Comp::reac(solver::reac_local_id lridx) const {
    assert(lridx.get() < pKProcs.size());
    return pKProcs[lridx.get()];
}

////////////////////////////////////////////////////////////////////////////////

void Comp::setupDeps() {
    for (auto const& k: pKProcs) {
        k->setupDeps();
    }
}

////////////////////////////////////////////////////////////////////////////////

void Comp::addIPatch(Patch* p) {
    AssertLog(std::find(pIPatches.begin(), pIPatches.end(), p) == pIPatches.end());
    pIPatches.push_back(p);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::addOPatch(Patch* p) {
    AssertLog(std::find(pOPatches.begin(), pOPatches.end(), p) == pOPatches.end());
    pOPatches.push_back(p);
}
////////////////////////////////////////////////////////////////////////////////

void Comp::setBounds(solver::spec_local_id i, int nc) {
    constexpr double delta = .05;
    if (nc > 3 / delta) {
        pPoolLB[i] = nc * (1 - delta);
        pPoolUB[i] = nc * (1 + delta);
    } else if (nc > 3) {
        pPoolLB[i] = nc - 3;
        pPoolUB[i] = nc + 3;
    } else if (nc > 0) {
        pPoolLB[i] = 1;
        pPoolUB[i] = 2 * nc;
    } else {
        pPoolLB[i] = 0;
        pPoolUB[i] = 0;
    }
    pPoolLB[i] -= delta;
    pPoolUB[i] += delta;
    // pPoolLB[i] = std::max(std::min(nc*(1 - delta), nc - 3.), 0.); // nc/(7 -
    // delta - 2*(3 - delta)/(1 + 2./(nc+1))); //*(1 - delta);// pPoolUB[i] = nc >
    // 0 ? nc*(7 - delta - 2*(3 - delta)/(1 + 2./(nc+1))) : 3; //std::max(nc*(1 +
    // delta), nc + 3.); // (1 + delta);//
}
////////////////////////////////////////////////////////////////////////////////

bool Comp::isOutOfBound(solver::spec_local_id i, int nc) {
    AssertLog(i < def()->countSpecs());
    if (nc > pPoolLB[i] && nc < pPoolUB[i]) {
        return false;
    }
    setBounds(i, nc);
    return true;
}

////////////////////////////////////////////////////////////////////////////////

const util::strongid_vector<solver::spec_local_id, double>& Comp::pools(
    wmrssa::PropensityRSSA prssa) const {
    switch (prssa) {
    case wmrssa::CURRENT:
        return def()->pools();
    case wmrssa::LOWERBOUND:
        return pPoolLB;
    case wmrssa::BOUNDS:
        return pPoolUB;
    default:
        AssertLog(false);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Comp::setupSpecDeps() {
    uint nspecs = def()->countSpecs();
    localSpecUpdKProcs.container().resize(nspecs);
    for (auto slidx: solver::spec_local_id::range(nspecs)) {
        solver::spec_global_id sgidx = def()->specL2G(slidx);
        for (auto const& k: pKProcs) {
            if (k->depSpecComp(sgidx, this)) {
                localSpecUpdKProcs[slidx].push_back(k);
            }
        }
        for (auto const& ip: pIPatches) {
            for (auto const& k: ip->kprocs()) {
                if (k->depSpecComp(sgidx, this)) {
                    localSpecUpdKProcs[slidx].push_back(k);
                }
            }
        }
        for (auto const& ip: pOPatches) {
            for (auto const& k: ip->kprocs()) {
                if (k->depSpecComp(sgidx, this)) {
                    localSpecUpdKProcs[slidx].push_back(k);
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<KProc*> const& Comp::getSpecUpdKProcs(solver::spec_local_id slidx) {
    return localSpecUpdKProcs[slidx];
}

}  // namespace steps::wmrssa
