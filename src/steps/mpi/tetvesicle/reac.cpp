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

#include "mpi/tetvesicle/reac.hpp"

// STEPS headers.
#include "math/constants.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/tetvesicle_rdef.hpp"
#include "solver/reacdef.hpp"
#include "util/checkpointing.hpp"

namespace steps::mpi::tetvesicle {

static inline double comp_ccst(double kcst, double vol, uint order) {
    double vscale = 1.0e3 * vol * math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;

    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;

    double ccst = kcst * pow(vscale, static_cast<double>(-o1));

    return ccst;
}

////////////////////////////////////////////////////////////////////////////////

Reac::Reac(solver::Reacdef* rdef, TetRDEF* tet)
    : pReacdef(rdef)
    , pTet(tet)
    , pCcst(0.0)
    , pKcst(0.0) {
    AssertLog(pReacdef != nullptr);
    AssertLog(pTet != nullptr);
    pType = KP_REAC;

    solver::reac_local_id lridx = pTet->compdef()->reacG2L(pReacdef->gidx());
    double kcst = pTet->compdef()->kcst(lridx);
    pKcst = kcst;
    pCcst = comp_ccst(kcst, pTet->vol(), pReacdef->order());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

Reac::~Reac() = default;

////////////////////////////////////////////////////////////////////////////////

void Reac::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pCcst);
    util::checkpoint(cp_file, pKcst);
    KProc::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void Reac::restore(std::fstream& cp_file) {
    util::restore(cp_file, pCcst);
    util::restore(cp_file, pKcst);
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void Reac::reset() {
    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    resetExtent();
    resetCcst();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void Reac::resetCcst() {
    solver::reac_local_id lridx = pTet->compdef()->reacG2L(pReacdef->gidx());
    double kcst = pTet->compdef()->kcst(lridx);
    // Also reset kcst
    pKcst = kcst;
    pCcst = comp_ccst(kcst, pTet->vol(), pReacdef->order());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void Reac::setKcst(double k) {
    AssertLog(k >= 0.0);
    pKcst = k;
    pCcst = comp_ccst(k, pTet->vol(), pReacdef->order());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void Reac::setupDeps() {
    AssertLog(pTet->getInHost());

    KProcPSet updset;

    // for (auto const & s : pReacdef->updColl()) pSpecChange.insert(s);

    // Search in local tetrahedron.
    uint nkprocs = pTet->countKProcs();

    for (uint k = 0; k < nkprocs; k++) {
        for (auto const& s: pReacdef->updColl()) {
            if (pTet->KProcDepSpecTet(k, pTet, s)) {
                updset.insert(pTet->getKProc(k));
                break;  // exit the species loop because we can only add this kp once
            }
        }
    }

    for (auto const& tri: pTet->nexttris()) {
        if (tri == nullptr) {
            continue;
        }

        if (pTet->getHost() != tri->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << tri->idx() << " and its compartment tetrahedron "
               << pTet->idx() << " belong to different hosts.\n";
            NotImplErrLog(os.str());
        }

        nkprocs = tri->countKProcs();
        for (uint sk = 0; sk < nkprocs; sk++) {
            for (auto const& s: pReacdef->updColl()) {
                if (tri->KProcDepSpecTet(sk, pTet, s)) {
                    updset.insert(tri->getKProc(sk));
                    break;  // exit the species loop because we can only add this kp once
                }
            }
        }
    }

    localUpdVec.assign(updset.begin(), updset.end());
}

////////////////////////////////////////////////////////////////////////////////

double Reac::rate(TetVesicleRDEF* /*solver*/) {
    if (inactive()) {
        return 0.0;
    }

    // This is necessary for TetVesicleRDEF solver since volume may have changed
    pCcst = comp_ccst(pKcst, pTet->vol(), pReacdef->order());

    // Prefetch some variables.
    solver::Compdef* cdef = pTet->compdef();
    const auto& lhs_vec = cdef->reac_lhs(cdef->reacG2L(pReacdef->gidx()));
    const auto& cnt_vec = pTet->pools();

    // Compute combinatorial part.
    double h_mu = 1.0;
    for (auto pool: lhs_vec.range()) {
        uint lhs = lhs_vec[pool];
        if (lhs == 0) {
            continue;
        }
        uint cnt = cnt_vec[pool];
        if (lhs > cnt) {
            h_mu = 0.0;
            break;
        }
        switch (lhs) {
        case 4: {
            h_mu *= static_cast<double>(cnt - 3);
        }
        case 3: {
            h_mu *= static_cast<double>(cnt - 2);
        }
        case 2: {
            h_mu *= static_cast<double>(cnt - 1);
        }
        case 1: {
            h_mu *= static_cast<double>(cnt);
            break;
        }
        default: {
            AssertLog(0);
            return 0.0;
        }
        }
    }

    // Multiply with scaled reaction constant.
    return h_mu * pCcst;
}

////////////////////////////////////////////////////////////////////////////////

void Reac::apply(const rng::RNGptr& /*rng*/,
                 double /*dt*/,
                 double /*simtime*/,
                 double period,
                 TetVesicleRDEF* /*solver*/) {
    const auto& local = pTet->pools();
    solver::Compdef* cdef = pTet->compdef();
    solver::reac_local_id l_ridx = cdef->reacG2L(pReacdef->gidx());
    const auto& upd_vec = cdef->reac_upd(l_ridx);
    for (auto i: upd_vec.range()) {
        if (pTet->clamped(i) == true) {
            continue;
        }
        int j = upd_vec[i];
        if (j == 0) {
            continue;
        }
        int nc = static_cast<int>(local[i]) + j;
        AssertLog(nc >= 0);
        pTet->setCount(i, static_cast<uint>(nc), period);
    }

    rExtent++;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::kproc_global_id> const& Reac::getRemoteUpdVec(int /*direction*/) const {
    return remoteUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<KProc*> const& Reac::getLocalUpdVec(int /*direction*/) const {
    return localUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

void Reac::resetOccupancies() {
    pTet->resetPoolOccupancy();
}

}  // namespace steps::mpi::tetvesicle
