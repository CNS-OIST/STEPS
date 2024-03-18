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

#include "mpi/tetvesicle/raftgen.hpp"

// STEPS headers.
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/tetvesicle_rdef.hpp"
#include "util/checkpointing.hpp"


namespace steps::mpi::tetvesicle {

RaftGen::RaftGen(solver::RaftGendef* rgdef, TriRDEF* tri)
    : pRaftGendef(rgdef)
    , pTri(tri) {
    AssertLog(pRaftGendef != nullptr);
    AssertLog(pTri != nullptr);

    pType = KP_RAFTGEN;

    double kcst = pRaftGendef->kcst();
    AssertLog(kcst >= 0.0);
    pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

RaftGen::~RaftGen() = default;

////////////////////////////////////////////////////////////////////////////////

void RaftGen::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pKcst);
    KProc::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void RaftGen::restore(std::fstream& cp_file) {
    util::restore(cp_file, pKcst);
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void RaftGen::reset() {
    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    resetExtent();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void RaftGen::setKcst(double k) {
    AssertLog(k >= 0.0);
    pKcst = k;
}

////////////////////////////////////////////////////////////////////////////////

void RaftGen::setupDeps() {
    // Need to do this because a raftgen removes species from triangle
    // during SSA

    AssertLog(pTri->getInHost());

    TetRDEF* itet = pTri->iTet();
    TetRDEF* otet = pTri->oTet();

    // Raftgendef doesn't have an 'update column' but there is no rhs, so
    // anything on lhs  will change upon application, held in dep_S

    uint nspecs_g = def()->countSpecs_global();
    // Anything in the triangles that depend on these
    // species have to be included for updates
    KProcPSet updset;

    uint nkprocs = pTri->countKProcs();
    // check if sk KProc in pTri depends on spec in pTri
    for (uint sk = 0; sk < nkprocs; sk++) {
        for (auto s: solver::spec_global_id::range(nspecs_g)) {
            if (def()->dep_S(s) != 0) {
                if (pTri->KProcDepSpecTri(sk, pTri, s)) {
                    updset.insert(pTri->getKProc(sk));
                }
            }
        }
    }

    if (itet != nullptr) {
        if (pTri->getHost() != itet->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron ";
            os << itet->idx() << " belong to different hosts.\n";
            NotImplErrLog(os.str());
        }

        // Checking Kprocs in inner tet that might depend on spec change in that tri
        nkprocs = itet->countKProcs();
        for (uint k = 0; k < nkprocs; k++) {
            // Changes in surface species affect reactions like vessreacs
            for (auto s: solver::spec_global_id::range(nspecs_g)) {
                if (def()->dep_S(s) != 0) {
                    if (itet->KProcDepSpecTri(k, pTri, s)) {
                        updset.insert(itet->getKProc(k));
                        break;  // exit the species loop because we can only add this kp once
                    }
                }
            }
        }
    }

    if (otet != nullptr) {
        if (pTri->getHost() != otet->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron ";
            os << otet->idx() << " belong to different hosts.\n";
            NotImplErrLog(os.str());
        }

        // Checking Kprocs in outer tet that might depend on spec change in that tri
        nkprocs = otet->countKProcs();
        for (uint k = 0; k < nkprocs; k++) {
            // Changes in surface species affect reactions like vessreacs
            for (auto s: solver::spec_global_id::range(nspecs_g)) {
                if (def()->dep_S(s) != 0) {
                    if (otet->KProcDepSpecTri(k, pTri, s)) {
                        updset.insert(otet->getKProc(k));
                        break;  // exit the species loop because we can only add this kp once
                    }
                }
            }
        }
    }

    localUpdVec.assign(updset.begin(), updset.end());

    return;
}

////////////////////////////////////////////////////////////////////////////////

double RaftGen::rate(TetVesicleRDEF* /*solver*/) {
    if (inactive()) {
        return 0.0;
    }

    solver::Patchdef* pdef = pTri->patchdef();

    // This uses global indices already

    const auto& lhs_s_vec = def()->lhs_S();

    for (auto sg: solver::spec_global_id::range(lhs_s_vec.size())) {
        uint lhs = lhs_s_vec[sg];
        if (lhs == 0) {
            continue;
        }

        // We need a spec or specs.
        solver::spec_local_id spec_lidx = pdef->specG2L(sg);

        // For development only. RaftGens should add their specs to patch
        AssertLog(spec_lidx.valid());

        uint cnt = pTri->pools()[spec_lidx];

        //  Compare to required lhs
        if (lhs > cnt) {
            //  The required species are not available
            return 0.0;
        }
    }

    return pKcst;
}

////////////////////////////////////////////////////////////////////////////////

void RaftGen::apply(const rng::RNGptr& /*rng*/,
                    double /*dt*/,
                    double /*simtime*/,
                    double period,
                    TetVesicleRDEF* /*solver*/) {
    // let the tri function do all the work- flag the raftgen for
    // true application during patch routine but also remove the species from
    // tri so they can't be removed by anything else in the mean time.
    tri()->applyRaftGen(def(), period);

    rExtent++;
}

////////////////////////////////////////////////////////////////////////////////

void RaftGen::resetOccupancies() {
    pTri->resetPoolOccupancy();
}

}  // namespace steps::mpi::tetvesicle
