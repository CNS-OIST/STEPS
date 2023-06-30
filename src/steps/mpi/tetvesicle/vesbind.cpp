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

#include "mpi/tetvesicle/vesbind.hpp"

// STEPS headers.
#include "math/constants.hpp"
#include "mpi/mpi_common.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/tetvesicle_rdef.hpp"
#include "mpi/tetvesicle/vesproxy.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"

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

VesBind::VesBind(solver::VesBinddef* vbdef, TetRDEF* tet)
    : pVesBinddef(vbdef)
    , pTet(tet)
    , pCcst(0.0)
    , pKcst(0.0) {
    AssertLog(pVesBinddef != nullptr);
    AssertLog(pTet != nullptr);
    pType = KP_VESBIND;

    double kcst = pVesBinddef->kcst();
    pKcst = kcst;

    pCcst = comp_ccst(kcst, pTet->vol(), pVesBinddef->order());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

VesBind::~VesBind() = default;

////////////////////////////////////////////////////////////////////////////////

void VesBind::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pCcst);
    util::checkpoint(cp_file, pKcst);
    KProc::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void VesBind::restore(std::fstream& cp_file) {
    util::restore(cp_file, pCcst);
    util::restore(cp_file, pKcst);
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void VesBind::reset() {
    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    resetExtent();
    resetCcst();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void VesBind::resetCcst() {
    double kcst = pVesBinddef->kcst();
    pKcst = kcst;
    pCcst = comp_ccst(kcst, pTet->vol(), pVesBinddef->order());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void VesBind::setKcst(double k) {
    AssertLog(k >= 0.0);
    pKcst = k;
    pCcst = comp_ccst(k, pTet->vol(), pVesBinddef->order());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void VesBind::setupDeps() {
    AssertLog(pTet->getInHost());

    // And the rest is for static kproc updates
    KProcPSet updset;

    // An application of vesbind can affect other
    // vesbind reactions and vesunbind reactions that belong to the tet.

    uint nkprocs = pTet->countKProcs();

    for (uint k = 0; k < nkprocs; k++) {
        if (pTet->KProcDepSpecTetVesSurface(k, pTet, def()->getSpec1gidx()) ||
            pTet->KProcDepSpecTetVesSurface(k, pTet, def()->getSpec2gidx()) ||
            pTet->KProcDepLinkSpecTetVesSurface(k, pTet, def()->getLinkSpec1gidx()) ||
            pTet->KProcDepLinkSpecTetVesSurface(k, pTet, def()->getLinkSpec2gidx())) {
            updset.insert(pTet->getKProc(k));
        }
    }

    localUpdVec.assign(updset.begin(), updset.end());
}

////////////////////////////////////////////////////////////////////////////////

double VesBind::rate(TetVesicleRDEF* /*solver*/) {
    if (inactive()) {
        return 0.0;
    }

    // This is necessary for TetVesicleRDEF solver since volume may have changed
    pCcst = comp_ccst(pKcst, pTet->vol(), pVesBinddef->order());

    // Need to check if we have an active maximum distance
    double max_distance = def()->max_distance();
    double min_distance = def()->min_distance();

    // 1 Need to find out which vesicles are in the tet. Compare to necessary
    // vesicles for this binding
    auto const& vesproxyrefs = pTet->getVesProxyrefs();

    // Need at least 2 vesicles
    if (vesproxyrefs.size() < 2) {
        return 0.0;
    }
    uint pairs = 0;

    solver::spec_global_id spec1_gidx = def()->getSpec1gidx();
    solver::spec_global_id spec2_gidx = def()->getSpec2gidx();
    solver::vesicle_global_id ves1_idx = def()->getVes1idx();
    solver::vesicle_global_id ves2_idx = def()->getVes2idx();

    uint ngspecs = def()->countSpecsGlobal();
    uint ngbspecs = def()->countLinkSpecsGlobal();

    // Write these out as separately in 'if else' at this level for clarity

    for (auto const& vesproxy1_mapit: vesproxyrefs) {
        auto const& vesproxy1 = vesproxy1_mapit.second;

        // Need to check if this is type Vesicle 1
        // Don't allow if vesicle has been earmarked for exocytosis at next update
        if (vesproxy1->idx() == ves1_idx && vesproxy1->exoApplied().unknown()) {
            // First need to check any species dependencies
            bool allspecdepsgood = true;

            if (def()->gotVDep1()) {
                for (auto s: solver::spec_global_id::range(ngspecs)) {
                    uint vdep1 = def()->vdep1(s);
                    if (vdep1 == 0) {
                        continue;
                    }
                    // We got a dependency species. Need to see if it's on the vesicle
                    // surface or not.
                    uint cnt = vesproxy1->getSpecCount_V(s);
                    if (cnt >= vdep1) {
                        continue;
                    } else {
                        allspecdepsgood = false;
                        break;
                    }
                }
            }

            if (allspecdepsgood and def()->gotLDep1()) {
                for (auto s: solver::linkspec_global_id::range(ngbspecs)) {
                    uint ldep1 = def()->ldep1(s);
                    if (ldep1 == 0) {
                        continue;
                    }
                    uint cnt = vesproxy1->getLinkSpecCount_V(s);
                    if (cnt >= ldep1) {
                        continue;
                    } else {
                        allspecdepsgood = false;
                        break;
                    }
                }
            }

            // Go to next vesicle if this one doesn't have the right dependencies
            if (!allspecdepsgood) {
                continue;
            }
            // Positions absolute, not relative to centre
            const auto& spec1_positions_abs = vesproxy1->getSurfaceSpecPos_absolute(spec1_gidx);

            for (auto const& spec1_pos_abs: spec1_positions_abs) {
                for (auto const& vesproxy2_mapit: vesproxyrefs) {
                    auto const& vesproxy2 = vesproxy2_mapit.second;

                    if (vesproxy1 == vesproxy2) {
                        continue;
                    }
                    if ((vesproxy2->idx() == ves2_idx) && (vesproxy2->exoApplied().unknown())) {
                        // Check dependencies first
                        if (def()->gotVDep2()) {
                            for (auto s: solver::spec_global_id::range(ngspecs)) {
                                uint vdep2 = def()->vdep2(s);
                                if (vdep2 == 0) {
                                    continue;
                                }
                                uint cnt = vesproxy2->getSpecCount_V(s);
                                if (cnt >= vdep2) {
                                    continue;
                                } else {
                                    allspecdepsgood = false;
                                    break;
                                }
                            }
                        }

                        if (allspecdepsgood and def()->gotLDep2()) {
                            for (auto s: solver::linkspec_global_id::range(ngbspecs)) {
                                uint ldep2 = def()->ldep2(s);
                                if (ldep2 == 0) {
                                    continue;
                                }
                                uint cnt = vesproxy2->getLinkSpecCount_V(s);
                                if (cnt >= ldep2) {
                                    continue;
                                } else {
                                    allspecdepsgood = false;
                                    break;
                                }
                            }
                        }

                        // Go to next vesicle if this one doesn't have the right
                        // dependencies
                        if (!allspecdepsgood) {
                            continue;
                        }
                        // Positions absolute, not relative to centre
                        const auto& spec2_positions_abs = vesproxy2->getSurfaceSpecPos_absolute(
                            spec2_gidx);

                        for (auto const& spec2_pos_abs: spec2_positions_abs) {
                            double dist = distance(spec1_pos_abs.second, spec2_pos_abs.second);
                            if (dist <= max_distance && dist >= min_distance) {
                                pairs += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    return pairs * pCcst;
}

////////////////////////////////////////////////////////////////////////////////

void VesBind::apply(const rng::RNGptr& /*rng*/,
                    double /*dt*/,
                    double /*simtime*/,
                    double /*period*/,
                    TetVesicleRDEF* /*solver*/) {
    // 1 Need to find out which vesicles are in the tet. Compare to necessary
    // vesicles for this binding
    auto const& vesproxyrefs = pTet->getVesProxyrefs();

    // Should be at least 2 vesicles
    AssertLog(vesproxyrefs.size() >= 2);

    // tetrahedron_global_id tet_gidx = pTet->idx();
    solver::spec_global_id spec1_gidx = def()->getSpec1gidx();
    solver::spec_global_id spec2_gidx = def()->getSpec2gidx();
    solver::vesicle_global_id ves1_idx = def()->getVes1idx();
    solver::vesicle_global_id ves2_idx = def()->getVes2idx();

    double max_distance = def()->max_distance();
    double min_distance = def()->min_distance();

    uint ngspecs = def()->countSpecsGlobal();
    uint ngbspecs = def()->countLinkSpecsGlobal();

    for (auto const& vesproxy1_mapit: vesproxyrefs) {
        auto const& vesproxy1 = vesproxy1_mapit.second;
        // Need to check if this is type Vesicle A
        if (vesproxy1->idx() == ves1_idx) {
            // Need to check any species dependencies
            bool allspecdepsgood = true;

            if (def()->gotVDep1()) {
                for (auto s: solver::spec_global_id::range(ngspecs)) {
                    uint vdep1 = def()->vdep1(s);
                    if (vdep1 == 0) {
                        continue;
                    }
                    // We got a dependency species. Need to see if it's on the vesicle
                    // surface or not.
                    uint cnt = vesproxy1->getSpecCount_V(s);
                    if (cnt >= vdep1) {
                        continue;
                    } else {
                        allspecdepsgood = false;
                        break;
                    }
                }
            }

            if (allspecdepsgood and def()->gotLDep1()) {
                for (auto s: solver::linkspec_global_id::range(ngbspecs)) {
                    uint ldep1 = def()->ldep1(s);
                    if (ldep1 == 0) {
                        continue;
                    }
                    uint cnt = vesproxy1->getLinkSpecCount_V(s);
                    if (cnt >= ldep1) {
                        continue;
                    } else {
                        allspecdepsgood = false;
                        break;
                    }
                }
            }

            // Go to next vesicle if this one doesn't have the right dependencies
            if (!allspecdepsgood) {
                continue;
            }
            // Not relative to vesicle centre
            auto spec1_positions = vesproxy1->getSurfaceSpecPos_absolute(spec1_gidx);

            for (uint spec1_idx_pos = 0; spec1_idx_pos < spec1_positions.size(); ++spec1_idx_pos) {
                math::position_abs spec1_pos = spec1_positions[spec1_idx_pos].second;

                for (auto const& vesproxy2_mapit: vesproxyrefs) {
                    auto const& vesproxy2 = vesproxy2_mapit.second;
                    if (vesproxy1 == vesproxy2) {
                        continue;
                    }
                    if (vesproxy2->idx() == ves2_idx) {
                        // Check dependencies first
                        if (def()->gotVDep2()) {
                            for (auto s: solver::spec_global_id::range(ngspecs)) {
                                uint vdep2 = def()->vdep2(s);
                                if (vdep2 == 0) {
                                    continue;
                                }
                                uint cnt = vesproxy2->getSpecCount_V(s);
                                if (cnt >= vdep2) {
                                    continue;
                                } else {
                                    allspecdepsgood = false;
                                    break;
                                }
                            }
                        }

                        if (allspecdepsgood && def()->gotLDep2()) {
                            for (auto s: solver::linkspec_global_id::range(ngbspecs)) {
                                uint ldep2 = def()->ldep2(s);
                                if (ldep2 == 0) {
                                    continue;
                                }
                                uint cnt = vesproxy2->getLinkSpecCount_V(s);
                                if (cnt >= ldep2) {
                                    continue;
                                } else {
                                    allspecdepsgood = false;
                                    break;
                                }
                            }
                        }

                        // Go to next vesicle if this one doesn't have the right
                        // dependencies
                        if (!allspecdepsgood) {
                            continue;
                        }
                        auto spec2_positions = vesproxy2->getSurfaceSpecPos_absolute(spec2_gidx);

                        for (uint spec2_idx_pos = 0; spec2_idx_pos < spec2_positions.size();
                             ++spec2_idx_pos) {
                            math::position_abs spec2_pos = spec2_positions[spec2_idx_pos].second;

                            // OK, these are absolute positions
                            double dist = distance(spec1_pos, spec2_pos);
                            if (dist <= max_distance && dist >= min_distance) {
                                // 1 Remove reactant species from the two vesicles
                                // 2 Create complex- add it to BOTH vesicles

                                // Need global species indices
                                solver::linkspec_global_id linkspec1_gidx =
                                    def()->getLinkSpec1gidx();
                                solver::linkspec_global_id linkspec2_gidx =
                                    def()->getLinkSpec2gidx();

                                // Remove specific spec, record positions relative to vesicle
                                // centre
                                vesproxy1->removeOneSurfSpec(spec1_gidx, spec1_idx_pos);
                                vesproxy2->removeOneSurfSpec(spec2_gidx, spec2_idx_pos);

                                // Link spec effectively has two positions. Both are relative to
                                // respective vesicle centres
                                // LinkSpec * link_spec = new
                                // LinkSpec(complexdef, tet_gidx, tet_gidx,
                                // spec1_pos_rel, spec2_pos_rel, vesproxy1, vesproxy2);
                                // vesproxy1->addLinkSpec(link_spec);
                                // vesproxy2->addLinkSpec(link_spec);
                                // ABOVE INFO NOW STORED IN TET AND MUST BE APPLIED AT NEXT
                                // UPDATE CHANCE
                                // pTet->addLinkSpec_proxy(complexdef, vesproxy1, vesproxy2,
                                // spec1_pos, spec2_pos);

                                // Add linkspec to a certain position
                                // These per vesicle indices need to be taken from a pool
                                // assigned to the vesproxy by comp
                                solver::linkspec_individual_id linkspec1_unique_idx =
                                    pTet->getNextLinkSpecUniqueIndex();
                                solver::linkspec_individual_id linkspec2_unique_idx =
                                    pTet->getNextLinkSpecUniqueIndex();

                                vesproxy1->addLinkSpec(linkspec1_unique_idx,
                                                       linkspec1_gidx,
                                                       spec1_pos);
                                vesproxy2->addLinkSpec(linkspec2_unique_idx,
                                                       linkspec2_gidx,
                                                       spec2_pos);

                                // Needs to be something like this to add to Comp so that pairs
                                // can be checked for vesicle unbinding
                                pTet->addNewLinkedSpecs(linkspec1_gidx,
                                                        linkspec2_gidx,
                                                        linkspec1_unique_idx,
                                                        linkspec2_unique_idx,
                                                        vesproxy1,
                                                        vesproxy2,
                                                        spec1_pos,
                                                        spec2_pos,
                                                        min_distance,
                                                        max_distance);

                                vesproxy1->updImmobility(def()->immobility());
                                vesproxy2->updImmobility(def()->immobility());

                                rExtent++;
                                return;
                            }
                        }
                    }
                }
            }
        }
    }

    // If we get here something went wrong- rate should be non-zero to try and
    // apply meaning vesicles and reactants should be available.
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

void VesBind::resetOccupancies() {
    pTet->resetPoolOccupancy();
}

}  // namespace steps::mpi::tetvesicle
