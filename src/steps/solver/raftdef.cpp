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

#include "solver/raftdef.hpp"

// STEPS headers.
#include "model/raft.hpp"
#include "model/raftsys.hpp"
#include "solver/compdef.hpp"
#include "solver/exocytosisdef.hpp"
#include "solver/patchdef.hpp"
#include "solver/raftdisdef.hpp"
#include "solver/raftendocytosisdef.hpp"
#include "solver/raftsreacdef.hpp"
#include "solver/statedef.hpp"
#include "solver/types.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

#include "solver/fwd.hpp"

namespace steps::solver {

Raftdef::Raftdef(Statedef* sd, raft_global_id idx, model::Raft* r)
    : pStatedef(sd)
    , pIdx(idx)
    , pDiameter(0)
    , pDcst(0)
    , pSetupRefsdone(false)
    , pSetupIndsdone(false)
    , pSpecsN_Rs(0)
    , pSpecsN_global(0)
    , pRaftSReacsN(0)
    , pRaftEndocytosisN(0)
    , pRaftDisN(0) {
    AssertLog(pStatedef != nullptr);
    AssertLog(r != nullptr);

    pName = r->getID();
    pDiameter = r->getDiameter();
    pDcst = r->getDcst();
    pRssys = r->getRaftsys();

    uint nrsreacs = pStatedef->countRaftSReacs();
    pRaftSReac_G2L.resize(nrsreacs);

    uint nrendos = pStatedef->countRaftEndocytosis();
    pRaftEndocytosis_G2L.resize(nrendos);

    uint nrdiss = pStatedef->countRaftDiss();
    pRaftDis_G2L.resize(nrdiss);

    pSpecsN_global = pStatedef->countSpecs();

    pSpec_G2L.resize(pSpecsN_global);
}

////////////////////////////////////////////////////////////////////////////////

Raftdef::~Raftdef() = default;

////////////////////////////////////////////////////////////////////////////////

void Raftdef::checkpoint(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void Raftdef::restore(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void Raftdef::reset() const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
}

////////////////////////////////////////////////////////////////////////////////

void Raftdef::setup_references() {
    AssertLog(pSetupRefsdone == false);
    AssertLog(pSetupIndsdone == false);

    const uint ngspecs = pStatedef->countSpecs();
    const uint ngrsreacs = pStatedef->countRaftSReacs();
    const uint ngrendos = pStatedef->countRaftEndocytosis();
    const uint ngraftdiss = pStatedef->countRaftDiss();

    //  rafts can appear in any compartment so any volume species have to be added
    //  to all comps

    // set up local sreac indices
    for (auto const& s: pRssys) {
        const auto& raftsreacs = pStatedef->model()->getRaftsys(s)->_getAllRaftSReacs();
        if (ngrsreacs == 0) {
            AssertLog(raftsreacs.empty() == true);
        }

        for (auto const& raftsr: raftsreacs) {
            raftsreac_global_id gidx = pStatedef->getRaftSReacIdx(raftsr.second);
            AssertLog(gidx < ngrsreacs);
            if (raftsreacG2L(gidx).valid()) {
                continue;
            }
            pRaftSReac_G2L[gidx.get()] = raftsreac_local_id(pRaftSReacsN++);
        }

        const auto& rendos = pStatedef->model()->getRaftsys(s)->_getAllRaftEndocytosiss();
        if (ngrendos == 0) {
            AssertLog(rendos.empty() == true);
        }

        for (auto const& endo: rendos) {
            raftendocytosis_global_id gidx = pStatedef->getRaftEndocytosisIdx(endo.second);
            AssertLog(gidx < ngrendos);
            if (raftendocytosisG2L(gidx).valid()) {
                continue;
            }
            pRaftEndocytosis_G2L[gidx.get()] = raftendocytosis_local_id(pRaftEndocytosisN++);
        }

        const auto& rdiss = pStatedef->model()->getRaftsys(s)->_getAllRaftDiss();
        if (ngraftdiss == 0) {
            AssertLog(rdiss.empty() == true);
        }

        for (auto const& rdis: rdiss) {
            raftdis_global_id gidx = pStatedef->getRaftDisIdx(rdis.second);
            AssertLog(gidx < ngraftdiss);
            if (raftdisG2L(gidx).valid()) {
                continue;
            }
            pRaftDis_G2L[gidx.get()] = raftdis_local_id(pRaftDisN++);
        }
    }

    for (auto sr: raftsreac_global_id::range(ngrsreacs)) {
        if (raftsreacG2L(sr).unknown()) {
            continue;
        }
        RaftSReacdef* raftsrdef = pStatedef->raftsreacdef(sr);
        AssertLog(raftsrdef != nullptr);

        for (auto s: spec_global_id::range(ngspecs)) {
            if (raftsrdef->reqspec_Rs(s) == true) {
                AssertLog(pStatedef->specdef(s) != nullptr);
                if (specG2L(s).unknown()) {
                    pSpec_G2L[s.get()] = spec_local_id(pSpecsN_Rs++);
                }
            }

            if (raftsrdef->reqspec_O(s) == true) {
                // Raft can exist in any compartment so add species to all compartments
                for (auto const& c: pStatedef->comps()) {
                    c->addSpec(s);
                }
            }

            if (raftsrdef->reqspec_S(s) == true) {
                // Difficult to predict which patches will be connected- potentially
                // all. Add to all
                for (auto const& p: pStatedef->patches()) {
                    p->addSpec(s);
                }
            }

            if (raftsrdef->reqspec_I(s) == true) {
                // Raft can exist in any compartment so add species to all compartments
                for (auto const& c: pStatedef->comps()) {
                    c->addSpec(s);
                }
            }
        }
    }

    for (auto rdis: raftdis_global_id::range(ngraftdiss)) {
        if (raftdisG2L(rdis).unknown()) {
            continue;
        }
        RaftDisdef* rdisdef = pStatedef->raftdisdef(rdis);
        AssertLog(rdisdef != nullptr);
        for (auto s: spec_global_id::range(ngspecs)) {
            if (rdisdef->reqspec_S(s) == true) {
                AssertLog(pStatedef->specdef(s) != nullptr);
                if (specG2L(s).unknown()) {
                    pSpec_G2L[s.get()] = spec_local_id(pSpecsN_Rs++);
                }
            }
        }
    }

    // Need to do this for endocytosis??
    for (auto rendo: raftendocytosis_global_id::range(ngrendos)) {
        if (raftendocytosisG2L(rendo).unknown()) {
            continue;
        }
        RaftEndocytosisdef* rendodef = pStatedef->raftendocytosisdef(rendo);
        AssertLog(rendodef != nullptr);
        for (auto s: spec_global_id::range(ngspecs)) {
            if (rendodef->reqspec_S(s) == true) {
                AssertLog(pStatedef->specdef(s) != nullptr);
                if (specG2L(s).unknown()) {
                    pSpec_G2L[s.get()] = spec_local_id(pSpecsN_Rs++);
                }
            }
        }
    }

    pSetupRefsdone = true;
}

////////////////////////////////////////////////////////////////////////////////

void Raftdef::setup_indices() {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == false);

    // DEAL WITH RAFT SPECIES
    uint ngspecs = pStatedef->countSpecs();

    if (countSpecs() != 0) {
        pSpec_L2G.resize(countSpecs());
        for (auto i: spec_global_id::range(ngspecs)) {
            spec_local_id lidx = specG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pSpec_L2G[lidx.get()] = i;
        }
    }

    // DEAL WITH RAFT SREACS
    if (pRaftSReacsN != 0) {
        // Set up local indices.
        pRaftSReac_L2G.resize(pRaftSReacsN);
        uint ngraftsreacs = pStatedef->countRaftSReacs();

        for (auto i: raftsreac_global_id::range(ngraftsreacs)) {
            raftsreac_local_id lidx = raftsreacG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pRaftSReac_L2G[lidx.get()] = i;
        }

        // Create _DEP, _LHS and _UPD vectors.
        uint arrsize_rs = countSpecs() * pRaftSReacsN;
        pRaftSReac_DEP_Rs_Spec.resize(arrsize_rs);
        pRaftSReac_LHS_Rs_Spec.resize(arrsize_rs);
        pRaftSReac_UPD_Rs_Spec.resize(arrsize_rs);

        // NOTE for compartment this has to be held as global indices since the same
        // raft can be present in any compartment in the system, in contrast to
        // surface reaction rules
        uint arrsize_o = countSpecs_O() * pRaftSReacsN;
        pRaftSReac_DEP_O_Spec.resize(arrsize_o);
        pRaftSReac_LHS_O_Spec.resize(arrsize_o);
        pRaftSReac_UPD_O_Spec.resize(arrsize_o);

        uint arrsize_i = countSpecs_I() * pRaftSReacsN;
        pRaftSReac_DEP_I_Spec.resize(arrsize_i);
        pRaftSReac_LHS_I_Spec.resize(arrsize_i);
        pRaftSReac_UPD_I_Spec.resize(arrsize_i);

        uint arrsize_s = countSpecs_S() * pRaftSReacsN;
        pRaftSReac_DEP_S_Spec.resize(arrsize_s);
        pRaftSReac_LHS_S_Spec.resize(arrsize_s);
        pRaftSReac_UPD_S_Spec.resize(arrsize_s);

        uint arrsize_rsdep = ngspecs * pRaftSReacsN;
        pRaftSReac_RsDEP_Spec.resize(arrsize_rsdep);
        uint arrsize_antirsdep = ngspecs * pRaftSReacsN;
        pRaftSReac_AntiRsDEP_Spec.resize(arrsize_antirsdep);

        // Fill the vectors with all kinds of useful information.
        for (auto ri: raftsreac_local_id::range(pRaftSReacsN)) {
            RaftSReacdef* raftsrdef = raftsreacdef(ri);

            // Handle surface stuff.
            for (auto si: spec_global_id::range(ngspecs)) {
                if (raftsrdef->reqspec_Rs(si) == false) {
                    continue;
                }

                // TODO: turn into error check?
                spec_local_id sil = pSpec_G2L[si.get()];
                AssertLog(sil.valid());

                uint aridx = _IDX_RaftSReac_Rs_Spec(ri, sil);
                pRaftSReac_DEP_Rs_Spec[aridx] = raftsrdef->dep_Rs(si);
                pRaftSReac_LHS_Rs_Spec[aridx] = raftsrdef->lhs_Rs(si);
                pRaftSReac_UPD_Rs_Spec[aridx] = raftsrdef->upd_Rs(si);
            }

            // Handle the inside comp stuff.
            if (raftsrdef->reqInside() == true) {
                for (auto si: spec_global_id::range(ngspecs)) {
                    if (raftsrdef->reqspec_I(si) == false) {
                        continue;
                    }

                    uint aridx = _IDX_RaftSReac_I_Spec_global(ri, si);
                    pRaftSReac_DEP_I_Spec[aridx] = raftsrdef->dep_I(si);
                    pRaftSReac_LHS_I_Spec[aridx] = raftsrdef->lhs_I(si);
                    pRaftSReac_UPD_I_Spec[aridx] = raftsrdef->upd_I(si);
                }
            }

            // Handle the outside comp stuff.
            if (raftsrdef->reqOutside() == true) {
                for (auto si: spec_global_id::range(ngspecs)) {
                    if (raftsrdef->reqspec_O(si) == false) {
                        continue;
                    }

                    uint aridx = _IDX_RaftSReac_O_Spec_global(ri, si);
                    pRaftSReac_DEP_O_Spec[aridx] = raftsrdef->dep_O(si);
                    pRaftSReac_LHS_O_Spec[aridx] = raftsrdef->lhs_O(si);
                    pRaftSReac_UPD_O_Spec[aridx] = raftsrdef->upd_O(si);
                }
            }

            // Handle the patch stuff.
            for (auto si: spec_global_id::range(ngspecs)) {
                if (raftsrdef->reqspec_S(si) == false) {
                    continue;
                }

                uint aridx = _IDX_RaftSReac_S_Spec_global(ri, si);
                pRaftSReac_DEP_S_Spec[aridx] = raftsrdef->dep_S(si);
                pRaftSReac_LHS_S_Spec[aridx] = raftsrdef->lhs_S(si);
                pRaftSReac_UPD_S_Spec[aridx] = raftsrdef->upd_S(si);
            }

            // Handle the raft surface dependency species stuff
            for (auto si: spec_global_id::range(ngspecs)) {
                uint aridx = _IDX_RaftSReac_Rsdep_Spec_global(ri, si);
                pRaftSReac_RsDEP_Spec[aridx] = raftsrdef->rsdep(si);
            }
            // Handle the raft surface anti-dependency species stuff
            for (auto si: spec_global_id::range(ngspecs)) {
                uint aridx = _IDX_RaftSReac_AntiRsdep_Spec_global(ri, si);
                pRaftSReac_AntiRsDEP_Spec[aridx] = raftsrdef->anti_rsdep(si);
            }
        }
    }

    if (countRaftEndocytosis() != 0) {
        // Set up local indices.
        pRaftEndocytosis_L2G.resize(countRaftEndocytosis());
        uint ngrendos = pStatedef->countRaftEndocytosis();

        for (auto i: raftendocytosis_global_id::range(ngrendos)) {
            raftendocytosis_local_id lidx = raftendocytosisG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pRaftEndocytosis_L2G[lidx.get()] = i;
        }
    }

    if (countRaftDiss() != 0) {
        // Set up local indices.
        pRaftDis_L2G.resize(countRaftDiss());
        uint ngraftdiss = pStatedef->countRaftDiss();

        for (auto i: raftdis_global_id::range(ngraftdiss)) {
            raftdis_local_id lidx = raftdisG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pRaftDis_L2G[lidx.get()] = i;
        }
    }

    pSetupIndsdone = true;
}

////////////////////////////////////////////////////////////////////////////////

RaftSReacdef* Raftdef::raftsreacdef(raftsreac_local_id lidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < pRaftSReacsN);
    return pStatedef->raftsreacdef(raftsreacL2G(lidx));
}

////////////////////////////////////////////////////////////////////////////////

RaftEndocytosisdef* Raftdef::raftendocytosisdef(raftendocytosis_local_id lidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < countRaftEndocytosis());
    return pStatedef->raftendocytosisdef(raftendocytosisL2G(lidx));
}

////////////////////////////////////////////////////////////////////////////////

RaftDisdef* Raftdef::raftdisdef(raftdis_local_id lidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < countRaftDiss());
    return pStatedef->raftdisdef(raftdisL2G(lidx));
}

////////////////////////////////////////////////////////////////////////////////

int Raftdef::raftsreac_dep_Rs(raftsreac_local_id srlidx, spec_local_id splidx) const noexcept {
    return pRaftSReac_DEP_Rs_Spec[splidx.get() + (srlidx.get() * countSpecs())];
}

////////////////////////////////////////////////////////////////////////////////

int Raftdef::raftsreac_dep_O(raftsreac_local_id srlidx, spec_global_id spgidx) const noexcept {
    return pRaftSReac_DEP_O_Spec[spgidx.get() + (srlidx.get() * countSpecs_O())];
}

////////////////////////////////////////////////////////////////////////////////

int Raftdef::raftsreac_dep_I(raftsreac_local_id srlidx, spec_global_id spgidx) const noexcept {
    return pRaftSReac_DEP_I_Spec[spgidx.get() + (srlidx.get() * countSpecs_I())];
}

////////////////////////////////////////////////////////////////////////////////

int Raftdef::raftsreac_dep_S(raftsreac_local_id srlidx, spec_global_id spgidx) const noexcept {
    return pRaftSReac_DEP_S_Spec[spgidx.get() + (srlidx.get() * countSpecs_S())];
}

}  // namespace steps::solver
