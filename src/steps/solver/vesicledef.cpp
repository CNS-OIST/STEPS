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

#include "solver/vesicledef.hpp"

// STEPS headers.
#include "math/constants.hpp"
#include "model/vesicle.hpp"
#include "model/vessurfsys.hpp"
#include "solver/compdef.hpp"
#include "solver/exocytosisdef.hpp"
#include "solver/patchdef.hpp"
#include "solver/statedef.hpp"
#include "solver/types.hpp"
#include "solver/vessdiffdef.hpp"
#include "solver/vessreacdef.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

namespace steps::solver {

Vesicledef::Vesicledef(Statedef* sd, vesicle_global_id idx, model::Vesicle* v)
    : pStatedef(sd)
    , pIdx(idx)
    , pDiameter(0)
    , pDcst(0)
    , pSetupRefsdone(false)
    , pSetupIndsdone(false)
    , pSpecsN_V(0)
    , pSpecsN_global(0)
    , pLinkSpecsN_global(0)
    , pVesSReacsN(0)
    , pExocytosisN(0)
    , pVesSDiffsN(0) {
    AssertLog(pStatedef != nullptr);
    AssertLog(v != nullptr);

    pName = v->getID();
    pDiameter = v->getDiameter();
    pDcst = v->getDcst();
    pVssys = v->getVesSurfsys();

    uint nvsreacs = pStatedef->countVesSReacs();
    pVesSReac_G2L.container().resize(nvsreacs);

    uint nexos = pStatedef->countExocytosis();
    pExocytosis_G2L.container().resize(nexos);

    uint nvsdiffs = pStatedef->countVesSDiffs();
    pVesSDiff_G2L.container().resize(nvsdiffs);

    pSpecsN_global = pStatedef->countSpecs();

    pSpec_G2L.container().resize(pSpecsN_global);

    pLinkSpecsN_global = pStatedef->countLinkSpecs();
}

////////////////////////////////////////////////////////////////////////////////

void Vesicledef::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pExocytosisKcst);
    util::checkpoint(cp_file, pExocytosisFlags);
}

////////////////////////////////////////////////////////////////////////////////

void Vesicledef::restore(std::fstream& cp_file) {
    util::restore(cp_file, pExocytosisKcst);
    util::restore(cp_file, pExocytosisFlags);
}

////////////////////////////////////////////////////////////////////////////////

double Vesicledef::vol() const noexcept {
    double radius = pDiameter / 2.0;
    return (4.0 / 3) * math::PI * radius * radius * radius;
}

////////////////////////////////////////////////////////////////////////////////

void Vesicledef::reset() {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);

    std::fill(pExocytosisFlags.begin(), pExocytosisFlags.end(), 0);

    for (auto i: exocytosis_local_id::range(countExocytosis())) {
        Exocytosisdef* exodef = exocytosisdef(i);
        pExocytosisKcst[i] = exodef->kcst();
    }
}

////////////////////////////////////////////////////////////////////////////////

void Vesicledef::setup_references() {
    AssertLog(pSetupRefsdone == false);
    AssertLog(pSetupIndsdone == false);

    const uint ngspecs = pStatedef->countSpecs();
    // const uint nglspecs = pStatedef->countLinkSpecs();
    const uint ngvsreacs = pStatedef->countVesSReacs();
    const uint ngvsdiffs = pStatedef->countVesSDiffs();
    const uint ngexos = pStatedef->countExocytosis();

    if (ngexos == 0) {
        AssertLog(pExocytosis_G2L.container().empty());
    }

    //  vesicles can appear in any compartment so any volume species have to be
    //  added to all comps

    // set up local sreac indices
    for (auto const& s: pVssys) {
        const auto& vessreacs = pStatedef->model()->getVesSurfsys(s)->_getAllVesSReacs();
        if (ngvsreacs == 0) {
            AssertLog(vessreacs.empty() == true);
        }

        for (auto const& vessr: vessreacs) {
            vessreac_global_id gidx = pStatedef->getVesSReacIdx((vessr.second));
            AssertLog(gidx < ngvsreacs);
            if (vessreacG2L(gidx).valid()) {
                continue;
            }
            pVesSReac_G2L[gidx] = vessreac_local_id(pVesSReacsN++);
        }

        const auto& sexos = pStatedef->model()->getVesSurfsys(s)->_getAllExocytosis();
        if (ngexos == 0) {
            AssertLog(sexos.empty() == true);
        }

        for (auto const& exo: sexos) {
            exocytosis_global_id gidx = pStatedef->getExocytosisIdx((exo.second));
            AssertLog(gidx < ngexos);
            if (exocytosisG2L(gidx).valid()) {
                continue;
            }
            pExocytosis_G2L[gidx] = exocytosis_local_id(pExocytosisN++);
        }

        const auto& vessdiffs = pStatedef->model()->getVesSurfsys(s)->_getAllVesSDiffs();
        if (ngvsdiffs == 0) {
            AssertLog(vessdiffs.empty() == true);
        }

        for (auto const& vessd: vessdiffs) {
            vessdiff_global_id gidx = pStatedef->getVesSDiffIdx((vessd.second));
            AssertLog(gidx < ngvsdiffs);
            if (vessurfdiffG2L(gidx).valid()) {
                continue;
            }
            pVesSDiff_G2L[gidx] = vessdiff_local_id(pVesSDiffsN++);
        }
    }

    for (auto sr: vessreac_global_id::range(ngvsreacs)) {
        if (vessreacG2L(sr).unknown()) {
            continue;
        }
        VesSReacdef* vessrdef = pStatedef->vessreacdef(sr);
        AssertLog(vessrdef != nullptr);

        for (auto s: spec_global_id::range(ngspecs)) {
            if (vessrdef->reqspec_V(s) == true) {
                AssertLog(pStatedef->specdef(s) != nullptr);
                if (specG2L(s).unknown()) {
                    pSpec_G2L[s] = spec_local_id(pSpecsN_V++);
                }
            }

            // Nothing to do for link species because we're using global ids

            if (vessrdef->reqspec_O(s) == true) {
                // Vesicle can exist in any compartment so add species to all
                // compartments
                for (auto const& c: pStatedef->comps()) {
                    c->addSpec(s);
                }
            }

            if (vessrdef->reqspec_S(s) == true) {
                // Difficult to predict which patches will be connected- potentially
                // all. Add to all
                for (auto const& p: pStatedef->patches()) {
                    p->addSpec(s);
                }
            }

            // Nothing to do for internal species because we'll just add all
            // somewhere??
        }
    }

    // Nothing TODO for exocytosis because it's very difficult to predict
    // which species will

    for (auto vsd: vessdiff_global_id::range(ngvsdiffs)) {
        if (vessurfdiffG2L(vsd).unknown()) {
            continue;
        }
        VesSDiffdef* vessddef = pStatedef->vessdiffdef(vsd);
        AssertLog(vessddef != nullptr);
        for (auto s: spec_global_id::range(ngspecs)) {
            if (vessddef->reqspec(s) == true) {
                if (specG2L(s).unknown()) {
                    pSpec_G2L[s] = spec_local_id(pSpecsN_V++);
                }
            }
        }
    }

    pSetupRefsdone = true;
}

////////////////////////////////////////////////////////////////////////////////

void Vesicledef::setup_indices() {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == false);

    // 1 -- DEAL WITH VESICLE SPECIES
    uint ngspecs = pStatedef->countSpecs();
    uint nglspecs = pStatedef->countLinkSpecs();

    if (countSpecs() != 0) {
        pSpec_L2G.container().resize(countSpecs());
        for (auto i: spec_global_id::range(ngspecs)) {
            spec_local_id lidx = specG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pSpec_L2G[lidx] = i;
        }
    }

    // 3 -- DEAL WITH VES SREACS
    if (countVesSReacs() != 0) {
        // Set up local indices.
        pVesSReac_L2G.container().resize(countVesSReacs());
        uint ngvessreacs = pStatedef->countVesSReacs();

        for (auto i: vessreac_global_id::range(ngvessreacs)) {
            vessreac_local_id lidx = vessreacG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pVesSReac_L2G[lidx] = i;
        }

        // Create _DEP, _LHS and _UPD vectors.
        uint arrsize_v = countSpecs() * countVesSReacs();
        pVesSReac_DEP_V_Spec.resize(arrsize_v);
        pVesSReac_LHS_V_Spec.resize(arrsize_v);
        pVesSReac_UPD_V_Spec.resize(arrsize_v);

        // Using global indices for B species
        uint arrsize_l = countLinkSpecs_V() * countVesSReacs();
        pVesSReac_DEP_L_Spec.resize(arrsize_l);
        pVesSReac_LHS_L_Spec.resize(arrsize_l);
        pVesSReac_UPD_L_Spec.resize(arrsize_l);

        // NOTE for compartment this has to be held as global indices since the same
        // vesicle can be present in any compartment in the system, in contrast to
        // surface reaction rules
        uint arrsize_o = countSpecs_O() * countVesSReacs();
        pVesSReac_DEP_O_Spec.resize(arrsize_o);
        pVesSReac_LHS_O_Spec.resize(arrsize_o);
        pVesSReac_UPD_O_Spec.resize(arrsize_o);

        uint arrsize_i = countSpecs_I() * countVesSReacs();
        pVesSReac_UPD_I_Spec.resize(arrsize_i);

        uint arrsize_s = countSpecs_S() * countVesSReacs();
        pVesSReac_DEP_S_Spec.resize(arrsize_s);
        pVesSReac_LHS_S_Spec.resize(arrsize_s);
        pVesSReac_UPD_S_Spec.resize(arrsize_s);

        uint arrsize_vdep = countSpecs_VDep() * countVesSReacs();
        pVesSReac_VDEP_Spec.resize(arrsize_vdep);

        // Fill the vectors with all kinds of useful information.
        for (auto ri: vessreac_local_id::range(countVesSReacs())) {
            VesSReacdef* vessrdef = vessreacdef(ri);

            // Handle surface stuff.
            for (auto si: spec_global_id::range(ngspecs)) {
                if (vessrdef->reqspec_V(si) == false) {
                    continue;
                }

                // TODO: turn into error check?
                spec_local_id sil = specG2L(si);
                AssertLog(sil.valid());

                uint aridx = _IDX_VesSReac_V_Spec(ri, sil);
                pVesSReac_DEP_V_Spec[aridx] = vessrdef->dep_V(si);
                pVesSReac_LHS_V_Spec[aridx] = vessrdef->lhs_V(si);
                pVesSReac_UPD_V_Spec[aridx] = vessrdef->upd_V(si);
            }

            for (auto si: linkspec_global_id::range(nglspecs)) {
                if (vessrdef->reqspec_L(si) == false) {
                    continue;
                }

                uint aridx = _IDX_VesSReac_L_Spec_global(ri, si);
                pVesSReac_DEP_L_Spec[aridx] = vessrdef->dep_L(si);
                pVesSReac_LHS_L_Spec[aridx] = vessrdef->lhs_L(si);
                pVesSReac_UPD_L_Spec[aridx] = vessrdef->upd_L(si);
            }

            // Handle the inside comp stuff.
            if (vessrdef->reqInside() == true) {
                for (auto si: spec_global_id::range(ngspecs)) {
                    if (vessrdef->reqspec_I(si) == false) {
                        continue;
                    }

                    uint aridx = _IDX_VesSReac_I_Spec_global(ri, si);
                    pVesSReac_UPD_I_Spec[aridx] = vessrdef->upd_I(si);
                }
            }

            // Handle the outside comp stuff.
            if (vessrdef->reqOutside() == true) {
                for (auto si: spec_global_id::range(ngspecs)) {
                    if (vessrdef->reqspec_O(si) == false) {
                        continue;
                    }

                    uint aridx = _IDX_VesSReac_O_Spec_global(ri, si);
                    pVesSReac_DEP_O_Spec[aridx] = vessrdef->dep_O(si);
                    pVesSReac_LHS_O_Spec[aridx] = vessrdef->lhs_O(si);
                    pVesSReac_UPD_O_Spec[aridx] = vessrdef->upd_O(si);
                }
            }

            // Handle the patch stuff.
            for (auto si: spec_global_id::range(ngspecs)) {
                if (vessrdef->reqspec_S(si) == false) {
                    continue;
                }

                uint aridx = _IDX_VesSReac_S_Spec_global(ri, si);
                pVesSReac_DEP_S_Spec[aridx] = vessrdef->dep_S(si);
                pVesSReac_LHS_S_Spec[aridx] = vessrdef->lhs_S(si);
                pVesSReac_UPD_S_Spec[aridx] = vessrdef->upd_S(si);
            }

            // Handle the vesicle surface dependency species stuff
            for (auto si: spec_global_id::range(ngspecs)) {
                uint aridx = _IDX_VesSReac_Vdep_Spec_global(ri, si);
                pVesSReac_VDEP_Spec[aridx] = vessrdef->vdep(si);
            }
        }
    }

    if (countExocytosis() != 0) {
        // Set up local indices.
        pExocytosis_L2G.container().resize(countExocytosis());
        uint ngexos = pStatedef->countExocytosis();

        for (auto i: exocytosis_global_id::range(ngexos)) {
            exocytosis_local_id lidx = exocytosisG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pExocytosis_L2G[lidx] = i;
        }
    }

    //  DEAL WITH SURFACE-DIFFUSION
    if (countVesSurfDiffs() != 0) {
        pVesSDiff_L2G.container().resize(countVesSurfDiffs());
        uint ngvessdiffs = pStatedef->countVesSDiffs();

        for (auto i: vessdiff_global_id::range(ngvessdiffs)) {
            vessdiff_local_id lidx = vessurfdiffG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pVesSDiff_L2G[lidx] = i;
        }

        uint arrsize = countSpecs() * countVesSurfDiffs();
        pVesSDiff_DEP_Spec.resize(arrsize);
        pVesSDiff_LIG.container().resize(countVesSurfDiffs());
        for (auto di: vessdiff_local_id::range(countVesSurfDiffs())) {
            VesSDiffdef* vessddef = vessurfdiffdef(di);
            pVesSDiff_LIG[di] = specG2L(vessddef->lig());
            for (auto si: spec_global_id::range(ngspecs)) {
                if (vessddef->reqspec(si) == false) {
                    continue;
                }
                spec_local_id sil = specG2L(si);
                AssertLog(sil.valid());
                uint aridx = _IDX_VesSDiff_Spec(di, sil);
                pVesSDiff_DEP_Spec[aridx] = vessddef->dep(si);
            }
        }
    }

    if (countExocytosis() != 0) {
        pExocytosisFlags.container().resize(countExocytosis());

        // Finally initialise Kcsts to user-supplied values
        pExocytosisKcst.container().resize(countExocytosis());
        for (auto i: exocytosis_local_id::range(countExocytosis())) {
            // exocytosisdef() returns global Reacdef by local index
            Exocytosisdef* exodef = exocytosisdef(i);
            pExocytosisKcst[i] = exodef->kcst();
        }
    }

    pSetupIndsdone = true;
}

////////////////////////////////////////////////////////////////////////////////

VesSReacdef* Vesicledef::vessreacdef(vessreac_local_id lidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < countVesSReacs());
    return pStatedef->vessreacdef(vessreacL2G(lidx));
}
////////////////////////////////////////////////////////////////////////////////

Exocytosisdef* Vesicledef::exocytosisdef(exocytosis_local_id lidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < countExocytosis());
    return pStatedef->exocytosisdef(exocytosisL2G(lidx));
}

////////////////////////////////////////////////////////////////////////////////

VesSDiffdef* Vesicledef::vessurfdiffdef(vessdiff_local_id dlidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(dlidx < countVesSurfDiffs());
    return pStatedef->vessdiffdef(vessurfdiffL2G(dlidx));
}

////////////////////////////////////////////////////////////////////////////////

void Vesicledef::setExoKcst(exocytosis_local_id exolidx, double kcst) {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    AssertLog(exolidx < countExocytosis());
    AssertLog(kcst >= 0.0);
    pExocytosisKcst[exolidx] = kcst;
}

////////////////////////////////////////////////////////////////////////////////

void Vesicledef::setExoActive(exocytosis_local_id exolidx, bool active) {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    if (active == true) {
        pExocytosisFlags.at(exolidx) &= ~INACTIVATED;
    } else {
        pExocytosisFlags.at(exolidx) |= INACTIVATED;
    }
}

////////////////////////////////////////////////////////////////////////////////

int Vesicledef::vessreac_dep_V(vessreac_local_id srlidx, spec_local_id splidx) const noexcept {
    return pVesSReac_DEP_V_Spec[splidx.get() + (srlidx.get() * countSpecs())];
}

////////////////////////////////////////////////////////////////////////////////

int Vesicledef::vessreac_dep_L(vessreac_local_id srlidx, linkspec_global_id spgidx) const noexcept {
    return pVesSReac_DEP_L_Spec[spgidx.get() + (srlidx.get() * countLinkSpecs_V())];
}

////////////////////////////////////////////////////////////////////////////////

int Vesicledef::vessreac_dep_O(vessreac_local_id srlidx, spec_global_id spgidx) const noexcept {
    return pVesSReac_DEP_O_Spec[spgidx.get() + (srlidx.get() * countSpecs_O())];
}

////////////////////////////////////////////////////////////////////////////////

int Vesicledef::vessreac_dep_S(vessreac_local_id srlidx, spec_global_id spgidx) const noexcept {
    return pVesSReac_DEP_S_Spec[spgidx.get() + (srlidx.get() * countSpecs_S())];
}

}  // namespace steps::solver
