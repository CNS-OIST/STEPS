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

// STL headers.
#include <sstream>
#include <string>

// STEPS headers.
#include "compdef.hpp"
#include "diffdef.hpp"
#include "ghkcurrdef.hpp"
#include "model/ghkcurr.hpp"
#include "model/ohmiccurr.hpp"
#include "model/sreac.hpp"
#include "model/vdepsreac.hpp"
#include "ohmiccurrdef.hpp"
#include "patchdef.hpp"
#include "solver/endocytosisdef.hpp"
#include "solver/raftgendef.hpp"
#include "sreacdef.hpp"
#include "types.hpp"
#include "vdepsreacdef.hpp"
// util
#include "util/checkpointing.hpp"
#include "util/error.hpp"
// logging
#include <easylogging++.h>

namespace steps::solver {

Patchdef::Patchdef(Statedef* sd, patch_global_id idx, wm::Patch* p)
    : pStatedef(sd)
    , pArea()
    , pIdx(idx)
    , pIcomp(nullptr)
    , pOcomp(nullptr)
    , pInner(nullptr)
    , pOuter(nullptr)
    , pSetupRefsdone(false)
    , pSetupIndsdone(false)
    , pSpecsN_I(0)
    , pSpecsN_S(0)
    , pSpecsN_O(0)
    , pSReacsN(0)
    , pSurfDiffsN(0)
    , pEndocytosisN(0)
    , pRaftGenN(0)
    , pVDepSReacsN(0)
    , pOhmicCurrsN(0)
    , pGHKcurrsN(0) {
    AssertLog(pStatedef != nullptr);
    AssertLog(p != nullptr);

    pName = p->getID();
    pArea = p->getArea();
    pPssys = p->getSurfsys();
    pIcomp = p->getIComp();
    pOcomp = p->getOComp();

    uint nspecs = pStatedef->countSpecs();
    pSpec_G2L.container().resize(nspecs);

    uint nsreacs = pStatedef->countSReacs();
    pSReac_G2L.container().resize(nsreacs);

    uint nsdiffs = pStatedef->countSurfDiffs();
    pSurfDiff_G2L.container().resize(nsdiffs);

    uint nendos = pStatedef->countEndocytosis();
    pEndocytosis_G2L.container().resize(nendos);

    uint nrgens = pStatedef->countRaftGens();
    pRaftGen_G2L.container().resize(nrgens);

    uint nohmiccurrs = pStatedef->countOhmicCurrs();
    pOhmicCurr_G2L.container().resize(nohmiccurrs);

    uint nghkcurrs = pStatedef->countGHKcurrs();
    pGHKcurr_G2L.container().resize(nghkcurrs);

    uint nvdepsreacs = pStatedef->countVDepSReacs();
    pVDepSReac_G2L.container().resize(nvdepsreacs);

    auto* tmp = dynamic_cast<tetmesh::TmPatch*>(p);
    if (tmp != nullptr) {
        for (auto* zone: tmp->getAllEndocyticZones()) {
            pEndocyticZonesdefs.push_back(new EndocyticZonedef(pStatedef, zone));
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

Patchdef::~Patchdef() {
    for (auto* zone: pEndocyticZonesdefs) {
        delete zone;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Patchdef::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pArea);
    util::checkpoint(cp_file, pPoolCount);
    util::checkpoint(cp_file, pPoolFlags);
    util::checkpoint(cp_file, pSReacKcst);
    util::checkpoint(cp_file, pSReacFlags);
    util::checkpoint(cp_file, pEndocytosisKcst);
    util::checkpoint(cp_file, pEndocytosisFlags);
    // pSurfDiffDcst not needed right now in contrast to Compdef pDiffDcst

    // Statedef does not own EndocyticZonesdefs
    for (auto const& ezone: pEndocyticZonesdefs) {
        ezone->checkpoint(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Patchdef::restore(std::fstream& cp_file) {
    util::restore(cp_file, pArea);
    util::restore(cp_file, pPoolCount);
    util::restore(cp_file, pPoolFlags);
    util::restore(cp_file, pSReacKcst);
    util::restore(cp_file, pSReacFlags);
    util::restore(cp_file, pEndocytosisKcst);
    util::restore(cp_file, pEndocytosisFlags);

    // Statedef does not own EndocyticZonesdefs
    for (auto const& ezone: pEndocyticZonesdefs) {
        ezone->restore(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Patchdef::setup_references() {
    AssertLog(pSetupRefsdone == false);
    AssertLog(pSetupIndsdone == false);

    // first find the inner and outer comps of this patch
    AssertLog(pIcomp != nullptr);
    comp_global_id icompidx = pStatedef->getCompIdx(pIcomp);
    pInner = pStatedef->compdef(icompidx);
    if (pOcomp != nullptr) {
        comp_global_id ocompidx = pStatedef->getCompIdx(pOcomp);
        pOuter = pStatedef->compdef(ocompidx);
    }

    const uint ngspecs = pStatedef->countSpecs();
    const uint ngsreacs = pStatedef->countSReacs();
    const uint ngsdiffs = pStatedef->countSurfDiffs();
    const uint ngendos = pStatedef->countEndocytosis();
    const uint ngraftgens = pStatedef->countRaftGens();
    const uint ngvdepsreacs = pStatedef->countVDepSReacs();
    const uint ngohmiccurrs = pStatedef->countOhmicCurrs();
    const uint ngghkcurrs = pStatedef->countGHKcurrs();

    if (ngspecs == 0) {
        AssertLog(pSpec_G2L.container().empty());
    }
    if (ngsreacs == 0) {
        AssertLog(pSReac_G2L.container().empty());
    }
    if (ngsdiffs == 0) {
        AssertLog(pSurfDiff_G2L.container().empty());
    }
    if (ngendos == 0) {
        AssertLog(pEndocytosis_G2L.container().empty());
    }
    if (ngraftgens == 0) {
        AssertLog(pRaftGen_G2L.container().empty());
    }
    if (ngvdepsreacs == 0) {
        AssertLog(pVDepSReac_G2L.container().empty());
    }
    if (ngohmiccurrs == 0) {
        AssertLog(pOhmicCurr_G2L.container().empty());
    }
    if (ngghkcurrs == 0) {
        AssertLog(pGHKcurr_G2L.container().empty());
    }

    // set up local sreac indices
    // auto s_end = pPssys.end();
    for (auto const& s: pPssys) {
        const auto& ssreacs = pStatedef->model()->getSurfsys(s)->_getAllSReacs();
        if (ngsreacs == 0) {
            AssertLog(ssreacs.empty() == true);
        }
        for (auto const& sr: ssreacs) {
            sreac_global_id gidx = pStatedef->getSReacIdx(sr.second);
            AssertLog(gidx < ngsreacs);
            if (sreacG2L(gidx).valid()) {
                continue;
            }
            pSReac_G2L[gidx] = sreac_local_id(pSReacsN++);
        }

        const auto& sdiffs = pStatedef->model()->getSurfsys(s)->_getAllDiffs();
        if (ngsdiffs == 0) {
            AssertLog(sdiffs.empty() == true);
        }
        for (auto const& sd: sdiffs) {
            surfdiff_global_id gidx = pStatedef->getSurfDiffIdx(sd.second);
            AssertLog(gidx < ngsdiffs);
            if (surfdiffG2L(gidx).valid()) {
                continue;
            }
            pSurfDiff_G2L[gidx] = surfdiff_local_id(pSurfDiffsN++);
        }

        const auto& sendos = pStatedef->model()->getSurfsys(s)->_getAllEndocytosis();
        if (ngendos == 0) {
            AssertLog(sendos.empty() == true);
        }
        for (auto const& endo: sendos) {
            endocytosis_global_id gidx = pStatedef->getEndocytosisIdx((endo.second));
            AssertLog(gidx < ngendos);
            if (endocytosisG2L(gidx).valid()) {
                continue;
            }
            pEndocytosis_G2L[gidx] = endocytosis_local_id(pEndocytosisN++);
        }

        const auto& rgens = pStatedef->model()->getSurfsys(s)->_getAllRaftGens();
        if (ngraftgens == 0) {
            AssertLog(rgens.empty() == true);
        }
        for (auto const& rgen: rgens) {
            raftgen_global_id gidx = pStatedef->getRaftGenIdx((rgen.second));
            AssertLog(gidx < ngraftgens);
            if (raftgenG2L(gidx).valid()) {
                continue;
            }
            pRaftGen_G2L[gidx] = raftgen_local_id(pRaftGenN++);
        }

        const auto& vdssreacs = pStatedef->model()->getSurfsys(s)->_getAllVDepSReacs();
        if (ngvdepsreacs == 0) {
            AssertLog(vdssreacs.empty() == true);
        }
        for (auto const& vdsr: vdssreacs) {
            vdepsreac_global_id gidx = pStatedef->getVDepSReacIdx(vdsr.second);
            AssertLog(gidx < ngvdepsreacs);
            if (vdepsreacG2L(gidx).valid()) {
                continue;
            }
            pVDepSReac_G2L[gidx] = vdepsreac_local_id(pVDepSReacsN++);
        }

        const auto& ocs = pStatedef->model()->getSurfsys(s)->_getAllOhmicCurrs();
        if (ngohmiccurrs == 0) {
            AssertLog(ocs.empty() == true);
        }
        for (auto const& oc: ocs) {
            ohmiccurr_global_id gidx = pStatedef->getOhmicCurrIdx(oc.second);
            AssertLog(gidx < ngohmiccurrs);
            if (ohmiccurrG2L(gidx).valid()) {
                continue;
            }
            pOhmicCurr_G2L[gidx] = ohmiccurr_local_id(pOhmicCurrsN++);
        }

        const auto& ghks = pStatedef->model()->getSurfsys(s)->_getAllGHKcurrs();
        if (ngghkcurrs == 0) {
            AssertLog(ghks.empty() == true);
        }
        for (auto const& ghk: ghks) {
            ghkcurr_global_id gidx = pStatedef->getGHKcurrIdx(ghk.second);
            AssertLog(gidx < ngghkcurrs);
            if (ghkcurrG2L(gidx).valid()) {
                continue;
            }
            pGHKcurr_G2L[gidx] = ghkcurr_local_id(pGHKcurrsN++);
        }
    }

    // Now add all species that appear in all surface reactions, ohmic currents,
    // ghk currents and voltage-dependent transitions/reactions that can occur
    // on this patch: to the patch, inner or outer compartment.
    for (auto sr: sreac_global_id::range(ngsreacs)) {
        if (sreacG2L(sr).unknown()) {
            continue;
        }

        SReacdef* srdef = pStatedef->sreacdef(sr);
        AssertLog(srdef != nullptr);
        for (auto s: spec_global_id::range(ngspecs)) {
            if (srdef->reqspec_S(s) == true) {
                AssertLog(pStatedef->specdef(s) != nullptr);
                if (specG2L(s).unknown()) {
                    pSpec_G2L[s] = spec_local_id(pSpecsN_S++);
                }
            }
            if (srdef->reqspec_I(s) == true) {
                AssertLog(pInner != nullptr);
                pInner->addSpec(s);
            }
            if (srdef->reqspec_O(s) == true) {
                if (pOuter == nullptr) {
                    std::ostringstream os;
                    os << "Can't add surface reaction '" << srdef->name() << "' to patch '";
                    os << name() << "'. Outer compartment not defined for this patch.";
                    ArgErrLog(os.str());
                }
                pOuter->addSpec(s);
            }
        }
    }

    for (auto sd: surfdiff_global_id::range(ngsdiffs)) {
        if (surfdiffG2L(sd).unknown()) {
            continue;
        }
        SurfDiffdef* sddef = pStatedef->surfdiffdef(sd);
        AssertLog(sddef != nullptr);
        for (auto s: spec_global_id::range(ngspecs)) {
            if (sddef->reqspec(s) == true) {
                if (specG2L(s).unknown()) {
                    pSpec_G2L[s] = spec_local_id(pSpecsN_S++);
                }
            }
        }
    }

    for (auto endo: endocytosis_global_id::range(ngendos)) {
        if (endocytosisG2L(endo).unknown()) {
            continue;
        }
        Endocytosisdef* endodef = pStatedef->endocytosisdef(endo);
        AssertLog(endodef != nullptr);
        for (auto s: spec_global_id::range(ngspecs)) {
            if (endodef->reqspec_S(s) == true) {
                AssertLog(pStatedef->specdef(s) != nullptr);
                if (specG2L(s).unknown()) {
                    pSpec_G2L[s] = spec_local_id(pSpecsN_S++);
                }
            }
        }
    }

    for (auto rgen: raftgen_global_id::range(ngraftgens)) {
        if (raftgenG2L(rgen).unknown()) {
            continue;
        }
        RaftGendef* rgendef = pStatedef->raftgendef(rgen);
        AssertLog(rgendef != nullptr);
        for (auto s: spec_global_id::range(ngspecs)) {
            if (rgendef->reqspec_S(s) == true) {
                AssertLog(pStatedef->specdef(s) != nullptr);
                if (specG2L(s).unknown()) {
                    pSpec_G2L[s] = spec_local_id(pSpecsN_S++);
                }
            }
        }
    }

    for (auto vdsr: vdepsreac_global_id::range(ngvdepsreacs)) {
        if (vdepsreacG2L(vdsr).unknown()) {
            continue;
        }
        VDepSReacdef* vdsrdef = pStatedef->vdepsreacdef(vdsr);
        AssertLog(vdsrdef != nullptr);
        for (auto s: spec_global_id::range(ngspecs)) {
            if (vdsrdef->reqspec_S(s) == true) {
                AssertLog(pStatedef->specdef(s) != nullptr);
                if (specG2L(s).unknown()) {
                    pSpec_G2L[s] = spec_local_id(pSpecsN_S++);
                }
            }
            if (vdsrdef->reqspec_I(s) == true) {
                AssertLog(pInner != nullptr);
                pInner->addSpec(s);
            }
            if (vdsrdef->reqspec_O(s) == true) {
                if (pOuter == nullptr) {
                    std::ostringstream os;
                    os << "Can't add voltage-dependent reaction '" << vdsrdef->name()
                       << "' to patch '";
                    os << name() << "'. Outer compartment not defined for this patch.";
                    ArgErrLog(os.str());
                }
                pOuter->addSpec(s);
            }
        }
    }

    for (auto oc: ohmiccurr_global_id::range(ngohmiccurrs)) {
        if (ohmiccurrG2L(oc).unknown()) {
            continue;
        }
        OhmicCurrdef* ocdef = pStatedef->ohmiccurrdef(oc);
        AssertLog(ocdef != nullptr);
        uint added = 0;
        for (auto s: spec_global_id::range(ngspecs)) {
            // Add the channel state
            if (ocdef->req(s) == true) {
                AssertLog(pStatedef->specdef(s) != nullptr);
                if (specG2L(s).unknown()) {
                    pSpec_G2L[s] = spec_local_id(pSpecsN_S++);
                }
                added += 1;
            }
        }
        // Only one channel state should be added per ohmic current
        AssertLog(added == 1);
    }

    for (auto ghk: ghkcurr_global_id::range(ngghkcurrs)) {
        if (ghkcurrG2L(ghk).unknown()) {
            continue;
        }
        GHKcurrdef* ghkdef = pStatedef->ghkcurrdef(ghk);
        AssertLog(ghkdef != nullptr);
        uint added = 0;
        for (auto s: spec_global_id::range(ngspecs)) {
            // Add the channel state
            if (ghkdef->req(s) == true) {
                AssertLog(pStatedef->specdef(s) != nullptr);
                // Only add the channel state, not the volume ion species (that affects
                // the GHK rate)
                if (ghkdef->req_v(s) == false) {
                    if (specG2L(s).unknown()) {
                        pSpec_G2L[s] = spec_local_id(pSpecsN_S++);
                    }
                    added += 1;
                }
            }
            // Add the volume ion species to the inner and outer compartment.
            if (ghkdef->req_v(s) == true) {
                AssertLog(pInner != nullptr);
                pInner->addSpec(s);
                if (pOuter == nullptr) {
                    if (ghkdef->voconc() < 0.0) {
                        std::ostringstream os;
                        os << "Can't add GHK current '" << ghkdef->name() << "' to patch '";
                        os << name() << "'. Outer compartment not defined for this patch ";
                        os << "and no virtual concentration has been defined.";
                        ArgErrLog(os.str());
                    }
                } else if (ghkdef->voconc() < 0.0) {
                    pOuter->addSpec(s);
                }
            }
        }
        // Only one channel state should be added per ghk current
        AssertLog(added == 1);
    }

    pSetupRefsdone = true;
}

////////////////////////////////////////////////////////////////////////////////

void Patchdef::addSpec(spec_global_id gidx) {
    AssertLog(pSetupIndsdone == false);
    AssertLog(pStatedef->specdef(gidx) != nullptr);
    if (specG2L(gidx).valid()) {
        return;
    }
    pSpec_G2L[gidx] = spec_local_id(pSpecsN_S++);
}

////////////////////////////////////////////////////////////////////////////////

void Patchdef::setup_indices() {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == false);

    // 1 -- DEAL WITH PATCH SPECIES
    // (Only if any species have been added to the patch)
    //   -> Setup local indices for all species.
    //
    // 2 -- COPY #SPECIES FOR INNER AND OUTER COMPS
    //   (These are required a lot during simulation, so it's sound
    //   to have them ready here to avoid an extra level of pointer
    //   lookup.)
    //
    // 3 -- DEAL WITH PATCH SREAC'S
    // (Only if any surface reactions have been added to the patch)
    //   -> Setup local indices for all surface reactions.
    //   -> The SReac objects have LHS, DEP and UPD vectors expressed in
    //      global species indices. Transform this to local indices:
    //      -> Pre-create pSReac_DEP, _LHS and _UPD vectors
    //          -> Always for surface (_S)
    //          -> For inner comp (_I)
    //          -> For outer comp (_O) if outer comp has been defined
    //      -> Loop over the SReacDef objects added to this patch:
    //          -> Fill out the newly created vectors by appropriate
    //             copying of the vectors defined the SReacDef object.
    //          -> While doing this, check whether everything can be
    //             resolved.
    //

    // 4 -- DEAL WITH PATCH SURFACE-DIFFUSION

    // 5 -- DEAL WITH ENDOCYTOSIS

    // 6 -- DEAL WITH RAFT GEN

    // 7 -- DEAL WITH PATCH VOLTAGE-DEPENDENT REACTIONS
    // (Only if any vdep reactions have been added to the patch)
    //   -> Setup local indices for all vdep reactions.
    //   -> The VDepSReac objects have LHS, DEP and UPD vectors expressed in
    //      global species indices. Transform this to local indices:
    //      -> Pre-create pVDepSReac_DEP, _LHS and _UPD vectors
    //          -> Always for surface (_S)
    //          -> For inner comp (_I)
    //          -> For outer comp (_O) if outer comp has been defined
    //      -> Loop over the VDepSReacDef objects added to this patch:
    //          -> Fill out the newly created vectors by appropriate
    //             copying of the vectors defined the VDepSReacDef object.
    //          -> While doing this, check whether everything can be
    //             resolved.
    // 5 -- DEAL WITH OHMIC CURRENTS
    // 6 -- DEAL WITH GHK CURRENTS
    // 7 -- DEAL WITH V-DEPENDENT TRANSITIONS

    // 1 -- DEAL WITH PATCH SPECIES
    uint ngspecs = pStatedef->countSpecs();
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

    // 2 -- COPY #SPECS FOR INNER AND OUTER COMPS
    if (pInner != nullptr) {
        pSpecsN_I = pInner->countSpecs();
    }
    if (pOuter != nullptr) {
        pSpecsN_O = pOuter->countSpecs();
    }

    // 3 -- DEAL WITH PATCH SREAC'S
    if (countSReacs() != 0) {
        // Set up local indices.
        pSReac_L2G.container().resize(countSReacs());
        uint ngsreacs = pStatedef->countSReacs();
        for (auto i: sreac_global_id::range(ngsreacs)) {
            sreac_local_id lidx = sreacG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pSReac_L2G[lidx] = i;
        }

        // Create _DEP, _LHS and _UPD vectors.
        uint arrsize_s = countSpecs() * countSReacs();
        pSReac_DEP_S_Spec.resize(arrsize_s);
        pSReac_LHS_S_Spec.resize(arrsize_s);
        pSReac_UPD_S_Spec.resize(arrsize_s);

        AssertLog(pInner != nullptr);  // Inner comp should exist
        {
            uint arrsize_i = countSpecs_I() * countSReacs();
            pSReac_DEP_I_Spec.resize(arrsize_i);
            pSReac_LHS_I_Spec.resize(arrsize_i);
            pSReac_UPD_I_Spec.resize(arrsize_i);
        }
        if (pOuter != nullptr)  // Only create if outer comp exists.
        {
            uint arrsize_o = countSpecs_O() * countSReacs();
            pSReac_DEP_O_Spec.resize(arrsize_o);
            pSReac_LHS_O_Spec.resize(arrsize_o);
            pSReac_UPD_O_Spec.resize(arrsize_o);
        }

        // Fill the vectors with all kinds of useful information.
        for (auto ri: sreac_local_id::range(countSReacs())) {
            SReacdef* srdef = sreacdef(ri);

            // Handle surface stuff.
            for (auto si: spec_global_id::range(ngspecs)) {
                if (srdef->reqspec_S(si) == false) {
                    continue;
                }

                // TODO: turn into error check?
                spec_local_id sil = pSpec_G2L[si];
                AssertLog(sil.valid());

                uint aridx = _IDX_SReac_S_Spec(ri, sil);
                pSReac_DEP_S_Spec[aridx] = srdef->dep_S(si);
                pSReac_LHS_S_Spec[aridx] = srdef->lhs_S(si);
                pSReac_UPD_S_Spec[aridx] = srdef->upd_S(si);
            }

            // Handle the inside comp stuff.
            if (srdef->reqInside() == true) {
                // TODO: turn into real error check?
                AssertLog(pInner != nullptr);

                for (auto si: spec_global_id::range(ngspecs)) {
                    if (srdef->reqspec_I(si) == false) {
                        continue;
                    }

                    // TODO: turn into error check?
                    spec_local_id sil = specG2L_I(si);
                    AssertLog(sil.valid());

                    uint aridx = _IDX_SReac_I_Spec(ri, sil);
                    pSReac_DEP_I_Spec[aridx] = srdef->dep_I(si);
                    pSReac_LHS_I_Spec[aridx] = srdef->lhs_I(si);
                    pSReac_UPD_I_Spec[aridx] = srdef->upd_I(si);
                }
            }

            // Handle the outside comp stuff.
            if (srdef->reqOutside() == true) {
                // TODO: turn into real error check?
                AssertLog(pOuter != nullptr);

                for (auto si: spec_global_id::range(ngspecs)) {
                    if (srdef->reqspec_O(si) == false) {
                        continue;
                    }

                    // TODO: turn into error check?
                    spec_local_id sil = specG2L_O(si);
                    AssertLog(sil.valid());

                    uint aridx = _IDX_SReac_O_Spec(ri, sil);
                    pSReac_DEP_O_Spec[aridx] = srdef->dep_O(si);
                    pSReac_LHS_O_Spec[aridx] = srdef->lhs_O(si);
                    pSReac_UPD_O_Spec[aridx] = srdef->upd_O(si);
                }
            }
        }
    }

    // 4 -- DEAL WITH PATCH SURFACE-DIFFUSION
    if (countSurfDiffs() != 0) {
        pSurfDiff_L2G.container().resize(countSurfDiffs());
        uint ngsdiffs = pStatedef->countSurfDiffs();

        for (auto i: surfdiff_global_id::range(ngsdiffs)) {
            surfdiff_local_id lidx = surfdiffG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pSurfDiff_L2G[lidx] = i;
        }

        uint arrsize = countSpecs() * countSurfDiffs();
        pSurfDiff_DEP_Spec.resize(arrsize);
        pSurfDiff_LIG.container().resize(countSurfDiffs());
        for (auto di: surfdiff_local_id::range(countSurfDiffs())) {
            SurfDiffdef* sddef = surfdiffdef(di);
            pSurfDiff_LIG[di] = specG2L(sddef->lig());
            for (auto si: spec_global_id::range(ngspecs)) {
                if (sddef->reqspec(si) == false) {
                    continue;
                }
                spec_local_id sil = specG2L(si);
                AssertLog(sil.valid());
                uint aridx = _IDX_SurfDiff_Spec(di, sil);
                pSurfDiff_DEP_Spec[aridx] = sddef->dep(si);
            }
        }
    }

    // 5 -- DEAL WITH ENDOCYTOSIS
    if (countEndocytosis() != 0) {
        // Set up local indices.
        pEndocytosis_L2G.container().resize(countEndocytosis());
        uint ngendos = pStatedef->countEndocytosis();
        for (auto i: endocytosis_global_id::range(ngendos)) {
            endocytosis_local_id lidx = endocytosisG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pEndocytosis_L2G[lidx] = i;
        }
    }

    // 6 -- DEAL WITH RAFT GEN
    if (countRaftGens() != 0) {
        // Set up local indices.
        pRaftGen_L2G.container().resize(countRaftGens());
        uint ngraftgens = pStatedef->countRaftGens();
        for (auto i: raftgen_global_id::range(ngraftgens)) {
            raftgen_local_id lidx = raftgenG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pRaftGen_L2G[lidx] = i;
        }
    }

    // 7 -- DEAL WITH PATCH VOLTAGE-DEPENDENT SURFACE REACTIONS
    if (countVDepSReacs() != 0) {
        // Set up local indices.
        pVDepSReac_L2G.container().resize(countVDepSReacs());
        uint ngvdsreacs = pStatedef->countVDepSReacs();
        for (auto i: vdepsreac_global_id::range(ngvdsreacs)) {
            vdepsreac_local_id lidx = vdepsreacG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pVDepSReac_L2G[lidx] = i;
        }

        // Create _DEP, _LHS and _UPD vectors.
        uint arrsize_s = countSpecs() * countVDepSReacs();
        pVDepSReac_DEP_S_Spec.resize(arrsize_s);
        pVDepSReac_LHS_S_Spec.resize(arrsize_s);
        pVDepSReac_UPD_S_Spec.resize(arrsize_s);

        AssertLog(pInner != nullptr);  // Inner comp should exist
        {
            uint arrsize_i = countSpecs_I() * countVDepSReacs();
            pVDepSReac_DEP_I_Spec.resize(arrsize_i);
            pVDepSReac_LHS_I_Spec.resize(arrsize_i);
            pVDepSReac_UPD_I_Spec.resize(arrsize_i);
        }
        if (pOuter != nullptr)  // Only create if outer comp exists.
        {
            uint arrsize_o = countSpecs_O() * countVDepSReacs();
            pVDepSReac_DEP_O_Spec.resize(arrsize_o);
            pVDepSReac_LHS_O_Spec.resize(arrsize_o);
            pVDepSReac_UPD_O_Spec.resize(arrsize_o);
        }

        // Fill the vectors with all kinds of useful information.
        for (auto ri: vdepsreac_local_id::range(countVDepSReacs())) {
            VDepSReacdef* vdsrdef = vdepsreacdef(ri);

            // Handle surface stuff.
            for (auto si: spec_global_id::range(ngspecs)) {
                if (vdsrdef->reqspec_S(si) == false) {
                    continue;
                }

                // TODO: turn into error check?
                spec_local_id sil = specG2L(si);
                AssertLog(sil.valid());

                uint aridx = _IDX_VDepSReac_S_Spec(ri, sil);
                pVDepSReac_DEP_S_Spec[aridx] = vdsrdef->dep_S(si);
                pVDepSReac_LHS_S_Spec[aridx] = vdsrdef->lhs_S(si);
                pVDepSReac_UPD_S_Spec[aridx] = vdsrdef->upd_S(si);
            }

            // Handle the inside comp stuff.
            if (vdsrdef->reqInside() == true) {
                // TODO: turn into real error check?
                AssertLog(pInner != nullptr);

                for (auto si: spec_global_id::range(ngspecs)) {
                    if (vdsrdef->reqspec_I(si) == false) {
                        continue;
                    }

                    // TODO: turn into error check?
                    spec_local_id sil = specG2L_I(si);
                    AssertLog(sil.valid());

                    uint aridx = _IDX_VDepSReac_I_Spec(ri, sil);
                    pVDepSReac_DEP_I_Spec[aridx] = vdsrdef->dep_I(si);
                    pVDepSReac_LHS_I_Spec[aridx] = vdsrdef->lhs_I(si);
                    pVDepSReac_UPD_I_Spec[aridx] = vdsrdef->upd_I(si);
                }
            }

            // Handle the outside comp stuff.
            if (vdsrdef->reqOutside() == true) {
                // TODO: turn into real error check?
                AssertLog(pOuter != nullptr);

                for (auto si: spec_global_id::range(ngspecs)) {
                    if (vdsrdef->reqspec_O(si) == false) {
                        continue;
                    }

                    // TODO: turn into error check?
                    spec_local_id sil = specG2L_O(si);
                    AssertLog(sil.valid());

                    uint aridx = _IDX_VDepSReac_O_Spec(ri, sil);
                    pVDepSReac_DEP_O_Spec[aridx] = vdsrdef->dep_O(si);
                    pVDepSReac_LHS_O_Spec[aridx] = vdsrdef->lhs_O(si);
                    pVDepSReac_UPD_O_Spec[aridx] = vdsrdef->upd_O(si);
                }
            }
        }
    }

    // 8 -- DEAL WITH OHMIC CURRENTS
    if (countOhmicCurrs() != 0) {
        // Set up local indices.
        pOhmicCurr_L2G.container().resize(countOhmicCurrs());
        uint ngohmiccurrs = pStatedef->countOhmicCurrs();
        for (auto i: ohmiccurr_global_id::range(ngohmiccurrs)) {
            ohmiccurr_local_id lidx = ohmiccurrG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pOhmicCurr_L2G[lidx] = i;
        }

        // Create local _DEP and _CHANSTATE vectors
        uint arrsize1 = countSpecs() * countOhmicCurrs();
        uint arrsize2 = countOhmicCurrs();
        pOhmicCurr_DEP_Spec.resize(arrsize1);
        pOhmicCurr_CHANSTATE.container().resize(arrsize2);

        // Fill the vectors with useful information
        for (auto ri: ohmiccurr_local_id::range(countOhmicCurrs())) {
            OhmicCurrdef* ocdef = ohmiccurrdef(ri);
            for (auto si: spec_global_id::range(ngspecs)) {
                if (ocdef->req(si) == false) {
                    continue;
                }

                spec_local_id slidx = specG2L(si);
                AssertLog(slidx.valid());

                uint aridx = _IDX_OhmicCurr_Spec(ri, slidx);
                pOhmicCurr_DEP_Spec[aridx] = ocdef->dep(si);
            }
            pOhmicCurr_CHANSTATE[ri] = specG2L(ocdef->chanstate());
        }
    }

    // 9 -- DEAL WITH GHK CURRENTS
    if (countGHKcurrs() != 0) {
        // Set up local indices.
        pGHKcurr_L2G.container().resize(countGHKcurrs());
        uint ngghkcurrs = pStatedef->countGHKcurrs();
        for (auto i: ghkcurr_global_id::range(ngghkcurrs)) {
            ghkcurr_local_id lidx = ghkcurrG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pGHKcurr_L2G[lidx] = i;
        }
        // Create local _DEP and _CHANSTATE vectors.
        uint arrsize1 = countSpecs() * countGHKcurrs();
        uint arrsize2 = countGHKcurrs();
        pGHKcurr_DEP_Spec.resize(arrsize1);
        pGHKcurr_CHANSTATE.container().resize(arrsize2);
        pGHKcurr_ION.container().resize(arrsize2);

        // Fill the vectors with useful information
        for (auto ri: ghkcurr_local_id::range(countGHKcurrs())) {
            GHKcurrdef* ghkdef = ghkcurrdef(ri);
            for (auto si: spec_global_id::range(ngspecs)) {
                if (ghkdef->req(si) == false) {
                    continue;
                }

                spec_local_id slidx = specG2L(si);
                // If not the volume ion species local index should be defined.
                if (ghkdef->req_v(si) == false) {
                    AssertLog(slidx.valid());
                }

                uint aridx = _IDX_GHKcurr_Spec(ri, slidx);
                // DEP information can be rate (value:2) not just stoichiometry
                // (value:1)
                pGHKcurr_DEP_Spec[aridx] = ghkdef->dep(si);
            }
            pGHKcurr_CHANSTATE[ri] = specG2L(ghkdef->chanstate());
            pGHKcurr_ION[ri] = specG2L(ghkdef->ion());
        }
    }

    // Initialise the pools and flags members to zeros.
    if (countSpecs() != 0) {
        pPoolCount.container().resize(countSpecs());
        pPoolFlags.container().resize(countSpecs());
    }

    if (countSReacs() != 0) {
        pSReacFlags.container().resize(countSReacs());

        // Finally initialise Kcsts to user-supplied values
        pSReacKcst.container().resize(countSReacs());

        for (auto i: sreac_local_id::range(countSReacs())) {
            // sreacdef() returns global Reacdef by local index
            SReacdef* sreac = sreacdef(i);
            pSReacKcst[i] = sreac->kcst();
        }
    }

    if (countSurfDiffs() != 0) {
        pSurfDiffDcst.container().resize(countSurfDiffs());

        for (auto i: surfdiff_local_id::range(countSurfDiffs())) {
            // sdiffdef() returns global SDiffdef by local index
            SurfDiffdef* sdiff = surfdiffdef(i);
            pSurfDiffDcst[i] = sdiff->dcst();
        }
    }

    if (countEndocytosis() != 0) {
        pEndocytosisFlags.container().resize(countEndocytosis());

        // Finally initialise Kcsts to user-supplied values
        pEndocytosisKcst.container().resize(countEndocytosis());
        for (auto i: endocytosis_local_id::range(countEndocytosis())) {
            // endocytosisdef() returns global Reacdef by local index
            Endocytosisdef* endodef = endocytosisdef(i);
            pEndocytosisKcst[i] = endodef->kcst();
        }
    }

    pSetupIndsdone = true;
}

////////////////////////////////////////////////////////////////////////////////

void Patchdef::setArea(double a) {
    AssertLog(a > 0.0);
    pArea = a;
}

////////////////////////////////////////////////////////////////////////////////

void Patchdef::reset() {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);

    std::fill(pPoolCount.begin(), pPoolCount.end(), 0.0);
    std::fill(pPoolFlags.begin(), pPoolFlags.end(), 0);
    std::fill(pSReacFlags.begin(), pSReacFlags.end(), 0);
    std::fill(pEndocytosisFlags.begin(), pEndocytosisFlags.end(), 0);

    for (auto i: sreac_local_id::range(countSReacs())) {
        // sreacdef() returns global Reacdef by local index
        SReacdef* sreac = sreacdef(i);
        pSReacKcst[i] = sreac->kcst();
    }

    for (auto i: surfdiff_local_id::range(countSurfDiffs())) {
        // sdiffdef() returns global SDiffdef by local index
        SurfDiffdef* sdiff = surfdiffdef(i);
        pSurfDiffDcst[i] = sdiff->dcst();
    }

    for (auto i: endocytosis_local_id::range(countEndocytosis())) {
        Endocytosisdef* endodef = endocytosisdef(i);
        pEndocytosisKcst[i] = endodef->kcst();
    }
}

////////////////////////////////////////////////////////////////////////////////

void Patchdef::setCount(spec_local_id slidx, double count) {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    AssertLog(slidx < countSpecs());
    AssertLog(count >= 0.0);
    pPoolCount[slidx] = count;
}

////////////////////////////////////////////////////////////////////////////////

void Patchdef::setClamped(spec_local_id slidx, bool clamp) {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    AssertLog(slidx < countSpecs());
    if (clamp == true) {
        pPoolFlags[slidx] |= CLAMPED;
    } else {
        pPoolFlags[slidx] &= ~CLAMPED;
    }
}

////////////////////////////////////////////////////////////////////////////////

SReacdef* Patchdef::sreacdef(sreac_local_id lidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < countSReacs());
    return pStatedef->sreacdef(sreacL2G(lidx));
}
////////////////////////////////////////////////////////////////////////////////

Endocytosisdef* Patchdef::endocytosisdef(endocytosis_local_id lidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < countEndocytosis());
    return pStatedef->endocytosisdef(endocytosisL2G(lidx));
}

////////////////////////////////////////////////////////////////////////////////

RaftGendef* Patchdef::raftgendef(raftgen_local_id lidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < countRaftGens());
    return pStatedef->raftgendef(raftgenL2G(lidx));
}

////////////////////////////////////////////////////////////////////////////////

SurfDiffdef* Patchdef::surfdiffdef(surfdiff_local_id dlidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(dlidx < countSurfDiffs());
    return pStatedef->surfdiffdef(surfdiffL2G(dlidx));
}

////////////////////////////////////////////////////////////////////////////////

VDepSReacdef* Patchdef::vdepsreacdef(vdepsreac_local_id lidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < countVDepSReacs());
    return pStatedef->vdepsreacdef(vdepsreacL2G(lidx));
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurrdef* Patchdef::ohmiccurrdef(ohmiccurr_local_id lidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < countOhmicCurrs());
    return pStatedef->ohmiccurrdef(ohmiccurrL2G(lidx));
}

////////////////////////////////////////////////////////////////////////////////

GHKcurrdef* Patchdef::ghkcurrdef(ghkcurr_local_id lidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < countGHKcurrs());
    return pStatedef->ghkcurrdef(ghkcurrL2G(lidx));
}

////////////////////////////////////////////////////////////////////////////////

void Patchdef::setKcst(sreac_local_id srlidx, double kcst) {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    AssertLog(srlidx < countSReacs());
    AssertLog(kcst >= 0.0);
    pSReacKcst[srlidx] = kcst;
}

////////////////////////////////////////////////////////////////////////////////

void Patchdef::setActive(sreac_local_id srlidx, bool active) {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    AssertLog(srlidx < countSReacs());
    if (active == true) {
        pSReacFlags[srlidx] &= ~INACTIVATED;
    } else {
        pSReacFlags[srlidx] |= INACTIVATED;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Patchdef::setEndoKcst(endocytosis_local_id endolidx, double kcst) {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    AssertLog(endolidx < pEndocytosisN);
    AssertLog(kcst >= 0.0);
    pEndocytosisKcst[endolidx] = kcst;
}

////////////////////////////////////////////////////////////////////////////////

void Patchdef::setEndoActive(endocytosis_local_id endolidx, bool active) {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    AssertLog(endolidx < pEndocytosisN);
    if (active == true) {
        pEndocytosisFlags[endolidx] &= ~INACTIVATED;
    } else {
        pEndocytosisFlags[endolidx] |= INACTIVATED;
    }
}

////////////////////////////////////////////////////////////////////////////////

spec_local_id Patchdef::specG2L_I(spec_global_id gidx) const {
    if (pInner == nullptr) {
        return {};
    }
    return pInner->specG2L(gidx);
}

////////////////////////////////////////////////////////////////////////////////

spec_local_id Patchdef::specG2L_O(spec_global_id gidx) const {
    if (pOuter == nullptr) {
        return {};
    }
    return pOuter->specG2L(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int Patchdef::ohmiccurr_dep_S(ohmiccurr_local_id oclidx, spec_local_id splidx) const {
    return pOhmicCurr_DEP_Spec[splidx.get() + (oclidx.get() * countOhmicCurrs())];
}

////////////////////////////////////////////////////////////////////////////////

spec_local_id Patchdef::ohmiccurr_chanstate(ohmiccurr_local_id oclidx) const {
    return pOhmicCurr_CHANSTATE[oclidx];
}

////////////////////////////////////////////////////////////////////////////////

int Patchdef::ghkcurr_dep_S(ghkcurr_local_id ghklidx, spec_local_id splidx) const {
    return pGHKcurr_DEP_Spec[splidx.get() + (ghklidx.get() * countGHKcurrs())];
}

////////////////////////////////////////////////////////////////////////////////

spec_local_id Patchdef::ghkcurr_chanstate(ghkcurr_local_id ghklidx) const {
    return pGHKcurr_CHANSTATE[ghklidx];
}

////////////////////////////////////////////////////////////////////////////////

spec_local_id Patchdef::ghkcurr_ion(ghkcurr_local_id ghklidx) const {
    return pGHKcurr_ION[ghklidx];
}

////////////////////////////////////////////////////////////////////////////////

uint Patchdef::surfdiff_dep(surfdiff_local_id dlidx, spec_local_id slidx) const {
    return pSurfDiff_DEP_Spec[slidx.get() + (dlidx.get() * countSpecs())];
}

}  // namespace steps::solver
