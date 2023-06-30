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

#include "statedef.hpp"

// STEPS headers.
#include "chandef.hpp"
#include "compdef.hpp"
#include "diffboundarydef.hpp"
#include "diffdef.hpp"
#include "endocytosisdef.hpp"
#include "exocytosisdef.hpp"
#include "ghkcurrdef.hpp"
#include "linkspecdef.hpp"
#include "model/diff.hpp"
#include "model/reac.hpp"
#include "model/sreac.hpp"
#include "ohmiccurrdef.hpp"
#include "patchdef.hpp"
#include "raftdef.hpp"
#include "raftdisdef.hpp"
#include "raftendocytosisdef.hpp"
#include "raftgendef.hpp"
#include "raftsreacdef.hpp"
#include "reacdef.hpp"
#include "sdiffboundarydef.hpp"
#include "specdef.hpp"
#include "sreacdef.hpp"
#include "vdepsreacdef.hpp"
#include "vesbinddef.hpp"
#include "vesicledef.hpp"
#include "vessdiffdef.hpp"
#include "vessreacdef.hpp"
#include "vesunbinddef.hpp"
// util
#include "util/checkpointing.hpp"
#include "util/error.hpp"
// logging
#include <easylogging++.h>

namespace steps::solver {

Statedef::Statedef(model::Model* m, wm::Geom* g, const rng::RNGptr& r)
    : pModel(m)
    , pGeom(g)
    , pRNG(r)
    , pTime(0.0)
    , pNSteps(0) {
    AssertLog(pModel != nullptr);
    AssertLog(pGeom != nullptr);

    // Create the def objects.
    // NOTE: The order is very important. For example all objects after SpecDef
    // need to know the number of species in the state so SpecDef must be first;
    // CompDef must know what reacs and diffs are in the system, so Reacdef and
    // Diffdef must be created before Compdef. Compdef MUST COME BEFORE Patchef,
    // etc.
    //
    uint nspecs = pModel->_countSpecs();
    AssertLog(nspecs > 0);
    for (auto sidx: spec_global_id::range(nspecs)) {
        auto* specdef = new Specdef(this, sidx, pModel->_getSpec(sidx));
        AssertLog(specdef != nullptr);
        pSpecdefs.container().push_back(specdef);
    }

    uint nlspecs = pModel->_countLinkSpecs();
    for (auto lsidx: linkspec_global_id::range(nlspecs)) {
        auto* linkspecdef = new LinkSpecdef(this, lsidx, pModel->_getLinkSpec(lsidx));
        AssertLog(linkspecdef != nullptr);
        pLinkSpecdefs.container().push_back(linkspecdef);
    }

    uint nchans = pModel->_countChans();
    for (auto chidx: chan_global_id::range(nchans)) {
        auto* chandef = new Chandef(this, chidx, pModel->_getChan(chidx));
        AssertLog(chandef != nullptr);
        pChandefs.container().push_back(chandef);
    }

    uint nreacs = pModel->_countReacs();
    for (auto ridx: reac_global_id::range(nreacs)) {
        auto* reacdef = new Reacdef(this, ridx, pModel->_getReac(ridx));
        AssertLog(reacdef != nullptr);
        pReacdefs.container().push_back(reacdef);
    }

    uint nvdiffs = pModel->_countVDiffs();
    for (auto didx: diff_global_id::range(nvdiffs)) {
        auto* diffdef = new Diffdef(this, didx, pModel->_getVDiff(didx));
        AssertLog(diffdef != nullptr);
        pDiffdefs.container().push_back(diffdef);
    }

    uint nsdiffs = pModel->_countSDiffs();
    for (auto didx: surfdiff_global_id::range(nsdiffs)) {
        auto* surfdiffdef = new SurfDiffdef(this, didx, pModel->_getSDiff(didx));
        AssertLog(surfdiffdef != nullptr);
        pSurfDiffdefs.container().push_back(surfdiffdef);
    }

    uint nsreacs = pModel->_countSReacs();
    for (auto sridx: sreac_global_id::range(nsreacs)) {
        auto* sreacdef = new SReacdef(this, sridx, pModel->_getSReac(sridx));
        AssertLog(sreacdef != nullptr);
        pSReacdefs.container().push_back(sreacdef);
    }

    uint nendos = pModel->_countEndocytosis();
    for (auto endoidx: endocytosis_global_id::range(nendos)) {
        auto* endocytosisdef = new Endocytosisdef(this, endoidx, pModel->_getEndocytosis(endoidx));
        AssertLog(endocytosisdef != nullptr);
        pEndocytosisdefs.container().push_back(endocytosisdef);
    }

    uint nexos = pModel->_countExocytosis();
    for (auto exoidx: exocytosis_global_id::range(nexos)) {
        auto* exocytosisdef = new Exocytosisdef(this, exoidx, pModel->_getExocytosis(exoidx));
        AssertLog(exocytosisdef != nullptr);
        pExocytosisdefs.container().push_back(exocytosisdef);
    }

    // Does this need vesicles set up first? If so an impossible triangle because
    // comps have to come before vesicles, ves-ves interactions have to come
    // before comp
    uint nvesbinds = pModel->_countVesBinds();
    for (auto vesbidx: vesbind_global_id::range(nvesbinds)) {
        auto* vesbinddef = new VesBinddef(this, vesbidx, pModel->_getVesBind(vesbidx));
        AssertLog(vesbinddef != nullptr);
        pVesBinddefs.container().push_back(vesbinddef);
    }
    uint nvesunbinds = pModel->_countVesUnbinds();
    for (auto vesubidx: vesunbind_global_id::range(nvesunbinds)) {
        auto* vesunbinddef = new VesUnbinddef(this, vesubidx, pModel->_getVesUnbind(vesubidx));
        AssertLog(vesunbinddef != nullptr);
        pVesUnbinddefs.container().push_back(vesunbinddef);
    }

    uint nvessreacs = pModel->_countVesSReacs();
    for (auto vessridx: vessreac_global_id::range(nvessreacs)) {
        auto* vessreacdef = new VesSReacdef(this, vessridx, pModel->_getVesSReac(vessridx));
        AssertLog(vessreacdef != nullptr);
        pVesSReacdefs.container().push_back(vessreacdef);
    }

    uint nvsdiffs = pModel->_countVesSDiffs();
    for (auto didx: vessdiff_global_id::range(nvsdiffs)) {
        auto* vessdiffdef = new VesSDiffdef(this, didx, pModel->_getVesSDiff(didx));
        AssertLog(vessdiffdef != nullptr);
        pVesSDiffdefs.container().push_back(vessdiffdef);
    }

    uint nraftsreacs = pModel->_countRaftSReacs();
    for (auto raftsridx: raftsreac_global_id::range(nraftsreacs)) {
        auto* raftsreacdef = new RaftSReacdef(this, raftsridx, pModel->_getRaftSReac(raftsridx));
        AssertLog(raftsreacdef != nullptr);
        pRaftSReacdefs.container().push_back(raftsreacdef);
    }

    uint nraftendos = pModel->_countRaftEndocytosis();
    for (auto raftendoidx: raftendocytosis_global_id::range(nraftendos)) {
        auto* raftendocytosisdef =
            new RaftEndocytosisdef(this, raftendoidx, pModel->_getRaftEndocytosis(raftendoidx));
        AssertLog(raftendocytosisdef != nullptr);
        pRaftEndocytosisdefs.container().push_back(raftendocytosisdef);
    }

    uint nraftgens = pModel->_countRaftGeneses();
    for (auto raftgenidx: raftgen_global_id::range(nraftgens)) {
        auto* raftgendef = new RaftGendef(this, raftgenidx, pModel->_getRaftGen(raftgenidx));
        AssertLog(raftgendef != nullptr);
        pRaftGendefs.container().push_back(raftgendef);
    }

    uint nraftdiss = pModel->_countRaftDiss();
    for (auto raftdisidx: raftdis_global_id::range(nraftdiss)) {
        auto* raftdisdef = new RaftDisdef(this, raftdisidx, pModel->_getRaftDis(raftdisidx));
        AssertLog(raftdisdef != nullptr);
        pRaftDisdefs.container().push_back(raftdisdef);
    }

    uint nvdsreacs = pModel->_countVDepSReacs();
    for (auto vdsridx: vdepsreac_global_id::range(nvdsreacs)) {
        auto* vdsrdef = new VDepSReacdef(this, vdsridx, pModel->_getVDepSReac(vdsridx));
        AssertLog(vdsrdef != nullptr);
        pVDepSReacdefs.container().push_back(vdsrdef);
    }

    uint nohmiccurrs = pModel->_countOhmicCurrs();
    for (auto ocidx: ohmiccurr_global_id::range(nohmiccurrs)) {
        auto* ocdef = new OhmicCurrdef(this, ocidx, pModel->_getOhmicCurr(ocidx));
        AssertLog(ocdef != nullptr);
        pOhmicCurrdefs.container().push_back(ocdef);
    }

    uint nghkcurrs = pModel->_countGHKcurrs();
    for (auto ghkidx: ghkcurr_global_id::range(nghkcurrs)) {
        auto* ghkdef = new GHKcurrdef(this, ghkidx, pModel->_getGHKcurr(ghkidx));
        AssertLog(ghkdef != nullptr);
        pGHKcurrdefs.container().push_back(ghkdef);
    }

    uint ncomps = pGeom->_countComps();
    AssertLog(ncomps > 0);
    for (auto cidx: comp_global_id::range(ncomps)) {
        auto* compdef = new Compdef(this, cidx, pGeom->_getComp(cidx));
        AssertLog(compdef != nullptr);
        pCompdefs.container().push_back(compdef);
    }

    uint npatches = pGeom->_countPatches();
    for (auto pidx: patch_global_id::range(npatches)) {
        auto* patchdef = new Patchdef(this, pidx, pGeom->_getPatch(pidx));
        AssertLog(patchdef != nullptr);
        pPatchdefs.container().push_back(patchdef);
    }

    if (auto* tetmesh = dynamic_cast<tetmesh::Tetmesh*>(pGeom)) {
        uint ndiffbs = tetmesh->_countDiffBoundaries();
        for (auto dbidx: diffboundary_global_id::range(ndiffbs)) {
            auto* diffboundarydef =
                new DiffBoundarydef(this, dbidx, tetmesh->_getDiffBoundary(dbidx));
            AssertLog(diffboundarydef != nullptr);
            pDiffBoundarydefs.container().push_back(diffboundarydef);
        }

        uint nsdiffbs = tetmesh->_countSDiffBoundaries();
        for (auto sdbidx: sdiffboundary_global_id::range(nsdiffbs)) {
            auto* sdiffboundarydef =
                new SDiffBoundarydef(this, sdbidx, tetmesh->_getSDiffBoundary(sdbidx));
            AssertLog(sdiffboundarydef != nullptr);
            pSDiffBoundarydefs.container().push_back(sdiffboundarydef);
        }
    }

    // These HAVE TO come after setting up vessreacs etc
    uint nvesicles = pModel->_countVesicles();
    for (auto vidx: vesicle_global_id::range(nvesicles)) {
        auto* vesicledef = new Vesicledef(this, vidx, pModel->_getVesicle(vidx));
        AssertLog(vesicledef != nullptr);
        pVesicledefs.container().push_back(vesicledef);
    }

    uint nrafts = pModel->_countRafts();
    for (auto ridx: raft_global_id::range(nrafts)) {
        auto* raftdef = new Raftdef(this, ridx, pModel->_getRaft(ridx));
        AssertLog(raftdef != nullptr);
        pRaftdefs.container().push_back(raftdef);
    }

    // Now setup all the def objects. This can't be achieved purely with
    // the constructors, e.g.  a patch may need to add species from its
    // surface reactions to inner, outer comp
    // NOTE: Again, order is important.
    //
    for (auto const& specdef: pSpecdefs) {
        specdef->setup();
    }
    for (auto const& linkspecdef: pLinkSpecdefs) {
        linkspecdef->setup();
    }
    for (auto const& chandef: pChandefs) {
        chandef->setup();
    }
    for (auto const& reacdef: pReacdefs) {
        reacdef->setup();
    }
    for (auto const& diffdef: pDiffdefs) {
        diffdef->setup();
    }
    for (auto const& surfdiffdef: pSurfDiffdefs) {
        surfdiffdef->setup();
    }
    for (auto const& sreacdef: pSReacdefs) {
        sreacdef->setup();
    }

    for (auto const& endodef: pEndocytosisdefs) {
        endodef->setup();
    }
    for (auto const& exodef: pExocytosisdefs) {
        exodef->setup();
    }
    for (auto const& vsrdef: pVesSReacdefs) {
        vsrdef->setup();
    }
    for (auto const& vsddef: pVesSDiffdefs) {
        vsddef->setup();
    }
    for (auto const& rsrdef: pRaftSReacdefs) {
        rsrdef->setup();
    }
    for (auto const& rendodef: pRaftEndocytosisdefs) {
        rendodef->setup();
    }
    for (auto const& rgendef: pRaftGendefs) {
        rgendef->setup();
    }
    for (auto const& rdisdef: pRaftDisdefs) {
        rdisdef->setup();
    }

    for (auto& vdepeacdef: pVDepSReacdefs) {
        vdepeacdef->setup();
    }
    for (auto const& ohmcurrdef: pOhmicCurrdefs) {
        ohmcurrdef->setup();
    }
    for (auto const& ghkcurrdef: pGHKcurrdefs) {
        ghkcurrdef->setup();
    }

    for (auto const& compdef: pCompdefs) {
        compdef->setup_references();
    }
    for (auto const& patchdef: pPatchdefs) {
        patchdef->setup_references();
    }

    // Needs to come after the vesicle interactions (vessdiff, vessreac, ) have
    // been setup.
    for (auto const& vesdef: pVesicledefs) {
        vesdef->setup_references();
    }
    for (auto const& raftdef: pRaftdefs) {
        raftdef->setup_references();
    }

    // These vesicle binding events are a special case- they belong to the
    // compartment but need to know vesicle indices before setup, so they need to
    // be setup after vesicles setup their indices. This is only going to work if
    // compdef doesn't actually need to do any setup on these objects.
    for (auto const& vbdef: pVesBinddefs) {
        vbdef->setup();
    }
    for (auto const& vbdef: pVesUnbinddefs) {
        vbdef->setup();
    }

    // Make local indices for species, (surface) reactions, diffusion rules
    // in compartments then patches. Separate method necessary since e.g.
    // Patchdef::setup_references can add species to Compdef
    for (auto const& compdef: pCompdefs) {
        compdef->setup_indices();
    }
    for (auto const& patchdef: pPatchdefs) {
        patchdef->setup_indices();
    }

    for (auto const& vesdef: pVesicledefs) {
        vesdef->setup_indices();
    }

    for (auto const& raftdef: pRaftdefs) {
        raftdef->setup_indices();
    }

    for (auto& diffboundarydef: pDiffBoundarydefs) {
        diffboundarydef->setup();
    }
    for (auto& sdiffboundarydef: pSDiffBoundarydefs) {
        sdiffboundarydef->setup();
    }
}

////////////////////////////////////////////////////////////////////////////////

Statedef::~Statedef() {
    auto c_end = pCompdefs.end();
    for (auto c = pCompdefs.begin(); c != c_end; ++c) {
        delete *c;
    }

    auto p_end = pPatchdefs.end();
    for (auto p = pPatchdefs.begin(); p != p_end; ++p) {
        delete *p;
    }

    auto db_end = pDiffBoundarydefs.end();
    for (auto db = pDiffBoundarydefs.begin(); db != db_end; ++db) {
        delete *db;
    }

    auto sdb_end = pSDiffBoundarydefs.end();
    for (auto sdb = pSDiffBoundarydefs.begin(); sdb != sdb_end; ++sdb) {
        delete *sdb;
    }

    auto r_end = pReacdefs.end();
    for (auto r = pReacdefs.begin(); r != r_end; ++r) {
        delete *r;
    }

    auto sr_end = pSReacdefs.end();
    for (auto sr = pSReacdefs.begin(); sr != sr_end; ++sr) {
        delete *sr;
    }

    auto endo_end = pEndocytosisdefs.end();
    for (auto endo = pEndocytosisdefs.begin(); endo != endo_end; ++endo) {
        delete *endo;
    }

    auto exo_end = pExocytosisdefs.end();
    for (auto exo = pExocytosisdefs.begin(); exo != exo_end; ++exo) {
        delete *exo;
    }

    auto vesb_end = pVesBinddefs.end();
    for (auto vesb = pVesBinddefs.begin(); vesb != vesb_end; ++vesb) {
        delete *vesb;
    }

    auto vesub_end = pVesUnbinddefs.end();
    for (auto vesub = pVesUnbinddefs.begin(); vesub != vesub_end; ++vesub) {
        delete *vesub;
    }

    auto vessr_end = pVesSReacdefs.end();
    for (auto vessr = pVesSReacdefs.begin(); vessr != vessr_end; ++vessr) {
        delete *vessr;
    }

    auto vsd_end = pVesSDiffdefs.end();
    for (auto vsd = pVesSDiffdefs.begin(); vsd != vsd_end; ++vsd) {
        delete *vsd;
    }

    auto raftsr_end = pRaftSReacdefs.end();
    for (auto raftsr = pRaftSReacdefs.begin(); raftsr != raftsr_end; ++raftsr) {
        delete *raftsr;
    }

    auto rendo_end = pRaftEndocytosisdefs.end();
    for (auto rendo = pRaftEndocytosisdefs.begin(); rendo != rendo_end; ++rendo) {
        delete *rendo;
    }

    auto rgen_end = pRaftGendefs.end();
    for (auto rgen = pRaftGendefs.begin(); rgen != rgen_end; ++rgen) {
        delete *rgen;
    }

    auto rdis_end = pRaftDisdefs.end();
    for (auto rdis = pRaftDisdefs.begin(); rdis != rdis_end; ++rdis) {
        delete *rdis;
    }

    auto sd_end = pSurfDiffdefs.end();
    for (auto sd = pSurfDiffdefs.begin(); sd != sd_end; ++sd) {
        delete *sd;
    }

    auto d_end = pDiffdefs.end();
    for (auto d = pDiffdefs.begin(); d != d_end; ++d) {
        delete *d;
    }

    auto ch_end = pChandefs.end();
    for (auto ch = pChandefs.begin(); ch != ch_end; ++ch) {
        delete *ch;
    }

    auto vsr_end = pVDepSReacdefs.end();
    for (auto vsr = pVDepSReacdefs.begin(); vsr != vsr_end; ++vsr) {
        delete *vsr;
    }

    auto oc_end = pOhmicCurrdefs.end();
    for (auto oc = pOhmicCurrdefs.begin(); oc != oc_end; ++oc) {
        delete *oc;
    }

    auto ghk_end = pGHKcurrdefs.end();
    for (auto ghk = pGHKcurrdefs.begin(); ghk != ghk_end; ++ghk) {
        delete *ghk;
    }

    auto s_end = pSpecdefs.end();
    for (auto s = pSpecdefs.begin(); s != s_end; ++s) {
        delete *s;
    }

    auto ls_end = pLinkSpecdefs.end();
    for (auto ls = pLinkSpecdefs.begin(); ls != ls_end; ++ls) {
        delete *ls;
    }

    auto v_end = pVesicledefs.end();
    for (auto v = pVesicledefs.begin(); v != v_end; ++v) {
        delete *v;
    }

    auto raft_end = pRaftdefs.end();
    for (auto raft = pRaftdefs.begin(); raft != raft_end; ++raft) {
        delete *raft;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Statedef::checkpoint(std::fstream& cp_file) {
    for (auto const& s: pSpecdefs) {
        s->checkpoint(cp_file);
    }

    for (auto const& ls: pLinkSpecdefs) {
        ls->checkpoint(cp_file);
    }

    for (auto const& ch: pChandefs) {
        ch->checkpoint(cp_file);
    }

    for (auto const& c: pCompdefs) {
        c->checkpoint(cp_file);
    }

    for (auto const& p: pPatchdefs) {
        p->checkpoint(cp_file);
    }

    for (auto const& v: pVesicledefs) {
        v->checkpoint(cp_file);
    }

    for (auto const& raft: pRaftdefs) {
        raft->checkpoint(cp_file);
    }

    for (auto const& r: pReacdefs) {
        r->checkpoint(cp_file);
    }

    for (auto const& sr: pSReacdefs) {
        sr->checkpoint(cp_file);
    }

    for (auto const& d: pDiffdefs) {
        d->checkpoint(cp_file);
    }

    for (auto const& sd: pSurfDiffdefs) {
        sd->checkpoint(cp_file);
    }

    for (auto const& endo: pEndocytosisdefs) {
        endo->checkpoint(cp_file);
    }

    for (auto const& exo: pExocytosisdefs) {
        exo->checkpoint(cp_file);
    }

    for (auto const& vesb: pVesBinddefs) {
        vesb->checkpoint(cp_file);
    }

    for (auto const& vesub: pVesUnbinddefs) {
        vesub->checkpoint(cp_file);
    }

    for (auto const& vsr: pVesSReacdefs) {
        vsr->checkpoint(cp_file);
    }

    for (auto const& vsd: pVesSDiffdefs) {
        vsd->checkpoint(cp_file);
    }

    for (auto const& raftsr: pRaftSReacdefs) {
        raftsr->checkpoint(cp_file);
    }

    for (auto const& rendo: pRaftEndocytosisdefs) {
        rendo->checkpoint(cp_file);
    }

    for (auto const& rgen: pRaftGendefs) {
        rgen->checkpoint(cp_file);
    }

    for (auto const& rdis: pRaftDisdefs) {
        rdis->checkpoint(cp_file);
    }

    for (auto const& vdsr: pVDepSReacdefs) {
        vdsr->checkpoint(cp_file);
    }

    for (auto const& oc: pOhmicCurrdefs) {
        oc->checkpoint(cp_file);
    }

    for (auto const& ghkc: pGHKcurrdefs) {
        ghkc->checkpoint(cp_file);
    }

    for (auto const& db: pDiffBoundarydefs) {
        db->checkpoint(cp_file);
    }

    for (auto const& sdb: pSDiffBoundarydefs) {
        sdb->checkpoint(cp_file);
    }

    util::checkpoint(cp_file, pTime);
    util::checkpoint(cp_file, pNSteps);
}

////////////////////////////////////////////////////////////////////////////////

void Statedef::restore(std::fstream& cp_file) {
    for (auto& s: pSpecdefs) {
        s->restore(cp_file);
    }

    for (auto& ls: pLinkSpecdefs) {
        ls->restore(cp_file);
    }

    for (auto& ch: pChandefs) {
        ch->restore(cp_file);
    }

    for (auto& c: pCompdefs) {
        c->restore(cp_file);
    }

    for (auto& p: pPatchdefs) {
        p->restore(cp_file);
    }

    for (auto& v: pVesicledefs) {
        v->restore(cp_file);
    }

    for (auto& raft: pRaftdefs) {
        raft->restore(cp_file);
    }

    for (auto& r: pReacdefs) {
        r->restore(cp_file);
    }

    for (auto& sr: pSReacdefs) {
        sr->restore(cp_file);
    }

    for (auto& d: pDiffdefs) {
        d->restore(cp_file);
    }

    for (auto& sd: pSurfDiffdefs) {
        sd->restore(cp_file);
    }

    for (auto& endo: pEndocytosisdefs) {
        endo->restore(cp_file);
    }

    for (auto& exo: pExocytosisdefs) {
        exo->restore(cp_file);
    }

    for (auto& vesb: pVesBinddefs) {
        vesb->restore(cp_file);
    }

    for (auto& vesub: pVesUnbinddefs) {
        vesub->restore(cp_file);
    }

    for (auto& vsr: pVesSReacdefs) {
        vsr->restore(cp_file);
    }

    for (auto& vsd: pVesSDiffdefs) {
        vsd->restore(cp_file);
    }

    for (auto& raftsr: pRaftSReacdefs) {
        raftsr->restore(cp_file);
    }

    for (auto& rendo: pRaftEndocytosisdefs) {
        rendo->restore(cp_file);
    }

    for (auto& rgen: pRaftGendefs) {
        rgen->restore(cp_file);
    }

    for (auto& rdis: pRaftDisdefs) {
        rdis->restore(cp_file);
    }

    for (auto& vdsr: pVDepSReacdefs) {
        vdsr->restore(cp_file);
    }

    for (auto& oc: pOhmicCurrdefs) {
        oc->restore(cp_file);
    }

    for (auto& ghkc: pGHKcurrdefs) {
        ghkc->restore(cp_file);
    }

    for (auto& db: pDiffBoundarydefs) {
        db->restore(cp_file);
    }

    for (auto& sdb: pSDiffBoundarydefs) {
        sdb->restore(cp_file);
    }

    util::restore(cp_file, pTime);
    util::restore(cp_file, pNSteps);
}

////////////////////////////////////////////////////////////////////////////////

Compdef* Statedef::compdef(comp_global_id gidx) const {
    return pCompdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

comp_global_id Statedef::getCompIdx(std::string const& c) const {
    uint maxcidx = pCompdefs.size();
    AssertLog(maxcidx > 0);
    AssertLog(maxcidx == pGeom->_countComps());
    comp_global_id cidx(0u);
    while (cidx < maxcidx) {
        if (c == pGeom->_getComp(cidx)->getID()) {
            return cidx;
        }
        ++cidx;
    }
    std::ostringstream os;
    os << "Geometry does not contain comp with string identifier '" << c << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

comp_global_id Statedef::getCompIdx(wm::Comp* comp) const {
    uint maxcidx = pCompdefs.size();
    AssertLog(maxcidx > 0);
    AssertLog(maxcidx == pGeom->_countComps());
    comp_global_id cidx(0u);
    while (cidx < maxcidx) {
        if (comp == pGeom->_getComp(cidx)) {
            return cidx;
        }
        ++cidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

Patchdef* Statedef::patchdef(patch_global_id gidx) const {
    return pPatchdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

patch_global_id Statedef::getPatchIdx(std::string const& p) const {
    uint maxpidx = pPatchdefs.size();
    AssertLog(maxpidx == pGeom->_countPatches());
    patch_global_id pidx(0u);
    while (pidx < maxpidx) {
        if (p == pGeom->_getPatch(pidx)->getID()) {
            return pidx;
        }
        ++pidx;
    }
    std::ostringstream os;
    os << "Geometry does not contain patch with string identifier '" << p << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

patch_global_id Statedef::getPatchIdx(wm::Patch* patch) const {
    uint maxpidx = pPatchdefs.size();
    AssertLog(maxpidx == pGeom->_countPatches());
    patch_global_id pidx(0u);
    while (pidx < maxpidx) {
        if (patch == pGeom->_getPatch(pidx)) {
            return pidx;
        }
        ++pidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

Specdef* Statedef::specdef(spec_global_id gidx) const {
    return pSpecdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

spec_global_id Statedef::getSpecIdx(std::string const& s) const {
    uint maxsidx = pSpecdefs.size();
    AssertLog(maxsidx > 0);
    AssertLog(maxsidx == pModel->_countSpecs());
    spec_global_id sidx(0u);
    while (sidx < maxsidx) {
        if (s == pModel->_getSpec(sidx)->getID()) {
            return sidx;
        }
        ++sidx;
    }
    std::ostringstream os;
    os << "Model does not contain species with string identifier '" << s << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

spec_global_id Statedef::getSpecIdx(model::Spec* spec) const {
    uint maxsidx = pSpecdefs.size();
    AssertLog(maxsidx > 0);
    AssertLog(maxsidx == pModel->_countSpecs());
    spec_global_id sidx(0u);
    while (sidx < maxsidx) {
        if (spec == pModel->_getSpec(sidx)) {
            return sidx;
        }
        ++sidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

LinkSpecdef* Statedef::linkspecdef(linkspec_global_id gidx) const {
    return pLinkSpecdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

linkspec_global_id Statedef::getLinkSpecIdx(std::string const& s) const {
    uint maxsidx = pLinkSpecdefs.size();
    AssertLog(maxsidx > 0);
    AssertLog(maxsidx == pModel->_countLinkSpecs());
    linkspec_global_id sidx(0u);
    while (sidx < maxsidx) {
        if (s == pModel->_getLinkSpec(sidx)->getID()) {
            return sidx;
        }
        ++sidx;
    }
    std::ostringstream os;
    os << "Model does not contain link species with string identifier '" << s << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

linkspec_global_id Statedef::getLinkSpecIdx(model::LinkSpec* lspec) const {
    uint maxsidx = pLinkSpecdefs.size();
    AssertLog(maxsidx > 0);
    AssertLog(maxsidx == pModel->_countLinkSpecs());
    linkspec_global_id sidx(0u);
    while (sidx < maxsidx) {
        if (lspec == pModel->_getLinkSpec(sidx)) {
            return sidx;
        }
        ++sidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

Reacdef* Statedef::reacdef(reac_global_id gidx) const {
    return pReacdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

reac_global_id Statedef::getReacIdx(std::string const& r) const {
    uint maxridx = pReacdefs.size();
    AssertLog(maxridx == pModel->_countReacs());
    reac_global_id ridx(0u);
    while (ridx < maxridx) {
        if (r == pModel->_getReac(ridx)->getID()) {
            return ridx;
        }
        ++ridx;
    }
    std::ostringstream os;
    os << "Model does not contain reac with string identifier '" << r << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

reac_global_id Statedef::getReacIdx(model::Reac* reac) const {
    uint maxridx = pReacdefs.size();
    AssertLog(maxridx == pModel->_countReacs());
    reac_global_id ridx(0u);
    while (ridx < maxridx) {
        if (reac == pModel->_getReac(ridx)) {
            return ridx;
        }
        ++ridx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

SReacdef* Statedef::sreacdef(sreac_global_id gidx) const {
    return pSReacdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

sreac_global_id Statedef::getSReacIdx(std::string const& sr) const {
    uint maxsridx = pSReacdefs.size();
    AssertLog(maxsridx == pModel->_countSReacs());
    sreac_global_id sridx(0u);
    while (sridx < maxsridx) {
        if (sr == pModel->_getSReac(sridx)->getID()) {
            return sridx;
        }
        ++sridx;
    }
    std::ostringstream os;
    os << "Model does not contain sreac with string identifier '" << sr << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

sreac_global_id Statedef::getSReacIdx(model::SReac* sreac) const {
    uint maxsridx = pSReacdefs.size();
    AssertLog(maxsridx == pModel->_countSReacs());
    sreac_global_id sridx(0u);
    while (sridx < maxsridx) {
        if (sreac == pModel->_getSReac(sridx)) {
            return sridx;
        }
        ++sridx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

Diffdef* Statedef::diffdef(diff_global_id gidx) const {
    return pDiffdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

diff_global_id Statedef::getDiffIdx(std::string const& d) const {
    uint maxdidx = pDiffdefs.size();
    AssertLog(maxdidx == pModel->_countVDiffs());
    diff_global_id didx(0u);
    while (didx < maxdidx) {
        if (d == pModel->_getVDiff(didx)->getID()) {
            return didx;
        }
        ++didx;
    }
    std::ostringstream os;
    os << "Model does not contain diff with string identifier '" << d << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

diff_global_id Statedef::getDiffIdx(model::Diff* diff) const {
    uint maxdidx = pDiffdefs.size();
    AssertLog(maxdidx == pModel->_countVDiffs());
    diff_global_id didx(0u);
    while (didx < maxdidx) {
        if (diff == pModel->_getVDiff(didx)) {
            return didx;
        }
        ++didx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

SurfDiffdef* Statedef::surfdiffdef(surfdiff_global_id gidx) const {
    return pSurfDiffdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

surfdiff_global_id Statedef::getSurfDiffIdx(std::string const& d) const {
    uint maxdidx = pSurfDiffdefs.size();
    AssertLog(maxdidx == pModel->_countSDiffs());
    surfdiff_global_id didx(0u);
    while (didx < maxdidx) {
        if (d == pModel->_getSDiff(didx)->getID()) {
            return didx;
        }
        ++didx;
    }
    std::ostringstream os;
    os << "Model does not contain diff with string identifier '" << d << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

surfdiff_global_id Statedef::getSurfDiffIdx(model::Diff* diff) const {
    uint maxdidx = pSurfDiffdefs.size();
    AssertLog(maxdidx == pModel->_countSDiffs());
    surfdiff_global_id didx(0u);
    while (didx < maxdidx) {
        if (diff == pModel->_getSDiff(didx)) {
            return didx;
        }
        ++didx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

Vesicledef* Statedef::vesicledef(vesicle_global_id gidx) const {
    return pVesicledefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

vesicle_global_id Statedef::getVesicleIdx(std::string const& v) const {
    uint maxvidx = pVesicledefs.size();
    AssertLog(maxvidx > 0);
    AssertLog(maxvidx == pModel->_countVesicles());
    vesicle_global_id vidx(0u);
    while (vidx < maxvidx) {
        if (v == pModel->_getVesicle(vidx)->getID()) {
            return vidx;
        }
        ++vidx;
    }
    std::ostringstream os;
    os << "Model does not contain vesicle with string identifier '" << v << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

vesicle_global_id Statedef::getVesicleIdx(model::Vesicle* vesicle) const {
    uint maxvidx = pVesicledefs.size();
    AssertLog(maxvidx > 0);
    AssertLog(maxvidx == pModel->_countVesicles());
    vesicle_global_id vidx(0u);
    while (vidx < maxvidx) {
        if (vesicle == pModel->_getVesicle(vidx)) {
            return vidx;
        }
        ++vidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

Raftdef* Statedef::raftdef(raft_global_id gidx) const {
    return pRaftdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

raft_global_id Statedef::getRaftIdx(std::string const& r) const {
    uint maxridx = pRaftdefs.size();
    AssertLog(maxridx > 0);
    AssertLog(maxridx == pModel->_countRafts());
    raft_global_id ridx(0u);
    while (ridx < maxridx) {
        if (r == pModel->_getRaft(ridx)->getID()) {
            return ridx;
        }
        ++ridx;
    }
    std::ostringstream os;
    os << "Model does not contain raft with string identifier '" << r << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

raft_global_id Statedef::getRaftIdx(model::Raft* raft) const {
    uint maxridx = pRaftdefs.size();
    AssertLog(maxridx > 0);
    AssertLog(maxridx == pModel->_countRafts());
    raft_global_id ridx(0u);
    while (ridx < maxridx) {
        if (raft == pModel->_getRaft(ridx)) {
            return ridx;
        }
        ++ridx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

Endocytosisdef* Statedef::endocytosisdef(endocytosis_global_id gidx) const {
    return pEndocytosisdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

endocytosis_global_id Statedef::getEndocytosisIdx(std::string const& endo) const {
    uint maxendoidx = pEndocytosisdefs.size();
    AssertLog(maxendoidx == pModel->_countEndocytosis());
    endocytosis_global_id endoidx(0u);
    while (endoidx < maxendoidx) {
        if (endo == pModel->_getEndocytosis(endoidx)->getID()) {
            return endoidx;
        }
        ++endoidx;
    }
    std::ostringstream os;
    os << "Model does not contain endocytosis with string identifier '" << endo << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

endocytosis_global_id Statedef::getEndocytosisIdx(model::Endocytosis* endocyt) const {
    uint maxendoidx = pEndocytosisdefs.size();
    AssertLog(maxendoidx == pModel->_countEndocytosis());
    endocytosis_global_id endoidx(0u);
    while (endoidx < maxendoidx) {
        if (endocyt == pModel->_getEndocytosis(endoidx)) {
            return endoidx;
        }
        ++endoidx;
    }

    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

Exocytosisdef* Statedef::exocytosisdef(exocytosis_global_id gidx) const {
    return pExocytosisdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

exocytosis_global_id Statedef::getExocytosisIdx(std::string const& exo) const {
    uint maxexoidx = pExocytosisdefs.size();
    AssertLog(maxexoidx == pModel->_countExocytosis());
    exocytosis_global_id exoidx(0u);
    while (exoidx < maxexoidx) {
        if (exo == pModel->_getExocytosis(exoidx)->getID()) {
            return exoidx;
        }
        ++exoidx;
    }
    std::ostringstream os;
    os << "Model does not contain exocytosis with string identifier '" << exo << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

exocytosis_global_id Statedef::getExocytosisIdx(model::Exocytosis* exocyt) const {
    uint maxexoidx = pExocytosisdefs.size();
    AssertLog(maxexoidx == pModel->_countExocytosis());
    exocytosis_global_id exoidx(0u);
    while (exoidx < maxexoidx) {
        if (exocyt == pModel->_getExocytosis(exoidx)) {
            return exoidx;
        }
        ++exoidx;
    }

    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

VesBinddef* Statedef::vesbinddef(vesbind_global_id gidx) const {
    return pVesBinddefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

vesbind_global_id Statedef::getVesBindIdx(std::string const& vb) const {
    uint maxidx = pVesBinddefs.size();
    AssertLog(maxidx == pModel->_countVesBinds());
    vesbind_global_id vbidx(0u);
    while (vbidx < maxidx) {
        if (vb == pModel->_getVesBind(vbidx)->getID()) {
            return vbidx;
        }
        ++vbidx;
    }
    std::ostringstream os;
    os << "Model does not contain vesicle binding reaction with string "
          "identifier '"
       << vb << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

vesbind_global_id Statedef::getVesBindIdx(model::VesBind* vesbind) const {
    uint maxidx = pVesBinddefs.size();
    AssertLog(maxidx == pModel->_countVesBinds());
    vesbind_global_id vbidx(0u);
    while (vbidx < maxidx) {
        if (vesbind == pModel->_getVesBind(vbidx)) {
            return vbidx;
        }
        ++vbidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

VesUnbinddef* Statedef::vesunbinddef(vesunbind_global_id gidx) const {
    return pVesUnbinddefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

vesunbind_global_id Statedef::getVesUnbindIdx(std::string const& vb) const {
    uint maxidx = pVesUnbinddefs.size();
    AssertLog(maxidx == pModel->_countVesUnbinds());
    vesunbind_global_id vbidx(0u);
    while (vbidx < maxidx) {
        if (vb == pModel->_getVesUnbind(vbidx)->getID()) {
            return vbidx;
        }
        ++vbidx;
    }
    std::ostringstream os;
    os << "Model does not contain vesicle unbinding reaction with string "
          "identifier '"
       << vb << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

vesunbind_global_id Statedef::getVesUnbindIdx(model::VesUnbind* vesunbind) const {
    uint maxidx = pVesUnbinddefs.size();
    AssertLog(maxidx == pModel->_countVesUnbinds());
    vesunbind_global_id vbidx(0u);
    while (vbidx < maxidx) {
        if (vesunbind == pModel->_getVesUnbind(vbidx)) {
            return vbidx;
        }
        ++vbidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

VesSDiffdef* Statedef::vessdiffdef(vessdiff_global_id gidx) const {
    return pVesSDiffdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

vessdiff_global_id Statedef::getVesSDiffIdx(std::string const& vsd) const {
    uint maxdidx = pVesSDiffdefs.size();
    AssertLog(maxdidx == pModel->_countVesSDiffs());
    vessdiff_global_id didx(0u);
    while (didx < maxdidx) {
        if (vsd == pModel->_getVesSDiff(didx)->getID()) {
            return didx;
        }
        ++didx;
    }
    std::ostringstream os;
    os << "Model does not contain vesicle surface diffusion with string "
          "identifier '"
       << vsd << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

vessdiff_global_id Statedef::getVesSDiffIdx(model::VesSDiff* vsdiff) const {
    uint maxdidx = pVesSDiffdefs.size();
    AssertLog(maxdidx == pModel->_countVesSDiffs());
    vessdiff_global_id didx(0u);
    while (didx < maxdidx) {
        if (vsdiff == pModel->_getVesSDiff(didx)) {
            return didx;
        }
        ++didx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

VesSReacdef* Statedef::vessreacdef(vessreac_global_id gidx) const {
    return pVesSReacdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

vessreac_global_id Statedef::getVesSReacIdx(std::string const& vsr) const {
    uint maxdidx = pVesSReacdefs.size();
    AssertLog(maxdidx == pModel->_countVesSReacs());
    vessreac_global_id didx(0u);
    while (didx < maxdidx) {
        if (vsr == pModel->_getVesSReac(didx)->getID()) {
            return didx;
        }
        ++didx;
    }
    std::ostringstream os;
    os << "Model does not contain vesicle surface reaction with string "
          "identifier '"
       << vsr << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

vessreac_global_id Statedef::getVesSReacIdx(model::VesSReac* vsreac) const {
    uint maxdidx = pVesSReacdefs.size();
    AssertLog(maxdidx == pModel->_countVesSReacs());
    vessreac_global_id didx(0u);
    while (didx < maxdidx) {
        if (vsreac == pModel->_getVesSReac(didx)) {
            return didx;
        }
        ++didx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

RaftSReacdef* Statedef::raftsreacdef(raftsreac_global_id gidx) const {
    return pRaftSReacdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

raftsreac_global_id Statedef::getRaftSReacIdx(std::string const& rsr) const {
    uint maxrsridx = pRaftSReacdefs.size();
    AssertLog(maxrsridx == pModel->_countRaftSReacs());
    raftsreac_global_id rsridx(0u);
    while (rsridx < maxrsridx) {
        if (rsr == pModel->_getRaftSReac(rsridx)->getID()) {
            return rsridx;
        }
        ++rsridx;
    }
    std::ostringstream os;
    os << "Model does not contain raft surface reaction with string identifier '" << rsr << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

raftsreac_global_id Statedef::getRaftSReacIdx(model::RaftSReac* rsreac) const {
    uint maxrsridx = pRaftSReacdefs.size();
    AssertLog(maxrsridx == pModel->_countRaftSReacs());
    raftsreac_global_id rsridx(0u);
    while (rsridx < maxrsridx) {
        if (rsreac == pModel->_getRaftSReac(rsridx)) {
            return rsridx;
        }
        ++rsridx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

RaftGendef* Statedef::raftgendef(raftgen_global_id gidx) const {
    return pRaftGendefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

raftgen_global_id Statedef::getRaftGenIdx(std::string const& rgen) const {
    uint maxrgenidx = pRaftGendefs.size();
    AssertLog(maxrgenidx == pModel->_countRaftGeneses());
    raftgen_global_id rgenidx(0u);
    while (rgenidx < maxrgenidx) {
        if (rgen == pModel->_getRaftGen(rgenidx)->getID()) {
            return rgenidx;
        }
        ++rgenidx;
    }
    std::ostringstream os;
    os << "Model does not contain raft genesis with string identifier '" << rgen << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

raftgen_global_id Statedef::getRaftGenIdx(model::RaftGen* rgen) const {
    uint maxrgenidx = pRaftGendefs.size();
    AssertLog(maxrgenidx == pModel->_countRaftGeneses());
    raftgen_global_id rgenidx(0u);
    while (rgenidx < maxrgenidx) {
        if (rgen == pModel->_getRaftGen(rgenidx)) {
            return rgenidx;
        }
        ++rgenidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

RaftDisdef* Statedef::raftdisdef(raftdis_global_id gidx) const {
    return pRaftDisdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

raftdis_global_id Statedef::getRaftDisIdx(std::string const& rdis) const {
    uint maxrdisidx = pRaftDisdefs.size();
    AssertLog(maxrdisidx == pModel->_countRaftDiss());
    raftdis_global_id rdisidx(0u);
    while (rdisidx < maxrdisidx) {
        if (rdis == pModel->_getRaftDis(rdisidx)->getID()) {
            return rdisidx;
        }
        ++rdisidx;
    }
    std::ostringstream os;
    os << "Model does not contain raft dissolution with string identifier '" << rdis << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

raftdis_global_id Statedef::getRaftDisIdx(model::RaftDis* rdis) const {
    uint maxrdisidx = pRaftDisdefs.size();
    AssertLog(maxrdisidx == pModel->_countRaftDiss());
    raftdis_global_id rdisidx(0u);
    while (rdisidx < maxrdisidx) {
        if (rdis == pModel->_getRaftDis(rdisidx)) {
            return rdisidx;
        }
        ++rdisidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}
////////////////////////////////////////////////////////////////////////////////

RaftEndocytosisdef* Statedef::raftendocytosisdef(raftendocytosis_global_id gidx) const {
    return pRaftEndocytosisdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

raftendocytosis_global_id Statedef::getRaftEndocytosisIdx(std::string const& rendo) const {
    uint maxendoidx = pRaftEndocytosisdefs.size();
    AssertLog(maxendoidx == pModel->_countRaftEndocytosis());
    raftendocytosis_global_id endoidx(0u);
    while (endoidx < maxendoidx) {
        if (rendo == pModel->_getRaftEndocytosis(endoidx)->getID()) {
            return endoidx;
        }
        ++endoidx;
    }
    std::ostringstream os;
    os << "Model does not contain raft endocytosis with string identifier '" << rendo << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

raftendocytosis_global_id Statedef::getRaftEndocytosisIdx(model::RaftEndocytosis* rendocyt) const {
    uint maxendoidx = pRaftEndocytosisdefs.size();
    AssertLog(maxendoidx == pModel->_countRaftEndocytosis());
    raftendocytosis_global_id endoidx(0u);
    while (endoidx < maxendoidx) {
        if (rendocyt == pModel->_getRaftEndocytosis(endoidx)) {
            return endoidx;
        }
        ++endoidx;
    }

    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurrdef* Statedef::ohmiccurrdef(ohmiccurr_global_id gidx) const {
    return pOhmicCurrdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

ohmiccurr_global_id Statedef::getOhmicCurrIdx(std::string const& oc) const {
    uint maxocidx = pOhmicCurrdefs.size();
    AssertLog(maxocidx == pModel->_countOhmicCurrs());
    ohmiccurr_global_id ocidx(0u);
    while (ocidx < maxocidx) {
        if (oc == pModel->_getOhmicCurr(ocidx)->getID()) {
            return ocidx;
        }
        ++ocidx;
    }
    std::ostringstream os;
    os << "Model does not contain ohmic current with string identifier '" << oc << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

ohmiccurr_global_id Statedef::getOhmicCurrIdx(model::OhmicCurr* ohmiccurr) const {
    uint maxocidx = pOhmicCurrdefs.size();
    AssertLog(maxocidx == pModel->_countOhmicCurrs());
    ohmiccurr_global_id ocidx(0u);
    while (ocidx < maxocidx) {
        if (ohmiccurr == pModel->_getOhmicCurr(ocidx)) {
            return ocidx;
        }
        ++ocidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

VDepSReacdef* Statedef::vdepsreacdef(vdepsreac_global_id gidx) const {
    return pVDepSReacdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

vdepsreac_global_id Statedef::getVDepSReacIdx(std::string const& vdsr) const {
    uint maxvdsridx = pVDepSReacdefs.size();
    AssertLog(maxvdsridx == pModel->_countVDepSReacs());
    vdepsreac_global_id vdsridx(0u);
    while (vdsridx < maxvdsridx) {
        if (vdsr == pModel->_getVDepSReac(vdsridx)->getID()) {
            return vdsridx;
        }
        ++vdsridx;
    }
    std::ostringstream os;
    os << "Model does not contain voltage-dependent reaction with string "
          "identifier '"
       << vdsr << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

vdepsreac_global_id Statedef::getVDepSReacIdx(model::VDepSReac* vdepsreac) const {
    uint maxvdsridx = pVDepSReacdefs.size();
    AssertLog(maxvdsridx == pModel->_countVDepSReacs());
    vdepsreac_global_id vdsridx(0u);
    while (vdsridx < maxvdsridx) {
        if (vdepsreac == pModel->_getVDepSReac(vdsridx)) {
            return vdsridx;
        }
        ++vdsridx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

GHKcurrdef* Statedef::ghkcurrdef(ghkcurr_global_id gidx) const {
    return pGHKcurrdefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

ghkcurr_global_id Statedef::getGHKcurrIdx(std::string const& ghk) const {
    uint maxghkidx = pGHKcurrdefs.size();
    AssertLog(maxghkidx == pModel->_countGHKcurrs());
    ghkcurr_global_id ghkidx(0u);
    while (ghkidx < maxghkidx) {
        if (ghk == pModel->_getGHKcurr(ghkidx)->getID()) {
            return ghkidx;
        }
        ++ghkidx;
    }
    std::ostringstream os;
    os << "Model does not contain ghk current with string identifier '" << ghk << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

ghkcurr_global_id Statedef::getGHKcurrIdx(model::GHKcurr* ghkcurr) const {
    uint maxghkidx = pGHKcurrdefs.size();
    AssertLog(maxghkidx == pModel->_countGHKcurrs());
    ghkcurr_global_id ghkidx(0u);
    while (ghkidx < maxghkidx) {
        if (ghkcurr == pModel->_getGHKcurr(ghkidx)) {
            return ghkidx;
        }
        ++ghkidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

membrane_global_id Statedef::getMembIdx(std::string const& m) const {
    if (auto* tetmesh = dynamic_cast<tetmesh::Tetmesh*>(pGeom)) {
        /* This isn't the right place for this check, rather the setup should check
        for
         * multiple membranes
        uint nmembs = tetmesh->_countMembs();
        if (nmembs != 1)
        {
            std::ostringstream os;
            os << "Only one Membrane may exist in simulation";
            ArgErrLog(os.str());
        }
        */
        tetmesh::Memb* membrane = tetmesh->_getMemb(0);
        if (m == membrane->getID()) {
            return membrane_global_id(0u);
        } else {
            std::ostringstream os;
            os << "Geometry does not contain membrane with string identifier '" << m << "'.";
            ArgErrLog(os.str());
        }

    } else {
        std::ostringstream os;
        os << "Membrane methods not available with well-mixed geometry";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

DiffBoundarydef* Statedef::diffboundarydef(diffboundary_global_id gidx) const {
    return pDiffBoundarydefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

diffboundary_global_id Statedef::getDiffBoundaryIdx(std::string const& d) const {
    uint maxdidx = pDiffBoundarydefs.size();
    if (auto* tetmesh = dynamic_cast<tetmesh::Tetmesh*>(pGeom)) {
        AssertLog(maxdidx == tetmesh->_countDiffBoundaries());
        diffboundary_global_id didx(0u);
        while (didx < maxdidx) {
            if (d == tetmesh->_getDiffBoundary(didx)->getID()) {
                return didx;
            }
            ++didx;
        }
        std::ostringstream os;
        os << "Geometry does not contain diff boundary with string identifier '" << d << "'.";
        ArgErrLog(os.str());
    } else {
        std::ostringstream os;
        os << "Diffusion boundary methods not available with well-mixed geometry";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

diffboundary_global_id Statedef::getDiffBoundaryIdx(tetmesh::DiffBoundary* diffb) const {
    uint maxdidx = pDiffBoundarydefs.size();
    if (auto* tetmesh = dynamic_cast<tetmesh::Tetmesh*>(pGeom)) {
        AssertLog(maxdidx == tetmesh->_countDiffBoundaries());
        diffboundary_global_id didx(0u);
        while (didx < maxdidx) {
            if (diffb == tetmesh->_getDiffBoundary(didx)) {
                return didx;
            }
            ++didx;
        }
        // Argument should be valid so we should not get here
        AssertLog(false);
    } else {
        std::ostringstream os;
        os << "Diffusion boundary methods not available with well-mixed geometry";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

SDiffBoundarydef* Statedef::sdiffboundarydef(sdiffboundary_global_id gidx) const {
    return pSDiffBoundarydefs.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

sdiffboundary_global_id Statedef::getSDiffBoundaryIdx(std::string const& sd) const {
    uint maxsdidx = pSDiffBoundarydefs.size();
    if (auto* tetmesh = dynamic_cast<tetmesh::Tetmesh*>(pGeom)) {
        AssertLog(maxsdidx == tetmesh->_countSDiffBoundaries());
        sdiffboundary_global_id sdidx(0u);
        while (sdidx < maxsdidx) {
            if (sd == tetmesh->_getSDiffBoundary(sdidx)->getID()) {
                return sdidx;
            }
            ++sdidx;
        }
        std::ostringstream os;
        os << "Geometry does not contain surface diffusion boundary with string "
              "identifier '"
           << sd << "'.";
        ArgErrLog(os.str());
    } else {
        std::ostringstream os;
        os << "Surface Diffusion Boundary methods not available with well-mixed "
              "geometry";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

sdiffboundary_global_id Statedef::getSDiffBoundaryIdx(tetmesh::SDiffBoundary* sdiffb) const {
    uint maxsdidx = pSDiffBoundarydefs.size();
    if (auto* tetmesh = dynamic_cast<tetmesh::Tetmesh*>(pGeom)) {
        AssertLog(maxsdidx == tetmesh->_countSDiffBoundaries());
        sdiffboundary_global_id sdidx(0u);
        while (sdidx < maxsdidx) {
            if (sdiffb == tetmesh->_getSDiffBoundary(sdidx)) {
                return sdidx;
            }
            ++sdidx;
        }
        // Argument should be valid so we should not get here
        AssertLog(false);
    } else {
        std::ostringstream os;
        os << "Surface Diffusion Boundary methods not available with well-mixed "
              "geometry";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Statedef::setTime(double t) {
    AssertLog(t >= 0.0);
    pTime = t;
}

////////////////////////////////////////////////////////////////////////////////

void Statedef::incTime(double dt) {
    AssertLog(dt >= 0.0);
    pTime += dt;
}

////////////////////////////////////////////////////////////////////////////////

void Statedef::incNSteps(uint i) {
    AssertLog(i != 0);
    pNSteps += i;
}

}  // namespace steps::solver
