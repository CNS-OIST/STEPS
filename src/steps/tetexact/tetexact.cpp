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

#include "tetexact.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>

#include "diff.hpp"
#include "ghkcurr.hpp"
#include "reac.hpp"
#include "sdiff.hpp"
#include "sreac.hpp"
#include "vdepsreac.hpp"
#include "wmvol.hpp"

#include "geom/comp.hpp"
#include "geom/memb.hpp"
#include "geom/tetmesh.hpp"
#include "geom/tmcomp.hpp"
#include "math/constants.hpp"
#include "math/point.hpp"
#include "solver/chandef.hpp"
#include "solver/compdef.hpp"
#include "solver/diffboundarydef.hpp"
#include "solver/efield/dVsolver.hpp"
#include "solver/ghkcurrdef.hpp"
#include "solver/ohmiccurrdef.hpp"
#include "solver/patchdef.hpp"
#include "solver/sdiffboundarydef.hpp"
#include "util/collections.hpp"
#include "util/distribute.hpp"
#include "util/error.hpp"

#include "util/profile/profiler_interface.hpp"

#include "util/checkpointing.hpp"

namespace steps::tetexact {

void schedIDXSet_To_Vec(SchedIDXSet const& s, SchedIDXVec& v) {
    v.resize(s.size());
    std::copy(s.begin(), s.end(), v.begin());
}

////////////////////////////////////////////////////////////////////////////////

Tetexact::Tetexact(model::Model* m, wm::Geom* g, const rng::RNGptr& r, int calcMembPot)
    : API(*m, *g, r)
    , pEFoption(static_cast<EF_solver>(calcMembPot)) {
    if (rng() == nullptr) {
        std::ostringstream os;
        os << "No RNG provided to solver initializer function";
        ArgErrLog(os.str());
    }

    // All initialization code now in _setup() to allow EField solver to be
    // derived and create EField local objects within the constructor
    _setup();
}

////////////////////////////////////////////////////////////////////////////////

Tetexact::~Tetexact() {
    for (auto const& c: pComps) {
        delete c;
    }
    for (auto const& p: pPatches) {
        delete p;
    }
    for (auto const& db: pDiffBoundaries) {
        delete db;
    }
    for (auto const& sdb: pSDiffBoundaries) {
        delete sdb;
    }
    for (auto const& wvol: pWmVols) {
        delete wvol;
    }
    for (auto const& t: pTets) {
        delete t;
    }
    for (auto const& t: pTris) {
        delete t;
    }
    for (auto const& g: nGroups) {
        g->free_indices();
        delete g;
    }
    for (auto const& g: pGroups) {
        g->free_indices();
        delete g;
    }
}

///////////////////////////////////////////////////////////////////////////////

void Tetexact::checkpoint(std::string const& file_name) {
    CLOG(INFO, "general_log") << "Checkpoint to " << file_name << "...";
    std::fstream cp_file;

    cp_file.open(file_name.c_str(), std::fstream::out | std::fstream::binary | std::fstream::trunc);

    API::checkpoint(cp_file);

    statedef().checkpoint(cp_file);

    for (auto c: pComps) {
        c->checkpoint(cp_file);
    }

    for (auto p: pPatches) {
        p->checkpoint(cp_file);
    }

    for (auto db: pDiffBoundaries) {
        db->checkpoint(cp_file);
    }

    for (auto sdb: pSDiffBoundaries) {
        sdb->checkpoint(cp_file);
    }

    for (auto wmv: pWmVols) {
        if (wmv != nullptr) {
            wmv->checkpoint(cp_file);
        }
    }

    for (auto t: pTets) {
        if (t != nullptr) {
            t->checkpoint(cp_file);
        }
    }

    for (auto t: pTris) {
        if (t != nullptr) {
            t->checkpoint(cp_file);
        }
    }

    for (auto k: pKProcs) {
        k->checkpoint(cp_file);
    }

    util::checkpoint(cp_file, nEntries);

    // checkpoint CR SSA
    util::checkpoint(cp_file, pSum);
    util::checkpoint(cp_file, nSum);
    util::checkpoint(cp_file, pA0);

    uint n_ngroups = nGroups.size();
    uint n_pgroups = pGroups.size();
    util::checkpoint(cp_file, n_ngroups);
    util::checkpoint(cp_file, n_pgroups);

    for (uint i = 0; i < n_ngroups; i++) {
        CRGroup* group = nGroups[i];
        util::checkpoint(cp_file, group->capacity);
        util::checkpoint(cp_file, group->size);
        util::checkpoint(cp_file, group->max);
        util::checkpoint(cp_file, group->sum);

        for (uint j = 0; j < group->size; j++) {
            uint idx = group->indices[j]->schedIDX().get();
            util::checkpoint(cp_file, idx);
        }
    }

    for (uint i = 0; i < n_pgroups; i++) {
        CRGroup* group = pGroups[i];
        util::checkpoint(cp_file, group->capacity);
        util::checkpoint(cp_file, group->size);
        util::checkpoint(cp_file, group->max);
        util::checkpoint(cp_file, group->sum);

        for (uint j = 0; j < group->size; j++) {
            uint idx = group->indices[j]->schedIDX().get();
            util::checkpoint(cp_file, idx);
        }
    }

    if (efflag()) {
        util::checkpoint(cp_file, pTemp);
        util::checkpoint(cp_file, pEFDT);
        pEField->checkpoint(cp_file);
    }

    cp_file.close();
    CLOG(INFO, "general_log") << "complete.\n";
}

///////////////////////////////////////////////////////////////////////////////

void Tetexact::restore(std::string const& file_name) {
    std::fstream cp_file;

    cp_file.open(file_name.c_str(), std::fstream::in | std::fstream::binary);

    cp_file.seekg(0);

    API::restore(cp_file);

    statedef().restore(cp_file);

    for (auto c: pComps) {
        c->restore(cp_file);
    }

    for (auto p: pPatches) {
        p->restore(cp_file);
    }

    for (auto db: pDiffBoundaries) {
        db->restore(cp_file);
    }

    for (auto sdb: pSDiffBoundaries) {
        sdb->restore(cp_file);
    }

    for (auto wmv: pWmVols) {
        if (wmv != nullptr) {
            wmv->restore(cp_file);
        }
    }

    for (auto t: pTets) {
        if (t != nullptr) {
            t->restore(cp_file);
        }
    }
    for (auto t: pTris) {
        if (t != nullptr) {
            t->restore(cp_file);
        }
    }

    for (auto k: pKProcs) {
        k->restore(cp_file);
    }

    util::compare(cp_file, nEntries, "Mismatched KProc count from restore.");

    // restore CR SSA
    util::restore(cp_file, pSum);
    util::restore(cp_file, nSum);
    util::restore(cp_file, pA0);

    uint n_ngroups;
    uint n_pgroups;

    util::restore(cp_file, n_ngroups);
    util::restore(cp_file, n_pgroups);
    nGroups.resize(n_ngroups);
    pGroups.resize(n_pgroups);

    for (uint i = 0; i < n_ngroups; i++) {
        unsigned capacity;
        unsigned size;
        double max;
        double sum;

        util::restore(cp_file, capacity);
        util::restore(cp_file, size);
        util::restore(cp_file, max);
        util::restore(cp_file, sum);

        nGroups[i] = new CRGroup(0, capacity);
        nGroups[i]->size = size;
        nGroups[i]->max = max;
        nGroups[i]->sum = sum;

        for (uint j = 0; j < size; j++) {
            uint idx;
            util::restore(cp_file, idx);
            nGroups[i]->indices[j] = pKProcs[idx];
        }
    }

    for (uint i = 0; i < n_pgroups; i++) {
        unsigned capacity;
        unsigned size;
        double max;
        double sum;

        util::restore(cp_file, capacity);
        util::restore(cp_file, size);
        util::restore(cp_file, max);
        util::restore(cp_file, sum);

        pGroups[i] = new CRGroup(0, capacity);
        pGroups[i]->size = size;
        pGroups[i]->max = max;
        pGroups[i]->sum = sum;

        for (uint j = 0; j < size; j++) {
            uint idx;
            util::restore(cp_file, idx);
            pGroups[i]->indices[j] = pKProcs[idx];
        }
    }

    if (efflag()) {
        util::restore(cp_file, pTemp);
        util::restore(cp_file, pEFDT);
        pEField->restore(cp_file);
    }

    cp_file.close();
}

////////////////////////////////////////////////////////////////////////////////


std::string Tetexact::getSolverName() const {
    return "tetexact";
}

////////////////////////////////////////////////////////////////////////////////

std::string Tetexact::getSolverDesc() const {
    return "SSA Composition and Rejection Exact Method in tetrahedral mesh";
}

////////////////////////////////////////////////////////////////////////////////

std::string Tetexact::getSolverAuthors() const {
    return "Stefan Wils, Iain Hepburn, Weiliang Chen";
}

////////////////////////////////////////////////////////////////////////////////

std::string Tetexact::getSolverEmail() const {
    return "steps.dev@gmail.com";
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setup() {
    // Perform upcast.
    pMesh = dynamic_cast<tetmesh::Tetmesh*>(&geom());
    if (pMesh == nullptr) {
        ArgErrLog(
            "Geometry description to solver::Tetexact solver "
            "constructor is not a valid tetmesh::Tetmesh object.");
    }

    // First initialise the pTets, pTris vector, because
    // want tets and tris to maintain indexing from Geometry
    auto ntets = pMesh->countTets();
    auto ntris = pMesh->countTris();
    auto ncomps = pMesh->_countComps();

    pTets.container().assign(ntets, nullptr);
    pTris.container().assign(ntris, nullptr);
    pWmVols.container().assign(ncomps, nullptr);

    // Now create the actual compartments.
    for (auto const& c: statedef().comps()) {
        const auto compdef_gidx = c->gidx();
        const auto comp_idx = _addComp(c.get());
        AssertLog(compdef_gidx == comp_idx);
    }
    // Create the actual patches.
    for (auto const& p: statedef().patches()) {
        const auto patchdef_gidx = p->gidx();
        const auto patch_idx = _addPatch(p.get());
        AssertLog(patchdef_gidx == patch_idx);
    }

    // Create the diffusion boundaries
    for (auto const& db: statedef().diffBoundaries()) {
        const auto diffboundary_gidx = db->gidx();
        const auto diffb_idx = _addDiffBoundary(db.get());
        AssertLog(diffboundary_gidx == diffb_idx);
    }

    // Create the surface diffusion boundaries
    for (auto const& sdb: statedef().sdiffBoundaries()) {
        const auto sdiffboundary_gidx = sdb->gidx();
        const auto sdiffb_idx = _addSDiffBoundary(sdb.get());
        AssertLog(sdiffboundary_gidx == sdiffb_idx);
    }

    const auto npatches = pPatches.size();
    AssertLog(pMesh->_countPatches() == npatches);
    for (auto p: solver::patch_global_id::range(npatches)) {
        // Add the tris for this patch
        // We have checked the indexing - p is the global index
        auto& wmpatch = pMesh->_getPatch(p);

        // Perform upcast
        auto tmpatch = dynamic_cast<tetmesh::TmPatch*>(&wmpatch);
        if (tmpatch == nullptr) {
            ArgErrLog("Well-mixed patches not supported in solver::Tetexact solver.");
        }
        auto* localpatch = pPatches[p];
        std::map<bar_id_t, std::vector<triangle_global_id>> bar2tri;

        // We need to go through all patches to record bar2tri mapping
        // for all connected triangle neighbors even they are in different
        // patches, because their information is needed for surface diffusion boundary

        for (auto bar_p: solver::patch_global_id::range(npatches)) {
            auto& patch = pMesh->_getPatch(bar_p);
            auto& bar_patch = dynamic_cast<tetmesh::TmPatch&>(patch);

            for (auto tri: bar_patch._getAllTriIndices()) {
                const auto bars = pMesh->_getTriBars(tri);
                for (uint i = 0; i < bars.size(); ++i) {
                    bar2tri[bars[i]].push_back(tri);
                }
            }
        }

        for (auto tri: tmpatch->_getAllTriIndices()) {
            AssertLog(pMesh->getTriPatch(tri) == tmpatch);

            auto area = pMesh->getTriArea(tri);

            // NB: Tri vertices may not be in consistent order, so use bar interface.
            const auto tri_bars = pMesh->_getTriBars(tri);
            double l[tri_bars.size()] = {0, 0, 0};

            for (uint j = 0; j < tri_bars.size(); ++j) {
                const auto v = pMesh->_getBar(tri_bars[j]);
                l[j] = distance(pMesh->_getVertex(v[0]), pMesh->_getVertex(v[1]));
            }


            // We need to get triangle neighbors in other patches as well
            // for surface diffution boundary

            // slow version as it loops over all triangles evern if they are not in a patch
            // std::vector<int> tris = pMesh->getTriTriNeighb(tri);


            std::array<triangle_global_id, 3> tris{{std::nullopt, std::nullopt, std::nullopt}};
            for (uint j = 0; j < tris.size(); ++j) {
                const auto& neighb_tris = bar2tri[tri_bars[j]];
                for (auto const& neighb_tri: neighb_tris) {
                    if (neighb_tri == tri || pMesh->getTriPatch(neighb_tri) == nullptr) {
                        continue;
                    }
                    tris[j] = neighb_tri;
                    break;
                }
            }

            const math::point3d baryc = pMesh->_getTriBarycenter(tri);

            double d[3] = {0, 0, 0};
            for (uint j = 0; j < 3; ++j) {
                if (tris[j].unknown()) {
                    continue;
                }
                d[j] = distance(baryc, pMesh->_getTriBarycenter(tris[j]));
            }

            const auto tri_tets = pMesh->_getTriTetNeighb(tri);
            _addTri(tri,
                    localpatch,
                    area,
                    l[0],
                    l[1],
                    l[2],
                    d[0],
                    d[1],
                    d[2],
                    tri_tets[0],
                    tri_tets[1],
                    tris[0],
                    tris[1],
                    tris[2]);
        }
    }

    ncomps = pComps.size();
    AssertLog(pMesh->_countComps() == ncomps);

    for (auto c: solver::comp_global_id::range(ncomps)) {
        // Now add the tets for this comp
        // We have checked the indexing- c is the global index
        auto& wmcomp = pMesh->_getComp(c);

        // Perform upcast
        auto tmcomp = dynamic_cast<tetmesh::TmComp*>(&wmcomp);
        if (tmcomp != nullptr) {
            auto* localcomp = pComps[c];

            for (auto tet: tmcomp->_getAllTetIndices()) {
                AssertLog(pMesh->getTetComp(tet) == tmcomp);

                double vol = pMesh->getTetVol(tet);

                const auto& tris = pMesh->_getTetTriNeighb(tet);

                const std::array<double, 4> a{{
                    pMesh->getTriArea(tris[0]),
                    pMesh->getTriArea(tris[1]),
                    pMesh->getTriArea(tris[2]),
                    pMesh->getTriArea(tris[3]),
                }};

                const auto tets = pMesh->_getTetTetNeighb(tet);
                math::point3d baryc = pMesh->_getTetBarycenter(tet);

                double d[4] = {0, 0, 0, 0};
                for (uint j = 0; j < 4; ++j) {
                    if (tets[j].unknown()) {
                        continue;
                    }
                    d[j] = distance(baryc, pMesh->_getTetBarycenter(tets[j]));
                }

                _addTet(tet,
                        localcomp,
                        vol,
                        a[0],
                        a[1],
                        a[2],
                        a[3],
                        d[0],
                        d[1],
                        d[2],
                        d[3],
                        tets[0],
                        tets[1],
                        tets[2],
                        tets[3]);
            }
        } else {
            // This means that this compartment is a well-mixed compartment
            // It will behave like a tetrahedral-based compartment, but
            // contain only one 'voxel' that is connected to all surface
            // triangles and has the same volume as the whole compartment
            steps::tetexact::Comp* localcomp = pComps[c];
            _addWmVol(c, localcomp, localcomp->def()->vol());
            AssertLog(pWmVols[c] != nullptr);

            // Now find all the triangles that reference this well-mixed volume
            // and set the inner or outer tetrahedron index accordingly.

            uint nopatches = wmcomp._countOPatches();
            for (uint i = 0; i < nopatches; ++i) {
                wm::Patch* op = wmcomp._getOPatch(i);
                //     Comp may have no outer patch
                if (op == nullptr) {
                    continue;
                }

                auto comp_opatch = dynamic_cast<tetmesh::TmPatch*>(op);
                if (comp_opatch == nullptr) {
                    ProgErrLog("Compartment outer patch is not a TmPatch.");
                }

                for (auto tri: comp_opatch->_getAllTriIndices()) {
                    pTris[tri]->setInnerTet(pWmVols[c]);
                    // Add triangle to WmVols' table of neighbouring triangles.
                    pWmVols[c]->setNextTri(pTris[tri]);
                }
            }

            auto nipatches = wmcomp._countIPatches();
            for (decltype(nipatches) i = 0; i < nipatches; ++i) {
                auto* ip = wmcomp._getIPatch(i);
                // Comp may not have an inner patch
                if (ip == nullptr) {
                    continue;
                }

                auto comp_ipatch = dynamic_cast<tetmesh::TmPatch*>(ip);
                if (comp_ipatch == nullptr) {
                    ProgErrLog("Compartment inner patch is not a TmPatch.");
                }

                for (auto tri: comp_ipatch->_getAllTriIndices()) {
                    pTris[tri]->setOuterTet(pWmVols[c]);
                    // Add triangle to WmVols' table of neighbouring triangles.
                    pWmVols[c]->setNextTri(pTris[tri]);
                }
            }
        }
    }

    // All tets and tris that belong to some comp or patch have been created
    // locally- now we can connect them locally
    // NOTE: currently if a tetrahedron's neighbour belongs to a different
    // comp they do not talk to each other (see Tet::setNextTet())
    //

    AssertLog(ntets == pTets.size());
    // pTets member size of all tets in geometry, but may not be filled with
    // local tets if they have not been added to a compartment
    for (auto t: tetrahedron_global_id::range(ntets)) {
        if (pTets[t] == nullptr) {
            continue;
        }

        for (uint j = 0; j < 4; ++j) {
            const auto tet = pTets[t]->tet(j);
            if (tet.valid() && pTets[tet] != nullptr) {
                pTets[t]->setNextTet(j, pTets[tet]);
            }
        }
        // Not setting Tet triangles at this point- only want to set
        // for surface triangles
    }
    AssertLog(ntris == pTris.size());

    for (auto t: triangle_global_id::range(ntris)) {
        // Looping over all possible tris, but only some have been added to a patch
        if (pTris[t] == nullptr) {
            continue;
        }

        for (uint j = 0; j < 3; ++j) {
            const auto tri = pTris[t]->tri(j);
            if (tri.valid() && pTris[tri] != nullptr) {
                pTris[t]->setNextTri(j, pTris[tri]);
            }
        }

        // By convention, triangles in a patch should have an inner tetrahedron defined
        // (neighbouring tets 'flipped' if necessary in Tetmesh)
        // but not necessarily an outer tet
        // 17/3/10- actually this is not the case any more with well-mixed compartments
        //
        const auto tetinner = pTris[t]->tet(0);
        const auto tetouter = pTris[t]->tet(1);

        // Now inside and outside tetrahedrons may be normal tetrahedrons, which
        // means compartment is a mesh compartment, or wmvols describing a
        // well-mixed compartment with multiple triangle connections.


        if (tetinner.valid()) {
            // NEW FOR THIS VERSION: Tris store index of inner and outer tet (outer may not exist if
            // on surface) but tets may not belong to a compartment, even inner tets now since they
            // may be well-mixed compartments
            //
            if (pTets[tetinner] != nullptr) {
                // A triangle may already have an inner tet defined as a well-mixed
                // volume, but that should not be the case here:
                AssertLog(pTris[t]->iTet() == nullptr);

                pTris[t]->setInnerTet(pTets[tetinner]);
                // Now add this triangle to inner tet's list of neighbours
                for (uint i = 0; i <= 4; ++i) {
                    // include assert for debugging purposes and remove
                    // once this is tested
                    AssertLog(i < 4);  //////////
                    // check if there is already a neighbouring tet or tri
                    // In theory if there is a tri to add, the tet should
                    // have less than 4 neighbouring tets added because
                    // a neighbouring tet(s) is in a different compartment

                    // THIS IS NOT THE CASE ANYMORE: tets in different compartments can be
                    // neighbours so as to allow for diffusion boundaries

                    //     Also check tris because in some cases a surface tet
                    // may have more than 1 neighbouring tri
                    // NOTE: The order here will end up being different to the
                    // neighbour order at the Tetmesh level

                    // Now with diffusion boundaries, meaning tets can have neighbours that
                    // are in different comps, we must check the compartment
                    auto* tet_in = pTets[tetinner];
                    if (tet_in->nextTet(i) != nullptr &&
                        tet_in->compdef() == tet_in->nextTet(i)->compdef()) {
                        continue;
                    }

                    if (tet_in->nextTri(i) != nullptr) {
                        continue;
                    }
                    tet_in->setNextTri(i, pTris[t]);
                    break;
                }
            }
        }

        // DEBUG 18/03/09:
        // Now correct check, previously didn't allow for tet index == 0
        if (tetouter.valid()) {
            if (pTets[tetouter] != nullptr) {
                // A triangle may already have an inner tet defined as a well-mixed
                // volume, but that should not be the case here:
                AssertLog(pTris[t]->oTet() == nullptr);

                pTris[t]->setOuterTet(pTets[tetouter]);
                // Add this triangle to outer tet's list of neighbours
                for (uint i = 0; i <= 4; ++i) {
                    AssertLog(i < 4);

                    // See above in that tets now store tets from different comps
                    steps::tetexact::Tet* tet_out = pTets[tetouter];

                    if (tet_out->nextTet(i) != nullptr &&
                        tet_out->compdef() == tet_out->nextTet(i)->compdef()) {
                        continue;
                    }

                    if (tet_out->nextTri(i) != nullptr) {
                        continue;
                    }
                    tet_out->setNextTri(i, pTris[t]);
                    break;
                }
            }
        }
    }

    // Now loop over the diffusion boundaries:
    // 1) get all the triangles and get the two tetrahedrons
    // 2) figure out which direction is the direction for a tetrahedron
    // 3) add the tetrahedron and the direction to local object

    // This is here because we need all tets to have been assigned correctly
    // to compartments. Check every one and set the compA and compB for the db
    const auto ndiffbnds = pDiffBoundaries.size();
    AssertLog(ndiffbnds == pMesh->_countDiffBoundaries());

    for (uint db = 0; db < ndiffbnds; ++db) {
        auto* localdiffb = pDiffBoundaries[db];

        auto compAidx = localdiffb->def()->compa();
        auto compBidx = localdiffb->def()->compb();
        auto& compAdef = statedef().compdef(compAidx);
        auto& compBdef = statedef().compdef(compBidx);

        for (auto dbtri: localdiffb->def()->tris()) {
            const auto tri_tets = pMesh->_getTriTetNeighb(dbtri);

            auto tetAidx = tri_tets[0];
            auto tetBidx = tri_tets[1];

            AssertLog(tetAidx.valid() && tetBidx.valid());

            const auto* tetA = pTets[tetAidx];
            const auto* tetB = pTets[tetBidx];
            AssertLog(tetA != nullptr && tetB != nullptr);

            const auto tetA_cdef = tetA->compdef();
            const auto tetB_cdef = tetB->compdef();

            if (tetA_cdef != &compAdef) {
                AssertLog(tetB_cdef == &compAdef);
                AssertLog(tetA_cdef == &compBdef);
            } else {
                AssertLog(tetB_cdef == &compBdef);
                AssertLog(tetA_cdef == &compAdef);
            }

            // Ok, checks over, lets get down to business
            int direction_idx_a = -1;
            int direction_idx_b = -1;

            const auto& tetA_tris = pMesh->_getTetTriNeighb(tetAidx);
            const auto& tetB_tris = pMesh->_getTetTriNeighb(tetBidx);

            for (uint i = 0; i < 4; ++i) {
                if (tetA_tris[i] == dbtri) {
                    AssertLog(direction_idx_a == -1);
                    direction_idx_a = i;
                }
                if (tetB_tris[i] == dbtri) {
                    AssertLog(direction_idx_b == -1);
                    direction_idx_b = i;
                }
            }
            AssertLog(direction_idx_a != -1);
            AssertLog(direction_idx_b != -1);

            // Set the tetrahedron and direction to the Diff Boundary object
            localdiffb->setTetDirection(tetAidx, static_cast<unsigned int>(direction_idx_a));
            localdiffb->setTetDirection(tetBidx, static_cast<unsigned int>(direction_idx_b));
        }

        localdiffb->setComps(_comp(compAidx), _comp(compBidx));

        // Before the kprocs are set up ( in _setup) the tetrahedrons need to know the diffusion
        // boundary direction, so let's do it here  - the diff bounday has had all
        // tetrahedrons added

        // Might as well copy the vectors because we need to index through
        auto const& tets = localdiffb->getTets();
        auto const& tets_direction = localdiffb->getTetDirection();

        ntets = tets.size();
        AssertLog(ntets <= pTets.size());
        AssertLog(tets_direction.size() == ntets);

        for (auto t = 0u; t < ntets; ++t) {
            pTets[tets[t]]->setDiffBndDirection(tets_direction[t]);
        }
    }


    // Now loop over the surface diffusion boundaries:
    // 1) get all the bars and get the two triangles
    // 2) figure out which direction is the direction for a triangle
    // 3) add the triangle and the direction to local object

    // This is here because we need all tris to have been assigned correctly
    // to patches. Check every one and set the patchA and patchB for the db
    auto nsdiffbnds = pSDiffBoundaries.size();
    AssertLog(nsdiffbnds == pMesh->_countSDiffBoundaries());

    for (uint sdb = 0; sdb < nsdiffbnds; ++sdb) {
        steps::tetexact::SDiffBoundary* localsdiffb = pSDiffBoundaries[sdb];

        solver::patch_global_id patchAidx = localsdiffb->def()->patcha();
        solver::patch_global_id patchBidx = localsdiffb->def()->patchb();
        const auto& patchAdef = statedef().patchdef(patchAidx);
        const auto& patchBdef = statedef().patchdef(patchBidx);

        for (auto sdbbar: localsdiffb->def()->bars()) {
            const auto bar_tris = pMesh->_getBarTriNeighb(bar_id_t(sdbbar));

            const auto triAidx = bar_tris[0];
            const auto triBidx = bar_tris[1];
            AssertLog(triAidx.valid() && triBidx.valid());

            const auto* triA = pTris[triAidx];
            const auto* triB = pTris[triBidx];
            AssertLog(triA != nullptr && triB != nullptr);

            const auto* triA_pdef = triA->patchdef();
            const auto* triB_pdef = triB->patchdef();
            AssertLog(triA_pdef != nullptr);
            AssertLog(triB_pdef != nullptr);

            if (triA_pdef != &patchAdef) {
                AssertLog(triB_pdef == &patchAdef);
                AssertLog(triA_pdef == &patchBdef);
            } else {
                AssertLog(triB_pdef == &patchBdef);
                AssertLog(triA_pdef == &patchAdef);
            }

            // Ok, checks over, lets get down to business
            int direction_idx_a = -1;
            int direction_idx_b = -1;

            const auto triA_bars = pMesh->_getTriBars(triAidx);
            const auto triB_bars = pMesh->_getTriBars(triBidx);

            for (uint i = 0; i < 3; ++i) {
                if (triA_bars[i] == sdbbar) {
                    AssertLog(direction_idx_a == -1);
                    direction_idx_a = i;
                }
                if (triB_bars[i] == sdbbar) {
                    AssertLog(direction_idx_b == -1);
                    direction_idx_b = i;
                }
            }
            AssertLog(direction_idx_a != -1);
            AssertLog(direction_idx_b != -1);

            // Set the tetrahedron and direction to the Diff Boundary object
            localsdiffb->setTriDirection(triAidx, static_cast<unsigned int>(direction_idx_a));
            localsdiffb->setTriDirection(triBidx, static_cast<unsigned int>(direction_idx_b));
        }

        localsdiffb->setPatches(_patch(patchAidx), _patch(patchBidx));


        // Might as well copy the vectors because we need to index through
        auto const& tris = localsdiffb->getTris();
        auto const& tris_direction = localsdiffb->getTriDirection();

        ntris = tris.size();
        AssertLog(ntris <= pTris.size());
        AssertLog(tris_direction.size() == ntris);

        for (uint t = 0; t < ntris; ++t) {
            pTris[tris[t]]->setSDiffBndDirection(tris_direction[t]);
        }
    }

    for (auto const& t: pTets) {
        if (t != nullptr) {
            t->setupKProcs(this);
        }
    }

    for (auto const& wmv: pWmVols) {
        if (wmv != nullptr) {
            wmv->setupKProcs(this);
        }
    }

    for (auto const& t: pTris) {
        if (t != nullptr) {
            t->setupKProcs(this, efflag());
        }
    }

    // Resolve all dependencies
    for (auto const& t: pTets) {
        // DEBUG: vector holds all possible tetrahedrons,
        // but they have not necessarily been added to a compartment.
        if (t == nullptr) {
            continue;
        }
        for (auto const& k: t->kprocs()) {
            k->setupDeps();
        }
    }

    for (auto const& wmv: pWmVols) {
        // Vector allows for all compartments to be well-mixed, so
        // hold null-pointer for mesh compartments
        if (wmv == nullptr) {
            continue;
        }
        for (auto const& k: wmv->kprocs()) {
            k->setupDeps();
        }
    }

    for (auto const& t: pTris) {
        // DEBUG: vector holds all possible triangles, but
        // only patch triangles are filled
        if (t == nullptr) {
            continue;
        }
        for (auto const& k: t->kprocs()) {
            k->setupDeps();
        }
    }

    // Create EField structures if EField is to be calculated
    if (efflag()) {
        _setupEField();
    }

    nEntries = pKProcs.size();

    // force update on zero order reactions
    _update();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setupEField() {
    using math::point3d;

    CLOG(INFO, "general_log") << "setupEfield" << std::endl;

    //// Note to self: for now roughly following flow from original code in sim/controller.py.
    //// code for setting up a mesh was in func_tetmesh constructor and called functions.

    AssertLog(efflag());

    switch (pEFoption) {
    case EF_DEFAULT:
    case EF_DV_BDSYS:
        pEField = solver::efield::make_EField<solver::efield::dVSolverBanded>();
        break;
    default:
        ArgErrLog("Unsupported E-Field solver.");
    }

    // Give temperature a default value of 20c
    pTemp = 293.15;

    auto nmembs = mesh()._countMembs();

    if (nmembs != 1) {
        std::ostringstream os;
        os << "Membrane potential solver currently supports only one ";
        os << "membrane description object.";
        ArgErrLog(os.str());
    }

    tetmesh::Memb* memb = mesh()._getMemb(0);
    AssertLog(memb != nullptr);

    // TODO: Decide what checks are needed for the membrane and implement them here

    pEFNTets = memb->countVolTets();
    pEFNTris = memb->countTris();
    pEFNVerts = memb->countVerts();

    pEFTets.resize(neftets() * 4);

    // All the triangles we will count here are the real membrane triangles,
    // virtual triangles will not require a capacitance.
    pEFTris.resize(neftris() * 3);

    pEFVerts.resize(nefverts() * 3);

    auto nverts = mesh().countVertices();
    auto ntris = mesh().countTris();
    auto ntets = mesh().countTets();

    pEFVert_GtoL.container().resize(nverts);
    pEFTri_GtoL.container().resize(ntris);
    pEFTet_GtoL.container().resize(ntets);
    pEFTri_LtoG.container().resize(neftris());

    // Copy the data to local structures.

    auto const& membverts = memb->_getAllVertIndices();
    AssertLog(membverts.size() == nefverts());
    for (uint efv = 0; efv < nefverts(); ++efv) {
        auto vertidx = membverts[efv];
        point3d verttemp = mesh()._getVertex(vertidx);
        uint efv2 = efv * 3;

        // CONVERTING TO MICRONS HERE. EFIELD OBJECT WILL NOT PERFORM THIS CONVERSION
        verttemp *= 1.0e6;
        pEFVerts[efv2] = verttemp[0];
        pEFVerts[efv2 + 1] = verttemp[1];
        pEFVerts[efv2 + 2] = verttemp[2];

        pEFVert_GtoL[vertidx] = vertex_id_t(efv);
    }

    const auto& membtets = memb->_getAllVolTetIndices();
    AssertLog(membtets.size() == neftets());
    for (uint eft = 0; eft < neftets(); ++eft) {
        auto tetidx = membtets[eft];
        const auto tettemp = mesh()._getTet(tetidx);
        uint eft2 = eft * 4;

        // Convert to indices used by EField object
        const auto tv0 = pEFVert_GtoL[tettemp[0]];
        const auto tv1 = pEFVert_GtoL[tettemp[1]];
        const auto tv2 = pEFVert_GtoL[tettemp[2]];
        const auto tv3 = pEFVert_GtoL[tettemp[3]];
        if (tv0.unknown() || tv1.unknown() || tv2.unknown() || tv3.unknown()) {
            std::ostringstream os;
            os << "Failed to create EField structures.";
            ProgErrLog(os.str());
        }

        pEFTets[eft2] = tv0;
        pEFTets[eft2 + 1] = tv1;
        pEFTets[eft2 + 2] = tv2;
        pEFTets[eft2 + 3] = tv3;

        pEFTet_GtoL[tetidx] = tetrahedron_local_id(eft);
    }

    auto const& membtris = memb->_getAllTriIndices();
    AssertLog(membtris.size() == neftris());

    pEFTris_vec.resize(neftris());

    for (uint eft = 0; eft < neftris(); ++eft) {
        auto triidx = membtris[eft];
        const auto tritemp = mesh()._getTri(triidx);
        uint eft2 = eft * 3;

        // Convert to indices used by EField object
        const auto tv0 = pEFVert_GtoL[tritemp[0]];
        const auto tv1 = pEFVert_GtoL[tritemp[1]];
        const auto tv2 = pEFVert_GtoL[tritemp[2]];
        if (tv0.unknown() || tv1.unknown() || tv2.unknown()) {
            std::ostringstream os;
            os << "Failed to create EField structures.";
            ProgErrLog(os.str());
        }

        pEFTris[eft2] = tv0;
        pEFTris[eft2 + 1] = tv1;
        pEFTris[eft2 + 2] = tv2;

        pEFTri_GtoL[triidx] = triangle_local_id(eft);
        pEFTri_LtoG[triangle_local_id(eft)] = triidx;

        // This is added now for quicker iteration during run()
        // Extremely important for larger meshes, orders of magnitude times faster
        pEFTris_vec[eft] = pTris[triidx];
    }

    CLOG(INFO, "general_log") << "Initting mesh with:" << std::endl;
    CLOG(INFO, "general_log") << "Number of EF verts:" << nefverts() << std::endl
                              << "Number of EF tris:" << neftris() << std::endl
                              << "Number of EF tets:" << neftets() << std::endl;


    pEField->initMesh(pEFVerts,
                      pEFTris,
                      pEFTets,
                      memb->_getOpt_method(),
                      memb->_getOpt_file_name(),
                      memb->_getSearch_percent());
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::saveMembOpt(std::string const& opt_file_name) {
    if (!efflag()) {
        std::ostringstream os;
        os << "saveMembOpt method only available if running EField ";
        ArgErrLog(os.str());
    }

    pEField->saveOptimal(opt_file_name);
}

////////////////////////////////////////////////////////////////////////////////

// 'Safe' global to local index translation methods that throw on error.

inline solver::spec_local_id Tetexact::specG2L_or_throw(Comp* comp,
                                                        solver::spec_global_id gidx) const {
    AssertLog(gidx < statedef().countSpecs());
    solver::spec_local_id lidx = comp->def()->specG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("species undefined in compartment");
    }
    return lidx;
}

inline solver::spec_local_id Tetexact::specG2L_or_throw(Patch* patch,
                                                        solver::spec_global_id gidx) const {
    AssertLog(gidx < statedef().countSpecs());
    solver::spec_local_id lidx = patch->def()->specG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("species undefined in patch");
    }
    return lidx;
}

inline solver::reac_local_id Tetexact::reacG2L_or_throw(Comp* comp,
                                                        solver::reac_global_id gidx) const {
    AssertLog(gidx < statedef().countReacs());
    solver::reac_local_id lidx = comp->def()->reacG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("reaction undefined in compartment");
    }
    return lidx;
}

inline solver::sreac_local_id Tetexact::sreacG2L_or_throw(Patch* patch,
                                                          solver::sreac_global_id gidx) const {
    AssertLog(gidx < statedef().countSReacs());
    solver::sreac_local_id lidx = patch->def()->sreacG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("surface reaction undefined in patch");
    }
    return lidx;
}

inline solver::diff_local_id Tetexact::diffG2L_or_throw(Comp* comp,
                                                        solver::diff_global_id gidx) const {
    AssertLog(gidx < statedef().countDiffs());
    solver::diff_local_id lidx = comp->def()->diffG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("diffusion rule undefined in compartment");
    }
    return lidx;
}

inline solver::surfdiff_local_id Tetexact::sdiffG2L_or_throw(
    Patch* patch,
    solver::surfdiff_global_id gidx) const {
    AssertLog(gidx < statedef().countSurfDiffs());
    solver::surfdiff_local_id lidx = patch->def()->surfdiffG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("diffusion rule undefined in patch");
    }
    return lidx;
}

inline solver::vdepsreac_local_id Tetexact::vdepsreacG2L_or_throw(
    Patch* patch,
    solver::vdepsreac_global_id gidx) const {
    AssertLog(gidx < statedef().countVDepSReacs());
    solver::vdepsreac_local_id lidx = patch->def()->vdepsreacG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("voltage-dependent surface reation undefined in patch");
    }
    return lidx;
}

////////////////////////////////////////////////////////////////////////////////

solver::comp_global_id Tetexact::_addComp(solver::Compdef* cdef) {
    auto comp = new Comp(cdef);
    AssertLog(comp != nullptr);
    auto compidx = pComps.size();
    pComps.container().push_back(comp);
    return solver::comp_global_id(compidx);
}

////////////////////////////////////////////////////////////////////////////////

solver::patch_global_id Tetexact::_addPatch(solver::Patchdef* pdef) {
    auto patch = new Patch(pdef);
    AssertLog(patch != nullptr);
    auto patchidx = pPatches.size();
    pPatches.container().push_back(patch);
    return solver::patch_global_id(patchidx);
}

////////////////////////////////////////////////////////////////////////////////

solver::diffboundary_global_id Tetexact::_addDiffBoundary(solver::DiffBoundarydef* dbdef) {
    auto diffb = new DiffBoundary(dbdef);
    AssertLog(diffb != nullptr);
    auto dbidx = static_cast<uint>(pDiffBoundaries.size());
    pDiffBoundaries.push_back(diffb);
    return solver::diffboundary_global_id(dbidx);
}

////////////////////////////////////////////////////////////////////////////////

solver::sdiffboundary_global_id Tetexact::_addSDiffBoundary(solver::SDiffBoundarydef* sdbdef) {
    auto sdiffb = new SDiffBoundary(sdbdef);
    AssertLog(sdiffb != nullptr);
    auto sdbidx = static_cast<uint>(pSDiffBoundaries.size());
    pSDiffBoundaries.push_back(sdiffb);
    return solver::sdiffboundary_global_id(sdbidx);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_addTet(tetrahedron_global_id tetidx,
                       steps::tetexact::Comp* comp,
                       double vol,
                       double a1,
                       double a2,
                       double a3,
                       double a4,
                       double d1,
                       double d2,
                       double d3,
                       double d4,
                       tetrahedron_global_id tet0,
                       tetrahedron_global_id tet1,
                       tetrahedron_global_id tet2,
                       tetrahedron_global_id tet3) {
    solver::Compdef* compdef = comp->def();
    auto localtet =
        new Tet(tetidx, compdef, vol, a1, a2, a3, a4, d1, d2, d3, d4, tet0, tet1, tet2, tet3);
    AssertLog(pTets.at(tetidx) == nullptr);
    pTets[tetidx] = localtet;
    comp->addTet(localtet);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_addWmVol(solver::comp_global_id cidx, steps::tetexact::Comp* comp, double vol) {
    solver::Compdef* compdef = comp->def();
    auto localtet = new WmVol(tetrahedron_global_id(cidx.get()), compdef, vol);
    AssertLog(cidx < pWmVols.size());
    pWmVols[cidx] = localtet;
    comp->addTet(localtet);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_addTri(triangle_global_id triidx,
                       steps::tetexact::Patch* patch,
                       double area,
                       double l0,
                       double l1,
                       double l2,
                       double d0,
                       double d1,
                       double d2,
                       tetrahedron_global_id tinner,
                       tetrahedron_global_id touter,
                       triangle_global_id tri0,
                       triangle_global_id tri1,
                       triangle_global_id tri2) {
    auto* patchdef = patch->def();
    const auto tri =
        new Tri(triidx, patchdef, area, l0, l1, l2, d0, d1, d2, tinner, touter, tri0, tri1, tri2);
    AssertLog(pTris.at(triidx) == nullptr);
    pTris[triidx] = tri;
    patch->addTri(tri);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::reset() {
    for (auto const& comp: pComps) {
        comp->reset();
    }

    for (auto const& patch: pPatches) {
        patch->reset();
    }

    if (efflag()) {
        pEField->setMembPotential(solver::membrane_global_id(0), DEFAULT_MEMB_POT);
    }

    for (auto const& tet: pTets) {
        if (tet != nullptr) {
            tet->reset();
        }
    }

    for (auto const& wmvol: pWmVols) {
        if (wmvol != nullptr) {
            wmvol->reset();
        }
    }

    for (auto const& t: pTris) {
        if (t != nullptr) {
            t->reset();
        }
    }

    for (auto const& group: nGroups) {
        free(group->indices);
        delete group;
    }
    nGroups.clear();

    for (auto const& group: pGroups) {
        free(group->indices);
        delete group;
    }
    pGroups.clear();

    pSum = 0.0;
    nSum = 0.0;
    pA0 = 0.0;

    _update();
    statedef().resetTime();
    statedef().resetNSteps();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::run(double endtime) {
    if (!efflag()) {
        if (endtime < statedef().time()) {
            std::ostringstream os;
            os << "Endtime is before current simulation time";
            ArgErrLog(os.str());
        }
        Instrumentor::phase_begin("run -> rd");
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
        Instrumentor::phase_end("run -> rd");
    } else if (efflag()) {
        // Run the simulation, including the EField calculation.
        // This loop will assume that the SSA dt is sufficiently small so
        // that a number of SSA events execute between every EField calculation.
        // A warning message will be printed if the SSA dt is larger than the EField dt.
        // The EField dt is actually a MAXIMUM dt- the actual time for the EField
        // calculation will be exact with respect to the last event time in the
        // SSA before reaching the EField dt.
        while (statedef().time() < endtime) {
            if (util::almost_equal(statedef().time(), endtime)) {
                statedef().setTime(endtime);
                break;
            }
            Instrumentor::phase_begin("run -> with efield -> rd");
            double starttime = statedef().time();
            // The zero propensity
            double a0 = getA0();
            // We need a bool to check if the SSA contains no possible events. In
            // this rare case (continue to) execute the EField calculation to the endtime.
            bool ssa_on = true;
            double ssa_dt = 0.0;
            if (a0 != 0.0) {
                ssa_dt = rng()->getExp(a0);
            } else {
                (ssa_on = false);
            }
            // Set the actual efield dt. This value will take a maximum pEFDT.
            double ef_dt = 0.0;

            double maxDt = std::min(endtime - starttime, pEFDT);

            while (ssa_on && (ef_dt + ssa_dt) < maxDt) {
                KProc* kp = _getNext();
                if (kp == nullptr) {
                    break;
                }
                _executeStep(kp, ssa_dt);
                ef_dt += ssa_dt;

                a0 = getA0();
                if (a0 != 0.0) {
                    ssa_dt = rng()->getExp(a0);
                } else {
                    (ssa_on = false);
                }
            }
            AssertLog(ef_dt < maxDt);

            // It's possible that ef_dt is zero here: ssa_dt is large, or has become large.
            // In that case print a warning but continue, running the EField simulation for EFDT
            if (!ssa_on || ef_dt == 0.0) {
                std::ostringstream os;
                // os << "\nWARNING: SSA tau is larger than EField dt.";
                // CLOG(INFO, "general_log") << os << std::endl;
            }

            // Align to efdt or endtime
            if (endtime - starttime > pEFDT) {
                statedef().setTime(starttime + maxDt);
            } else {
                statedef().setTime(endtime);
            }

            // Now to perform the EField calculation. This means finding ohmic and GHK
            // currents from triangles during the ef_dt and applying these to the EField
            // object.

            triangle_local_id tlidx(0);
            double sttime = statedef().time();
            Instrumentor::phase_end("run -> with efield -> rd");
            Instrumentor::phase_begin("run -> with efield -> efield");
            for (auto const& eft: pEFTris_vec) {
                double v = pEField->getTriV(tlidx);
                double cur = eft->computeI(v, maxDt, sttime, efdt());
                pEField->setTriI(tlidx, cur);
                tlidx++;
            }

            pEField->advance(maxDt);
            Instrumentor::phase_end("run -> with efield -> efield");
            Instrumentor::phase_begin("run -> with efield -> kproc update");
            _update(pVdepKProcs.begin(), pVdepKProcs.end());
            Instrumentor::phase_end("run -> with efield -> kproc update");
        }
    }

    else {
        AssertLog(false);
    }
}

////////////////////////////////////////////////////////////////////////

void Tetexact::advance(double adv) {
    if (adv < 0.0) {
        std::ostringstream os;
        os << "Time to advance cannot be negative";
        ArgErrLog(os.str());
    }

    double endtime = statedef().time() + adv;
    run(endtime);
}

////////////////////////////////////////////////////////////////////////

void Tetexact::step() {
    if (efflag()) {
        std::ostringstream os;
        os << "Method not available with EField calculation.";
        ArgErrLog(os.str());
    }

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

////////////////////////////////////////////////////////////////////////////////

double Tetexact::getTime() const {
    return statedef().time();
}

////////////////////////////////////////////////////////////////////////

uint Tetexact::getNSteps() const {
    return statedef().nsteps();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::setTime(double time) {
    statedef().setTime(time);
}

////////////////////////////////////////////////////////////////////////

void Tetexact::setNSteps(uint nsteps) {
    statedef().setNSteps(nsteps);
}

////////////////////////////////////////////////////////////////////////

void Tetexact::setTemp(double t) {
    if (!efflag()) {
        std::ostringstream os;
        os << "\nWARNING: Temperature set in simulation without membrane ";
        os << "potential calculation will be ignored.\n";
        CLOG(INFO, "general_log") << os.str() << std::endl;
    }
    AssertLog(t >= 0.0);
    pTemp = t;
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getCompVol(solver::comp_global_id cidx) const {
    return _comp(cidx)->vol();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getCompSpecCount(solver::comp_global_id cidx, solver::spec_global_id sidx) const {
    Comp* comp = _comp(cidx);
    solver::spec_local_id slidx = specG2L_or_throw(comp, sidx);

    uint count = 0;
    for (auto& tet: comp->tets()) {
        count += tet->pools()[slidx];
    }

    return count;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setCompSpecCount(solver::comp_global_id cidx,
                                 solver::spec_global_id sidx,
                                 double n) {
    Comp* comp = _comp(cidx);
    solver::spec_local_id slidx = specG2L_or_throw(comp, sidx);

    // functions for distribution:
    auto set_count = [slidx](WmVol* tet, uint c) { tet->setCount(slidx, c); };
    auto inc_count = [slidx](WmVol* tet) { tet->incCount(slidx, 1); };
    auto weight = [](const WmVolPVecCI& tet) { return (*tet)->vol(); };

    util::distribute_quantity<uint>(n,
                                    comp->bgnTet(),
                                    comp->endTet(),
                                    weight,
                                    set_count,
                                    inc_count,
                                    *rng(),
                                    comp->def()->vol());

    for (auto& tet: comp->tets()) {
        _updateSpec(*tet);
    }
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getCompSpecAmount(solver::comp_global_id cidx,
                                    solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getCompSpecCount(cidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setCompSpecAmount(solver::comp_global_id cidx,
                                  solver::spec_global_id sidx,
                                  double a) {
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setCompSpecCount(cidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getCompSpecConc(solver::comp_global_id cidx, solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    return _getCompSpecCount(cidx, sidx) / (1.0e3 * _comp(cidx)->vol() * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setCompSpecConc(solver::comp_global_id cidx,
                                solver::spec_global_id sidx,
                                double c) {
    AssertLog(c >= 0.0);
    double count = c * (1.0e3 * _comp(cidx)->vol() * math::AVOGADRO);
    // the following method does all the necessary argument checking
    _setCompSpecCount(cidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getCompSpecClamped(solver::comp_global_id cidx, solver::spec_global_id sidx) const {
    Comp* comp = _comp(cidx);
    solver::spec_local_id slidx = specG2L_or_throw(comp, sidx);

    for (auto const& tet: comp->tets()) {
        if (!tet->clamped(slidx)) {
            return false;
        }
    }

    return true;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setCompSpecClamped(solver::comp_global_id cidx,
                                   solver::spec_global_id sidx,
                                   bool b) {
    Comp* comp = _comp(cidx);
    solver::spec_local_id slidx = specG2L_or_throw(comp, sidx);

    // Set the flag in def object, though this may not be necessary
    comp->def()->setClamped(slidx, b);
    for (auto const& tet: comp->tets()) {
        tet->setClamped(slidx, b);
    }
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getCompReacK(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
    Comp* comp = _comp(cidx);
    solver::reac_local_id lridx = reacG2L_or_throw(comp, ridx);

    // We're just returning the default value for this comp, individual
    // tets may have different Kcsts set individually
    return comp->def()->kcst(lridx);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setCompReacK(solver::comp_global_id cidx, solver::reac_global_id ridx, double kf) {
    AssertLog(kf >= 0.0);
    Comp* comp = _comp(cidx);
    solver::reac_local_id lridx = reacG2L_or_throw(comp, ridx);

    // First set the default value for the comp
    comp->def()->setKcst(lridx, kf);

    // Now update all tetrahedra in this comp
    for (auto& tet: comp->tets()) {
        tet->reac(lridx).setKcst(kf);
    }

    // Rates have changed
    _update();
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getCompReacActive(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
    Comp* comp = _comp(cidx);
    solver::reac_local_id lridx = reacG2L_or_throw(comp, ridx);

    for (auto const& tet: comp->tets()) {
        if (tet->reac(lridx).inactive()) {
            return false;
        }
    }

    return true;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setCompReacActive(solver::comp_global_id cidx,
                                  solver::reac_global_id ridx,
                                  bool a) {
    Comp* comp = _comp(cidx);
    solver::reac_local_id lridx = reacG2L_or_throw(comp, ridx);

    // Set the default value for the comp, though this is not entirely
    // necessary
    comp->def()->setActive(lridx, a);

    for (auto const& tet: comp->tets()) {
        tet->reac(lridx).setActive(a);
    }

    // It's cheaper to just recompute everything.
    _update();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getCompDiffD(solver::comp_global_id cidx, solver::diff_global_id didx) const {
    Comp* comp = _comp(cidx);
    solver::diff_local_id ldidx = diffG2L_or_throw(comp, didx);

    // We're just returning the default value for this comp, individual
    // tets may have different Dcsts set individually
    return comp->def()->dcst(ldidx);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setCompDiffD(solver::comp_global_id cidx, solver::diff_global_id didx, double dk) {
    AssertLog(dk >= 0.0);
    Comp* comp = _comp(cidx);
    solver::diff_local_id ldidx = diffG2L_or_throw(comp, didx);

    // First set the default value for the comp
    comp->def()->setDcst(ldidx, dk);

    // Now update all tets in this comp
    for (auto& wmvol: comp->tets()) {
        auto tet = dynamic_cast<Tet*>(wmvol);
        if (tet == nullptr) {
            ArgErrLog("cannot change diffusion constant in well-mixed compartment");
        }

        tet->diff(ldidx).setDcst(dk);
    }

    // Rates have changed
    _update();
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getCompDiffActive(solver::comp_global_id cidx, solver::diff_global_id didx) const {
    Comp* comp = _comp(cidx);
    solver::diff_local_id ldidx = diffG2L_or_throw(comp, didx);

    for (auto const& wmvol: comp->tets()) {
        auto tet = dynamic_cast<Tet*>(wmvol);
        if (tet == nullptr) {
            ArgErrLog("diffusion activation not defined in well-mixed compartment");
        }

        if (tet->diff(ldidx).inactive()) {
            return false;
        }
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setCompDiffActive(solver::comp_global_id cidx,
                                  solver::diff_global_id didx,
                                  bool act) {
    Comp* comp = _comp(cidx);
    solver::diff_local_id ldidx = diffG2L_or_throw(comp, didx);

    for (auto const& wmvol: comp->tets()) {
        auto tet = dynamic_cast<Tet*>(wmvol);
        if (tet == nullptr) {
            ArgErrLog("diffusion activation not defined in well-mixed compartment");
        }

        tet->diff(ldidx).setActive(act);
    }

    // It's cheaper to just recompute everything.
    _update();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getPatchArea(solver::patch_global_id pidx) const {
    return _patch(pidx)->area();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getPatchSpecCount(solver::patch_global_id pidx,
                                    solver::spec_global_id sidx) const {
    Patch* patch = _patch(pidx);
    solver::spec_local_id slidx = specG2L_or_throw(patch, sidx);

    uint count = 0;
    for (auto& tri: patch->tris()) {
        count += tri->pools()[slidx];
    }
    return count;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setPatchSpecCount(solver::patch_global_id pidx,
                                  solver::spec_global_id sidx,
                                  double n) {
    //////
    Patch* patch = _patch(pidx);
    solver::spec_local_id slidx = specG2L_or_throw(patch, sidx);

    // functions for distribution:
    auto set_count = [slidx](Tri* tri, uint c) { tri->setCount(slidx, c); };
    auto inc_count = [slidx](Tri* tri) { tri->incCount(slidx, 1); };
    auto weight = [](const TriPVecCI& tri) { return (*tri)->area(); };

    util::distribute_quantity<uint>(n,
                                    patch->bgnTri(),
                                    patch->endTri(),
                                    weight,
                                    set_count,
                                    inc_count,
                                    *rng(),
                                    patch->def()->area());

    for (auto& tri: patch->tris()) {
        _updateSpec(*tri);
    }
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getPatchSpecAmount(solver::patch_global_id pidx,
                                     solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getPatchSpecCount(pidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setPatchSpecAmount(solver::patch_global_id pidx,
                                   solver::spec_global_id sidx,
                                   double a) {
    AssertLog(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setPatchSpecCount(pidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getPatchSpecClamped(solver::patch_global_id pidx,
                                    solver::spec_global_id sidx) const {
    Patch* patch = _patch(pidx);
    solver::spec_local_id slidx = specG2L_or_throw(patch, sidx);

    for (auto& tri: patch->tris()) {
        if (!tri->clamped(slidx)) {
            return false;
        }
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setPatchSpecClamped(solver::patch_global_id pidx,
                                    solver::spec_global_id sidx,
                                    bool buf) {
    Patch* patch = _patch(pidx);
    solver::spec_local_id slidx = specG2L_or_throw(patch, sidx);

    // Set the flag in def object for consistency, though this is not
    // entirely necessary
    patch->def()->setClamped(slidx, buf);

    for (auto& tri: patch->tris()) {
        tri->setClamped(slidx, buf);
    }
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getPatchSReacK(solver::patch_global_id pidx, solver::sreac_global_id ridx) const {
    Patch* patch = _patch(pidx);
    solver::sreac_local_id lsridx = sreacG2L_or_throw(patch, ridx);

    // We're just returning the default value for this patch, individual
    // triangles may have different Kcsts set
    return patch->def()->kcst(lsridx);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setPatchSReacK(solver::patch_global_id pidx,
                               solver::sreac_global_id ridx,
                               double kf) {
    AssertLog(kf >= 0.0);
    Patch* patch = _patch(pidx);
    solver::sreac_local_id lsridx = sreacG2L_or_throw(patch, ridx);

    // First set the default values for this patch
    patch->def()->setKcst(lsridx, kf);

    // Now update all triangles in this patch
    for (auto& tri: patch->tris()) {
        tri->sreac(lsridx).setKcst(kf);
    }

    // Rates have changed
    _update();
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getPatchSReacActive(solver::patch_global_id pidx,
                                    solver::sreac_global_id ridx) const {
    Patch* patch = _patch(pidx);
    solver::sreac_local_id lsridx = sreacG2L_or_throw(patch, ridx);

    for (auto& tri: patch->tris()) {
        if (tri->sreac(lsridx).inactive()) {
            return false;
        }
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setDiffBoundarySpecDiffusionActive(solver::diffboundary_global_id dbidx,
                                                   solver::spec_global_id sidx,
                                                   bool act) {
    // Need to do two things:
    // 1) check if the species is defined in both compartments conencted
    // by the diffusion boundary
    // 2) loop over all tetrahedrons around the diff boundary and then the
    // diffusion rules and activate diffusion if the diffusion rule
    // relates to this species
    DiffBoundary* diffb = _diffboundary(dbidx);
    specG2L_or_throw(diffb->compA(), sidx);
    specG2L_or_throw(diffb->compB(), sidx);

    const auto& bdtets = diffb->getTets();
    const auto& bdtetsdir = diffb->getTetDirection();

    // Have to use indices rather than iterator because need access to the
    // tet direction
    const auto ntets = bdtets.size();

    for (auto bdt = 0u; bdt != ntets; ++bdt) {
        Tet* tet = pTets[bdtets[bdt]];
        auto direction = bdtetsdir[bdt];
        AssertLog(direction < 4);

        // Each diff kproc then has access to the species through it's defined parent
        uint ndiffs = tet->compdef()->countDiffs();
        for (auto d: solver::diff_local_id::range(ndiffs)) {
            Diff& diff = tet->diff(d);
            // sidx is the global species index; so is the lig() return from diffdef
            solver::spec_global_id specgidx = diff.def()->lig();
            if (specgidx == sidx) {
                diff.setDiffBndActive(direction, act);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getDiffBoundarySpecDiffusionActive(solver::diffboundary_global_id dbidx,
                                                   solver::spec_global_id sidx) const {
    DiffBoundary* diffb = _diffboundary(dbidx);
    specG2L_or_throw(diffb->compA(), sidx);
    specG2L_or_throw(diffb->compB(), sidx);

    const auto& bdtets = diffb->getTets();
    const auto& bdtetsdir = diffb->getTetDirection();

    // Have to use indices rather than iterator because need access to the
    // tet direction

    const auto ntets = static_cast<index_t>(bdtets.size());

    for (auto bdt: tetrahedron_global_id::range(ntets)) {
        Tet* tet = pTets[bdtets[bdt.get()]];
        auto direction = bdtetsdir[bdt.get()];
        AssertLog(direction < 4);

        // Each diff kproc then has access to the species through it's defined parent
        auto ndiffs = tet->compdef()->countDiffs();
        for (auto d: solver::diff_local_id::range(ndiffs)) {
            Diff& diff = tet->diff(d);
            // sidx is the global species index; so is the lig() return from diffdef
            solver::spec_global_id specgidx = diff.def()->lig();
            if (specgidx == sidx) {
                // Just need to check the first one
                return diff.getDiffBndActive(direction);
            }
        }
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setDiffBoundarySpecDcst(solver::diffboundary_global_id dbidx,
                                        solver::spec_global_id sidx,
                                        double dcst,
                                        solver::comp_global_id direction_comp) {
    auto* diffb = _diffboundary(dbidx);
    specG2L_or_throw(diffb->compA(), sidx);
    specG2L_or_throw(diffb->compB(), sidx);

    solver::Compdef* dirc_compdef = nullptr;
    if (direction_comp != std::numeric_limits<uint>::max()) {
        dirc_compdef = _comp(direction_comp)->def();
    }

    auto const& bdtets = diffb->getTets();
    auto const& bdtetsdir = diffb->getTetDirection();

    auto const ntets = bdtets.size();

    for (auto bdt = 0u; bdt != ntets; ++bdt) {
        Tet* tet = pTets[bdtets[bdt]];
        // if tet compdef equals to dirc_compdef,
        // it is the desination tet so diff should not be changed
        // NULL (bidirection) and source tet are both different
        // fromdirc_compdef
        if (dirc_compdef == tet->compdef()) {
            continue;
        }
        auto direction = bdtetsdir[bdt];
        AssertLog(direction < 4);

        // Each diff kproc then has access to the species through it's defined parent
        auto ndiffs = tet->compdef()->countDiffs();
        for (auto d: solver::diff_local_id::range(ndiffs)) {
            Diff& diff = tet->diff(d);
            // sidx is the global species index; so is the lig() return from diffdef
            auto specgidx = diff.def()->lig();
            if (specgidx == sidx) {
                // The following function will automatically activate diffusion
                // in this direction if necessary
                diff.setDirectionDcst(direction, dcst);
                _updateElement(diff);
            }
        }
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setSDiffBoundarySpecDiffusionActive(solver::sdiffboundary_global_id sdbidx,
                                                    solver::spec_global_id sidx,
                                                    bool act) {
    // Need to do two things:
    // 1) check if the species is defined in both patches connected
    // by the surface diffusion boundary
    // 2) loop over all triangles around the sdiff boundary and then the
    // diffusion rules and activate diffusion if the diffusion rule
    // relates to this species

    SDiffBoundary* sdiffb = _sdiffboundary(sdbidx);
    specG2L_or_throw(sdiffb->patchA(), sidx);
    specG2L_or_throw(sdiffb->patchB(), sidx);

    const auto& sbdtris = sdiffb->getTris();
    const auto& sbdtrisdir = sdiffb->getTriDirection();

    // Have to use indices rather than iterator because need access to the
    // tri direction
    const auto ntris = sbdtris.size();

    for (auto sbdt = 0u; sbdt != ntris; ++sbdt) {
        Tri* tri = pTris[sbdtris[sbdt]];
        auto direction = sbdtrisdir[sbdt];
        AssertLog(direction < 3);

        // Each sdiff kproc then has access to the species through its defined parent
        auto nsdiffs = tri->patchdef()->countSurfDiffs();
        for (auto sd: solver::surfdiff_local_id::range(nsdiffs)) {
            SDiff& sdiff = tri->sdiff(sd);
            // sidx is the global species index; so is the lig() return from diffdef
            auto specgidx = sdiff.def()->lig();
            if (specgidx == sidx) {
                sdiff.setSDiffBndActive(direction, act);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getSDiffBoundarySpecDiffusionActive(solver::sdiffboundary_global_id sdbidx,
                                                    solver::spec_global_id sidx) const {
    SDiffBoundary* sdiffb = _sdiffboundary(sdbidx);
    specG2L_or_throw(sdiffb->patchA(), sidx);
    specG2L_or_throw(sdiffb->patchB(), sidx);

    const auto& sbdtris = sdiffb->getTris();
    const auto& sbdtrisdir = sdiffb->getTriDirection();

    // Have to use indices rather than iterator because need access to the
    // tet direction

    const auto ntris = sbdtris.size();

    for (auto sbdt = 0u; sbdt != ntris; ++sbdt) {
        Tri* tri = pTris[sbdtris[sbdt]];
        auto direction = sbdtrisdir[sbdt];
        AssertLog(direction < 3);

        // Each sdiff kproc then has access to the species through its defined parent
        auto nsdiffs = tri->patchdef()->countSurfDiffs();
        for (auto sd: solver::surfdiff_local_id::range(nsdiffs)) {
            SDiff& sdiff = tri->sdiff(sd);
            // sidx is the global species index; so is the lig() return from diffdef
            auto specgidx = sdiff.def()->lig();
            if (specgidx == sidx) {
                // Just need to check the first one
                return sdiff.getSDiffBndActive(direction);
            }
        }
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setSDiffBoundarySpecDcst(solver::sdiffboundary_global_id sdbidx,
                                         solver::spec_global_id sidx,
                                         double dcst,
                                         solver::patch_global_id direction_patch) {
    SDiffBoundary* sdiffb = _sdiffboundary(sdbidx);
    specG2L_or_throw(sdiffb->patchA(), sidx);
    specG2L_or_throw(sdiffb->patchB(), sidx);

    solver::Patchdef* dirp_patchdef = nullptr;
    if (direction_patch != std::numeric_limits<uint>::max()) {
        dirp_patchdef = _patch(direction_patch)->def();
    }

    const auto& sbdtris = sdiffb->getTris();
    const auto& sbdtrisdir = sdiffb->getTriDirection();

    const auto ntris = sbdtris.size();

    for (auto sbdt = 0u; sbdt != ntris; ++sbdt) {
        Tri* tri = pTris[sbdtris[sbdt]];

        if (dirp_patchdef == tri->patchdef()) {
            continue;
        }
        auto direction = sbdtrisdir[sbdt];
        AssertLog(direction < 3);

        // Each diff kproc then has access to the species through its defined parent
        auto nsdiffs = tri->patchdef()->countSurfDiffs();
        for (auto sd: solver::surfdiff_local_id::range(nsdiffs)) {
            SDiff& sdiff = tri->sdiff(sd);
            // sidx is the global species index; so is the lig() return from diffdef
            auto specgidx = sdiff.def()->lig();
            if (specgidx == sidx) {
                // The following function will automatically activate diffusion
                // in this direction if necessary
                sdiff.setDirectionDcst(direction, dcst);
                _updateElement(sdiff);
            }
        }
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setPatchSReacActive(solver::patch_global_id pidx,
                                    solver::sreac_global_id ridx,
                                    bool a) {
    Patch* patch = _patch(pidx);
    solver::sreac_local_id lsridx = sreacG2L_or_throw(patch, ridx);

    // First set the flags in def object for consistency, though this is
    // not entirely necessary for this solver
    patch->def()->setActive(lsridx, a);

    for (auto const& tri: patch->tris()) {
        tri->sreac(lsridx).setActive(a);
    }

    // It's cheaper to just recompute everything.
    _update();
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getPatchVDepSReacActive(solver::patch_global_id pidx,
                                        solver::vdepsreac_global_id vsridx) const {
    Patch* patch = _patch(pidx);
    solver::vdepsreac_local_id lvsridx = vdepsreacG2L_or_throw(patch, vsridx);

    for (auto& tri: patch->tris()) {
        if (tri->vdepsreac(lvsridx).inactive()) {
            return false;
        }
    }

    return true;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setPatchVDepSReacActive(solver::patch_global_id pidx,
                                        solver::vdepsreac_global_id vsridx,
                                        bool a) {
    Patch* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::vdepsreac_local_id lvsridx = vdepsreacG2L_or_throw(patch, vsridx);


    // Not necessary and not possible to set the flags in def object
    for (auto& tri: patch->tris()) {
        tri->vdepsreac(lvsridx).setActive(a);
    }

    // It's cheaper to just recompute everything.
    _update();
}
////////////////////////////////////////////////////////////////////////////////

void Tetexact::addKProc(steps::tetexact::KProc* kp, bool Vdep) {
    AssertLog(kp != nullptr);

    SchedIDX nidx(static_cast<uint>(pKProcs.size()));  // because pKProcs.size() is ulong
    pKProcs.push_back(kp);
    if (Vdep) {
        pVdepKProcs.push_back(kp);
    }
    kp->setSchedIDX(nidx);
}

////////////////////////////////////////////////////////////////////////////////

steps::tetexact::KProc* Tetexact::_getNext() const {
    AssertLog(pA0 >= 0.0);
    // Quick check to see whether nothing is there.
    if (pA0 == 0.0) {
        return nullptr;
    }

    double selector = pA0 * rng()->getUnfII();

    double partial_sum = 0.0;

    const auto n_neg_groups = nGroups.size();
    const auto n_pos_groups = pGroups.size();

    for (uint i = 0; i < n_neg_groups; i++) {
        CRGroup* group = nGroups[i];
        if (group->size == 0) {
            continue;
        }

        if (selector > partial_sum + group->sum) {
            partial_sum += group->sum;

            continue;
        }

        double g_max = group->max;
        double random_rate = g_max * rng()->getUnfII();
        ;
        uint group_size = group->size;
        uint random_pos = rng()->get() % group_size;
        KProc* random_kp = group->indices[random_pos];

        while (random_kp->crData.rate <= random_rate) {
            random_rate = g_max * rng()->getUnfII();
            random_pos = rng()->get() % group_size;
            random_kp = group->indices[random_pos];
        }

        return random_kp;
    }


    for (uint i = 0; i < n_pos_groups; i++) {
        CRGroup* group = pGroups[i];
        if (group->size == 0) {
            continue;
        }

        if (selector > partial_sum + group->sum) {
            partial_sum += group->sum;
            continue;
        }

        double g_max = group->max;
        double random_rate = g_max * rng()->getUnfII();
        ;
        uint group_size = group->size;
        uint random_pos = rng()->get() % group_size;
        KProc* random_kp = group->indices[random_pos];


        while (random_kp->crData.rate <= random_rate) {
            random_rate = g_max * rng()->getUnfII();
            random_pos = rng()->get() % group_size;
            random_kp = group->indices[random_pos];
        }

        return random_kp;
    }

    // Precision rounding error force clean up
    // Force the search in the last non-empty group
    for (int i = n_pos_groups - 1; i >= 0; i--) {
        CRGroup* group = pGroups[i];
        if (group->size == 0) {
            continue;
        }

        double g_max = group->max;
        double random_rate = g_max * rng()->getUnfII();
        ;
        uint group_size = group->size;
        uint random_pos = rng()->get() % group_size;
        KProc* random_kp = group->indices[random_pos];


        while (random_kp->crData.rate <= random_rate) {
            random_rate = g_max * rng()->getUnfII();
            random_pos = rng()->get() % group_size;
            random_kp = group->indices[random_pos];
        }

        return random_kp;
    }

    for (int i = n_neg_groups - 1; i >= 0; i--) {
        CRGroup* group = nGroups[i];
        if (group->size == 0) {
            continue;
        }

        double g_max = group->max;
        double random_rate = g_max * rng()->getUnfII();
        ;
        uint group_size = group->size;
        uint random_pos = rng()->get() % group_size;
        KProc* random_kp = group->indices[random_pos];


        while (random_kp->crData.rate <= random_rate) {
            random_rate = g_max * rng()->getUnfII();
            random_pos = rng()->get() % group_size;
            random_kp = group->indices[random_pos];
        }

        return random_kp;
    }

    // Precision rounding error force clean up - Complete

    std::ostringstream os;
    os << "Cannot find any suitable entry.\n";
    os << "A0: " << std::setprecision(15) << pA0 << "\n";
    os << "Selector: " << std::setprecision(15) << selector << "\n";
    os << "Current Partial Sum: " << std::setprecision(15) << partial_sum << "\n";

    os << "Distribution of group sums\n";
    os << "Negative groups\n";

    for (uint i = 0; i < n_neg_groups; i++) {
        os << i << ": " << std::setprecision(15) << nGroups[i]->sum << "\n";
    }
    os << "Positive groups\n";
    for (uint i = 0; i < n_pos_groups; i++) {
        os << i << ": " << std::setprecision(15) << pGroups[i]->sum << "\n";
    }

    ProgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////
/*
void Tetexact::_reset()
{

}
*/
////////////////////////////////////////////////////////////////////////////////

void Tetexact::_executeStep(steps::tetexact::KProc* kp, double dt) {
    std::vector<KProc*> const& upd = kp->apply(rng(), dt, statedef().time());
    _update(upd.begin(), upd.end());
    statedef().incTime(dt);
    statedef().incNSteps(1);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_updateSpec(WmVol& tet) {
    KProcPSet updset;

    // Loop over tet.
    for (auto const& kproc: tet.kprocs()) {
        updset.insert(kproc);
    }

    for (auto const& tri: tet.nexttris()) {
        if (tri == nullptr) {
            continue;
        }
        for (auto const& kproc: tri->kprocs()) {
            updset.insert(kproc);
        }
    }

    // Send the list of kprocs that need to be updated to the schedule.
    _update(updset.begin(), updset.end());
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_updateSpec(Tri& tri) {
    _update(tri.kprocs().begin(), tri.kprocs().end());
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getCompReacH(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
    Comp* comp = _comp(cidx);
    solver::reac_local_id lridx = reacG2L_or_throw(comp, ridx);

    double h = 0.0;
    for (auto& tet: comp->tets()) {
        h += tet->reac(lridx).h();
    }
    return h;
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getCompReacC(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
    Comp* comp = _comp(cidx);
    solver::reac_local_id lridx = reacG2L_or_throw(comp, ridx);

    if (comp->tets().empty()) {
        return 0.0;
    }

    double c = 0.0;
    double v = 0.0;
    for (auto& tet: comp->tets()) {
        double v2 = tet->vol();
        c += tet->reac(lridx).c() * v2;
        v += v2;
    }
    AssertLog(v > 0.0);
    return c / v;
}

////////////////////////////////////////////////////////////////////////////////

long double Tetexact::_getCompReacA(solver::comp_global_id cidx,
                                    solver::reac_global_id ridx) const {
    Comp* comp = _comp(cidx);
    solver::reac_local_id lridx = reacG2L_or_throw(comp, ridx);

    long double a = 0.0L;
    for (auto& tet: comp->tets()) {
        a += static_cast<long double>(tet->reac(lridx).rate());
    }
    return a;
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long Tetexact::_getCompReacExtent(solver::comp_global_id cidx,
                                                solver::reac_global_id ridx) const {
    Comp* comp = _comp(cidx);
    solver::reac_local_id lridx = reacG2L_or_throw(comp, ridx);

    unsigned long long x = 0;
    for (auto& tet: comp->tets()) {
        x += tet->reac(lridx).getExtent();
    }
    return x;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_resetCompReacExtent(solver::comp_global_id cidx, solver::reac_global_id ridx) {
    Comp* comp = _comp(cidx);
    solver::reac_local_id lridx = reacG2L_or_throw(comp, ridx);

    // The 'local' Comp object has same index as solver::Compdef object
    for (auto& tet: comp->tets()) {
        tet->reac(lridx).resetExtent();
    }
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getPatchSReacH(solver::patch_global_id pidx, solver::sreac_global_id ridx) const {
    Patch* patch = _patch(pidx);
    solver::sreac_local_id lsridx = sreacG2L_or_throw(patch, ridx);

    double h = 0.0;
    for (auto& tri: patch->tris()) {
        h += tri->sreac(lsridx).h();
    }
    return h;
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getPatchSReacC(solver::patch_global_id pidx, solver::sreac_global_id ridx) const {
    Patch* patch = _patch(pidx);
    solver::sreac_local_id lsridx = sreacG2L_or_throw(patch, ridx);

    if (patch->tris().empty()) {
        return 0.0;
    }

    double c = 0.0;
    double a = 0.0;
    for (auto& tri: patch->tris()) {
        double a2 = tri->area();
        c += tri->sreac(lsridx).c() * a2;
        a += a2;
    }
    AssertLog(a > 0.0);
    return c / a;
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getPatchSReacA(solver::patch_global_id pidx, solver::sreac_global_id ridx) const {
    Patch* patch = _patch(pidx);
    solver::sreac_local_id lsridx = sreacG2L_or_throw(patch, ridx);

    double a = 0.0;
    for (auto& tri: patch->tris()) {
        a += tri->sreac(lsridx).rate();
    }
    return a;
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long Tetexact::_getPatchSReacExtent(solver::patch_global_id pidx,
                                                  solver::sreac_global_id ridx) const {
    Patch* patch = _patch(pidx);
    solver::sreac_local_id lsridx = sreacG2L_or_throw(patch, ridx);

    unsigned long long x = 0;
    for (auto& tri: patch->tris()) {
        x += tri->sreac(lsridx).getExtent();
    }
    return x;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_resetPatchSReacExtent(solver::patch_global_id pidx, solver::sreac_global_id ridx) {
    Patch* patch = _patch(pidx);
    solver::sreac_local_id lsridx = sreacG2L_or_throw(patch, ridx);

    for (auto& tri: patch->tris()) {
        tri->sreac(lsridx).resetExtent();
    }
}

////////////////////////////////////////////////////////////////////////////////

const Tet& Tetexact::_getTet(tetrahedron_global_id tidx) const {
    auto tet = pTets.at(tidx);
    if (tet == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
        ArgErrLog(os.str());
    }
    return *tet;
}

////////////////////////////////////////////////////////////////////////////////

Tet& Tetexact::_getTet(tetrahedron_global_id tidx) {
    auto tet = pTets.at(tidx);
    if (tet == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
        ArgErrLog(os.str());
    }
    return *tet;
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTetVol(tetrahedron_global_id tidx) const {
    return _getTet(tidx).vol();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTetVol(tetrahedron_global_id /*tidx*/, double /*vol*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getTetSpecDefined(tetrahedron_global_id tidx, solver::spec_global_id sidx) const {
    AssertLog(sidx < statedef().countSpecs());
    auto tet = pTets.at(tidx);

    if (tet == nullptr) {
        return false;
    }

    solver::spec_local_id slidx = tet->compdef()->specG2L(sidx);
    return slidx.valid();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTetSpecCount(tetrahedron_global_id tidx, solver::spec_global_id sidx) const {
    AssertLog(sidx < statedef().countSpecs());

    const auto& tet = _getTet(tidx);
    solver::spec_local_id slidx = tet.compdef()->specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    return tet.pools()[slidx];
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTetSpecCount(tetrahedron_global_id tidx, solver::spec_global_id sidx, double n) {
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(n >= 0.0);
    auto& tet = _getTet(tidx);

    if (n > std::numeric_limits<unsigned int>::max()) {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << std::numeric_limits<unsigned int>::max() << ").\n";
        ArgErrLog(os.str());
    }

    solver::spec_local_id slidx = tet.compdef()->specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    double n_int = std::floor(n);
    double n_frc = n - n_int;
    uint c = static_cast<uint>(n_int);
    if (n_frc > 0.0) {
        double rand01 = rng()->getUnfIE();
        if (rand01 < n_frc) {
            c++;
        }
    }

    // Tet object updates def level Comp object counts
    tet.setCount(slidx, c);
    _updateSpec(tet);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTetSpecAmount(tetrahedron_global_id tidx, solver::spec_global_id sidx) const {
    // following method does all necessary argument checking
    double count = _getTetSpecCount(tidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTetSpecAmount(tetrahedron_global_id tidx,
                                 solver::spec_global_id sidx,
                                 double m) {
    // convert amount in mols to number of molecules
    double m2 = m * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setTetSpecCount(tidx, sidx, m2);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTetSpecConc(tetrahedron_global_id tidx, solver::spec_global_id sidx) const {
    // following method does all necessary argument checking
    double count = _getTetSpecCount(tidx, sidx);
    Tet* tet = pTets[tidx];
    double vol = tet->vol();
    return count / (1.0e3 * vol * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTetSpecConc(tetrahedron_global_id tidx, solver::spec_global_id sidx, double c) {
    AssertLog(c >= 0.0);
    auto& tet = _getTet(tidx);

    double count = c * (1.0e3 * tet.vol() * math::AVOGADRO);
    // the following method does all the necessary argument checking
    _setTetSpecCount(tidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getTetSpecClamped(tetrahedron_global_id tidx, solver::spec_global_id sidx) const {
    AssertLog(sidx < statedef().countSpecs());
    auto& tet = _getTet(tidx);

    solver::spec_local_id slidx = tet.compdef()->specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    return tet.clamped(slidx);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTetSpecClamped(tetrahedron_global_id tidx,
                                  solver::spec_global_id sidx,
                                  bool buf) {
    AssertLog(sidx < statedef().countSpecs());
    auto& tet = _getTet(tidx);

    solver::spec_local_id slidx = tet.compdef()->specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    tet.setClamped(slidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTetReacK(tetrahedron_global_id tidx, solver::reac_global_id ridx) const {
    AssertLog(ridx < statedef().countReacs());
    auto& tet = _getTet(tidx);

    solver::reac_local_id lridx = tet.compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    return tet.reac(lridx).kcst();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTetReacK(tetrahedron_global_id tidx, solver::reac_global_id ridx, double kf) {
    AssertLog(ridx < statedef().countReacs());
    AssertLog(kf >= 0.0);
    auto& tet = _getTet(tidx);

    solver::reac_local_id lridx = tet.compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "\nReaction undefined in tetrahedron.";
        ArgErrLog(os.str());
    }

    tet.reac(lridx).setKcst(kf);

    _updateElement(tet.reac(lridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getTetReacActive(tetrahedron_global_id tidx, solver::reac_global_id ridx) const {
    AssertLog(ridx < statedef().countReacs());
    auto& tet = _getTet(tidx);

    solver::reac_local_id lridx = tet.compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    return !tet.reac(lridx).inactive();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTetReacActive(tetrahedron_global_id tidx,
                                 solver::reac_global_id ridx,
                                 bool act) {
    AssertLog(ridx < statedef().countReacs());
    auto& tet = _getTet(tidx);

    solver::reac_local_id lridx = tet.compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    tet.reac(lridx).setActive(act);

    _updateElement(tet.reac(lridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTetDiffD(tetrahedron_global_id tidx,
                              solver::diff_global_id didx,
                              tetrahedron_global_id direction_tet) const {
    AssertLog(didx < statedef().countDiffs());
    auto& tet = _getTet(tidx);

    solver::diff_local_id ldidx = tet.compdef()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    if (direction_tet.unknown()) {
        return tet.diff(ldidx).dcst();
    } else {
        int direction = tet.getTetDirection(direction_tet);
        if (direction == -1) {
            std::ostringstream os;
            os << "Tetrahedron " << direction_tet << " is not a neighbor of tetrahedron " << tidx
               << ".\n";
            ArgErrLog(os.str());
        }

        return tet.diff(ldidx).dcst(direction);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTetDiffD(tetrahedron_global_id tidx,
                            solver::diff_global_id didx,
                            double dk,
                            tetrahedron_global_id direction_tet) {
    AssertLog(didx < statedef().countDiffs());
    auto& tet = _getTet(tidx);

    solver::diff_local_id ldidx = tet.compdef()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    if (direction_tet.unknown()) {
        tet.diff(ldidx).setDcst(dk);
    } else {
        int direction = tet.getTetDirection(direction_tet);
        if (direction == -1) {
            std::ostringstream os;
            os << "Tetrahedron " << direction_tet << " is not a neighbor of tetrahedron " << tidx
               << ".\n";
            ArgErrLog(os.str());
        }

        tet.diff(ldidx).setDirectionDcst(direction, dk);
    }
    _updateElement(tet.diff(ldidx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getTetDiffActive(tetrahedron_global_id tidx, solver::diff_global_id didx) const {
    AssertLog(didx < statedef().countDiffs());
    auto& tet = _getTet(tidx);

    solver::diff_local_id ldidx = tet.compdef()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    return !tet.diff(ldidx).inactive();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTetDiffActive(tetrahedron_global_id tidx,
                                 solver::diff_global_id didx,
                                 bool act) {
    AssertLog(didx < statedef().countDiffs());

    auto& tet = _getTet(tidx);
    solver::diff_local_id ldidx = tet.compdef()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    tet.diff(ldidx).setActive(act);

    _updateElement(tet.diff(ldidx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTetReacH(tetrahedron_global_id tidx, solver::reac_global_id ridx) const {
    AssertLog(ridx < statedef().countReacs());

    auto& tet = _getTet(tidx);
    solver::reac_local_id lridx = tet.compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    return tet.reac(lridx).h();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTetReacC(tetrahedron_global_id tidx, solver::reac_global_id ridx) const {
    AssertLog(ridx < statedef().countReacs());

    auto& tet = _getTet(tidx);
    solver::reac_local_id lridx = tet.compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    return tet.reac(lridx).c();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTetReacA(tetrahedron_global_id tidx, solver::reac_global_id ridx) const {
    AssertLog(ridx < statedef().countReacs());

    auto& tet = _getTet(tidx);
    solver::reac_local_id lridx = tet.compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    return tet.reac(lridx).rate();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTetDiffA(tetrahedron_global_id tidx, solver::diff_global_id didx) const {
    AssertLog(didx < statedef().countDiffs());

    auto& tet = _getTet(tidx);
    solver::diff_local_id ldidx = tet.compdef()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    return tet.diff(ldidx).rate();
}

////////////////////////////////////////////////////////////////////////////////

const Tri& Tetexact::_getTri(triangle_global_id tidx) const {
    auto tri = pTris.at(tidx);
    if (tri == nullptr) {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        ArgErrLog(os.str());
    }
    return *tri;
}

////////////////////////////////////////////////////////////////////////////////

Tri& Tetexact::_getTri(triangle_global_id tidx) {
    auto tri = pTris.at(tidx);
    if (tri == nullptr) {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        ArgErrLog(os.str());
    }
    return *tri;
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriArea(triangle_global_id tidx) const {
    return _getTri(tidx).area();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTriArea(triangle_global_id /*tidx*/, double /*area*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getTriSpecDefined(triangle_global_id tidx, solver::spec_global_id sidx) const {
    AssertLog(sidx < statedef().countSpecs());
    auto tri = pTris.at(tidx);

    if (tri == nullptr) {
        return false;
    }

    solver::spec_local_id slidx = tri->patchdef()->specG2L(sidx);
    return slidx.valid();
}

////////////////////////////////////////////////////////////////////////////////


double Tetexact::_getTriSpecCount(triangle_global_id tidx, solver::spec_global_id sidx) const {
    AssertLog(sidx < statedef().countSpecs());

    auto& tri = _getTri(tidx);
    solver::spec_local_id slidx = tri.patchdef()->specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    return tri.pools()[slidx];
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTriSpecCount(triangle_global_id tidx, solver::spec_global_id sidx, double n) {
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(n >= 0.0);

    if (n > std::numeric_limits<unsigned int>::max()) {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << std::numeric_limits<unsigned int>::max() << ").\n";
        ArgErrLog(os.str());
    }

    auto& tri = _getTri(tidx);
    solver::spec_local_id slidx = tri.patchdef()->specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    double n_int = std::floor(n);
    double n_frc = n - n_int;
    uint c = static_cast<uint>(n_int);
    if (n_frc > 0.0) {
        double rand01 = rng()->getUnfIE();
        if (rand01 < n_frc) {
            c++;
        }
    }

    // Tri object updates counts in def level Comp object
    tri.setCount(slidx, c);
    _updateSpec(tri);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriSpecAmount(triangle_global_id tidx, solver::spec_global_id sidx) const {
    // following method does all necessary argument checking
    return _getTriSpecCount(tidx, sidx) / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTriSpecAmount(triangle_global_id tidx, solver::spec_global_id sidx, double m) {
    // the following method does all the necessary argument checking
    _setTriSpecCount(tidx, sidx, m * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////


bool Tetexact::_getTriSpecClamped(triangle_global_id tidx, solver::spec_global_id sidx) const {
    AssertLog(sidx < statedef().countSpecs());

    auto& tri = _getTri(tidx);
    solver::spec_local_id slidx = tri.patchdef()->specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    return tri.clamped(slidx);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTriSpecClamped(triangle_global_id tidx, solver::spec_global_id sidx, bool buf) {
    AssertLog(sidx < statedef().countSpecs());

    auto& tri = _getTri(tidx);
    solver::spec_local_id slidx = tri.patchdef()->specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    tri.setClamped(slidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriSReacK(triangle_global_id tidx, solver::sreac_global_id ridx) const {
    AssertLog(ridx < statedef().countSReacs());

    auto& tri = _getTri(tidx);
    solver::sreac_local_id lsridx = tri.patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    return tri.sreac(lsridx).kcst();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTriSReacK(triangle_global_id tidx, solver::sreac_global_id ridx, double kf) {
    AssertLog(ridx < statedef().countSReacs());

    auto& tri = _getTri(tidx);
    solver::sreac_local_id lsridx = tri.patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    tri.sreac(lsridx).setKcst(kf);

    _updateElement(tri.sreac(lsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getTriSReacActive(triangle_global_id tidx, solver::sreac_global_id ridx) const {
    AssertLog(ridx < statedef().countSReacs());

    auto& tri = _getTri(tidx);
    solver::sreac_local_id lsridx = tri.patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    return !tri.sreac(lsridx).inactive();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTriSReacActive(triangle_global_id tidx, solver::sreac_global_id ridx, bool act) {
    AssertLog(ridx < statedef().countSReacs());

    auto& tri = _getTri(tidx);
    solver::sreac_local_id lsridx = tri.patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    tri.sreac(lsridx).setActive(act);

    _updateElement(tri.sreac(lsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriSDiffD(triangle_global_id tidx,
                               solver::surfdiff_global_id didx,
                               triangle_global_id direction_tri) const {
    AssertLog(didx < statedef().countSurfDiffs());

    auto& tri = _getTri(tidx);
    solver::surfdiff_local_id ldidx = tri.patchdef()->surfdiffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    if (direction_tri.unknown()) {
        return tri.sdiff(ldidx).dcst();

    } else {
        int direction = tri.getTriDirection(direction_tri);
        if (direction == -1) {
            std::ostringstream os;
            os << "Triangle " << direction_tri << " is not a neighbor of triangle " << tidx
               << ".\n";
            ArgErrLog(os.str());
        }

        return tri.sdiff(ldidx).dcst(direction);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTriSDiffD(triangle_global_id tidx,
                             solver::surfdiff_global_id didx,
                             double dk,
                             triangle_global_id direction_tri) {
    AssertLog(didx < statedef().countSurfDiffs());

    auto& tri = _getTri(tidx);
    solver::surfdiff_local_id ldidx = tri.patchdef()->surfdiffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    if (direction_tri.unknown()) {
        tri.sdiff(ldidx).setDcst(dk);

    } else {
        int direction = tri.getTriDirection(direction_tri);
        if (direction == -1) {
            std::ostringstream os;
            os << "Triangle " << direction_tri << " is not a neighbor of triangle " << tidx
               << ".\n";
            ArgErrLog(os.str());
        }


        tri.sdiff(ldidx).setDirectionDcst(direction, dk);
    }
    _updateElement(tri.sdiff(ldidx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getTriVDepSReacActive(triangle_global_id tidx,
                                      solver::vdepsreac_global_id vsridx) const {
    AssertLog(vsridx < statedef().countVDepSReacs());

    auto& tri = _getTri(tidx);
    solver::vdepsreac_local_id lvsridx = tri.patchdef()->vdepsreacG2L(vsridx);
    if (lvsridx.unknown()) {
        std::ostringstream os;
        os << "Voltage-dependent surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    return !tri.vdepsreac(lvsridx).inactive();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTriVDepSReacActive(triangle_global_id tidx,
                                      solver::vdepsreac_global_id vsridx,
                                      bool act) {
    AssertLog(vsridx < statedef().countVDepSReacs());

    auto& tri = _getTri(tidx);
    solver::vdepsreac_local_id lvsridx = tri.patchdef()->vdepsreacG2L(vsridx);
    if (lvsridx.unknown()) {
        std::ostringstream os;
        os << "Voltage-dependent surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    tri.vdepsreac(lvsridx).setActive(act);

    _updateElement(tri.vdepsreac(lvsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriSReacH(triangle_global_id tidx, solver::sreac_global_id ridx) const {
    AssertLog(ridx < statedef().countSReacs());

    auto& tri = _getTri(tidx);
    solver::sreac_local_id lsridx = tri.patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    return tri.sreac(lsridx).h();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriSReacC(triangle_global_id tidx, solver::sreac_global_id ridx) const {
    AssertLog(ridx < statedef().countSReacs());

    auto& tri = _getTri(tidx);
    solver::sreac_local_id lsridx = tri.patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    return tri.sreac(lsridx).c();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriSReacA(triangle_global_id tidx, solver::sreac_global_id ridx) const {
    AssertLog(ridx < statedef().countSReacs());

    auto& tri = _getTri(tidx);
    solver::sreac_local_id lsridx = tri.patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    return tri.sreac(lsridx).rate();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::setEfieldDT(double efdt) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    if (efdt <= 0.0) {
        std::ostringstream os;
        os << "EField dt must be graeter than zero.";
        ArgErrLog(os.str());
    }
    pEFDT = efdt;
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTetV(tetrahedron_global_id tidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto loctidx = pEFTet_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    // EField object should convert value to base s.i. units
    return pEField->getTetV(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTetV(tetrahedron_global_id tidx, double v) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto loctidx = pEFTet_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    // EField object should convert to millivolts
    pEField->setTetV(loctidx, v);

    // Voltage-dependent rates may have changed
    _update();
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getTetVClamped(tetrahedron_global_id tidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto loctidx = pEFTet_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    return pEField->getTetVClamped(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTetVClamped(tetrahedron_global_id tidx, bool cl) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto loctidx = pEFTet_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    pEField->setTetVClamped(loctidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriV(triangle_global_id tidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }
    // EField object should convert value to base s.i. units
    return pEField->getTriV(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTriV(triangle_global_id tidx, double v) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    // EField object should convert to millivolts
    pEField->setTriV(loctidx, v);

    // Voltage-dependent rates may have changed
    _update();
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getTriVClamped(triangle_global_id tidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    return pEField->getTriVClamped(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTriVClamped(triangle_global_id tidx, bool cl) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    pEField->setTriVClamped(loctidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTriOhmicErev(triangle_global_id tidx,
                                solver::ohmiccurr_global_id ocgidx,
                                double erev) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }

    auto loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    Tri* tri = pTris[tidx];

    solver::ohmiccurr_local_id locidx = tri->patchdef()->ohmiccurrG2L(ocgidx);
    if (locidx.unknown()) {
        std::ostringstream os;
        os << "Ohmic current undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    tri->setOCerev(locidx, erev);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriOhmicErev(triangle_global_id tidx,
                                  solver::ohmiccurr_global_id ocgidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }

    auto loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    Tri* tri = pTris[tidx];

    solver::ohmiccurr_local_id locidx = tri->patchdef()->ohmiccurrG2L(ocgidx);
    if (locidx.unknown()) {
        std::ostringstream os;
        os << "Ohmic current undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    return tri->getOCerev(locidx);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriOhmicI(triangle_global_id tidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }

    auto& tri = _getTri(tidx);

    const auto loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    return tri.getOhmicI(pEField->getTriV(loctidx), efdt());
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriOhmicI(triangle_global_id tidx, solver::ohmiccurr_global_id ocidx) const {
    AssertLog(ocidx < statedef().countOhmicCurrs());

    auto& tri = _getTri(tidx);
    const auto loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    solver::ohmiccurr_local_id locidx = tri.patchdef()->ohmiccurrG2L(ocidx);
    if (locidx.unknown()) {
        std::ostringstream os;
        os << "Ohmic current undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    return tri.getOhmicI(locidx, pEField->getTriV(loctidx), efdt());
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriGHKI(triangle_global_id tidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }

    return _getTri(tidx).getGHKI();
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriGHKI(triangle_global_id tidx, solver::ghkcurr_global_id ghkidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }

    auto& tri = _getTri(tidx);

    solver::ghkcurr_local_id locidx = tri.patchdef()->ghkcurrG2L(ghkidx);
    if (locidx.unknown()) {
        std::ostringstream os;
        os << "GHK current undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    return tri.getGHKI(locidx);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriI(triangle_global_id tidx) const {
    return _getTriGHKI(tidx) + _getTriOhmicI(tidx);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getVertIClamp(vertex_id_t vidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto locvidx = pEFVert_GtoL[vidx];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }

    return pEField->getVertIClamp(locvidx);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setVertIClamp(vertex_id_t vidx, double cur) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto locvidx = pEFVert_GtoL[vidx];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }

    // EField object should convert to required units
    pEField->setVertIClamp(locvidx, cur);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getTriIClamp(triangle_global_id tidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    return pEField->getTriIClamp(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTriIClamp(triangle_global_id tidx, double cur) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    // EField object should convert to required units
    pEField->setTriIClamp(loctidx, cur);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setTriCapac(triangle_global_id tidx, double cap) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    // EField object should convert to required units
    pEField->setTriCapac(loctidx, cap);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::_getVertV(vertex_id_t vidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto locvidx = pEFVert_GtoL[vidx];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }
    // EField object should convert value to base s.i. units
    return pEField->getVertV(locvidx);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setVertV(vertex_id_t vidx, double v) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto locvidx = pEFVert_GtoL[vidx];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }
    // EField object should convert to millivolts
    pEField->setVertV(locvidx, v);

    // Voltage-dependent rates may have changed
    _update();
}

////////////////////////////////////////////////////////////////////////////////

bool Tetexact::_getVertVClamped(vertex_id_t vidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto locvidx = pEFVert_GtoL[vidx];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }

    return pEField->getVertVClamped(locvidx);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setVertVClamped(vertex_id_t vidx, bool cl) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    const auto locvidx = pEFVert_GtoL[vidx];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }
    // EField object should convert to millivolts
    pEField->setVertVClamped(locvidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setMembRes(solver::membrane_global_id midx, double ro, double vrev) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    if (ro <= 0.0) {
        std::ostringstream os;
        os << "Resistivity must be greater than zero.";
        ArgErrLog(os.str());
    }
    // EField object should convert to required units
    AssertLog(midx.get() == 0);
    pEField->setSurfaceResistivity(midx, ro, vrev);
}

std::pair<double, double> Tetexact::_getMembRes(solver::membrane_global_id midx) const {
    if (!efflag()) {
        ArgErrLog("Method not available: EField calculation not included in simulation.");
    }
    // EField object should convert to required units
    AssertLog(midx.get() == 0);
    return pEField->getSurfaceResistivity();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setMembPotential(solver::membrane_global_id midx, double v) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    // EField object should convert to millivolts
    AssertLog(midx.get() == 0);
    pEField->setMembPotential(midx, v);

    _update();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setMembCapac(solver::membrane_global_id midx, double cm) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    if (cm < 0.0) {
        std::ostringstream os;
        os << "Capacitance must be greater than or equal to zero.";
        ArgErrLog(os.str());
    }


    // EField object should convert to required units
    AssertLog(midx.get() == 0);
    pEField->setMembCapac(midx, cm);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_setMembVolRes(solver::membrane_global_id midx, double ro) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    if (ro < 0.0) {
        std::ostringstream os;
        os << "Resistivity must be greater than or equal to zero.";
        ArgErrLog(os.str());
    }
    // EField object should convert to required units
    AssertLog(midx.get() == 0);
    pEField->setMembVolRes(midx, ro);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::_updateElement(KProc& kp) {
    double new_rate = kp.rate(this);

    CRKProcData& data = kp.crData;
    double old_rate = data.rate;

    data.rate = new_rate;

    if (old_rate == new_rate) {
        return;
    }


    // new rate in positive groups
    if (new_rate >= 0.5) {
        // pow is the same
        int old_pow = data.pow;
        int new_pow;
        std::frexp(new_rate, &new_pow);

        if (old_pow == new_pow && data.recorded) {
            CRGroup* old_group = _getGroup(old_pow);

            old_group->sum += (new_rate - old_rate);
        }
        // pow is not the same
        else {
            data.pow = new_pow;


            if (data.recorded) {
                // remove old
                CRGroup* old_group = _getGroup(old_pow);
                (old_group->size)--;

                if (old_group->size == 0) {
                    old_group->sum = 0.0;
                } else {
                    old_group->sum -= old_rate;

                    KProc* last = old_group->indices[old_group->size];
                    old_group->indices[data.pos] = last;
                    last->crData.pos = data.pos;
                }
            }

            // add new
            if (static_cast<int>(pGroups.size()) <= new_pow) {
                _extendPGroups(new_pow + 1);
            }

            CRGroup* new_group = pGroups[new_pow];

            AssertLog(new_group != nullptr);
            if (new_group->size == new_group->capacity) {
                _extendGroup(new_group);
            }
            uint pos = new_group->size;
            new_group->indices[pos] = &kp;
            new_group->size++;
            new_group->sum += new_rate;
            data.pos = pos;
        }
        data.recorded = true;

    }
    // new rate in negative group
    else if (new_rate < 0.5 && new_rate > 1e-20) {
        int old_pow = data.pow;
        int new_pow;
        frexp(new_rate, &new_pow);

        if (old_pow == new_pow && data.recorded) {
            CRGroup* old_group = _getGroup(old_pow);

            old_group->sum += (new_rate - old_rate);
        }
        // pow is not the same
        else {
            data.pow = new_pow;

            if (data.recorded) {
                CRGroup* old_group = _getGroup(old_pow);
                (old_group->size)--;

                if (old_group->size == 0) {
                    old_group->sum = 0.0;
                } else {
                    old_group->sum -= old_rate;

                    KProc* last = old_group->indices[old_group->size];
                    old_group->indices[data.pos] = last;
                    last->crData.pos = data.pos;
                }
            }

            // add new

            if (static_cast<int>(nGroups.size()) <= -new_pow) {
                _extendNGroups(-new_pow + 1);
            }

            CRGroup* new_group = nGroups[-new_pow];

            if (new_group->size == new_group->capacity) {
                _extendGroup(new_group);
            }
            uint pos = new_group->size;
            new_group->indices[pos] = &kp;
            new_group->size++;
            new_group->sum += new_rate;
            data.pos = pos;
        }
        data.recorded = true;
    }

    else {
        if (data.recorded) {
            CRGroup* old_group = _getGroup(data.pow);

            // remove old
            old_group->size--;

            if (old_group->size == 0) {
                old_group->sum = 0.0;
            } else {
                old_group->sum -= old_rate;

                KProc* last = old_group->indices[old_group->size];
                old_group->indices[data.pos] = last;
                last->crData.pos = data.pos;
            }
        }
        data.recorded = false;
    }
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// ROI Data Access
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

std::vector<double> Tetexact::getBatchTetSpecCounts(const std::vector<index_t>& tets,
                                                    std::string const& s) const {
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;

    auto ntets = tets.size();
    solver::spec_global_id sgidx = statedef().getSpecIdx(s);
    std::vector<double> data(ntets, 0.0);

    for (uint t = 0; t < ntets; t++) {
        const tetrahedron_global_id tidx(tets[t]);

        if (tidx >= pTets.size()) {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            ArgErrLog(os.str());
        }

        if (pTets[tidx] == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        Tet* tet = pTets[tidx];
        solver::spec_local_id slidx = tet->compdef()->specG2L(sgidx);
        if (slidx.unknown()) {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }

        data[t] = tet->pools()[slidx];
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log") << "Species " << s
                                     << " has not been defined in the following tetrahedrons, fill "
                                        "in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
    return data;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> Tetexact::getBatchTriSpecCounts(const std::vector<index_t>& tris,
                                                    std::string const& s) const {
    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;


    auto ntris = tris.size();
    solver::spec_global_id sgidx = statedef().getSpecIdx(s);
    std::vector<double> data(ntris, 0.0);

    for (uint t = 0; t < ntris; t++) {
        const triangle_global_id tidx(tris[t]);
        auto tri = pTris.at(tidx);
        if (tri == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        solver::spec_local_id slidx = tri->patchdef()->specG2L(sgidx);
        if (slidx.unknown()) {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }

        data[t] = tri->pools()[slidx];
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log") << "The following triangles have not been assigned to a "
                                        "patch, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log") << "Species " << s
                                     << " has not been defined in the following triangles, fill in "
                                        "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::getBatchTetSpecCountsNP(const index_t* indices,
                                       size_t input_size,
                                       std::string const& s,
                                       double* counts,
                                       size_t output_size) const {
    if (input_size != output_size) {
        std::ostringstream os;
        os << "Error: output array (counts) size should be the same as input array (indices) "
              "size.\n";
        ArgErrLog(os.str());
    }

    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (uint t = 0; t < input_size; t++) {
        const tetrahedron_global_id tidx(indices[t]);

        if (tidx >= pTets.size()) {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            ArgErrLog(os.str());
        }

        if (pTets[tidx] == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        Tet* tet = pTets[tidx];
        solver::spec_local_id slidx = tet->compdef()->specG2L(sgidx);
        if (slidx.unknown()) {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }

        counts[t] = tet->pools()[slidx];
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log") << "Species " << s
                                     << " has not been defined in the following tetrahedrons, fill "
                                        "in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::getBatchTriSpecCountsNP(const index_t* indices,
                                       size_t input_size,
                                       std::string const& s,
                                       double* counts,
                                       size_t output_size) const {
    if (input_size != output_size) {
        std::ostringstream os;
        os << "Error: output array (counts) size should be the same as input array (indices) "
              "size.\n";
        ArgErrLog(os.str());
    }

    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;


    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (uint t = 0; t < input_size; t++) {
        const triangle_global_id tidx(indices[t]);
        auto tri = pTris.at(tidx);
        if (tri == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        solver::spec_local_id slidx = tri->patchdef()->specG2L(sgidx);
        if (slidx.unknown()) {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }

        counts[t] = tri->pools()[slidx];
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log") << "The following triangles have not been assigned to a "
                                        "patch, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log") << "Species " << s
                                     << " has not been defined in the following triangles, fill in "
                                        "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// ROI Data Access
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

std::vector<double> Tetexact::getROITetSpecCounts(const std::string& ROI_id,
                                                  std::string const& s) const {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    const auto size = roi->second.size();
    std::vector<double> data(size);

    getBatchTetSpecCountsNP(
        reinterpret_cast<const index_t*>(roi->second.data()), size, s, &data.front(), data.size());
    return data;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> Tetexact::getROITriSpecCounts(const std::string& ROI_id,
                                                  std::string const& s) const {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    const auto size = roi->second.size();
    std::vector<double> data(size);
    getBatchTriSpecCountsNP(
        reinterpret_cast<const index_t*>(roi->second.data()), size, s, &data.front(), data.size());
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::getROITetSpecCountsNP(const std::string& ROI_id,
                                     std::string const& s,
                                     double* counts,
                                     size_t output_size) const {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    getBatchTetSpecCountsNP(reinterpret_cast<const index_t*>(roi->second.data()),
                            roi->second.size(),
                            s,
                            counts,
                            output_size);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::getROITriSpecCountsNP(const std::string& ROI_id,
                                     std::string const& s,
                                     double* counts,
                                     size_t output_size) const {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    getBatchTriSpecCountsNP(reinterpret_cast<const index_t*>(roi->second.data()),
                            roi->second.size(),
                            s,
                            counts,
                            output_size);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::getROIVol(const std::vector<tetrahedron_global_id>& tets) const {
    double sum = 0.0;
    for (const auto& tet: tets) {
        sum += pTets[tet]->vol();
    }
    return sum;
}

double Tetexact::getROIVol(const std::string& ROI_id) const {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    return getROIVol(roi->second);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::getROIArea(const std::string& ROI_id) const {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    double sum = 0.0;
    for (auto const& tri: roi->second) {
        sum += pTris[tri]->area();
    }
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::getROITriSpecCount(const std::vector<triangle_global_id>& indices,
                                    const std::string& s) const {
    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;
    double sum = 0.0;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: indices) {
        auto tri = pTris.at(tidx);
        if (tri == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        solver::spec_local_id slidx = tri->patchdef()->specG2L(sgidx);
        if (slidx.unknown()) {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }

        sum += tri->pools()[slidx];
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log") << "The following triangles have not been assigned to a "
                                        "patch, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log") << "Species " << s
                                     << " has not been defined in the following triangles, fill in "
                                        "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
    return sum;
}

double Tetexact::getROITetSpecCount(const std::vector<tetrahedron_global_id>& indices,
                                    const std::string& s) const {
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;
    double sum = 0.0;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: indices) {
        auto tet = pTets.at(tidx);
        if (tet == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        solver::spec_local_id slidx = tet->compdef()->specG2L(sgidx);
        if (slidx.unknown()) {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }

        sum += tet->pools()[slidx];
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log") << "Species " << s
                                     << " has not been defined in the following tetrahedrons, fill "
                                        "in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
    return sum;
}


double Tetexact::getROISpecCount(const std::string& ROI_id, std::string const& s) const {
    {
        auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id, 0, false);
        if (roi != mesh().rois.end<tetmesh::ROI_TET>()) {
            return getROITetSpecCount(roi->second, s);
        }
    }
    {
        auto const& roi = mesh().rois.get<tetmesh::ROI_TRI>(ROI_id, 0, false);
        if (roi != mesh().rois.end<tetmesh::ROI_TRI>()) {
            return getROITriSpecCount(roi->second, s);
        }
    }
    std::ostringstream os;
    os << "Error: Cannot find suitable ROI for the function call getROICount.\n";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////


void Tetexact::setROITriSpecCount(const std::vector<triangle_global_id>& indices,
                                  std::string const& s,
                                  double count) {
    double total_weight = 0.0;
    const solver::spec_global_id sgidx = statedef().getSpecIdx(s);
    std::vector<Tri*> apply;

    for (auto const& tidx: indices) {
        if (tidx >= pTris.size()) {
            ArgErrLog("ROI refers to nonexistent triangle " + std::to_string(tidx));
        }

        Tri* tri = pTris[tidx];
        if (tri == nullptr || tri->patchdef()->specG2L(sgidx).unknown()) {
            continue;
        }

        apply.push_back(tri);
        total_weight += tri->area();
    }

    util::distribute_quantity<uint>(
        count,
        apply.begin(),
        apply.end(),
        [](const std::vector<Tri*>::iterator& tri) { return (*tri)->area(); },
        [sgidx](Tri* tri, uint c) { tri->setCount(tri->patchdef()->specG2L(sgidx), c); },
        [sgidx](Tri* tri) { tri->incCount(tri->patchdef()->specG2L(sgidx), 1); },
        *rng(),
        total_weight);

    for (auto& tri: apply) {
        _updateSpec(*tri);
    }
}

void Tetexact::setROITetSpecCount(const std::vector<tetrahedron_global_id>& indices,
                                  std::string const& s,
                                  double count) {
    double total_weight = 0.0;
    const solver::spec_global_id sgidx = statedef().getSpecIdx(s);
    std::vector<Tet*> apply;

    for (auto const& tidx: indices) {
        if (tidx >= pTets.size()) {
            ArgErrLog("ROI refers to nonexistent tetrahedron " + std::to_string(tidx));
        }

        Tet* tet = pTets[tidx];
        if (tet == nullptr || tet->compdef()->specG2L(sgidx).unknown()) {
            continue;
        }

        apply.push_back(tet);
        total_weight += tet->vol();
    }

    util::distribute_quantity<uint>(
        count,
        apply.begin(),
        apply.end(),
        [](const std::vector<Tet*>::iterator& tet) { return (*tet)->vol(); },
        [sgidx](Tet* tet, uint c) { tet->setCount(tet->compdef()->specG2L(sgidx), c); },
        [sgidx](Tet* tet) { tet->incCount(tet->compdef()->specG2L(sgidx), 1); },
        *rng(),
        total_weight);

    for (auto& tet: apply) {
        _updateSpec(*tet);
    }
}

void Tetexact::setROISpecCount(const std::string& ROI_id, std::string const& s, double count) {
    {
        auto const& roi = mesh().rois.get<tetmesh::ROI_TRI>(ROI_id, 0, false);
        if (roi != mesh().rois.end<tetmesh::ROI_TRI>()) {
            setROITriSpecCount(roi->second, s, count);
            return;
        }
    }
    {
        auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id, 0, false);
        if (roi != mesh().rois.end<tetmesh::ROI_TET>()) {
            setROITetSpecCount(roi->second, s, count);
            return;
        }
    }
    ArgErrLog("can only set counts in tetrahedra or triangle ROIs");
}

////////////////////////////////////////////////////////////////////////////////


double Tetexact::getROISpecAmount(const std::string& ROI_id, std::string const& s) const {
    double count = getROISpecCount(ROI_id, s);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////


void Tetexact::setROISpecAmount(const std::string& ROI_id, std::string const& s, double amount) {
    setROISpecCount(ROI_id, s, amount * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

double Tetexact::getROISpecConc(const std::string& ROI_id, std::string const& s) const {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    const double count = getROITetSpecCount(roi->second, s);
    const double vol = getROIVol(roi->second);
    return count / (1.0e3 * vol * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::setROISpecConc(const std::string& ROI_id, std::string const& s, double conc) {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id, 0, false);
    if (roi == mesh().rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("can only set concentrations in tetrahedra ROIs");
    }

    double total_weight = 0.0;
    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    std::vector<Tet*> apply;
    for (auto const& tidx: roi->second) {
        if (tidx >= pTets.size()) {
            ArgErrLog("ROI refers to nonexistent tetrahedron " + std::to_string(tidx));
        }

        Tet* tet = pTets[tidx];
        if (tet == nullptr || tet->compdef()->specG2L(sgidx).unknown()) {
            continue;
        }

        apply.push_back(tet);
        total_weight += tet->vol();
    }

    double count = conc * (1.0e3 * total_weight * math::AVOGADRO);

    util::distribute_quantity<uint>(
        count,
        apply.begin(),
        apply.end(),
        [](const std::vector<Tet*>::iterator& tet) { return (*tet)->vol(); },
        [sgidx](Tet* tet, uint c) { tet->setCount(tet->compdef()->specG2L(sgidx), c); },
        [sgidx](Tet* tet) { tet->incCount(tet->compdef()->specG2L(sgidx), 1); },
        *rng(),
        total_weight);

    for (auto& tet: apply) {
        _updateSpec(*tet);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::setROITriSpecClamped(const std::vector<triangle_global_id>& indices,
                                    std::string const& s,
                                    bool b) {
    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: indices) {
        auto tri = pTris.at(tidx);
        if (tri == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        solver::spec_local_id slidx = tri->patchdef()->specG2L(sgidx);
        if (slidx.unknown()) {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }

        tri->setClamped(slidx, b);
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log") << "The following triangles have not been assigned to a "
                                        "patch, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log") << "Species " << s
                                     << " has not been defined in the following triangles, fill in "
                                        "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
}

void Tetexact::setROITetSpecClamped(const std::vector<tetrahedron_global_id>& indices,
                                    std::string const& s,
                                    bool b) {
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: indices) {
        auto tet = pTets.at(tidx);
        if (tet == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        solver::spec_local_id slidx = tet->compdef()->specG2L(sgidx);
        if (slidx.unknown()) {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }

        tet->setClamped(slidx, b);
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log") << "Species " << s
                                     << " has not been defined in the following tetrahedrons, fill "
                                        "in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
}

void Tetexact::setROISpecClamped(const std::string& ROI_id, std::string const& s, bool b) {
    {
        auto const& roi = mesh().rois.get<tetmesh::ROI_TRI>(ROI_id, 0, false);
        if (roi != mesh().rois.end<tetmesh::ROI_TRI>()) {
            setROITriSpecClamped(roi->second, s, b);
            return;
        }
    }
    {
        auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id, 0, false);
        if (roi != mesh().rois.end<tetmesh::ROI_TET>()) {
            setROITetSpecClamped(roi->second, s, b);
            return;
        }
    }
    std::ostringstream os;
    os << "Error: Cannot find suitable ROI for the function call setROICount.\n";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::setROIReacK(const std::string& ROI_id, std::string const& r, double kf) {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_reac_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream reac_undefined;

    const solver::reac_global_id rgidx = statedef().getReacIdx(r);

    for (auto const& tidx: roi->second) {
        auto tet = pTets.at(tidx);
        if (tet == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        solver::reac_local_id rlidx = tet->compdef()->reacG2L(rgidx);
        if (rlidx.unknown()) {
            reac_undefined << tidx << " ";
            has_reac_warning = true;
            continue;
        }

        tet->reac(rlidx).setKcst(kf);
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_reac_warning) {
        CLOG(WARNING, "general_log") << "Reac " << r
                                     << " has not been defined in the following tetrahedrons, no "
                                        "change is applied to them:\n";
        CLOG(WARNING, "general_log") << reac_undefined.str() << "\n";
    }

    _update();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::setROISReacK(const std::string& ROI_id, std::string const& sr, double kf) {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tri_warning = false;
    bool has_sreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream sreac_undefined;

    solver::sreac_global_id srgidx = statedef().getSReacIdx(sr);

    for (auto const& tidx: roi->second) {
        auto tri = pTris.at(tidx);
        if (tri == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        solver::sreac_local_id srlidx = tri->patchdef()->sreacG2L(srgidx);
        if (srlidx.unknown()) {
            sreac_undefined << tidx << " ";
            has_sreac_warning = true;
            continue;
        }

        tri->sreac(srlidx).setKcst(kf);
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log") << "The following triangles have not been assigned to a "
                                        "patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_sreac_warning) {
        CLOG(WARNING, "general_log")
            << "SReac " << sr
            << " has not been defined in the following patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << sreac_undefined.str() << "\n";
    }

    _update();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::setROIDiffD(const std::string& ROI_id, std::string const& d, double dk) {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_diff_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream diff_undefined;

    solver::diff_global_id dgidx = statedef().getDiffIdx(d);

    for (auto const& tidx: roi->second) {
        auto tet = pTets.at(tidx);
        if (tet == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        solver::diff_local_id dlidx = tet->compdef()->diffG2L(dgidx);
        if (dlidx.unknown()) {
            diff_undefined << tidx << " ";
            has_diff_warning = true;
            continue;
        }

        tet->diff(dlidx).setDcst(dk);
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_diff_warning) {
        CLOG(WARNING, "general_log") << "Diff " << d
                                     << " has not been defined in the following tetrahedrons, no "
                                        "change is applied to them:\n";
        CLOG(WARNING, "general_log") << diff_undefined.str() << "\n";
    }

    _update();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::setROIReacActive(const std::string& ROI_id, std::string const& r, bool a) {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_reac_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream reac_undefined;

    solver::reac_global_id rgidx = statedef().getReacIdx(r);

    for (auto const& tidx: roi->second) {
        auto tet = pTets.at(tidx);
        if (tet == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        solver::reac_local_id rlidx = tet->compdef()->reacG2L(rgidx);
        if (rlidx.unknown()) {
            reac_undefined << tidx << " ";
            has_reac_warning = true;
            continue;
        }

        tet->reac(rlidx).setActive(a);
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_reac_warning) {
        CLOG(WARNING, "general_log") << "Reac " << r
                                     << " has not been defined in the following tetrahedrons, no "
                                        "change is applied to them:\n";
        CLOG(WARNING, "general_log") << reac_undefined.str() << "\n";
    }

    _update();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::setROISReacActive(const std::string& ROI_id, std::string const& sr, bool a) {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tri_warning = false;
    bool has_sreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream sreac_undefined;

    solver::sreac_global_id srgidx = statedef().getSReacIdx(sr);

    for (auto const& tidx: roi->second) {
        auto tri = pTris.at(tidx);
        if (tri == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        solver::sreac_local_id srlidx = tri->patchdef()->sreacG2L(srgidx);
        if (srlidx.unknown()) {
            sreac_undefined << tidx << " ";
            has_sreac_warning = true;
            continue;
        }

        tri->sreac(srlidx).setActive(a);
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log") << "The following triangles have not been assigned to a "
                                        "patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_sreac_warning) {
        CLOG(WARNING, "general_log")
            << "SReac " << sr
            << " has not been defined in the following patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << sreac_undefined.str() << "\n";
    }

    _update();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::setROIDiffActive(const std::string& ROI_id, std::string const& d, bool a) {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_diff_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream diff_undefined;

    solver::diff_global_id dgidx = statedef().getDiffIdx(d);

    for (auto const& tidx: roi->second) {
        auto tet = pTets.at(tidx);
        if (tet == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        solver::diff_local_id dlidx = tet->compdef()->diffG2L(dgidx);
        if (dlidx.unknown()) {
            diff_undefined << tidx << " ";
            has_diff_warning = true;
            continue;
        }

        tet->diff(dlidx).setActive(a);
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_diff_warning) {
        CLOG(WARNING, "general_log") << "Diff " << d
                                     << " has not been defined in the following tetrahedrons, no "
                                        "change is applied to them:\n";
        CLOG(WARNING, "general_log") << diff_undefined.str() << "\n";
    }

    _update();
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::setROIVDepSReacActive(const std::string& ROI_id, std::string const& vsr, bool a) {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tri_warning = false;
    bool has_vsreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream vsreac_undefined;

    solver::vdepsreac_global_id vsrgidx = statedef().getVDepSReacIdx(vsr);

    for (auto const& tidx: roi->second) {
        auto tri = pTris.at(tidx);
        if (tri == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        solver::vdepsreac_local_id vsrlidx = tri->patchdef()->vdepsreacG2L(vsrgidx);
        if (vsrlidx.unknown()) {
            vsreac_undefined << tidx << " ";
            has_vsreac_warning = true;
            continue;
        }

        tri->vdepsreac(vsrlidx).setActive(a);
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log") << "The following triangles have not been assigned to a "
                                        "patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_vsreac_warning) {
        CLOG(WARNING, "general_log")
            << "VDepSReac " << vsr
            << " has not been defined in the following patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << vsreac_undefined.str() << "\n";
    }
    _update();
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long Tetexact::getROIReacExtent(const std::string& ROI_id,
                                              std::string const& r) const {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_reac_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream reac_undefined;

    solver::reac_global_id rgidx = statedef().getReacIdx(r);

    unsigned long long sum = 0;

    for (auto const& tidx: roi->second) {
        auto tet = pTets.at(tidx);
        if (tet == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        solver::reac_local_id rlidx = tet->compdef()->reacG2L(rgidx);
        if (rlidx.unknown()) {
            reac_undefined << tidx << " ";
            has_reac_warning = true;
            continue;
        }

        sum += tet->reac(rlidx).getExtent();
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_reac_warning) {
        CLOG(WARNING, "general_log") << "Reac " << r
                                     << " has not been defined in the following tetrahedrons, no "
                                        "change is applied to them:\n";
        CLOG(WARNING, "general_log") << reac_undefined.str() << "\n";
    }

    return sum;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::resetROIReacExtent(const std::string& ROI_id, std::string const& r) {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_reac_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream reac_undefined;

    solver::reac_global_id rgidx = statedef().getReacIdx(r);

    for (auto const& tidx: roi->second) {
        auto tet = pTets.at(tidx);
        if (tet == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        solver::reac_local_id rlidx = tet->compdef()->reacG2L(rgidx);
        if (rlidx.unknown()) {
            reac_undefined << tidx << " ";
            has_reac_warning = true;
            continue;
        }

        tet->reac(rlidx).resetExtent();
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_reac_warning) {
        CLOG(WARNING, "general_log") << "Reac " << r
                                     << " has not been defined in the following tetrahedrons, no "
                                        "change is applied to them:\n";
        CLOG(WARNING, "general_log") << reac_undefined.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long Tetexact::getROISReacExtent(const std::string& ROI_id,
                                               std::string const& sr) const {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tri_warning = false;
    bool has_sreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream sreac_undefined;

    solver::sreac_global_id srgidx = statedef().getSReacIdx(sr);

    unsigned long long sum = 0;

    for (auto const& tidx: roi->second) {
        auto tri = pTris.at(tidx);
        if (tri == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        solver::sreac_local_id srlidx = tri->patchdef()->sreacG2L(srgidx);
        if (srlidx.unknown()) {
            sreac_undefined << tidx << " ";
            has_sreac_warning = true;
            continue;
        }

        sum += tri->sreac(srlidx).getExtent();
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log") << "The following triangles have not been assigned to a "
                                        "patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_sreac_warning) {
        CLOG(WARNING, "general_log")
            << "SReac " << sr
            << " has not been defined in the following patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << sreac_undefined.str() << "\n";
    }

    return sum;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::resetROISReacExtent(const std::string& ROI_id, std::string const& sr) {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tri_warning = false;
    bool has_sreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream sreac_undefined;

    solver::sreac_global_id srgidx = statedef().getSReacIdx(sr);

    for (auto const& tidx: roi->second) {
        auto tri = pTris.at(tidx);
        if (tri == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        solver::sreac_local_id srlidx = tri->patchdef()->sreacG2L(srgidx);
        if (srlidx.unknown()) {
            sreac_undefined << tidx << " ";
            has_sreac_warning = true;
            continue;
        }

        tri->sreac(srlidx).resetExtent();
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log") << "The following triangles have not been assigned to a "
                                        "patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_sreac_warning) {
        CLOG(WARNING, "general_log")
            << "SReac " << sr
            << " has not been defined in the following patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << sreac_undefined.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long Tetexact::getROIDiffExtent(const std::string& ROI_id,
                                              std::string const& d) const {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_diff_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream diff_undefined;

    solver::diff_global_id dgidx = statedef().getDiffIdx(d);

    unsigned long long sum = 0;

    for (auto const& tidx: roi->second) {
        auto tet = pTets.at(tidx);
        if (tet == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        solver::diff_local_id dlidx = tet->compdef()->diffG2L(dgidx);
        if (dlidx.unknown()) {
            diff_undefined << tidx << " ";
            has_diff_warning = true;
            continue;
        }

        sum += tet->diff(dlidx).getExtent();
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_diff_warning) {
        CLOG(WARNING, "general_log") << "Diff " << d
                                     << " has not been defined in the following tetrahedrons, no "
                                        "change is applied to them:\n";
        CLOG(WARNING, "general_log") << diff_undefined.str() << "\n";
    }

    return sum;
}

////////////////////////////////////////////////////////////////////////////////

void Tetexact::resetROIDiffExtent(const std::string& ROI_id, std::string const& d) {
    auto const& roi = mesh().rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh().rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_diff_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream diff_undefined;

    solver::diff_global_id dgidx = statedef().getDiffIdx(d);

    for (auto const& tidx: roi->second) {
        auto tet = pTets.at(tidx);
        if (tet == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        solver::diff_local_id dlidx = tet->compdef()->diffG2L(dgidx);
        if (dlidx.unknown()) {
            diff_undefined << tidx << " ";
            has_diff_warning = true;
            continue;
        }

        tet->diff(dlidx).resetExtent();
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_diff_warning) {
        CLOG(WARNING, "general_log") << "Diff " << d
                                     << " has not been defined in the following tetrahedrons, no "
                                        "change is applied to them:\n";
        CLOG(WARNING, "general_log") << diff_undefined.str() << "\n";
    }
}

}  // namespace steps::tetexact
