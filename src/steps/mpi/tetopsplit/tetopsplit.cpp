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

#include "tetopsplit.hpp"

// Standard library headers.
#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <iomanip>
#include <iterator>
#include <numeric>

#include <mpi.h>

#include "diff.hpp"
#include "ghkcurr.hpp"
#include "reac.hpp"
#include "sdiff.hpp"
#include "sreac.hpp"
#include "vdepsreac.hpp"

#include "math/constants.hpp"
#include "math/point.hpp"
#include "mpi/mpi_common.hpp"
#include "solver/chandef.hpp"
#include "solver/compdef.hpp"
#include "solver/diffboundarydef.hpp"
#include "solver/ohmiccurrdef.hpp"
#include "solver/patchdef.hpp"
#include "solver/sdiffboundarydef.hpp"
#include "solver/types.hpp"
#include "solver/vdepsreacdef.hpp"

#include "solver/efield/dVsolver.hpp"
#include "solver/efield/efield.hpp"
#ifdef USE_PETSC
#include "solver/efield/dVsolver_petsc.hpp"
#endif

#include "util/checkpointing.hpp"
#include "util/collections.hpp"
#include "util/distribute.hpp"
#include "util/error.hpp"
#include "util/profile/profiler_interface.hpp"

namespace steps::mpi::tetopsplit {

////////////////////////////////////////////////////////////////////////////////

void schedIDXSet_To_Vec(SchedIDXSet const& s, SchedIDXVec& v) {
    v.resize(s.size());
    std::copy(s.begin(), s.end(), v.begin());
}

////////////////////////////////////////////////////////////////////////////////

TetOpSplitP::TetOpSplitP(model::Model* m,
                         wm::Geom* g,
                         const rng::RNGptr& r,
                         int calcMembPot,
                         std::vector<int> const& tet_hosts,
                         const std::map<triangle_global_id, int>& tri_hosts,
                         std::vector<int> const& wm_hosts)
    : API(m, g, r)
    , pEFoption(static_cast<EF_solver>(calcMembPot))
    , gen(rd()) {
    if (rng() == nullptr) {
        std::ostringstream os;
        os << "No RNG provided to solver initializer function";
        ArgErrLog(os.str());
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nHosts);

    // Perform upcast.
    pMesh = dynamic_cast<tetmesh::Tetmesh*>(geom());
    if (pMesh == nullptr) {
        ArgErrLog(
            "Geometry description to solver::Tetexact solver "
            "constructor is not a valid tetmesh::Tetmesh object.");
    }

    if (tet_hosts.empty()) {
        // partition the mesh using metis
        _partition();
    } else {
        for (uint t = 0; t < tet_hosts.size(); t++) {
            tetHosts.emplace(tetrahedron_global_id(t), tet_hosts[t]);
        }
        triHosts = tri_hosts;
        wmHosts = wm_hosts;
    }

    // All initialization code now in _setup() to allow EField solver to be
    // derived and create EField local objects within the constructor
    _setup();
    _updateLocal();
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

TetOpSplitP::~TetOpSplitP() {
    for (auto& c: pComps) {
        delete c;
    }
    for (auto& p: pPatches) {
        delete p;
    }
    for (auto& db: pDiffBoundaries) {
        delete db;
    }
    for (auto& sdb: pSDiffBoundaries) {
        delete sdb;
    }
    for (auto& wvol: pWmVols) {
        delete wvol;
    }
    for (auto& t: pTets) {
        delete t;
    }
    for (auto& t: pTris) {
        delete t;
    }
    for (auto& g: nGroups) {
        g->free_indices();
        delete g;
    }
    for (auto& g: pGroups) {
        g->free_indices();
        delete g;
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_partition() {
    uint n_tets_tris[2];
    std::vector<tetrahedron_global_id> tet_ids;
    std::vector<int> tet_hosts;

    std::vector<triangle_global_id> tri_ids;
    std::vector<int> tri_hosts;

    if (myRank == 0) {
        pMesh->_autoPartition(1, 0, nHosts, tetHosts, triHosts);
        n_tets_tris[0] = tetHosts.size();
        tet_ids.reserve(tetHosts.size());
        tet_hosts.reserve(tetHosts.size());
        for (auto [tet, host]: tetHosts) {
            tet_ids.push_back(tet);
            tet_hosts.push_back(host);
        }
        n_tets_tris[1] = triHosts.size();
        tri_ids.reserve(triHosts.size());
        tri_hosts.reserve(triHosts.size());
        for (auto [tri, host]: triHosts) {
            tri_ids.push_back(tri);
            tri_hosts.push_back(host);
        }
    }

    MPI_Bcast(n_tets_tris, 2, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    if (myRank != 0) {
        tet_ids.resize(n_tets_tris[0]);
        tet_hosts.resize(n_tets_tris[0]);
        tri_ids.resize(n_tets_tris[1]);
        tri_hosts.resize(n_tets_tris[1]);
    }

    MPI_Bcast(tet_ids.data(), n_tets_tris[0], MPI_STEPS_INDEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(tet_hosts.data(), n_tets_tris[0], MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(tri_ids.data(), n_tets_tris[1], MPI_STEPS_INDEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(tri_hosts.data(), n_tets_tris[1], MPI_INT, 0, MPI_COMM_WORLD);

    if (myRank != 0) {
        for (uint t = 0; t < n_tets_tris[0]; t++) {
            tetHosts[tet_ids[t]] = tet_hosts[t];
        }
        for (uint t = 0; t < n_tets_tris[1]; t++) {
            triHosts[tri_ids[t]] = tri_hosts[t];
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::checkpoint(std::string const& file_name) {
    std::fstream cp_file;

    cp_file.open(file_name.c_str(), std::fstream::out | std::fstream::binary | std::fstream::trunc);
    util::checkpoint(cp_file, myRank);

    API::checkpoint(cp_file);

    statedef().checkpoint(cp_file);

    for (auto const& c: pComps) {
        c->checkpoint(cp_file);
    }

    for (auto const& p: pPatches) {
        p->checkpoint(cp_file);
    }

    for (auto const& db: pDiffBoundaries) {
        db->checkpoint(cp_file);
    }

    for (auto const& sdb: pSDiffBoundaries) {
        sdb->checkpoint(cp_file);
    }

    for (auto const& t: pTets) {
        if (t != nullptr) {
            t->checkpoint(cp_file);
        }
    }

    for (auto const& t: pTris) {
        if (t != nullptr) {
            t->checkpoint(cp_file);
        }
    }

    util::checkpoint(cp_file, pSum);
    util::checkpoint(cp_file, nSum);
    util::checkpoint(cp_file, pA0);

    util::checkpoint(cp_file, diffSep);
    util::checkpoint(cp_file, sdiffSep);

    for (auto const& k: pKProcs) {
        if (k != nullptr) {
            k->checkpoint(cp_file);
        }
    }

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

    util::checkpoint(cp_file, recomputeUpdPeriod);
    util::checkpoint(cp_file, reacExtent);
    util::checkpoint(cp_file, diffExtent);
    util::checkpoint(cp_file, nIteration);
    util::checkpoint(cp_file, updPeriod);
    util::checkpoint(cp_file, diffApplyThreshold);

    if (efflag()) {
        util::checkpoint(cp_file, pTemp);
        util::checkpoint(cp_file, pEFDT);
        pEField->checkpoint(cp_file);
    }
}

///////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::restore(std::string const& file_name) {
    std::fstream cp_file;

    cp_file.open(file_name.c_str(), std::fstream::in | std::fstream::binary);

    cp_file.seekg(0);
    util::compare(cp_file, myRank);

    API::restore(cp_file);

    statedef().restore(cp_file);

    for (auto const& c: pComps) {
        c->restore(cp_file);
    }

    for (auto const& p: pPatches) {
        p->restore(cp_file);
    }

    for (auto const& db: pDiffBoundaries) {
        db->restore(cp_file);
    }

    for (auto const& sdb: pSDiffBoundaries) {
        sdb->restore(cp_file);
    }

    for (auto const& t: pTets) {
        if (t != nullptr) {
            t->restore(cp_file);
        }
    }

    for (auto const& t: pTris) {
        if (t != nullptr) {
            t->restore(cp_file);
        }
    }

    util::restore(cp_file, pSum);
    util::restore(cp_file, nSum);
    util::restore(cp_file, pA0);

    util::restore(cp_file, diffSep);
    util::restore(cp_file, sdiffSep);

    for (auto const& k: pKProcs) {
        if (k != nullptr) {
            k->restore(cp_file);
        }
    }

    std::vector<Diff*> diffs_temp(pDiffs);
    for (auto const& d: diffs_temp) {
        auto pos = d->crData.pos;
        AssertLog(pos < pDiffs.size());
        pDiffs[pos] = d;
    }

    std::vector<SDiff*> sdiffs_temp(pSDiffs);
    for (auto const& sd: sdiffs_temp) {
        auto pos = sd->crData.pos;
        AssertLog(pos < pSDiffs.size());
        pSDiffs[pos] = sd;
    }

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

    util::restore(cp_file, recomputeUpdPeriod);
    util::restore(cp_file, reacExtent);
    util::restore(cp_file, diffExtent);
    util::restore(cp_file, nIteration);
    util::restore(cp_file, updPeriod);
    util::restore(cp_file, diffApplyThreshold);

    if (efflag()) {
        util::restore(cp_file, pTemp);
        util::restore(cp_file, pEFDT);
        pEField->restore(cp_file);
        _refreshEFTrisV();
    }

    cp_file.close();
}

////////////////////////////////////////////////////////////////////////////////

std::string TetOpSplitP::getSolverName() const {
    return "Parallel TetOpSplit";
}

////////////////////////////////////////////////////////////////////////////////

std::string TetOpSplitP::getSolverDesc() const {
    return "Parallel approximate stochastic method in tetrahedral mesh";
}

////////////////////////////////////////////////////////////////////////////////

std::string TetOpSplitP::getSolverAuthors() const {
    return "Iain Hepburn, Weiliang Chen, Stefan Wils, Sam Yates";
}

////////////////////////////////////////////////////////////////////////////////

std::string TetOpSplitP::getSolverEmail() const {
    return "steps.dev@gmail.com";
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setup() {
    // First initialise the pTets, pTris vector, because
    // want tets and tris to maintain indexing from Geometry
    uint ntets = mesh()->countTets();
    uint ntris = mesh()->countTris();
    uint ncomps = mesh()->_countComps();

    pTets.container().assign(ntets, nullptr);
    pTris.container().assign(ntris, nullptr);
    pWmVols.container().assign(ncomps, nullptr);
    diffSep = 0;
    sdiffSep = 0;
    // Now create the actual compartments.
    for (auto const& c: statedef().comps()) {
        solver::comp_global_id compdef_gidx = c->gidx();
        uint comp_idx = _addComp(c);
        AssertLog(compdef_gidx == comp_idx);
    }
    // Create the actual patches.
    for (auto const& p: statedef().patches()) {
        solver::patch_global_id patchdef_gidx = p->gidx();
        uint patch_idx = _addPatch(p);
        AssertLog(patchdef_gidx == patch_idx);
    }

    // Create the diffusion boundaries
    for (auto const& db: statedef().diffBoundaries()) {
        solver::diffboundary_global_id diffboundary_gidx = db->gidx();
        uint diffb_idx = _addDiffBoundary(db);
        AssertLog(diffboundary_gidx == diffb_idx);
    }

    // Create the surface diffusion boundaries
    for (auto const& sdb: statedef().sdiffBoundaries()) {
        solver::sdiffboundary_global_id sdiffboundary_gidx = sdb->gidx();
        uint sdiffb_idx = _addSDiffBoundary(sdb);
        AssertLog(sdiffboundary_gidx == sdiffb_idx);
    }

    auto npatches = pPatches.size();
    AssertLog(mesh()->_countPatches() == npatches);
    int rk;
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);

    for (auto p: solver::patch_global_id::range(npatches)) {
        // Add the tris for this patch
        // We have checked the indexing - p is the global index
        wm::Patch* wmpatch = mesh()->_getPatch(p);

        // Perform upcast
        auto* tmpatch = dynamic_cast<tetmesh::TmPatch*>(wmpatch);
        if (tmpatch == nullptr) {
            ArgErrLog(
                "Well-mixed patches not supported in "
                "solver::TetOpSplitP solver.");
        }

        Patch* localpatch = pPatches[p];

        /*** Previous impl
        // Create a map between edges and adjacent tris in this patch
        std::map<uint, std::vector<uint> > bar2tri;
        for (uint tri: tmpatch->_getAllTriIndices())
        {
            const uint *bars = pMesh->_getTriBars(tri);
            for (int i = 0; i < 3; ++i)
                bar2tri[bars[i]].push_back(tri);
        }

        for (uint tri: tmpatch->_getAllTriIndices())
        ***/

        // Create a map between edges and adjacent tris in this patch
        std::map<bar_id_t, std::vector<triangle_global_id>> bar2tri;

        // We need to go through all patches to record bar2tri mapping
        // for all connected triangle neighbors even they are in different
        // patches, because their information is needed for surface diffusion
        // boundary

        for (auto bar_p: solver::patch_global_id::range(npatches)) {
            auto patch = pMesh->_getPatch(bar_p);
            AssertLog(patch != nullptr);
            auto* bar_patch = dynamic_cast<tetmesh::TmPatch*>(patch);
            if (bar_patch == nullptr) {
                ProgErrLog("Unable to cast the Patch object to TmPatch.");
            } else {
                for (const auto& tri: bar_patch->_getAllTriIndices()) {
                    for (const auto& bar: pMesh->_getTriBars(tri)) {
                        bar2tri[bar].push_back(tri);
                    }
                }
            }
        }

        auto const& tri_idxs = tmpatch->_getAllTriIndices();

        for (auto i = 0u; i < tri_idxs.size(); ++i)
        //*** end of new impl
        {
            auto tri = tri_idxs[i];
            AssertLog(pMesh->getTriPatch(tri) == tmpatch);

            double area = pMesh->getTriArea(tri);

            // NB: Tri vertices may not be in consistent order, so use bar interface.
            const auto& tri_bars = pMesh->_getTriBars(tri);
            double l[3] = {0, 0, 0};

            for (auto j = 0u; j < tri_bars.size(); ++j) {
                const auto* v = pMesh->_getBar(tri_bars[j]);
                l[j] = distance(pMesh->_getVertex(v[0]), pMesh->_getVertex(v[1]));
            }

            // Get neighboring tris
            std::array<triangle_global_id, 3> tris{{std::nullopt, std::nullopt, std::nullopt}};
            for (auto j = 0u; j < tri_bars.size(); ++j) {
                const std::vector<triangle_global_id>& neighb_tris = bar2tri[tri_bars[j]];
                for (const auto& neighb_tri: neighb_tris) {
                    if (neighb_tri == tri || pMesh->getTriPatch(neighb_tri) == nullptr) {
                        continue;
                    }
                    tris[j] = neighb_tri;
                    break;
                }
            }

            const math::point3d& baryc = pMesh->_getTriBarycenter(tri);
            math::point3d d;
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
    AssertLog(mesh()->_countComps() == ncomps);

    for (auto c: solver::comp_global_id::range(ncomps)) {
        // Add the tets for this comp.
        // We have checked the indexing - c is the global index.
        auto wmcomp = mesh()->_getComp(c);
        // Perform upcast
        auto tmcomp = dynamic_cast<tetmesh::TmComp*>(wmcomp);
        if (tmcomp != nullptr) {
            Comp* localcomp = pComps[c];

            for (const auto& tet: tmcomp->_getAllTetIndices()) {
                AssertLog(pMesh->getTetComp(tet) == tmcomp);

                double vol = pMesh->getTetVol(tet);

                const auto& tris = pMesh->_getTetTriNeighb(tet);

                double a[4] = {0, 0, 0, 0};
                for (uint j = 0; j < 4; ++j) {
                    a[j] = pMesh->getTriArea(tris[j]);
                }

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
            Comp* localcomp = pComps[c];
            _addWmVol(c, localcomp, localcomp->def()->vol());
            AssertLog(pWmVols[c] != nullptr);

            // Now find all the triangles that reference this well-mixed volume
            // and set the inner or outer tetrahedron index accordingly.

            uint nopatches = wmcomp->_countOPatches();
            for (uint i = 0; i < nopatches; ++i) {
                wm::Patch* op = wmcomp->_getOPatch(i);
                //     Comp may have no outer patch
                if (op == nullptr) {
                    continue;
                }
                auto* comp_opatch = dynamic_cast<tetmesh::TmPatch*>(op);
                if (comp_opatch == nullptr) {
                    ProgErrLog("Compartment outer patch is not a TmPatch.");
                } else {
                    for (auto tri: comp_opatch->_getAllTriIndices()) {
                        pTris[tri]->setInnerTet(pWmVols[c]);
                        // Add triangle to WmVols' table of neighbouring triangles.
                        pWmVols[c]->setNextTri(pTris[tri]);
                    }
                }
            }

            uint nipatches = wmcomp->_countIPatches();
            for (uint i = 0; i < nipatches; ++i) {
                wm::Patch* ip = wmcomp->_getIPatch(i);
                // Comp may not have an inner patch
                if (ip == nullptr) {
                    continue;
                }
                auto* comp_ipatch = dynamic_cast<tetmesh::TmPatch*>(ip);
                if (comp_ipatch == nullptr) {
                    ProgErrLog("Compartment inner patch is not a TmPatch.");
                } else {
                    for (auto tri: comp_ipatch->_getAllTriIndices()) {
                        pTris[tri]->setOuterTet(pWmVols[c]);
                        // Add triangle to WmVols' table of neighbouring triangles.
                        pWmVols[c]->setNextTri(pTris[tri]);
                    }
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
            auto tet = pTets[t]->tet(j);
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
            auto tri = pTris[t]->tri(j);
            if (tri.valid() && pTris[tri] != nullptr) {
                pTris[t]->setNextTri(j, pTris[tri]);
            }
        }

        // By convention, triangles in a patch should have an inner tetrahedron
        // defined (neighbouring tets 'flipped' if necessary in Tetmesh) but not
        // necessarily an outer tet 17/3/10- actually this is not the case any more
        // with well-mixed compartments
        //
        auto tetinner = pTris[t]->tet(0);
        auto tetouter = pTris[t]->tet(1);

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

                    // Now with diffusion boundaries, meaning tets can have neighbours
                    // that are in different comps, we must check the compartment
                    Tet* tet_in = pTets[tetinner];
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
                    Tet* tet_out = pTets[tetouter];

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
    uint ndiffbnds = pDiffBoundaries.size();
    AssertLog(ndiffbnds == mesh()->_countDiffBoundaries());

    for (uint db = 0; db < ndiffbnds; ++db) {
        DiffBoundary* localdiffb = pDiffBoundaries[db];

        solver::comp_global_id compAidx = localdiffb->def()->compa();
        solver::comp_global_id compBidx = localdiffb->def()->compb();
        solver::Compdef* compAdef = statedef().compdef(compAidx);
        solver::Compdef* compBdef = statedef().compdef(compBidx);

        for (auto dbtri: localdiffb->def()->tris()) {
            const auto tri_tets = pMesh->_getTriTetNeighb(dbtri);

            auto tetAidx = tri_tets[0];
            auto tetBidx = tri_tets[1];
            AssertLog(tetAidx.valid() && tetBidx.valid());

            Tet* tetA = _tet(tetAidx);
            Tet* tetB = _tet(tetBidx);
            AssertLog(tetA != nullptr && tetB != nullptr);

            solver::Compdef* tetA_cdef = tetA->compdef();
            solver::Compdef* tetB_cdef = tetB->compdef();
            AssertLog(tetA_cdef != nullptr);
            AssertLog(tetB_cdef != nullptr);

            if (tetA_cdef != compAdef) {
                AssertLog(tetB_cdef == compAdef);
                AssertLog(tetA_cdef == compBdef);
            } else {
                AssertLog(tetB_cdef == compBdef);
                AssertLog(tetA_cdef == compAdef);
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
            localdiffb->setTetDirection(tetAidx, direction_idx_a);
            localdiffb->setTetDirection(tetBidx, direction_idx_b);
        }
        localdiffb->setComps(_comp(compAidx), _comp(compBidx));

        // Before the kprocs are set up ( in _setup) the tetrahedrons need to know
        // the diffusion boundary direction, so let's do it here  - the diff bounday
        // has had all tetrahedrons added

        // Might as well copy the vectors because we need to index through
        auto const& tets = localdiffb->getTets();
        auto const& tets_direction = localdiffb->getTetDirection();

        ntets = tets.size();
        AssertLog(ntets <= pTets.size());
        AssertLog(tets_direction.size() == ntets);

        for (uint t = 0; t < ntets; ++t) {
            _tet(tets[t])->setDiffBndDirection(tets_direction[t]);
        }
    }

    // Now loop over the surface diffusion boundaries:
    // 1) get all the bars and get the two triangles
    // 2) figure out which direction is the direction for a triangle
    // 3) add the triangle and the direction to local object

    // This is here because we need all tris to have been assigned correctly
    // to patches. Check every one and set the patchA and patchB for the db
    auto nsdiffbnds = pSDiffBoundaries.size();
    AssertLog(nsdiffbnds == mesh()->_countSDiffBoundaries());

    for (auto sdb = 0u; sdb < nsdiffbnds; ++sdb) {
        SDiffBoundary* localsdiffb = pSDiffBoundaries[sdb];

        solver::patch_global_id patchAidx = localsdiffb->def()->patcha();
        solver::patch_global_id patchBidx = localsdiffb->def()->patchb();
        solver::Patchdef* patchAdef = statedef().patchdef(patchAidx);
        solver::Patchdef* patchBdef = statedef().patchdef(patchBidx);

        for (auto sdbbar: localsdiffb->def()->bars()) {
            const auto* bar_tris = pMesh->_getBarTriNeighb(bar_id_t(sdbbar));

            auto triAidx = bar_tris[0];
            auto triBidx = bar_tris[1];
            AssertLog(triAidx.valid() && triBidx.valid());

            Tri* triA = _tri(triAidx);
            Tri* triB = _tri(triBidx);
            AssertLog(triA != nullptr && triB != nullptr);

            solver::Patchdef* triA_pdef = triA->patchdef();
            solver::Patchdef* triB_pdef = triB->patchdef();
            AssertLog(triA_pdef != nullptr);
            AssertLog(triB_pdef != nullptr);

            if (triA_pdef != patchAdef) {
                AssertLog(triB_pdef == patchAdef);
                AssertLog(triA_pdef == patchBdef);
            } else {
                AssertLog(triB_pdef == patchBdef);
                AssertLog(triA_pdef == patchAdef);
            }

            // Ok, checks over, lets get down to business
            int direction_idx_a = -1;
            int direction_idx_b = -1;

            const auto& triA_bars = pMesh->_getTriBars(triAidx);
            const auto& triB_bars = pMesh->_getTriBars(triBidx);

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
            localsdiffb->setTriDirection(triAidx, direction_idx_a);
            localsdiffb->setTriDirection(triBidx, direction_idx_b);
        }
        localsdiffb->setPatches(_patch(patchAidx), _patch(patchBidx));

        // Might as well copy the vectors because we need to index through
        auto const& tris = localsdiffb->getTris();
        auto const& tris_direction = localsdiffb->getTriDirection();

        ntris = tris.size();
        AssertLog(ntris <= pTris.size());
        AssertLog(tris_direction.size() == ntris);

        for (auto t = 0u; t < ntris; ++t) {
            _tri(tris[t])->setSDiffBndDirection(tris_direction[t]);
        }
    }

    for (auto& t: pTets) {
        if (t != nullptr) {
            t->setupKProcs(this);
        }
    }

    for (auto& wmv: pWmVols) {
        if (wmv != nullptr) {
            wmv->setupKProcs(this);
        }
    }

    for (auto& t: pTris) {
        if (t != nullptr) {
            t->setupKProcs(this, efflag());
        }
    }
    // Resolve all dependencies

    // DEBUG: vector holds all possible tetrahedrons,
    // but they have not necessarily been added to a compartment.
    for (auto& t: pTets) {
        if ((t != nullptr) && t->getInHost()) {
            t->setupDeps();
        }
    }

    // Vector allows for all compartments to be well-mixed, so
    // hold null-pointer for mesh compartments
    for (auto& wmv: pWmVols) {
        if ((wmv != nullptr) && wmv->getInHost()) {
            wmv->setupDeps();
        }
    }

    // DEBUG: vector holds all possible triangles, but
    // only patch triangles are filled
    for (auto& t: pTris) {
        if ((t != nullptr) && t->getInHost()) {
            t->setupDeps();
        }
    }
    // Create EField structures if EField is to be calculated
    if (efflag()) {
        _setupEField();
    }

    for (auto& tet: boundaryTets) {
        tet->setupBufferLocations();
    }
    for (auto& tri: boundaryTris) {
        tri->setupBufferLocations();
    }
    // just in case
    neighbHosts.erase(myRank);
    nNeighbHosts = neighbHosts.size();

    // construct remote molecule change buffers
    remoteChanges.clear();
    for (auto& neighbor: neighbHosts) {
        remoteChanges[neighbor] = {};
    }

    nEntries = pKProcs.size();
    diffSep = pDiffs.size();
    sdiffSep = pSDiffs.size();
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setupEField() {
    using math::point3d;
    using solver::efield::dVSolverBanded;
    using solver::efield::make_EField;

    //// Note to self: for now roughly following flow from original code in
    /// sim/controller.py. / code for setting up a mesh was in func_tetmesh
    /// constructor and called functions.

    AssertLog(efflag());

    switch (pEFoption) {
    case EF_DEFAULT:
    case EF_DV_BDSYS:
        pEField = make_EField<dVSolverBanded>();
        break;
#ifdef USE_PETSC
    case EF_DV_PETSC:
        using solver::efield::dVSolverPETSC;
        pEField = make_EField<dVSolverPETSC>(MPI_COMM_WORLD);
        break;
#endif
    default:
        ArgErrLog("Unsupported E-Field solver.");
    }

    // Give temperature a default value of 20c
    pTemp = 293.15;

    uint nmembs = mesh()->_countMembs();

    if (nmembs != 1) {
        std::ostringstream os;
        os << "Membrane potential solver requires only one ";
        os << "membrane description object per membrane but " << nmembs << " were found.";
        if (nmembs == 0) {
            os << "\nConsider using EF_NONE in the Simulation constructor instead.";
        }

        ArgErrLog(os.str());
    }

    tetmesh::Memb* memb = mesh()->_getMemb(0);
    AssertLog(memb != nullptr);

    // TODO: Decide what checks are needed for the membrane and implement them
    // here

    pEFNTets = memb->countVolTets();
    pEFNTris = memb->countTris();
    pEFNVerts = memb->countVerts();

    std::vector<vertex_id_t> EFTets(pEFNTets * 4);

    // All the triangles we will count here are the real membrane triangles,
    // virtual triangles will not require a capacitance.
    std::vector<vertex_id_t> EFTris(pEFNTris * 3);

    std::vector<double> EFVerts(pEFNVerts * 3);

    auto nverts = mesh()->countVertices();
    auto ntris = mesh()->countTris();
    auto ntets = mesh()->countTets();

    pEFVert_GtoL.container().resize(nverts);
    pEFTri_GtoL.container().resize(ntris);
    pEFTet_GtoL.container().resize(ntets);
    pEFTri_LtoG.container().resize(neftris());

    // Copy the data to local structures.

    auto const& membverts = memb->_getAllVertIndices();
    AssertLog(membverts.size() == nefverts());
    for (uint efv = 0; efv < nefverts(); ++efv) {
        const auto vertidx = membverts[efv];
        auto verttemp = mesh()->_getVertex(vertidx);
        uint efv2 = efv * 3;

        // CONVERTING TO MICRONS HERE. EFIELD OBJECT WILL NOT PERFORM THIS
        // CONVERSION
        verttemp *= 1.0e6;
        EFVerts[efv2] = verttemp[0];
        EFVerts[efv2 + 1] = verttemp[1];
        EFVerts[efv2 + 2] = verttemp[2];

        pEFVert_GtoL[vertidx] = vertex_id_t(efv);
    }

    const auto& membtets = memb->_getAllVolTetIndices();
    AssertLog(membtets.size() == neftets());
    for (uint eft = 0; eft < neftets(); ++eft) {
        auto const& tetidx = membtets[eft];
        const auto& tettemp = mesh()->_getTet(tetidx);
        uint eft2 = eft * 4;

        // Convert to indices used by EField object
        auto tv0 = pEFVert_GtoL[tettemp[0]];
        auto tv1 = pEFVert_GtoL[tettemp[1]];
        auto tv2 = pEFVert_GtoL[tettemp[2]];
        auto tv3 = pEFVert_GtoL[tettemp[3]];
        if (tv0.unknown() || tv1.unknown() || tv2.unknown() || tv3.unknown()) {
            std::ostringstream os;
            os << "Failed to create EField structures.";
            ProgErrLog(os.str());
        }

        EFTets[eft2] = tv0;
        EFTets[eft2 + 1] = tv1;
        EFTets[eft2 + 2] = tv2;
        EFTets[eft2 + 3] = tv3;

        pEFTet_GtoL[tetidx] = tetrahedron_local_id(eft);
    }

    const auto& membtris = memb->_getAllTriIndices();
    AssertLog(membtris.size() == neftris());

    pEFTris_vec.resize(pEFNTris);
    EFTrisV.container().resize(pEFNTris);

    EFTrisI_permuted.resize(pEFNTris);
    EFTrisI_idx.resize(pEFNTris);

    EFTrisI_offset.assign(nHosts, 0);
    EFTrisI_count.assign(nHosts, 0);

    std::vector<index_t> local_eftri_indices;

    for (uint eft = 0; eft < pEFNTris; ++eft) {
        auto triidx = membtris[eft];
        const auto& tritemp = mesh()->_getTri(triidx);
        uint eft2 = eft * 3;

        // Convert to indices used by EField object
        auto tv0 = pEFVert_GtoL[tritemp[0]];
        auto tv1 = pEFVert_GtoL[tritemp[1]];
        auto tv2 = pEFVert_GtoL[tritemp[2]];
        if (tv0.unknown() || tv1.unknown() || tv2.unknown()) {
            std::ostringstream os;
            os << "Failed to create EField structures.";
            ProgErrLog(os.str());
        }

        EFTris[eft2] = tv0;
        EFTris[eft2 + 1] = tv1;
        EFTris[eft2 + 2] = tv2;

        pEFTri_GtoL[triidx] = triangle_local_id(eft);
        pEFTri_LtoG[triangle_local_id(eft)] = triidx;

        // This is added now for quicker iteration during run()
        // Extremely important for larger meshes, orders of magnitude times faster
        Tri* tri_p = pTris[triidx];
        pEFTris_vec[eft] = tri_p;

        int tri_host = tri_p->getHost();
        ++EFTrisI_count[tri_host];
        if (myRank == tri_host) {
            local_eftri_indices.push_back(eft);
        }
    }

    const int* count_begin = EFTrisI_count.data();
    std::partial_sum(count_begin, count_begin + (nHosts - 1), EFTrisI_offset.data() + 1);

    AssertLog(local_eftri_indices.size() == static_cast<uint>(EFTrisI_count[myRank]));

    MPI_Allgatherv(local_eftri_indices.data(),
                   static_cast<int>(local_eftri_indices.size()),
                   MPI_STEPS_INDEX,
                   EFTrisI_idx.data(),
                   EFTrisI_count.data(),
                   EFTrisI_offset.data(),
                   MPI_STEPS_INDEX,
                   MPI_COMM_WORLD);

    pEField->initMesh(EFVerts,
                      EFTris,
                      EFTets,
                      memb->_getOpt_method(),
                      memb->_getOpt_file_name(),
                      memb->_getSearch_percent());

    // Triangles need to be set to some initial voltage, which they can read from
    // the Efield pointer. _setup() will read those voltages to initialise the
    // voltage-dependent reactions
    _refreshEFTrisV();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::saveMembOpt(std::string const& opt_file_name) {
    if (myRank != 0) {
        return;
    }
    if (!efflag()) {
        std::ostringstream os;
        os << "saveMembOpt method only available if running EField ";
        ArgErrLog(os.str());
    }

    pEField->saveOptimal(opt_file_name);
}

////////////////////////////////////////////////////////////////////////////////

uint TetOpSplitP::_addComp(solver::Compdef* cdef) {
    auto comp = new Comp(cdef);
    auto compidx = pComps.size();
    pComps.container().push_back(comp);
    return compidx;
}

////////////////////////////////////////////////////////////////////////////////

uint TetOpSplitP::_addPatch(solver::Patchdef* pdef) {
    auto patch = new Patch(pdef);
    auto patchidx = pPatches.size();
    pPatches.container().push_back(patch);
    return patchidx;
}

////////////////////////////////////////////////////////////////////////////////

uint TetOpSplitP::_addDiffBoundary(solver::DiffBoundarydef* dbdef) {
    auto diffb = new DiffBoundary(dbdef);
    auto dbidx = pDiffBoundaries.size();
    pDiffBoundaries.push_back(diffb);
    return dbidx;
}

////////////////////////////////////////////////////////////////////////////////

uint TetOpSplitP::_addSDiffBoundary(solver::SDiffBoundarydef* sdbdef) {
    auto sdiffb = new SDiffBoundary(sdbdef);
    auto sdbidx = pSDiffBoundaries.size();
    pSDiffBoundaries.push_back(sdiffb);
    return sdbidx;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_addTet(tetrahedron_global_id tetidx,
                          Comp* comp,
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
    auto localtet = new Tet(tetidx,
                            compdef,
                            vol,
                            a1,
                            a2,
                            a3,
                            a4,
                            d1,
                            d2,
                            d3,
                            d4,
                            tet0,
                            tet1,
                            tet2,
                            tet3,
                            myRank,
                            tetHosts.at(tetidx));
    AssertLog(pTets.at(tetidx) == nullptr);
    pTets[tetidx] = localtet;
    comp->addTet(localtet);

    // MPISTEPS
    localtet->setSolver(this);
    // MPISTEPS
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_addWmVol(solver::comp_global_id cidx, Comp* comp, double vol) {
    solver::Compdef* compdef = comp->def();
    auto* localtet =
        new WmVol(tetrahedron_global_id(cidx.get()), compdef, vol, myRank, wmHosts[cidx]);
    AssertLog(localtet != nullptr);
    AssertLog(cidx < pWmVols.size());
    pWmVols[cidx] = localtet;
    comp->addTet(localtet);

    // MPISTEPS
    localtet->setSolver(this);
    // MPISTEPS
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_addTri(triangle_global_id triidx,
                          Patch* patch,
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
    auto tri = new Tri(triidx,
                       patchdef,
                       area,
                       l0,
                       l1,
                       l2,
                       d0,
                       d1,
                       d2,
                       tinner,
                       touter,
                       tri0,
                       tri1,
                       tri2,
                       myRank,
                       triHosts.at(triidx));
    AssertLog(pTris.at(triidx) == nullptr);
    pTris[triidx] = tri;
    patch->addTri(tri);

    // MPISTEPS
    tri->setSolver(this);
    // MPISTEPS
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::reset() {
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
        if (tet == nullptr or !tet->getInHost()) {
            continue;
        }
        tet->reset();
    }

    for (auto const& wmvol: pWmVols) {
        if (wmvol == nullptr or !wmvol->getInHost()) {
            continue;
        }
        wmvol->reset();
    }

    for (auto const& tri: pTris) {
        if (tri == nullptr or !tri->getInHost()) {
            continue;
        }
        tri->reset();
    }

    auto ngroups = nGroups.size();
    for (uint i = 0; i < ngroups; i++) {
        free(nGroups[i]->indices);
        delete nGroups[i];
    }
    nGroups.clear();

    ngroups = pGroups.size();
    for (uint i = 0; i < ngroups; i++) {
        free(pGroups[i]->indices);
        delete pGroups[i];
    }
    pGroups.clear();

    pSum = 0.0;
    nSum = 0.0;
    pA0 = 0.0;

    reacExtent = 0.0;
    diffExtent = 0.0;
    nIteration = 0.0;

    statedef().resetTime();
    statedef().resetNSteps();
    _updateLocal();

    compTime = 0.0;
    syncTime = 0.0;
    idleTime = 0.0;

    efieldTime = 0.0;
    rdTime = 0.0;
    dataExchangeTime = 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::run(double endtime) {
    if (endtime < statedef().time()) {
        std::ostringstream os;
        os << "Endtime is before current simulation time";
        ArgErrLog(os.str());
    } else if (util::almost_equal(statedef().time(), endtime)) {
        statedef().setTime(endtime);
    } else {
        if (recomputeUpdPeriod) {
            _computeUpdPeriod();
        }
        if (efflag()) {
            _runWithEField(endtime);
        } else {
            _runWithoutEField(endtime);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_runWithoutEField(double endtime) {
    MPI_Barrier(MPI_COMM_WORLD);

    // This is the time (in seconds) to the next diffusion update. The upper
    // limit, so as to avoid systematic slowing of diffusion, is the inverse of
    // the highest local single-molecule diffusion rate of occupied tets.

    // This bool tracks if the update time has been aligned to the endtime. Trying
    // to do this by instead comparing two doubles can cause infinite loops.

    // decide globally
    bool aligned = false;

    double update_period = updPeriod;

    MPI_Request* requests = nullptr;

    // here we assume that all molecule counts have been updated so the rates are
    // accurate
    while (statedef().time() < endtime and not aligned) {
#ifdef MPI_PROFILING
        double timing_start = MPI_Wtime();
#endif
        Instrumentor::phase_begin("runWithoutEField -> Operator Split: SSA");

        double pre_ssa_time = statedef().time();
        // Update period may take us past the endtime- adjust if so
        if (pre_ssa_time + updPeriod > endtime) {
            update_period = endtime - pre_ssa_time;
            aligned = true;
        }

        // *********************** Operator Split: SSA
        // *********************************

        // Run SSA for the update period

        double cumulative_dt = 0.0;

        // Store a sequence of kprocs actually applied
        KProcPSet applied_ssa_kprocs;
        while (true) {
            KProc* kp = _getNext();
            if (kp == nullptr) {
                break;
            }
            double a0 = getA0();
            if (a0 == 0.0) {
                break;
            }
            double dt = rng()->getExp(a0);
            if (cumulative_dt + dt > update_period || statedef().time() + dt > endtime) {
                break;
            }
            cumulative_dt += dt;

            _executeStep(kp, dt, cumulative_dt);
            reacExtent += 1;

            applied_ssa_kprocs.insert(kp);
        }

        // Now for each process advance to jump time and get jump randomly
        // ****************************************************************************

        // *********************** Operator Split: Diffusion
        // ***************************

        // Apply diffusion after the update period

        Instrumentor::phase_end("runWithoutEField -> Operator Split: SSA");
#ifdef MPI_PROFILING
        double timing_end = MPI_Wtime();
        compTime += (timing_end - timing_start);
        timing_start = MPI_Wtime();
#endif

        // wait until previous loop finishes sending diffusion data
        if (requests != nullptr) {
            MPI_Waitall(nNeighbHosts, requests, MPI_STATUSES_IGNORE);
            delete[] requests;
        }

        // create new requests for this loop
        requests = new MPI_Request[nNeighbHosts];

        for (auto& neighbor: neighbHosts) {
            remoteChanges[neighbor].clear();
        }

#ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        idleTime += (timing_end - timing_start);
        timing_start = MPI_Wtime();
#endif
        Instrumentor::phase_begin("runWithoutEField -> Operator Split: Diffusion");

        // Track how many diffusion 'steps' we do, simply for bookkeeping
        uint nsteps = 0;

        // to reduce memory cost we use directions to retrieve the update list in
        // upd process
        std::vector<KProc*> applied_diffs;
        std::vector<int> directions;

        for (uint pos = 0; pos < diffSep; pos++) {
            Diff* d = pDiffs[pos];
            double rate = d->crData.rate;
            if (rate == 0) {
                continue;
            }
            // rate is the rate (scaled_dcst * population)
            double scaleddcst = d->getScaledDcst();

            // The number of molecules available for diffusion for this diffusion rule.
            // We use the rate as a way to get the population before any diffusion is applied
            // since it will only be updated after all the diffusions are processed.
            double population = std::round(rate / scaleddcst);

            // t1, AKA 'X', is a fractional number between 0 and 1: the update period
            // divided by the local mean single-molecule dwellperiod. This fraction
            // gives the mean proportion of molecules to diffuse.
            double t1 = update_period * scaleddcst;

            if (t1 >= 1.0) {
                t1 = 1.0;
            }

            // Calculate the occupancy, that is the integrated molecules over the
            // period (units s)
            double occupancy = d->getTet()->getPoolOccupancy(d->getLigLidx()) +
                               population *
                                   (update_period - d->getTet()->getLastUpdate(d->getLigLidx()));

            // n is, correctly, a binomial, but the binomial function requires
            // rounding to an integer.

            // occupancy/update_period gives the mean number of molecules during the
            // period
            double n_double = occupancy / update_period;

            // could be higher than those available - a source of error
            if (n_double > population) {
                n_double = population;
            }
            double n_int = std::floor(n_double);
            double n_frc = n_double - n_int;
            uint mean_n = static_cast<uint>(n_int);

            // deal linearly with the fraction
            if (n_frc > 0.0) {
                double rand01 = rng()->getUnfIE();
                if (rand01 < n_frc) {
                    mean_n++;
                }
            }

            // Find the binomial n
            uint nmolcs = rng()->getBinom(mean_n, t1);

            if (nmolcs == 0) {
                continue;
            }

            // we apply here
            if (nmolcs > diffApplyThreshold) {
                int direction = d->apply(rng(), nmolcs);
                if (applied_diffs.empty() or applied_diffs.back() != d or
                    directions.back() != direction) {
                    applied_diffs.push_back(d);
                    directions.push_back(direction);
                }
            } else {
                for (uint ai = 0; ai < nmolcs; ++ai) {
                    int direction = d->apply(rng());
                    if (applied_diffs.empty() or applied_diffs.back() != d or
                        directions.back() != direction) {
                        applied_diffs.push_back(d);
                        directions.push_back(direction);
                    }
                }
            }
            nsteps += nmolcs;
            diffExtent += nmolcs;
        }

        // surface diffusion

        for (uint pos = 0; pos < sdiffSep; pos++) {
            SDiff* d = pSDiffs[pos];
            double rate = d->crData.rate;

            if (rate == 0) {
                continue;
            }
            // rate is the rate (scaled_dcst * population)
            double scaleddcst = d->getScaledDcst();

            // The number of molecules available for diffusion for this diffusion rule.
            // We use the rate as a way to get the population before any diffusion is applied
            // since it will only be updated after all the diffusions are processed.
            double population = std::round(rate / scaleddcst);

            // t1, AKA 'X', is a fractional number between 0 and 1: the update period
            // divided by the local mean single-molecule dwellperiod. This fraction
            // gives the mean proportion of molecules to diffuse.
            double t1 = update_period * scaleddcst;

            if (t1 >= 1.0) {
                t1 = 1.0;
            }

            double occupancy = d->getTri()->getPoolOccupancy(d->getLigLidx()) +
                               population *
                                   (update_period - d->getTri()->getLastUpdate(d->getLigLidx()));

            // n is, correctly, a binomial, but the binomial function requires
            // rounding to an integer.

            // occupancy/update_period gives the mean number of molecules during the
            // period
            double n_double = occupancy / update_period;

            // could be higher than those available - a source of error
            if (n_double > population) {
                n_double = population;
            }
            double n_int = std::floor(n_double);
            double n_frc = n_double - n_int;
            uint mean_n = static_cast<uint>(n_int);

            // deal linearly with the fraction
            if (n_frc > 0.0) {
                double rand01 = rng()->getUnfIE();
                if (rand01 < n_frc) {
                    mean_n++;
                }
            }

            // Find the binomial n
            uint nmolcs = rng()->getBinom(mean_n, t1);

            if (nmolcs == 0) {
                continue;
            }
            // we apply here
            if (nmolcs > diffApplyThreshold) {
                int direction = d->apply(rng(), nmolcs);
                if (applied_diffs.empty() or applied_diffs.back() != d or
                    directions.back() != direction) {
                    applied_diffs.push_back(d);
                    directions.push_back(direction);
                }
            } else {
                for (uint ai = 0; ai < nmolcs; ++ai) {
                    int direction = d->apply(rng());
                    if (applied_diffs.empty() or applied_diffs.back() != d or
                        directions.back() != direction) {
                        applied_diffs.push_back(d);
                        directions.push_back(direction);
                    }
                }
            }
            nsteps += nmolcs;
            diffExtent += nmolcs;
        }

        Instrumentor::phase_end("runWithoutEField -> Operator Split: Diffusion");
#ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        compTime += (timing_end - timing_start);
#endif

        Instrumentor::phase_begin("runWithoutEField -> _remoteSyncAndUpdate");
        _remoteSyncAndUpdate(requests, applied_diffs, directions);
        Instrumentor::phase_end("runWithoutEField -> _remoteSyncAndUpdate");

// *********************** Operator Split: SSA *********************************
#ifdef MPI_PROFILING
        starttime = MPI_Wtime();
#endif
        Instrumentor::phase_begin("runWithoutEField -> Operator Split: SSA");

        for (auto const& akp: applied_ssa_kprocs) {
            akp->resetOccupancies();
        }

        // by the end of the ssa iteration,
        // force the state time to be aligned with the user-defined endtime
        if (aligned) {
            statedef().setTime(endtime);
        } else {
            statedef().setTime(pre_ssa_time + update_period);
        }
        if (nsteps > 0) {
            statedef().incNSteps(nsteps);
        }

        nIteration += 1;

        Instrumentor::phase_end("runWithoutEField -> Operator Split: SSA");
#ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        compTime += (timing_end - timing_start);
#endif
    }
    if (requests != nullptr) {
        MPI_Waitall(nNeighbHosts, requests, MPI_STATUSES_IGNORE);
        delete[] requests;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_refreshEFTrisV() {
    for (auto tlidx: triangle_local_id::range(pEFNTris)) {
        EFTrisV[tlidx] = pEField->getTriV(tlidx);
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_runWithEField(double endtime) {
#ifdef MPI_PROFILING
    double timing_start = MPI_Wtime();
#endif

#ifdef MPI_PROFILING
    double timing_end = MPI_Wtime();
    efieldTime += (timing_end - timing_start);
#endif

    while (statedef().time() < endtime) {
#ifdef MPI_PROFILING
        timing_start = MPI_Wtime();
#endif
        Instrumentor::phase_begin("runWithEField -> rd");

        double t0 = statedef().time();

        double t1 = std::min(t0 + pEFDT, endtime);
        if (util::almost_equal(t1, endtime)) {
            t1 = endtime;
        }
        _runWithoutEField(t1);

        Instrumentor::phase_end("runWithEField -> rd");
#ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        rdTime += (timing_end - timing_start);
#endif

#ifdef MPI_PROFILING
        timing_start = MPI_Wtime();
#endif
        Instrumentor::phase_begin("runWithEField -> efield");
        // update host-local currents
        int i_begin = EFTrisI_offset[myRank];
        int i_end = i_begin + EFTrisI_count[myRank];

        double sttime = statedef().time();
        double real_ef_dt = sttime - t0;
        for (int i = i_begin; i < i_end; ++i) {
            auto tlidx = EFTrisI_idx[i];
            EFTrisI_permuted[i] =
                pEFTris_vec[tlidx.get()]->computeI(EFTrisV[tlidx], real_ef_dt, sttime, efdt());
        }

        Instrumentor::phase_end("runWithEField -> efield");
#ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        efieldTime += (timing_end - timing_start);
#endif

#ifdef MPI_PROFILING
        timing_start = MPI_Wtime();
#endif
        MPI_Allgatherv(MPI_IN_PLACE,
                       0,
                       MPI_DATATYPE_NULL,
                       EFTrisI_permuted.data(),
                       EFTrisI_count.data(),
                       EFTrisI_offset.data(),
                       MPI_DOUBLE,
                       MPI_COMM_WORLD);

#ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        dataExchangeTime += (timing_end - timing_start);
#endif

#ifdef MPI_PROFILING
        timing_start = MPI_Wtime();
#endif
        Instrumentor::phase_begin("runWithEField -> efield");

        for (uint i = 0; i < pEFNTris; i++) {
            pEField->setTriI(EFTrisI_idx[i], EFTrisI_permuted[i]);
        }

        pEField->advance(real_ef_dt);
        _refreshEFTrisV();

        Instrumentor::phase_end("runWithEField -> efield");
#ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        efieldTime += (timing_end - timing_start);
#endif

#ifdef MPI_PROFILING
        timing_start = MPI_Wtime();
#endif
        Instrumentor::phase_begin("runWithEField -> kproc update");

        _updateLocal(pVdepKProcs);

        Instrumentor::phase_end("runWithEField -> kproc update");
#ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        rdTime += (timing_end - timing_start);
#endif
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////

void TetOpSplitP::advance(double adv) {
    if (adv < 0.0) {
        std::ostringstream os;
        os << "Time to advance cannot be negative";
        ArgErrLog(os.str());
    }

    double endtime = statedef().time() + adv;
    run(endtime);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::step() {
    std::ostringstream os;
    os << "This function is not available for this solver!";
    NotImplErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::getTime() const {
    return statedef().time();
}

////////////////////////////////////////////////////////////////////////

uint TetOpSplitP::getNSteps() const {
    return statedef().nsteps();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setTime(double time) {
    statedef().setTime(time);
}

////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setNSteps(uint nsteps) {
    statedef().setNSteps(nsteps);
}

////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setTemp(double t) {
    if (efflag() == false) {
        std::ostringstream os;
        os << "\nWARNING: Temperature set in simulation without membrane ";
        os << "potential calculation will be ignored.\n";
        // CLOG(INFO, "general_log") << os << std::endl;
    }
    AssertLog(t >= 0.0);
    pTemp = t;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getCompVol(solver::comp_global_id cidx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    Comp* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    return comp->vol();
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getCompSpecCount(solver::comp_global_id cidx,
                                      solver::spec_global_id sidx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(statedef().countComps() == pComps.size());
    Comp* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::spec_local_id slidx = comp->def()->specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    uint local_total_count = 0;
    uint total_count = 0;
    for (auto& t: comp->tets()) {
        if (t->getInHost()) {
            local_total_count += t->pools()[slidx];
        }
    }
    MPI_Allreduce(&local_total_count, &total_count, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    return total_count;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setCompSpecCount(solver::comp_global_id cidx,
                                    solver::spec_global_id sidx,
                                    double n) {
    MPI_Barrier(MPI_COMM_WORLD);
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(n >= 0.0);
    AssertLog(statedef().countComps() == pComps.size());
    Comp* comp = _comp(cidx);
    const solver::spec_local_id slidx = comp->def()->specG2L(sidx);

    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // only do the distribution in rank 0
    // then bcast to other ranks

    if (myRank == 0) {
        // functions for distribution:
        auto set_count = [slidx](WmVol* tet, uint c) { tet->setCount(slidx, c); };
        auto inc_count = [slidx](WmVol* tet) { tet->incCount(slidx, 1, 0.0, true); };
        auto weight = [](const WmVolPVecCI& tet) { return (*tet)->vol(); };

        util::distribute_quantity<uint>(n,
                                        comp->tets().begin(),
                                        comp->tets().end(),
                                        weight,
                                        set_count,
                                        inc_count,
                                        *rng(),
                                        comp->def()->vol());
    }

    std::vector<uint> counts(comp->countTets());

    if (myRank == 0) {
        std::transform(comp->tets().begin(),
                       comp->tets().end(),
                       counts.begin(),
                       [slidx](const WmVolP& t) { return t->pools()[slidx]; });
    }

    MPI_Bcast(counts.data(), counts.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    uint curr_pos = 0;
    for (auto const& t: comp->tets()) {
        if (myRank != 0) {
            t->setCount(slidx, counts[curr_pos]);
        }
        _updateSpec(*t, sidx);
        curr_pos++;
    }
    _updateSum();
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getCompSpecAmount(solver::comp_global_id cidx,
                                       solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getCompSpecCount(cidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setCompSpecAmount(solver::comp_global_id cidx,
                                     solver::spec_global_id sidx,
                                     double a) {
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setCompSpecCount(cidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getCompSpecConc(solver::comp_global_id cidx,
                                     solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getCompSpecCount(cidx, sidx);
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    double vol = comp->vol();
    return count / (1.0e3 * vol * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setCompSpecConc(solver::comp_global_id cidx,
                                   solver::spec_global_id sidx,
                                   double c) {
    AssertLog(c >= 0.0);
    AssertLog(cidx < statedef().countComps());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    double count = c * (1.0e3 * comp->vol() * math::AVOGADRO);
    // the following method does all the necessary argument checking
    _setCompSpecCount(cidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getCompSpecClamped(solver::comp_global_id cidx,
                                      solver::spec_global_id sidx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(statedef().countComps() == pComps.size());
    Comp* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::spec_local_id lsidx = comp->def()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    bool local_clamped = true;
    for (auto const& t: comp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        if (!t->clamped(lsidx)) {
            local_clamped = false;
        }
    }

    bool global_clamped = false;

    MPI_Allreduce(&local_clamped, &global_clamped, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);

    return global_clamped;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setCompSpecClamped(solver::comp_global_id cidx,
                                      solver::spec_global_id sidx,
                                      bool b) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(statedef().countComps() == pComps.size());
    Comp* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::spec_local_id lsidx = comp->def()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // Set the flag in def object, though this may not be necessary
    comp->def()->setClamped(lsidx, b);

    for (auto const& t: comp->tets()) {
        t->setClamped(lsidx, b);
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getCompReacK(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    AssertLog(statedef().countComps() == pComps.size());
    Comp* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = comp->def()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // We're just returning the default value for this comp, individual
    // tets may have different Kcsts set individually
    return (comp->def()->kcst(lridx));
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setCompReacK(solver::comp_global_id cidx,
                                solver::reac_global_id ridx,
                                double kf) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    AssertLog(statedef().countComps() == pComps.size());
    AssertLog(kf >= 0.0);
    Comp* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = comp->def()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // First set the default value for the comp
    comp->def()->setKcst(lridx, kf);

    // Now update all tetrahedra in this comp
    for (auto const& t: comp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        t->reac(lridx).setKcst(kf);
    }

    // Rates have changed
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getCompReacActive(solver::comp_global_id cidx,
                                     solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    AssertLog(statedef().countComps() == pComps.size());
    Comp* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = comp->def()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    bool local_active = true;

    for (auto const& t: comp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        if (t->reac(lridx).inactive()) {
            local_active = false;
        }
    }

    bool global_active = false;

    MPI_Allreduce(&local_active, &global_active, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
    return global_active;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setCompReacActive(solver::comp_global_id cidx,
                                     solver::reac_global_id ridx,
                                     bool a) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    AssertLog(statedef().countComps() == pComps.size());
    Comp* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = comp->def()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // Set the default value for the comp, though this is not entirely
    // necessary
    comp->def()->setActive(lridx, a);

    for (auto const& t: comp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        t->reac(lridx).setActive(a);
    }
    // It's cheaper to just recompute everything.
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getCompDiffD(solver::comp_global_id cidx, solver::diff_global_id didx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(didx < statedef().countDiffs());
    AssertLog(statedef().countComps() == pComps.size());
    Comp* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::diff_local_id ldidx = comp->def()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // We're just returning the default value for this comp, individual
    // tets may have different Dcsts set individually
    return (comp->def()->dcst(ldidx));
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setCompDiffD(solver::comp_global_id cidx,
                                solver::diff_global_id didx,
                                double dk) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(didx < statedef().countDiffs());
    AssertLog(statedef().countComps() == pComps.size());
    AssertLog(dk >= 0.0);
    Comp* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::diff_local_id ldidx = comp->def()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in compartment.\n";
        ArgErrLog(os.str());
    }
    recomputeUpdPeriod = true;
    // First set the default value for the comp
    comp->def()->setDcst(ldidx, dk);

    // Now update all tets in this comp
    for (auto const& t: comp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        if (auto tet = dynamic_cast<Tet*>(t)) {
            tet->diff(ldidx).setDcst(dk);
        } else {
            std::ostringstream os;
            os << "Cannot change diffusion constant in well-mixed compartment.";
            ArgErrLog(os.str());
        }
    }

    // Rates have changed
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getCompDiffActive(solver::comp_global_id cidx,
                                     solver::diff_global_id didx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(didx < statedef().countDiffs());
    AssertLog(statedef().countComps() == pComps.size());
    Comp* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::diff_local_id ldidx = comp->def()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    bool local_active = true;

    for (auto const& t: comp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        // WmVol * wmcomp = (*t);
        if (auto tet = dynamic_cast<Tet*>(t)) {
            if (tet->diff(ldidx).inactive()) {
                local_active = false;
            }
        } else {
            std::ostringstream os;
            os << "Diffusion activation not defined in well-mixed compartment.\n";
            ArgErrLog(os.str());
        }
    }
    bool global_active = false;
    MPI_Allreduce(&local_active, &global_active, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
    return global_active;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setCompDiffActive(solver::comp_global_id cidx,
                                     solver::diff_global_id didx,
                                     bool act) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(didx < statedef().countDiffs());
    AssertLog(statedef().countComps() == pComps.size());
    Comp* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::diff_local_id ldidx = comp->def()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    for (auto const& t: comp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        // WmVol * wmcomp = (*t);
        if (auto tet = dynamic_cast<Tet*>(t)) {
            tet->diff(ldidx).setActive(act);
        } else {
            std::ostringstream os;
            os << "Cannot change diffusion constant in well-mixed compartment.\n";
            ArgErrLog(os.str());
        }
    }

    recomputeUpdPeriod = true;

    // It's cheaper to just recompute everything.
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getPatchArea(solver::patch_global_id pidx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    Patch* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    return patch->area();
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getPatchSpecCount(solver::patch_global_id pidx,
                                       solver::spec_global_id sidx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(statedef().countPatches() == pPatches.size());
    Patch* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::spec_local_id slidx = patch->def()->specG2L(sidx);
    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    uint local_total_count = 0;
    uint total_count = 0;
    for (auto const& t: patch->tris()) {
        if (t->getInHost()) {
            local_total_count += t->pools()[slidx];
        }
    }
    MPI_Allreduce(&local_total_count, &total_count, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

    return total_count;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setPatchSpecCount(solver::patch_global_id pidx,
                                     solver::spec_global_id sidx,
                                     double n) {
    MPI_Barrier(MPI_COMM_WORLD);
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(statedef().countPatches() == pPatches.size());
    AssertLog(n >= 0.0);
    Patch* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::spec_local_id slidx = patch->def()->specG2L(sidx);

    if (slidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // only do the distribution in rank 0
    // then bcast to other ranks

    if (myRank == 0) {
        // functions for distribution:
        auto set_count = [slidx](Tri* tri, uint c) { tri->setCount(slidx, c); };
        auto inc_count = [slidx](Tri* tri) { tri->incCount(slidx, 1, 0.0, true); };
        auto weight = [](const TriPVecCI& tri) { return (*tri)->area(); };

        util::distribute_quantity<uint>(n,
                                        patch->tris().begin(),
                                        patch->tris().end(),
                                        weight,
                                        set_count,
                                        inc_count,
                                        *rng(),
                                        patch->def()->area());
    }
    std::vector<uint> counts(patch->countTris());
    if (myRank == 0) {
        std::transform(patch->tris().begin(),
                       patch->tris().end(),
                       counts.begin(),
                       [slidx](const TriP& t) { return t->pools()[slidx]; });
    }

    MPI_Bcast(counts.data(), counts.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    uint curr_pos = 0;
    for (auto const& t: patch->tris()) {
        if (myRank != 0) {
            // set count only don't need sync
            t->setCount(slidx, counts[curr_pos]);
        }
        _updateSpec(*t, sidx);
        curr_pos++;
    }
    _updateSum();
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getPatchSpecAmount(solver::patch_global_id pidx,
                                        solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getPatchSpecCount(pidx, sidx);
    return (count / math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setPatchSpecAmount(solver::patch_global_id pidx,
                                      solver::spec_global_id sidx,
                                      double a) {
    AssertLog(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setPatchSpecCount(pidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getPatchSpecClamped(solver::patch_global_id pidx,
                                       solver::spec_global_id sidx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(statedef().countPatches() == pPatches.size());
    Patch* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::spec_local_id lsidx = patch->def()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    bool local_clamped = true;

    for (auto const& t: patch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        if (t->clamped(lsidx) == false) {
            local_clamped = false;
        }
    }
    bool global_clamped = false;

    MPI_Allreduce(&local_clamped, &global_clamped, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);

    return global_clamped;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setPatchSpecClamped(solver::patch_global_id pidx,
                                       solver::spec_global_id sidx,
                                       bool buf) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(statedef().countPatches() == pPatches.size());
    Patch* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::spec_local_id lsidx = patch->def()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // Set the flag in def object for consistency, though this is not
    // entirely necessary
    patch->def()->setClamped(lsidx, buf);

    for (auto const& t: patch->tris()) {
        t->setClamped(lsidx, buf);
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getPatchSReacK(solver::patch_global_id pidx,
                                    solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    AssertLog(statedef().countPatches() == pPatches.size());
    Patch* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // We're just returning the default value for this patch, individual
    // triangles may have different Kcsts set
    return (patch->def()->kcst(lsridx));
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setPatchSReacK(solver::patch_global_id pidx,
                                  solver::sreac_global_id ridx,
                                  double kf) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    AssertLog(statedef().countPatches() == pPatches.size());
    AssertLog(kf >= 0.0);
    Patch* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // First set the default values for this patch
    patch->def()->setKcst(lsridx, kf);

    // Now update all triangles in this patch
    for (auto const& t: patch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        t->sreac(lsridx).setKcst(kf);
    }

    // Rates have changed
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getPatchSReacActive(solver::patch_global_id pidx,
                                       solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    AssertLog(statedef().countPatches() == pPatches.size());
    Patch* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }
    bool local_active = true;
    for (auto const& t: patch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        if (t->sreac(lsridx).inactive()) {
            local_active = false;
        }
    }
    bool global_active = false;
    MPI_Allreduce(&local_active, &global_active, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
    return global_active;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setDiffBoundarySpecDiffusionActive(solver::diffboundary_global_id dbidx,
                                                      solver::spec_global_id sidx,
                                                      bool act) {
    AssertLog(dbidx < statedef().countDiffBoundaries());
    AssertLog(sidx < statedef().countSpecs());

    // Need to do two things:
    // 1) check if the species is defined in both compartments conencted
    // by the diffusion boundary
    // 2) loop over all tetrahedrons around the diff boundary and then the
    // diffusion rules and activate diffusion if the diffusion rule
    // relates to this species

    DiffBoundary* diffb = _diffboundary(dbidx);
    Comp* compA = diffb->compA();
    Comp* compB = diffb->compB();

    solver::spec_local_id lsidxA = compA->def()->specG2L(sidx);
    solver::spec_local_id lsidxB = compB->def()->specG2L(sidx);

    if (lsidxA.unknown() || lsidxB.unknown()) {
        std::ostringstream os;
        os << "Species undefined in compartments connected by diffusion "
              "boundary.\n";
        ArgErrLog(os.str());
    }

    auto const& bdtets = diffb->getTets();
    auto const& bdtetsdir = diffb->getTetDirection();

    // Have to use indices rather than iterator because need access to the
    // tet direction
    auto ntets = bdtets.size();

    for (auto bdt = 0u; bdt != ntets; ++bdt) {
        Tet* tet = _tet(bdtets[bdt]);
        if (!tet->getInHost()) {
            continue;
        }
        auto direction = bdtetsdir[bdt];
        AssertLog(direction < 4);

        // Each diff kproc then has access to the species through it's defined
        // parent
        auto ndiffs = tet->compdef()->countDiffs();
        for (auto d: solver::diff_local_id::range(ndiffs)) {
            Diff& diff = tet->diff(d);
            // sidx is the global species index; so is the lig() return from diffdef
            auto specgidx = diff.def()->lig();
            if (specgidx == sidx) {
                diff.setDiffBndActive(direction, act);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getDiffBoundarySpecDiffusionActive(solver::diffboundary_global_id dbidx,
                                                      solver::spec_global_id sidx) const {
    AssertLog(dbidx < statedef().countDiffBoundaries());
    AssertLog(sidx < statedef().countSpecs());

    DiffBoundary* diffb = _diffboundary(dbidx);
    Comp* compA = diffb->compA();
    Comp* compB = diffb->compB();

    solver::spec_local_id lsidxA = compA->def()->specG2L(sidx);
    solver::spec_local_id lsidxB = compB->def()->specG2L(sidx);

    if (lsidxA.unknown() or lsidxB.unknown()) {
        std::ostringstream os;
        os << "Species undefined in compartments connected by diffusion "
              "boundary.\n";
        ArgErrLog(os.str());
    }

    auto const& bdtets = diffb->getTets();
    auto const& bdtetsdir = diffb->getTetDirection();

    // Have to use indices rather than iterator because need access to the
    // tet direction

    bool local_active = true;

    const auto ntets = bdtets.size();

    for (auto bdt = 0u; bdt != ntets; ++bdt) {
        Tet* tet = _tet(bdtets[bdt]);
        if (!tet->getInHost()) {
            continue;
        }
        uint direction = bdtetsdir[bdt];
        AssertLog(direction < 4);

        // Each diff kproc then has access to the species through it's defined
        // parent
        auto ndiffs = tet->compdef()->countDiffs();
        for (auto d: solver::diff_local_id::range(ndiffs)) {
            Diff& diff = tet->diff(d);
            // sidx is the global species index; so is the lig() return from diffdef
            auto specgidx = diff.def()->lig();
            if (specgidx == sidx) {
                // Just need to check the first one
                if (diff.getDiffBndActive(direction)) {
                    local_active = true;
                } else {
                    local_active = false;
                }
                break;
            }
        }
    }
    bool global_active = false;
    MPI_Allreduce(&local_active, &global_active, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
    return global_active;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setDiffBoundarySpecDcst(solver::diffboundary_global_id dbidx,
                                           solver::spec_global_id sidx,
                                           double dcst,
                                           solver::comp_global_id direction_comp) {
    AssertLog(dbidx < statedef().countDiffBoundaries());
    AssertLog(sidx < statedef().countSpecs());

    DiffBoundary* diffb = _diffboundary(dbidx);
    Comp* compA = diffb->compA();
    Comp* compB = diffb->compB();

    solver::spec_local_id lsidxA = compA->def()->specG2L(sidx);
    solver::spec_local_id lsidxB = compB->def()->specG2L(sidx);

    if (lsidxA.unknown() || lsidxB.unknown()) {
        std::ostringstream os;
        os << "Species undefined in compartments connected by diffusion "
              "boundary.\n";
        ArgErrLog(os.str());
    }
    recomputeUpdPeriod = true;

    solver::Compdef* dirc_compdef = nullptr;
    if (direction_comp.valid()) {
        dirc_compdef = _comp(direction_comp)->def();
    }

    auto const& bdtets = diffb->getTets();
    auto const& bdtetsdir = diffb->getTetDirection();

    const auto ntets = bdtets.size();

    for (auto bdt = 0u; bdt != ntets; ++bdt) {
        Tet* tet = _tet(bdtets[bdt]);
        if (!tet->getInHost()) {
            continue;
        }
        // if tet compdef equals to dirc_compdef,
        // it is the desination tet so diff should not be changed
        // nullptr (bidirection) and source tet are both different
        // fromdirc_compdef
        if (dirc_compdef == tet->compdef()) {
            continue;
        }
        auto direction = bdtetsdir[bdt];
        AssertLog(direction < 4);

        // Each diff kproc then has access to the species through it's defined
        // parent
        auto ndiffs = tet->compdef()->countDiffs();
        for (auto d: solver::diff_local_id::range(ndiffs)) {
            Diff& diff = tet->diff(d);
            // sidx is the global species index; so is the lig() return from diffdef
            auto specgidx = diff.def()->lig();
            if (specgidx == sidx) {
                diff.setDiffBndActive(direction, true);
                diff.setDirectionDcst(direction, dcst);
                _updateElement(&diff);
            }
        }
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setPatchSReacActive(solver::patch_global_id pidx,
                                       solver::sreac_global_id ridx,
                                       bool a) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    AssertLog(statedef().countPatches() == pPatches.size());
    Patch* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // First set the flags in def object for consistency, though this is
    // not entirely necessary for this solver
    patch->def()->setActive(lsridx, a);

    for (auto const& t: patch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        t->sreac(lsridx).setActive(a);
    }
    // It's cheaper to just recompute everything.
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getPatchVDepSReacActive(solver::patch_global_id pidx,
                                           solver::vdepsreac_global_id vsridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(vsridx < statedef().countVDepSReacs());
    AssertLog(statedef().countPatches() == pPatches.size());
    Patch* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::vdepsreac_local_id lvsridx = patch->def()->vdepsreacG2L(vsridx);
    if (lvsridx.unknown()) {
        std::ostringstream os;
        os << "Voltage-dependent surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    bool local_active = true;

    for (auto const& t: patch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        if (t->vdepsreac(lvsridx).inactive()) {
            local_active = false;
        }
    }
    bool global_active = false;
    MPI_Allreduce(&local_active, &global_active, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
    return global_active;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setPatchVDepSReacActive(solver::patch_global_id pidx,
                                           solver::vdepsreac_global_id vsridx,
                                           bool a) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(vsridx < statedef().countVDepSReacs());
    AssertLog(statedef().countPatches() == pPatches.size());
    Patch* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::vdepsreac_local_id lvsridx = patch->def()->vdepsreacG2L(vsridx);
    if (lvsridx.unknown()) {
        std::ostringstream os;
        os << "Voltage-dependent surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // Not necessary and not possible to set the flags in def object

    for (auto const& t: patch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        t->vdepsreac(lvsridx).setActive(a);
    }
    // It's cheaper to just recompute everything.
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setSDiffBoundarySpecDiffusionActive(solver::sdiffboundary_global_id sdbidx,
                                                       solver::spec_global_id sidx,
                                                       bool act) {
    AssertLog(sdbidx < statedef().countSDiffBoundaries());
    AssertLog(sidx < statedef().countSpecs());

    // Need to do two things:
    // 1) check if the species is defined in both patches connected
    // by the surface diffusion boundary
    // 2) loop over all triangles around the sdiff boundary and then the
    // diffusion rules and activate diffusion if the diffusion rule
    // relates to this species

    SDiffBoundary* sdiffb = _sdiffboundary(sdbidx);
    Patch* patchA = sdiffb->patchA();
    Patch* patchB = sdiffb->patchB();

    solver::spec_local_id lsidxA = patchA->def()->specG2L(sidx);
    solver::spec_local_id lsidxB = patchB->def()->specG2L(sidx);

    if (lsidxA.unknown() || lsidxB.unknown()) {
        std::ostringstream os;
        os << "Species undefined in patches connected by surface diffusion "
              "boundary.\n";
        ArgErrLog(os.str());
    }

    auto const& sbdtris = sdiffb->getTris();
    auto const& sbdtrisdir = sdiffb->getTriDirection();

    // Have to use indices rather than iterator because need access to the
    // tri direction
    auto ntris = sbdtris.size();

    for (uint sbdt = 0; sbdt != ntris; ++sbdt) {
        Tri* tri = _tri(sbdtris[sbdt]);
        if (!tri->getInHost()) {
            continue;
        }
        uint direction = sbdtrisdir[sbdt];
        AssertLog(direction < 3);

        // Each sdiff kproc then has access to the species through its defined
        // parent
        uint nsdiffs = tri->patchdef()->countSurfDiffs();
        for (auto sd: solver::surfdiff_local_id::range(nsdiffs)) {
            SDiff& sdiff = tri->sdiff(sd);
            // sidx is the global species index; so is the lig() return from diffdef
            solver::spec_global_id specgidx = sdiff.def()->lig();
            if (specgidx == sidx) {
                sdiff.setSDiffBndActive(direction, act);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getSDiffBoundarySpecDiffusionActive(solver::sdiffboundary_global_id sdbidx,
                                                       solver::spec_global_id sidx) const {
    AssertLog(sdbidx < statedef().countSDiffBoundaries());
    AssertLog(sidx < statedef().countSpecs());

    SDiffBoundary* sdiffb = _sdiffboundary(sdbidx);
    Patch* patchA = sdiffb->patchA();
    Patch* patchB = sdiffb->patchB();

    solver::spec_local_id lsidxA = patchA->def()->specG2L(sidx);
    solver::spec_local_id lsidxB = patchB->def()->specG2L(sidx);

    if (lsidxA.unknown() || lsidxB.unknown()) {
        std::ostringstream os;
        os << "Species undefined in patches connected by surface diffusion "
              "boundary.\n";
        ArgErrLog(os.str());
    }

    auto const& sbdtris = sdiffb->getTris();
    auto const& sbdtrisdir = sdiffb->getTriDirection();

    // Have to use indices rather than iterator because need access to the
    // tet direction

    short local_active = 1;  // true

    uint ntris = sbdtris.size();

    for (uint sbdt = 0; sbdt != ntris; ++sbdt) {
        Tri* tri = _tri(sbdtris[sbdt]);
        if (!tri->getInHost()) {
            continue;
        }
        uint direction = sbdtrisdir[sbdt];
        AssertLog(direction < 3);

        // Each sdiff kproc then has access to the species through its defined
        // parent
        uint nsdiffs = tri->patchdef()->countSurfDiffs();
        for (auto sd: solver::surfdiff_local_id::range(nsdiffs)) {
            SDiff& sdiff = tri->sdiff(sd);
            // sidx is the global species index; so is the lig() return from diffdef
            solver::spec_global_id specgidx = sdiff.def()->lig();
            if (specgidx == sidx) {
                // Just need to check the first one
                if (sdiff.getSDiffBndActive(direction)) {
                    local_active = 1;
                } else {
                    local_active = 0;
                }
                break;
            }
        }
    }
    short global_active = 0;
    MPI_Allreduce(&local_active, &global_active, 1, MPI_SHORT, MPI_LAND, MPI_COMM_WORLD);
    return global_active != 0;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setSDiffBoundarySpecDcst(solver::sdiffboundary_global_id sdbidx,
                                            solver::spec_global_id sidx,
                                            double dcst,
                                            solver::patch_global_id direction_patch) {
    SDiffBoundary* sdiffb = _sdiffboundary(sdbidx);
    Patch* patchA = sdiffb->patchA();
    Patch* patchB = sdiffb->patchB();

    solver::spec_local_id lsidxA = patchA->def()->specG2L(sidx);
    solver::spec_local_id lsidxB = patchB->def()->specG2L(sidx);

    if (lsidxA.unknown() || lsidxB.unknown()) {
        std::ostringstream os;
        os << "Species undefined in patches connected by surface diffusion "
              "boundary.\n";
        ArgErrLog(os.str());
    }

    recomputeUpdPeriod = true;

    solver::Patchdef* dirp_patchdef = nullptr;
    if (direction_patch.valid()) {
        dirp_patchdef = _patch(direction_patch)->def();
    }

    auto const& sbdtris = sdiffb->getTris();
    auto const& sbdtrisdir = sdiffb->getTriDirection();

    uint ntris = sbdtris.size();

    for (uint sbdt = 0; sbdt != ntris; ++sbdt) {
        Tri* tri = _tri(sbdtris[sbdt]);

        if (!tri->getInHost()) {
            continue;
        }

        if (dirp_patchdef == tri->patchdef()) {
            continue;
        }
        uint direction = sbdtrisdir[sbdt];
        AssertLog(direction < 3);

        // Each diff kproc then has access to the species through its defined parent
        uint nsdiffs = tri->patchdef()->countSurfDiffs();
        for (auto sd: solver::surfdiff_local_id::range(nsdiffs)) {
            SDiff& sdiff = tri->sdiff(sd);
            // sidx is the global species index; so is the lig() return from diffdef
            solver::spec_global_id specgidx = sdiff.def()->lig();
            if (specgidx == sidx) {
                // The following function will automatically activate diffusion
                // in this direction if necessary
                sdiff.setDirectionDcst(direction, dcst);
                _updateElement(&sdiff);
            }
        }
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

solver::kproc_global_id TetOpSplitP::addKProc(KProc* kp, bool Vdep) {
    SchedIDX nidx(static_cast<uint>(pKProcs.size()));  // because pKProcs.size() is ulong
    pKProcs.push_back(kp);
    if (Vdep) {
        pVdepKProcs.push_back(kp);
    }
    return nidx;
}

////////////////////////////////////////////////////////////////////////////////

KProc* TetOpSplitP::_getNext() const {
    AssertLog(pA0 >= 0.0);
    // Quick check to see whether nothing is there.
    if (pA0 == 0.0) {
        return nullptr;
    }

    double selector = pA0 * rng()->getUnfII();

    double partial_sum = 0.0;

    auto n_neg_groups = nGroups.size();
    auto n_pos_groups = pGroups.size();

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
    for (auto i = n_pos_groups - 1; i != UINT_MAX; i--) {
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

    for (auto i = n_neg_groups - 1; i != UINT_MAX; i--) {
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

void TetOpSplitP::_executeStep(KProc* kp, double dt, double period) {
    kp->apply(rng(), dt, statedef().time(), period);
    statedef().incTime(dt);

    // as in 0.6.1 reaction and surface reaction only require updates of local
    // KProcs, it may change if VDepSurface reaction is added in the future
    std::vector<KProc*> const& upd = kp->getLocalUpdVec();
    _updateLocal(upd);
    statedef().incNSteps(1);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_updateSpec(WmVol& tet, solver::spec_global_id spec_gidx) {
    // NOTE: this function does not update the Sum of popensity, _updateSum() is
    // required after calling it.
    AssertLog(_getTetSpecDefined(tet.idx(), spec_gidx));

    if (!tet.getInHost()) {
        return;
    }

    KProcPSet updset;

    // Loop over tet.
    uint nkprocs = tet.countKProcs();

    for (uint k = 0; k < nkprocs; k++) {
        if (tet.KProcDepSpecTet(k, &tet, spec_gidx)) {
            updset.insert(tet.getKProc(k));
        }
    }

    for (auto const& tri: tet.nexttris()) {
        if (tri == nullptr) {
            continue;
        }
        nkprocs = tri->countKProcs();
        for (uint sk = 0; sk < nkprocs; sk++) {
            if (tri->KProcDepSpecTet(sk, &tet, spec_gidx)) {
                updset.insert(tri->getKProc(sk));
            }
        }
    }

    _updateLocal(updset);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_updateSpec(Tri& tri, solver::spec_global_id spec_gidx) {
    // NOTE: this function does not update the Sum of popensity, _updateSum() is
    // required after calling it.
    AssertLog(_getTriSpecDefined(tri.idx(), spec_gidx));

    if (!tri.getInHost()) {
        return;
    }

    KProcPSet updset;

    uint nkprocs = tri.countKProcs();

    for (uint sk = 0; sk < nkprocs; sk++) {
        if (tri.KProcDepSpecTri(sk, &tri, spec_gidx)) {
            updset.insert(tri.getKProc(sk));
        }
    }
    _updateLocal(updset);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getCompReacH(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = comp->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx];
    AssertLog(lcomp->def() == comp);

    if (lcomp->tets().empty()) {
        return 0.0;
    }

    double local_h = 0.0;
    for (auto t: lcomp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        Reac& reac = t->reac(lridx);
        local_h += reac.h();
    }
    double global_h = 0.0;
    MPI_Allreduce(&local_h, &global_h, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_h;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getCompReacC(solver::comp_global_id cidx, solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = comp->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx];
    AssertLog(lcomp->def() == comp);

    if (lcomp->tets().empty()) {
        return 0.0;
    }
    double local_c = 0.0;
    double local_v = 0.0;
    for (auto t: lcomp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        double v = t->vol();
        Reac& reac = t->reac(lridx);
        local_c += reac.c() * v;
        local_v += v;
    }
    double global_c = 0.0;
    double global_v = 0.0;
    MPI_Allreduce(&local_c, &global_c, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_v, &global_v, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_c / global_v;
}

////////////////////////////////////////////////////////////////////////////////

long double TetOpSplitP::_getCompReacA(solver::comp_global_id cidx,
                                       solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = comp->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx];
    AssertLog(lcomp->def() == comp);

    if (lcomp->tets().empty()) {
        return 0.0L;
    }

    long double local_a = 0.0L;
    for (auto t: lcomp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        Reac& reac = t->reac(lridx);
        local_a += static_cast<long double>(reac.rate());
    }
    long double global_a = 0.0L;
    MPI_Allreduce(&local_a, &global_a, 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_a;
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetOpSplitP::_getCompReacExtent(solver::comp_global_id cidx,
                                                   solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = comp->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx];
    AssertLog(lcomp->def() == comp);

    if (lcomp->tets().empty()) {
        return 0;
    }

    unsigned long long local_x = 0;
    for (auto t: lcomp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        Reac& reac = t->reac(lridx);
        local_x += reac.getExtent();
    }

    unsigned long long global_x = 0;
    MPI_Allreduce(&local_x, &global_x, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return global_x;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_resetCompReacExtent(solver::comp_global_id cidx, solver::reac_global_id ridx) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(ridx < statedef().countReacs());
    solver::Compdef* comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = comp->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    Comp* lcomp = pComps[cidx];
    AssertLog(lcomp->def() == comp);

    if (lcomp->tets().empty()) {
        return;
    }

    for (auto t: lcomp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        t->reac(lridx).resetExtent();
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getPatchSReacH(solver::patch_global_id pidx,
                                    solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    solver::sreac_local_id lsridx = patch->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx];
    AssertLog(lpatch->def() == patch);

    if (lpatch->tris().empty()) {
        return 0.0;
    }

    double local_h = 0.0;
    for (auto t: lpatch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        SReac& sreac = t->sreac(lsridx);
        local_h += sreac.h();
    }

    double global_h = 0.0;
    MPI_Allreduce(&local_h, &global_h, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_h;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getPatchSReacC(solver::patch_global_id pidx,
                                    solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    solver::sreac_local_id lsridx = patch->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx];
    AssertLog(lpatch->def() == patch);

    if (lpatch->tris().empty()) {
        return 0.0;
    }

    double local_c = 0.0;
    double local_a = 0.0;
    for (auto t: lpatch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        double a = t->area();
        SReac& sreac = t->sreac(lsridx);
        local_c += sreac.c() * a;
        local_a += a;
    }
    double global_c = 0.0;
    double global_a = 0.0;
    MPI_Allreduce(&local_c, &global_c, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_a, &global_a, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_c / global_a;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getPatchSReacA(solver::patch_global_id pidx,
                                    solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    solver::sreac_local_id lsridx = patch->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx];
    AssertLog(lpatch->def() == patch);

    if (lpatch->tris().empty()) {
        return 0.0;
    }

    double local_a = 0.0;
    for (auto t: lpatch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        SReac& sreac = t->sreac(lsridx);
        local_a += sreac.rate();
    }
    double global_a = 0.0;
    MPI_Allreduce(&local_a, &global_a, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_a;
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetOpSplitP::_getPatchSReacExtent(solver::patch_global_id pidx,
                                                     solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    solver::sreac_local_id lsridx = patch->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx];
    AssertLog(lpatch->def() == patch);

    if (lpatch->tris().empty()) {
        return 0;
    }

    unsigned long long local_x = 0;
    for (auto t: lpatch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        SReac& sreac = t->sreac(lsridx);
        local_x += sreac.getExtent();
    }
    unsigned long long global_x = 0.0;
    MPI_Allreduce(&local_x, &global_x, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return global_x;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_resetPatchSReacExtent(solver::patch_global_id pidx,
                                         solver::sreac_global_id ridx) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(ridx < statedef().countSReacs());
    solver::Patchdef* patch = statedef().patchdef(pidx);
    solver::sreac_local_id lsridx = patch->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        ArgErrLog(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    Patch* lpatch = pPatches[pidx];
    AssertLog(lpatch->def() == patch);

    if (lpatch->tris().empty()) {
        return;
    }

    for (auto t: lpatch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        t->sreac(lsridx).resetExtent();
    }
}

////////////////////////////////////////////////////////////////////////////////

const Tet& TetOpSplitP::_getTet(tetrahedron_global_id tidx) const {
    auto tet = pTets.at(tidx);
    if (tet == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
        ArgErrLog(os.str());
    }
    return *tet;
}

////////////////////////////////////////////////////////////////////////////////

Tet& TetOpSplitP::_getTet(tetrahedron_global_id tidx) {
    auto tet = pTets.at(tidx);
    if (tet == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
        ArgErrLog(os.str());
    }
    return *tet;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTetVol(tetrahedron_global_id tidx) const {
    return _getTet(tidx).vol();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTetVol(tetrahedron_global_id /*tidx*/, double /*vol*/) {
    std::ostringstream os;
    os << "Can not change tetrahedron volume in a mesh based solver.\n";
    NotImplErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getTetSpecDefined(tetrahedron_global_id tidx,
                                     solver::spec_global_id sidx) const {
    AssertLog(sidx < statedef().countSpecs());
    auto tet = pTets.at(tidx);

    if (tet == nullptr) {
        return false;
    }

    solver::spec_local_id lsidx = tet->compdef()->specG2L(sidx);
    return lsidx.valid();
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTetSpecCount(tetrahedron_global_id tidx,
                                     solver::spec_global_id sidx) const {
    MPI_Barrier(MPI_COMM_WORLD);
    AssertLog(sidx < statedef().countSpecs());
    auto& tet = _getTet(tidx);
    solver::spec_local_id lsidx = tet.compdef()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    uint count = tet.pools()[lsidx];
    MPI_Bcast(&count, 1, MPI_UNSIGNED, tetHosts.at(tidx), MPI_COMM_WORLD);
    return count;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTetSpecCount(tetrahedron_global_id tidx,
                                   solver::spec_global_id sidx,
                                   double n) {
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(n >= 0.0);
    auto& tet = _getTet(tidx);

    if (n > UINT_MAX) {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << UINT_MAX << ").\n";
        ArgErrLog(os.str());
    }

    if (tet.getInHost()) {
        solver::spec_local_id lsidx = tet.compdef()->specG2L(sidx);
        if (lsidx.unknown()) {
            std::ostringstream os;
            os << "Species undefined in tetrahedron.\n";
            ArgErrLog(os.str());
        }

        uint count = 0;

        double n_int = std::floor(n);
        double n_frc = n - n_int;
        count = static_cast<uint>(n_int);
        if (n_frc > 0.0) {
            double rand01 = rng()->getUnfIE();
            if (rand01 < n_frc) {
                count++;
            }
        }

        // don't need sync
        tet.setCount(lsidx, count);
        _updateSpec(tet, sidx);
        _updateSum();
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTetSpecAmount(tetrahedron_global_id tidx,
                                      solver::spec_global_id sidx) const {
    // following method does all necessary argument checking
    double count = _getTetSpecCount(tidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTetSpecAmount(tetrahedron_global_id tidx,
                                    solver::spec_global_id sidx,
                                    double m) {
    // convert amount in mols to number of molecules
    double m2 = m * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setTetSpecCount(tidx, sidx, m2);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTetSpecConc(tetrahedron_global_id tidx, solver::spec_global_id sidx) const {
    // following method does all necessary argument checking
    double count = _getTetSpecCount(tidx, sidx);
    Tet* tet = pTets[tidx];
    double vol = tet->vol();
    return (count / (1.0e3 * vol * math::AVOGADRO));
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTetSpecConc(tetrahedron_global_id tidx,
                                  solver::spec_global_id sidx,
                                  double c) {
    AssertLog(c >= 0.0);
    auto& tet = _getTet(tidx);
    double count = c * (1.0e3 * tet.vol() * math::AVOGADRO);
    // the following method does all the necessary argument checking
    _setTetSpecCount(tidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getTetSpecClamped(tetrahedron_global_id tidx,
                                     solver::spec_global_id sidx) const {
    AssertLog(sidx < statedef().countSpecs());
    auto& tet = _getTet(tidx);

    solver::spec_local_id lsidx = tet.compdef()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    return tet.clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTetSpecClamped(tetrahedron_global_id tidx,
                                     solver::spec_global_id sidx,
                                     bool buf) {
    AssertLog(sidx < statedef().countSpecs());
    auto tet = pTets.at(tidx);

    if (tet == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }

    solver::spec_local_id lsidx = tet->compdef()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    tet->setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTetReacK(tetrahedron_global_id tidx, solver::reac_global_id ridx) const {
    AssertLog(ridx < statedef().countReacs());
    auto& tet = _getTet(tidx);
    auto host = _getHost(tidx);

    solver::reac_local_id lridx = tet.compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }
    double kcst = 0;
    if (tet.getInHost()) {
        kcst = tet.reac(lridx).kcst();
    }
    MPI_Bcast(&kcst, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return kcst;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTetReacK(tetrahedron_global_id tidx, solver::reac_global_id ridx, double kf) {
    AssertLog(ridx < statedef().countReacs());
    AssertLog(kf >= 0.0);
    auto& tet = _getTet(tidx);

    solver::reac_local_id lridx = tet.compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "\nReaction undefined in tetrahedron.";
        ArgErrLog(os.str());
    }

    if (!tet.getInHost()) {
        return;
    }
    tet.reac(lridx).setKcst(kf);
    _updateElement(&tet.reac(lridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getTetReacActive(tetrahedron_global_id tidx, solver::reac_global_id ridx) const {
    AssertLog(ridx < statedef().countReacs());

    auto& tet = _getTet(tidx);
    auto host = _getHost(tidx);

    solver::reac_local_id lridx = tet.compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    bool active = false;
    if (tet.getInHost()) {
        if (tet.reac(lridx).inactive()) {
            active = false;
        } else {
            active = true;
        }
    }
    MPI_Bcast(&active, 1, MPI_C_BOOL, host, MPI_COMM_WORLD);
    return active;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTetReacActive(tetrahedron_global_id tidx,
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
    if (!tet.getInHost()) {
        return;
    }
    tet.reac(lridx).setActive(act);
    _updateElement(&tet.reac(lridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTetDiffD(tetrahedron_global_id tidx,
                                 solver::diff_global_id didx,
                                 tetrahedron_global_id direction_tet) const {
    AssertLog(didx < statedef().countDiffs());

    auto& tet = _getTet(tidx);
    auto host = _getHost(tidx);

    solver::diff_local_id ldidx = tet.compdef()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }
    double dcst = 0.0;
    if (tet.getInHost()) {
        if (direction_tet.unknown()) {
            dcst = tet.diff(ldidx).dcst();
        } else {
            int direction = tet.getTetDirection(direction_tet);
            if (direction == -1) {
                std::ostringstream os;
                os << "Tetrahedron " << direction_tet << " is not a neighbor of tetrahedron "
                   << tidx << ".\n";
                ArgErrLog(os.str());
            }

            dcst = tet.diff(ldidx).dcst(direction);
        }
    }
    MPI_Bcast(&dcst, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return dcst;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTetDiffD(tetrahedron_global_id tidx,
                               solver::diff_global_id didx,
                               double dk,
                               tetrahedron_global_id direction_tet) {
    AssertLog(didx < statedef().countDiffs());
    auto& tet = _getTet(tidx);
    recomputeUpdPeriod = true;

    solver::diff_local_id ldidx = tet.compdef()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    if (!tet.getInHost()) {
        return;
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
    _updateElement(&tet.diff(ldidx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getTetDiffActive(tetrahedron_global_id tidx, solver::diff_global_id didx) const {
    AssertLog(didx < statedef().countDiffs());
    auto& tet = _getTet(tidx);
    auto host = _getHost(tidx);

    solver::diff_local_id ldidx = tet.compdef()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }
    bool active = false;
    if (tet.getInHost()) {
        if (tet.diff(ldidx).inactive()) {
            active = false;
        } else {
            active = true;
        }
    }
    MPI_Bcast(&active, 1, MPI_C_BOOL, host, MPI_COMM_WORLD);
    return active;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTetDiffActive(tetrahedron_global_id tidx,
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
    if (!tet.getInHost()) {
        return;
    }
    tet.diff(ldidx).setActive(act);

    recomputeUpdPeriod = true;

    _updateElement(&tet.diff(ldidx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTetReacH(tetrahedron_global_id tidx, solver::reac_global_id ridx) const {
    AssertLog(ridx < statedef().countReacs());
    auto& tet = _getTet(tidx);
    auto host = _getHost(tidx);

    solver::reac_local_id lridx = tet.compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }
    double h = 0;
    if (tet.getInHost()) {
        h = tet.reac(lridx).h();
    }
    MPI_Bcast(&h, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return h;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTetReacC(tetrahedron_global_id tidx, solver::reac_global_id ridx) const {
    AssertLog(ridx < statedef().countReacs());
    auto& tet = _getTet(tidx);
    auto host = _getHost(tidx);

    solver::reac_local_id lridx = tet.compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    double c = 0;
    if (tet.getInHost()) {
        c = tet.reac(lridx).c();
    }
    MPI_Bcast(&c, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return c;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTetReacA(tetrahedron_global_id tidx, solver::reac_global_id ridx) const {
    AssertLog(ridx < statedef().countReacs());
    auto& tet = _getTet(tidx);
    auto host = _getHost(tidx);

    solver::reac_local_id lridx = tet.compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }
    double a = 0;
    if (tet.getInHost()) {
        a = tet.reac(lridx).rate();
    }
    MPI_Bcast(&a, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return a;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTetDiffA(tetrahedron_global_id tidx, solver::diff_global_id didx) const {
    AssertLog(didx < statedef().countDiffs());
    auto& tet = _getTet(tidx);
    auto host = _getHost(tidx);

    solver::diff_local_id ldidx = tet.compdef()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }
    double a = 0;
    if (tet.getInHost()) {
        a = tet.diff(ldidx).rate();
    }
    MPI_Bcast(&a, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return a;
}

////////////////////////////////////////////////////////////////////////////////

const Tri& TetOpSplitP::_getTri(triangle_global_id tidx) const {
    auto tri = pTris.at(tidx);
    if (tri == nullptr) {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.";
        ArgErrLog(os.str());
    }
    return *tri;
}

////////////////////////////////////////////////////////////////////////////////

Tri& TetOpSplitP::_getTri(triangle_global_id tidx) {
    auto tri = pTris.at(tidx);
    if (tri == nullptr) {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.";
        ArgErrLog(os.str());
    }
    return *tri;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriArea(triangle_global_id tidx) const {
    return _getTri(tidx).area();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTriArea(triangle_global_id /*tidx*/, double /*area*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getTriSpecDefined(triangle_global_id tidx, solver::spec_global_id sidx) const {
    AssertLog(sidx < statedef().countSpecs());
    auto tri = pTris.at(tidx);

    if (tri == nullptr) {
        return false;
    }

    solver::spec_local_id lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx.unknown()) {
        return false;
    } else {
        return true;
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriSpecCount(triangle_global_id tidx, solver::spec_global_id sidx) const {
    MPI_Barrier(MPI_COMM_WORLD);
    AssertLog(sidx < statedef().countSpecs());
    auto& tri = _getTri(tidx);
    auto host = _getHost(tidx);

    solver::spec_local_id lsidx = tri.patchdef()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    uint count = tri.pools()[lsidx];

    MPI_Bcast(&count, 1, MPI_UNSIGNED, host, MPI_COMM_WORLD);
    return count;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTriSpecCount(triangle_global_id tidx, solver::spec_global_id sidx, double n) {
    MPI_Barrier(MPI_COMM_WORLD);
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(n >= 0.0);
    auto& tri = _getTri(tidx);

    if (n > UINT_MAX) {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << UINT_MAX << ").\n";
        ArgErrLog(os.str());
    }

    if (tri.getInHost()) {
        solver::spec_local_id lsidx = tri.patchdef()->specG2L(sidx);
        if (lsidx.unknown()) {
            std::ostringstream os;
            os << "Species undefined in triangle.\n";
            ArgErrLog(os.str());
        }

        uint count = 0;

        double n_int = std::floor(n);
        double n_frc = n - n_int;
        count = static_cast<uint>(n_int);
        if (n_frc > 0.0) {
            double rand01 = rng()->getUnfIE();
            if (rand01 < n_frc) {
                count++;
            }
        }

        tri.setCount(lsidx, count);
        _updateSpec(tri, sidx);
        _updateSum();
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriSpecAmount(triangle_global_id tidx, solver::spec_global_id sidx) const {
    // following method does all necessary argument checking
    double count = _getTriSpecCount(tidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTriSpecAmount(triangle_global_id tidx,
                                    solver::spec_global_id sidx,
                                    double m) {
    // convert amount in mols to number of molecules
    double m2 = m * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setTriSpecCount(tidx, sidx, m2);
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getTriSpecClamped(triangle_global_id tidx, solver::spec_global_id sidx) const {
    AssertLog(sidx < statedef().countSpecs());

    auto& tri = _getTri(tidx);
    solver::spec_local_id lsidx = tri.patchdef()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    return tri.clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTriSpecClamped(triangle_global_id tidx,
                                     solver::spec_global_id sidx,
                                     bool buf) {
    AssertLog(sidx < statedef().countSpecs());
    auto& tri = _getTri(tidx);
    solver::spec_local_id lsidx = tri.patchdef()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    tri.setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriSReacK(triangle_global_id tidx, solver::sreac_global_id ridx) const {
    AssertLog(ridx < statedef().countSReacs());
    auto& tri = _getTri(tidx);
    auto host = _getHost(tidx);
    solver::sreac_local_id lsridx = tri.patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }
    double kcst = 0;
    if (tri.getInHost()) {
        kcst = tri.sreac(lsridx).kcst();
    }
    MPI_Bcast(&kcst, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return kcst;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTriSReacK(triangle_global_id tidx, solver::sreac_global_id ridx, double kf) {
    AssertLog(ridx < statedef().countSReacs());
    auto& tri = _getTri(tidx);

    solver::sreac_local_id lsridx = tri.patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }
    if (!tri.getInHost()) {
        return;
    }

    tri.sreac(lsridx).setKcst(kf);
    _updateElement(&tri.sreac(lsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getTriSReacActive(triangle_global_id tidx, solver::sreac_global_id ridx) const {
    AssertLog(ridx < statedef().countSReacs());

    auto& tri = _getTri(tidx);
    auto host = _getHost(tidx);

    solver::sreac_local_id lsridx = tri.patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }
    bool active = false;
    if (tri.getInHost()) {
        if (tri.sreac(lsridx).inactive()) {
            active = false;
        } else {
            active = true;
        }
    }
    MPI_Bcast(&active, 1, MPI_C_BOOL, host, MPI_COMM_WORLD);
    return active;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTriSReacActive(triangle_global_id tidx,
                                     solver::sreac_global_id ridx,
                                     bool act) {
    AssertLog(ridx < statedef().countSReacs());

    auto& tri = _getTri(tidx);

    solver::sreac_local_id lsridx = tri.patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }
    if (!tri.getInHost()) {
        return;
    }
    tri.sreac(lsridx).setActive(act);
    _updateElement(&tri.sreac(lsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriSDiffD(triangle_global_id tidx,
                                  solver::surfdiff_global_id didx,
                                  triangle_global_id direction_tri) const {
    AssertLog(didx < statedef().countSurfDiffs());
    auto& tri = _getTri(tidx);
    auto host = _getHost(tidx);

    solver::surfdiff_local_id ldidx = tri.patchdef()->surfdiffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in triangle.\n";
        ArgErrLog(os.str());
    }
    double dcst = 0.0;
    if (tri.getInHost()) {
        if (direction_tri.unknown()) {
            dcst = tri.sdiff(ldidx).dcst();

        } else {
            int direction = tri.getTriDirection(direction_tri);
            if (direction == -1) {
                std::ostringstream os;
                os << "Triangle " << direction_tri << " is not a neighbor of triangle " << tidx
                   << ".\n";
                ArgErrLog(os.str());
            }

            dcst = tri.sdiff(ldidx).dcst(direction);
        }
    }
    MPI_Bcast(&dcst, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return dcst;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTriSDiffD(triangle_global_id tidx,
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
    recomputeUpdPeriod = true;
    if (!tri.getInHost()) {
        return;
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
    _updateElement(&tri.sdiff(ldidx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getTriVDepSReacActive(triangle_global_id tidx,
                                         solver::vdepsreac_global_id vsridx) const {
    AssertLog(vsridx < statedef().countVDepSReacs());

    auto& tri = _getTri(tidx);
    auto host = _getHost(tidx);

    solver::vdepsreac_local_id lvsridx = tri.patchdef()->vdepsreacG2L(vsridx);
    if (lvsridx.unknown()) {
        std::ostringstream os;
        os << "Voltage-dependent surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }
    bool active = false;
    if (tri.getInHost()) {
        if (tri.vdepsreac(lvsridx).inactive()) {
            active = false;
        } else {
            active = true;
        }
    }
    MPI_Bcast(&active, 1, MPI_C_BOOL, host, MPI_COMM_WORLD);
    return active;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTriVDepSReacActive(triangle_global_id tidx,
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
    if (!tri.getInHost()) {
        return;
    }
    tri.vdepsreac(lvsridx).setActive(act);
    _updateElement(&tri.vdepsreac(lvsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriSReacH(triangle_global_id tidx, solver::sreac_global_id ridx) const {
    AssertLog(ridx < statedef().countSReacs());

    auto& tri = _getTri(tidx);
    auto host = _getHost(tidx);

    solver::sreac_local_id lsridx = tri.patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }
    double h = 0;
    if (tri.getInHost()) {
        h = tri.sreac(lsridx).h();
    }
    MPI_Bcast(&h, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return h;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriSReacC(triangle_global_id tidx, solver::sreac_global_id ridx) const {
    AssertLog(ridx < statedef().countSReacs());

    auto& tri = _getTri(tidx);
    auto host = _getHost(tidx);

    solver::sreac_local_id lsridx = tri.patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    double c = 0;
    if (tri.getInHost()) {
        c = tri.sreac(lsridx).c();
    }
    MPI_Bcast(&c, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return c;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriSReacA(triangle_global_id tidx, solver::sreac_global_id ridx) const {
    AssertLog(ridx < statedef().countSReacs());

    auto& tri = _getTri(tidx);
    auto host = _getHost(tidx);

    solver::sreac_local_id lsridx = tri.patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    double a = 0;
    if (tri.getInHost()) {
        a = tri.sreac(lsridx).rate();
    }
    MPI_Bcast(&a, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return a;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setEfieldDT(double efdt) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
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

double TetOpSplitP::_getTetV(tetrahedron_global_id tidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    auto loctidx = pEFTet_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    return pEField->getTetV(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTetV(tetrahedron_global_id tidx, double v) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    auto loctidx = pEFTet_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    // EField object should convert to millivolts
    pEField->setTetV(loctidx, v);

    // separate structure to store the EField triangle voltage, may need
    // refreshing.
    _refreshEFTrisV();

    // Voltage-dependent reactions may have changed
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getTetVClamped(tetrahedron_global_id tidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    auto loctidx = pEFTet_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    return pEField->getTetVClamped(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTetVClamped(tetrahedron_global_id tidx, bool cl) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    auto loctidx = pEFTet_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    pEField->setTetVClamped(loctidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriV(triangle_global_id tidx) const {
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
    // EField object should convert value to base s.i. units
    return EFTrisV[loctidx];
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTriV(triangle_global_id tidx, double v) {
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

    EFTrisV[loctidx] = v;
    // EField object should convert to millivolts
    pEField->setTriV(loctidx, v);

    // Voltage-dependent reactions may have changed
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getTriVClamped(triangle_global_id tidx) const {
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

    return pEField->getTriVClamped(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTriVClamped(triangle_global_id tidx, bool cl) {
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

    pEField->setTriVClamped(loctidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTriOhmicErev(triangle_global_id tidx,
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

    if (tri->getInHost()) {
        solver::ohmiccurr_local_id locidx = tri->patchdef()->ohmiccurrG2L(ocgidx);
        if (locidx.unknown()) {
            ArgErrLog("Ohmic current undefined in triangle.\n");
        }
        tri->setOCerev(locidx, erev);
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriOhmicErev(triangle_global_id tidx,
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

    auto& tri = _getTri(tidx);
    auto host = _getHost(tidx);

    double erev = 0.0;
    if (tri.getInHost()) {
        solver::ohmiccurr_local_id locidx = tri.patchdef()->ohmiccurrG2L(ocgidx);
        if (locidx.unknown()) {
            ArgErrLog("Ohmic current undefined in triangle.\n");
        }
        erev = tri.getOCerev(locidx);
    }

    MPI_Bcast(&erev, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return erev;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriOhmicI(triangle_global_id tidx) const {
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
    auto& tri = _getTri(tidx);
    auto host = _getHost(tidx);
    double cur = 0.0;
    if (tri.getInHost()) {
        cur = tri.getOhmicI(EFTrisV[loctidx], efdt());
    }
    MPI_Bcast(&cur, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return cur;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriOhmicI(triangle_global_id tidx,
                                  solver::ohmiccurr_global_id ocidx) const {
    AssertLog(ocidx < statedef().countOhmicCurrs());
    auto& tri = _getTri(tidx);
    auto host = _getHost(tidx);

    auto loctidx = pEFTri_GtoL[tidx];
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
    double cur = 0.0;
    if (tri.getInHost()) {
        cur = tri.getOhmicI(locidx, EFTrisV[loctidx], efdt());
    }
    MPI_Bcast(&cur, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return cur;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriGHKI(triangle_global_id tidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }

    auto& tri = _getTri(tidx);
    auto host = _getHost(tidx);

    double cur = 0.0;
    if (tri.getInHost()) {
        cur = tri.getGHKI();
    }
    MPI_Bcast(&cur, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return cur;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriGHKI(triangle_global_id tidx, solver::ghkcurr_global_id ghkidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }

    auto& tri = _getTri(tidx);
    auto host = _getHost(tidx);

    solver::ghkcurr_local_id locidx = tri.patchdef()->ghkcurrG2L(ghkidx);
    if (locidx.unknown()) {
        std::ostringstream os;
        os << "GHK current undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    double cur = 0.0;
    if (tri.getInHost()) {
        cur = tri.getGHKI(locidx);
    }
    MPI_Bcast(&cur, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return cur;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriI(triangle_global_id tidx) const {
    return _getTriGHKI(tidx) + _getTriOhmicI(tidx);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getVertIClamp(vertex_id_t vidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    auto locvidx = pEFVert_GtoL[vidx];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }

    // EField object should convert to required units
    return pEField->getVertIClamp(locvidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setVertIClamp(vertex_id_t vidx, double cur) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    auto locvidx = pEFVert_GtoL[vidx];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }

    // EField object should convert to required units
    pEField->setVertIClamp(locvidx, cur);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getTriIClamp(triangle_global_id tidx) const {
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

    // EField object should convert to required units
    return pEField->getTriIClamp(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTriIClamp(triangle_global_id tidx, double cur) {
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

    // EField object should convert to required units
    pEField->setTriIClamp(loctidx, cur);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setTriCapac(triangle_global_id tidx, double cap) {
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

    // EField object should convert to required units
    pEField->setTriCapac(loctidx, cap);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::_getVertV(vertex_id_t vidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    auto locvidx = pEFVert_GtoL[vidx];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }
    return pEField->getVertV(locvidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setVertV(vertex_id_t vidx, double v) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    auto locvidx = pEFVert_GtoL[vidx];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }
    // EField object should convert to millivolts
    pEField->setVertV(locvidx, v);

    // separate structure to store the EField triangle voltage, needs refreshing.
    _refreshEFTrisV();
    // Voltage-dependent reactions may have changed
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

bool TetOpSplitP::_getVertVClamped(vertex_id_t vidx) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    auto locvidx = pEFVert_GtoL[vidx];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }

    return pEField->getVertVClamped(locvidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setVertVClamped(vertex_id_t vidx, bool cl) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    auto locvidx = pEFVert_GtoL[vidx];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }

    // EField object should convert to millivolts
    pEField->setVertVClamped(locvidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setMembRes(solver::membrane_global_id midx, double ro, double vrev) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
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

std::pair<double, double> TetOpSplitP::_getMembRes(solver::membrane_global_id midx) const {
    if (!efflag()) {
        ArgErrLog("Method not available: EField calculation not included in simulation.");
    }
    // EField object should convert to required units
    AssertLog(midx.get() == 0);
    return pEField->getSurfaceResistivity();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setMembPotential(solver::membrane_global_id midx, double v) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    // EField object should convert to millivolts
    AssertLog(midx.get() == 0);
    pEField->setMembPotential(midx, v);

    // separate structure to store the EField triangle voltage, needs refreshing.
    _refreshEFTrisV();
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setMembCapac(solver::membrane_global_id midx, double cm) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    if (cm < 0.0) {
        std::ostringstream os;
        os << "Capacitance must be greater than or equal to zero.";
        ArgErrLog(os.str());
    }

    AssertLog(midx.get() == 0);
    // EField object should convert to required units
    pEField->setMembCapac(midx, cm);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_setMembVolRes(solver::membrane_global_id midx, double ro) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    if (ro < 0.0) {
        std::ostringstream os;
        os << "Resistivity must be greater than or equal to zero.";
        ArgErrLog(os.str());
    }
    AssertLog(midx.get() == 0);
    // EField object should convert to required units
    pEField->setMembVolRes(midx, ro);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_computeUpdPeriod() {
    double local_max_rate = 0.0;

    for (uint pos = 0; pos < diffSep; pos++) {
        Diff* d = pDiffs[pos];
        // Now ignoring inactive diffusion
        double scaleddcst = 0.0;
        if (d->active()) {
            scaleddcst = d->getScaledDcst();
        }
        if (scaleddcst > local_max_rate) {
            local_max_rate = scaleddcst;
        }
    }

    for (uint pos = 0; pos < sdiffSep; pos++) {
        SDiff* d = pSDiffs[pos];
        // Now ignoring inactive diffusion
        double scaleddcst = 0.0;
        if (d->active()) {
            scaleddcst = d->getScaledDcst();
        }
        if (scaleddcst > local_max_rate) {
            local_max_rate = scaleddcst;
        }
    }
    // get global max rate
    double global_max_rate = 0;
    MPI_Allreduce(&local_max_rate, &global_max_rate, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (global_max_rate < 0.0) {
        std::ostringstream os;
        os << "Maximum scaled diffusion constant is " << global_max_rate
           << ". This should not happen in this solver.\n";
        ArgErrLog(os.str());
    }
    updPeriod = 1.0 / global_max_rate;
    recomputeUpdPeriod = false;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_updateLocal(KProcPSet const& upd_entries) {
    for (auto& kp: upd_entries) {
        AssertLog(kp != nullptr);
        _updateElement(kp);
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_updateLocal(std::vector<KProc*> const& upd_entries) {
    for (auto& kp: upd_entries) {
        AssertLog(kp != nullptr);
        _updateElement(kp);
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_updateLocal(std::vector<uint> const& upd_entries) {
    for (auto& upd_idx: upd_entries) {
        KProc* kp = pKProcs[upd_idx];
        if (kp != nullptr) {
            _updateElement(kp);
        }
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_updateLocal() {
    for (uint i = 0; i < nEntries; i++) {
        KProc* kp = pKProcs[i];
        if (kp != nullptr) {
            _updateElement(kp);
        }
    }
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

CRGroup* TetOpSplitP::_getGroup(int pow) {
    if (pow >= 0) {
        return pGroups[pow];
    } else {
        return nGroups[-pow];
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_extendPGroups(uint new_size) {
    auto curr_size = pGroups.size();

    while (curr_size < new_size) {
        pGroups.push_back(new CRGroup(curr_size));
        curr_size++;
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_extendNGroups(uint new_size) {
    uint curr_size = nGroups.size();

    while (curr_size < new_size) {
        nGroups.push_back(new CRGroup(-curr_size));
        curr_size++;
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_extendGroup(CRGroup* group, uint size) {
    group->capacity += size;
    group->indices = static_cast<KProc**>(
        realloc(group->indices, sizeof(KProc*) * group->capacity));
    if (group->indices == nullptr) {
        SysErrLog("DirectCR: unable to allocate memory for SSA group.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_updateSum() {
    pA0 = 0.0;

    auto n_neg_groups = nGroups.size();
    auto n_pos_groups = pGroups.size();

    for (uint i = 0; i < n_neg_groups; i++) {
        pA0 += nGroups[i]->sum;
    }

    for (uint i = 0; i < n_pos_groups; i++) {
        pA0 += pGroups[i]->sum;
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_updateElement(KProc* kp) {
    if (kp->getType() == KP_DIFF || kp->getType() == KP_SDIFF) {
        kp->crData.rate = kp->rate(this);
        return;
    }

    double new_rate = kp->rate(this);

    CRKProcData& data = kp->crData;
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
            if (pGroups.size() <= static_cast<unsigned>(new_pow)) {
                _extendPGroups(new_pow + 1);
            }

            CRGroup* new_group = pGroups[new_pow];

            AssertLog(new_group != nullptr);
            if (new_group->size == new_group->capacity) {
                _extendGroup(new_group);
            }
            uint pos = new_group->size;
            new_group->indices[pos] = kp;
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

            if (nGroups.size() <= static_cast<uint>(-new_pow)) {
                _extendNGroups(-new_pow + 1);
            }

            CRGroup* new_group = nGroups[-new_pow];

            if (new_group->size == new_group->capacity) {
                _extendGroup(new_group);
            }

            uint pos = new_group->size;
            new_group->indices[pos] = kp;
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
// Batch data recording
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetOpSplitP::getBatchTetSpecCounts(const std::vector<index_t>& tets,
                                                       std::string const& s) const {
    auto ntets = tets.size();
    std::vector<double> global_counts(ntets, 0.0);

    TetOpSplitP::getBatchTetSpecCountsNP(tets.data(), ntets, s, global_counts.data(), ntets);

    return global_counts;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetOpSplitP::getBatchTriSpecCounts(const std::vector<index_t>& tris,
                                                       std::string const& s) const {
    auto ntris = tris.size();
    std::vector<double> global_counts(ntris, 0.0);

    TetOpSplitP::getBatchTriSpecCountsNP(tris.data(), ntris, s, global_counts.data(), ntris);

    return global_counts;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setBatchTetSpecConcs(const std::vector<index_t>& tets,
                                       std::string const& s,
                                       const std::vector<double>& concs) {
    TetOpSplitP::setBatchTetSpecConcsNP(tets.data(), tets.size(), s, concs.data(), concs.size());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetOpSplitP::getBatchTetSpecConcs(const std::vector<index_t>& tets,
                                                      std::string const& s) const {
    auto ntets = tets.size();
    std::vector<double> global_concs(ntets, 0.0);

    TetOpSplitP::getBatchTetSpecConcsNP(tets.data(), ntets, s, global_concs.data(), ntets);

    return global_concs;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::getBatchTetSpecCountsNP(const index_t* indices,
                                          size_t input_size,
                                          std::string const& s,
                                          double* counts,
                                          size_t output_size) const {
    if (input_size != output_size) {
        std::ostringstream os;
        os << "Error: output array (counts) size should be the same as input array "
              "(indices) size.\n";
        ArgErrLog(os.str());
    }
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;
    std::vector<double> local_counts(input_size, 0.0);

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (uint t = 0; t < input_size; t++) {
        tetrahedron_global_id tidx(indices[t]);

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

        if (tet->getInHost()) {
            local_counts[t] = tet->pools()[slidx];
        }
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "Species " << s
            << " has not been defined in the following tetrahedrons, fill in zeros "
               "at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
    MPI_Allreduce(local_counts.data(), counts, input_size, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::getBatchTriSpecCountsNP(const index_t* indices,
                                          size_t input_size,
                                          std::string const& s,
                                          double* counts,
                                          size_t output_size) const {
    if (input_size != output_size) {
        std::ostringstream os;
        os << "Error: output array (counts) size should be the same as input array "
              "(indices) size.\n";
        ArgErrLog(os.str());
    }
    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);
    std::vector<double> local_counts(input_size, 0.0);
    for (uint t = 0; t < input_size; t++) {
        triangle_global_id tidx(indices[t]);

        if (tidx >= pTris.size()) {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            ArgErrLog(os.str());
        }

        if (pTris[tidx] == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        Tri* tri = pTris[tidx];
        solver::spec_local_id slidx = tri->patchdef()->specG2L(sgidx);
        if (slidx.unknown()) {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }
        if (tri->getInHost()) {
            local_counts[t] = tri->pools()[slidx];
        }
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a patch, fill in "
               "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "Species " << s
            << " has not been defined in the following triangles, fill in zeros at "
               "target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
    MPI_Allreduce(local_counts.data(), counts, input_size, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setBatchTetSpecConcsNP(const index_t* indices,
                                         size_t ntets,
                                         std::string const& s,
                                         const double* concs,
                                         size_t output_size) {
    if (ntets != output_size) {
        std::ostringstream os;
        os << "Error: output array (concs) size should be the same as input array "
              "(indices) size.\n";
        ArgErrLog(os.str());
    }
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (size_t t = 0; t < ntets; t++) {
        tetrahedron_global_id tidx(indices[t]);

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
        if (tet->getInHost()) {
            _setTetSpecConc(tidx, sgidx, concs[t]);
        }
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "Species " << s
            << " has not been defined in the following tetrahedrons, fill in zeros "
               "at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::getBatchTetSpecConcsNP(const index_t* indices,
                                         size_t ntets,
                                         std::string const& s,
                                         double* global_concs,
                                         size_t output_size) const {
    if (ntets != output_size) {
        std::ostringstream os;
        os << "Error: output array (concs) size should be the same as input array "
              "(indices) size.\n";
        ArgErrLog(os.str());
    }
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);
    std::vector<double> local_concs(ntets, 0.0);

    for (uint t = 0; t < ntets; t++) {
        tetrahedron_global_id tidx(indices[t]);

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
        if (tet->getInHost()) {
            double count = tet->pools()[slidx];
            double vol = tet->vol();
            local_concs[t] = (count / (1.0e3 * vol * math::AVOGADRO));
        }
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "Species " << s
            << " has not been defined in the following tetrahedrons, fill in zeros "
               "at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
    MPI_Allreduce(local_concs.data(), global_concs, ntets, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::sumBatchTetCountsNP(const index_t* indices,
                                        size_t input_size,
                                        std::string const& s) {
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);
    double partial_sum = 0.0;
    for (uint t = 0; t < input_size; t++) {
        tetrahedron_global_id tidx(indices[t]);

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

        if (tet->getInHost()) {
            partial_sum += tet->pools()[slidx];
        }
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "Species " << s
            << " has not been defined in the following tetrahedrons, fill in zeros "
               "at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
    double global_sum = 0.0;
    MPI_Allreduce(&partial_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_sum;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::sumBatchTriCountsNP(const index_t* indices,
                                        size_t input_size,
                                        std::string const& s) {
    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;

    double partial_sum = 0.0;
    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (uint t = 0; t < input_size; t++) {
        triangle_global_id tidx(indices[t]);

        if (tidx >= pTris.size()) {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            ArgErrLog(os.str());
        }

        if (pTris[tidx] == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        Tri* tri = pTris[tidx];
        solver::spec_local_id slidx = tri->patchdef()->specG2L(sgidx);
        if (slidx.unknown()) {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }
        if (tri->getInHost()) {
            partial_sum += tri->pools()[slidx];
        }
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a patch, fill in "
               "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "Species " << s
            << " has not been defined in the following triangles, fill in zeros at "
               "target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }

    double global_sum = 0.0;
    MPI_Allreduce(&partial_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_sum;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::sumBatchTriGHKIsNP(const index_t* indices,
                                       size_t input_size,
                                       std::string const& ghk) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    // the following may raise exception if string is unknown
    solver::ghkcurr_global_id ghkidx = statedef().getGHKcurrIdx(ghk);

    double partial_sum = 0.0;

    for (auto t = 0u; t < input_size; t++) {
        triangle_global_id tidx(indices[t]);

        if (tidx >= mesh()->countTris()) {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }

        Tri* tri = pTris[tidx];

        solver::ghkcurr_local_id locidx = tri->patchdef()->ghkcurrG2L(ghkidx);
        if (locidx.unknown()) {
            std::ostringstream os;
            os << "GHK current undefined in triangle.\n";
            ArgErrLog(os.str());
        }

        if (tri->getInHost()) {
            partial_sum += tri->getGHKI(locidx);
        }
    }
    double global_sum = 0.0;
    MPI_Allreduce(&partial_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_sum;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::sumBatchTriOhmicIsNP(const index_t* indices,
                                         size_t input_size,
                                         std::string const& oc) {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }

    // the following may raise exception if string is unknown
    solver::ohmiccurr_global_id ocidx = statedef().getOhmicCurrIdx(oc);

    double partial_sum = 0.0;

    for (uint t = 0; t < input_size; t++) {
        triangle_global_id tidx(indices[t]);
        auto loctidx = pEFTri_GtoL[tidx];
        if (tidx >= mesh()->countTris()) {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }

        Tri* tri = pTris[tidx];

        solver::ohmiccurr_local_id locidx = tri->patchdef()->ohmiccurrG2L(ocidx);
        if (locidx.unknown()) {
            std::ostringstream os;
            os << "Ohmic current undefined in triangle.\n";
            ArgErrLog(os.str());
        }

        if (tri->getInHost()) {
            partial_sum += tri->getOhmicI(locidx, EFTrisV[loctidx], efdt());
        }
    }
    double global_sum = 0.0;
    MPI_Allreduce(&partial_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_sum;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::getBatchTriOhmicIsNP(const index_t* indices,
                                       size_t input_size,
                                       std::string const& oc,
                                       double* counts,
                                       size_t output_size) const {
    if (input_size != output_size) {
        std::ostringstream os;
        os << "Error: output array (counts) size should be the same as input array "
              "(indices) size.\n";
        ArgErrLog(os.str());
    }

    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;

    solver::ohmiccurr_global_id ocidx = statedef().getOhmicCurrIdx(oc);
    std::vector<double> local_counts(input_size, 0.0);
    for (uint t = 0; t < input_size; t++) {
        triangle_global_id tidx(indices[t]);

        if (tidx >= pTris.size()) {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            ArgErrLog(os.str());
        }

        if (pTris[tidx] == nullptr) {
            tri_not_assign << tidx << ' ';
            has_tri_warning = true;
            continue;
        }

        Tri* tri = pTris[tidx];
        solver::ohmiccurr_local_id locidx = tri->patchdef()->ohmiccurrG2L(ocidx);
        if (locidx.unknown()) {
            spec_undefined << tidx << ' ';
            has_spec_warning = true;
            continue;
        }
        if (tri->getInHost()) {
            auto loctidx = pEFTri_GtoL[tidx];
            local_counts[t] = tri->getOhmicI(locidx, EFTrisV[loctidx], efdt());
        }
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a patch, fill in "
               "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << '\n';
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "Ohmic Current " << oc
            << " has not been defined in the following triangles, fill in zeros at "
               "target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << '\n';
    }
    MPI_Allreduce(local_counts.data(), counts, input_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::getBatchTriGHKIsNP(const index_t* indices,
                                     size_t input_size,
                                     std::string const& ghk,
                                     double* counts,
                                     size_t output_size) const {
    if (input_size != output_size) {
        std::ostringstream os;
        os << "Error: output array (counts) size should be the same as input array "
              "(indices) size.\n";
        ArgErrLog(os.str());
    }

    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;

    solver::ghkcurr_global_id ghkidx = statedef().getGHKcurrIdx(ghk);
    std::vector<double> local_counts(input_size, 0.0);
    for (uint t = 0; t < input_size; t++) {
        triangle_global_id tidx(indices[t]);

        if (tidx >= pTris.size()) {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            ArgErrLog(os.str());
        }

        if (pTris[tidx] == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        Tri* tri = pTris[tidx];
        solver::ghkcurr_local_id locidx = tri->patchdef()->ghkcurrG2L(ghkidx);
        if (locidx.unknown()) {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }
        if (tri->getInHost()) {
            local_counts[t] = tri->getGHKI(locidx);
        }
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a patch, fill in "
               "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "GHk Current " << ghk
            << " has not been defined in the following triangles, fill in zeros at "
               "target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
    MPI_Allreduce(local_counts.data(), counts, input_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::getBatchTriVsNP(const index_t* indices,
                                  size_t input_size,
                                  double* counts,
                                  size_t output_size) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }

    if (input_size != output_size) {
        std::ostringstream os;
        os << "Error: output array (counts) size should be the same as input array "
              "(indices) size.\n";
        ArgErrLog(os.str());
    }

    bool has_tri_warning = false;
    std::ostringstream tri_not_assign;

    for (uint t = 0; t < input_size; t++) {
        triangle_global_id tidx(indices[t]);
        if (tidx >= pTris.size()) {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            ArgErrLog(os.str());
        }
        auto loctidx = pEFTri_GtoL[tidx];
        if (loctidx.unknown()) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }
        counts[t] = EFTrisV[loctidx];
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a membrane, fill "
               "in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::getBatchTetVsNP(const index_t* indices,
                                  size_t input_size,
                                  double* counts,
                                  size_t output_size) const {
    if (!efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }

    if (input_size != output_size) {
        std::ostringstream os;
        os << "Error: output array (counts) size should be the same as input array "
              "(indices) size.\n";
        ArgErrLog(os.str());
    }

    bool has_tet_warning = false;
    std::ostringstream tet_not_assign;

    for (uint t = 0; t < input_size; t++) {
        tetrahedron_global_id tidx(indices[t]);

        if (tidx >= pTets.size()) {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            ArgErrLog(os.str());
        }

        auto loctidx = pEFTet_GtoL[tidx];

        if (loctidx.unknown()) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }
        counts[t] = pEField->getTetV(loctidx);
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log")
            << "The following tetrahedrons have not been assigned to a conduction "
               "volume, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::getBatchTriBatchOhmicIsNP(const index_t* indices,
                                            size_t input_size,
                                            std::vector<std::string>& ocs,
                                            double* counts,
                                            size_t output_size) const {
    auto n_ocs = ocs.size();
    if (input_size * n_ocs != output_size) {
        std::ostringstream os;
        os << "Error: output array (counts) size should be the same as input array "
              "(indices) size times the number of ocs. \n";
        ArgErrLog(os.str());
    }

    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;

    std::vector<solver::ohmiccurr_global_id> ocidxs;
    ocidxs.reserve(n_ocs);

    for (auto& oc: ocs) {
        ocidxs.push_back(statedef().getOhmicCurrIdx(oc));
    }

    std::vector<double> local_counts(output_size, 0.0);
    for (uint t = 0; t < input_size; t++) {
        triangle_global_id tidx(indices[t]);

        if (tidx >= pTris.size()) {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            ArgErrLog(os.str());
        }

        if (pTris[tidx] == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        Tri* tri = pTris[tidx];

        for (auto oc_counter = 0u; oc_counter < ocidxs.size(); oc_counter++) {
            solver::ohmiccurr_local_id locidx = tri->patchdef()->ohmiccurrG2L(ocidxs[oc_counter]);
            if (locidx.unknown()) {
                spec_undefined << tidx << ":" << ocs[oc_counter] << " ";
                has_spec_warning = true;
                continue;
            }
            if (tri->getInHost()) {
                auto loctidx = pEFTri_GtoL[tidx];
                local_counts[t * n_ocs + oc_counter] =
                    tri->getOhmicI(locidx, EFTrisV[loctidx], efdt());
            }
        }
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a patch, fill in "
               "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log") << "Ohmic Current:Triangle undefined, fill in "
                                        "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
    MPI_Allreduce(local_counts.data(), counts, output_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::getBatchTriBatchGHKIsNP(const index_t* indices,
                                          size_t input_size,
                                          std::vector<std::string>& ghks,
                                          double* counts,
                                          size_t output_size) const {
    auto n_ghks = ghks.size();
    if (input_size * n_ghks != output_size) {
        std::ostringstream os;
        os << "Error: output array (counts) size should be the same as input array "
              "(indices) size times the number of ghks. \n";
        ArgErrLog(os.str());
    }

    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;

    std::vector<solver::ghkcurr_global_id> ghkidxs;
    ghkidxs.reserve(n_ghks);

    for (auto& ghk: ghks) {
        ghkidxs.push_back(statedef().getGHKcurrIdx(ghk));
    }

    std::vector<double> local_counts(output_size, 0.0);

    for (uint t = 0; t < input_size; t++) {
        triangle_global_id tidx(indices[t]);

        if (tidx >= pTris.size()) {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            ArgErrLog(os.str());
        }

        if (pTris[tidx] == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        Tri* tri = pTris[tidx];
        for (auto ghk_counter = 0u; ghk_counter < ghkidxs.size(); ghk_counter++) {
            solver::ghkcurr_local_id locidx = tri->patchdef()->ghkcurrG2L(ghkidxs[ghk_counter]);
            if (locidx.unknown()) {
                spec_undefined << tidx << ":" << ghks[ghk_counter] << " ";
                has_spec_warning = true;
                continue;
            }
            if (tri->getInHost()) {
                local_counts[t * n_ghks + ghk_counter] = tri->getGHKI(locidx);
            }
        }
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a patch, fill in "
               "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log") << "GHK Current:Triangle undefined, fill in "
                                        "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
    MPI_Allreduce(local_counts.data(), counts, output_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////
// ROI Data Access
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetOpSplitP::getROITetSpecCounts(const std::string& ROI_id,
                                                     std::string const& s) const {
    auto const& roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    const auto size = roi->second.size();
    std::vector<double> data(size);
    getBatchTetSpecCountsNP(
        reinterpret_cast<const index_t*>(roi->second.data()), size, s, &data.front(), data.size());
    return data;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetOpSplitP::getROITriSpecCounts(const std::string& ROI_id,
                                                     std::string const& s) const {
    auto const& roi = mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    const auto size = roi->second.size();
    std::vector<double> data(size);
    getBatchTriSpecCountsNP(reinterpret_cast<const index_t*>(roi->second.data()),
                            roi->second.size(),
                            s,
                            &data.front(),
                            data.size());
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::getROITetSpecCountsNP(const std::string& ROI_id,
                                        std::string const& s,
                                        double* counts,
                                        size_t output_size) const {
    auto const& roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    getBatchTetSpecCountsNP(reinterpret_cast<const index_t*>(roi->second.data()),
                            roi->second.size(),
                            s,
                            counts,
                            output_size);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::getROITriSpecCountsNP(const std::string& ROI_id,
                                        std::string const& s,
                                        double* counts,
                                        size_t output_size) const {
    auto const& roi = mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    getBatchTriSpecCountsNP(reinterpret_cast<const index_t*>(roi->second.data()),
                            roi->second.size(),
                            s,
                            counts,
                            output_size);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::getROIVol(const std::string& ROI_id) const {
    auto const& roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    double sum = 0.0;
    for (auto const& tidx: roi->second) {
        sum += pTets[tidx]->vol();
    }
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::getROIArea(const std::string& ROI_id) const {
    auto const& roi = mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    double sum = 0.0;
    for (auto const& tidx: roi->second) {
        sum += pTris[tidx]->area();
    }
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::getROITetSpecCount(const std::vector<tetrahedron_global_id>& tetrahedrons,
                                       const std::string& s) const {
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;
    double local_sum = 0.0;
    double global_sum = 0.0;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: tetrahedrons) {
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

        // compute local sum for each process
        if (tet->getInHost()) {
            local_sum += tet->pools()[slidx];
        }
    }

    // gather global sum
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "Species " << s
            << " has not been defined in the following tetrahedrons, fill in zeros "
               "at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
    return global_sum;
}

double TetOpSplitP::getROITriSpecCount(const std::vector<triangle_global_id>& triangles,
                                       const std::string& s) const {
    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;
    double local_sum = 0.0;
    double global_sum = 0.0;

    const solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: triangles) {
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

        // compute local sum for each process
        if (tri->getInHost()) {
            local_sum += tri->pools()[slidx];
        }
    }

    // gather global sum
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (has_tri_warning) {
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a patch, fill in "
               "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "Species " << s
            << " has not been defined in the following triangles, fill in zeros at "
               "target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
    return global_sum;
}

double TetOpSplitP::getROISpecCount(const std::string& ROI_id, std::string const& s) const {
    {
        auto const& roi =
            mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id, 0 /* count */, false /* warning */);
        if (roi != mesh()->rois.end<tetmesh::ROI_TRI>()) {
            return getROITriSpecCount(roi->second, s);
        }
    }
    {
        auto const& roi =
            mesh()->rois.get<tetmesh::ROI_TET>(ROI_id, 0 /* count */, false /* warning */);
        if (roi != mesh()->rois.end<tetmesh::ROI_TET>()) {
            return getROITetSpecCount(roi->second, s);
        }
    }
    ArgErrLog("Error: Cannot find suitable ROI for the function call getROICount.\n");
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setROITriSpecCount(const std::vector<triangle_global_id>& triangles,
                                     const std::string& s,
                                     double count) {
    double totalarea = 0.0;
    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;
    std::vector<triangle_global_id> apply_indices;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: triangles) {
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

        apply_indices.push_back(tidx);
        totalarea += tri->area();
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a patch, fill in "
               "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "Species " << s
            << " has not been defined in the following triangles, fill in zeros at "
               "target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }

    auto ind_size = apply_indices.size();

    std::vector<double> apply_count(ind_size, 0.0);

    // Create distribution at node 0
    if (myRank == 0) {
        uint c = static_cast<uint>(count);
        uint nremoved = 0;

        for (uint t = 0; t < ind_size; t++) {
            auto tidx = apply_indices[t];
            Tri* tri = pTris[tidx];

            if ((count == 0.0) || (nremoved == c)) {
                break;
            }

            double fract = static_cast<double>(c) * (tri->area() / totalarea);
            uint n3 = static_cast<uint>(std::floor(fract));

            double n3_frac = fract - static_cast<double>(n3);
            if (n3_frac > 0.0) {
                double rand01 = rng()->getUnfIE();
                if (rand01 < n3_frac) {
                    n3++;
                }
            }

            nremoved += n3;

            if (nremoved >= c) {
                n3 -= (nremoved - c);
                nremoved = c;
            }

            apply_count[t] = n3;
        }
        AssertLog(nremoved <= c);
        c -= nremoved;
        while (c != 0) {
            double accum = 0.0;
            double selector = rng()->getUnfIE() * totalarea;
            for (uint t = 0; t < ind_size; t++) {
                auto tidx = apply_indices[t];
                Tri* tri = pTris[tidx];
                accum += tri->area();
                if (selector < accum) {
                    apply_count[t] += 1.0;
                    break;
                }
            }
            c--;
        }
    }

    MPI_Bcast(apply_count.data(), ind_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // counts need to be set globally so that sync can be avoided
    for (uint t = 0; t < ind_size; t++) {
        auto tidx = apply_indices[t];
        auto& tri = _getTri(tidx);

        solver::spec_local_id slidx = tri.patchdef()->specG2L(sgidx);
        tri.setCount(slidx, apply_count[t]);
        _updateSpec(tri, sgidx);
    }
    _updateSum();
}

void TetOpSplitP::setROITetSpecCount(const std::vector<tetrahedron_global_id>& tetrahedrons,
                                     const std::string& s,
                                     double count) {
    double totalvol = 0.0;
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;
    std::vector<tetrahedron_global_id> apply_indices;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: tetrahedrons) {
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

        apply_indices.push_back(tidx);
        totalvol += tet->vol();
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "Species " << s
            << " has not been defined in the following tetrahedrons, fill in zeros "
               "at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }

    auto ind_size = apply_indices.size();

    std::vector<double> apply_count(ind_size, 0.0);

    // Create distribution at node 0
    if (myRank == 0) {
        uint c = static_cast<uint>(count);
        uint nremoved = 0;

        for (uint t = 0; t < ind_size; t++) {
            auto tidx = apply_indices[t];
            Tet* tet = pTets[tidx];

            if ((count == 0.0) || (nremoved == c)) {
                break;
            }

            double fract = static_cast<double>(c) * (tet->vol() / totalvol);
            uint n3 = static_cast<uint>(std::floor(fract));

            double n3_frac = fract - static_cast<double>(n3);
            if (n3_frac > 0.0) {
                double rand01 = rng()->getUnfIE();
                if (rand01 < n3_frac) {
                    n3++;
                }
            }

            nremoved += n3;

            if (nremoved >= c) {
                n3 -= (nremoved - c);
                nremoved = c;
            }

            apply_count[t] = n3;
        }
        AssertLog(nremoved <= c);
        c -= nremoved;
        while (c != 0) {
            double accum = 0.0;
            double selector = rng()->getUnfIE() * totalvol;
            for (uint t = 0; t < ind_size; t++) {
                auto tidx = apply_indices[t];
                Tet* tet = pTets[tidx];
                accum += tet->vol();
                if (selector < accum) {
                    apply_count[t] += 1.0;
                    break;
                }
            }
            c--;
        }
    }

    MPI_Bcast(apply_count.data(), ind_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // set the counts golbally and update local KProcs
    for (uint t = 0; t < ind_size; t++) {
        auto tidx = apply_indices[t];
        auto& tet = _getTet(tidx);
        solver::spec_local_id slidx = tet.compdef()->specG2L(sgidx);
        tet.setCount(slidx, apply_count[t]);
        _updateSpec(tet, sgidx);
    }

    _updateSum();
}

// WEILIANG: Can we apply Sam's set count method (as in e.g. setCompCount) here?
void TetOpSplitP::setROISpecCount(const std::string& ROI_id, std::string const& s, double count) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (count > UINT_MAX) {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << UINT_MAX << ").\n";
        ArgErrLog(os.str());
    }

    auto const& tri_roi = mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id, 0, false);
    if (tri_roi != mesh()->rois.end<tetmesh::ROI_TRI>()) {
        setROITriSpecCount(tri_roi->second, s, count);
    } else {
        auto const& tet_roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id, 0, false);
        if (tet_roi != mesh()->rois.end<tetmesh::ROI_TET>()) {
            setROITetSpecCount(tet_roi->second, s, count);
        } else {
            std::ostringstream os;
            os << "Error: Cannot find suitable ROI for the function call "
                  "getROICount.\n";
            ArgErrLog(os.str());
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::getROISpecAmount(const std::string& ROI_id, std::string const& s) const {
    double count = getROISpecCount(ROI_id, s);
    return (count / math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setROISpecAmount(const std::string& ROI_id, std::string const& s, double amount) {
    setROISpecCount(ROI_id, s, amount * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////
// WEILIANG: Can we apply Sam's set count method (as in e.g. setCompCount) here?
void TetOpSplitP::setROISpecConc(const std::string& ROI_id, std::string const& s, double conc) {
    if (conc > UINT_MAX) {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << UINT_MAX << ").\n";
        ArgErrLog(os.str());
    }

    auto const& roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    std::vector<tetrahedron_global_id> apply_indices;

    double totalvol = 0.0;
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: roi->second) {
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

        apply_indices.push_back(tidx);
        totalvol += tet->vol();
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "Species " << s
            << " has not been defined in the following tetrahedrons, fill in zeros "
               "at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }

    int count = conc * (1.0e3 * totalvol * math::AVOGADRO);

    setROISpecCount(ROI_id, s, count);
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::getROISpecConc(const std::string& ROI_id, const std::string& s) const {
    auto const& roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    const double count = getROITetSpecCount(roi->second, s);
    double vol = getROIVol(ROI_id);
    return count / (1.0e3 * vol * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setROITriSpecClamped(const std::vector<triangle_global_id>& triangles,
                                       const std::string& s,
                                       bool b) {
    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: triangles) {
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
        if (tri->getInHost()) {
            tri->setClamped(slidx, b);
        }
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a patch, fill in "
               "zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "Species " << s
            << " has not been defined in the following triangles, fill in zeros at "
               "target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
}

void TetOpSplitP::setROITetSpecClamped(const std::vector<tetrahedron_global_id>& tetrahedrons,
                                       std::string const& s,
                                       bool b) {
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: tetrahedrons) {
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

        if (tet->getInHost()) {
            tet->setClamped(slidx, b);
        }
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, fill in zeros at target positions:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        CLOG(WARNING, "general_log")
            << "Species " << s
            << " has not been defined in the following tetrahedrons, fill in zeros "
               "at target positions:\n";
        CLOG(WARNING, "general_log") << spec_undefined.str() << "\n";
    }
}

void TetOpSplitP::setROISpecClamped(const std::string& ROI_id, std::string const& s, bool b) {
    {
        auto const& roi = mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id, 0, false);
        if (roi != mesh()->rois.end<tetmesh::ROI_TRI>()) {
            setROITriSpecClamped(roi->second, s, b);
            return;
        }
    }
    {
        auto const& roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id, 0, false);
        if (roi != mesh()->rois.end<tetmesh::ROI_TET>()) {
            setROITetSpecClamped(roi->second, s, b);
            return;
        }
    }
    std::ostringstream os;
    os << "Error: Cannot find suitable ROI for the function call getROICount.\n";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setROIReacK(const std::string& ROI_id, std::string const& r, double kf) {
    auto const& roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TET>()) {
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

        if (tet->getInHost()) {
            tet->reac(rlidx).setKcst(kf);
        }
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_reac_warning) {
        CLOG(WARNING, "general_log")
            << "Reac " << r
            << " has not been defined in the following tetrahedrons, no change is "
               "applied to them:\n";
        CLOG(WARNING, "general_log") << reac_undefined.str() << "\n";
    }

    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setROISReacK(const std::string& ROI_id, std::string const& sr, double kf) {
    auto const& roi = mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TRI>()) {
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

        if (tri->getInHost()) {
            tri->sreac(srlidx).setKcst(kf);
        }
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a patch, no "
               "change is applied to them:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_sreac_warning) {
        CLOG(WARNING, "general_log") << "SReac " << sr
                                     << " has not been defined in the following "
                                        "patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << sreac_undefined.str() << "\n";
    }

    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setROIDiffD(const std::string& ROI_id, std::string const& d, double dk) {
    auto const& roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TET>()) {
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

        if (tet->getInHost()) {
            tet->diff(dlidx).setDcst(dk);
        }
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_diff_warning) {
        CLOG(WARNING, "general_log")
            << "Diff " << d
            << " has not been defined in the following tetrahedrons, no change is "
               "applied to them:\n";
        CLOG(WARNING, "general_log") << diff_undefined.str() << "\n";
    }

    recomputeUpdPeriod = true;

    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setROIReacActive(const std::string& ROI_id, std::string const& r, bool a) {
    auto const& roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TET>()) {
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

        if (tet->getInHost()) {
            tet->reac(rlidx).setActive(a);
        }
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_reac_warning) {
        CLOG(WARNING, "general_log")
            << "Reac " << r
            << " has not been defined in the following tetrahedrons, no change is "
               "applied to them:\n";
        CLOG(WARNING, "general_log") << reac_undefined.str() << "\n";
    }

    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setROISReacActive(const std::string& ROI_id, std::string const& sr, bool a) {
    auto const& roi = mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TRI>()) {
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

        if (tri->getInHost()) {
            tri->sreac(srlidx).setActive(a);
        }
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a patch, no "
               "change is applied to them:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_sreac_warning) {
        CLOG(WARNING, "general_log") << "SReac " << sr
                                     << " has not been defined in the following "
                                        "patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << sreac_undefined.str() << "\n";
    }

    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setROIDiffActive(const std::string& ROI_id, std::string const& d, bool a) {
    auto const& roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TET>()) {
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

        if (tet->getInHost()) {
            tet->diff(dlidx).setActive(a);
        }
    }

    if (has_tet_warning) {
        CLOG(WARNING, "general_log") << "The following tetrahedrons have not been assigned to a "
                                        "compartment, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << tet_not_assign.str() << "\n";
    }

    if (has_diff_warning) {
        CLOG(WARNING, "general_log")
            << "Diff " << d
            << " has not been defined in the following tetrahedrons, no change is "
               "applied to them:\n";
        CLOG(WARNING, "general_log") << diff_undefined.str() << "\n";
    }

    recomputeUpdPeriod = true;

    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setROIVDepSReacActive(const std::string& ROI_id, std::string const& vsr, bool a) {
    auto const& roi = mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TRI>()) {
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

        if (tri->getInHost()) {
            tri->vdepsreac(vsrlidx).setActive(a);
        }
    }

    if (has_tri_warning) {
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a patch, no "
               "change is applied to them:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_vsreac_warning) {
        CLOG(WARNING, "general_log") << "VDepSReac " << vsr
                                     << " has not been defined in the following "
                                        "patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << vsreac_undefined.str() << "\n";
    }
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetOpSplitP::getROIReacExtent(const std::string& ROI_id,
                                                 std::string const& r) const {
    std::ostringstream os;
    os << "This function has not been implemented!";
    NotImplErrLog(os.str());

    auto const& roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TET>()) {
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
        CLOG(WARNING, "general_log")
            << "Reac " << r
            << " has not been defined in the following tetrahedrons, no change is "
               "applied to them:\n";
        CLOG(WARNING, "general_log") << reac_undefined.str() << "\n";
    }

    return sum;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::resetROIReacExtent(const std::string& ROI_id, std::string const& r) {
    std::ostringstream os;
    os << "This function has not been implemented!";
    NotImplErrLog(os.str());

    auto const& roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TET>()) {
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
        CLOG(WARNING, "general_log")
            << "Reac " << r
            << " has not been defined in the following tetrahedrons, no change is "
               "applied to them:\n";
        CLOG(WARNING, "general_log") << reac_undefined.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetOpSplitP::getROISReacExtent(const std::string& ROI_id,
                                                  std::string const& sr) const {
    std::ostringstream os;
    os << "This function has not been implemented!";
    NotImplErrLog(os.str());

    auto const& roi = mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TRI>()) {
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
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a patch, no "
               "change is applied to them:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_sreac_warning) {
        CLOG(WARNING, "general_log") << "SReac " << sr
                                     << " has not been defined in the following "
                                        "patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << sreac_undefined.str() << "\n";
    }

    return sum;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::resetROISReacExtent(const std::string& ROI_id, std::string const& sr) {
    std::ostringstream os;
    os << "This function has not been implemented!";
    NotImplErrLog(os.str());

    auto const& roi = mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TRI>()) {
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
        CLOG(WARNING, "general_log")
            << "The following triangles have not been assigned to a patch, no "
               "change is applied to them:\n";
        CLOG(WARNING, "general_log") << tri_not_assign.str() << "\n";
    }

    if (has_sreac_warning) {
        CLOG(WARNING, "general_log") << "SReac " << sr
                                     << " has not been defined in the following "
                                        "patch, no change is applied to them:\n";
        CLOG(WARNING, "general_log") << sreac_undefined.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetOpSplitP::getROIDiffExtent(const std::string& ROI_id,
                                                 std::string const& d) const {
    std::ostringstream os;
    os << "This function has not been implemented!";
    NotImplErrLog(os.str());

    auto const& roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TET>()) {
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
        CLOG(WARNING, "general_log")
            << "Diff " << d
            << " has not been defined in the following tetrahedrons, no change is "
               "applied to them:\n";
        CLOG(WARNING, "general_log") << diff_undefined.str() << "\n";
    }

    return sum;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::resetROIDiffExtent(const std::string& ROI_id, std::string const& d) {
    std::ostringstream os;
    os << "This function has not been implemented!";
    NotImplErrLog(os.str());

    auto const& roi = mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == mesh()->rois.end<tetmesh::ROI_TET>()) {
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
        CLOG(WARNING, "general_log")
            << "Diff " << d
            << " has not been defined in the following tetrahedrons, no change is "
               "applied to them:\n";
        CLOG(WARNING, "general_log") << diff_undefined.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

int TetOpSplitP::getTetHostRank(tetrahedron_global_id tidx) const {
    return _getHost(tidx);
}

////////////////////////////////////////////////////////////////////////////////

int TetOpSplitP::getTriHostRank(triangle_global_id tidx) const {
    return _getHost(tidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_updateLocal(uint* upd_entries, uint buffer_size) {
    for (uint i = 0; i < buffer_size; i++) {
        if (pKProcs[upd_entries[i]] != nullptr) {
            _updateElement(pKProcs[upd_entries[i]]);
        }
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::addDiff(Diff* diff) {
    diff->crData.pos = pDiffs.size();
    pDiffs.push_back(diff);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_updateDiff(Diff* diff) {
    double new_rate = diff->rate(this);

    diff->crData.rate = new_rate;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::addSDiff(SDiff* sdiff) {
    sdiff->crData.pos = pSDiffs.size();
    pSDiffs.push_back(sdiff);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_updateSDiff(SDiff* sdiff) {
    double new_rate = sdiff->rate(this);

    sdiff->crData.rate = new_rate;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::addNeighHost(int host) {
    neighbHosts.insert(host);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::registerBoundaryTet(Tet* tet) {
    boundaryTets.insert(tet);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::registerBoundaryTri(Tri* tri) {
    boundaryTris.insert(tri);
}

////////////////////////////////////////////////////////////////////////////////

uint TetOpSplitP::registerRemoteMoleculeChange(int svol_host,
                                               uint loc,
                                               SubVolType svol_type,
                                               unsigned long idx,
                                               solver::spec_local_id slidx,
                                               uint change) {
    uint new_loc = remoteChanges[svol_host].size();

    if (new_loc == 0 || new_loc - 4 < loc) {
        remoteChanges[svol_host].push_back(svol_type);
        remoteChanges[svol_host].push_back(idx);
        remoteChanges[svol_host].push_back(slidx.get());
        remoteChanges[svol_host].push_back(change);
    } else {
        uint stored_type = remoteChanges[svol_host][loc];
        uint stored_idx = remoteChanges[svol_host][loc + 1];
        uint stored_slidx = remoteChanges[svol_host][loc + 2];

        if (stored_type == svol_type && stored_idx == idx && stored_slidx == slidx) {
            remoteChanges[svol_host][loc + 3] += change;
            new_loc = loc;
        } else {
            remoteChanges[svol_host].push_back(svol_type);
            remoteChanges[svol_host].push_back(idx);
            remoteChanges[svol_host].push_back(slidx.get());
            remoteChanges[svol_host].push_back(change);
        }
    }
    return new_loc;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::_remoteSyncAndUpdate(void* requests,
                                       std::vector<KProc*>& applied_diffs,
                                       std::vector<int>& directions) {
#ifdef MPI_PROFILING
    double timing_start = MPI_Wtime();
#endif

    auto requestsPtr = static_cast<MPI_Request*>(requests);

    uint request_count = 0;
    for (auto& dest: neighbHosts) {
        MPI_Isend(remoteChanges[dest].data(),
                  remoteChanges[dest].size(),
                  MPI_UNSIGNED,
                  dest,
                  OPSPLIT_MOLECULE_CHANGE,
                  MPI_COMM_WORLD,
                  &(requestsPtr[request_count]));
        request_count++;
    }

    MPI_Status status;
    KProcPSet upd_kprocs;

    std::set<int> await_neighbors(neighbHosts);

#ifdef MPI_PROFILING
    double timing_end = MPI_Wtime();
    syncTime += (timing_end - timing_start);
#endif

    while (!await_neighbors.empty()) {
#ifdef MPI_PROFILING
        timing_start = MPI_Wtime();
#endif
        int flag = 0;
        int data_source = 0;
        for (auto& neighbor: await_neighbors) {
            MPI_Iprobe(neighbor, OPSPLIT_MOLECULE_CHANGE, MPI_COMM_WORLD, &flag, &status);
            if (flag != 0) {
                data_source = neighbor;
                break;
            }
        }
#ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        idleTime += (timing_end - timing_start);
#endif
        if (flag == 0) {
            continue;
        }

#ifdef MPI_PROFILING
        starttime = MPI_Wtime();
#endif
        // receive data
        int change_size = 0;
        MPI_Get_count(&status, MPI_UNSIGNED, &change_size);
        std::vector<uint> changes(change_size);
        MPI_Recv(changes.data(),
                 change_size,
                 MPI_UNSIGNED,
                 status.MPI_SOURCE,
                 OPSPLIT_MOLECULE_CHANGE,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        // apply changes
        uint nchanges = change_size / 4;
        for (uint c = 0; c < nchanges; c++) {
            uint type = changes[c * 4];
            uint idx = changes[c * 4 + 1];
            solver::spec_local_id slidx(changes[c * 4 + 2]);
            uint value = changes[c * 4 + 3];

            if (type == SUB_WM) {
                pWmVols[solver::comp_global_id(idx)]->incCount(slidx, value);
            }
            if (type == SUB_TET) {
                tetrahedron_global_id tetidx(idx);
                pTets[tetidx]->incCount(slidx, value);
                std::vector<KProc*> const& remote_upd = pTets[tetidx]->getSpecUpdKProcs(slidx);
                upd_kprocs.insert(remote_upd.begin(), remote_upd.end());
            }
            if (type == SUB_TRI) {
                triangle_global_id triidx(idx);
                pTris[triidx]->incCount(slidx, value);
                std::vector<KProc*> const& remote_upd = pTris[triidx]->getSpecUpdKProcs(slidx);
                upd_kprocs.insert(remote_upd.begin(), remote_upd.end());
            }
        }

        await_neighbors.erase(data_source);
#ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        syncTime += (timing_end - timing_start);
#endif
    }

#ifdef MPI_PROFILING
    timing_start = MPI_Wtime();
#endif

    auto napply = applied_diffs.size();
    for (uint i = 0; i < napply; i++) {
        KProc* kp = applied_diffs[i];
        int direction = directions[i];

        std::vector<KProc*> const& local_upd = kp->getLocalUpdVec(direction);

        for (auto& upd_kp: local_upd) {
            _updateElement(upd_kp);
        }
    }

    // update kprocs caused by remote molecule changes
    for (auto& upd_kp: upd_kprocs) {
        _updateElement(upd_kp);
    }
    _updateSum();

#ifdef MPI_PROFILING
    timing_end = MPI_Wtime();
    compTime += (timing_end - timing_start);
#endif
}

////////////////////////////////////////////////////////////////////////////

int TetOpSplitP::_getHost(tetrahedron_global_id tgidx) const {
    auto result = tetHosts.find(tgidx);
    if (result == tetHosts.end()) {
        std::ostringstream os;
        os << "Tetrahedron " << tgidx << " has not been assigned to a host.\n";
        ArgErrLog(os.str());
    }
    return result->second;
}

////////////////////////////////////////////////////////////////////////////////

int TetOpSplitP::_getHost(triangle_global_id tgidx) const {
    auto result = triHosts.find(tgidx);
    if (result == triHosts.end()) {
        std::ostringstream os;
        os << "Triangle " << tgidx << " has not been assigned to a host.\n";
        ArgErrLog(os.str());
    }
    return result->second;
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::repartitionAndReset(std::vector<int> const& tet_hosts,
                                      std::map<uint, int> const& tri_hosts,
                                      std::vector<int> const& wm_hosts) {
    pKProcs.clear();
    pVdepKProcs.clear();
    pDiffs.clear();
    pSDiffs.clear();
    neighbHosts.clear();
    boundaryTets.clear();
    boundaryTris.clear();

    auto ngroups = nGroups.size();
    for (uint i = 0; i < ngroups; i++) {
        free(nGroups[i]->indices);
        delete nGroups[i];
    }
    nGroups.clear();

    ngroups = pGroups.size();
    for (uint i = 0; i < ngroups; i++) {
        free(pGroups[i]->indices);
        delete pGroups[i];
    }
    pGroups.clear();

    tetHosts.clear();
    for (uint t = 0; t < tet_hosts.size(); t++) {
        tetHosts.emplace(tetrahedron_global_id(t), tet_hosts[t]);
    }

    triHosts.clear();
    for (auto [tri, host]: tri_hosts) {
        triHosts.emplace(triangle_global_id(tri), host);
    }
    wmHosts.container().assign(wm_hosts.begin(), wm_hosts.end());

    auto ntets = pTets.size();
    for (auto t: tetrahedron_global_id::range(ntets)) {
        if (pTets[t] == nullptr) {
            continue;
        }
        pTets[t]->repartition(this, myRank, tetHosts[t]);
    }

    for (auto wm: pWmVols.range()) {
        if (pWmVols[wm] == nullptr) {
            continue;
        }
        pWmVols[wm]->repartition(this, myRank, wmHosts[wm]);
    }

    for (auto& t: pTris) {
        if (t == nullptr) {
            continue;
        }
        t->repartition(this, myRank, triHosts[t->idx()]);
    }

    for (auto& t: pTets) {
        if ((t != nullptr) && t->getInHost()) {
            t->setupDeps();
        }
    }

    // Vector allows for all compartments to be well-mixed, so
    // hold null-pointer for mesh compartments
    for (auto& wmv: pWmVols) {
        if ((wmv != nullptr) && wmv->getInHost()) {
            wmv->setupDeps();
        }
    }
    // DEBUG: vector holds all possible triangles, but
    // only patch triangles are filled
    for (auto& t: pTris) {
        if ((t != nullptr) && t->getInHost()) {
            t->setupDeps();
        }
    }

    for (auto& tet: boundaryTets) {
        tet->setupBufferLocations();
    }
    for (auto& tri: boundaryTris) {
        tri->setupBufferLocations();
    }

    if (efflag()) {
        std::ostringstream os;
        os << "Repartition of EField is not implemented:\n";
        ArgErrLog(os.str());
    }

    neighbHosts.erase(myRank);
    nNeighbHosts = neighbHosts.size();

    // construct remote molecule change buffers
    remoteChanges.clear();
    for (auto neighbor: neighbHosts) {
        remoteChanges[neighbor] = {};
    }

    nEntries = pKProcs.size();
    diffSep = pDiffs.size();
    sdiffSep = pSDiffs.size();
    reset();
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

void TetOpSplitP::setDiffApplyThreshold(int threshold) {
    diffApplyThreshold = threshold;
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetOpSplitP::getReacExtent(bool local) {
    if (local) {
        return reacExtent;
    }
    unsigned long long sum = 0;
    MPI_Allreduce(&reacExtent, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetOpSplitP::getDiffExtent(bool local) {
    if (local) {
        return diffExtent;
    }
    unsigned long long sum = 0;
    MPI_Allreduce(&diffExtent, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::getNIteration() const {
    return nIteration;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::getCompTime() const {
    return compTime;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::getSyncTime() const {
    return syncTime;
}

////////////////////////////////////////////////////////////////////////////////
double TetOpSplitP::getIdleTime() const {
    return idleTime;
}

////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::getEFieldTime() const {
    return efieldTime;
}
////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::getRDTime() const {
    return rdTime;
}
////////////////////////////////////////////////////////////////////////////////

double TetOpSplitP::getDataExchangeTime() const {
    return dataExchangeTime;
}

}  // namespace steps::mpi::tetopsplit
