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

#include "mpi/tetvesicle/tetvesicle_rdef.hpp"

// Standard library & STL headers.
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mpi.h>
#include <vector>

// STEPS headers.

#include "mpi/tetvesicle/comp_rdef.hpp"
#include "mpi/tetvesicle/diff.hpp"
#include "mpi/tetvesicle/diffboundary.hpp"
#include "mpi/tetvesicle/endocytosis.hpp"
#include "mpi/tetvesicle/exocytosis.hpp"
#include "mpi/tetvesicle/ghkcurr.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/reac.hpp"
#include "mpi/tetvesicle/sdiff.hpp"
#include "mpi/tetvesicle/sreac.hpp"
#include "mpi/tetvesicle/tet_rdef.hpp"
#include "mpi/tetvesicle/tri_rdef.hpp"
#include "mpi/tetvesicle/vdepsreac.hpp"

#include "solver/chandef.hpp"
#include "solver/compdef.hpp"
#include "solver/ghkcurrdef.hpp"
#include "solver/linkspecdef.hpp"
#include "solver/ohmiccurrdef.hpp"
#include "solver/patchdef.hpp"
#include "solver/raftendocytosisdef.hpp"
#include "solver/statedef.hpp"
#include "solver/vdepsreacdef.hpp"
#include "solver/vessdiffdef.hpp"
#include "solver/vessreacdef.hpp"

#include "solver/efield/dVsolver.hpp"
#include "solver/efield/efield.hpp"
#ifdef USE_PETSC
#include "solver/efield/dVsolver_petsc.hpp"
#endif

#include "geom/memb.hpp"
#include "geom/tetmesh.hpp"
#include "geom/tmcomp.hpp"
#include "math/constants.hpp"

#include "math/point.hpp"
#include "math/tetrahedron.hpp"
#include "util/checkpointing.hpp"
#include "util/collections.hpp"
#include "util/distribute.hpp"
#include "util/error.hpp"

#include <fau.de/overlap.hpp>

namespace steps::mpi::tetvesicle {

TetVesicleRDEF::TetVesicleRDEF(model::Model* m, wm::Geom* g, const rng::RNGptr& r, int calcMembPot)
    : API(*m, *g, r)
    , pMesh(nullptr)
    , pA0(0.0)
    , pEFoption(static_cast<EF_solver>(calcMembPot))
    , pTemp(0.0)
    , pEField(nullptr)
    , pEFDT(1.0e-5)
    , pEFNVerts(0)
    , pEFNTris(0)
    , pEFTris_vec(0)
    , pEFNTets(0) {
    if (rng() == nullptr) {
        std::ostringstream os;
        os << "No RNG provided to solver initializer function";
        ArgErrLog(os.str());
    }

    // All initialization code now in _setup() to allow EField solver to be
    // derived and create EField local objects within the constructor

    std::cout.setf(std::ios::unitbuf);

    pMesh = dynamic_cast<tetmesh::Tetmesh*>(&geom());
    if (pMesh == nullptr) {
        ArgErrLog(
            "Geometry description to solver::TetVesicleRDEF solver "
            "constructor is not a valid tetmesh::Tetmesh object.");
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &myRank_World);
    MPI_Comm_size(MPI_COMM_WORLD, &nHosts_World);
    nHosts_RDEF = nHosts_World - 1;
    _partition();
    MPI_Barrier(MPI_COMM_WORLD);
    dataTypeUtil.commitAllDataTypes();
    _setup();
    // may remove this barrier later.
    MPI_Barrier(MPI_COMM_WORLD);
    _updateLocal();
    MPI_Barrier(MPI_COMM_WORLD);

    // assign the start of the LinkSpec individual IDs uniquely for this solver.
    // VesRaft gets the first block, half the total, so there is no overlap with RDEF cores.
    index_t frac = std::numeric_limits<index_t>::max() / (2 * nHosts_RDEF);
    pNextLinkSpecUniqueID.set((std::numeric_limits<index_t>::max() / 2) + myRank_RDEF * frac);
    pMaxLinkSpecUniqueID.set((std::numeric_limits<index_t>::max() / 2) + (myRank_RDEF + 1) * frac -
                             1);
    // PS IDs don't have to be unique across RDEF ranks, and they will be reassigned by VesRaft.
    // Since these are smaller changes that can be reset, we only use 10% of the available numbers.
    pNextPointSpecUniqueID.set(std::numeric_limits<index_t>::max() * 0.9);
    pMaxPointSpecUniqueID.set(std::numeric_limits<index_t>::max());
}

////////////////////////////////////////////////////////////////////////////////

TetVesicleRDEF::~TetVesicleRDEF() {
    for (auto c: pComps) {
        delete c;
    }
    for (auto p: pPatches) {
        delete p;
    }
    for (auto db: pDiffBoundaries) {
        delete db;
    }
    for (auto t: pTets) {
        delete t;
    }
    for (auto t: pTris) {
        delete t;
    }
    for (auto g: nGroups) {
        g->free_indices();
        delete g;
    }
    for (auto g: pGroups) {
        g->free_indices();
        delete g;
    }

    dataTypeUtil.freeAllDataTypes();
}

///////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::checkpoint(std::string const& file_name) {
    std::fstream cp_file;

    cp_file.open(file_name.c_str(), std::fstream::out | std::fstream::binary | std::fstream::trunc);
    util::checkpoint(cp_file, myRank_World);

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

    // Kprocs cannot change order so simple c & r
    for (auto const& k: pKProcs) {
        if (k != nullptr) {
            k->checkpoint(cp_file);
        }
    }
    // NOTE: pDiffs and pSDiffs are encapsulated in pKProcs.

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

    if (_efflag()) {
        util::checkpoint(cp_file, pTemp);
        pEField->checkpoint(cp_file);
        util::checkpoint(cp_file, pEFDT);
    }

    util::checkpoint(cp_file, pNextLinkSpecUniqueID);
    util::checkpoint(cp_file, pNextPointSpecUniqueID);

    util::checkpoint(cp_file, recomputeUpdPeriod);
    util::checkpoint(cp_file, reacExtent);
    util::checkpoint(cp_file, diffExtent);
    util::checkpoint(cp_file, nIteration);
    util::checkpoint(cp_file, updPeriod);
    util::checkpoint(cp_file, diffApplyThreshold);
}

///////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::restore(std::string const& file_name) {
    std::fstream cp_file;

    cp_file.open(file_name.c_str(), std::fstream::in | std::fstream::binary);

    cp_file.seekg(0);
    util::compare(cp_file, myRank_World);

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

    // Kprocs cannot change order so simple c & r
    for (auto const& k: pKProcs) {
        if (k != nullptr) {
            k->restore(cp_file);
        }
    }

    // NOTE pDiffs and pSDiffs have changed order during simulation from time of setup
    // but they have been restored in pKProcs. Order in pDiffs, pSDiffs needs to be updated

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

    if (_efflag()) {
        util::restore(cp_file, pTemp);
        pEField->restore(cp_file);
        util::restore(cp_file, pEFDT);
        _refreshEFTrisV();
    }

    util::restore(cp_file, pNextLinkSpecUniqueID);
    util::restore(cp_file, pNextPointSpecUniqueID);

    util::restore(cp_file, recomputeUpdPeriod);
    util::restore(cp_file, reacExtent);
    util::restore(cp_file, diffExtent);
    util::restore(cp_file, nIteration);
    util::restore(cp_file, updPeriod);
    util::restore(cp_file, diffApplyThreshold);

    cp_file.close();

    // Sync with vesraft and do a global update
    _useVesV2R();
    _useRaftV2R();
    _updateLocal();
    recomputeUpdPeriod = true;
}

////////////////////////////////////////////////////////////////////////////////

std::string TetVesicleRDEF::getSolverName() const {
    return "TetVesicleRDEF";
}

////////////////////////////////////////////////////////////////////////////////

std::string TetVesicleRDEF::getSolverDesc() const {
    return "Parallel approximate stochastic solver with added Vesicle-related "
           "functionality.";
}

////////////////////////////////////////////////////////////////////////////////

std::string TetVesicleRDEF::getSolverAuthors() const {
    return "Iain Hepburn, Weiliang Chen, Stefan Wils, Blue Brain Project.";
}

////////////////////////////////////////////////////////////////////////////////

std::string TetVesicleRDEF::getSolverEmail() const {
    return "steps.dev@gmail.com";
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_partition() {
    if (myRank_World == 0) {
        ProgErr("A TetVesicleRDEF solver is created in the VesRaft rank.");
    } else {
        vesraftRank_World = 0;
        RDEFmasterRank_World = 1;
        RDEFmasterRank_RDEF = 0;
        myRank_RDEF = myRank_World - 1;
        MPI_Comm_split(MPI_COMM_WORLD, 1, myRank_RDEF, &RDEFComm);
    }

    uint n_tets_tris[2];
    std::vector<tetrahedron_global_id> tet_ids;
    std::vector<int> tet_hosts;

    std::vector<triangle_global_id> tri_ids;
    std::vector<int> tri_hosts;

    MPI_Bcast(n_tets_tris, 2, MPI_UNSIGNED, vesraftRank_World, MPI_COMM_WORLD);

    tet_ids.resize(n_tets_tris[0]);
    tet_hosts.resize(n_tets_tris[0]);
    tri_ids.resize(n_tets_tris[1]);
    tri_hosts.resize(n_tets_tris[1]);

    MPI_Bcast(tet_ids.data(), n_tets_tris[0], MPI_STEPS_INDEX, vesraftRank_World, MPI_COMM_WORLD);
    MPI_Bcast(tet_hosts.data(), n_tets_tris[0], MPI_INT, vesraftRank_World, MPI_COMM_WORLD);
    MPI_Bcast(tri_ids.data(), n_tets_tris[1], MPI_STEPS_INDEX, vesraftRank_World, MPI_COMM_WORLD);
    MPI_Bcast(tri_hosts.data(), n_tets_tris[1], MPI_INT, vesraftRank_World, MPI_COMM_WORLD);

    for (uint t = 0; t < n_tets_tris[0]; t++) {
        tetHosts[tet_ids[t]] = tet_hosts[t];
    }
    for (uint t = 0; t < n_tets_tris[1]; t++) {
        triHosts[tri_ids[t]] = tri_hosts[t];
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setup() {
    // First initialise the pTets, pTris vector, because
    // want tets and tris to maintain indexing from Geometry
    uint ntets = pMesh->countTets();
    uint ntris = pMesh->countTris();

    pTets.container().assign(ntets, nullptr);
    pTris.container().assign(ntris, nullptr);

    // Now create the actual compartments.
    for (auto const& c: statedef().comps()) {
        solver::comp_global_id compdef_gidx = c->gidx();
        uint comp_idx = _addComp(c.get(), pMesh);
        AssertLog(compdef_gidx == comp_idx);
    }
    // Create the actual patches.
    for (auto const& p: statedef().patches()) {
        solver::patch_global_id patchdef_gidx = p->gidx();
        uint patch_idx = _addPatch(p.get(), pMesh);
        AssertLog(patchdef_gidx == patch_idx);
    }

    // Create the diffusion boundaries
    for (auto const& db: statedef().diffBoundaries()) {
        solver::diffboundary_global_id diffboundary_gidx = db->gidx();
        uint diffb_idx = _addDiffBoundary(db.get());
        AssertLog(diffboundary_gidx == diffb_idx);
    }

    // Create the surface diffusion boundaries
    for (auto const& sdb: statedef().sdiffBoundaries()) {
        solver::sdiffboundary_global_id sdiffboundary_gidx = sdb->gidx();
        uint sdiffb_idx = _addSDiffBoundary(sdb.get());
        AssertLog(sdiffboundary_gidx == sdiffb_idx);
    }

    uint npatches = pPatches.size();
    AssertLog(pMesh->_countPatches() == npatches);

    for (auto p: solver::patch_global_id::range(npatches)) {
        // Add the tris for this patch
        // We have checked the indexing - p is the global index
        wm::Patch& wmpatch = pMesh->_getPatch(p);

        // Perform upcast
        auto tmpatch = dynamic_cast<tetmesh::TmPatch*>(&wmpatch);
        if (tmpatch == nullptr) {
            ArgErrLog(
                "Well-mixed patches not supported in solver::TetVesicleRDEF "
                "solver.");
        }
        PatchRDEF* localpatch = pPatches[p];

        // Create a map between edges and adjacent tris in this patch
        std::map<bar_id_t, std::vector<triangle_global_id>> bar2tri;

        // We need to go through all patches to record bar2tri mapping
        // for all connected triangle neighbors even they are in different
        // patches, because their information is needed for surface diffusion
        // boundary

        for (auto bar_p: solver::patch_global_id::range(npatches)) {
            auto& patch = pMesh->_getPatch(bar_p);
            auto* bar_patch = dynamic_cast<tetmesh::TmPatch*>(&patch);
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

        auto const& tri_ids = tmpatch->_getAllTriIndices();
        for (auto i = 0u; i < tri_ids.size(); ++i) {
            auto tri = tri_ids[i];

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
            const math::point3d& trinorm = pMesh->_getTriNorm(tri);

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
                    tris[2],
                    baryc,
                    trinorm);
        }
    }

    const auto ncomps = pComps.size();
    AssertLog(pMesh->_countComps() == ncomps);

    for (auto c: solver::comp_global_id::range(ncomps)) {
        // Now add the tets for this comp
        // We have checked the indexing- c is the global index
        wm::Comp& wmcomp = pMesh->_getComp(c);

        // Perform upcast
        auto* tmcomp = dynamic_cast<tetmesh::TmComp*>(&wmcomp);
        if (tmcomp != nullptr) {
            CompRDEF* localcomp = pComps[c];

            for (auto const& tet: tmcomp->_getAllTetIndices()) {
                AssertLog(pMesh->getTetComp(tet) == tmcomp);

                double vol = pMesh->getTetVol(tet);

                const auto& tris = pMesh->_getTetTriNeighb(tet);

                std::array<double, 4> areas{0, 0, 0, 0};
                for (uint j = 0; j < areas.size(); ++j) {
                    areas[j] = pMesh->getTriArea(tris[j]);
                }

                const auto tets = pMesh->_getTetTetNeighb(tet);
                const auto& baryc = pMesh->_getTetBarycenter(tet);

                std::array<double, 4> distances = {0, 0, 0, 0};
                for (uint j = 0; j < distances.size(); ++j) {
                    if (tets[j] == tetrahedron_global_id::unknown_value()) {
                        continue;
                    }
                    distances[j] = math::distance(baryc, pMesh->_getTetBarycenter(tets[j]));
                }

                _addTet(tet,
                        localcomp,
                        vol,
                        areas[0],
                        areas[1],
                        areas[2],
                        areas[3],
                        distances[0],
                        distances[1],
                        distances[2],
                        distances[3],
                        tets[0],
                        tets[1],
                        tets[2],
                        tets[3],
                        baryc);
            }

            /// TODO?
            // SETUP VESICLE Proxies

        } else {
            ProgErrLog("Well-mixed compartments not supported for vesicle simulations.");
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
            tetrahedron_global_id tet_neighbor = pTets[t]->tet(j);
            if (tet_neighbor.valid() && pTets[tet_neighbor] != nullptr) {
                pTets[t]->setNextTet(j, pTets[tet_neighbor]);
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
            triangle_global_id tri_neighb = pTris[t]->tri(j);
            if (tri_neighb.valid() && pTris[tri_neighb] != nullptr) {
                pTris[t]->setNextTri(j, pTris[tri_neighb]);
            }
        }

        // By convention, triangles in a patch should have an inner tetrahedron
        // defined (neighbouring tets 'flipped' if necessary in Tetmesh) but not
        // necessarily an outer tet 17/3/10- actually this is not the case any more
        // with well-mixed compartments
        //
        tetrahedron_global_id tetinner = pTris[t]->tet(0);
        tetrahedron_global_id tetouter = pTris[t]->tet(1);

        // Now inside and outside tetrahedrons must be normal tetrahedrons

        if (tetinner.valid()) {
            // NEW FOR THIS VERSION: Tris store index of inner and outer tet (outer
            // may not exist if on surface) but tets may not belong to a compartment,
            // even inner tets now since they may be well-mixed compartments
            //
            if (pTets[tetinner] != nullptr) {
                // A triangle may already have an inner tet defined as a well-mixed
                // volume, but that should not be the case here:
                AssertLog(pTris[t]->iTet() == nullptr);

                pTris[t]->setInnerTet(pTets[tetinner]);
                // Now add this triangle to inner tet's list of neighbours
                for (uint i = 0; i <= 4; ++i) {
                    // include AssertLog for debugging purposes and remove
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
                    TetRDEF* tet_in = pTets[tetinner];
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
                    TetRDEF* tet_out = pTets[tetouter];

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

    /// TODO? setup raft proxies

    // Need to setup Rafts here after the triangles have been setup, because
    // Patches are going to create lookup tables depending on the connectivity
    // for (auto const &patch : pPatches) {
    // All patches can contain all rafts, for now
    // patch->setupRafts();
    //}

    // Now loop over the diffusion boundaries:
    // 1) get all the triangles and get the two tetrahedrons
    // 2) figure out which direction is the direction for a tetrahedron
    // 3) add the tetrahedron and the direction to local object

    // This is here because we need all tets to have been assigned correctly
    // to compartments. Check every one and set the compA and compB for the db
    uint ndiffbnds = pDiffBoundaries.size();
    AssertLog(ndiffbnds == pMesh->_countDiffBoundaries());

    for (auto db: solver::diffboundary_global_id::range(ndiffbnds)) {
        DiffBoundary* localdiffb = pDiffBoundaries[db];

        solver::comp_global_id compAidx = localdiffb->def()->compa();
        solver::comp_global_id compBidx = localdiffb->def()->compb();
        const auto& compAdef = statedef().compdef(compAidx);
        const auto& compBdef = statedef().compdef(compBidx);

        for (auto dbtri: localdiffb->def()->tris()) {
            const auto tri_tets = pMesh->_getTriTetNeighb(dbtri);

            auto tetAidx = tri_tets[0];
            auto tetBidx = tri_tets[1];
            AssertLog(tetAidx != tetrahedron_global_id::unknown_value() &&
                      tetBidx != tetrahedron_global_id::unknown_value());

            TetRDEF* tetA = _tet(tetAidx);
            TetRDEF* tetB = _tet(tetBidx);
            AssertLog(tetA != nullptr && tetB != nullptr);

            solver::Compdef* tetA_cdef = tetA->compdef();
            solver::Compdef* tetB_cdef = tetB->compdef();
            AssertLog(tetA_cdef != nullptr);
            AssertLog(tetB_cdef != nullptr);

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
    AssertLog(nsdiffbnds == _mesh()->_countSDiffBoundaries());

    for (auto sdb: solver::sdiffboundary_global_id::range(nsdiffbnds)) {
        SDiffBoundary* localsdiffb = pSDiffBoundaries[sdb];

        solver::patch_global_id patchAidx = localsdiffb->def()->patcha();
        solver::patch_global_id patchBidx = localsdiffb->def()->patchb();
        const auto& patchAdef = statedef().patchdef(patchAidx);
        const auto& patchBdef = statedef().patchdef(patchBidx);

        for (auto sdbbar: localsdiffb->def()->bars()) {
            const auto* bar_tris = pMesh->_getBarTriNeighb(bar_id_t(sdbbar));

            auto triAidx = bar_tris[0];
            auto triBidx = bar_tris[1];
            AssertLog(triAidx != triangle_global_id::unknown_value() &&
                      triBidx != triangle_global_id::unknown_value());

            TriRDEF* triA = _tri(triAidx);
            TriRDEF* triB = _tri(triBidx);
            AssertLog(triA != nullptr && triB != nullptr);

            solver::Patchdef* triA_pdef = triA->patchdef();
            solver::Patchdef* triB_pdef = triB->patchdef();
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

    for (auto& t: pTris) {
        if (t != nullptr) {
            t->setupKProcs(this, _efflag());
        }
    }

    // Resolve all dependencies
    for (auto t: pTets) {
        if ((t != nullptr) && t->getInHost()) {
            t->setupDeps();
        }
    }

    for (auto t: pTris) {
        if ((t != nullptr) && t->getInHost()) {
            t->setupDeps();
        }
    }

    // Create EField structures if EField is to be calculated
    if (_efflag() == true) {
        _setupEField();
    }

    for (auto& tet: boundaryTets) {
        tet->setupBufferLocations();
    }
    for (auto& tri: boundaryTris) {
        tri->setupBufferLocations();
    }
    // just in case
    neighbHosts.erase(myRank_World);
    nNeighbHosts = neighbHosts.size();

    // construct remote molecule change buffers
    remoteChanges.clear();
    for (auto& neighbor: neighbHosts) {
        remoteChanges[neighbor] = {};
    }

    nEntries = pKProcs.size();

    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setupEField() {
    using solver::efield::dVSolverBanded;
    using solver::efield::make_EField;

    //// Note to self: for now roughly following flow from original code in
    /// sim/controller.py. / code for setting up a mesh was in func_tetmesh
    /// constructor and called functions.

    AssertLog(_efflag() == true);

    switch (pEFoption) {
    case EF_DEFAULT:
    case EF_DV_BDSYS:
        pEField = make_EField<dVSolverBanded>();
        break;
#ifdef USE_PETSC
    case EF_DV_PETSC:
        using solver::efield::dVSolverPETSC;
        pEField = make_EField<dVSolverPETSC>(RDEFComm);
        break;
#endif
    default:
        ArgErrLog("Unsupported E-Field solver.");
    }

    // Give temperature a default value of 20c
    pTemp = 293.15;

    uint nmembs = _mesh()->_countMembs();

    if (nmembs != 1) {
        std::ostringstream os;
        os << "Membrane potential solver currently supports only one ";
        os << "membrane description object.";
        ArgErrLog(os.str());
    }

    tetmesh::Memb* memb = _mesh()->_getMemb(0);
    AssertLog(memb != nullptr);

    pEFNTets = memb->countVolTets();
    pEFNTris = memb->countTris();
    pEFNVerts = memb->countVerts();

    std::vector<vertex_id_t> EFTets(pEFNTets * 4);

    // All the triangles we will count here are the real membrane triangles,
    // virtual triangles will not require a capacitance.
    std::vector<vertex_id_t> EFTris(pEFNTris * 3);

    std::vector<double> EFVerts(pEFNVerts * 3);

    auto nverts = _mesh()->countVertices();
    auto ntris = _mesh()->countTris();
    auto ntets = _mesh()->countTets();

    pEFVert_GtoL.resize(nverts);
    pEFTri_GtoL.container().resize(ntris);
    pEFTet_GtoL.container().resize(ntets);
    pEFTri_LtoG.container().resize(_neftris());

    // Copy the data to local structures.

    auto const& membverts = memb->_getAllVertIndices();
    AssertLog(membverts.size() == _nefverts());
    for (uint efv = 0; efv < _nefverts(); ++efv) {
        const auto vertidx = membverts[efv];
        math::point3d verttemp = _mesh()->_getVertex(vertidx);
        uint efv2 = efv * 3;

        // CONVERTING TO MICRONS HERE. EFIELD OBJECT WILL NOT PERFORM THIS
        // CONVERSION
        verttemp *= 1.0e6;
        EFVerts[efv2] = verttemp[0];
        EFVerts[efv2 + 1] = verttemp[1];
        EFVerts[efv2 + 2] = verttemp[2];

        pEFVert_GtoL[vertidx.get()].set(efv);
    }

    const auto& membtets = memb->_getAllVolTetIndices();
    AssertLog(membtets.size() == _neftets());
    for (uint eft = 0; eft < _neftets(); ++eft) {
        auto const& tetidx = membtets[eft];
        const auto& tettemp = _mesh()->_getTet(tetidx);
        uint eft2 = eft * 4;

        // Convert to indices used by EField object
        auto tv0 = pEFVert_GtoL[tettemp[0].get()];
        auto tv1 = pEFVert_GtoL[tettemp[1].get()];
        auto tv2 = pEFVert_GtoL[tettemp[2].get()];
        auto tv3 = pEFVert_GtoL[tettemp[3].get()];
        if (tv0.unknown() || tv1.unknown() || tv2.unknown() || tv3.unknown()) {
            std::ostringstream os;
            os << "Failed to create EField structures.";
            ProgErrLog(os.str());
        }

        EFTets[eft2] = tv0;
        EFTets[eft2 + 1] = tv1;
        EFTets[eft2 + 2] = tv2;
        EFTets[eft2 + 3] = tv3;

        pEFTet_GtoL[tetidx].set(eft);
    }

    const auto& membtris = memb->_getAllTriIndices();
    AssertLog(membtris.size() == _neftris());

    pEFTris_vec.container().resize(_neftris());

    EFTrisV.container().resize(pEFNTris);

    EFTrisI_permuted.resize(pEFNTris);
    EFTrisI_idx.resize(pEFNTris);

    EFTrisI_offset.assign(nHosts_RDEF, 0);
    EFTrisI_count.assign(nHosts_RDEF, 0);

    std::vector<index_t> local_eftri_indices;

    for (uint eft = 0; eft < _neftris(); ++eft) {
        auto triidx = membtris[eft];
        const auto& tritemp = _mesh()->_getTri(triidx);
        uint eft2 = eft * 3;

        // Convert to indices used by EField object
        auto tv0 = pEFVert_GtoL[tritemp[0].get()];
        auto tv1 = pEFVert_GtoL[tritemp[1].get()];
        auto tv2 = pEFVert_GtoL[tritemp[2].get()];
        if (tv0.unknown() || tv1.unknown() || tv2.unknown()) {
            std::ostringstream os;
            os << "Failed to create EField structures.";
            ProgErrLog(os.str());
        }

        EFTris[eft2] = tv0;
        EFTris[eft2 + 1] = tv1;
        EFTris[eft2 + 2] = tv2;

        pEFTri_GtoL[triidx].set(eft);
        pEFTri_LtoG[triangle_local_id(eft)] = triidx;

        // This is added now for quicker iteration during run()
        // Extremely important for larger meshes, orders of magnitude times faster
        TriRDEF* tri_p = pTris[triidx];
        pEFTris_vec[triangle_local_id(eft)] = tri_p;

        int tri_host = tri_p->getHost();
        ++EFTrisI_count[tri_host - 1];
        if (myRank_World == tri_host) {
            local_eftri_indices.push_back(eft);
        }
    }

    const int* count_begin = EFTrisI_count.data();
    std::partial_sum(count_begin, count_begin + (nHosts_RDEF - 1), EFTrisI_offset.data() + 1);

    AssertLog(local_eftri_indices.size() == static_cast<uint>(EFTrisI_count[myRank_RDEF]));

    MPI_Allgatherv(local_eftri_indices.data(),
                   static_cast<int>(local_eftri_indices.size()),
                   MPI_STEPS_INDEX,
                   EFTrisI_idx.data(),
                   EFTrisI_count.data(),
                   EFTrisI_offset.data(),
                   MPI_STEPS_INDEX,
                   RDEFComm);

    // pEField->initMesh(_nefverts(), EFVerts, _neftris(), pEFTris, _neftets(),
    // pEFTets, memb->_getOpt_method(), memb->_getOpt_file_name(),
    // memb->_getSearch_percent());
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

// 'Safe' global to local index translation methods that throw on error.

inline solver::spec_local_id TetVesicleRDEF::_specG2L_or_throw(CompRDEF* comp,
                                                               solver::spec_global_id gidx) const {
    AssertLog(gidx < statedef().countSpecs());
    solver::spec_local_id lidx = comp->def()->specG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("species undefined in compartment");
    }
    return lidx;
}

inline solver::spec_local_id TetVesicleRDEF::_specG2L_or_throw(PatchRDEF* patch,
                                                               solver::spec_global_id gidx) const {
    AssertLog(gidx < statedef().countSpecs());
    solver::spec_local_id lidx = patch->def()->specG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("species undefined in patch");
    }
    return lidx;
}

#if 0
inline solver::spec_local_id TetVesicleRDEF::_specG2L_or_throw(TetRDEF *tet, solver::spec_global_id gidx) const {
    AssertLog(gidx < statedef().countSpecs());
    solver::spec_local_id lidx = tet->compdef()->specG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("species undefined in tetrahedron");
    }
    return lidx;
}

inline solver::spec_local_id TetVesicleRDEF::_specG2L_or_throw(TriRDEF *tri, solver::spec_global_id gidx) const {
    AssertLog(gidx < statedef().countSpecs());
    solver::spec_local_id lidx = tri->patchdef()->SPEC_LIDX_UNDEFINED(gidx);

    if (lidx.unknown()) {
        ArgErrLog("species undefined in triangle");
    }
    return lidx;
}
#endif

inline solver::reac_local_id TetVesicleRDEF::_reacG2L_or_throw(CompRDEF* comp,
                                                               solver::reac_global_id gidx) const {
    AssertLog(gidx < statedef().countReacs());
    solver::reac_local_id lidx = comp->def()->reacG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("reaction undefined in compartment");
    }
    return lidx;
}

inline solver::sreac_local_id TetVesicleRDEF::_sreacG2L_or_throw(
    PatchRDEF* patch,
    solver::sreac_global_id gidx) const {
    AssertLog(gidx < statedef().countSReacs());
    solver::sreac_local_id lidx = patch->def()->sreacG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("surface reaction undefined in patch");
    }
    return lidx;
}

inline solver::diff_local_id TetVesicleRDEF::_diffG2L_or_throw(CompRDEF* comp,
                                                               solver::diff_global_id gidx) const {
    AssertLog(gidx < statedef().countDiffs());
    solver::diff_local_id lidx = comp->def()->diffG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("diffusion rule undefined in compartment");
    }
    return lidx;
}

inline solver::vdepsreac_local_id TetVesicleRDEF::_vdepsreacG2L_or_throw(
    PatchRDEF* patch,
    solver::vdepsreac_global_id gidx) const {
    AssertLog(gidx < statedef().countVDepSReacs());
    solver::vdepsreac_local_id lidx = patch->def()->vdepsreacG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("voltage-dependent surface reation undefined in patch");
    }
    return lidx;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_addComp(solver::Compdef* cdef, tetmesh::Tetmesh* mesh) {
    auto* comp = new CompRDEF(cdef, mesh, this);
    AssertLog(comp != nullptr);
    uint compidx = pComps.size();
    pComps.container().push_back(comp);
    return compidx;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_addPatch(solver::Patchdef* pdef, tetmesh::Tetmesh* mesh) {
    auto* patch = new PatchRDEF(pdef, mesh, this);
    AssertLog(patch != nullptr);
    uint patchidx = pPatches.size();
    pPatches.container().push_back(patch);
    return patchidx;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_addDiffBoundary(solver::DiffBoundarydef* dbdef) {
    auto diffb = new DiffBoundary(dbdef);
    auto dbidx = pDiffBoundaries.size();
    pDiffBoundaries.container().push_back(diffb);
    return dbidx;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_addSDiffBoundary(solver::SDiffBoundarydef* sdbdef) {
    auto sdiffb = new SDiffBoundary(sdbdef);
    auto sdbidx = pSDiffBoundaries.size();
    pSDiffBoundaries.container().push_back(sdiffb);
    return sdbidx;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_addTet(tetrahedron_global_id tetidx,
                             CompRDEF* comp,
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
                             tetrahedron_global_id tet3,
                             const math::point3d& baryc) {
    solver::Compdef* compdef = comp->def();
    auto* localtet = new TetRDEF(tetidx,
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
                                 baryc,
                                 myRank_World,
                                 tetHosts[tetidx]);
    AssertLog(localtet != nullptr);
    AssertLog(tetidx < pTets.size());
    AssertLog(pTets[tetidx] == nullptr);
    pTets[tetidx] = localtet;
    comp->addTet(localtet);

    // MPISTEPS
    localtet->setSolverRDEF(this);
    // MPISTEPS
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_addTri(triangle_global_id triidx,
                             PatchRDEF* patch,
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
                             triangle_global_id tri2,
                             const math::point3d& baryc,
                             const math::point3d& trinorm) {
    solver::Patchdef* patchdef = patch->def();
    auto* tri = new TriRDEF(triidx,
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
                            baryc,
                            trinorm,
                            myRank_World,
                            triHosts[triidx]);
    AssertLog(tri != nullptr);
    AssertLog(triidx.get() < pTris.size());
    AssertLog(pTris[triidx] == nullptr);
    pTris[triidx] = tri;
    patch->addTri(tri);

    // MPISTEPS
    tri->setSolverRDEF(this);
    // MPISTEPS
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::reset() {
    for (auto const& comp: pComps) {
        comp->reset();
    }

    for (auto const& patch: pPatches) {
        patch->reset();
    }

    if (_efflag() == true) {
        pEField->setMembPotential(solver::membrane_global_id(0), DEFAULT_MEMB_POT);
    }

    for (auto const& tet: pTets) {
        if (tet == nullptr or !tet->getInHost()) {
            continue;
        }
        tet->reset();
    }

    for (auto const& tri: pTris) {
        if (tri == nullptr or !tri->getInHost()) {
            continue;
        }
        tri->reset();
    }

    uint ngroups = nGroups.size();
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

    diffSep = 0;
    sdiffSep = 0;

    statedef().resetTime();
    statedef().resetNSteps();
    _updateLocal();  // why?
    // for zero order reactions the rate is non-zero from the start
    // so need to add to the CR search

    // Have to reset the kcsts in VesSreacdef objects. Can't do it via the
    // local VesSReacs because they are dynamic, there may be none in existance.
    // HAVE to go through statedef
    for (auto const& vssr: statedef().vessreacs()) {
        vssr->reset();
    }

    for (auto const& exo: statedef().exocytosiss()) {
        exo->reset();
    }

    for (auto const& raftendo: statedef().raftendocytosiss()) {
        raftendo->reset();
    }

    // MPI STUFF
    reacExtent = 0.0;
    diffExtent = 0.0;
    nIteration = 0.0;

    compTime = 0.0;
    syncTime = 0.0;
    idleTime = 0.0;

    efieldTime = 0.0;
    rdTime = 0.0;
    dataExchangeTime = 0.0;

    recomputeUpdPeriod = true;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::run(double endtime) {
    ArgErrLogIf(endtime < statedef().time(), "Endtime is before current simulation time ");

    // At this point, all the tet and tri counts are up do date because calls like
    // setTetCount are dealt with by this solver, but VesRaft may have updates for
    // vesicles and rafts.
    //  only do these dramatic updates if necessary

    bool require_vesicle_comm;
    MPI_Bcast(&require_vesicle_comm, 1, MPI_C_BOOL, vesraftRank_World, MPI_COMM_WORLD);

    if (require_vesicle_comm) {
        _useVesV2R();
        _useRaftV2R();
        _updateLocal();  // necessary?
        recomputeUpdPeriod = true;
    }

    while (true) {
        double ves_endtime;
        MPI_Bcast(&ves_endtime, 1, MPI_DOUBLE, vesraftRank_World, MPI_COMM_WORLD);

        if (ves_endtime < 0.0) {
            // negative number from VesRaft means we are finished
            statedef().setTime(endtime);
            return;
        } else {
            AssertLog(!(ves_endtime > endtime));
            _runOneVesicleStep(ves_endtime);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_runOneVesicleStep(double ves_endtime) {
    if (recomputeUpdPeriod) {
        _computeUpdPeriod();
    }

    if (_efflag()) {
        _runWithEField(ves_endtime);
    } else {
        _runWithoutEField(ves_endtime);
    }

    MPI_Barrier(RDEFComm);  // Necessary??

    // >>>>>>>>>>>>> COMMUNICATION, RDEF -> VESRAFT
    _constructVesR2V();
    _constructRaftR2V();

    for (auto const& tri: pTris) {
        if (tri == nullptr or !tri->getInHost()) {
            continue;
        }
        for (auto spec_lidx: solver::spec_local_id::range(tri->patchdef()->countSpecs())) {
            _regTriPoolSync(tri->idx(),
                            tri->patchdef()->specL2G(spec_lidx),
                            tri->pools()[spec_lidx]);
        }
    }

    _syncPools(RDEF_TO_VESRAFT);
    // <<<<<<<<<<<<<< COMMUNICATION, RDEF -> VESRAFT

    // >>>>>>>>>>>>> COMMUNICATION, VESRAFT -> RDEF
    // Now need to recieve data from master (after it has done its business)
    _useVesV2R();
    _useRaftV2R();
    _syncPools(VESRAFT_TO_RDEF);
    // <<<<<<<<<<<<<<  COMMUNICATION, VESRAFT -> RDEF

    _updateLocal();  // for good measure
    recomputeUpdPeriod = true;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_runWithoutEField(double endtime) {
    MPI_Barrier(RDEFComm);

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
    while ((static_cast<int>(statedef().time() < endtime) & static_cast<int>(!aligned)) != 0) {
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
            double a0 = pA0;
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

        // Track how many diffusion 'steps' we do, simply for bookkeeping
        uint nsteps = 0;

        // to reduce memory cost we use directions to retrieve the update list in
        // upd process
        std::vector<KProc*> applied_diffs;
        std::vector<int> directions;

        for (size_t d_it = 0; d_it < diffSep; d_it++) {
            Diff* d = pDiffs[d_it];
            double rate = d->crData.rate;

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

            if (t1 >= 1.0) {  // this automatically allows for the minimum update period.
                t1 = 1.0;
            }

            // Calculate the occupancy, that is the integrated molecules over the
            // period (units s)
            double occupancy = d->getTet()->getPoolOccupancy(d->getLigLidx()) +
                               population *
                                   (update_period - d->getTet()->getLastUpdate(d->getLigLidx()));

            // n is, correctly, a double, but the binomial function requires
            // rounding to an integer.

            // occupancy/update_period gives the mean number of molecules during the
            // period
            double n_double = occupancy / update_period;

            // could be higher than those available with minimum update period
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

            if (mean_n == 0) {
                continue;
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
        for (size_t sd_it = 0; sd_it < sdiffSep; sd_it++) {
            SDiff* d = pSDiffs[sd_it];
            double rate = d->crData.rate;

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

            // n is, correctly, a double, but the binomial function requires
            // rounding to an integer.

            // occupancy/update_period gives the mean number of molecules during the
            // period
            double n_double = occupancy / update_period;

            // could be higher than those available with minimum update period
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

            if (mean_n == 0) {
                continue;
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

        _remoteSyncAndUpdate(requests, applied_diffs, directions);

        // *********************** Operator Split: SSA
        // *********************************

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
    }
    if (requests != nullptr) {
        MPI_Waitall(nNeighbHosts, requests, MPI_STATUSES_IGNORE);
        delete[] requests;
    }
    MPI_Barrier(RDEFComm);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_refreshEFTrisV() {
    for (auto tlidx: triangle_local_id::range(pEFNTris)) {
        EFTrisV[tlidx] = pEField->getTriV(tlidx);
    }
    pEFTrisVStale = false;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_runWithEField(double endtime) {
    while (statedef().time() < endtime) {
        double t0 = statedef().time();
        double t1 = std::min(t0 + pEFDT, endtime);
        if (util::almost_equal(t1, endtime)) {
            t1 = endtime;
        }
        _runWithoutEField(t1);

        // update host-local currents
        int i_begin = EFTrisI_offset[myRank_RDEF];
        int i_end = i_begin + EFTrisI_count[myRank_RDEF];

        double sttime = statedef().time();
        double real_ef_dt = sttime - t0;
        for (int i = i_begin; i < i_end; ++i) {
            auto tlidx = EFTrisI_idx[i];
            EFTrisI_permuted[i] =
                pEFTris_vec[tlidx]->computeI(EFTrisV[tlidx], real_ef_dt, sttime, getEfieldDT());
        }

        MPI_Allgatherv(MPI_IN_PLACE,
                       0,
                       MPI_DATATYPE_NULL,
                       EFTrisI_permuted.data(),
                       EFTrisI_count.data(),
                       EFTrisI_offset.data(),
                       MPI_DOUBLE,
                       RDEFComm);

        for (uint i = 0; i < pEFNTris; i++) {
            pEField->setTriI(EFTrisI_idx[i], EFTrisI_permuted[i]);
        }

        pEField->advance(real_ef_dt);
        _refreshEFTrisV();

        _updateLocal(pVdepKProcs);
    }
    MPI_Barrier(RDEFComm);
}

////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::advance(double adv) {
    if (adv < 0.0) {
        std::ostringstream os;
        os << "Time to advance cannot be negative";
        ArgErrLog(os.str());
    }

    double endtime = statedef().time() + adv;
    run(endtime);
}

////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::step() {
    NotImplErrLog("This function is not implemented.");
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getTime() const {
    return statedef().time();
}

////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getA0() const {
    return MPI_ConditionalReduce<double>(pA0, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::getNSteps() const {
    return statedef().nsteps();
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::setTime(double time) {
    statedef().setTime(time);
}

////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::setNSteps(uint nsteps) {
    statedef().setNSteps(nsteps);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getVesicleDT() const {
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::setVesicleDT(double /*dt*/) {}

////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::setTemp(double t) {
    if (_efflag() == false) {
        CLOG(WARNING, "general_log") << "Temperature set in simulation without membrane potential "
                                        "calculation will be ignored.\n";
    }
    AssertLog(t >= 0.0);
    pTemp = t;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getCompVol(solver::comp_global_id cidx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    CompRDEF* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    return comp->vol();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getCompSpecCount(solver::comp_global_id cidx,
                                         solver::spec_global_id sidx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    CompRDEF* comp = _comp(cidx);
    AssertLog(comp != nullptr);

    const solver::spec_local_id slidx = _specG2L_or_throw(comp, sidx);

    uint local_count = 0;
    for (auto const& tet: comp->tets()) {
        if (tet->getInHost()) {
            local_count += tet->pools()[slidx];
        }
    }

    return MPI_ConditionalReduce<uint>(local_count, MPI_UNSIGNED, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setCompSpecCount(solver::comp_global_id cidx,
                                       solver::spec_global_id sidx,
                                       double n) {
    MPI_Barrier(MPI_COMM_WORLD);
    AssertLog(cidx < statedef().countComps());
    AssertLog(n >= 0.0);
    AssertLog(statedef().countComps() == pComps.size());

    CompRDEF* comp = _comp(cidx);
    const solver::spec_local_id slidx = _specG2L_or_throw(comp, sidx);

    std::vector<uint> counts(comp->countTets());

    if (myRank_RDEF == RDEFmasterRank_RDEF) {
        // functions for distribution:
        auto set_count = [slidx](TetRDEF* tet, uint c) { tet->setCount(slidx, c); };
        auto inc_count = [slidx](TetRDEF* tet) { tet->incCount(slidx, 1, 0.0, true); };
        auto weight = [](const TetRDEFPVecCI& tet) { return (*tet)->vol(); };

        util::distribute_quantity<uint>(
            n, comp->tets().begin(), comp->tets().end(), weight, set_count, inc_count, *rng());
        std::transform(comp->tets().begin(),
                       comp->tets().end(),
                       counts.begin(),
                       [slidx](const TetRDEFP& t) { return t->pools()[slidx]; });
    }

    if (nHosts_RDEF != 1) {
        MPI_Bcast(counts.data(), counts.size(), MPI_UNSIGNED, RDEFmasterRank_RDEF, RDEFComm);
    }

    uint curr_pos = 0;
    for (auto const& t: comp->tets()) {
        if (myRank_RDEF != RDEFmasterRank_RDEF) {
            t->setCount(slidx, counts[curr_pos]);
        }
        _updateSpec(t, sidx);
        curr_pos++;
    }
    _updateSum();
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getCompSpecAmount(solver::comp_global_id cidx,
                                          solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getCompSpecCount(cidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setCompSpecAmount(solver::comp_global_id cidx,
                                        solver::spec_global_id sidx,
                                        double a) {
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setCompSpecCount(cidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getCompSpecConc(solver::comp_global_id cidx,
                                        solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    return _getCompSpecCount(cidx, sidx) / (1.0e3 * _comp(cidx)->vol() * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setCompSpecConc(solver::comp_global_id cidx,
                                      solver::spec_global_id sidx,
                                      double c) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    AssertLog(c >= 0.0);
    CompRDEF* comp = _comp(cidx);
    AssertLog(comp != nullptr);

    double count = c * (1.0e3 * comp->vol() * math::AVOGADRO);
    // the following method does argument checking on sidx
    _setCompSpecCount(cidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleRDEF::_getCompSpecClamped(solver::comp_global_id cidx,
                                         solver::spec_global_id sidx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());

    CompRDEF* comp = _comp(cidx);
    AssertLog(comp != nullptr);

    solver::spec_local_id lsidx = _specG2L_or_throw(comp, sidx);

    bool local_clamped = true;
    for (auto const& t: comp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        if (!t->clamped(lsidx)) {
            local_clamped = false;
        }
    }

    return MPI_ConditionalReduce<bool>(local_clamped, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setCompSpecClamped(solver::comp_global_id cidx,
                                         solver::spec_global_id sidx,
                                         bool b) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    CompRDEF* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::spec_local_id lsidx = _specG2L_or_throw(comp, sidx);

    // Set the flag in def object, though this may not be necessary
    comp->def()->setClamped(lsidx, b);

    for (auto const& t: comp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        t->setClamped(lsidx, b);
    }

    for (auto const& t: boundaryTets) {
        if (t->compdef() == comp->def()) {
            t->setClamped(lsidx, b);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getCompReacK(solver::comp_global_id cidx,
                                     solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    CompRDEF* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = _reacG2L_or_throw(comp, ridx);

    // We're just returning the default value for this comp, individual
    // tets may have different Kcsts set individually
    return comp->def()->kcst(lridx);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setCompReacK(solver::comp_global_id cidx,
                                   solver::reac_global_id ridx,
                                   double kf) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    AssertLog(kf >= 0.0);
    CompRDEF* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = _reacG2L_or_throw(comp, ridx);

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

bool TetVesicleRDEF::_getCompReacActive(solver::comp_global_id cidx,
                                        solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    CompRDEF* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = _reacG2L_or_throw(comp, ridx);

    bool local_active = true;

    for (auto const& t: comp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        if (t->reac(lridx).inactive()) {
            local_active = false;
        }
    }
    return MPI_ConditionalReduce<bool>(local_active, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setCompReacActive(solver::comp_global_id cidx,
                                        solver::reac_global_id ridx,
                                        bool a) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    CompRDEF* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = _reacG2L_or_throw(comp, ridx);

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

double TetVesicleRDEF::_getCompDiffD(solver::comp_global_id cidx,
                                     solver::diff_global_id didx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    CompRDEF* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::diff_local_id ldidx = _diffG2L_or_throw(comp, didx);

    // We're just returning the default value for this comp, individual
    // tets may have different Dcsts set individually
    return comp->def()->dcst(ldidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setCompDiffD(solver::comp_global_id cidx,
                                   solver::diff_global_id didx,
                                   double dk) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    AssertLog(dk >= 0.0);
    CompRDEF* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::diff_local_id ldidx = _diffG2L_or_throw(comp, didx);

    recomputeUpdPeriod = true;

    // First set the default value for the comp
    comp->def()->setDcst(ldidx, dk);

    // Now update all tets in this comp
    for (auto const& t: comp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        t->diff(ldidx).setDcst(dk);
    }

    // Rates have changed
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleRDEF::_getCompDiffActive(solver::comp_global_id cidx,
                                        solver::diff_global_id didx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    CompRDEF* comp = _comp(cidx);
    solver::diff_local_id ldidx = _diffG2L_or_throw(comp, didx);

    bool local_active = true;

    for (auto const& t: comp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        if (t->diff(ldidx).inactive()) {
            local_active = false;
        }
    }

    return MPI_ConditionalReduce<bool>(local_active, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setCompDiffActive(solver::comp_global_id cidx,
                                        solver::diff_global_id didx,
                                        bool act) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    CompRDEF* comp = _comp(cidx);
    AssertLog(comp != nullptr);
    solver::diff_local_id ldidx = _diffG2L_or_throw(comp, didx);

    for (auto const& t: comp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        t->diff(ldidx).setActive(act);
    }

    recomputeUpdPeriod = true;

    // It's cheaper to just recompute everything.
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

//////////////////// VESICULAR, Comp-specific ///////////////////////////

void TetVesicleRDEF::_setCompVesicleCount(solver::comp_global_id /*cidx*/,
                                          solver::vesicle_global_id /*vidx*/,
                                          uint /*n*/) {}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getCompVesicleCount(solver::comp_global_id cidx,
                                          solver::vesicle_global_id vidx) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTetVesicleDcst(tetrahedron_global_id tidx,
                                        solver::vesicle_global_id vidx,
                                        double dcst) {}

////////////////////////////////////////////////////////////////////////////////

solver::vesicle_individual_id TetVesicleRDEF::_addCompVesicle(solver::comp_global_id cidx,
                                                              solver::vesicle_global_id vidx) {
    return solver::vesicle_individual_id(MPI_ConditionalBcast<index_t>(
        0, MPI_STEPS_INDEX, vesraftRank_World, myRank_World, syncOutput, outputRank));
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getCompVesicleSurfaceSpecCount(solver::comp_global_id cidx,
                                                     solver::vesicle_global_id vidx,
                                                     solver::spec_global_id sidx) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getCompVesicleInnerSpecCount(solver::comp_global_id cidx,
                                                   solver::vesicle_global_id vidx,
                                                   solver::spec_global_id sidx) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::map<solver::vesicle_individual_id, uint> TetVesicleRDEF::_getCompVesicleSurfaceSpecCountMap(
    solver::comp_global_id cidx,
    solver::vesicle_global_id vidx,
    solver::spec_global_id sidx) const {
    std::map<solver::vesicle_individual_id, uint> return_map;
    MPI_ConditionalBcast<solver::vesicle_individual_id, uint>(
        return_map, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_map;
}

////////////////////////////////////////////////////////////////////////////////

std::map<solver::vesicle_individual_id, uint>
TetVesicleRDEF::_getCompVesicleSurfaceLinkSpecCountMap(solver::comp_global_id cidx,
                                                       solver::vesicle_global_id vidx,
                                                       solver::linkspec_global_id bsidx) const {
    std::map<solver::vesicle_individual_id, uint> return_map;
    MPI_ConditionalBcast<solver::vesicle_individual_id, uint>(
        return_map, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_map;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getCompVesicleSurfaceLinkSpecCount(solver::comp_global_id cidx,
                                                         solver::vesicle_global_id vidx,
                                                         solver::linkspec_global_id bsidx) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::vesicle_individual_id> TetVesicleRDEF::_getAllVesicleIndices() const {
    std::vector<solver::vesicle_individual_id> ves_indices;
    MPI_ConditionalBcast(ves_indices, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return ves_indices;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::vesicle_individual_id> TetVesicleRDEF::_getAllVesicleIndicesOnPath() const {
    std::vector<solver::vesicle_individual_id> ves_indices;
    MPI_ConditionalBcast(ves_indices, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return ves_indices;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::vesicle_individual_id> TetVesicleRDEF::_getCompVesicleIndices(
    solver::comp_global_id cidx,
    solver::vesicle_global_id vidx) const {
    std::vector<solver::vesicle_individual_id> return_vec;
    MPI_ConditionalBcast(return_vec, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

solver::comp_global_id TetVesicleRDEF::_getSingleVesicleCompartment(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index) const {
    return solver::comp_global_id(
        MPI_ConditionalBcast<index_t>(solver::comp_global_id::unknown_value(),
                                      MPI_STEPS_INDEX,
                                      vesraftRank_World,
                                      myRank_World,
                                      syncOutput,
                                      outputRank));
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setSingleVesicleInnerSpecCount(solver::vesicle_global_id vidx,
                                                     solver::vesicle_individual_id ves_unique_index,
                                                     solver::spec_global_id sidx,
                                                     uint c) {}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getSingleVesicleSurfaceSpecCount(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::spec_global_id sidx) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getSingleVesicleInnerSpecCount(solver::vesicle_global_id vidx,
                                                     solver::vesicle_individual_id ves_unique_index,
                                                     solver::spec_global_id sidx) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setSingleVesicleSurfaceSpecCount(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::spec_global_id sidx,
    uint c) {}


////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleRDEF::_getSingleSpecPosSpherical(
    solver::spec_global_id sidx,
    solver::pointspec_individual_id ps_unique_id) const {
    std::vector<double> return_vec;

    MPI_ConditionalBcast(
        return_vec, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);

    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_deleteSingleVesicle(solver::vesicle_global_id vidx,
                                          solver::vesicle_individual_id ves_unique_index) {}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getSingleVesicleSurfaceLinkSpecCount(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::linkspec_global_id lsidx) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::linkspec_individual_id> TetVesicleRDEF::_getSingleVesicleSurfaceLinkSpecIndices(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::linkspec_global_id lsidx) const {
    std::vector<solver::linkspec_individual_id> return_vec;
    MPI_ConditionalBcast(return_vec, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::pointspec_individual_id> TetVesicleRDEF::_getSingleVesicleSurfaceSpecIndices(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::spec_global_id sidx) const {
    std::vector<solver::pointspec_individual_id> return_vec;

    MPI_ConditionalBcast(return_vec, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setSingleVesiclePos(solver::vesicle_global_id vidx,
                                          solver::vesicle_individual_id ves_unique_index,
                                          const std::vector<double>& pos,
                                          bool force) {}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleRDEF::_getSingleVesiclePos(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index) const {
    std::vector<double> return_vec;
    auto nentries = MPI_ConditionalBcast<std::size_t>(
        return_vec.size(), MPI_STD_SIZE_T, vesraftRank_World, myRank_World, syncOutput, outputRank);
    MPI_ConditionalBcast<double>(
        return_vec, nentries, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> TetVesicleRDEF::_getSingleVesicleSurfaceLinkSpecPos(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::linkspec_global_id bsidx) {
    std::vector<std::vector<double>> return_mv;
    MPI_ConditionalBcast<double>(
        return_mv, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_mv;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleRDEF::_getSingleLinkSpecPos(
    solver::linkspec_individual_id ls_unique_id) const {
    std::vector<double> pos_vec;
    auto nentries = MPI_ConditionalBcast<std::size_t>(
        pos_vec.size(), MPI_STD_SIZE_T, vesraftRank_World, myRank_World, syncOutput, outputRank);
    MPI_ConditionalBcast<double>(
        pos_vec, nentries, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return pos_vec;
}

////////////////////////////////////////////////////////////////////////////////

solver::linkspec_individual_id TetVesicleRDEF::_getSingleLinkSpecLinkedTo(
    solver::linkspec_individual_id ls_unique_id) const {
    solver::linkspec_individual_id ls_id;
    return MPI_ConditionalBcast(
        ls_id, MPI_STEPS_INDEX, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

solver::vesicle_individual_id TetVesicleRDEF::_getSingleLinkSpecVes(
    solver::linkspec_individual_id ls_unique_id) const {
    solver::vesicle_individual_id ves_id;
    return MPI_ConditionalBcast(
        ves_id, MPI_STEPS_INDEX, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> TetVesicleRDEF::_getSingleVesicleSurfaceSpecPos(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::spec_global_id sidx) {
    std::vector<std::vector<double>> return_mv;
    MPI_ConditionalBcast<>(
        return_mv, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_mv;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> TetVesicleRDEF::_getSingleVesicleSurfaceSpecPosSpherical(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::spec_global_id sidx) const {
    std::vector<std::vector<double>> return_vec;
    MPI_ConditionalBcast(
        return_vec, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setSingleVesicleSurfaceSpecPosSpherical(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::spec_global_id sidx,
    const std::vector<std::vector<double>>& pos_spherical) {}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getSingleVesicleImmobility(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index) const {
    return MPI_ConditionalBcast<uint>(
        0u, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::tetrahedron_global_id> TetVesicleRDEF::_getSingleVesicleOverlapTets(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index) const {
    std::vector<steps::tetrahedron_global_id> return_vec;
    MPI_ConditionalBcast(return_vec, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getPatchArea(solver::patch_global_id pidx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    PatchRDEF* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    return patch->area();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getPatchSpecCount(solver::patch_global_id pidx,
                                          solver::spec_global_id sidx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    PatchRDEF* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::spec_local_id slidx = _specG2L_or_throw(patch, sidx);

    uint local_count = 0;
    for (auto const& t: patch->tris()) {
        if (t->getInHost()) {
            local_count += t->pools()[slidx];
        }
    }

    return MPI_ConditionalReduce<uint>(local_count, MPI_UNSIGNED, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setPatchSpecCount(solver::patch_global_id pidx,
                                        solver::spec_global_id sidx,
                                        double n) {
    MPI_Barrier(MPI_COMM_WORLD);
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    AssertLog(n >= 0.0);
    PatchRDEF* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::spec_local_id slidx = _specG2L_or_throw(patch, sidx);

    std::vector<uint> counts(patch->countTris());

    if (myRank_RDEF == RDEFmasterRank_RDEF) {
        // functions for distribution:
        auto set_count = [slidx](TriRDEF* tri, uint c) { tri->setCount(slidx, c); };
        auto inc_count = [slidx](TriRDEF* tri) { tri->incCount(slidx, 1, 0.0, true); };
        auto weight = [](const TriRDEFPVecCI& tri) { return (*tri)->area(); };

        util::distribute_quantity<uint>(n,
                                        patch->tris().begin(),
                                        patch->tris().end(),
                                        weight,
                                        set_count,
                                        inc_count,
                                        *rng(),
                                        patch->def()->area());
        std::transform(patch->tris().begin(),
                       patch->tris().end(),
                       counts.begin(),
                       [slidx](TriRDEF* const& t) { return t->pools()[slidx]; });
    }
    if (nHosts_RDEF != 1) {
        MPI_Bcast(counts.data(), counts.size(), MPI_UNSIGNED, RDEFmasterRank_RDEF, RDEFComm);
    }
    uint curr_pos = 0;
    for (auto const& t: patch->tris()) {
        if (myRank_RDEF != RDEFmasterRank_RDEF) {
            // set count only don't need sync
            t->setCount(slidx, counts[curr_pos]);
        }
        _updateSpec(t, sidx);
        curr_pos++;
    }

    _updateSum();
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getPatchSpecAmount(solver::patch_global_id pidx,
                                           solver::spec_global_id sidx) const {
    // the following method does all the necessary argument checking
    double count = _getPatchSpecCount(pidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setPatchSpecAmount(solver::patch_global_id pidx,
                                         solver::spec_global_id sidx,
                                         double a) {
    AssertLog(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setPatchSpecCount(pidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleRDEF::_getPatchSpecClamped(solver::patch_global_id pidx,
                                          solver::spec_global_id sidx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    PatchRDEF* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::spec_local_id lsidx = _specG2L_or_throw(patch, sidx);

    bool local_clamped = true;

    for (auto const& t: patch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        if (t->clamped(lsidx) == false) {
            local_clamped = false;
        }
    }
    return MPI_ConditionalReduce<bool>(local_clamped, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setPatchSpecClamped(solver::patch_global_id pidx,
                                          solver::spec_global_id sidx,
                                          bool buf) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    PatchRDEF* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::spec_local_id lsidx = _specG2L_or_throw(patch, sidx);

    // Set the flag in def object for consistency, though this is not
    // entirely necessary
    patch->def()->setClamped(lsidx, buf);


    for (auto const& t: patch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        t->setClamped(lsidx, buf);
    }

    for (auto const& t: boundaryTris) {
        if (t->patchdef() == patch->def()) {
            t->setClamped(lsidx, buf);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getPatchSReacK(solver::patch_global_id pidx,
                                       solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    PatchRDEF* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lsridx = _sreacG2L_or_throw(patch, ridx);

    // We're just returning the default value for this patch, individual
    // triangles may have different Kcsts set
    return patch->def()->kcst(lsridx);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setPatchSReacK(solver::patch_global_id pidx,
                                     solver::sreac_global_id ridx,
                                     double kf) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    AssertLog(kf >= 0.0);
    PatchRDEF* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lsridx = _sreacG2L_or_throw(patch, ridx);

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

bool TetVesicleRDEF::_getPatchSReacActive(solver::patch_global_id pidx,
                                          solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    PatchRDEF* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lsridx = _sreacG2L_or_throw(patch, ridx);

    bool local_active = true;
    for (auto const& t: patch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        if (t->sreac(lsridx).inactive()) {
            local_active = false;
        }
    }
    return MPI_ConditionalReduce<bool>(local_active, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setDiffBoundarySpecDiffusionActive(solver::diffboundary_global_id dbidx,
                                                         solver::spec_global_id sidx,
                                                         bool act) {
    AssertLog(dbidx < statedef().countDiffBoundaries());

    // Need to do two things:
    // 1) check if the species is defined in both compartments conencted
    // by the diffusion boundary
    // 2) loop over all tetrahedrons around the diff boundary and then the
    // diffusion rules and activate diffusion if the diffusion rule
    // relates to this species

    DiffBoundary* diffb = _diffboundary(dbidx);
    CompRDEF* compA = diffb->compA();
    CompRDEF* compB = diffb->compB();

    solver::spec_local_id lsidxA = _specG2L_or_throw(compA, sidx);
    solver::spec_local_id lsidxB = _specG2L_or_throw(compB, sidx);

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
    uint ntets = bdtets.size();

    for (uint bdt = 0; bdt != ntets; ++bdt) {
        TetRDEF* tet = _tet(bdtets[bdt]);
        if (!tet->getInHost()) {
            continue;
        }
        uint direction = bdtetsdir[bdt];
        AssertLog(direction < 4);

        // Each diff kproc then has access to the species through it's defined
        // parent
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

bool TetVesicleRDEF::_getDiffBoundarySpecDiffusionActive(solver::diffboundary_global_id dbidx,
                                                         solver::spec_global_id sidx) const {
    AssertLog(dbidx < statedef().countDiffBoundaries());

    DiffBoundary* diffb = _diffboundary(dbidx);
    CompRDEF* compA = diffb->compA();
    CompRDEF* compB = diffb->compB();

    solver::spec_local_id lsidxA = _specG2L_or_throw(compA, sidx);
    solver::spec_local_id lsidxB = _specG2L_or_throw(compB, sidx);

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

    uint ntets = bdtets.size();

    for (uint bdt = 0u; bdt != ntets; ++bdt) {
        TetRDEF* tet = _tet(bdtets[bdt]);
        if (!tet->getInHost()) {
            continue;
        }
        uint direction = bdtetsdir[bdt];
        AssertLog(direction < 4);

        // Each diff kproc then has access to the species through it's defined
        // parent
        uint ndiffs = tet->compdef()->countDiffs();
        for (auto d: solver::diff_local_id::range(ndiffs)) {
            Diff& diff = tet->diff(d);
            // sidx is the global species index; so is the lig() return from diffdef
            solver::spec_global_id specgidx = diff.def()->lig();
            if (specgidx == sidx) {
                local_active = diff.getDiffBndActive(direction);
                break;
            }
        }
    }
    return MPI_ConditionalReduce<bool>(local_active, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setDiffBoundarySpecDcst(solver::diffboundary_global_id dbidx,
                                              solver::spec_global_id sidx,
                                              double dcst,
                                              solver::comp_global_id direction_comp) {
    AssertLog(dbidx < statedef().countDiffBoundaries());

    DiffBoundary* diffb = _diffboundary(dbidx);
    CompRDEF* compA = diffb->compA();
    CompRDEF* compB = diffb->compB();

    solver::spec_local_id lsidxA = _specG2L_or_throw(compA, sidx);
    solver::spec_local_id lsidxB = _specG2L_or_throw(compB, sidx);

    if (lsidxA.unknown() or lsidxB.unknown()) {
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

    uint ntets = bdtets.size();

    for (uint bdt = 0u; bdt != ntets; ++bdt) {
        TetRDEF* tet = _tet(bdtets[bdt]);
        if (!tet->getInHost()) {
            continue;
        }
        // if tet compdef equals to dirc_compdef,
        // it is the desination tet so diff should not be changed
        // NULL (bidirection) and source tet are both different
        // fromdirc_compdef
        if (dirc_compdef == tet->compdef()) {
            continue;
        }
        uint direction = bdtetsdir[bdt];
        AssertLog(direction < 4);

        // Each diff kproc then has access to the species through it's defined
        // parent
        uint ndiffs = tet->compdef()->countDiffs();
        for (auto d: solver::diff_local_id::range(ndiffs)) {
            Diff& diff = tet->diff(d);
            // sidx is the global species index; so is the lig() return from diffdef
            solver::spec_global_id specgidx = diff.def()->lig();
            if (specgidx == sidx) {
#ifdef DIRECTIONAL_DCST_DEBUG
                CLOG(DEBUG, "steps_debug")
                    << "direction: " << direction << " dcst: " << dcst << "\n";
#endif

                diff.setDiffBndActive(direction, true);
                diff.setDirectionDcst(direction, dcst);
                _updateElement(&diff);
            }
        }
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setPatchSReacActive(solver::patch_global_id pidx,
                                          solver::sreac_global_id ridx,
                                          bool a) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    PatchRDEF* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lsridx = _sreacG2L_or_throw(patch, ridx);

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

bool TetVesicleRDEF::_getPatchVDepSReacActive(solver::patch_global_id pidx,
                                              solver::vdepsreac_global_id vsridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    PatchRDEF* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::vdepsreac_local_id lvsridx = _vdepsreacG2L_or_throw(patch, vsridx);

    bool local_active = true;

    for (auto const& t: patch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        if (t->vdepsreac(lvsridx).inactive()) {
            local_active = false;
        }
    }
    return MPI_ConditionalReduce<bool>(local_active, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setPatchVDepSReacActive(solver::patch_global_id pidx,
                                              solver::vdepsreac_global_id vsridx,
                                              bool a) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    PatchRDEF* patch = _patch(pidx);
    AssertLog(patch != nullptr);
    solver::vdepsreac_local_id lvsridx = _vdepsreacG2L_or_throw(patch, vsridx);

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

void TetVesicleRDEF::_setSDiffBoundarySpecDiffusionActive(solver::sdiffboundary_global_id sdbidx,
                                                          solver::spec_global_id sidx,
                                                          bool act) {
    AssertLog(sdbidx < statedef().countSDiffBoundaries());

    // Need to do two things:
    // 1) check if the species is defined in both patches connected
    // by the surface diffusion boundary
    // 2) loop over all triangles around the sdiff boundary and then the
    // diffusion rules and activate diffusion if the diffusion rule
    // relates to this species

    SDiffBoundary* sdiffb = _sdiffboundary(sdbidx);
    PatchRDEF* patchA = sdiffb->patchA();
    PatchRDEF* patchB = sdiffb->patchB();

    solver::spec_local_id lsidxA = patchA->def()->specG2L(sidx);
    solver::spec_local_id lsidxB = patchB->def()->specG2L(sidx);

    if (lsidxA.unknown() or lsidxB.unknown()) {
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
        TriRDEF* tri = _tri(sbdtris[sbdt]);
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
            solver::spec_global_id specgidx = sdiff.sdef()->lig();
            if (specgidx == sidx) {
                sdiff.setSDiffBndActive(direction, act);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleRDEF::_getSDiffBoundarySpecDiffusionActive(solver::sdiffboundary_global_id sdbidx,
                                                          solver::spec_global_id sidx) const {
    AssertLog(sdbidx < statedef().countSDiffBoundaries());

    SDiffBoundary* sdiffb = _sdiffboundary(sdbidx);
    PatchRDEF* patchA = sdiffb->patchA();
    PatchRDEF* patchB = sdiffb->patchB();

    solver::spec_local_id lsidxA = patchA->def()->specG2L(sidx);
    solver::spec_local_id lsidxB = patchB->def()->specG2L(sidx);

    if (lsidxA.unknown() or lsidxB.unknown()) {
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
        TriRDEF* tri = _tri(sbdtris[sbdt]);
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
            solver::spec_global_id specgidx = sdiff.sdef()->lig();
            if (specgidx == sidx) {
                local_active = static_cast<short>(sdiff.getSDiffBndActive(direction));
                break;
            }
        }
    }

    return MPI_ConditionalReduce<bool>(
        local_active != 0, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setSDiffBoundarySpecDcst(solver::sdiffboundary_global_id sdbidx,
                                               solver::spec_global_id sidx,
                                               double dcst,
                                               solver::patch_global_id direction_patch) {
    SDiffBoundary* sdiffb = _sdiffboundary(sdbidx);
    PatchRDEF* patchA = sdiffb->patchA();
    PatchRDEF* patchB = sdiffb->patchB();

    solver::spec_local_id lsidxA = patchA->def()->specG2L(sidx);
    solver::spec_local_id lsidxB = patchB->def()->specG2L(sidx);

    if (lsidxA.unknown() or lsidxB.unknown()) {
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
        TriRDEF* tri = _tri(sbdtris[sbdt]);

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
            solver::spec_global_id specgidx = sdiff.sdef()->lig();
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
///////////////////// RAFTS, ENDOCYTOSIS, Patch-specific  //////////////////////

void TetVesicleRDEF::_setPatchRaftCount(solver::patch_global_id pidx,
                                        solver::raft_global_id ridx,
                                        uint n) {}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getPatchRaftCount(solver::patch_global_id pidx,
                                        solver::raft_global_id ridx) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getSingleRaftImmobility(
    solver::raft_global_id ridx,
    solver::raft_individual_id raft_individual_index) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::raft_individual_id> TetVesicleRDEF::_getPatchRaftIndices(
    solver::patch_global_id pidx,
    solver::raft_global_id ridx) const {
    std::vector<solver::raft_individual_id> return_vec;
    MPI_ConditionalBcast(return_vec, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

solver::patch_global_id TetVesicleRDEF::_getSingleRaftPatch(
    solver::raft_global_id ridx,
    solver::raft_individual_id raft_unique_index) const {
    return solver::patch_global_id(
        MPI_ConditionalBcast<index_t>(solver::patch_global_id::unknown_value(),
                                      MPI_STEPS_INDEX,
                                      vesraftRank_World,
                                      myRank_World,
                                      syncOutput,
                                      outputRank));
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleRDEF::_getSingleRaftPos(
    solver::raft_global_id ridx,
    solver::raft_individual_id raft_individual_index) const {
    std::vector<double> return_vec;
    auto nentries = MPI_ConditionalBcast<std::size_t>(
        return_vec.size(), MPI_STD_SIZE_T, vesraftRank_World, myRank_World, syncOutput, outputRank);
    MPI_ConditionalBcast<double>(
        return_vec, nentries, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

std::map<solver::raft_individual_id, uint> TetVesicleRDEF::_getPatchRaftSpecCountMap(
    solver::patch_global_id pidx,
    solver::raft_global_id ridx,
    solver::spec_global_id sidx) const {
    std::map<solver::raft_individual_id, uint> return_map;
    MPI_ConditionalBcast<solver::raft_individual_id, uint>(
        return_map, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_map;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getPatchRaftSpecCount(solver::patch_global_id pidx,
                                            solver::raft_global_id ridx,
                                            solver::spec_global_id sidx) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getSingleRaftSpecCount(solver::raft_global_id ridx,
                                             solver::raft_individual_id raft_individual_index,
                                             solver::spec_global_id sidx) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setSingleRaftSpecCount(solver::raft_global_id ridx,
                                             solver::raft_individual_id raft_individual_index,
                                             solver::spec_global_id sidx,
                                             uint c) {}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setPatchEndocyticZoneEndocytosisActive(solver::patch_global_id pidx,
                                                             std::string const& zone,
                                                             solver::endocytosis_global_id endogidx,
                                                             bool active) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setPatchEndocyticZoneEndocytosisK(solver::patch_global_id pidx,
                                                        std::string const& zone,
                                                        solver::endocytosis_global_id endogidx,
                                                        double k) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getPatchEndocyticZoneEndocytosisExtent(
    solver::patch_global_id pidx,
    std::string const& zone,
    solver::endocytosis_global_id endogidx) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::EndocytosisEvent> TetVesicleRDEF::_getPatchEndocyticZoneEndocytosisEvents(
    solver::patch_global_id pidx,
    std::string const& zone,
    solver::endocytosis_global_id endogidx) const {
    std::vector<solver::EndocytosisEvent> events;
    MPI_ConditionalBcast(events,
                         dataTypeUtil.MPI_ExocytosisEventSync,
                         vesraftRank_World,
                         myRank_World,
                         syncOutput,
                         outputRank);
    return events;
}

////////////////////////////////////////////////////////////////////////////////

solver::kproc_global_id TetVesicleRDEF::addKProc_(KProc* kp, bool Vdep) {
    // AssertLog (kp != nullptr); removed because nullptr can now be used

    SchedIDX nidx(static_cast<uint>(pKProcs.size()));  // because pKProcs.size() is ulong
    pKProcs.push_back(kp);
    if (Vdep) {
        pVdepKProcs.push_back(kp);
    }
    return nidx;  // calling func now has the responsiblity of calling
                  // kproc::setSchedIDX
}

////////////////////////////////////////////////////////////////////////////////

KProc* TetVesicleRDEF::_getNext() const {
    AssertLog(pA0 >= 0.0);
    // Quick check to see whether nothing is there.
    if (pA0 == 0.0) {
        return nullptr;
    }

    double selector = pA0 * rng()->getUnfII();

    double partial_sum = 0.0;

    uint n_neg_groups = nGroups.size();
    uint n_pos_groups = pGroups.size();

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

void TetVesicleRDEF::_executeStep(KProc* kp, double dt, double period) {
    kp->apply(rng(), dt, statedef().time(), period, this);
    statedef().incTime(dt);

    std::vector<KProc*> const& upd = kp->getLocalUpdVec();
    _updateLocal(upd);

    statedef().incNSteps(1);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_updateSpec(TetRDEF* tet, solver::spec_global_id spec_gidx) {
    // NOTE: this function does not update the Sum of popensity, _updateSum() is
    // required after calling it.
    AssertLog(_getTetSpecDefined(tet->idx(), spec_gidx));

    if (!tet->getInHost()) {
        return;
    }
    KProcPSet updset;

    // Loop over tet.
    uint nkprocs = tet->countKProcs();

    for (uint k = 0; k < nkprocs; k++) {
        if (tet->KProcDepSpecTet(k, tet, spec_gidx)) {
            updset.insert(tet->getKProc(k));
        }
    }

    for (auto const& tri: tet->nexttris()) {
        if (tri == nullptr) {
            continue;
        }
        nkprocs = tri->countKProcs();
        for (uint sk = 0; sk < nkprocs; sk++) {
            if (tri->KProcDepSpecTet(sk, tet, spec_gidx)) {
                updset.insert(tri->getKProc(sk));
            }
        }
    }

    /* TetOpSplit uses this- is it necessary?
    for (auto & kp : updset) {
        _updateElement(kp);
    }
    */

    _updateLocal(updset);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_updateSpec(TriRDEF* tri, solver::spec_global_id spec_gidx) {
    // NOTE: this function does not update the Sum of popensity, _updateSum() is
    // required after calling it.
    AssertLog(_getTriSpecDefined(tri->idx(), spec_gidx));

    if (!tri->getInHost()) {
        return;
    }
    KProcPSet updset;

    uint nkprocs = tri->countKProcs();

    for (uint sk = 0; sk < nkprocs; sk++) {
        if (tri->KProcDepSpecTri(sk, tri, spec_gidx)) {
            updset.insert(tri->getKProc(sk));
        }
    }

    /* TetOpSplit uses this- is it necessary?
    for (auto & kp : updset) {
        _updateElement(kp);
    }
    */

    _updateLocal(updset);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getCompReacH(solver::comp_global_id cidx,
                                     solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    CompRDEF* lcomp = _comp(cidx);
    AssertLog(lcomp != nullptr);
    solver::reac_local_id lridx = _reacG2L_or_throw(lcomp, ridx);

    double local_h = 0.0;
    for (auto t: lcomp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        Reac& reac = t->reac(lridx);
        local_h += reac.h();
    }
    return MPI_ConditionalReduce<double>(local_h, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getCompReacC(solver::comp_global_id cidx,
                                     solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    CompRDEF* lcomp = _comp(cidx);
    AssertLog(lcomp != nullptr);
    solver::reac_local_id lridx = _reacG2L_or_throw(lcomp, ridx);

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
    auto global_c =
        MPI_ConditionalReduce<double>(local_c, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
    auto global_v =
        MPI_ConditionalReduce<double>(local_v, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
    return global_c / global_v;
}

////////////////////////////////////////////////////////////////////////////////

long double TetVesicleRDEF::_getCompReacA(solver::comp_global_id cidx,
                                          solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    CompRDEF* lcomp = _comp(cidx);
    AssertLog(lcomp != nullptr);
    solver::reac_local_id lridx = _reacG2L_or_throw(lcomp, ridx);

    long double local_a = 0.0L;
    for (auto t: lcomp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        Reac& reac = t->reac(lridx);
        local_a += static_cast<long double>(reac.rate());
    }
    return MPI_ConditionalReduce<long double>(
        local_a, MPI_LONG_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetVesicleRDEF::_getCompReacExtent(solver::comp_global_id cidx,
                                                      solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    CompRDEF* lcomp = _comp(cidx);
    AssertLog(lcomp != nullptr);
    solver::reac_local_id lridx = _reacG2L_or_throw(lcomp, ridx);

    unsigned long long local_x = 0;
    for (auto t: lcomp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        Reac& reac = t->reac(lridx);
        local_x += reac.getExtent();
    }

    return MPI_ConditionalReduce<unsigned long long>(
        local_x, MPI_UNSIGNED_LONG_LONG, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_resetCompReacExtent(solver::comp_global_id cidx,
                                          solver::reac_global_id ridx) {
    AssertLog(cidx < statedef().countComps());
    CompRDEF* lcomp = _comp(cidx);
    AssertLog(lcomp != nullptr);
    solver::reac_local_id lridx = _reacG2L_or_throw(lcomp, ridx);

    for (auto t: lcomp->tets()) {
        if (!t->getInHost()) {
            continue;
        }
        Reac& reac = t->reac(lridx);
        reac.resetExtent();
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getPatchSReacH(solver::patch_global_id pidx,
                                       solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    PatchRDEF* lpatch = _patch(pidx);
    AssertLog(lpatch != nullptr);
    solver::sreac_local_id lsridx = _sreacG2L_or_throw(lpatch, ridx);

    double local_h = 0.0;
    for (auto t: lpatch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        SReac& sreac = t->sreac(lsridx);
        local_h += sreac.h();
    }

    return MPI_ConditionalReduce<double>(local_h, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getPatchSReacC(solver::patch_global_id pidx,
                                       solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    PatchRDEF* lpatch = _patch(pidx);
    AssertLog(lpatch != nullptr);
    solver::sreac_local_id lsridx = _sreacG2L_or_throw(lpatch, ridx);

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
    auto global_c =
        MPI_ConditionalReduce<double>(local_c, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
    auto global_a =
        MPI_ConditionalReduce<double>(local_a, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
    return global_c / global_a;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getPatchSReacA(solver::patch_global_id pidx,
                                       solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    PatchRDEF* lpatch = _patch(pidx);
    AssertLog(lpatch != nullptr);
    solver::sreac_local_id lsridx = _sreacG2L_or_throw(lpatch, ridx);

    double local_a = 0.0;
    for (auto t: lpatch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        SReac& sreac = t->sreac(lsridx);
        local_a += sreac.rate();
    }
    return MPI_ConditionalReduce<double>(local_a, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetVesicleRDEF::_getPatchSReacExtent(solver::patch_global_id pidx,
                                                        solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    PatchRDEF* lpatch = _patch(pidx);
    AssertLog(lpatch != nullptr);
    solver::sreac_local_id lsridx = _sreacG2L_or_throw(lpatch, ridx);

    unsigned long long local_x = 0;
    for (auto t: lpatch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        SReac& sreac = t->sreac(lsridx);
        local_x += sreac.getExtent();
    }
    return MPI_ConditionalReduce<unsigned long long>(
        local_x, MPI_UNSIGNED_LONG_LONG, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_resetPatchSReacExtent(solver::patch_global_id pidx,
                                            solver::sreac_global_id ridx) {
    AssertLog(pidx < statedef().countPatches());
    PatchRDEF* lpatch = _patch(pidx);
    AssertLog(lpatch != nullptr);
    solver::sreac_local_id lsridx = _sreacG2L_or_throw(lpatch, ridx);

    for (auto t: lpatch->tris()) {
        if (!t->getInHost()) {
            continue;
        }
        SReac& sreac = t->sreac(lsridx);
        sreac.resetExtent();
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTetVol(tetrahedron_global_id tidx) const {
    AssertLog(tidx < pTets.size());
    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
        ArgErrLog(os.str());
    }
    return pTets[tidx]->staticVol();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getTetReducedVol(tetrahedron_global_id tidx) const {
    ArgErrLogIf(tidx >= pTets.size(), "Tetrahedron index out of range.");

    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
        ArgErrLog(os.str());
    }
    return MPI_ConditionalBcast<double>(
        0, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTetVol(tetrahedron_global_id /*tidx*/, double /*vol*/) {
    std::ostringstream os;
    os << "Can not change tetrahedron volume in a mesh based solver.\n";
    NotImplErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleRDEF::_getTetSpecDefined(tetrahedron_global_id tidx,
                                        solver::spec_global_id sidx) const {
    AssertLog(tidx < pTets.size());
    AssertLog(sidx < statedef().countSpecs());

    TetRDEF* tet = _getTet(tidx);
    solver::spec_local_id lsidx = tet->compdef()->specG2L(sidx);
    return lsidx.valid();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTetSpecCount(tetrahedron_global_id tidx,
                                        solver::spec_global_id sidx) const {
    AssertLog(tidx < pTets.size());
    AssertLog(sidx < statedef().countSpecs());
    int host = _getTetHost(tidx);
    TetRDEF* tet = _getTet(tidx);
    solver::spec_local_id lsidx = tet->compdef()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    double count = tet->pools()[lsidx];
    return MPI_ConditionalBcast<double>(
        count, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTetSpecCount(tetrahedron_global_id tidx,
                                      solver::spec_global_id sidx,
                                      double n) {
    AssertLog(tidx < pTets.size());
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(n >= 0.0);

    _getTetHost(tidx);

    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }
    if (n > UINT_MAX) {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << UINT_MAX << ").\n";
        ArgErrLog(os.str());
    }

    TetRDEF* tet = pTets[tidx];

    if (tet->getInHost()) {
        solver::spec_local_id lsidx = tet->compdef()->specG2L(sidx);
        if (lsidx.unknown()) {
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
        tet->setCount(lsidx, c);
        _updateSpec(tet, sidx);
        _updateSum();
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTetSpecAmount(tetrahedron_global_id tidx,
                                         solver::spec_global_id sidx) const {
    // following method does all necessary argument checking
    double count = _getTetSpecCount(tidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTetSpecAmount(tetrahedron_global_id tidx,
                                       solver::spec_global_id sidx,
                                       double m) {
    // convert amount in mols to number of molecules
    double m2 = m * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setTetSpecCount(tidx, sidx, m2);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTetSpecConc(tetrahedron_global_id tidx,
                                       solver::spec_global_id sidx) const {
    // following method does all necessary argument checking
    double count = _getTetSpecCount(tidx, sidx);
    TetRDEF* tet = pTets[tidx];
    double vol = tet->staticVol();
    return count / (1.0e3 * vol * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTetSpecConc(tetrahedron_global_id tidx,
                                     solver::spec_global_id sidx,
                                     double c) {
    AssertLog(c >= 0.0);
    AssertLog(tidx < pTets.size());

    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
        ArgErrLog(os.str());
    }

    TetRDEF* tet = pTets[tidx];
    double count = c * (1.0e3 * tet->staticVol() * math::AVOGADRO);
    // the following method does all the necessary argument checking
    _setTetSpecCount(tidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleRDEF::_getTetSpecClamped(tetrahedron_global_id tidx,
                                        solver::spec_global_id sidx) const {
    // TODO
    // is clamped check stored in vesraft, or need to fetch frm rdef?
    AssertLog(tidx < pTets.size());
    AssertLog(sidx < statedef().countSpecs());

    TetRDEF* tet = _getTet(tidx);
    int host = _getTetHost(tidx);

    solver::spec_local_id lsidx = tet->compdef()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    bool clamped = tet->clamped(lsidx);
    return MPI_ConditionalBcast<bool>(
        clamped, MPI_C_BOOL, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTetSpecClamped(tetrahedron_global_id tidx,
                                        solver::spec_global_id sidx,
                                        bool buf) {
    AssertLog(tidx < pTets.size());
    AssertLog(sidx < statedef().countSpecs());

    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }

    TetRDEF* tet = pTets[tidx];

    solver::spec_local_id lsidx = tet->compdef()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    tet->setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTetReacK(tetrahedron_global_id tidx, solver::reac_global_id ridx) const {
    AssertLog(tidx < pTets.size());
    AssertLog(ridx < statedef().countReacs());

    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }
    int host = _getTetHost(tidx);

    TetRDEF* tet = pTets[tidx];

    solver::reac_local_id lridx = tet->compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    double kcst = 0;
    if (host == myRank_World) {
        kcst = tet->reac(lridx).kcst();
    }
    return MPI_ConditionalBcast<double>(
        kcst, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTetReacK(tetrahedron_global_id tidx,
                                  solver::reac_global_id ridx,
                                  double kf) {
    AssertLog(tidx < pTets.size());
    AssertLog(ridx < statedef().countReacs());
    AssertLog(kf >= 0.0);
    int host = _getTetHost(tidx);
    if (host != myRank_World) {
        return;
    }

    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }

    TetRDEF* tet = pTets[tidx];

    solver::reac_local_id lridx = tet->compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "\nReaction undefined in tetrahedron.";
        ArgErrLog(os.str());
    }

    tet->reac(lridx).setKcst(kf);

    _updateElement(&tet->reac(lridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleRDEF::_getTetReacActive(tetrahedron_global_id tidx,
                                       solver::reac_global_id ridx) const {
    AssertLog(tidx < pTets.size());
    AssertLog(ridx < statedef().countReacs());

    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }
    int host = _getTetHost(tidx);

    TetRDEF* tet = pTets[tidx];

    solver::reac_local_id lridx = tet->compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    bool active = false;
    if (host == myRank_World) {
        active = tet->reac(lridx).active();
    }
    return MPI_ConditionalBcast<bool>(
        active, MPI_C_BOOL, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTetReacActive(tetrahedron_global_id tidx,
                                       solver::reac_global_id ridx,
                                       bool act) {
    AssertLog(tidx < pTets.size());
    AssertLog(ridx < statedef().countReacs());

    int host = _getTetHost(tidx);
    if (host != myRank_World) {
        return;
    }

    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }

    TetRDEF* tet = pTets[tidx];

    solver::reac_local_id lridx = tet->compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    tet->reac(lridx).setActive(act);

    _updateElement(&tet->reac(lridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTetDiffD(tetrahedron_global_id tidx,
                                    solver::diff_global_id didx,
                                    tetrahedron_global_id direction_tet) const {
    AssertLog(tidx < pTets.size());
    AssertLog(didx < statedef().countDiffs());

    int host = _getTetHost(tidx);

    TetRDEF* tet = _getTet(tidx);

    solver::diff_local_id ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    double dcst = 0.0;

    if (host == myRank_World) {
        if (direction_tet.unknown()) {
            dcst = tet->diff(ldidx).dcst();
        } else {
            int direction = tet->getTetDirection(direction_tet);
            if (direction == -1) {
                std::ostringstream os;
                os << "Tetrahedron " << direction_tet << " is not a neighbor of tetrahedron "
                   << tidx << ".\n";
                ArgErrLog(os.str());
            }

            dcst = tet->diff(ldidx).dcst(direction);
        }
    }

    return MPI_ConditionalBcast<double>(
        dcst, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTetDiffD(tetrahedron_global_id tidx,
                                  solver::diff_global_id didx,
                                  double dk,
                                  tetrahedron_global_id direction_tet) {
#ifdef DIRECTIONAL_DCST_DEBUG
    CLOG_IF(direction_tet != std::numeric_limits<uint>::max(), DEBUG, "steps_debug")
        << "tidx: " << tidx << " didx: " << didx << " dk: " << dk
        << " direction tet: " << direction_tet << "\n";

    CLOG_IF(direction_tet == std::numeric_limits<uint>::max(), DEBUG, "steps_debug")
        << "tidx: " << tidx << " didx: " << didx << " dk: " << dk << " for all directions.\n";
#endif

    AssertLog(tidx < pTets.size());
    AssertLog(didx < statedef().countDiffs());

    int host = _getTetHost(tidx);

    recomputeUpdPeriod = true;
    if (host != myRank_World) {
        return;
    }
    TetRDEF* tet = _getTet(tidx);

    solver::diff_local_id ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    if (direction_tet.unknown()) {
        tet->diff(ldidx).setDcst(dk);
    } else {
        int direction = tet->getTetDirection(direction_tet);
        if (direction == -1) {
            std::ostringstream os;
            os << "Tetrahedron " << direction_tet << " is not a neighbor of tetrahedron " << tidx
               << ".\n";
            ArgErrLog(os.str());
        }

#ifdef DIRECTIONAL_DCST_DEBUG
        CLOG(DEBUG, "steps_debug")
            << "use tet " << direction_tet << " to set direction " << direction << ".\n";
#endif

        tet->diff(ldidx).setDirectionDcst(direction, dk);
    }
    _updateElement(&tet->diff(ldidx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleRDEF::_getTetDiffActive(tetrahedron_global_id tidx,
                                       solver::diff_global_id didx) const {
    AssertLog(tidx < pTets.size());
    AssertLog(didx < statedef().countDiffs());

    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }

    int host = _getTetHost(tidx);

    TetRDEF* tet = pTets[tidx];

    solver::diff_local_id ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    bool active = false;
    if (host == myRank_World) {
        active = tet->diff(ldidx).active();
    }
    return MPI_ConditionalBcast<bool>(
        active, MPI_C_BOOL, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTetDiffActive(tetrahedron_global_id tidx,
                                       solver::diff_global_id didx,
                                       bool act) {
    AssertLog(tidx < pTets.size());
    AssertLog(didx < statedef().countDiffs());

    int host = _getTetHost(tidx);

    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }

    TetRDEF* tet = pTets[tidx];

    solver::diff_local_id ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    recomputeUpdPeriod = true;
    if (host == myRank_World) {
        tet->diff(ldidx).setActive(act);

        _updateElement(&tet->diff(ldidx));
        _updateSum();
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTetReacH(tetrahedron_global_id tidx, solver::reac_global_id ridx) const {
    AssertLog(tidx < pTets.size());
    AssertLog(ridx < statedef().countReacs());

    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }

    int host = _getTetHost(tidx);

    TetRDEF* tet = pTets[tidx];

    solver::reac_local_id lridx = tet->compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    double h = 0;
    if (host == myRank_World) {
        h = tet->reac(lridx).h();
    }
    return MPI_ConditionalBcast<double>(h, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTetReacC(tetrahedron_global_id tidx, solver::reac_global_id ridx) const {
    AssertLog(tidx < pTets.size());
    AssertLog(ridx < statedef().countReacs());

    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }

    int host = _getTetHost(tidx);

    TetRDEF* tet = pTets[tidx];

    solver::reac_local_id lridx = tet->compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    double c = 0;
    if (host == myRank_World) {
        c = tet->reac(lridx).c();
    }
    return MPI_ConditionalBcast<double>(c, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTetReacA(tetrahedron_global_id tidx, solver::reac_global_id ridx) const {
    AssertLog(tidx < pTets.size());
    AssertLog(ridx < statedef().countReacs());

    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }

    int host = _getTetHost(tidx);

    TetRDEF* tet = pTets[tidx];

    solver::reac_local_id lridx = tet->compdef()->reacG2L(ridx);
    if (lridx.unknown()) {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    double a = 0;
    if (host == myRank_World) {
        a = tet->reac(lridx).rate();
    }
    return MPI_ConditionalBcast<double>(a, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTetDiffA(tetrahedron_global_id tidx, solver::diff_global_id didx) const {
    AssertLog(tidx < pTets.size());
    AssertLog(didx < statedef().countDiffs());

    if (pTets[tidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }

    int host = _getTetHost(tidx);

    TetRDEF* tet = pTets[tidx];

    solver::diff_local_id ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    double a = 0;
    if (host == myRank_World) {
        a = tet->diff(ldidx).rate();
    }
    return MPI_ConditionalBcast<double>(a, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriArea(triangle_global_id tidx) const {
    AssertLog(tidx < pTris.size());
    return _getTri(tidx)->area();
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTriArea(triangle_global_id /*tidx*/, double /*area*/) {
    NotImplErrLog("This function is not implemented.");
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleRDEF::_getTriSpecDefined(triangle_global_id tidx,
                                        solver::spec_global_id sidx) const {
    AssertLog(tidx < pTris.size());
    AssertLog(sidx < statedef().countSpecs());

    TriRDEF* tri = _getTri(tidx);
    solver::spec_local_id lsidx = tri->patchdef()->specG2L(sidx);
    return lsidx.valid();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriSpecCount(triangle_global_id tidx,
                                        solver::spec_global_id sidx) const {
    AssertLog(tidx < pTris.size());
    AssertLog(sidx < statedef().countSpecs());
    int host = _getTriHost(tidx);
    TriRDEF* tri = _getTri(tidx);
    solver::spec_local_id lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    uint count = tri->pools()[lsidx];
    const auto it = triHosts.find(tidx);
    if (it == triHosts.end()) {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a host.\n";
        ArgErrLog(os.str());
    }
    return MPI_ConditionalBcast<double>(
        count, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTriSpecCount(triangle_global_id tidx,
                                      solver::spec_global_id sidx,
                                      double n) {
    MPI_Barrier(RDEFComm);
    AssertLog(tidx < pTris.size());
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(n >= 0.0);

    if (n > UINT_MAX) {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << UINT_MAX << ").\n";
        ArgErrLog(os.str());
    }

    TriRDEF* tri = _getTri(tidx);

    if (tri->getInHost()) {
        solver::spec_local_id lsidx = tri->patchdef()->specG2L(sidx);
        if (lsidx.unknown()) {
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
        tri->setCount(lsidx, c);
        _updateSpec(tri, sidx);
        _updateSum();
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriSpecAmount(triangle_global_id tidx,
                                         solver::spec_global_id sidx) const {
    // following method does all necessary argument checking
    double count = _getTriSpecCount(tidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTriSpecAmount(triangle_global_id tidx,
                                       solver::spec_global_id sidx,
                                       double m) {
    // convert amount in mols to number of molecules
    double m2 = m * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setTriSpecCount(tidx, sidx, m2);
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleRDEF::_getTriSpecClamped(triangle_global_id tidx,
                                        solver::spec_global_id sidx) const {
    AssertLog(tidx < pTris.size());
    AssertLog(sidx < statedef().countSpecs());

    TriRDEF* tri = _getTri(tidx);
    int host = _getTriHost(tidx);
    solver::spec_local_id lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    bool clamped = tri->clamped(lsidx);
    return MPI_ConditionalBcast<bool>(
        clamped, MPI_C_BOOL, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTriSpecClamped(triangle_global_id tidx,
                                        solver::spec_global_id sidx,
                                        bool buf) {
    AssertLog(tidx < pTris.size());
    AssertLog(sidx < statedef().countSpecs());

    TriRDEF* tri = _getTri(tidx);

    solver::spec_local_id lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    tri->setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriSReacK(triangle_global_id tidx, solver::sreac_global_id ridx) const {
    AssertLog(tidx < pTris.size());
    AssertLog(ridx < statedef().countSReacs());
    int host = _getTriHost(tidx);
    TriRDEF* tri = _getTri(tidx);

    solver::sreac_local_id lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    double kcst = 0;
    if (host == myRank_World) {
        kcst = tri->sreac(lsridx).kcst();
    }
    return MPI_ConditionalBcast<double>(
        kcst, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTriSReacK(triangle_global_id tidx,
                                   solver::sreac_global_id ridx,
                                   double kf) {
    AssertLog(tidx < pTris.size());
    AssertLog(ridx < statedef().countSReacs());

    int host = _getTriHost(tidx);
    if (host != myRank_World) {
        return;
    }

    TriRDEF* tri = _getTri(tidx);

    solver::sreac_local_id lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    tri->sreac(lsridx).setKcst(kf);

    _updateElement(&tri->sreac(lsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleRDEF::_getTriSReacActive(triangle_global_id tidx,
                                        solver::sreac_global_id ridx) const {
    AssertLog(tidx < pTris.size());
    AssertLog(ridx < statedef().countSReacs());

    int host = _getTriHost(tidx);
    TriRDEF* tri = _getTri(tidx);

    solver::sreac_local_id lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    bool active = false;
    if (host == myRank_World) {
        active = tri->sreac(lsridx).active();
    }
    return MPI_ConditionalBcast<bool>(
        active, MPI_C_BOOL, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTriSReacActive(triangle_global_id tidx,
                                        solver::sreac_global_id ridx,
                                        bool act) {
    AssertLog(tidx < pTris.size());
    AssertLog(ridx < statedef().countSReacs());

    int host = _getTriHost(tidx);
    if (host != myRank_World) {
        return;
    }
    TriRDEF* tri = _getTri(tidx);

    solver::sreac_local_id lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    tri->sreac(lsridx).setActive(act);
    _updateElement(&tri->sreac(lsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriSDiffD(triangle_global_id tidx,
                                     solver::surfdiff_global_id didx,
                                     triangle_global_id direction_tri) const {
    AssertLog(tidx < pTris.size());
    AssertLog(didx < statedef().countSurfDiffs());

    int host = _getTriHost(tidx);

    TriRDEF* tri = _getTri(tidx);

    solver::surfdiff_local_id ldidx = tri->patchdef()->surfdiffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    double dcst = 0.0;

    if (host == myRank_World) {
        if (direction_tri.unknown()) {
            dcst = tri->sdiff(ldidx).dcst();
        } else {
            int direction = tri->getTriDirection(direction_tri);
            if (direction == -1) {
                std::ostringstream os;
                os << "Triangle " << direction_tri << " is not a neighbor of triangle " << tidx
                   << ".\n";
                ArgErrLog(os.str());
            }

            dcst = tri->sdiff(ldidx).dcst(direction);
        }
    }

    return MPI_ConditionalBcast<double>(
        dcst, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTriSDiffD(triangle_global_id tidx,
                                   solver::surfdiff_global_id didx,
                                   double dk,
                                   triangle_global_id direction_tri) {
    AssertLog(tidx < pTris.size());
    AssertLog(didx < statedef().countSurfDiffs());
    int host = _getTriHost(tidx);

    TriRDEF* tri = _getTri(tidx);

    solver::surfdiff_local_id ldidx = tri->patchdef()->surfdiffG2L(didx);
    if (ldidx.unknown()) {
        std::ostringstream os;
        os << "Diffusion rule undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    recomputeUpdPeriod = true;
    if (host != myRank_World) {
        return;
    }
    if (direction_tri.unknown()) {
        tri->sdiff(ldidx).setDcst(dk);
    } else {
        int direction = tri->getTriDirection(direction_tri);
        if (direction == -1) {
            std::ostringstream os;
            os << "Triangle " << direction_tri << " is not a neighbor of triangle " << tidx
               << ".\n";
            ArgErrLog(os.str());
        }
        tri->sdiff(ldidx).setDirectionDcst(direction, dk);
    }

    _updateElement(&tri->sdiff(ldidx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleRDEF::_getTriVDepSReacActive(triangle_global_id tidx,
                                            solver::vdepsreac_global_id vsridx) const {
    AssertLog(tidx < pTris.size());
    AssertLog(vsridx < statedef().countVDepSReacs());

    int host = _getTriHost(tidx);

    TriRDEF* tri = _getTri(tidx);

    solver::vdepsreac_local_id lvsridx = tri->patchdef()->vdepsreacG2L(vsridx);
    if (lvsridx.unknown()) {
        std::ostringstream os;
        os << "Voltage-dependent surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    bool active = false;
    if (host == myRank_World) {
        active = tri->vdepsreac(lvsridx).active();
    }
    return MPI_ConditionalBcast<bool>(
        active, MPI_C_BOOL, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTriVDepSReacActive(triangle_global_id tidx,
                                            solver::vdepsreac_global_id vsridx,
                                            bool act) {
    AssertLog(tidx < pTris.size());
    AssertLog(vsridx < statedef().countVDepSReacs());

    int host = _getTriHost(tidx);

    TriRDEF* tri = _getTri(tidx);

    solver::vdepsreac_local_id lvsridx = tri->patchdef()->vdepsreacG2L(vsridx);
    if (lvsridx.unknown()) {
        std::ostringstream os;
        os << "Voltage-dependent surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }
    if (host != myRank_World) {
        return;
    }

    tri->vdepsreac(lvsridx).setActive(act);
    _updateElement(&tri->vdepsreac(lvsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriSReacH(triangle_global_id tidx, solver::sreac_global_id ridx) const {
    AssertLog(tidx < pTris.size());
    AssertLog(ridx < statedef().countSReacs());

    int host = _getTriHost(tidx);

    TriRDEF* tri = _getTri(tidx);

    solver::sreac_local_id lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    double h = 0;
    if (host == myRank_World) {
        h = tri->sreac(lsridx).h();
    }
    return MPI_ConditionalBcast<double>(h, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriSReacC(triangle_global_id tidx, solver::sreac_global_id ridx) const {
    AssertLog(tidx < pTris.size());
    AssertLog(ridx < statedef().countSReacs());
    int host = _getTriHost(tidx);

    TriRDEF* tri = _getTri(tidx);

    solver::sreac_local_id lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    double c = 0;
    if (host == myRank_World) {
        c = tri->sreac(lsridx).c();
    }
    return MPI_ConditionalBcast<double>(c, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriSReacA(triangle_global_id tidx, solver::sreac_global_id ridx) const {
    AssertLog(tidx < pTris.size());
    AssertLog(ridx < statedef().countSReacs());

    int host = _getTriHost(tidx);

    TriRDEF* tri = _getTri(tidx);

    solver::sreac_local_id lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx.unknown()) {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    double a = 0;
    if (host == myRank_World) {
        a = tri->sreac(lsridx).rate();
    }
    return MPI_ConditionalBcast<double>(a, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getTriRaftCount(triangle_global_id tidx, solver::raft_global_id ridx) const {
    return MPI_ConditionalBcast<uint>(
        0.0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTriRaftCount(triangle_global_id /*tidx*/,
                                      solver::raft_global_id /*ridx*/,
                                      uint /*n*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

solver::raft_individual_id TetVesicleRDEF::_addTriRaft(triangle_global_id tidx,
                                                       solver::raft_global_id ridx) {
    return solver::raft_individual_id(MPI_ConditionalBcast<index_t>(
        0, MPI_STEPS_INDEX, vesraftRank_World, myRank_World, syncOutput, outputRank));
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::setEfieldDT(double efdt) {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    if (efdt <= 0.0) {
        std::ostringstream os;
        os << "EField dt must be greater than zero.";
        ArgErrLog(os.str());
    }
    pEFDT = efdt;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTetV(tetrahedron_global_id tidx) const {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    tetrahedron_local_id loctidx = pEFTet_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    // EField object should convert value to base s.i. units
    double v = pEField->getTetV(loctidx);
    return MPI_ConditionalBcast<double>(
        v, MPI_DOUBLE, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTetV(tetrahedron_global_id tidx, double v) {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    tetrahedron_local_id loctidx = pEFTet_GtoL[tidx];
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

bool TetVesicleRDEF::_getTetVClamped(tetrahedron_global_id tidx) const {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    tetrahedron_local_id loctidx = pEFTet_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    bool clamped = pEField->getTetVClamped(loctidx);
    return MPI_ConditionalBcast<bool>(
        clamped, MPI_C_BOOL, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTetVClamped(tetrahedron_global_id tidx, bool cl) {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    tetrahedron_local_id loctidx = pEFTet_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    pEField->setTetVClamped(loctidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getTriV_(triangle_global_id tidx) const {
    triangle_local_id loctidx = pEFTri_GtoL[tidx];
    AssertLog(loctidx.valid());

    return EFTrisV[loctidx];
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriV(triangle_global_id tidx) const {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    triangle_local_id loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    double v = EFTrisV[loctidx];
    return MPI_ConditionalBcast<double>(
        v, MPI_DOUBLE, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTriV(triangle_global_id tidx, double v) {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    triangle_local_id loctidx = pEFTri_GtoL[tidx];
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

bool TetVesicleRDEF::_getTriVClamped(triangle_global_id tidx) const {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    triangle_local_id loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    bool clamped = pEField->getTriVClamped(loctidx);
    return MPI_ConditionalBcast<bool>(
        clamped, MPI_C_BOOL, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTriVClamped(triangle_global_id tidx, bool cl) {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    triangle_local_id loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    pEField->setTriVClamped(loctidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTriOhmicErev(triangle_global_id tidx,
                                      solver::ohmiccurr_global_id ocgidx,
                                      double erev) {
    if (!_efflag()) {
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

    int host = _getTriHost(tidx);
    TriRDEF* tri = _getTri(tidx);

    if (host == myRank_World) {
        solver::ohmiccurr_local_id locidx = tri->patchdef()->ohmiccurrG2L(ocgidx);
        if (locidx.unknown()) {
            ArgErrLog("Ohmic current undefined in triangle.\n");
        }
        tri->setOCerev(locidx, erev);
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriOhmicErev(triangle_global_id tidx,
                                        solver::ohmiccurr_global_id ocgidx) const {
    if (!_efflag()) {
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

    int host = _getTriHost(tidx);
    TriRDEF* tri = _getTri(tidx);

    double erev = 0.0;
    if (host == myRank_World) {
        solver::ohmiccurr_local_id locidx = tri->patchdef()->ohmiccurrG2L(ocgidx);
        if (locidx.unknown()) {
            ArgErrLog("Ohmic current undefined in triangle.\n");
        }
        erev = tri->getOCerev(locidx);
    }

    return MPI_ConditionalBcast<double>(
        erev, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriOhmicI(triangle_global_id tidx) const {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    int host = _getTriHost(tidx);
    TriRDEF* tri = _getTri(tidx);

    triangle_local_id loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    double cur = 0.0;
    if (host == myRank_World) {
        cur = tri->getOhmicI(EFTrisV[loctidx], getEfieldDT());
    }
    return MPI_ConditionalBcast<double>(
        cur, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriOhmicI(triangle_global_id tidx,
                                     solver::ohmiccurr_global_id ocidx) const {
    AssertLog(tidx < pTris.size());
    AssertLog(ocidx < statedef().countOhmicCurrs());

    triangle_local_id loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    int host = _getTriHost(tidx);
    TriRDEF* tri = _getTri(tidx);

    solver::ohmiccurr_local_id locidx = tri->patchdef()->ohmiccurrG2L(ocidx);
    if (locidx.unknown()) {
        std::ostringstream os;
        os << "Ohmic current undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    double cur = 0.0;
    if (host == myRank_World) {
        cur = tri->getOhmicI(locidx, EFTrisV[loctidx], getEfieldDT());
    }
    return MPI_ConditionalBcast<double>(
        cur, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriGHKI(triangle_global_id tidx) const {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }

    int host = _getTriHost(tidx);
    TriRDEF* tri = _getTri(tidx);

    double cur = 0.0;
    if (host == myRank_World) {
        cur = tri->getGHKI(getEfieldDT());
    }
    return MPI_ConditionalBcast<double>(
        cur, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriGHKI(triangle_global_id tidx,
                                   solver::ghkcurr_global_id ghkidx) const {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }

    int host = _getTriHost(tidx);
    TriRDEF* tri = _getTri(tidx);

    solver::ghkcurr_local_id locidx = tri->patchdef()->ghkcurrG2L(ghkidx);
    if (locidx.unknown()) {
        std::ostringstream os;
        os << "GHK current undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    double cur = 0.0;
    if (host == myRank_World) {
        cur = tri->getGHKI(locidx, getEfieldDT());
    }
    return MPI_ConditionalBcast<double>(
        cur, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriI(triangle_global_id tidx) const {
    return _getTriGHKI(tidx) + _getTriOhmicI(tidx);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getVertIClamp(vertex_id_t vidx) const {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    auto locvidx = pEFVert_GtoL[vidx.get()];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }

    double value = pEField->getVertIClamp(locvidx);
    return MPI_ConditionalBcast<double>(
        value, MPI_DOUBLE, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}
////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setVertIClamp(vertex_id_t vidx, double cur) {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    vertex_id_t locvidx = pEFVert_GtoL[vidx.get()];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }

    // EField object should convert to required units
    pEField->setVertIClamp(locvidx, cur);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getTriIClamp(triangle_global_id tidx) const {
    if (!_efflag()) {
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

    double value = pEField->getTriIClamp(loctidx);
    return MPI_ConditionalBcast<double>(
        value, MPI_DOUBLE, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTriIClamp(triangle_global_id tidx, double cur) {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    triangle_local_id loctidx = pEFTri_GtoL[tidx];
    if (loctidx.unknown()) {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    // EField object should convert to required units
    pEField->setTriIClamp(loctidx, cur);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setTriCapac(triangle_global_id tidx, double cap) {
    if (!_efflag()) {
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

double TetVesicleRDEF::_getVertV(vertex_id_t vidx) const {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    vertex_id_t locvidx = pEFVert_GtoL[vidx.get()];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }
    // EField object should return base s.i. units
    double value = pEField->getVertV(locvidx);
    return MPI_ConditionalBcast<double>(
        value, MPI_DOUBLE, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setVertV(vertex_id_t vidx, double v) {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    vertex_id_t locvidx = pEFVert_GtoL[vidx.get()];
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

bool TetVesicleRDEF::_getVertVClamped(vertex_id_t vidx) const {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    vertex_id_t locvidx = pEFVert_GtoL[vidx.get()];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }

    bool clamped = pEField->getVertVClamped(locvidx);
    return MPI_ConditionalBcast<bool>(
        clamped, MPI_C_BOOL, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setVertVClamped(vertex_id_t vidx, bool cl) {
    if (!_efflag()) {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in "
              "simulation.";
        ArgErrLog(os.str());
    }
    vertex_id_t locvidx = pEFVert_GtoL[vidx.get()];
    if (locvidx.unknown()) {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }
    // EField object should convert to millivolts
    pEField->setVertVClamped(locvidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setMembRes(solver::membrane_global_id midx, double ro, double vrev) {
    if (!_efflag()) {
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

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setMembPotential(solver::membrane_global_id midx, double v) {
    if (!_efflag()) {
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

void TetVesicleRDEF::_setMembCapac(solver::membrane_global_id midx, double cm) {
    if (!_efflag()) {
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

void TetVesicleRDEF::_setMembVolRes(solver::membrane_global_id midx, double ro) {
    if (!_efflag()) {
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

void TetVesicleRDEF::createPath(std::string const& /*path*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::addPathPoint(std::string const& /*path*/,
                                  uint /*point_idx*/,
                                  const std::vector<double>& /*position*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::addPathBranch(std::string const& /*path*/,
                                   uint /*sourcepoint_idx*/,
                                   const std::map<uint, double>& /*destpoints_indxs*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_addPathVesicle(std::string const& /*path_name*/,
                                     solver::vesicle_global_id /*vidx*/,
                                     double /*speed*/,
                                     const std::map<solver::spec_global_id, uint>& /*spec_deps*/,
                                     const std::vector<double>& stoch_stepsize /*stoch_stepsize*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_addVesicleDiffusionGroup(
    solver::vesicle_global_id /*vidx*/,
    const std::vector<solver::comp_global_id>& /*comp_indices*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setVesSReacK(solver::vessreac_global_id vsridx, double kf) {
    AssertLog(kf >= 0.0);

    AssertLog(vsridx < statedef().countVesSReacs());

    statedef().vessreacdef(vsridx).setKcst(kf);

    // Rates have changed
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setExocytosisK(solver::exocytosis_global_id exoidx, double kf) {
    AssertLog(kf >= 0.0);

    AssertLog(exoidx < statedef().countExocytosis());

    statedef().exocytosisdef(exoidx).setKcst(kf);

    // Rates have changed
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getVesSReacExtent(solver::vessreac_global_id vsridx) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getExocytosisExtent(solver::exocytosis_global_id exoidx) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::ExocytosisEvent> TetVesicleRDEF::_getExocytosisEvents(
    solver::exocytosis_global_id exoidx) {
    std::vector<solver::ExocytosisEvent> events;
    MPI_ConditionalBcast(events,
                         dataTypeUtil.MPI_ExocytosisEventSync,
                         vesraftRank_World,
                         myRank_World,
                         syncOutput,
                         outputRank);
    return events;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::_getRaftEndocytosisExtent(solver::raftendocytosis_global_id rendoidx) const {
    return MPI_ConditionalBcast<uint>(
        0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::RaftEndocytosisEvent> TetVesicleRDEF::_getRaftEndocytosisEvents(
    solver::raftendocytosis_global_id rendoidx) {
    std::vector<solver::RaftEndocytosisEvent> events;
    MPI_ConditionalBcast(events,
                         dataTypeUtil.MPI_RaftEndocytosisEventSync,
                         vesraftRank_World,
                         myRank_World,
                         syncOutput,
                         outputRank);
    return events;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setRaftEndocytosisK(solver::raftendocytosis_global_id rendoidx, double kcst) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::_getSingleRaftRaftEndocytosisK(
    solver::raft_global_id ridx,
    solver::raft_individual_id raft_unique_index,
    solver::raftendocytosis_global_id rendoidx) const {
    return MPI_ConditionalBcast<double>(
        0, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setSingleRaftRaftEndocytosisK(solver::raft_global_id ridx,
                                                    solver::raft_individual_id raft_unique_index,
                                                    solver::raftendocytosis_global_id rendoidx,
                                                    double k) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setSingleRaftSReacActive(solver::raft_global_id ridx,
                                               solver::raft_individual_id raft_unique_index,
                                               solver::raftsreac_global_id rsreacidx,
                                               bool active) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleRDEF::_getSingleRaftSReacActive(solver::raft_global_id ridx,
                                               solver::raft_individual_id raft_unique_index,
                                               solver::raftsreac_global_id rsreacidx) const {
    return MPI_ConditionalBcast<bool>(
        true, MPI_C_BOOL, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_setVesicleSurfaceLinkSpecSDiffD(solver::vesicle_global_id /*vidx*/,
                                                      solver::linkspec_global_id /*lsidx*/,
                                                      double /*dcst*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_computeUpdPeriod() {
    double local_max_rate = 0.0;

    for (uint pos = 0; pos < pDiffs.size(); pos++) {
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

    for (uint pos = 0; pos < pSDiffs.size(); pos++) {
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
    MPI_Allreduce(&local_max_rate, &global_max_rate, 1, MPI_DOUBLE, MPI_MAX, RDEFComm);

    if (global_max_rate < 0.0) {
        std::ostringstream os;
        os << "Maximum scaled diffusion constant is " << global_max_rate
           << ". This should not happen in this solver.\n";
        ArgErrLog(os.str());
    }
    updPeriod = 1.0 / global_max_rate;

    if (updPeriod < minUpdPeriod) {
        updPeriod = minUpdPeriod;
    }
    recomputeUpdPeriod = false;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_updateLocal(KProcPSet const& upd_entries) {
    for (auto& kp: upd_entries) {
        AssertLog(kp != nullptr);
        _updateElement(kp);
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_updateLocal(std::vector<KProc*> const& upd_entries) {
    for (auto& kp: upd_entries) {
        AssertLog(kp != nullptr);
        _updateElement(kp);
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_updateLocal(std::vector<uint> const& upd_entries) {
    for (auto& upd_idx: upd_entries) {
        KProc* kp = pKProcs[upd_idx];
        if (kp != nullptr) {
            _updateElement(kp);
        }
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_updateLocal() {
    for (uint i = 0; i < nEntries; i++) {
        KProc* kp = pKProcs[i];
        if (kp != nullptr) {
            _updateElement(kp);
        }
    }
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

CRGroup* TetVesicleRDEF::_getGroup(int pow) {
    if (pow >= 0) {
        return pGroups[pow];
    } else {
        return nGroups[-pow];
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_extendPGroups(uint new_size) {
    auto curr_size = pGroups.size();

    while (curr_size < new_size) {
        pGroups.push_back(new CRGroup(curr_size));
        curr_size++;
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_extendNGroups(uint new_size) {
    uint curr_size = nGroups.size();

    while (curr_size < new_size) {
        nGroups.push_back(new CRGroup(-curr_size));
        curr_size++;
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_extendGroup(CRGroup* group, uint size) {
    group->capacity += size;
    group->indices = static_cast<KProc**>(
        realloc(group->indices, sizeof(KProc*) * group->capacity));
    if (group->indices == nullptr) {
        SysErrLog("DirectCR: unable to allocate memory for SSA group.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_updateSum() {
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

void TetVesicleRDEF::_updateElement(KProc* kp)  // OK, function checked
{
    double new_rate = kp->rate(this);
    CRKProcData& data = kp->crData;

    if (kp->getType() == KP_DIFF) {
        if (new_rate <= 0.0 && data.rate > 0.0) {
            diffSep--;
            Diff* temp = pDiffs[diffSep];
            uint swap_pos = data.pos;
            data.pos = diffSep;
            pDiffs[diffSep] = dynamic_cast<Diff*>(kp);
            pDiffs[swap_pos] = temp;
            temp->crData.pos = swap_pos;
        } else if (new_rate > 0.0 && data.rate <= 0.0) {
            Diff* temp = pDiffs[diffSep];
            uint swap_pos = data.pos;
            data.pos = diffSep;
            pDiffs[diffSep] = dynamic_cast<Diff*>(kp);
            pDiffs[swap_pos] = temp;
            temp->crData.pos = swap_pos;
            diffSep++;
        }
        data.rate = new_rate;
        return;
    } else if (kp->getType() == KP_SDIFF) {
        if (new_rate <= 0.0 && data.rate > 0.0) {
            sdiffSep--;
            SDiff* temp = pSDiffs[sdiffSep];
            uint swap_pos = data.pos;
            data.pos = sdiffSep;
            pSDiffs[sdiffSep] = dynamic_cast<SDiff*>(kp);
            pSDiffs[swap_pos] = temp;
            temp->crData.pos = swap_pos;
        } else if (new_rate > 0.0 && data.rate <= 0.0) {
            SDiff* temp = pSDiffs[sdiffSep];
            uint swap_pos = data.pos;
            data.pos = sdiffSep;
            pSDiffs[sdiffSep] = dynamic_cast<SDiff*>(kp);
            pSDiffs[swap_pos] = temp;
            temp->crData.pos = swap_pos;
            sdiffSep++;
        }
        data.rate = new_rate;
        return;
    }

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
        frexp(new_rate, &new_pow);

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

            if (nGroups.size() <= static_cast<unsigned>(-new_pow)) {
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
// ROI Data Access
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleRDEF::getBatchTetSpecCounts(std::vector<index_t> const& tets,
                                                          std::string const& s) const {
    size_t ntets = tets.size();
    std::vector<double> counts(ntets, 0.0);
    getBatchTetSpecCountsNP(tets.data(), ntets, s, counts.data(), ntets);
    return counts;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleRDEF::getBatchTriSpecCounts(std::vector<index_t> const& tris,
                                                          std::string const& s) const {
    size_t ntris = tris.size();
    std::vector<double> counts(ntris, 0.0);
    getBatchTriSpecCountsNP(tris.data(), ntris, s, counts.data(), ntris);
    return counts;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::getBatchTetSpecCountsNP(const index_t* indices,
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

        TetRDEF* tet = pTets[tidx];
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
    MPI_ConditionalReduce<double>(
        local_counts.data(), counts, input_size, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::getBatchTriSpecCountsNP(const index_t* indices,
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
    std::vector<double> local_counts(input_size, 0.0);

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

        TriRDEF* tri = pTris[tidx];
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
    MPI_ConditionalReduce<double>(
        local_counts.data(), counts, input_size, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// ROI Data Access
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleRDEF::getROITetSpecCounts(const std::string& ROI_id,
                                                        std::string const& s) const {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    const auto size = roi->second.size();
    std::vector<double> data(size);
    getBatchTetSpecCountsNP(
        reinterpret_cast<const index_t*>(roi->second.data()), size, s, &data.front(), data.size());
    return data;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleRDEF::getROITriSpecCounts(const std::string& ROI_id,
                                                        std::string const& s) const {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TRI>()) {
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

void TetVesicleRDEF::getROITetSpecCountsNP(const std::string& ROI_id,
                                           std::string const& s,
                                           double* counts,
                                           size_t output_size) const {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    getBatchTetSpecCountsNP(reinterpret_cast<const index_t*>(roi->second.data()),
                            roi->second.size(),
                            s,
                            counts,
                            output_size);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::getROITriSpecCountsNP(const std::string& ROI_id,
                                           std::string const& s,
                                           double* counts,
                                           size_t output_size) const {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    getBatchTriSpecCountsNP(reinterpret_cast<const index_t*>(roi->second.data()),
                            roi->second.size(),
                            s,
                            counts,
                            output_size);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getROIVol(const std::string& ROI_id) const {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    double sum = 0.0;
    for (auto const& tidx: roi->second) {
        sum += pTets[tidx]->staticVol();
    }
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getROIArea(const std::string& ROI_id) const {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    double sum = 0.0;
    for (auto const& tidx: roi->second) {
        sum += pTris[tidx]->area();
    }
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getROITetSpecCount(const std::vector<tetrahedron_global_id>& tetrahedrons,
                                          const std::string& s) const {
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;
    double local_sum = 0.0;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: tetrahedrons) {
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

        TetRDEF* tet = pTets[tidx];
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
    return MPI_ConditionalReduce<double>(local_sum, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getROITriSpecCount(const std::vector<triangle_global_id>& triangles,
                                          const std::string& s) const {
    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;
    double local_sum = 0.0;

    const solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: triangles) {
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

        TriRDEF* tri = pTris[tidx];
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
    return MPI_ConditionalReduce<double>(local_sum, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getROISpecCount(const std::string& ROI_id, std::string const& s) const {
    {
        auto const& roi =
            _mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id, 0 /* count */, false /* warning */);
        if (roi != _mesh()->rois.end<tetmesh::ROI_TRI>()) {
            return getROITriSpecCount(roi->second, s);
        }
    }
    {
        auto const& roi =
            _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id, 0 /* count */, false /* warning */);
        if (roi != _mesh()->rois.end<tetmesh::ROI_TET>()) {
            return getROITetSpecCount(roi->second, s);
        }
    }
    ArgErrLog("Error: Cannot find suitable ROI for the function call getROICount.\n");
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::setROITriSpecCount(const std::vector<triangle_global_id>& triangles,
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

        TriRDEF* tri = pTris[tidx];
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
    if (myRank_RDEF == RDEFmasterRank_RDEF) {
        uint c = static_cast<uint>(count);
        uint nremoved = 0;

        for (uint t = 0; t < ind_size; t++) {
            auto tidx = apply_indices[t];
            TriRDEF* tri = pTris[tidx];

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
                TriRDEF* tri = pTris[tidx];
                accum += tri->area();
                if (selector < accum) {
                    apply_count[t] += 1.0;
                    break;
                }
            }
            c--;
        }
    }
    if (nHosts_RDEF != 1) {
        MPI_Bcast(apply_count.data(), ind_size, MPI_DOUBLE, RDEFmasterRank_RDEF, RDEFComm);
    }

    // counts need to be set globally so that sync can be avoided
    for (uint t = 0; t < ind_size; t++) {
        auto tidx = apply_indices[t];
        TriRDEF* tri = pTris[tidx];
        solver::spec_local_id slidx = tri->patchdef()->specG2L(sgidx);
        tri->setCount(slidx, apply_count[t]);
        _updateSpec(tri, sgidx);
    }
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::setROITetSpecCount(const std::vector<tetrahedron_global_id>& tetrahedrons,
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

        TetRDEF* tet = pTets[tidx];
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
    if (myRank_RDEF == RDEFmasterRank_RDEF) {
        uint c = static_cast<uint>(count);
        uint nremoved = 0;

        for (uint t = 0; t < ind_size; t++) {
            auto tidx = apply_indices[t];
            TetRDEF* tet = pTets[tidx];

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
                TetRDEF* tet = pTets[tidx];
                accum += tet->vol();
                if (selector < accum) {
                    apply_count[t] += 1.0;
                    break;
                }
            }
            c--;
        }
    }

    if (nHosts_RDEF != 1) {
        MPI_Bcast(apply_count.data(), ind_size, MPI_DOUBLE, RDEFmasterRank_RDEF, RDEFComm);
    }
    // set the counts golbally and update local KProcs
    for (uint t = 0; t < ind_size; t++) {
        auto tidx = apply_indices[t];
        TetRDEF* tet = pTets[tidx];
        solver::spec_local_id slidx = tet->compdef()->specG2L(sgidx);
        tet->setCount(slidx, apply_count[t]);
        _updateSpec(tet, sgidx);
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

// WEILIANG: Can we apply Sam's set count method (as in e.g. setCompCount) here?
void TetVesicleRDEF::setROISpecCount(const std::string& ROI_id,
                                     std::string const& s,
                                     double count) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (count > UINT_MAX) {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << UINT_MAX << ").\n";
        ArgErrLog(os.str());
    }

    auto const& tri_roi = _mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id, 0, false);
    if (tri_roi != _mesh()->rois.end<tetmesh::ROI_TRI>()) {
        setROITriSpecCount(tri_roi->second, s, count);
    } else {
        auto const& tet_roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id, 0, false);
        if (tet_roi != _mesh()->rois.end<tetmesh::ROI_TET>()) {
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

double TetVesicleRDEF::getROISpecAmount(const std::string& ROI_id, std::string const& s) const {
    double count = getROISpecCount(ROI_id, s);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::setROISpecAmount(const std::string& ROI_id,
                                      std::string const& s,
                                      double amount) {
    setROISpecCount(ROI_id, s, amount * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////
// WEILIANG: Can we apply Sam's set count method (as in e.g. setCompCount) here?
void TetVesicleRDEF::setROISpecConc(const std::string& ROI_id, std::string const& s, double conc) {
    double vol = getROIVol(ROI_id);
    int count = conc * (1.0e3 * vol * math::AVOGADRO);
    setROISpecCount(ROI_id, s, count);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getROISpecConc(const std::string& ROI_id, const std::string& s) const {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    const double count = getROITetSpecCount(roi->second, s);
    double vol = getROIVol(ROI_id);
    return count / (1.0e3 * vol * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::setROITriSpecClamped(const std::vector<triangle_global_id>& triangles,
                                          const std::string& s,
                                          bool b) {
    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: triangles) {
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

        TriRDEF* tri = pTris[tidx];
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

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::setROITetSpecClamped(const std::vector<tetrahedron_global_id>& tetrahedrons,
                                          std::string const& s,
                                          bool b) {
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;

    solver::spec_global_id sgidx = statedef().getSpecIdx(s);

    for (auto const& tidx: tetrahedrons) {
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

        TetRDEF* tet = pTets[tidx];
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

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::setROISpecClamped(const std::string& ROI_id, std::string const& s, bool b) {
    {
        auto const& roi = _mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id, 0, false);
        if (roi != _mesh()->rois.end<tetmesh::ROI_TRI>()) {
            setROITriSpecClamped(roi->second, s, b);
            return;
        }
    }
    {
        auto const& roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id, 0, false);
        if (roi != _mesh()->rois.end<tetmesh::ROI_TET>()) {
            setROITetSpecClamped(roi->second, s, b);
            return;
        }
    }
    std::ostringstream os;
    os << "Error: Cannot find suitable ROI for the function call getROICount.\n";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::setROIReacK(const std::string& ROI_id, std::string const& r, double kf) {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_reac_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream reac_undefined;

    solver::reac_global_id rgidx = statedef().getReacIdx(r);

    for (auto const& tidx: roi->second) {
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

        TetRDEF* tet = pTets[tidx];
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

void TetVesicleRDEF::setROISReacK(const std::string& ROI_id, std::string const& sr, double kf) {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tri_warning = false;
    bool has_sreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream sreac_undefined;

    solver::sreac_global_id srgidx = statedef().getSReacIdx(sr);

    for (auto const& tidx: roi->second) {
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

        TriRDEF* tri = pTris[tidx];
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

void TetVesicleRDEF::setROIDiffD(const std::string& ROI_id, std::string const& d, double dk) {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_diff_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream diff_undefined;

    solver::diff_global_id dgidx = statedef().getDiffIdx(d);

    for (auto const& tidx: roi->second) {
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

        TetRDEF* tet = pTets[tidx];
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

void TetVesicleRDEF::setROIReacActive(const std::string& ROI_id, std::string const& r, bool a) {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_reac_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream reac_undefined;

    solver::reac_global_id rgidx = statedef().getReacIdx(r);

    for (auto const& tidx: roi->second) {
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

        TetRDEF* tet = pTets[tidx];
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

void TetVesicleRDEF::setROISReacActive(const std::string& ROI_id, std::string const& sr, bool a) {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tri_warning = false;
    bool has_sreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream sreac_undefined;

    solver::sreac_global_id srgidx = statedef().getSReacIdx(sr);

    for (auto const& tidx: roi->second) {
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

        TriRDEF* tri = pTris[tidx];
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

void TetVesicleRDEF::setROIDiffActive(const std::string& ROI_id, std::string const& d, bool a) {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_diff_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream diff_undefined;

    solver::diff_global_id dgidx = statedef().getDiffIdx(d);

    for (auto const& tidx: roi->second) {
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

        TetRDEF* tet = pTets[tidx];
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

void TetVesicleRDEF::setROIVDepSReacActive(const std::string& ROI_id,
                                           std::string const& vsr,
                                           bool a) {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tri_warning = false;
    bool has_vsreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream vsreac_undefined;

    solver::vdepsreac_global_id vsrgidx = statedef().getVDepSReacIdx(vsr);

    for (auto const& tidx: roi->second) {
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

        TriRDEF* tri = pTris[tidx];
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

unsigned long long TetVesicleRDEF::getROIReacExtent(const std::string& ROI_id,
                                                    std::string const& r) const {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_reac_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream reac_undefined;

    solver::reac_global_id rgidx = statedef().getReacIdx(r);

    unsigned long long local_sum = 0;

    for (auto const& tidx: roi->second) {
        if (tidx >= pTets.size()) {
            std::ostringstream oss;
            oss << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            ArgErrLog(oss.str());
        }

        if (pTets[tidx] == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        TetRDEF* tet = pTets[tidx];
        solver::reac_local_id rlidx = tet->compdef()->reacG2L(rgidx);
        if (rlidx.unknown()) {
            reac_undefined << tidx << " ";
            has_reac_warning = true;
            continue;
        }

        if (tet->getInHost()) {
            local_sum += tet->reac(rlidx).getExtent();
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

    return MPI_ConditionalReduce<unsigned long long>(
        local_sum, MPI_UNSIGNED_LONG_LONG, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::resetROIReacExtent(const std::string& ROI_id, std::string const& r) {
    std::ostringstream os;
    os << "This function has not been implemented!";
    NotImplErrLog(os.str());

    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_reac_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream reac_undefined;

    solver::reac_global_id rgidx = statedef().getReacIdx(r);

    for (auto const& tidx: roi->second) {
        if (tidx >= pTets.size()) {
            std::ostringstream oss;
            oss << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            ArgErrLog(oss.str());
        }

        if (pTets[tidx] == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        TetRDEF* tet = pTets[tidx];
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

unsigned long long TetVesicleRDEF::getROISReacExtent(const std::string& ROI_id,
                                                     std::string const& sr) const {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tri_warning = false;
    bool has_sreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream sreac_undefined;

    solver::sreac_global_id srgidx = statedef().getSReacIdx(sr);

    unsigned long long local_sum = 0;

    for (auto const& tidx: roi->second) {
        if (tidx >= pTris.size()) {
            std::ostringstream oss;
            oss << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            ArgErrLog(oss.str());
        }

        if (pTris[tidx] == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        TriRDEF* tri = pTris[tidx];
        solver::sreac_local_id srlidx = tri->patchdef()->sreacG2L(srgidx);
        if (srlidx.unknown()) {
            sreac_undefined << tidx << " ";
            has_sreac_warning = true;
            continue;
        }

        if (tri->getInHost()) {
            local_sum += tri->sreac(srlidx).getExtent();
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

    return MPI_ConditionalReduce<unsigned long long>(
        local_sum, MPI_UNSIGNED_LONG_LONG, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::resetROISReacExtent(const std::string& ROI_id, std::string const& sr) {
    std::ostringstream os;
    os << "This function has not been implemented!";
    NotImplErrLog(os.str());

    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TRI>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TRI>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tri_warning = false;
    bool has_sreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream sreac_undefined;

    solver::sreac_global_id srgidx = statedef().getSReacIdx(sr);

    for (auto const& tidx: roi->second) {
        if (tidx >= pTris.size()) {
            std::ostringstream oss;
            oss << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            ArgErrLog(oss.str());
        }

        if (pTris[tidx] == nullptr) {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        TriRDEF* tri = pTris[tidx];
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

unsigned long long TetVesicleRDEF::getROIDiffExtent(const std::string& ROI_id,
                                                    std::string const& d) const {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_diff_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream diff_undefined;

    solver::diff_global_id dgidx = statedef().getDiffIdx(d);

    unsigned long long local_sum = 0;

    for (auto const& tidx: roi->second) {
        if (tidx >= pTets.size()) {
            std::ostringstream oss;
            oss << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            ArgErrLog(oss.str());
        }

        if (pTets[tidx] == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        TetRDEF* tet = pTets[tidx];
        solver::diff_local_id dlidx = tet->compdef()->diffG2L(dgidx);
        if (dlidx.unknown()) {
            diff_undefined << tidx << " ";
            has_diff_warning = true;
            continue;
        }

        if (tet->getInHost()) {
            local_sum += tet->diff(dlidx).getExtent();
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

    return MPI_ConditionalReduce<unsigned long long>(
        local_sum, MPI_UNSIGNED_LONG_LONG, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::resetROIDiffExtent(const std::string& ROI_id, std::string const& d) {
    std::ostringstream os;
    os << "This function has not been implemented!";
    NotImplErrLog(os.str());

    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }

    bool has_tet_warning = false;
    bool has_diff_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream diff_undefined;

    solver::diff_global_id dgidx = statedef().getDiffIdx(d);

    for (auto const& tidx: roi->second) {
        if (tidx >= pTets.size()) {
            std::ostringstream oss;
            oss << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            ArgErrLog(oss.str());
        }

        if (pTets[tidx] == nullptr) {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        TetRDEF* tet = pTets[tidx];
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

void TetVesicleRDEF::_updateLocal(uint* upd_entries, uint buffer_size) {
    for (uint i = 0; i < buffer_size; i++) {
        if (pKProcs[upd_entries[i]] != nullptr) {
            _updateElement(pKProcs[upd_entries[i]]);
        }
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::addDiff_(Diff* diff) {
    diff->crData.pos = pDiffs.size();
    pDiffs.push_back(diff);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::addSDiff_(SDiff* sdiff) {
    sdiff->crData.pos = pSDiffs.size();
    pSDiffs.push_back(sdiff);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::addNeighHost_(int host) {
    neighbHosts.insert(host);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::registerBoundaryTet_(TetRDEF* tet) {
    boundaryTets.insert(tet);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::registerBoundaryTri_(TriRDEF* tri) {
    boundaryTris.insert(tri);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleRDEF::registerRemoteMoleculeChange_(int svol_host,
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

void TetVesicleRDEF::_remoteSyncAndUpdate(void* requests,
                                          std::vector<KProc*>& applied_diffs,
                                          std::vector<int>& directions) {
#ifdef MPI_PROFILING
    double starttime = MPI_Wtime();
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
    double endtime = MPI_Wtime();
    syncTime += (endtime - starttime);
#endif

    while (!await_neighbors.empty()) {
#ifdef MPI_PROFILING
        starttime = MPI_Wtime();
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
        endtime = MPI_Wtime();
        idleTime += (endtime - starttime);
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

            if (type == SUB_TET) {
                auto tet = pTets[tetrahedron_global_id(idx)];
                tet->incCount(slidx, value);
                std::vector<KProc*> const& remote_upd = tet->getSpecUpdKProcs(slidx);
                upd_kprocs.insert(remote_upd.begin(), remote_upd.end());
            }
            if (type == SUB_TRI) {
                auto tri = pTris[triangle_global_id(idx)];
                tri->incCount(slidx, value);
                std::vector<KProc*> const& remote_upd = tri->getSpecUpdKProcs(slidx);
                upd_kprocs.insert(remote_upd.begin(), remote_upd.end());
            }
        }

        await_neighbors.erase(data_source);

#ifdef MPI_PROFILING
        endtime = MPI_Wtime();
        syncTime += (endtime - starttime);
#endif
    }

#ifdef MPI_PROFILING
    starttime = MPI_Wtime();
#endif

    auto napply = applied_diffs.size();
    for (uint i = 0; i < napply; i++) {
        KProc* kp = applied_diffs[i];
        _updateElement(kp);  // First update the applied kproc itself

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
    endtime = MPI_Wtime();
    compTime += (endtime - starttime);
#endif
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::setDiffApplyThreshold(int threshold) {
    diffApplyThreshold = threshold;
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetVesicleRDEF::getReacExtent(bool local) {
    if (local) {
        return reacExtent;
    }

    unsigned long long sum;
    MPI_Reduce(&reacExtent, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, RDEFComm);
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetVesicleRDEF::getDiffExtent(bool local) {
    if (local) {
        return diffExtent;
    }
    unsigned long long sum;
    MPI_Reduce(&diffExtent, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, RDEFComm);
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getNIteration() const {
    return nIteration;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getCompTime() const {
    return compTime;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getSyncTime() const {
    return syncTime;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getIdleTime() const {
    return idleTime;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getEFieldTime() const {
    return efieldTime;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getRDTime() const {
    return rdTime;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleRDEF::getDataExchangeTime() const {
    return dataExchangeTime;
}

////////////////////////////////////////////////////////////////////////////////

TetRDEF* TetVesicleRDEF::_getTet(tetrahedron_global_id tgidx) const {
    if (pTets[tgidx] == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tgidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }
    return pTets[tgidx];
}

////////////////////////////////////////////////////////////////////////////////

TriRDEF* TetVesicleRDEF::_getTri(triangle_global_id tgidx) const {
    if (pTris[tgidx] == nullptr) {
        std::ostringstream os;
        os << "Triangle " << tgidx << " has not been assigned to a patch.\n";
        ArgErrLog(os.str());
    }
    return pTris[tgidx];
}

////////////////////////////////////////////////////////////////////////////////

int TetVesicleRDEF::_getTetHost(tetrahedron_global_id tgidx) const {
    auto host_result = tetHosts.find(tgidx);
    if (host_result == tetHosts.end()) {
        std::ostringstream os;
        os << "Tetrahedron " << tgidx << " has not been assigned to a host.\n";
        ArgErrLog(os.str());
    }
    return host_result->second;
}

////////////////////////////////////////////////////////////////////////////////

int TetVesicleRDEF::_getTriHost(triangle_global_id tgidx) const {
    auto host_result = triHosts.find(tgidx);
    if (host_result == triHosts.end()) {
        std::ostringstream os;
        os << "Triangle " << tgidx << " has not been assigned to a host.\n";
        ArgErrLog(os.str());
    }
    return host_result->second;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_regTetPoolSync(tetrahedron_global_id tet_gidx,
                                     solver::spec_global_id spec_gidx,
                                     uint count) {
    tetPoolCountSyncs_Vec.push_back({tet_gidx.get(), spec_gidx, count});
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_regTriPoolSync(triangle_global_id tri_gidx,
                                     solver::spec_global_id spec_gidx,
                                     uint count) {
    triPoolCountSyncs_Vec.push_back({tri_gidx.get(), spec_gidx, count});
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_syncPools(SyncDirection direction) {
    switch (direction) {
    case VESRAFT_TO_RDEF: {
        MPI_BcastVec<PoolCountSync>(tetPoolCountSyncs_Vec,
                                    dataTypeUtil.MPI_PoolCountSync,
                                    vesraftRank_World,
                                    myRank_World,
                                    MPI_COMM_WORLD);
        MPI_BcastVec<PoolCountSync>(triPoolCountSyncs_Vec,
                                    dataTypeUtil.MPI_PoolCountSync,
                                    vesraftRank_World,
                                    myRank_World,
                                    MPI_COMM_WORLD);

        // The tet pool data received from VesRaft is an increase (not the absolute number)
        for (auto const& entry: tetPoolCountSyncs_Vec) {
            TetRDEF* tet = pTets[tetrahedron_global_id(entry.container_global_index)];
            AssertLog(tet != nullptr);  // TODO remove after testing
            solver::spec_local_id lsidx(
                tet->compdef()->specG2L(solver::spec_global_id(entry.spec_global_index)));
            AssertLog(entry.count > 0);
            if (!tet->clamped(lsidx)) {
                tet->incCount(lsidx, entry.count, 0.0, true);
            }
        }

        // This data is the count, not an inc, so setCount is appropriate here
        for (auto const& entry: triPoolCountSyncs_Vec) {
            TriRDEF* tri = pTris[triangle_global_id(entry.container_global_index)];
            AssertLog(tri != nullptr);  // TODO remove after testing
            solver::spec_local_id lsidx(
                tri->patchdef()->specG2L(solver::spec_global_id(entry.spec_global_index)));
            if (!tri->clamped(lsidx)) {
                tri->setCount(lsidx, entry.count);
            }
        }

        tetPoolCountSyncs_Vec.clear();
        triPoolCountSyncs_Vec.clear();

        break;
    }
    case RDEF_TO_VESRAFT: {
        MPI_GatherVec<PoolCountSync>(tetPoolCountSyncs_Vec,
                                     dataTypeUtil.MPI_PoolCountSync,
                                     vesraftRank_World,
                                     myRank_World,
                                     nHosts_World,
                                     MPI_COMM_WORLD);
        MPI_GatherVec<PoolCountSync>(triPoolCountSyncs_Vec,
                                     dataTypeUtil.MPI_PoolCountSync,
                                     vesraftRank_World,
                                     myRank_World,
                                     nHosts_World,
                                     MPI_COMM_WORLD);

        tetPoolCountSyncs_Vec.clear();
        triPoolCountSyncs_Vec.clear();

        break;
    }
    default:
        break;
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_useVesV2R() {
    // TODO RECEIVE FLAG FROM MASTER AND ONLY DO ALL THIS IF NECESSARY
    // which flag?

    // receive vesProxyV2R_Vec, vesSurfSpecV2R_Vec, vesLinkSpecV2R_Vec data
    // from VesRaft
    MPI_BcastVec<TetV2R>(
        tetV2R_Vec, dataTypeUtil.MPI_TetV2R, vesraftRank_World, myRank_World, MPI_COMM_WORLD);
    MPI_BcastVec<VesProxyV2R>(vesProxyV2R_Vec,
                              dataTypeUtil.MPI_VesProxyV2R,
                              vesraftRank_World,
                              myRank_World,
                              MPI_COMM_WORLD);
    MPI_BcastVec<VesSurfSpecV2R>(vesSurfSpecV2R_Vec,
                                 dataTypeUtil.MPI_VesSurfSpecV2R,
                                 vesraftRank_World,
                                 myRank_World,
                                 MPI_COMM_WORLD);
    MPI_BcastVec<VesLinkSpecV2R>(vesLinkSpecV2R_Vec,
                                 dataTypeUtil.MPI_VesLinkSpecV2R,
                                 vesraftRank_World,
                                 myRank_World,
                                 MPI_COMM_WORLD);

    for (auto const& tet: pTets) {
        if (tet == nullptr or !tet->getInHost()) {
            continue;
        }
        tet->clearVesProxyrefs();
    }

    for (auto const& tet_v2r: tetV2R_Vec) {
        if (tetHosts.find(tet_v2r.tetrahedron_global_index) == tetHosts.end()) {
            ProgErrLog("Tetrahedron not assigned host.\n");
        }

        TetRDEF* tet = _tet(tet_v2r.tetrahedron_global_index);

        if (tet == nullptr) {
            ProgErrLog("Tetrahedron not assigned to a compartment.\n");
        }

        if (!tet->getInHost()) {
            continue;
        }

        tet->setOverlap(tet_v2r.overlap);
    }

    for (auto const& ves_proxy_v2r: vesProxyV2R_Vec) {
        if (tetHosts.find(ves_proxy_v2r.tetrahedron_global_index) == tetHosts.end()) {
            ProgErrLog("Tetrahedron not assigned host.\n");
        }

        TetRDEF* tet = _tet(ves_proxy_v2r.tetrahedron_global_index);

        if (tet == nullptr) {
            ProgErrLog("Tetrahedron not assigned to a compartment.\n");
        }

        if (!tet->getInHost()) {
            continue;
        }

        const solver::vesicle_global_id ves_gidx(ves_proxy_v2r.vesicle_global_index);
        auto& ves_def = statedef().vesicledef(ves_gidx);

        const solver::vesicle_individual_id ves_individual_id(
            ves_proxy_v2r.vesicle_individual_index);

        bool contains_link = ves_proxy_v2r.contains_link;

        math::position_abs ves_pos{ves_proxy_v2r.vesicle_central_position[0],
                                   ves_proxy_v2r.vesicle_central_position[1],
                                   ves_proxy_v2r.vesicle_central_position[2]};

        tet->createVesProxyref(&ves_def, ves_individual_id, ves_pos, contains_link);
    }

    for (auto const& ves_surfspec_v2r: vesSurfSpecV2R_Vec) {
        const tetrahedron_global_id tet_gidx(ves_surfspec_v2r.tetrahedron_global_index);

        if (tetHosts.find(tet_gidx) == tetHosts.end()) {
            ProgErrLog("Tetrahedron not assigned host.\n");
        }

        TetRDEF* tet = _tet(tet_gidx);

        if (tet == nullptr) {
            ProgErrLog("Tetrahedron not assigned to a compartment.\n");
        }

        if (!tet->getInHost()) {
            continue;
        }

        const solver::vesicle_individual_id ves_individual_id(
            ves_surfspec_v2r.vesicle_individual_index);

        // Following will error check the ID, which should exist. So no need to
        // check for a nullptr or anything
        VesProxy* vesproxy = tet->getVesProxyref(ves_individual_id);

        const solver::spec_global_id spec_gidx(ves_surfspec_v2r.surface_spec_global_index);
        const solver::pointspec_individual_id spec_idx(
            ves_surfspec_v2r.surface_spec_individual_index);

        math::position_abs spec_pos_abs{ves_surfspec_v2r.surface_spec_position_abs[0],
                                        ves_surfspec_v2r.surface_spec_position_abs[1],
                                        ves_surfspec_v2r.surface_spec_position_abs[2]};

        vesproxy->addSurfSpec(spec_gidx, spec_idx, spec_pos_abs);
    }

    for (auto const& ves_linkspec_v2r: vesLinkSpecV2R_Vec) {
        const tetrahedron_global_id tet_gidx(ves_linkspec_v2r.tetrahedron_global_index);

        if (tetHosts.find(tet_gidx) == tetHosts.end()) {
            ProgErrLog("Tetrahedron not assigned host.\n");
        }

        TetRDEF* tet = _tet(tet_gidx);

        if (tet == nullptr) {
            ProgErrLog("Tetrahedron not assigned to a compartment.\n");
        }

        if (!tet->getInHost()) {
            continue;
        }

        const solver::vesicle_individual_id ves_individual_id(
            ves_linkspec_v2r.vesicle_individual_index);

        const solver::linkspec_global_id linkspec_global_id(ves_linkspec_v2r.linkspec_global_index);

        const solver::linkspec_individual_id linkspec_individual_id(
            ves_linkspec_v2r.linkspec_individual_index);

        math::position_abs linkspec_pos_abs{ves_linkspec_v2r.linkspec_position_abs[0],
                                            ves_linkspec_v2r.linkspec_position_abs[1],
                                            ves_linkspec_v2r.linkspec_position_abs[2]};

        VesProxy* vesproxy = tet->getVesProxyref(ves_individual_id);

        vesproxy->addLinkSpec(linkspec_individual_id, linkspec_global_id, linkspec_pos_abs);
        // TODO makesure this is enough. Is position needed??
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_constructVesR2V() {
    vesProxyR2V_Vec.clear();
    vesSurfSpecR2V_Vec.clear();
    vesInnerSpecR2V_Vec.clear();
    vesLinkSpecPairR2V_Vec.clear();
    vesLinkSpecR2V_Vec.clear();

    for (auto const& tet: pTets) {
        if (tet == nullptr or !tet->getInHost()) {
            continue;
        }

        tetrahedron_global_id tet_gidx = tet->idx();

        for (auto const& linkspec_newpair: tet->getNewLinkedSpecs()) {
            VesLinkSpecPairR2V ves_linkspec_pair{};
            ves_linkspec_pair.tetrahedron_global_index = tet_gidx;
            ves_linkspec_pair.vesicle1_individual_index = linkspec_newpair.vesicle1_individual_id;
            ves_linkspec_pair.vesicle2_individual_index = linkspec_newpair.vesicle2_individual_id;
            ves_linkspec_pair.linkspec1_global_index = linkspec_newpair.linkspec1_global_id;
            ves_linkspec_pair.linkspec2_global_index = linkspec_newpair.linkspec2_global_id;
            ves_linkspec_pair.linkspec1_individual_index = linkspec_newpair.linkspec1_individual_id;
            ves_linkspec_pair.linkspec2_individual_index = linkspec_newpair.linkspec2_individual_id;
            ves_linkspec_pair.linkspec1_position_abs = linkspec_newpair.linkspec1_pos_absolute;
            ves_linkspec_pair.linkspec2_position_abs = linkspec_newpair.linkspec2_pos_absolute;
            ves_linkspec_pair.min_length = linkspec_newpair.min_length;
            ves_linkspec_pair.max_length = linkspec_newpair.max_length;

            vesLinkSpecPairR2V_Vec.emplace_back(ves_linkspec_pair);
        }

        tet->clearNewLinkedSpecs();


        for (auto const& ves: tet->getVesProxyrefs()) {
            VesProxyR2V ves_proxy_r2v{};

            // solver::vesicle_global_id ves_gidx = ves.second->idx();
            solver::vesicle_individual_id ves_uniqueidx = ves.second->getUniqueIndex();

            ves_proxy_r2v.tetrahedron_global_index = tet_gidx;
            // ves_proxy_r2v.vesicle_global_index = ves_gidx.get();
            ves_proxy_r2v.vesicle_individual_index = ves_uniqueidx;
            ves_proxy_r2v.immobility_update = ves.second->getImmobilityUpdate();
            ves_proxy_r2v.exo_applied_global_index = ves.second->exoApplied();

            vesProxyR2V_Vec.emplace_back(ves_proxy_r2v);

            for (auto const& surf_specs: ves.second->getSurfSpecs()) {
                // surf_specs map .second holds vector of all species of this global idx
                // (global idx is surf_specs.first)
                for (auto const& surf_spec: surf_specs.second) {
                    VesSurfSpecR2V ves_surfspec_r2v{};

                    ves_surfspec_r2v.tetrahedron_global_index = tet_gidx;
                    ves_surfspec_r2v.vesicle_individual_index = ves_uniqueidx;
                    ves_surfspec_r2v.surface_spec_global_index = surf_specs.first;
                    ves_surfspec_r2v.surface_spec_individual_index = surf_spec.first;
                    // OK- vesproxy holds surface species absolute positions
                    ves_surfspec_r2v.surface_spec_position_abs = surf_spec.second;

                    vesSurfSpecR2V_Vec.emplace_back(ves_surfspec_r2v);
                }
            }

            for (auto const& inner_specs: ves.second->getSpecCounts_I()) {
                VesInnerSpecR2V ves_innerspec_r2v{};

                ves_innerspec_r2v.tetrahedron_global_index = tet_gidx;
                ves_innerspec_r2v.vesicle_individual_index = ves_uniqueidx;
                ves_innerspec_r2v.inner_spec_global_index = inner_specs.first;
                ves_innerspec_r2v.count = inner_specs.second;

                vesInnerSpecR2V_Vec.emplace_back(ves_innerspec_r2v);
            }

            for (auto const& link_spec_uniqueid: ves.second->getLinkSpecUpd()) {
                VesLinkSpecR2V ves_linkspec_r2v{};

                ves_linkspec_r2v.tetrahedron_global_index = tet_gidx;
                ves_linkspec_r2v.linkspec_individual_index = link_spec_uniqueid;
                ves_linkspec_r2v.linkspec_global_index = ves.second->getLinkSpecGidx(
                    link_spec_uniqueid);
                ves_linkspec_r2v.vesicle_individual_index = ves_uniqueidx;

                vesLinkSpecR2V_Vec.emplace_back(ves_linkspec_r2v);
            }
            ves.second->clearLinkSpecUpd();
        }
    }

    // send vesProxyR2V_Vec, vesSurfSpecR2V_Vec, vesInnerSpecR2V_Vec,
    // vesLinkSpecPairR2V_Vec, vesLinkSpecR2V_Vec data to VesRaft
    MPI_GatherVec<VesProxyR2V>(vesProxyR2V_Vec,
                               dataTypeUtil.MPI_VesProxyR2V,
                               vesraftRank_World,
                               myRank_World,
                               nHosts_World,
                               MPI_COMM_WORLD);
    MPI_GatherVec<VesSurfSpecR2V>(vesSurfSpecR2V_Vec,
                                  dataTypeUtil.MPI_VesSurfSpecR2V,
                                  vesraftRank_World,
                                  myRank_World,
                                  nHosts_World,
                                  MPI_COMM_WORLD);
    MPI_GatherVec<VesInnerSpecR2V>(vesInnerSpecR2V_Vec,
                                   dataTypeUtil.MPI_VesInnerSpecR2V,
                                   vesraftRank_World,
                                   myRank_World,
                                   nHosts_World,
                                   MPI_COMM_WORLD);
    MPI_GatherVec<VesLinkSpecPairR2V>(vesLinkSpecPairR2V_Vec,
                                      dataTypeUtil.MPI_VesLinkSpecPairR2V,
                                      vesraftRank_World,
                                      myRank_World,
                                      nHosts_World,
                                      MPI_COMM_WORLD);
    MPI_GatherVec<VesLinkSpecR2V>(vesLinkSpecR2V_Vec,
                                  dataTypeUtil.MPI_VesLinkSpecR2V,
                                  vesraftRank_World,
                                  myRank_World,
                                  nHosts_World,
                                  MPI_COMM_WORLD);

    // PointSpec IDs will be reassigned in VesRaft, so we can reset here
    pNextPointSpecUniqueID.set(std::numeric_limits<index_t>::max() * 0.9);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_useRaftV2R() {
    // TODO RECEIVE FLAG FROM MASTER AND ONLY DO ALL THIS IF NECESSARY ()

    // receive raftProxyV2R_Vec and raftSpecV2R_Vec data from VesRaft
    MPI_BcastVec<RaftProxyV2R>(raftProxyV2R_Vec,
                               dataTypeUtil.MPI_RaftProxyV2R,
                               vesraftRank_World,
                               myRank_World,
                               MPI_COMM_WORLD);
    MPI_BcastVec<RaftSpecV2R>(raftSpecV2R_Vec,
                              dataTypeUtil.MPI_RaftSpecV2R,
                              vesraftRank_World,
                              myRank_World,
                              MPI_COMM_WORLD);

    MPI_BcastVec<RaftSReacInactiveV2R>(raftSReacInactiveV2R_Vec,
                                       dataTypeUtil.MPI_RaftSReacInactiveV2R,
                                       vesraftRank_World,
                                       myRank_World,
                                       MPI_COMM_WORLD);

    for (auto const& tri: pTris) {
        if (tri == nullptr or !tri->getInHost()) {
            continue;
        }
        tri->clearRaftProxyrefs();
    }

    for (auto const& raft_proxy_v2r: raftProxyV2R_Vec) {
        const triangle_global_id tri_gidx(raft_proxy_v2r.triangle_global_index);

        if (triHosts.find(tri_gidx) == triHosts.end()) {
            ProgErrLog("Triangle not assigned host.\n");
        }

        TriRDEF* tri = _tri(tri_gidx);

        if (tri == nullptr) {
            ProgErrLog("Triangle not assigned to a patch.\n");
        }

        if (!tri->getInHost()) {
            continue;
        }

        const solver::raft_global_id raft_gidx(raft_proxy_v2r.raft_global_index);
        auto& raft_def = statedef().raftdef(raft_gidx);

        const solver::raft_individual_id raft_individual_id(raft_proxy_v2r.raft_individual_index);

        tri->createRaftProxyref(&raft_def, raft_individual_id);
    }

    for (auto const& raft_spec_v2r: raftSpecV2R_Vec) {
        const triangle_global_id tri_gidx(raft_spec_v2r.triangle_global_index);

        if (triHosts.find(tri_gidx) == triHosts.end()) {
            ProgErrLog("Triangle not assigned host.\n");
        }

        TriRDEF* tri = _tri(tri_gidx);

        if (tri == nullptr) {
            ProgErrLog("Triangle not assigned to a patch.\n");
        }

        if (!tri->getInHost()) {
            continue;
        }

        const solver::raft_individual_id raft_individual_id(raft_spec_v2r.raft_individual_index);

        const solver::spec_global_id spec_gidx(raft_spec_v2r.spec_global_index);
        RaftProxy* raftproxy = tri->getRaftProxyref(raft_individual_id);

        raftproxy->setSpecCountByGidx(spec_gidx, raft_spec_v2r.count);
    }


    for (auto const& raft_srinactive_v2r: raftSReacInactiveV2R_Vec) {
        const triangle_global_id tri_gidx(raft_srinactive_v2r.triangle_global_index);

        if (triHosts.find(tri_gidx) == triHosts.end()) {
            ProgErrLog("Triangle not assigned host.\n");
        }

        TriRDEF* tri = _tri(tri_gidx);

        if (tri == nullptr) {
            ProgErrLog("Triangle not assigned to a patch.\n");
        }

        if (!tri->getInHost()) {
            continue;
        }

        const solver::raft_individual_id raft_individual_id(
            raft_srinactive_v2r.raft_individual_index);

        const solver::raftsreac_global_id raftsreac_gidx(
            raft_srinactive_v2r.raftsreac_global_index);
        RaftProxy* raftproxy = tri->getRaftProxyref(raft_individual_id);

        raftproxy->setRaftSReacInActive(raftsreac_gidx);
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::_constructRaftR2V() {
    raftProxyR2V_Vec.clear();
    raftSpecR2V_Vec.clear();
    triRaftGenR2V_Vec.clear();

    for (auto const& tri: pTris) {
        if (tri == nullptr or !tri->getInHost()) {
            continue;
        }

        for (auto const& raftgen: tri->getAppliedRaftGens()) {
            triRaftGenR2V_Vec.push_back({tri->idx(), raftgen.first, raftgen.second});
        }
        tri->clearAppliedRaftGens();


        for (auto const& raft: tri->getRaftProxyrefs()) {
            RaftProxyR2V raft_proxy_r2v{};

            // solver::raft_global_id raft_gidx = raft->idx();

            raft_proxy_r2v.triangle_global_index = tri->idx();
            raft_proxy_r2v.raft_individual_index = raft.first;
            raft_proxy_r2v.immobility_update = raft.second->getImmobilityUpdate();

            raftProxyR2V_Vec.emplace_back(raft_proxy_r2v);

            for (auto const& specs: raft.second->getSpecs()) {
                RaftSpecR2V raft_spec_r2v{};
                raft_spec_r2v.triangle_global_index = tri->idx();
                raft_spec_r2v.raft_individual_index = raft.first;
                raft_spec_r2v.spec_global_index.set(specs.first);
                raft_spec_r2v.count = specs.second;

                raftSpecR2V_Vec.emplace_back(raft_spec_r2v);
            }
        }
    }

    MPI_GatherVec<RaftGenCountR2V>(triRaftGenR2V_Vec,
                                   dataTypeUtil.MPI_RaftGenCountR2V,
                                   vesraftRank_World,
                                   myRank_World,
                                   nHosts_World,
                                   MPI_COMM_WORLD);

    // send raftProxyR2V_Vec and raftSpecR2V_Vec data to VesRaft
    MPI_GatherVec<RaftProxyR2V>(raftProxyR2V_Vec,
                                dataTypeUtil.MPI_RaftProxyR2V,
                                vesraftRank_World,
                                myRank_World,
                                nHosts_World,
                                MPI_COMM_WORLD);
    MPI_GatherVec<RaftSpecR2V>(raftSpecR2V_Vec,
                               dataTypeUtil.MPI_RaftSpecR2V,
                               vesraftRank_World,
                               myRank_World,
                               nHosts_World,
                               MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////

void TetVesicleRDEF::setOutputSync(bool enable_sync, int output_rank) noexcept {
    syncOutput = enable_sync;
    outputRank = output_rank;
}

////////////////////////////////////////////////////////////////////////////

bool TetVesicleRDEF::getOutputSyncStatus() const noexcept {
    return syncOutput;
}

////////////////////////////////////////////////////////////////////////////

int TetVesicleRDEF::getOutputSyncRank() const noexcept {
    return outputRank;
}

}  // namespace steps::mpi::tetvesicle
