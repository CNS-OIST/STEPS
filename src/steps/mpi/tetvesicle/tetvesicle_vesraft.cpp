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

#include "mpi/tetvesicle/tetvesicle_vesraft.hpp"

// Standard library & STL headers.
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <vector>

#include <mpi.h>

// STEPS headers.
#include "mpi/mpi_common.hpp"
#include "util/common.hpp"

#include "mpi/tetvesicle/comp_vesraft.hpp"
#include "mpi/tetvesicle/diff.hpp"
#include "mpi/tetvesicle/diffboundary.hpp"
#include "mpi/tetvesicle/endocytosis.hpp"
#include "mpi/tetvesicle/exocytosis.hpp"
#include "mpi/tetvesicle/ghkcurr.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/linkspec.hpp"
#include "mpi/tetvesicle/patch_vesraft.hpp"
#include "mpi/tetvesicle/raftendocytosis.hpp"
#include "mpi/tetvesicle/reac.hpp"
#include "mpi/tetvesicle/sdiff.hpp"
#include "mpi/tetvesicle/sreac.hpp"
#include "mpi/tetvesicle/tet_vesraft.hpp"
#include "mpi/tetvesicle/tri_vesraft.hpp"
#include "mpi/tetvesicle/vdepsreac.hpp"

#include "solver/chandef.hpp"
#include "solver/compdef.hpp"
#include "solver/diffdef.hpp"
#include "solver/exocytosisdef.hpp"
#include "solver/ghkcurrdef.hpp"
#include "solver/linkspecdef.hpp"
#include "solver/ohmiccurrdef.hpp"
#include "solver/patchdef.hpp"
#include "solver/raftendocytosisdef.hpp"
#include "solver/reacdef.hpp"
#include "solver/sreacdef.hpp"
#include "solver/statedef.hpp"
#include "solver/types.hpp"
#include "solver/vdepsreacdef.hpp"
#include "solver/vessdiffdef.hpp"
#include "solver/vessreacdef.hpp"
#include "solver/vesunbinddef.hpp"

#include "solver/efield/dVsolver.hpp"
#include "solver/efield/efield.hpp"
#ifdef USE_PETSC
#include "solver/efield/dVsolver_petsc.hpp"
#endif

#include "geom/tetmesh.hpp"
#include "math/constants.hpp"
#include "math/point.hpp"
#include "math/tetrahedron.hpp"
#include "util/checkpointing.hpp"
#include "util/collections.hpp"
#include "util/distribute.hpp"
#include "util/error.hpp"

#include <fau.de/overlap.hpp>

#include <time.h>

namespace steps::mpi::tetvesicle {

using math::point3d;

TetVesicleVesRaft::TetVesicleVesRaft(model::Model* m,
                                     wm::Geom* g,
                                     const rng::RNGptr& r,
                                     int /*calcMembPot*/)
    : API(m, g, r)
    , pMesh(nullptr)
    , pRequireVesicleCommunication(false)
    , pQtablesize_spec(1000)
    , pQtablesize_linkspec(1000) {
    if (rng() == nullptr) {
        std::ostringstream os;
        os << "No RNG provided to solver initializer function";
        ArgErrLog(os.str());
    }

    // All initialization code now in _setup() to allow EField solver to be
    // derived and create EField local objects within the constructor

    std::cout.setf(std::ios::unitbuf);

    pMesh = dynamic_cast<tetmesh::Tetmesh*>(geom());
    if (pMesh == nullptr) {
        ArgErrLog(
            "Geometry description to solver::TetVesicleVesRaft solver "
            "constructor is not a valid tetmesh::Tetmesh object.");
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank_World);
    MPI_Comm_size(MPI_COMM_WORLD, &nHosts_World);

    _partition();
    MPI_Barrier(MPI_COMM_WORLD);
    dataTypeUtil.commitAllDataTypes();
    _setup();
    // may remove this barrier later.
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

TetVesicleVesRaft::~TetVesicleVesRaft() {
    for (auto c: pComps) {
        delete c;
    }
    for (auto p: pPatches) {
        delete p;
    }

    for (auto t: pTets) {
        delete t;
    }
    for (auto t: pTet_ext) {
        delete t;
    }
    for (auto t: pTris) {
        delete t;
    }

    for (auto& qtit: pQtables_spec) {
        for (auto& q: qtit.second) {
            if (q != nullptr) {
                delete q;
            }
        }
        qtit.second.container().clear();
    }

    for (auto& qtit: pQtables_linkspec) {
        for (auto& q: qtit.second) {
            if (q != nullptr) {
                delete q;
            }
        }
        qtit.second.container().clear();
    }

    for (auto lsp: pLinkSpecPairs) {
        delete lsp;
    }

    dataTypeUtil.freeAllDataTypes();
}

///////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::checkpoint(std::string const& file_name) {
    std::fstream cp_file;

    cp_file.open(file_name.c_str(), std::fstream::out | std::fstream::binary | std::fstream::trunc);
    util::checkpoint(cp_file, myRank_World);

    API::checkpoint(cp_file);

    // Basic data
    util::checkpoint(cp_file, pTemp);
    util::checkpoint(cp_file, pEFDT);
    util::checkpoint(cp_file, pVesicledt);
    util::checkpoint(cp_file, pVesicleDefaultdt);
    util::checkpoint(cp_file, pVesicles_count);
    util::checkpoint(cp_file, pRafts_count);
    util::checkpoint(cp_file, pNextPointSpecUniqueID);
    util::checkpoint(cp_file, maxWalkDistSqFact);
    util::checkpoint(cp_file, minNbTetVisited);
    util::checkpoint(cp_file, pNextLinkSpecUniqueID);
    util::checkpoint(cp_file, pVesSpecD);
    util::checkpoint(cp_file, pVesLinkSpecD);

    statedef().checkpoint(cp_file);

    for (auto const& c: pComps) {
        c->checkpoint(cp_file);
    }

    for (auto const& p: pPatches) {
        p->checkpoint(cp_file);
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

    uint lsp_size = pLinkSpecPairs.size();
    util::checkpoint(cp_file, lsp_size);
    for (auto const& lsp: pLinkSpecPairs) {
        lsp->checkpoint(cp_file);
    }

    // Note: Don't need to do anything for pVesicles and pRafts because Comp and Patch own them

    for (auto const& ves_to_qt: pQtables_spec) {
        std::map<solver::spec_global_id, double> spec_tau;
        for (auto spec_gidx: solver::spec_global_id::range(ves_to_qt.second.size())) {
            auto qt = ves_to_qt.second[spec_gidx];
            if (qt != nullptr) {
                spec_tau[spec_gidx] = qt->getTau();
            }
        }
        util::checkpoint(cp_file, spec_tau);
    }

    for (auto const& ves_to_qt: pQtables_linkspec) {
        std::map<solver::linkspec_global_id, double> linkspec_tau;
        for (auto linkspec_gidx: solver::linkspec_global_id::range(ves_to_qt.second.size())) {
            auto qt = ves_to_qt.second[linkspec_gidx];
            if (qt != nullptr) {
                linkspec_tau[linkspec_gidx] = qt->getTau();
            }
        }
        util::checkpoint(cp_file, linkspec_tau);
    }

    std::vector<std::string> path_ids;
    for (auto const& path: pPaths) {
        path_ids.emplace_back(path.first);
    }
    util::checkpoint(cp_file, path_ids);

    for (auto const& path: pPaths) {
        path.second->checkpoint(cp_file);
    }

    cp_file.close();
}

///////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::restore(std::string const& file_name) {
    std::fstream cp_file;

    cp_file.open(file_name.c_str(), std::fstream::in | std::fstream::binary);

    cp_file.seekg(0);
    util::compare(cp_file, myRank_World);

    API::restore(cp_file);

    util::restore(cp_file, pTemp);
    util::restore(cp_file, pEFDT);
    util::restore(cp_file, pVesicledt);
    util::compare(cp_file, pVesicleDefaultdt);
    util::restore(cp_file, pVesicles_count);
    util::restore(cp_file, pRafts_count);
    util::restore(cp_file, pNextPointSpecUniqueID);
    util::restore(cp_file, maxWalkDistSqFact);
    util::restore(cp_file, minNbTetVisited);
    util::restore(cp_file, pNextLinkSpecUniqueID);
    util::restore(cp_file, pVesSpecD);
    util::restore(cp_file, pVesLinkSpecD);

    statedef().restore(cp_file);

    for (auto const& c: pComps) {
        c->restore(cp_file);
    }

    for (auto const& p: pPatches) {
        p->restore(cp_file);
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

    uint lsp_size;
    util::restore(cp_file, lsp_size);

    for (uint i = 0; i < lsp_size; ++i) {
        auto linkspecpair = new LinkSpecPair(cp_file, this);  // constructor does the restore
        pLinkSpecPairs.insert(linkspecpair);
    }

    for (auto ves_gidx: solver::vesicle_global_id::range(statedef().countVesicles())) {
        std::map<solver::spec_global_id, double> spec_tau;
        util::restore(cp_file, spec_tau);
        for (auto const& st: spec_tau) {
            pQtables_spec[ves_gidx][st.first] = new Qtable(pQtablesize_spec, st.second, rng());
        }
    }

    for (auto ves_gidx: solver::vesicle_global_id::range(statedef().countVesicles())) {
        std::map<solver::linkspec_global_id, double> linkspec_tau;
        util::restore(cp_file, linkspec_tau);
        for (auto const& lst: linkspec_tau) {
            pQtables_linkspec[ves_gidx][lst.first] =
                new Qtable(pQtablesize_linkspec, lst.second, rng());
        }
    }

    std::vector<std::string> path_ids;
    util::restore(cp_file, path_ids);
    for (auto const& path_id: path_ids) {
        Path* path = new Path(path_id);
        pPaths[path_id] = path;
        path->restore(cp_file);
    }

    cp_file.close();

    // Sync with RDEF cores
    _constructVesV2R();
    _constructRaftV2R();
    pRequireVesicleCommunication = false;
}

////////////////////////////////////////////////////////////////////////////////

std::string TetVesicleVesRaft::getSolverName() const {
    return "TetVesicleVesRaft";
}

////////////////////////////////////////////////////////////////////////////////

std::string TetVesicleVesRaft::getSolverDesc() const {
    return "Parallel approximate stochastic solver with added Vesicle-related "
           "functionality.";
}

////////////////////////////////////////////////////////////////////////////////

std::string TetVesicleVesRaft::getSolverAuthors() const {
    return "Iain Hepburn, Weiliang Chen, Stefan Wils, Blue Brain Project.";
}

////////////////////////////////////////////////////////////////////////////////

std::string TetVesicleVesRaft::getSolverEmail() const {
    return "steps.dev@gmail.com";
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_partition() {
    if (myRank_World != 0) {
        ProgErr("A TetVesicleVesRaft solver is created in a VesRDEF rank.");
    } else {
        vesraftRank_World = 0;
        RDEFmasterRank_World = 1;
        MPI_Comm single_node_comm;
        MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &single_node_comm);
    }

    uint n_tets_tris[2];
    std::vector<tetrahedron_global_id> tet_ids;
    std::vector<int> tet_hosts;

    std::vector<triangle_global_id> tri_ids;
    std::vector<int> tri_hosts;

    pMesh->_autoPartition(1, 1, nHosts_World - 1, tetHosts, triHosts);
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

    MPI_Bcast(n_tets_tris, 2, MPI_UNSIGNED, vesraftRank_World, MPI_COMM_WORLD);

    MPI_Bcast(tet_ids.data(), n_tets_tris[0], MPI_STEPS_INDEX, vesraftRank_World, MPI_COMM_WORLD);
    MPI_Bcast(tet_hosts.data(), n_tets_tris[0], MPI_INT, vesraftRank_World, MPI_COMM_WORLD);
    MPI_Bcast(tri_ids.data(), n_tets_tris[1], MPI_STEPS_INDEX, vesraftRank_World, MPI_COMM_WORLD);
    MPI_Bcast(tri_hosts.data(), n_tets_tris[1], MPI_INT, vesraftRank_World, MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setup() {
    // First initialise the pTets, pTris vector, because
    // want tets and tris to maintain indexing from Geometry

    uint ntets = pMesh->countTets();
    uint ntris = pMesh->countTris();

    pTets.container().assign(ntets, nullptr);
    pTris.container().assign(ntris, nullptr);

    pTet_ext.container().assign(ntets, nullptr);

    // Now create the actual compartments.
    for (auto const& c: statedef().comps()) {
        solver::comp_global_id compdef_gidx = c->gidx();
        uint comp_idx = _addComp(c, pMesh);
        AssertLog(compdef_gidx == comp_idx);
    }
    // Create the actual patches.
    for (auto const& p: statedef().patches()) {
        solver::patch_global_id patchdef_gidx = p->gidx();
        uint patch_idx = _addPatch(p, pMesh);
        AssertLog(patchdef_gidx == patch_idx);
    }

    uint npatches = pPatches.size();
    AssertLog(pMesh->_countPatches() == npatches);

    for (auto p: solver::patch_global_id::range(npatches)) {
        // Add the tris for this patch
        // We have checked the indexing - p is the global index
        wm::Patch* wmpatch = pMesh->_getPatch(p);

        // Perform upcast
        auto* tmpatch = dynamic_cast<tetmesh::TmPatch*>(wmpatch);
        if (tmpatch == nullptr) {
            ArgErrLog(
                "Well-mixed patches not supported in "
                "solver::TetVesicleVesRaft solver.");
        }
        PatchVesRaft* localpatch = pPatches[p];

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

            const point3d& baryc = pMesh->_getTriBarycenter(tri);
            const point3d& trinorm = pMesh->_getTriNorm(tri);

            point3d d;
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
        wm::Comp* wmcomp = pMesh->_getComp(c);

        // Perform upcast
        auto* tmcomp = dynamic_cast<tetmesh::TmComp*>(wmcomp);
        if (tmcomp != nullptr) {
            CompVesRaft* localcomp = pComps[c];

            for (auto const& tet: tmcomp->_getAllTetIndices()) {
                AssertLog(pMesh->getTetComp(tet) == tmcomp);

                double vol = pMesh->getTetVol(tet);

                const auto& tris = pMesh->_getTetTriNeighb(tet);

                double a[4] = {0, 0, 0, 0};
                for (uint j = 0; j < 4; ++j) {
                    a[j] = pMesh->getTriArea(tris[j]);
                }

                const auto tets = pMesh->_getTetTetNeighb(tet);
                point3d baryc = pMesh->_getTetBarycenter(tet);

                double d[4] = {0, 0, 0, 0};
                for (uint j = 0; j < 4; ++j) {
                    if (tets[j] == tetrahedron_global_id::unknown_value()) {
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
                        tets[3],
                        baryc);
            }

            // SETUP VESICLES

            // Should be safe here to add the vesicles
            // All compartments can hold all types of vesicle??
            localcomp->setupVesicles();

            // for Qtables.
            uint nSpecs_global = statedef().countSpecs();
            uint nLinkSpecs_global = statedef().countLinkSpecs();

            for (auto const& ves: statedef().vesicles()) {
                solver::vesicle_global_id vidx = ves->gidx();

                // Set up the surface diffusion stuff
                pQtables_spec[vidx] = std::vector<Qtable*>(nSpecs_global, nullptr);
                pVesSpecD[vidx] = std::vector<double>(nSpecs_global, 0.0);

                pQtables_linkspec[vidx] = std::vector<Qtable*>(nLinkSpecs_global, nullptr);
                pVesLinkSpecD[vidx] = std::vector<double>(nLinkSpecs_global, 0.0);

                for (auto vsd_idx: solver::vessdiff_local_id::range(ves->countVesSurfDiffs())) {
                    solver::VesSDiffdef* vsddef = ves->vessurfdiffdef(vsd_idx);
                    AssertLog(vsddef != nullptr);
                    setVesicleSpecDiffD_(vidx, vsddef->lig(), vsddef->dcst());
                }

                for (auto const& bs: statedef().linkspecs()) {
                    double dcst = bs->dcst();
                    if (dcst > 0.0) {
                        _setVesicleSurfaceLinkSpecSDiffD(vidx, bs->gidx(), dcst);
                    }
                }
            }
        } else {
            ProgErrLog("Well-mixed compartments not supported for vesicle simulations.");
        }
    }

    // All tets and tris that belong to some comp or patch have been created
    // locally- now we can connect them locally
    // NOTE: currently if a tetrahedron's neighbour belongs to a different
    // comp they do not talk to each other (see Tet::setNextTet())

    AssertLog(ntets == pTets.size());

    // pTets member size of all tets in geometry, but may not be filled with
    // local tets if they have not been added to a compartment
    for (auto t: tetrahedron_global_id::range(ntets)) {
        if (pTets[t] == nullptr) {
            continue;
        }

        for (uint j = 0; j < 4; ++j) {
            tetrahedron_global_id tet = pTets[t]->tet(j);
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
            triangle_global_id tri = pTris[t]->tri(j);
            if (tri.valid() && pTris[tri] != nullptr) {
                pTris[t]->setNextTri(j, pTris[tri]);
            }
        }

        // By convention, triangles in a patch should have an inner tetrahedron
        // defined (neighbouring tets 'flipped' if necessary in Tetmesh) but not
        // necessarily an outer tet 17/3/10- actually this is not the case any more
        // with well-mixed compartments

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
                    TetVesRaft* tet_in = pTets[tetinner];
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
                    TetVesRaft* tet_out = pTets[tetouter];

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

    // Need to setup Rafts here after the triangles have been setup, because
    // Patches are going to create lookup tables depending on the connectivity
    for (auto const& patch: pPatches) {
        // All patches can contain all rafts, for now
        patch->setupRafts();

        patch->setupEndocyticZones();
    }

    // Fill linkspecpair -> vesunbind map
    for (auto const& vub: statedef().vesunbinds()) {
        auto& val = pLinkSpecPair2VesUnbinds[std::make_tuple(vub->getVes1idx(),
                                                             vub->getVes2idx(),
                                                             vub->getLinkSpec1gidx(),
                                                             vub->getLinkSpec2gidx())];
        val.first.push_back(vub);
        val.second += vub->kcst();
    }

    // Note, this must match up with the intitial pNextPointSpecUniqueID on RDEF ranks
    pRDEFminPointSpecUniqueID.set(std::numeric_limits<index_t>::max() * 0.9);
}

////////////////////////////////////////////////////////////////////////////////

// 'Safe' global to local index translation methods that throw on error.

inline solver::spec_local_id TetVesicleVesRaft::_specG2L_or_throw(
    CompVesRaft* comp,
    solver::spec_global_id gidx) const {
    AssertLog(gidx < statedef().countSpecs());
    solver::spec_local_id lidx = comp->def()->specG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("species undefined in compartment");
    }
    return lidx;
}

inline solver::spec_local_id TetVesicleVesRaft::_specG2L_or_throw(
    PatchVesRaft* patch,
    solver::spec_global_id gidx) const {
    AssertLog(gidx < statedef().countSpecs());
    solver::spec_local_id lidx = patch->def()->specG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("species undefined in patch");
    }
    return lidx;
}

#if 0
inline solver::spec_local_id TetVesicleVesRaft::_specG2L_or_throw(TetVesRaft *tet, solver::spec_global_id gidx) const {
    AssertLog(gidx < statedef().countSpecs());
    solver::spec_local_id lidx = tet->compdef()->specG2L(gidx);

    if (lidx.unknown()) {
       ArgErrLog("species undefined in tetrahedron");
    }
    return lidx;
}

inline solver::spec_local_id TetVesicleVesRaft::_specG2L_or_throw(Tri *tri, solver::spec_global_id gidx) const {
    AssertLog(gidx < statedef().countSpecs());
    solver::spec_local_id lidx = tri->patchdef()->SPEC_LIDX_UNDEFINED(gidx);

    if (lidx.unknown()) {
        ArgErrLog("species undefined in triangle");
    }
    return lidx;
}
#endif

inline solver::reac_local_id TetVesicleVesRaft::_reacG2L_or_throw(
    CompVesRaft* comp,
    solver::reac_global_id gidx) const {
    AssertLog(gidx < statedef().countReacs());
    solver::reac_local_id lidx = comp->def()->reacG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("reaction undefined in compartment");
    }
    return lidx;
}

inline solver::sreac_local_id TetVesicleVesRaft::_sreacG2L_or_throw(
    PatchVesRaft* patch,
    solver::sreac_global_id gidx) const {
    AssertLog(gidx < statedef().countSReacs());
    solver::sreac_local_id lidx = patch->def()->sreacG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("surface reaction undefined in patch");
    }
    return lidx;
}

inline solver::diff_local_id TetVesicleVesRaft::_diffG2L_or_throw(
    CompVesRaft* comp,
    solver::diff_global_id gidx) const {
    AssertLog(gidx < statedef().countDiffs());
    solver::diff_local_id lidx = comp->def()->diffG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("diffusion rule undefined in compartment");
    }
    return lidx;
}

inline solver::endocytosis_local_id TetVesicleVesRaft::_endoG2L_or_throw(
    PatchVesRaft* patch,
    solver::endocytosis_global_id gidx) const {
    AssertLog(gidx < statedef().countEndocytosis());

    solver::endocytosis_local_id lidx = patch->def()->endocytosisG2L(gidx);

    if (lidx.unknown()) {
        ArgErrLog("endocytosis undefined in patch");
    }
    return lidx;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_addComp(solver::Compdef* cdef, tetmesh::Tetmesh* mesh) {
    auto* comp = new CompVesRaft(cdef, mesh, this);
    AssertLog(comp != nullptr);
    uint compidx = pComps.size();
    pComps.container().push_back(comp);
    return compidx;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_addPatch(solver::Patchdef* pdef, tetmesh::Tetmesh* mesh) {
    auto* patch = new PatchVesRaft(pdef, mesh, this);
    AssertLog(patch != nullptr);
    uint patchidx = pPatches.size();
    pPatches.container().push_back(patch);
    return patchidx;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_addTet(tetrahedron_global_id tetidx,
                                CompVesRaft* comp,
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
                                math::point3d baryc) {
    solver::Compdef* compdef = comp->def();
    auto* localtet = new TetVesRaft(
        tetidx, compdef, vol, a1, a2, a3, a4, d1, d2, d3, d4, tet0, tet1, tet2, tet3, baryc);
    AssertLog(localtet != nullptr);
    AssertLog(tetidx.get() < pTets.size());
    AssertLog(pTets[tetidx] == nullptr);
    pTets[tetidx] = localtet;
    comp->addTet(localtet);

    // Setup tetrahedrons for overlap library

    const auto tetverts = _mesh()->getTet(tetidx);

    std::vector<double> vert0 = _mesh()->getVertex(vertex_id_t(tetverts[0]));
    vector_t vertex0(vert0[0], vert0[1], vert0[2]);

    std::vector<double> vert1 = _mesh()->getVertex(vertex_id_t(tetverts[1]));
    vector_t vertex1(vert1[0], vert1[1], vert1[2]);

    std::vector<double> vert2 = _mesh()->getVertex(vertex_id_t(tetverts[2]));
    vector_t vertex2(vert2[0], vert2[1], vert2[2]);

    std::vector<double> vert3 = _mesh()->getVertex(vertex_id_t(tetverts[3]));
    vector_t vertex3(vert3[0], vert3[1], vert3[2]);

    // Buggy if I replace the following line with the stuff commented out after
    // it. Problem with tet_anticlockwise??
    auto* tex = new Tetrahedron(vertex0, vertex1, vertex2, vertex3);
    /*
    math::point3d v0 = math::point3d(vert0[0], vert0[1], vert0[2]);
    math::point3d v1 = math::point3d(vert1[0], vert1[1], vert1[2]);
    math::point3d v2 = math::point3d(vert2[0], vert2[1], vert2[2]);
    math::point3d v3 = math::point3d(vert3[0], vert3[1], vert3[2]);

    Tetrahedron * tex;

    if (tet_anticlockwise(v0, v1, v2,  v3))
    {
        tex = new Tetrahedron(vertex0, vertex1, vertex2,  vertex3);
    }
    else
    {
        AssertLog(tet_anticlockwise(v0, v2, v1, v3));
        tex = new Tetrahedron(vertex0, vertex2, vertex1,  vertex3);
        }
        */

    pTet_ext[tetidx] = tex;

    // MPISTEPS

    localtet->setSolverVesRaft(this);
    // MPISTEPS
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_addTri(triangle_global_id triidx,
                                PatchVesRaft* patch,
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
                                point3d baryc,
                                point3d trinorm) {
    solver::Patchdef* patchdef = patch->def();
    auto* tri = new TriVesRaft(triidx,
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
                               trinorm);
    AssertLog(tri != nullptr);
    AssertLog(triidx.get() < pTris.size());
    AssertLog(pTris[triidx] == nullptr);
    pTris[triidx] = tri;
    patch->addTri(tri);

    // MPISTEPS
    tri->setSolverVesRaft(this);
    // MPISTEPS
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::reset() {
    for (auto const& comp: pComps) {
        comp->reset();
    }

    for (auto const& patch: pPatches) {
        patch->reset();
    }

    for (auto const& tet: pTets) {
        if (tet == nullptr) {
            continue;
        }
        tet->reset();
    }

    for (auto const& t: pTris) {
        if (t == nullptr) {
            continue;
        }
        t->reset();
    }

    statedef().resetTime();
    statedef().resetNSteps();

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

    for (auto const& lsp: pLinkSpecPairs) {
        delete lsp;
    }
    pLinkSpecPairs.clear();

    pVesicles_count = 0;
    pRafts_count = 0;

    pVesicledt = pVesicleDefaultdt;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::run(double endtime) {
    ArgErrLogIf(endtime < statedef().time(), "Endtime is before current simulation time ");

    //  only do these dramatic updates if necessary
    MPI_Bcast(&pRequireVesicleCommunication, 1, MPI_C_BOOL, vesraftRank_World, MPI_COMM_WORLD);

    if (pRequireVesicleCommunication) {
        _constructVesV2R();
        _constructRaftV2R();
        pRequireVesicleCommunication = false;
    }

    while (true) {
        if (util::almost_equal(statedef().time(), endtime)) {
            // We're done. Need to tell RDEF with a negative number
            double endtime_reached = -1.0;
            MPI_Bcast(&endtime_reached, 1, MPI_DOUBLE, vesraftRank_World, MPI_COMM_WORLD);

            statedef().setTime(endtime);
            return;
        } else {
            double ves_endtime = statedef().time() + pVesicledt;
            double vesicle_dt = pVesicledt;

            if (ves_endtime > endtime) {
                ves_endtime = endtime;
                vesicle_dt = endtime - statedef().time();
            }

            MPI_Bcast(&ves_endtime, 1, MPI_DOUBLE, vesraftRank_World, MPI_COMM_WORLD);

            // RDEF does its reaction-diffusion business, then sends updated vesicle
            // info
            _useVesR2V();
            _useRaftR2V();

            _syncPools(RDEF_TO_VESRAFT);  // _syncPools take care of clearing vectors

            _runVesicle(vesicle_dt);  // May add to tetPoolCountSyncs_Vec
            _runRaft(vesicle_dt);

            statedef().setTime(statedef().time() + vesicle_dt);

            _constructVesV2R();
            _constructRaftV2R();

            // Keep all tri sync for now.
            for (auto const& tri: pTris) {
                if (tri == nullptr) {
                    continue;
                }
                for (auto spec_lidx: solver::spec_local_id::range(tri->patchdef()->countSpecs())) {
                    regTriPoolSync_(tri->idx(),
                                    tri->patchdef()->specL2G(spec_lidx),
                                    tri->pools()[spec_lidx]);
                }
            }

            _syncPools(VESRAFT_TO_RDEF);  // contains a subset for tets set by any exocytosis events
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_runVesicle(double dt) {
    // New routine to apply vesicle unbinding events. Do first because a)
    // unbinding can be fast and b) bnding has already had a chance

    // First store any linkspecpairs to go
    std::vector<std::set<const LinkSpecPair*, util::DerefPtrLess<LinkSpecPair>>::iterator>
        linkspecpair_del;

    for (auto lspit = pLinkSpecPairs.begin(); lspit != pLinkSpecPairs.end(); ++lspit) {
        auto const* linkspecpair = *lspit;
        auto linkspec1 = linkspecpair->getLinkSpec1();
        auto linkspec2 = linkspecpair->getLinkSpec2();

        auto pair_linkspec1_gidx = linkspec1->getGidx();
        auto pair_linkspec2_gidx = linkspec2->getGidx();

        auto vesicle1 = linkspecpair->getVesicle1();
        auto vesicle2 = linkspecpair->getVesicle2();

        auto pair_vesicle1_gidx = vesicle1->idx();
        auto pair_vesicle2_gidx = vesicle2->idx();

        solver::VesUnbinddef* selectedUnbind = nullptr;
        double total_rate = 0.0;

        // First check if (A,B) reactions exist
        auto linkSpecKey1 = std::make_tuple(pair_vesicle1_gidx,
                                            pair_vesicle2_gidx,
                                            pair_linkspec1_gidx,
                                            pair_linkspec2_gidx);
        auto it1 = pLinkSpecPair2VesUnbinds.find(linkSpecKey1);
        if (it1 != pLinkSpecPair2VesUnbinds.end()) {
            total_rate += it1->second.second;
        }

        auto it2 = pLinkSpecPair2VesUnbinds.end();

        bool isSymmetric = linkspecpair->isSymmetric();

        if (not isSymmetric) {
            // If the linkspec pair is asymmetrical, we need to check the existence of (B,A)
            // reactions
            auto linkSpecKey2 = std::make_tuple(linkspecpair->getVesicle2()->idx(),
                                                linkspecpair->getVesicle1()->idx(),
                                                linkspecpair->getLinkSpec2()->getGidx(),
                                                linkspecpair->getLinkSpec1()->getGidx());
            it2 = pLinkSpecPair2VesUnbinds.find(linkSpecKey2);
            if (it2 != pLinkSpecPair2VesUnbinds.end()) {
                total_rate += it2->second.second;
            }
        }

        if (total_rate > 0.0) {
            double vub_dt = rng()->getExp(total_rate);
            if (vub_dt < dt) {
                double rand = rng()->getUnfIE() * total_rate;
                double acc = 0;
                if (it1 != pLinkSpecPair2VesUnbinds.end()) {
                    for (auto& vub: it1->second.first) {
                        acc += vub->kcst();
                        if (rand < acc) {
                            selectedUnbind = vub;
                            if (isSymmetric and ((rng()->get() % 2) != 0)) {
                                // Although the LHS is symmetric, the RHS might not be, so we need
                                // to swap with 50% prob.
                                std::swap(vesicle1, vesicle2);
                                std::swap(linkspec1, linkspec2);
                            }
                            break;
                        }
                    }
                }
                if (selectedUnbind == nullptr) {
                    // If the vesunbind reaction was not selected among (A,B) reactions
                    AssertLog(it2 != pLinkSpecPair2VesUnbinds.end());
                    for (auto& vub: it2->second.first) {
                        acc += vub->kcst();
                        if (rand < acc) {
                            selectedUnbind = vub;
                            // If the selected vesunbind reaction corresponds to the reversed case,
                            // swap vesicles and linkspecs.
                            std::swap(vesicle1, vesicle2);
                            std::swap(linkspec1, linkspec2);
                            break;
                        }
                    }
                    AssertLog(selectedUnbind != nullptr);
                }

                // Apply vesunbind reaction
                auto vub_spec1_gidx = selectedUnbind->getSpec1gidx();
                auto vub_spec2_gidx = selectedUnbind->getSpec2gidx();

                math::position_rel_to_ves linkspec1_pos = linkspec1->getPosCartesian_rel();
                math::position_rel_to_ves linkspec2_pos = linkspec2->getPosCartesian_rel();
                linkspec1->removeLinkSpecPair();
                linkspec2->removeLinkSpecPair();
                vesicle1->remLinkSpec(linkspec1->getUniqueID(), linkspec1);
                vesicle2->remLinkSpec(linkspec2->getUniqueID(), linkspec2);
                vesicle1->addOneSurfSpec(vub_spec1_gidx, linkspec1_pos);
                vesicle2->addOneSurfSpec(vub_spec2_gidx, linkspec2_pos);
                vesicle1->updImmobility(selectedUnbind->immobility());
                vesicle2->updImmobility(selectedUnbind->immobility());

                linkspecpair_del.emplace_back(lspit);
            }
        }
    }
    // Delete the linkspecpairs that need deleting.
    for (auto const& lspit: linkspecpair_del) {
        delete *lspit;
        pLinkSpecPairs.erase(lspit);
    }
    // Finally, do the vesicle diffusion
    for (auto const& comp: pComps) {
        comp->runVesicle(dt);
    }
}

////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::advance(double adv) {
    if (adv < 0.0) {
        std::ostringstream os;
        os << "Time to advance cannot be negative";
        ArgErrLog(os.str());
    }

    double endtime = statedef().time() + adv;
    run(endtime);
}

////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::step() {
    NotImplErrLog("This function is not implemented.");
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::getTime() const {
    return statedef().time();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::getA0() const {
    return MPI_ConditionalReduce<double>(0, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::getNSteps() const {
    return statedef().nsteps();
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setTime(double time) {
    statedef().setTime(time);
}

////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setNSteps(uint nsteps) {
    statedef().setNSteps(nsteps);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::getVesicleDT() const {
    return MPI_ConditionalBcast<double>(
        pVesicledt, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setVesicleDT(double dt) {
    if (dt <= 0.0) {
        std::ostringstream os;
        os << "Vesicle dt must be greater than zero.";
        ArgErrLog(os.str());
    }
    pVesicledt = dt;
    _recalcQtables();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getCompVol(solver::comp_global_id cidx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    CompVesRaft* comp = getComp_(cidx);
    AssertLog(comp != nullptr);
    return comp->vol();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getCompSpecCount(solver::comp_global_id cidx,
                                            solver::spec_global_id sidx) const {
    return MPI_ConditionalReduce<uint>(0, MPI_UNSIGNED, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setCompSpecCount(solver::comp_global_id cidx,
                                          solver::spec_global_id sidx,
                                          double n) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getCompSpecAmount(solver::comp_global_id cidx,
                                             solver::spec_global_id sidx) const {
    double count = _getCompSpecCount(cidx, sidx);
    return (count / math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setCompSpecAmount(solver::comp_global_id cidx,
                                           solver::spec_global_id sidx,
                                           double a) {
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setCompSpecCount(cidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getCompSpecConc(solver::comp_global_id cidx,
                                           solver::spec_global_id sidx) const {
    return _getCompSpecCount(cidx, sidx) / (1.0e3 * getComp_(cidx)->vol() * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setCompSpecConc(solver::comp_global_id cidx,
                                         solver::spec_global_id sidx,
                                         double c) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    AssertLog(c >= 0.0);
    CompVesRaft* comp = getComp_(cidx);
    AssertLog(comp != nullptr);

    double count = c * (1.0e3 * comp->vol() * math::AVOGADRO);
    // the following method does argument checking on sidx
    _setCompSpecCount(cidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getCompSpecClamped(solver::comp_global_id cidx,
                                            solver::spec_global_id sidx) const {
    return MPI_ConditionalReduce<bool>(true, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setCompSpecClamped(solver::comp_global_id /*cidx*/,
                                            solver::spec_global_id /*sidx*/,
                                            bool /*b*/) { /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getCompReacK(solver::comp_global_id cidx,
                                        solver::reac_global_id ridx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    CompVesRaft* comp = getComp_(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = _reacG2L_or_throw(comp, ridx);

    // We're just returning the default value for this comp, individual
    // tets may have different Kcsts set individually
    return comp->def()->kcst(lridx);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setCompReacK(solver::comp_global_id cidx,
                                      solver::reac_global_id ridx,
                                      double kf) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    AssertLog(kf >= 0.0);
    CompVesRaft* comp = getComp_(cidx);
    AssertLog(comp != nullptr);
    solver::reac_local_id lridx = _reacG2L_or_throw(comp, ridx);

    // First set the default value for the comp
    comp->def()->setKcst(lridx, kf);
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getCompReacActive(solver::comp_global_id cidx,
                                           solver::reac_global_id ridx) const {
    return MPI_ConditionalReduce<bool>(true, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setCompReacActive(solver::comp_global_id /*cidx*/,
                                           solver::reac_global_id /*ridx*/,
                                           bool /*a*/) { /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getCompDiffD(solver::comp_global_id cidx,
                                        solver::diff_global_id didx) const {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    CompVesRaft* comp = getComp_(cidx);
    AssertLog(comp != nullptr);
    solver::diff_local_id ldidx = _diffG2L_or_throw(comp, didx);

    // We're just returning the default value for this comp, individual
    // tets may have different Dcsts set individually
    return comp->def()->dcst(ldidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setCompDiffD(solver::comp_global_id cidx,
                                      solver::diff_global_id didx,
                                      double dk) {
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    AssertLog(dk >= 0.0);
    CompVesRaft* comp = getComp_(cidx);
    AssertLog(comp != nullptr);
    solver::diff_local_id ldidx = _diffG2L_or_throw(comp, didx);

    // First set the default value for the comp
    comp->def()->setDcst(ldidx, dk);
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getCompDiffActive(solver::comp_global_id cidx,
                                           solver::diff_global_id didx) const {
    return MPI_ConditionalReduce<bool>(true, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setCompDiffActive(solver::comp_global_id cidx,
                                           solver::diff_global_id didx,
                                           bool act) {}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setCompVesicleCount(solver::comp_global_id cidx,
                                             solver::vesicle_global_id vidx,
                                             uint n) {
    CompVesRaft* comp = getComp_(cidx);

    comp->setVesicleCount(vidx, n);

    pRequireVesicleCommunication = true;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getCompVesicleCount(solver::comp_global_id cidx,
                                             solver::vesicle_global_id vidx) const {
    CompVesRaft* comp = getComp_(cidx);

    uint vesicle_count = comp->getVesicleCount(vidx);
    return MPI_ConditionalBcast<uint>(
        vesicle_count, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTetVesicleDcst(tetrahedron_global_id tidx,
                                           solver::vesicle_global_id vidx,
                                           double dcst) {
    auto* comp = _getTet(tidx)->getCompVesRaft();

    comp->setVesicleTetDcst(vidx, tidx, dcst);
}

////////////////////////////////////////////////////////////////////////////////

solver::vesicle_individual_id TetVesicleVesRaft::_addCompVesicle(solver::comp_global_id cidx,
                                                                 solver::vesicle_global_id vidx) {
    CompVesRaft* comp = getComp_(cidx);

    // This is done in this layer due to behavior of the comp::addVesicle
    // function The function expects a position to test
    solver::Vesicledef* vesdef = statedef().vesicledef(vidx);
    math::position_abs pos;
    solver::vesicle_individual_id added_vesicle;
    uint attempts = 0;
    // added_vesicle will return -1 if unsuccesful
    while (added_vesicle.unknown()) {
        attempts++;
        if (attempts > 10000) {
            CLOG(WARNING, "general_log") << "Unable to add a vesicle: too many iterations.\n";
            return {};
        }

        tetrahedron_global_id tet_gidx = comp->getRandPosByTetStaticVols(&pos);
        added_vesicle = comp->addVesicle(vesdef, pos, tet_gidx);
    }

    pRequireVesicleCommunication = true;

    return MPI_ConditionalBcast(
        added_vesicle, MPI_STEPS_INDEX, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getCompVesicleSurfaceSpecCount(solver::comp_global_id cidx,
                                                        solver::vesicle_global_id vidx,
                                                        solver::spec_global_id sidx) const {
    CompVesRaft* comp = getComp_(cidx);
    uint return_count = comp->getVesicleSurfaceSpecCount(vidx, sidx);
    return MPI_ConditionalBcast<uint>(
        return_count, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getCompVesicleInnerSpecCount(solver::comp_global_id cidx,
                                                      solver::vesicle_global_id vidx,
                                                      solver::spec_global_id sidx) const {
    CompVesRaft* comp = getComp_(cidx);
    uint return_count = comp->getVesicleInnerSpecCount(vidx, sidx);
    return MPI_ConditionalBcast<uint>(
        return_count, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::map<solver::vesicle_individual_id, uint> TetVesicleVesRaft::_getCompVesicleSurfaceSpecCountMap(
    solver::comp_global_id cidx,
    solver::vesicle_global_id vidx,
    solver::spec_global_id sidx) const {
    CompVesRaft* comp = getComp_(cidx);
    auto return_map = comp->getVesicleSurfaceSpecCountMap(vidx, sidx);
    MPI_ConditionalBcast(
        return_map, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_map;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getCompVesicleSurfaceLinkSpecCount(
    solver::comp_global_id cidx,
    solver::vesicle_global_id vidx,
    solver::linkspec_global_id lsidx) const {
    CompVesRaft* comp = getComp_(cidx);
    uint return_count = comp->getVesicleLinkSpecCount(vidx, lsidx);
    return MPI_ConditionalBcast<uint>(
        return_count, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::map<solver::vesicle_individual_id, uint>
TetVesicleVesRaft::_getCompVesicleSurfaceLinkSpecCountMap(solver::comp_global_id cidx,
                                                          solver::vesicle_global_id vidx,
                                                          solver::linkspec_global_id lsidx) const {
    CompVesRaft* comp = getComp_(cidx);
    auto return_map = comp->getVesicleLinkSpecCountMap(vidx, lsidx);
    MPI_ConditionalBcast(
        return_map, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_map;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::vesicle_individual_id> TetVesicleVesRaft::_getCompVesicleIndices(
    solver::comp_global_id cidx,
    solver::vesicle_global_id vidx) const {
    CompVesRaft* comp = getComp_(cidx);
    std::vector<solver::vesicle_individual_id> return_vec(comp->getVesicleIndices(vidx));
    MPI_ConditionalBcast(return_vec, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

solver::comp_global_id TetVesicleVesRaft::_getSingleVesicleCompartment(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index) const {
    for (auto* comp: pComps) {
        auto& vesicles = comp->getAllVesicles();
        auto vesit = vesicles.find(vidx);
        if (vesit != vesicles.end()) {
            for (auto& ves: vesit->second) {
                if (ves->getUniqueIndex() == ves_unique_index) {
                    return MPI_ConditionalBcast(comp->def()->gidx(),
                                                MPI_STEPS_INDEX,
                                                vesraftRank_World,
                                                myRank_World,
                                                syncOutput,
                                                outputRank);
                }
            }
        }
    }
    return MPI_ConditionalBcast(solver::comp_global_id(),
                                MPI_STEPS_INDEX,
                                vesraftRank_World,
                                myRank_World,
                                syncOutput,
                                outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setSingleVesicleInnerSpecCount(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::spec_global_id sidx,
    uint c) {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    if (vesicle_it == pVesicles.end()) {
        CLOG(WARNING, "general_log") << "Vesicle unique id " << ves_unique_index << "unknown.\n";
        return;
    }
    if (vesicle_it->second->idx() != vidx) {
        CLOG(WARNING, "general_log")
            << "Incorrect vesicle type for vesicle unique ID " << ves_unique_index << ".\n ";
        return;
    }

    vesicle_it->second->setInnerSpecCount(sidx, c);

    pRequireVesicleCommunication = true;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getSingleVesicleSurfaceSpecCount(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::spec_global_id sidx) const {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    if (vesicle_it == pVesicles.end()) {
        CLOG(WARNING, "general_log") << "Vesicle unique id " << ves_unique_index << "unknown.\n";
        return MPI_ConditionalBcast<uint>(
            0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
    }

    if (vesicle_it->second->idx() != vidx) {
        CLOG(WARNING, "general_log")
            << "Incorrect vesicle type for vesicle unique ID " << ves_unique_index << ".\n ";
        return MPI_ConditionalBcast<uint>(
            0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
    }

    uint spec_count = vesicle_it->second->getSurfSpecCount(sidx);

    return MPI_ConditionalBcast<uint>(
        spec_count, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getSingleVesicleInnerSpecCount(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::spec_global_id sidx) const {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    if (vesicle_it == pVesicles.end()) {
        CLOG(WARNING, "general_log") << "Vesicle unique id " << ves_unique_index << "unknown.\n";
        return MPI_ConditionalBcast<uint>(
            0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
    }

    if (vesicle_it->second->idx() != vidx) {
        CLOG(WARNING, "general_log")
            << "Incorrect vesicle type for vesicle unique ID " << ves_unique_index << ".\n ";
        return MPI_ConditionalBcast<uint>(
            0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
    }

    uint spec_count = vesicle_it->second->getInnerSpecCount(sidx);

    return MPI_ConditionalBcast<uint>(
        spec_count, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setSingleVesicleSurfaceSpecCount(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::spec_global_id sidx,
    uint c) {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    if (vesicle_it == pVesicles.end()) {
        CLOG(WARNING, "general_log") << "Vesicle unique id " << ves_unique_index << "unknown.\n";
        return;
    }

    if (vesicle_it->second->idx() != vidx) {
        CLOG(WARNING, "general_log")
            << "Incorrect vesicle type for vesicle unique ID " << ves_unique_index << ".\n ";
        return;
    }

    vesicle_it->second->setSurfSpecCount(sidx, c);

    pRequireVesicleCommunication = true;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_deleteSingleVesicle(solver::vesicle_global_id vidx,
                                             solver::vesicle_individual_id ves_unique_index) {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    if (vesicle_it == pVesicles.end()) {
        CLOG(WARNING, "general_log") << "Vesicle unique id " << ves_unique_index << "unknown.\n";
        return;
    }
    if (vesicle_it->second->idx() != vidx) {
        CLOG(WARNING, "general_log")
            << "Incorrect vesicle type for vesicle unique ID " << ves_unique_index << ".\n ";
        return;
    }

    vesicle_it->second->comp()->deleteSingleVesicle(vesicle_it->second);

    pRequireVesicleCommunication = true;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getSingleVesicleSurfaceLinkSpecCount(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::linkspec_global_id lsidx) const {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    if (vesicle_it == pVesicles.end()) {
        CLOG(WARNING, "general_log") << "Vesicle unique id " << ves_unique_index << "unknown.\n";
        return MPI_ConditionalBcast<uint>(
            0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
    }

    if (vesicle_it->second->idx() != vidx) {
        CLOG(WARNING, "general_log")
            << "Incorrect vesicle type for vesicle unique ID " << ves_unique_index << ".\n ";
        return MPI_ConditionalBcast<uint>(
            0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
    }

    uint spec_count = vesicle_it->second->getLinkSpecCount(lsidx);

    return MPI_ConditionalBcast<uint>(
        spec_count, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::linkspec_individual_id>
TetVesicleVesRaft::_getSingleVesicleSurfaceLinkSpecIndices(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::linkspec_global_id lsidx) const {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    std::vector<solver::linkspec_individual_id> empty_vec;
    if (vesicle_it == pVesicles.end()) {
        CLOG(WARNING, "general_log") << "Vesicle unique id " << ves_unique_index << "unknown.\n";
        MPI_ConditionalBcast(empty_vec, vesraftRank_World, myRank_World, syncOutput, outputRank);
        return empty_vec;
    }

    if (vesicle_it->second->idx() != vidx) {
        CLOG(WARNING, "general_log")
            << "Incorrect vesicle type for vesicle unique ID " << ves_unique_index << ".\n ";
        MPI_ConditionalBcast(empty_vec, vesraftRank_World, myRank_World, syncOutput, outputRank);
        return empty_vec;
    }

    std::vector<solver::linkspec_individual_id> return_vec = vesicle_it->second->getLinkSpecIndices(
        lsidx);

    MPI_ConditionalBcast(return_vec, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::pointspec_individual_id> TetVesicleVesRaft::_getSingleVesicleSurfaceSpecIndices(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::spec_global_id sidx) const {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    std::vector<solver::pointspec_individual_id> empty_vec;
    if (vesicle_it == pVesicles.end()) {
        CLOG(WARNING, "general_log") << "Vesicle unique id " << ves_unique_index << "unknown.\n";
        MPI_ConditionalBcast(empty_vec, vesraftRank_World, myRank_World, syncOutput, outputRank);
        return empty_vec;
    }

    if (vesicle_it->second->idx() != vidx) {
        CLOG(WARNING, "general_log")
            << "Incorrect vesicle type for vesicle unique ID " << ves_unique_index << ".\n ";
        MPI_ConditionalBcast(empty_vec, vesraftRank_World, myRank_World, syncOutput, outputRank);
        return empty_vec;
    }

    std::vector<solver::pointspec_individual_id> return_vec =
        vesicle_it->second->getSurfaceSpecIndices(sidx);

    MPI_ConditionalBcast(return_vec, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setCompSingleVesiclePos(solver::comp_global_id cidx,
                                                 solver::vesicle_global_id vidx,
                                                 solver::vesicle_individual_id ves_unique_index,
                                                 const std::vector<double>& pos,
                                                 bool force) {
    CompVesRaft* comp = getComp_(cidx);
    comp->setVesiclePos(vidx, ves_unique_index, pos, force);

    pRequireVesicleCommunication = true;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleVesRaft::_getSingleVesiclePos(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index) const {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    ProgErrLogIf(vesicle_it == pVesicles.end(), "Vesicle unique id unknown.\n");
    ProgErrLogIf(vesicle_it->second->idx() != vidx, "Incorrect vesicle type.\n ");

    math::position_abs v_pos_p3d = vesicle_it->second->getPosition();
    std::vector<double> return_vec = {v_pos_p3d[0], v_pos_p3d[1], v_pos_p3d[2]};

    auto nentries = MPI_ConditionalBcast<std::size_t>(
        return_vec.size(), MPI_STD_SIZE_T, vesraftRank_World, myRank_World, syncOutput, outputRank);
    MPI_ConditionalBcast<double>(
        return_vec, nentries, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> TetVesicleVesRaft::_getSingleVesicleSurfaceLinkSpecPos(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::linkspec_global_id lsidx) {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    ProgErrLogIf(vesicle_it == pVesicles.end(), "Vesicle unique id unknown.\n");
    ProgErrLogIf(vesicle_it->second->idx() != vidx, "Incorrect vesicle type.\n ");

    std::vector<std::vector<double>> return_mv = vesicle_it->second->getLinkSpecPos(lsidx);

    MPI_ConditionalBcast<double>(
        return_mv, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_mv;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleVesRaft::_getSingleLinkSpecPos(
    solver::linkspec_individual_id ls_unique_id) const {
    auto ls_it = pLinkSpecs.find(ls_unique_id);
    ProgErrLogIf(ls_it == pLinkSpecs.end(), "Link species unique id unknown.\n");

    math::position_abs pos = ls_it->second->getPosCartesian_abs();
    std::vector<double> pos_vec(std::initializer_list<double>({pos[0], pos[1], pos[2]}));

    auto nentries = MPI_ConditionalBcast<std::size_t>(
        pos_vec.size(), MPI_STD_SIZE_T, vesraftRank_World, myRank_World, syncOutput, outputRank);
    MPI_ConditionalBcast(
        pos_vec, nentries, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return pos_vec;
}

////////////////////////////////////////////////////////////////////////////////

solver::linkspec_individual_id TetVesicleVesRaft::_getSingleLinkSpecLinkedTo(
    solver::linkspec_individual_id ls_unique_id) const {
    auto ls_it = pLinkSpecs.find(ls_unique_id);
    ProgErrLogIf(ls_it == pLinkSpecs.end(), "Link species unique id unknown.\n");

    auto ls_id = ls_it->second->getLinkedSpec()->getUniqueID();

    return MPI_ConditionalBcast(
        ls_id, MPI_STEPS_INDEX, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

solver::vesicle_individual_id TetVesicleVesRaft::_getSingleLinkSpecVes(
    solver::linkspec_individual_id ls_unique_id) const {
    auto ls_it = pLinkSpecs.find(ls_unique_id);

    if (ls_it == pLinkSpecs.end()) {
        return MPI_ConditionalBcast(solver::vesicle_individual_id(),
                                    MPI_STEPS_INDEX,
                                    vesraftRank_World,
                                    myRank_World,
                                    syncOutput,
                                    outputRank);
    } else {
        auto ves_id = ls_it->second->getVesicle()->getUniqueIndex();
        return MPI_ConditionalBcast(
            ves_id, MPI_STEPS_INDEX, vesraftRank_World, myRank_World, syncOutput, outputRank);
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> TetVesicleVesRaft::_getSingleVesicleSurfaceSpecPos(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::spec_global_id sidx) {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    ProgErrLogIf(vesicle_it == pVesicles.end(), "Vesicle unique id unknown.\n");
    ProgErrLogIf(vesicle_it->second->idx() != vidx, "Incorrect vesicle type.\n ");

    auto return_mv = vesicle_it->second->getSurfaceSpecPos(sidx);

    MPI_ConditionalBcast(
        return_mv, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_mv;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleVesRaft::_getSingleSpecPosSpherical(
    solver::spec_global_id sidx,
    solver::pointspec_individual_id ps_unique_id) const {
    std::vector<double> return_vec;

    for (auto const& ves: pVesicles) {
        for (auto const& ps: ves.second->getSurfSpecs()[sidx]) {
            if (ps->getUniqueIndex() == ps_unique_id) {
                auto pos = ps->getPosSpherical();
                return_vec.emplace_back(pos.getTheta());
                return_vec.emplace_back(pos.getPhi());
            }
        }
    }

    MPI_ConditionalBcast(
        return_vec, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);

    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> TetVesicleVesRaft::_getSingleVesicleSurfaceSpecPosSpherical(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::spec_global_id sidx) const {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    if (vesicle_it == pVesicles.end()) {
        CLOG(WARNING, "general_log") << "Vesicle unique id " << ves_unique_index << "unknown.\n";
        return {};
    }
    if (vesicle_it->second->idx() != vidx) {
        CLOG(WARNING, "general_log")
            << "Incorrect vesicle type for vesicle unique ID " << ves_unique_index << ".\n ";
        return {};
    }

    auto return_vec = vesicle_it->second->getSurfaceSpecPosSpherical(sidx);
    MPI_ConditionalBcast(
        return_vec, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);

    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setSingleVesicleSurfaceSpecPosSpherical(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index,
    solver::spec_global_id sidx,
    const std::vector<std::vector<double>>& pos_spherical) {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    if (vesicle_it == pVesicles.end()) {
        CLOG(WARNING, "general_log") << "Vesicle unique id " << ves_unique_index << "unknown.\n";
        return;
    }
    if (vesicle_it->second->idx() != vidx) {
        CLOG(WARNING, "general_log")
            << "Incorrect vesicle type for vesicle unique ID " << ves_unique_index << ".\n ";
        return;
    }

    vesicle_it->second->setSurfaceSpecPosSpherical(sidx, pos_spherical);

    pRequireVesicleCommunication = true;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getSingleVesicleImmobility(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index) const {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    ProgErrLogIf(vesicle_it == pVesicles.end(), "Vesicle unique id unknown.\n");
    ProgErrLogIf(vesicle_it->second->idx() != vidx, "Incorrect vesicle type.\n ");

    return MPI_ConditionalBcast<uint>(vesicle_it->second->getImmobility(),
                                      MPI_UNSIGNED,
                                      vesraftRank_World,
                                      myRank_World,
                                      syncOutput,
                                      outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::tetrahedron_global_id> TetVesicleVesRaft::_getSingleVesicleOverlapTets(
    solver::vesicle_global_id vidx,
    solver::vesicle_individual_id ves_unique_index) const {
    auto const& vesicle_it = pVesicles.find(ves_unique_index);
    ProgErrLogIf(vesicle_it == pVesicles.end(), "Vesicle unique id unknown.\n");
    ProgErrLogIf(vesicle_it->second->idx() != vidx, "Incorrect vesicle type.\n ");

    std::vector<steps::tetrahedron_global_id> return_vec = vesicle_it->second->getOverlapVec_gidx();

    MPI_ConditionalBcast(return_vec, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getPatchArea(solver::patch_global_id pidx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    PatchVesRaft* patch = getPatch_(pidx);
    AssertLog(patch != nullptr);
    return patch->area();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getPatchSpecCount(solver::patch_global_id pidx,
                                             solver::spec_global_id sidx) const {
    return MPI_ConditionalReduce<uint>(0, MPI_UNSIGNED, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setPatchSpecCount(solver::patch_global_id pidx,
                                           solver::spec_global_id sidx,
                                           double n) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getPatchSpecAmount(solver::patch_global_id pidx,
                                              solver::spec_global_id sidx) const {
    double count = _getPatchSpecCount(pidx, sidx);
    return (count / math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setPatchSpecAmount(solver::patch_global_id pidx,
                                            solver::spec_global_id sidx,
                                            double a) {
    AssertLog(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setPatchSpecCount(pidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getPatchSpecClamped(solver::patch_global_id pidx,
                                             solver::spec_global_id sidx) const {
    return MPI_ConditionalReduce<bool>(true, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setPatchSpecClamped(solver::patch_global_id /*pidx*/,
                                             solver::spec_global_id /*sidx*/,
                                             bool /*buf*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getPatchSReacK(solver::patch_global_id pidx,
                                          solver::sreac_global_id ridx) const {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    PatchVesRaft* patch = getPatch_(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lsridx = _sreacG2L_or_throw(patch, ridx);

    // We're just returning the default value for this patch, individual
    // triangles may have different Kcsts set
    return patch->def()->kcst(lsridx);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setPatchSReacK(solver::patch_global_id pidx,
                                        solver::sreac_global_id ridx,
                                        double kf) {
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    AssertLog(kf >= 0.0);
    PatchVesRaft* patch = getPatch_(pidx);
    AssertLog(patch != nullptr);
    solver::sreac_local_id lsridx = _sreacG2L_or_throw(patch, ridx);

    // First set the default values for this patch
    patch->def()->setKcst(lsridx, kf);
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getPatchSReacActive(solver::patch_global_id pidx,
                                             solver::sreac_global_id ridx) const {
    return MPI_ConditionalReduce<bool>(true, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setDiffBoundarySpecDiffusionActive(
    solver::diffboundary_global_id /*dbidx*/,
    solver::spec_global_id /*sidx*/,
    bool /*act*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getDiffBoundarySpecDiffusionActive(solver::diffboundary_global_id dbidx,
                                                            solver::spec_global_id sidx) const {
    return MPI_ConditionalReduce<bool>(true, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setDiffBoundarySpecDcst(solver::diffboundary_global_id /*dbidx*/,
                                                 solver::spec_global_id /*sidx*/,
                                                 double /*dcst*/,
                                                 solver::comp_global_id /*direction_comp*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setPatchSReacActive(solver::patch_global_id /*pidx*/,
                                             solver::sreac_global_id /*ridx*/,
                                             bool /*a*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getPatchVDepSReacActive(solver::patch_global_id pidx,
                                                 solver::vdepsreac_global_id vsridx) const {
    return MPI_ConditionalReduce<bool>(true, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setPatchVDepSReacActive(solver::patch_global_id /*pidx*/,
                                                 solver::vdepsreac_global_id /*vsridx*/,
                                                 bool /*a*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setSDiffBoundarySpecDiffusionActive(
    solver::sdiffboundary_global_id /*sdbidx*/,
    solver::spec_global_id /*sidx*/,
    bool /*act*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getSDiffBoundarySpecDiffusionActive(solver::sdiffboundary_global_id sdbidx,
                                                             solver::spec_global_id sidx) const {
    return MPI_ConditionalReduce<bool>(true, MPI_C_BOOL, MPI_LAND, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setSDiffBoundarySpecDcst(solver::sdiffboundary_global_id /*sdbidx*/,
                                                  solver::spec_global_id /*sidx*/,
                                                  double /*dcst*/,
                                                  solver::patch_global_id /*direction_patch*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////
///////////////////// RAFTS, ENDOCYTOSIS, Patch-specific
/////////////////////////

void TetVesicleVesRaft::_setPatchRaftCount(solver::patch_global_id pidx,
                                           solver::raft_global_id ridx,
                                           uint n) {
    PatchVesRaft* patch = getPatch_(pidx);

    patch->setRaftCount(ridx, n);

    pRequireVesicleCommunication = true;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getPatchRaftCount(solver::patch_global_id pidx,
                                           solver::raft_global_id ridx) const {
    PatchVesRaft* patch = getPatch_(pidx);
    return MPI_ConditionalBcast<uint>(patch->getRaftCount(ridx),
                                      MPI_UNSIGNED,
                                      vesraftRank_World,
                                      myRank_World,
                                      syncOutput,
                                      outputRank);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getSingleRaftImmobility(
    solver::raft_global_id ridx,
    solver::raft_individual_id raft_unique_index) const {
    auto const& raft_it = pRafts.find(raft_unique_index);
    ProgErrLogIf(raft_it == pRafts.end(), "Raft unique id unknown.\n");
    ProgErrLogIf(raft_it->second->idx() != ridx, "Incorrect Raft type.\n ");

    uint immobility = raft_it->second->getImmobility();

    return MPI_ConditionalBcast<uint>(
        immobility, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::raft_individual_id> TetVesicleVesRaft::_getPatchRaftIndices(
    solver::patch_global_id pidx,
    solver::raft_global_id ridx) const {
    PatchVesRaft* patch = getPatch_(pidx);
    std::vector<solver::raft_individual_id> return_vec = patch->getRaftIndices(ridx);
    MPI_ConditionalBcast(return_vec, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

solver::patch_global_id TetVesicleVesRaft::_getSingleRaftPatch(
    solver::raft_global_id ridx,
    solver::raft_individual_id raft_unique_index) const {
    for (auto* patch: pPatches) {
        auto& rafts = patch->getAllRafts();
        auto raftit = rafts.find(ridx);
        if (raftit != rafts.end()) {
            for (auto& raft: raftit->second) {
                if (raft->getUniqueIndex() == raft_unique_index) {
                    return MPI_ConditionalBcast<index_t>(patch->def()->gidx(),
                                                         MPI_STEPS_INDEX,
                                                         vesraftRank_World,
                                                         myRank_World,
                                                         syncOutput,
                                                         outputRank);
                }
            }
        }
    }
    return MPI_ConditionalBcast(solver::patch_global_id(),
                                MPI_STEPS_INDEX,
                                vesraftRank_World,
                                myRank_World,
                                syncOutput,
                                outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleVesRaft::_getSingleRaftPos(
    solver::raft_global_id ridx,
    solver::raft_individual_id raft_unique_index) const {
    auto const& raft_it = pRafts.find(raft_unique_index);
    ProgErrLogIf(raft_it == pRafts.end(), "Raft unique id unknown.\n");
    ProgErrLogIf(raft_it->second->idx() != ridx, "Incorrect Raft type.\n ");

    math::position_abs r_pos_p3d = raft_it->second->getPosition();
    std::vector<double> return_vec = {r_pos_p3d[0], r_pos_p3d[1], r_pos_p3d[2]};

    auto nentries = MPI_ConditionalBcast<std::size_t>(
        return_vec.size(), MPI_STD_SIZE_T, vesraftRank_World, myRank_World, syncOutput, outputRank);
    MPI_ConditionalBcast<double>(
        return_vec, nentries, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_vec;
}

////////////////////////////////////////////////////////////////////////////////

std::map<solver::raft_individual_id, uint> TetVesicleVesRaft::_getPatchRaftSpecCountMap(
    solver::patch_global_id pidx,
    solver::raft_global_id ridx,
    solver::spec_global_id sidx) const {
    PatchVesRaft* patch = getPatch_(pidx);
    std::map<solver::raft_individual_id, uint> return_map = patch->getRaftSpecCountMap(ridx, sidx);
    MPI_ConditionalBcast(
        return_map, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
    return return_map;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getPatchRaftSpecCount(solver::patch_global_id pidx,
                                               solver::raft_global_id ridx,
                                               solver::spec_global_id sidx) const {
    PatchVesRaft* patch = getPatch_(pidx);
    uint return_count = patch->getRaftSpecCount(ridx, sidx);
    return MPI_ConditionalBcast<uint>(
        return_count, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getSingleRaftSpecCount(solver::raft_global_id ridx,
                                                solver::raft_individual_id raft_unique_index,
                                                solver::spec_global_id sidx) const {
    auto const& raft_it = pRafts.find(raft_unique_index);
    if (raft_it == pRafts.end()) {
        CLOG(WARNING, "general_log") << "Raft unique id " << raft_unique_index << "unknown.\n";
        return MPI_ConditionalBcast<uint>(
            0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
    }
    if (raft_it->second->idx() != ridx) {
        CLOG(WARNING, "general_log")
            << "Incorrect Raft type for raft unique ID " << raft_unique_index << ".\n ";
        return MPI_ConditionalBcast<uint>(
            0, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
    }

    uint spec_count = raft_it->second->getSpecCountByGidx(sidx);

    return MPI_ConditionalBcast<uint>(
        spec_count, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setSingleRaftSpecCount(solver::raft_global_id ridx,
                                                solver::raft_individual_id raft_unique_index,
                                                solver::spec_global_id sidx,
                                                uint c) {
    setSingleRaftSpecCount_(ridx, raft_unique_index, sidx, c);

    pRequireVesicleCommunication = true;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setSingleRaftSpecCount_(solver::raft_global_id ridx,
                                                solver::raft_individual_id raft_unique_index,
                                                solver::spec_global_id sidx,
                                                uint c) {
    auto const& raft_it = pRafts.find(raft_unique_index);
    if (raft_it == pRafts.end()) {
        CLOG(WARNING, "general_log") << "Raft unique id " << raft_unique_index << "unknown.\n";
        return;
    }
    if (raft_it->second->idx() != ridx) {
        CLOG(WARNING, "general_log")
            << "Incorrect Raft type for raft unique ID " << raft_unique_index << ".\n ";
        return;
    }

    raft_it->second->setSpecCountByGidx(sidx, c);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setPatchEndocyticZoneEndocytosisActive(
    solver::patch_global_id pidx,
    std::string const& zone,
    solver::endocytosis_global_id endogidx,
    bool active) {
    PatchVesRaft* patch = getPatch_(pidx);
    solver::endocytosis_local_id endolidx = _endoG2L_or_throw(patch, endogidx);

    patch->getEndocytosis(endolidx, zone)->setActive(active);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setPatchEndocyticZoneEndocytosisK(solver::patch_global_id pidx,
                                                           std::string const& zone,
                                                           solver::endocytosis_global_id endogidx,
                                                           double k) {
    PatchVesRaft* patch = getPatch_(pidx);
    solver::endocytosis_local_id endolidx = _endoG2L_or_throw(patch, endogidx);

    patch->getEndocytosis(endolidx, zone)->setKcst(k);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getPatchEndocyticZoneEndocytosisExtent(
    solver::patch_global_id pidx,
    std::string const& zone,
    solver::endocytosis_global_id endogidx) const {
    PatchVesRaft* patch = getPatch_(pidx);
    solver::endocytosis_local_id endolidx = _endoG2L_or_throw(patch, endogidx);

    return MPI_ConditionalBcast<uint>(patch->getEndocytosis(endolidx, zone)->getExtent(),
                                      MPI_UNSIGNED,
                                      vesraftRank_World,
                                      myRank_World,
                                      syncOutput,
                                      outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::EndocytosisEvent> TetVesicleVesRaft::_getPatchEndocyticZoneEndocytosisEvents(
    solver::patch_global_id pidx,
    std::string const& zone,
    solver::endocytosis_global_id endogidx) const {
    PatchVesRaft* patch = getPatch_(pidx);
    solver::endocytosis_local_id endolidx = _endoG2L_or_throw(patch, endogidx);

    auto events = patch->getEndocytosis(endolidx, zone)->getEvents();
    MPI_ConditionalBcast(events,
                         dataTypeUtil.MPI_EndocytosisEventSync,
                         vesraftRank_World,
                         myRank_World,
                         syncOutput,
                         outputRank);
    return events;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getCompReacH(solver::comp_global_id cidx,
                                        solver::reac_global_id ridx) const {
    return MPI_ConditionalReduce<double>(0.0, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getCompReacC(solver::comp_global_id cidx,
                                        solver::reac_global_id ridx) const {
    auto global_c = MPI_ConditionalReduce<double>(0.0, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
    auto global_v = MPI_ConditionalReduce<double>(0.0, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
    return global_c / global_v;
}

////////////////////////////////////////////////////////////////////////////////

long double TetVesicleVesRaft::_getCompReacA(solver::comp_global_id cidx,
                                             solver::reac_global_id ridx) const {
    return MPI_ConditionalReduce<long double>(
        0.0L, MPI_LONG_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetVesicleVesRaft::_getCompReacExtent(solver::comp_global_id cidx,
                                                         solver::reac_global_id ridx) const {
    return MPI_ConditionalReduce<unsigned long long>(
        0.0L, MPI_UNSIGNED_LONG_LONG, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_resetCompReacExtent(solver::comp_global_id /*cidx*/,
                                             solver::reac_global_id /*ridx*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getPatchSReacH(solver::patch_global_id pidx,
                                          solver::sreac_global_id ridx) const {
    return MPI_ConditionalReduce<double>(0.0, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getPatchSReacC(solver::patch_global_id pidx,
                                          solver::sreac_global_id ridx) const {
    auto global_c = MPI_ConditionalReduce<double>(0.0, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
    auto global_a = MPI_ConditionalReduce<double>(0.0, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
    return global_c / global_a;
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getPatchSReacA(solver::patch_global_id pidx,
                                          solver::sreac_global_id ridx) const {
    return MPI_ConditionalReduce<double>(0.0, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetVesicleVesRaft::_getPatchSReacExtent(solver::patch_global_id pidx,
                                                           solver::sreac_global_id ridx) const {
    return MPI_ConditionalReduce<unsigned long long>(
        0L, MPI_UNSIGNED_LONG_LONG, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_resetPatchSReacExtent(solver::patch_global_id /*pidx*/,
                                               solver::sreac_global_id /*ridx*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTetVol(tetrahedron_global_id tidx) const {
    auto tet = pTets.at(tidx);
    if (tet == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
        ArgErrLog(os.str());
    }
    return tet->staticVol();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::getTetReducedVol(tetrahedron_global_id tidx) const {
    auto tet = pTets.at(tidx);
    if (tet == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
        ArgErrLog(os.str());
    }
    double vol = tet->vol();
    return MPI_ConditionalBcast<double>(
        vol, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTetVol(tetrahedron_global_id /*tidx*/, double /*vol*/) {
    std::ostringstream os;
    os << "Can not change tetrahedron volume in a mesh based solver.\n";
    NotImplErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getTetSpecDefined(tetrahedron_global_id tidx,
                                           solver::spec_global_id sidx) const {
    AssertLog(tidx < pTets.size());
    AssertLog(sidx < statedef().countSpecs());

    TetVesRaft* tet = _getTet(tidx);
    solver::spec_local_id lsidx = tet->compdef()->specG2L(sidx);
    return lsidx.valid();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTetSpecCount(tetrahedron_global_id tidx,
                                           solver::spec_global_id sidx) const {
    int host = _getTetHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTetSpecCount(tetrahedron_global_id /*tidx*/,
                                         solver::spec_global_id /*sidx*/,
                                         double /*n*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTetSpecAmount(tetrahedron_global_id tidx,
                                            solver::spec_global_id sidx) const {
    double count = _getTetSpecCount(tidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTetSpecAmount(tetrahedron_global_id /*tidx*/,
                                          solver::spec_global_id /*sidx*/,
                                          double /*m*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTetSpecConc(tetrahedron_global_id tidx,
                                          solver::spec_global_id sidx) const {
    // following method does all necessary argument checking
    double count = _getTetSpecCount(tidx, sidx);
    TetVesRaft* tet = _getTet(tidx);
    double vol = tet->staticVol();
    return (count / (1.0e3 * vol * math::AVOGADRO));
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTetSpecConc(tetrahedron_global_id /*tidx*/,
                                        solver::spec_global_id /*sidx*/,
                                        double /*c*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getTetSpecClamped(tetrahedron_global_id tidx,
                                           solver::spec_global_id sidx) const {
    int host = _getTetHost(tidx);
    return MPI_ConditionalBcast<bool>(
        false, MPI_C_BOOL, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTetSpecClamped(tetrahedron_global_id /*tidx*/,
                                           solver::spec_global_id /*sidx*/,
                                           bool /*buf*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTetReacK(tetrahedron_global_id tidx,
                                       solver::reac_global_id ridx) const {
    int host = _getTetHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTetReacK(tetrahedron_global_id /*tidx*/,
                                     solver::reac_global_id /*ridx*/,
                                     double /*kf*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getTetReacActive(tetrahedron_global_id tidx,
                                          solver::reac_global_id ridx) const {
    int host = _getTetHost(tidx);
    return MPI_ConditionalBcast<bool>(
        false, MPI_C_BOOL, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTetReacActive(tetrahedron_global_id /*tidx*/,
                                          solver::reac_global_id /*ridx*/,
                                          bool /*act*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTetDiffD(tetrahedron_global_id tidx,
                                       solver::diff_global_id didx,
                                       tetrahedron_global_id direction_tet) const {
    int host = _getTetHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTetDiffD(tetrahedron_global_id /*tidx*/,
                                     solver::diff_global_id /*didx*/,
                                     double /*dk*/,
                                     tetrahedron_global_id /*direction_tet*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getTetDiffActive(tetrahedron_global_id tidx,
                                          solver::diff_global_id didx) const {
    int host = _getTetHost(tidx);
    return MPI_ConditionalBcast<bool>(
        false, MPI_C_BOOL, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTetDiffActive(tetrahedron_global_id /*tidx*/,
                                          solver::diff_global_id /*didx*/,
                                          bool /*act*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTetReacH(tetrahedron_global_id tidx,
                                       solver::reac_global_id ridx) const {
    int host = _getTetHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTetReacC(tetrahedron_global_id tidx,
                                       solver::reac_global_id ridx) const {
    int host = _getTetHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTetReacA(tetrahedron_global_id tidx,
                                       solver::reac_global_id ridx) const {
    int host = _getTetHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTetDiffA(tetrahedron_global_id tidx,
                                       solver::diff_global_id didx) const {
    int host = _getTetHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriArea(triangle_global_id tidx) const {
    return _getTri(tidx)->area();
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTriArea(triangle_global_id /*tidx*/, double /*area*/) {
    NotImplErrLog("This function is not implemented.");
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getTriSpecDefined(triangle_global_id tidx,
                                           solver::spec_global_id sidx) const {
    AssertLog(sidx < statedef().countSpecs());

    TriVesRaft* tri = _getTri(tidx);
    solver::spec_local_id lsidx = tri->patchdef()->specG2L(sidx);
    return lsidx.valid();
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriSpecCount(triangle_global_id tidx,
                                           solver::spec_global_id sidx) const {
    int host = _getTriHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::getTriSpecCount_(triangle_global_id tidx,
                                           solver::spec_global_id sidx) const {
    AssertLog(sidx < statedef().countSpecs());

    TriVesRaft* tri = _getTri(tidx);
    solver::spec_local_id lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx.unknown()) {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        ArgErrLog(os.str());
    }
    return tri->pools()[lsidx];
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTriSpecCount(triangle_global_id /*tidx*/,
                                         solver::spec_global_id /*sidx*/,
                                         double /*n*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setTriSpecCount_(triangle_global_id tidx,
                                         solver::spec_global_id sidx,
                                         double n) {
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(n >= 0.0);

    if (n > UINT_MAX) {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << UINT_MAX << ").\n";
        ArgErrLog(os.str());
    }

    TriVesRaft* tri = _getTri(tidx);

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
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriSpecAmount(triangle_global_id tidx,
                                            solver::spec_global_id sidx) const {
    double count = _getTriSpecCount(tidx, sidx);
    return count / math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTriSpecAmount(triangle_global_id /*tidx*/,
                                          solver::spec_global_id /*sidx*/,
                                          double /*m*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getTriSpecClamped(triangle_global_id tidx,
                                           solver::spec_global_id sidx) const {
    int host = _getTriHost(tidx);
    return MPI_ConditionalBcast<bool>(
        false, MPI_C_BOOL, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTriSpecClamped(triangle_global_id /*tidx*/,
                                           solver::spec_global_id /*sidx*/,
                                           bool /*buf*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriSReacK(triangle_global_id tidx,
                                        solver::sreac_global_id ridx) const {
    int host = _getTriHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTriSReacK(triangle_global_id /*tidx*/,
                                      solver::sreac_global_id /*ridx*/,
                                      double /*kf*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getTriSReacActive(triangle_global_id tidx,
                                           solver::sreac_global_id ridx) const {
    int host = _getTriHost(tidx);
    return MPI_ConditionalBcast<bool>(
        false, MPI_C_BOOL, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTriSReacActive(triangle_global_id /*tidx*/,
                                           solver::sreac_global_id /*ridx*/,
                                           bool /*act*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriSDiffD(triangle_global_id tidx,
                                        solver::surfdiff_global_id didx,
                                        triangle_global_id direction_tri) const {
    int host = _getTriHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTriSDiffD(triangle_global_id /*unused*/,
                                      solver::surfdiff_global_id /*unused*/,
                                      double /*unused*/,
                                      triangle_global_id /*unused*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getTriVDepSReacActive(triangle_global_id tidx,
                                               solver::vdepsreac_global_id vsridx) const {
    int host = _getTriHost(tidx);
    return MPI_ConditionalBcast<bool>(
        false, MPI_C_BOOL, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTriVDepSReacActive(triangle_global_id /*tidx*/,
                                               solver::vdepsreac_global_id /*vsridx*/,
                                               bool /*act*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriSReacH(triangle_global_id tidx,
                                        solver::sreac_global_id ridx) const {
    int host = _getTriHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriSReacC(triangle_global_id tidx,
                                        solver::sreac_global_id ridx) const {
    int host = _getTriHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriSReacA(triangle_global_id tidx,
                                        solver::sreac_global_id ridx) const {
    int host = _getTriHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getTriRaftCount(triangle_global_id tidx,
                                         solver::raft_global_id ridx) const {
    AssertLog(ridx < statedef().countRafts());

    TriVesRaft* tri = _getTri(tidx);
    PatchVesRaft* patch = tri->patchVesRaft();

    return MPI_ConditionalBcast<uint>(patch->getRaftCount(ridx, tri),
                                      MPI_UNSIGNED,
                                      vesraftRank_World,
                                      myRank_World,
                                      syncOutput,
                                      outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTriRaftCount(triangle_global_id tidx,
                                         solver::raft_global_id ridx,
                                         uint n) {
    AssertLog(ridx < statedef().countRafts());

    TriVesRaft* tri = _getTri(tidx);
    PatchVesRaft* patch = tri->patchVesRaft();

    patch->setRaftCount(ridx, tri, n);

    pRequireVesicleCommunication = true;
}

////////////////////////////////////////////////////////////////////////////////

solver::raft_individual_id TetVesicleVesRaft::_addTriRaft(triangle_global_id tidx,
                                                          solver::raft_global_id ridx) {
    AssertLog(ridx < statedef().countRafts());

    TriVesRaft* tri = _getTri(tidx);
    PatchVesRaft* patch = tri->patchVesRaft();

    auto added_raft = patch->addRaft(ridx, tri);

    pRequireVesicleCommunication = true;

    return MPI_ConditionalBcast(
        added_raft, MPI_STEPS_INDEX, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTetV(tetrahedron_global_id /*tidx*/) const {
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTetV(tetrahedron_global_id /*tidx*/, double /*v*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getTetVClamped(tetrahedron_global_id tidx) const {
    return MPI_ConditionalBcast<bool>(
        false, MPI_C_BOOL, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTetVClamped(tetrahedron_global_id /*tidx*/, bool /*cl*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriV(triangle_global_id tidx) const {
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTriV(triangle_global_id /*tidx*/, double /*v*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getTriVClamped(triangle_global_id tidx) const {
    return MPI_ConditionalBcast<bool>(
        false, MPI_C_BOOL, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}
////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTriVClamped(triangle_global_id /*tidx*/, bool /*cl*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriOhmicErev(triangle_global_id tidx,
                                           solver::ohmiccurr_global_id ocgidx) const {
    int host = _getTriHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTriOhmicErev(triangle_global_id tidx,
                                         solver::ohmiccurr_global_id ocgidx,
                                         double erev) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriOhmicI(triangle_global_id tidx) const {
    int host = _getTriHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriOhmicI(triangle_global_id tidx,
                                        solver::ohmiccurr_global_id ocidx) const {
    int host = _getTriHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriGHKI(triangle_global_id tidx) const {
    int host = _getTriHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriGHKI(triangle_global_id tidx,
                                      solver::ghkcurr_global_id ghkidx) const {
    int host = _getTriHost(tidx);
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, host, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriI(triangle_global_id tidx) const {
    return _getTriGHKI(tidx) + _getTriOhmicI(tidx);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getVertIClamp(vertex_id_t vidx) const {
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}
////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setVertIClamp(vertex_id_t /*vidx*/, double /*i*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getTriIClamp(triangle_global_id tidx) const {
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTriIClamp(triangle_global_id /*tidx*/, double /*i*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setTriCapac(triangle_global_id /*tidx*/, double /*cm*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getVertV(vertex_id_t vidx) const {
    return MPI_ConditionalBcast<double>(
        0.0, MPI_DOUBLE, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setVertV(vertex_id_t /*vidx*/, double /*v*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getVertVClamped(vertex_id_t vidx) const {
    return MPI_ConditionalBcast<bool>(
        false, MPI_C_BOOL, RDEFmasterRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setVertVClamped(vertex_id_t /*vidx*/, bool /*cl*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setMembRes(solver::membrane_global_id /*midx*/,
                                    double /*ro*/,
                                    double /*vrev*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setMembPotential(solver::membrane_global_id /*midx*/, double /*v*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setMembCapac(solver::membrane_global_id /*midx*/, double /*cm*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setMembVolRes(solver::membrane_global_id /*midx*/, double /*ro*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

solver::vesicle_individual_id TetVesicleVesRaft::getVesicleNextIndex_() {
    solver::vesicle_individual_id next_idx(pVesicles_count);

    pVesicles_count += 1;

    return next_idx;
}

////////////////////////////////////////////////////////////////////////////////

solver::raft_individual_id TetVesicleVesRaft::getRaftNextIndex_() {
    solver::raft_individual_id next_idx(pRafts_count);

    pRafts_count += 1;

    return next_idx;
}

////////////////////////////////////////////////////////////////////////////////

solver::pointspec_individual_id TetVesicleVesRaft::getPointSpecNextIndex_() {
    ++pNextPointSpecUniqueID;
    ProgErrLogIf(pNextPointSpecUniqueID >= pRDEFminPointSpecUniqueID,
                 "VesRaft has run out of PointSpec IDs!\n ");
    return pNextPointSpecUniqueID;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_runRaft(double dt) {
    for (auto const& patch: pPatches) {
        patch->runRaft(dt);
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setVesicleSpecDiffD_(solver::vesicle_global_id vidx,
                                             solver::spec_global_id spec_gidx,
                                             double d) {
    // This should be at this level- stored by the comp object. Then
    // Whenever the vesicles need to add their qtables they check
    // whether the diff d has been defined in this comp or not.

    AssertLog(pVesSpecD.count(vidx) == 1);
    AssertLog(pVesSpecD[vidx].size() > spec_gidx.get());

    pVesSpecD[vidx][spec_gidx] = d;

    _recalcQtable_spec(vidx, spec_gidx, d);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_recalcQtable_spec(solver::vesicle_global_id vidx,
                                           solver::spec_global_id spec_gidx,
                                           double d) {
    if (d <= 0.0) {
        delete pQtables_spec[vidx][spec_gidx];
        pQtables_spec[vidx][spec_gidx] = nullptr;
        return;
    }

    double radius = statedef().vesicledef(vidx)->diameter() / 2.0;
    double tau = (2.0 * d * pVesicledt) / (radius * radius);

    if (pQtables_spec[vidx][spec_gidx] != nullptr) {
        pQtables_spec[vidx][spec_gidx]->reinit(pQtablesize_spec, tau);
    } else {
        pQtables_spec[vidx][spec_gidx] = new Qtable(pQtablesize_spec, tau, rng());
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setVesicleSurfaceLinkSpecSDiffD(solver::vesicle_global_id vidx,
                                                         solver::linkspec_global_id linkspec_gidx,
                                                         double d) {
    // This should be at this level- stored by the comp object. Then
    // Whenever the vesicles need to add their qtables they check
    // whether the diff d has been defined in this comp or not.
    AssertLog(pVesLinkSpecD.count(vidx) == 1);
    AssertLog(pVesLinkSpecD[vidx].size() > linkspec_gidx.get());

    pVesLinkSpecD[vidx][linkspec_gidx] = d;

    _recalcQtable_linkspec(vidx, linkspec_gidx, d);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_recalcQtable_linkspec(solver::vesicle_global_id vidx,
                                               solver::linkspec_global_id linkspec_gidx,
                                               double d) {
    if (d <= 0.0) {
        delete pQtables_linkspec[vidx][linkspec_gidx];
        pQtables_linkspec[vidx][linkspec_gidx] = nullptr;
        return;
    }

    double radius = statedef().vesicledef(vidx)->diameter() / 2.0;
    double tau = (2.0 * d * pVesicledt) / (radius * radius);

    if (pQtables_linkspec[vidx][linkspec_gidx] != nullptr) {
        pQtables_linkspec[vidx][linkspec_gidx]->reinit(pQtablesize_spec, tau);
    } else {
        pQtables_linkspec[vidx][linkspec_gidx] = new Qtable(pQtablesize_linkspec, tau, rng());
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_recalcQtables() {
    for (auto const& ves_to_qt: pQtables_spec) {
        solver::vesicle_global_id vidx = ves_to_qt.first;  // for clarity
        for (auto spec_gidx: solver::spec_global_id::range(ves_to_qt.second.size())) {
            double d = pVesSpecD[vidx][spec_gidx];
            _recalcQtable_spec(vidx, spec_gidx, d);
        }
    }

    for (auto const& ves_to_qt: pQtables_linkspec) {
        solver::vesicle_global_id vidx = ves_to_qt.first;  // for clarity
        for (auto linkspec_gidx: solver::linkspec_global_id::range(ves_to_qt.second.size())) {
            double d = pVesLinkSpecD[vidx][linkspec_gidx];
            _recalcQtable_linkspec(vidx, linkspec_gidx, d);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::getQPhiSpec_(solver::vesicle_global_id ves_gidx,
                                       solver::spec_global_id spec_gidx) {
    if (pQtables_spec[ves_gidx][spec_gidx] == nullptr) {
        return 0.0;
    } else {
        return pQtables_spec[ves_gidx][spec_gidx]->getPhi();
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::getQPhiLinkspec_(solver::vesicle_global_id ves_gidx,
                                           solver::linkspec_global_id linkspec_gidx) {
    if (pQtables_linkspec[ves_gidx][linkspec_gidx] == nullptr) {
        return 0.0;
    } else {
        return pQtables_linkspec[ves_gidx][linkspec_gidx]->getPhi();
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::addOverlap_(std::map<tetrahedron_global_id, double>& tets_overlap,
                                    Vesicle* ves) const {
    for (auto const& tet: tets_overlap) {
        tetrahedron_global_id tet_gidx = tet.first;
        tet_(tet_gidx)->changeOverlap(tet.second);
        tet_(tet_gidx)->addVesref(ves);
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::removeOverlap_(std::map<tetrahedron_global_id, double>& tets_overlap,
                                       Vesicle* ves) const {
    for (auto const& tet: tets_overlap) {
        tetrahedron_global_id tet_gidx = tet.first;
        tet_(tet_gidx)->changeOverlap(-(tet.second));
        tet_(tet_gidx)->removeVesref(ves);
    }
}

////////////////////////////////////////////////////////////////////////////////

LinkSpec* TetVesicleVesRaft::getLinkSpec_(solver::linkspec_individual_id linkspec_id) {
    for (auto const& vesicle: pVesicles) {
        for (auto const& linkspec: vesicle.second->getLinkSpecs()) {
            if (linkspec.first == linkspec_id) {
                return linkspec.second;
            }
        }
    }

    ProgErrLog("LinkSpec unknown!");
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::createPath(std::string const& id) {
    // Check id first
    if (pPaths.find(id) != pPaths.end()) {
        ArgErrLog("Path already exists with this ID.");
    }

    Path* path = new Path(id);

    pPaths[id] = path;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::addPathPoint(std::string const& path_name,
                                     uint point_id,
                                     const std::vector<double>& position) {
    if (pPaths.find(path_name) == pPaths.end()) {
        ArgErrLog("Path ID unknown.");
    }

    pPaths[path_name]->addPoint(point_id, {position[0], position[1], position[2]});
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::addPathBranch(std::string const& path_name,
                                      uint point_id,
                                      const std::map<uint, double>& dest_points) {
    if (pPaths.find(path_name) == pPaths.end()) {
        ArgErrLog("Path ID unknown.");
    }

    if (dest_points.size() == 0) {
        ArgErrLog("There must be at least one terminal point.");
    }

    pPaths[path_name]->addBranch(point_id, dest_points);
}

////////////////////////////////////////////////////////////////////////////////

std::map<std::string, std::map<uint, std::pair<std::vector<double>, std::map<uint, double>>>>
TetVesicleVesRaft::getAllPaths() const {
    std::map<std::string, std::map<uint, std::pair<std::vector<double>, std::map<uint, double>>>>
        ret;

    for (const auto& p: pPaths) {
        ret.emplace(p.first, p.second->getPathMap());
    }

    return ret;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_addPathVesicle(std::string const& path_name,
                                        solver::vesicle_global_id ves_idx,
                                        double speed,
                                        const std::map<solver::spec_global_id, uint>& spec_deps,
                                        const std::vector<double>& stoch_stepsize) {
    if (pPaths.find(path_name) == pPaths.end()) {
        ArgErrLog("Path ID unknown.");
    }

    if (speed <= 0.0) {
        ArgErrLog("Speed must be non-zero and positive.");
    }

    pPaths[path_name]->addVesicle(ves_idx, speed, spec_deps, stoch_stepsize);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::removeLinkSpecPair_(const LinkSpecPair* lsp) {
    // This method should only be called if the link species associated to the linkSpecPair
    // have not yet been deleted. If they have, calling this method could lead to segfaults
    // or incorrect result from `pLinkSpecPairs.find(lsp)`.
    auto it = pLinkSpecPairs.find(lsp);
    ProgErrLogIf(it == pLinkSpecPairs.end(), "Link spec pair could not be removed");
    delete lsp;
    pLinkSpecPairs.erase(it);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Path*> TetVesicleVesRaft::vesicleCrossedPaths_(const math::position_abs& ves_pos,
                                                           solver::vesicle_global_id ves_gidx,
                                                           double ves_rad) {
    std::vector<Path*> crossed_paths;
    for (auto const& path: pPaths) {
        if (path.second->crossedPath(ves_pos, ves_gidx, ves_rad)) {
            uint random_idx = rng()->get() % (crossed_paths.size() + 1);
            crossed_paths.emplace(crossed_paths.begin() + random_idx, path.second);
        }
    }

    return crossed_paths;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_addVesicleDiffusionGroup(
    solver::vesicle_global_id vidx,
    const std::vector<solver::comp_global_id>& comp_indices) {
    std::vector<CompVesRaft*> comps;

    for (auto const& cidx: comp_indices) {
        CompVesRaft* comp = getComp_(cidx);
        comps.emplace_back(comp);
    }

    for (auto const& comp: comps) {
        comp->addVesiclePermittedComps(vidx, comps);
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setVesSReacK(solver::vessreac_global_id vsridx, double kf) {
    AssertLog(kf >= 0.0);

    AssertLog(vsridx < statedef().countVesSReacs());

    statedef().vessreacdef(vsridx)->setKcst(kf);

    // Rates have changed
    //_updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setExocytosisK(solver::exocytosis_global_id exoidx, double kf) {
    AssertLog(kf >= 0.0);

    AssertLog(exoidx < statedef().countExocytosis());

    statedef().exocytosisdef(exoidx)->setKcst(kf);

    // Rates have changed
    //_updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getVesSReacExtent(solver::vessreac_global_id vsridx) const {
    AssertLog(vsridx < statedef().countVesSReacs());

    uint extent = statedef().vessreacdef(vsridx)->getExtent();
    return MPI_ConditionalBcast<uint>(
        extent, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getExocytosisExtent(solver::exocytosis_global_id exoidx) const {
    AssertLog(exoidx < statedef().countExocytosis());

    uint extent = statedef().exocytosisdef(exoidx)->getExtent();
    return MPI_ConditionalBcast<uint>(
        extent, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::ExocytosisEvent> TetVesicleVesRaft::_getExocytosisEvents(
    solver::exocytosis_global_id exoidx) {
    AssertLog(exoidx < statedef().countExocytosis());

    auto events = statedef().exocytosisdef(exoidx)->getEvents();
    MPI_ConditionalBcast(events,
                         dataTypeUtil.MPI_ExocytosisEventSync,
                         vesraftRank_World,
                         myRank_World,
                         syncOutput,
                         outputRank);
    return events;
}

////////////////////////////////////////////////////////////////////////////////

uint TetVesicleVesRaft::_getRaftEndocytosisExtent(
    solver::raftendocytosis_global_id rendoidx) const {
    AssertLog(rendoidx < statedef().countRaftEndocytosis());

    uint extent = statedef().raftendocytosisdef(rendoidx)->getExtent();
    return MPI_ConditionalBcast<uint>(
        extent, MPI_UNSIGNED, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::RaftEndocytosisEvent> TetVesicleVesRaft::_getRaftEndocytosisEvents(
    solver::raftendocytosis_global_id rendoidx) {
    AssertLog(rendoidx < statedef().countRaftEndocytosis());

    auto events = statedef().raftendocytosisdef(rendoidx)->getEvents();
    MPI_ConditionalBcast(events,
                         dataTypeUtil.MPI_RaftEndocytosisEventSync,
                         vesraftRank_World,
                         myRank_World,
                         syncOutput,
                         outputRank);
    return events;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setRaftEndocytosisK(solver::raftendocytosis_global_id rendoidx,
                                             double kcst) {
    AssertLog(rendoidx < statedef().countRaftEndocytosis());

    auto raftendo_def = statedef().raftendocytosisdef(rendoidx);
    raftendo_def->setKcst(kcst);

    // Also need to change all the RaftEndocytosis currently in existence
    for (auto const& raft_it: pRafts) {
        for (auto const& raftendo: raft_it.second->raftendos()) {
            if (raftendo->endodef() == raftendo_def) {
                raftendo->setKcst(kcst);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::_getSingleRaftRaftEndocytosisK(
    solver::raft_global_id ridx,
    solver::raft_individual_id raft_unique_index,
    solver::raftendocytosis_global_id rendoidx) const {
    double kcst{-1};
    AssertLog(rendoidx < statedef().countRaftEndocytosis());

    auto raftendo_def = statedef().raftendocytosisdef(rendoidx);

    auto const& raft_it = pRafts.find(raft_unique_index);
    if (raft_it == pRafts.end()) {
        CLOG(WARNING, "general_log") << "Raft unique id " << raft_unique_index << "unknown.\n";
    } else if (raft_it->second->idx() != ridx) {
        CLOG(WARNING, "general_log")
            << "Incorrect Raft type for raft unique ID " << raft_unique_index << ".\n ";
    } else {
        for (auto const& raftendo: raft_it->second->raftendos()) {
            if (raftendo->endodef() == raftendo_def) {
                kcst = raftendo->kcst();
            }
        }
        if (kcst == -1) {
            CLOG(WARNING, "general_log") << "Could not find raft endocytosis in raft unique ID "
                                         << raft_unique_index << ".\n ";
        }
    }

    return MPI_ConditionalBcast<double>(
        kcst, MPI_DOUBLE, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setSingleRaftRaftEndocytosisK(solver::raft_global_id ridx,
                                                       solver::raft_individual_id raft_unique_index,
                                                       solver::raftendocytosis_global_id rendoidx,
                                                       double k) {
    AssertLog(rendoidx < statedef().countRaftEndocytosis());

    auto raftendo_def = statedef().raftendocytosisdef(rendoidx);

    auto const& raft_it = pRafts.find(raft_unique_index);
    if (raft_it == pRafts.end()) {
        CLOG(WARNING, "general_log") << "Raft unique id " << raft_unique_index << "unknown.\n";
        return;
    }
    if (raft_it->second->idx() != ridx) {
        CLOG(WARNING, "general_log")
            << "Incorrect Raft type for raft unique ID " << raft_unique_index << ".\n ";
        return;
    }

    for (auto const& raftendo: raft_it->second->raftendos()) {
        if (raftendo->endodef() == raftendo_def) {
            raftendo->setKcst(k);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_setSingleRaftSReacActive(solver::raft_global_id ridx,
                                                  solver::raft_individual_id raft_unique_index,
                                                  solver::raftsreac_global_id rsreacidx,
                                                  bool active) {
    AssertLog(rsreacidx < statedef().countRaftSReacs());

    auto const& raft_it = pRafts.find(raft_unique_index);
    if (raft_it == pRafts.end()) {
        CLOG(WARNING, "general_log") << "Raft unique id " << raft_unique_index << "unknown.\n";
        return;
    }
    if (raft_it->second->idx() != ridx) {
        CLOG(WARNING, "general_log")
            << "Incorrect Raft type for raft unique ID " << raft_unique_index << ".\n ";
        return;
    }

    raft_it->second->setRaftSReacActive(rsreacidx, active);
}

////////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::_getSingleRaftSReacActive(solver::raft_global_id ridx,
                                                  solver::raft_individual_id raft_unique_index,
                                                  solver::raftsreac_global_id rsreacidx) const {
    AssertLog(rsreacidx < statedef().countRaftSReacs());

    auto const& raft_it = pRafts.find(raft_unique_index);
    if (raft_it == pRafts.end()) {
        CLOG(WARNING, "general_log") << "Raft unique id " << raft_unique_index << "unknown.\n";
        return true;
    }
    if (raft_it->second->idx() != ridx) {
        CLOG(WARNING, "general_log")
            << "Incorrect Raft type for raft unique ID " << raft_unique_index << ".\n ";
        return true;
    }

    bool active = raft_it->second->getRaftSReacActive(rsreacidx);

    return MPI_ConditionalBcast<bool>(
        active, MPI_C_BOOL, vesraftRank_World, myRank_World, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////
// ROI Data Access
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleVesRaft::getBatchTetSpecCounts(std::vector<index_t> const& tets,
                                                             std::string const& s) const {
    size_t ntets = tets.size();
    std::vector<double> counts(ntets, 0.0);
    getBatchTriSpecCountsNP(tets.data(), ntets, s, counts.data(), ntets);
    return counts;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleVesRaft::getBatchTriSpecCounts(std::vector<index_t> const& tris,
                                                             std::string const& s) const {
    size_t ntris = tris.size();
    std::vector<double> counts(ntris, 0.0);
    getBatchTriSpecCountsNP(tris.data(), ntris, s, counts.data(), ntris);
    return counts;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::getBatchTetSpecCountsNP(const index_t* indices,
                                                size_t input_size,
                                                std::string const& s,
                                                double* counts,
                                                size_t output_size) const {
    std::vector<double> local_counts(input_size, 0.0);
    MPI_ConditionalReduce<double>(
        local_counts.data(), counts, input_size, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::getBatchTriSpecCountsNP(const index_t* indices,
                                                size_t input_size,
                                                std::string const& s,
                                                double* counts,
                                                size_t output_size) const {
    std::vector<double> local_counts(input_size, 0.0);
    MPI_ConditionalReduce<double>(
        local_counts.data(), counts, input_size, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////
// ROI Data Access
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TetVesicleVesRaft::getROITetSpecCounts(const std::string& ROI_id,
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

std::vector<double> TetVesicleVesRaft::getROITriSpecCounts(const std::string& ROI_id,
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

void TetVesicleVesRaft::getROITetSpecCountsNP(const std::string& ROI_id,
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

void TetVesicleVesRaft::getROITriSpecCountsNP(const std::string& ROI_id,
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

double TetVesicleVesRaft::getROIVol(const std::string& ROI_id) const {
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

double TetVesicleVesRaft::getROIArea(const std::string& ROI_id) const {
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

double TetVesicleVesRaft::getROITetSpecCount(const std::vector<tetrahedron_global_id>& tetrahedrons,
                                             const std::string& s) const {
    return MPI_ConditionalReduce<double>(0.0, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::getROITriSpecCount(const std::vector<triangle_global_id>& triangles,
                                             const std::string& s) const {
    return MPI_ConditionalReduce<double>(0.0, MPI_DOUBLE, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::getROISpecCount(const std::string& ROI_id, std::string const& s) const {
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

void TetVesicleVesRaft::setROITriSpecCount(const std::vector<triangle_global_id>& /*unused*/,
                                           const std::string& /*unused*/,
                                           double /*unused*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setROITetSpecCount(const std::vector<tetrahedron_global_id>& s,
                                           const std::string& /*unused*/,
                                           double /*unused*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

// WEILIANG: Can we apply Sam's set count method (as in e.g. setCompCount)
// here?
void TetVesicleVesRaft::setROISpecCount(const std::string& /*ROI_id*/,
                                        std::string const& /*s*/,
                                        double /*count*/) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::getROISpecAmount(const std::string& ROI_id, std::string const& s) const {
    double count = getROISpecCount(ROI_id, s);
    return (count / math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setROISpecAmount(const std::string& /*ROI_id*/,
                                         std::string const& /*s*/,
                                         double /*amount*/) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////
// WEILIANG: Can we apply Sam's set count method (as in e.g. setCompCount)
// here?
void TetVesicleVesRaft::setROISpecConc(const std::string& /*ROI_id*/,
                                       std::string const& /*s*/,
                                       double /*conc*/) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double TetVesicleVesRaft::getROISpecConc(const std::string& ROI_id, const std::string& s) const {
    auto const& roi = _mesh()->rois.get<tetmesh::ROI_TET>(ROI_id);
    if (roi == _mesh()->rois.end<tetmesh::ROI_TET>()) {
        ArgErrLog("ROI check fail, please make sure the ROI stores correct elements.");
    }
    const double count = getROITetSpecCount(roi->second, s);
    double vol = getROIVol(ROI_id);
    return count / (1.0e3 * vol * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setROITriSpecClamped(const std::vector<triangle_global_id>& /*unused*/,
                                             const std::string& /*unused*/,
                                             bool /*unused*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setROITetSpecClamped(const std::vector<tetrahedron_global_id>& /*unused*/,
                                             std::string const& /*unused*/,
                                             bool /*unused*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setROISpecClamped(const std::string& /*ROI_id*/,
                                          std::string const& /*s*/,
                                          bool /*b*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setROIReacK(const std::string& /*ROI_id*/,
                                    std::string const& /*r*/,
                                    double /*kf*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setROISReacK(const std::string& /*ROI_id*/,
                                     std::string const& /*sr*/,
                                     double /*kf*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setROIDiffD(const std::string& /*ROI_id*/,
                                    std::string const& /*d*/,
                                    double /*dk*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setROIReacActive(const std::string& /*ROI_id*/,
                                         std::string const& /*r*/,
                                         bool /*a*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setROISReacActive(const std::string& /*ROI_id*/,
                                          std::string const& /*sr*/,
                                          bool /*a*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setROIDiffActive(const std::string& /*ROI_id*/,
                                         std::string const& /*d*/,
                                         bool /*act*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setROIVDepSReacActive(const std::string& /*ROI_id*/,
                                              std::string const& /*vsr*/,
                                              bool /*a*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetVesicleVesRaft::getROIReacExtent(const std::string& /*ROI_id*/,
                                                       std::string const& /*r*/) const {
    return MPI_ConditionalReduce<unsigned long long>(
        0, MPI_UNSIGNED_LONG_LONG, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::resetROIReacExtent(const std::string& /*ROI_id*/,
                                           std::string const& /*r*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetVesicleVesRaft::getROISReacExtent(const std::string& /*ROI_id*/,
                                                        std::string const& /*sr*/) const {
    return MPI_ConditionalReduce<unsigned long long>(
        0, MPI_UNSIGNED_LONG_LONG, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::resetROISReacExtent(const std::string& /*ROI_id*/,
                                            std::string const& /*sr*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long TetVesicleVesRaft::getROIDiffExtent(const std::string& /*ROI_id*/,
                                                       std::string const& /*d*/) const {
    return MPI_ConditionalReduce<unsigned long long>(
        0, MPI_UNSIGNED_LONG_LONG, MPI_SUM, syncOutput, outputRank);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::resetROIDiffExtent(const std::string& /*ROI_id*/,
                                           std::string const& /*s*/) {
    /* empty function */
}

////////////////////////////////////////////////////////////////////////////////

TetVesRaft* TetVesicleVesRaft::_getTet(tetrahedron_global_id tgidx) const {
    auto tet = pTets.at(tgidx);
    if (tet == nullptr) {
        std::ostringstream os;
        os << "Tetrahedron " << tgidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }
    return tet;
}

////////////////////////////////////////////////////////////////////////////////

TriVesRaft* TetVesicleVesRaft::_getTri(triangle_global_id tgidx) const {
    auto tri = pTris.at(tgidx);
    if (tri == nullptr) {
        std::ostringstream os;
        os << "Triangle " << tgidx << " has not been assigned to a patch.\n";
        ArgErrLog(os.str());
    }
    return tri;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::regTetPoolSync_(tetrahedron_global_id tet_gidx,
                                        solver::spec_global_id spec_gidx,
                                        uint count) {
    if (tet_(tet_gidx)->compdef()->specG2L(spec_gidx).unknown()) {
        std::ostringstream os;
        os << "Species is undefined in Tetrahedron: cannot register for sync.\n";
        ArgErrLog(os.str());
    }

    tetPoolCountSyncs_Vec.emplace_back(PoolCountSync{tet_gidx.get(), spec_gidx, count});
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::regTriPoolSync_(triangle_global_id tri_gidx,
                                        solver::spec_global_id spec_gidx,
                                        uint count) {
    if (tri_(tri_gidx)->patchdef()->specG2L(spec_gidx).unknown()) {
        std::ostringstream os;
        os << "Species is undefined in Triangle: cannot register for sync.\n";
        ArgErrLog(os.str());
    }

    triPoolCountSyncs_Vec.emplace_back(PoolCountSync{tri_gidx.get(), spec_gidx, count});
}

////////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_syncPools(SyncDirection direction) {
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

        for (auto const& entry: tetPoolCountSyncs_Vec) {
            TetVesRaft* tet = pTets[tetrahedron_global_id(entry.container_global_index)];
            AssertLog(tet != nullptr);  // TODO remove after testing
            solver::spec_local_id lsidx = tet->compdef()->specG2L(entry.spec_global_index);
            tet->setCount(lsidx, entry.count);
        }

        for (auto const& entry: triPoolCountSyncs_Vec) {
            TriVesRaft* tri = pTris[triangle_global_id(entry.container_global_index)];
            AssertLog(tri != nullptr);  // TODO remove after testing
            solver::spec_local_id lsidx = tri->patchdef()->specG2L(entry.spec_global_index);
            tri->setCount(lsidx, entry.count);
        }

        tetPoolCountSyncs_Vec.clear();
        triPoolCountSyncs_Vec.clear();

        break;
    }
    default:
        break;
    }
}

////////////////////////////////////////////////////////////////////////////

int TetVesicleVesRaft::_getTetHost(tetrahedron_global_id tgidx) const {
    auto host_result = tetHosts.find(tgidx);
    if (host_result == tetHosts.end()) {
        std::ostringstream os;
        os << "Tetrahedron " << tgidx << " has not been assigned to a host.\n";
        ArgErrLog(os.str());
    }
    return host_result->second;
}

////////////////////////////////////////////////////////////////////////////////

int TetVesicleVesRaft::_getTriHost(triangle_global_id tgidx) const {
    auto host_result = triHosts.find(tgidx);
    if (host_result == triHosts.end()) {
        std::ostringstream os;
        os << "Triangle " << tgidx << " has not been assigned to a host.\n";
        ArgErrLog(os.str());
    }
    return host_result->second;
}

////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_constructVesV2R() {
    tetV2R_Vec.clear();
    vesProxyV2R_Vec.clear();
    vesSurfSpecV2R_Vec.clear();
    vesLinkSpecV2R_Vec.clear();

    for (auto const& tet: pTets) {
        if (tet == nullptr) {
            continue;
        }

        tetrahedron_global_id tet_gidx = tet->idx();

        TetV2R tet_v2r{};

        double overlap = tet->getOverlap();

        tet_v2r.tetrahedron_global_index = tet_gidx;
        tet_v2r.overlap = overlap;

        tetV2R_Vec.emplace_back(tet_v2r);

        // No need to do vesproxies if full overlap- no vesicle surface
        if (tet->isFullOverlap()) {
            continue;
        }

        for (auto const& ves: tet->getVesrefs()) {
            VesProxyV2R ves_proxy_v2r{};

            solver::vesicle_global_id ves_gidx = ves->idx();
            solver::vesicle_individual_id ves_uniqueidx = ves->getUniqueIndex();
            const math::position_abs& ves_pos = ves->getPosition();

            // ves proxy structure has all the information it needs so can be
            // completed here
            ves_proxy_v2r.tetrahedron_global_index = tet_gidx;
            ves_proxy_v2r.vesicle_global_index = ves_gidx;
            ves_proxy_v2r.vesicle_individual_index = ves_uniqueidx;
            ves_proxy_v2r.contains_link = ves->containsLink();
            ves_proxy_v2r.vesicle_central_position = ves_pos;

            vesProxyV2R_Vec.emplace_back(ves_proxy_v2r);

            // ves surfspec structure, one per surfspec
            for (auto const& surf_specs: ves->getSurfSpecs()) {
                // surf_specs map .second holds vector of all species of this global
                // idx (global idx is surf_specs.first)
                for (auto const& surf_spec: surf_specs.second) {
                    if (surf_spec->getOverlapTet_gidx() != tet_gidx) {
                        continue;
                    }

                    VesSurfSpecV2R ves_surfspec_v2r{};

                    ves_surfspec_v2r.tetrahedron_global_index = tet_gidx;
                    ves_surfspec_v2r.vesicle_individual_index = ves_uniqueidx;
                    ves_surfspec_v2r.surface_spec_global_index = surf_specs.first;
                    ves_surfspec_v2r.surface_spec_individual_index = surf_spec->getUniqueIndex();
                    // The positions held in pointspec object are relative to vesicle
                    // centre
                    const auto& surf_spec_pos_rel = surf_spec->getPosCartesian();

                    ves_surfspec_v2r.surface_spec_position_abs[0] = surf_spec_pos_rel[0] +
                                                                    ves_pos[0];
                    ves_surfspec_v2r.surface_spec_position_abs[1] = surf_spec_pos_rel[1] +
                                                                    ves_pos[1];
                    ves_surfspec_v2r.surface_spec_position_abs[2] = surf_spec_pos_rel[2] +
                                                                    ves_pos[2];
                    vesSurfSpecV2R_Vec.emplace_back(ves_surfspec_v2r);
                }
            }

            // Inner species not required in this direction, because they don't
            // affect rates of anything within SSA

            // ves linkspec structure, one per linkspec
            for (auto const& link_spec: ves->getLinkSpecs()) {
                if (link_spec.second->getOverlapTet_gidx() != tet_gidx) {
                    continue;
                }

                VesLinkSpecV2R ves_linkspec_v2r{};

                ves_linkspec_v2r.tetrahedron_global_index = tet_gidx;
                ves_linkspec_v2r.vesicle_individual_index = ves_uniqueidx;
                ves_linkspec_v2r.linkspec_global_index = link_spec.second->getGidx();
                ves_linkspec_v2r.linkspec_individual_index = link_spec.first;
                ves_linkspec_v2r.linkspec_position_abs = link_spec.second->getPosCartesian_abs();

                vesLinkSpecV2R_Vec.emplace_back(ves_linkspec_v2r);
            }
        }
    }

    // send vesProxyV2R_Vec, vesSurfSpecV2R_Vec, vesLinkSpecV2R_Vec data to RDEF
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
}

////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_useVesR2V() {
    // receive data from RDEF
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

    // Drastic first step of clearing all the vesicle surface species. Does this
    // belong here??
    for (auto const& comp: pComps) {
        for (auto& ves_map: comp->getAllVesicles()) {
            for (auto const& ves: ves_map.second) {
                ves->clearSurfSpecs();
            }
        }
    }

    for (auto const& ves_proxy_r2v: vesProxyR2V_Vec) {
        const auto tet_gidx = ves_proxy_r2v.tetrahedron_global_index;

        TetVesRaft* tet = tet_(tet_gidx);

        if (tet == nullptr) {
            ProgErrLog("Tetrahedron not assigned to a compartment.\n");
        }

        const auto ves_individual_id = ves_proxy_r2v.vesicle_individual_index;

        int const immobility_upd = ves_proxy_r2v.immobility_update;
        const auto exo_applied_gidx = ves_proxy_r2v.exo_applied_global_index;

        bool found_ves = false;
        for (auto const& ves: tet->getVesrefs()) {
            if (ves->getUniqueIndex() == ves_individual_id) {
                found_ves = true;
                ves->updImmobility(immobility_upd);
                if (exo_applied_gidx.valid()) {
                    ves->scheduleExocytosis(exo_applied_gidx);
                }
                break;
            }
        }
        if (found_ves == false) {
            ProgErrLog("Vesicle not found.\n");
        }
    }

    for (auto const& ves_surfspec_r2v: vesSurfSpecR2V_Vec) {
        const auto tet_gidx = ves_surfspec_r2v.tetrahedron_global_index;

        if (tetHosts.find(tet_gidx) == tetHosts.end()) {
            ProgErrLog("Tetrahedron not assigned host.\n");
        }

        TetVesRaft* tet = tet_(tet_gidx);

        if (tet == nullptr) {
            ProgErrLog("Tetrahedron not assigned to a compartment.\n");
        }

        const auto ves_individual_id = ves_surfspec_r2v.vesicle_individual_index;

        const auto spec_gidx = ves_surfspec_r2v.surface_spec_global_index;
        auto spec_idx = ves_surfspec_r2v.surface_spec_individual_index;

        math::position_abs spec_pos_abs{ves_surfspec_r2v.surface_spec_position_abs[0],
                                        ves_surfspec_r2v.surface_spec_position_abs[1],
                                        ves_surfspec_r2v.surface_spec_position_abs[2]};

        // Re-assign point spec individual id attributed in RDEF ranks
        if (spec_idx >= pRDEFminPointSpecUniqueID) {
            spec_idx.set(getPointSpecNextIndex_().get());
        }

        bool found_ves = false;
        for (auto const& ves: tet->getVesrefs()) {
            if (ves->getUniqueIndex() == ves_individual_id) {
                found_ves = true;
                ves->addOneSurfSpec(spec_gidx, spec_idx, tet_gidx, spec_pos_abs);
                break;
            }
        }
        if (found_ves == false) {
            ProgErrLog("Vesicle not found.\n");
        }
    }

    for (auto const& ves_innerspec_r2v: vesInnerSpecR2V_Vec) {
        const auto tet_gidx = ves_innerspec_r2v.tetrahedron_global_index;

        if (tetHosts.find(tet_gidx) == tetHosts.end()) {
            ProgErrLog("Tetrahedron not assigned host.\n");
        }

        TetVesRaft* tet = tet_(tet_gidx);

        if (tet == nullptr) {
            ProgErrLog("Tetrahedron not assigned to a compartment.\n");
        }

        const auto ves_individual_id = ves_innerspec_r2v.vesicle_individual_index;

        const auto spec_gidx = ves_innerspec_r2v.inner_spec_global_index;
        uint count = ves_innerspec_r2v.count;

        bool found_ves = false;
        for (auto const& ves: tet->getVesrefs()) {
            if (ves->getUniqueIndex() == ves_individual_id) {
                found_ves = true;
                ves->incInnerSpecCount(spec_gidx, count);
                break;
            }
        }
        if (found_ves == false) {
            ProgErrLog("Vesicle not found.\n");
        }
    }

    //  vesLinkSpecPairR2V_Vec contains NEW LinkSpecs, created by VesBind. These must be added
    //  first, before going through all LinkSpecs in vesLinkSpecR2V_Vec
    // At this stage the IDs unique within RDEF cores, but not sequential. Ids will be
    // reassigned, with a map from the RDEF ID to the new, unique, sequential ID on the VesRaft
    // core
    std::map<solver::linkspec_individual_id, solver::linkspec_individual_id>
        ls_ind_id_rdef_to_vesraft;

    for (auto const& ves_linkspecpair_r2v: vesLinkSpecPairR2V_Vec) {
        const auto tet_gidx = ves_linkspecpair_r2v.tetrahedron_global_index;

        if (tetHosts.find(tet_gidx) == tetHosts.end()) {
            ProgErrLog("Tetrahedron not assigned host.\n");
        }

        TetVesRaft* tet = tet_(tet_gidx);

        if (tet == nullptr) {
            ProgErrLog("Tetrahedron not assigned to a compartment.\n");
        }

        const auto ves1_individual_id = ves_linkspecpair_r2v.vesicle1_individual_index;
        const auto ves2_individual_id = ves_linkspecpair_r2v.vesicle2_individual_index;

        const auto linkspec1_global_id = ves_linkspecpair_r2v.linkspec1_global_index;
        const auto linkspec2_global_id = ves_linkspecpair_r2v.linkspec2_global_index;

        solver::LinkSpecdef* linkspec1_def = statedef().linkspecdef(linkspec1_global_id);
        solver::LinkSpecdef* linkspec2_def = statedef().linkspecdef(linkspec2_global_id);

        const auto linkspec1_individual_id_rdef = ves_linkspecpair_r2v.linkspec1_individual_index;
        const auto linkspec2_individual_id_rdef = ves_linkspecpair_r2v.linkspec2_individual_index;

        const auto linkspec1_individual_id = pNextLinkSpecUniqueID++;
        const auto linkspec2_individual_id = pNextLinkSpecUniqueID++;

        AssertLog(pNextLinkSpecUniqueID.get() < std::numeric_limits<index_t>::max() / 2);

        AssertLog(ls_ind_id_rdef_to_vesraft.count(linkspec1_individual_id_rdef) == 0);
        AssertLog(ls_ind_id_rdef_to_vesraft.count(linkspec2_individual_id_rdef) == 0);

        ls_ind_id_rdef_to_vesraft[linkspec1_individual_id_rdef] = linkspec1_individual_id;
        ls_ind_id_rdef_to_vesraft[linkspec2_individual_id_rdef] = linkspec2_individual_id;

        double min_length = ves_linkspecpair_r2v.min_length;
        double max_length = ves_linkspecpair_r2v.max_length;

        // At this stage, the vesicles should both be present in the one
        // tetrahedron
        Vesicle* ves1 = nullptr;
        Vesicle* ves2 = nullptr;

        for (auto const& ves: tet->getVesrefs()) {
            if (ves->getUniqueIndex() == ves1_individual_id) {
                ves1 = ves;
                continue;
            }
            if (ves->getUniqueIndex() == ves2_individual_id) {
                ves2 = ves;
                continue;
            }
        }
        if (ves1 == nullptr || ves2 == nullptr) {
            ProgErrLog("Failed to find linkspecs in vesicle");
        }

        math::position_abs ves1_pos = ves1->getPosition();
        math::position_rel_to_ves linkspec1_pos_rel{
            ves_linkspecpair_r2v.linkspec1_position_abs[0] - ves1_pos[0],
            ves_linkspecpair_r2v.linkspec1_position_abs[1] - ves1_pos[1],
            ves_linkspecpair_r2v.linkspec1_position_abs[2] - ves1_pos[2]};

        auto linkspec1 =
            new LinkSpec(linkspec1_def, linkspec1_individual_id, ves1, tet_gidx, linkspec1_pos_rel);

        ves1->addLinkSpec(linkspec1_individual_id, linkspec1);

        math::position_abs ves2_pos = ves2->getPosition();
        math::position_rel_to_ves linkspec2_pos_rel{
            ves_linkspecpair_r2v.linkspec2_position_abs[0] - ves2_pos[0],
            ves_linkspecpair_r2v.linkspec2_position_abs[1] - ves2_pos[1],
            ves_linkspecpair_r2v.linkspec2_position_abs[2] - ves2_pos[2]};

        auto linkspec2 =
            new LinkSpec(linkspec2_def, linkspec2_individual_id, ves2, tet_gidx, linkspec2_pos_rel);

        ves2->addLinkSpec(linkspec2_individual_id, linkspec2);

        auto linkspecpair =
            new LinkSpecPair(linkspec1, linkspec2, ves1, ves2, min_length, max_length);
        // NOTE LinkSpecPairs are cleaned up in runvesicle as vesunbinding occurs

        linkspec1->addLinkSpecPair(linkspecpair);
        linkspec2->addLinkSpecPair(linkspecpair);

        pLinkSpecPairs.insert(linkspecpair);
    }

    //  vesLinkSpecR2V_Vec contains LinkSpecs that should already exist on this solver (RDEF
    // can't destroy them) however, the ID may have changed, so we update them all.

    for (auto const& ves_linkspec_r2v: vesLinkSpecR2V_Vec) {
        const auto tet_gidx = ves_linkspec_r2v.tetrahedron_global_index;

        if (tetHosts.find(tet_gidx) == tetHosts.end()) {
            ProgErrLog("Tetrahedron not assigned host.\n");
        }

        TetVesRaft* tet = tet_(tet_gidx);

        if (tet == nullptr) {
            ProgErrLog("Tetrahedron not assigned to a compartment.\n");
        }

        const auto linkspec_global_id = ves_linkspec_r2v.linkspec_global_index;

        solver::LinkSpecdef* linkspec_def = statedef().linkspecdef(linkspec_global_id);

        const auto linkspec_individual_id_rdef = ves_linkspec_r2v.linkspec_individual_index;
        solver::linkspec_individual_id linkspec_individual_id = linkspec_individual_id_rdef;
        auto it = ls_ind_id_rdef_to_vesraft.find(linkspec_individual_id_rdef);
        if (it != ls_ind_id_rdef_to_vesraft.end()) {
            linkspec_individual_id = it->second;
        }

        const auto ves_individual_id = ves_linkspec_r2v.vesicle_individual_index;

        bool found_ves = false;
        for (auto const& ves: tet->getVesrefs()) {
            if (ves->getUniqueIndex() == ves_individual_id) {
                found_ves = true;

                ves->updateLinkSpec(linkspec_individual_id, linkspec_def);
                break;
            }
        }
        if (found_ves == false) {
            ProgErrLog("Vesicle not found.\n");
        }
    }
}

////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_constructRaftV2R() {
    // TODO CHECK FLAG, AND ONLY DO ALL THIS IF NECESSARY. IF NOT, SIMPLY SEND
    // FLAG TO CLIENTS

    raftProxyV2R_Vec.clear();
    raftSpecV2R_Vec.clear();
    raftSReacInactiveV2R_Vec.clear();

    for (auto const& patch: pPatches) {
        for (auto const& raft_type: patch->getAllRafts()) {
            solver::raft_global_id raft_gidx = raft_type.first;

            for (auto const& raft: raft_type.second) {
                solver::raft_individual_id raft_uniqueidx = raft->getUniqueIndex();

                for (triangle_global_id tri_gidx: raft->getOverlapVec()) {
                    RaftProxyV2R raft_proxy_v2r{};

                    raft_proxy_v2r.triangle_global_index = tri_gidx;
                    raft_proxy_v2r.raft_global_index = raft_gidx;
                    raft_proxy_v2r.raft_individual_index = raft_uniqueidx;

                    raftProxyV2R_Vec.emplace_back(raft_proxy_v2r);

                    for (solver::raftsreac_global_id rsreac_gidx: raft->raftsreacsinactive()) {
                        RaftSReacInactiveV2R raft_srinactive_v2r{};

                        raft_srinactive_v2r.triangle_global_index = tri_gidx;
                        raft_srinactive_v2r.raft_individual_index = raft_uniqueidx;
                        raft_srinactive_v2r.raftsreac_global_index = rsreac_gidx;

                        raftSReacInactiveV2R_Vec.emplace_back(raft_srinactive_v2r);
                    }
                }

                // Separate loop for species, because tris are returned, not given
                for (auto spec_gidx: solver::spec_global_id::range(statedef().countSpecs())) {
                    for (auto const& tri_counts: raft->getTriSpecCounts(spec_gidx, rng())) {
                        RaftSpecV2R raft_spec_v2r{};

                        raft_spec_v2r.triangle_global_index = tri_counts.first;
                        raft_spec_v2r.raft_individual_index = raft_uniqueidx;
                        raft_spec_v2r.spec_global_index = spec_gidx;
                        raft_spec_v2r.count = tri_counts.second;

                        raftSpecV2R_Vec.emplace_back(raft_spec_v2r);
                    }
                }
            }
        }
    }

    //  send raftProxyV2R_Vec and raftSpecV2R_Vec data to RDEF
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
}

////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::_useRaftR2V() {
    // receive data from RDEF

    MPI_GatherVec<RaftGenCountR2V>(triRaftGenR2V_Vec,
                                   dataTypeUtil.MPI_RaftGenCountR2V,
                                   vesraftRank_World,
                                   myRank_World,
                                   nHosts_World,
                                   MPI_COMM_WORLD);

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

    // Clear all the raft specs
    for (auto const& patch: pPatches) {
        for (auto& raft_map: patch->getAllRafts()) {
            for (auto const& raft: raft_map.second) {
                raft->clearSpecs();
            }
        }
    }

    for (auto const& raft_proxy_r2v: raftProxyR2V_Vec) {
        const auto tri_gidx = raft_proxy_r2v.triangle_global_index;

        TriVesRaft* tri = tri_(tri_gidx);

        if (tri == nullptr) {
            ProgErrLog("Triangle not assigned to a patch.\n");
        }

        const auto raft_individual_id = raft_proxy_r2v.raft_individual_index;
        int const& mobility_update = raft_proxy_r2v.immobility_update;

        bool found_raft = false;
        for (auto const& raft: tri->getRaftrefs()) {
            if (raft->getUniqueIndex() == raft_individual_id) {
                found_raft = true;
                raft->updImmobility(mobility_update);
                break;
            }
        }
        if (found_raft == false) {
            ProgErrLog("Raft not found.\n");
        }
    }

    for (auto const& raft_spec_r2v: raftSpecR2V_Vec) {
        const auto tri_gidx = raft_spec_r2v.triangle_global_index;

        TriVesRaft* tri = tri_(tri_gidx);

        if (tri == nullptr) {
            ProgErrLog("Triangle not assigned to a patch.\n");
        }

        const auto raft_individual_id = raft_spec_r2v.raft_individual_index;
        const auto spec_gidx = raft_spec_r2v.spec_global_index;

        uint const& count = raft_spec_r2v.count;

        bool found_raft = false;
        for (auto const& raft: tri->getRaftrefs()) {
            if (raft->getUniqueIndex() == raft_individual_id) {
                found_raft = true;
                raft->incSpecCountByGidx(spec_gidx, count);
                break;
            }
        }
        if (found_raft == false) {
            ProgErrLog("Raft not found.\n");
        }
    }

    for (auto const& raftgen_count_r2v: triRaftGenR2V_Vec) {
        const auto tri_gidx = raftgen_count_r2v.triangle_global_index;
        TriVesRaft* tri = tri_(tri_gidx);
        if (tri == nullptr) {
            ProgErrLog("Triangle not assigned to a patch.\n");
        }

        tri->addRaftGen(raftgen_count_r2v.raftgen_global_index, raftgen_count_r2v.count);
    }
}

////////////////////////////////////////////////////////////////////////////

void TetVesicleVesRaft::setOutputSync(bool enable_sync, int output_rank) {
    syncOutput = enable_sync;
    outputRank = output_rank;
}

////////////////////////////////////////////////////////////////////////////

bool TetVesicleVesRaft::getOutputSyncStatus() const {
    return syncOutput;
}

////////////////////////////////////////////////////////////////////////////

int TetVesicleVesRaft::getOutputSyncRank() const {
    return outputRank;
}

}  // namespace steps::mpi::tetvesicle
