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

#include "mpi/tetvesicle/patch_vesraft.hpp"

// STEPS headers.
#include "math/point.hpp"
#include "math/triangle.hpp"
#include "mpi/tetvesicle/endocytosis.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/reac.hpp"
#include "mpi/tetvesicle/tetvesicle_vesraft.hpp"
#include "mpi/tetvesicle/tri_vesraft.hpp"
#include "solver/raftdef.hpp"
#include "solver/raftgendef.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"
#include "util/distribute.hpp"

namespace steps::mpi::tetvesicle {

////////////////////////////////////////////////////////////////////////////////

PatchVesRaft::PatchVesRaft(solver::Patchdef* patchdef,
                           tetmesh::Tetmesh* mesh,
                           TetVesicleVesRaft* vesraft)
    : pPatchdef(patchdef)
    , pArea(0.0)
    , pMesh(mesh)
    , pVesRaft(vesraft) {
    AssertLog(pPatchdef != nullptr);
    pRNG = pPatchdef->statedef()->rng();
}

////////////////////////////////////////////////////////////////////////////////

PatchVesRaft::~PatchVesRaft() {
    for (auto* endo: pEndosVec) {
        delete endo;
    }
}

////////////////////////////////////////////////////////////////////////////////

void PatchVesRaft::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pArea);
    std::map<solver::raft_global_id, uint> rafts;
    for (auto const& raft: pRafts) {
        rafts[raft.first] = raft.second.size();
    }

    util::checkpoint(cp_file, rafts);

    for (auto const& raft: pRafts) {
        for (auto const& raft_p: raft.second) {
            raft_p->checkpoint(cp_file);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void PatchVesRaft::restore(std::fstream& cp_file) {
    util::compare(cp_file, pArea);
    std::map<solver::raft_global_id, uint> rafts;
    util::restore(cp_file, rafts);

    for (auto const& raft: rafts) {
        for (uint i = 0; i < raft.second; ++i) {
            // Constructor does the restore
            solver::Raftdef* raftdef = def()->statedef()->raftdef(raft.first);
            Raft* r = new Raft(raftdef, this, cp_file);

            pRafts[raft.first].push_back(r);

            triangle_global_id tri_central = r->tri_central()->idx();
            addTriRaftsRefs(tri_central, r);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void PatchVesRaft::reset() {
    def()->reset();

    for (auto& rit: pRafts) {
        for (auto r: rit.second) {
            delete r;
        }
        rit.second.clear();
    }

    for (auto const& endo: pEndosMap) {
        for (auto const& zone: endo.second) {
            zone.second->resetExtent();
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void PatchVesRaft::setupRafts() {
    for (auto const& r: def()->statedef()->rafts()) {
        solver::raft_global_id raftdef_gidx = r->gidx();

        if (pRafts.find(raftdef_gidx) != pRafts.end()) {
            std::ostringstream os;
            os << "Raft initialisation error. Index " << raftdef_gidx << " is not unique.\n";
            ProgErrLog(os.str());
        }

        std::list<Raft*> raft;

        pRafts[raftdef_gidx] = raft;

        // Set up the overlap tables
        for (auto const& tri: pTris) {
            triangle_global_id tri_globalidx = tri->idx();
            std::vector<triangle_global_id> overlap;

            // Start with the actual triangle itself
            overlap.push_back(tri_globalidx);

            // Call it overlap if it overlaps with one of the vertices. This will
            // miss a certain overlap types where overlap is with 'bar' but no
            // vertices. Potential to improve in the future

            math::point3d tri_baryc = mesh()->_getTriBarycenter(tri_globalidx);

            for (auto const& trineighb: pTris) {
                if (tri == trineighb) {
                    continue;
                }
                triangle_global_id trineighb_globalidx = trineighb->idx();

                const auto& trineighb_verts = mesh()->getTri(trineighb_globalidx);

                math::point3d tnv0(mesh()->_getVertex(vertex_id_t(trineighb_verts[0])));
                math::point3d tnv1(mesh()->_getVertex(vertex_id_t(trineighb_verts[1])));
                math::point3d tnv2(mesh()->_getVertex(vertex_id_t(trineighb_verts[2])));

                double min_distance = std::min({math::distance(tri_baryc, tnv0),
                                                math::distance(tri_baryc, tnv1),
                                                math::distance(tri_baryc, tnv2)});

                if (min_distance < r->diameter() / 2.0) {
                    overlap.push_back(trineighb_globalidx);
                }
            }
            pTriOverlap[raftdef_gidx][tri_globalidx] = overlap;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void PatchVesRaft::setupEndocyticZones() {
    for (auto* zone: def()->endocyticZones()) {
        std::vector<TriVesRaft*> tris_p;

        for (auto const& tidx: zone->tris()) {
            auto local = pTriidcs_G_to_L.find(tidx);
            if (local == pTriidcs_G_to_L.end()) {
                ArgErrLog("Triangle does not belong to this patch.");
            }
            tris_p.push_back(pTris[local->second]);
        }

        if (pZones.count(zone->name()) != 0) {
            std::ostringstream os;
            os << "Zone ID " << zone->name() << " is already in use.\n";
            ProgErrLog(os.str());
        }

        pZones[zone->name()] = tris_p;

        // Add all endocytosis reactions that are in the patch
        for (auto endoidx: solver::endocytosis_local_id::range(def()->countEndocytosis())) {
            if (pEndosMap.count(endoidx) == 0) {
                pEndosMap[endoidx] = std::map<std::string, Endocytosis*>();
            }
            // Create the endocytosis rule
            solver::Endocytosisdef* endodef = def()->endocytosisdef(endoidx);
            auto* endo = new Endocytosis(endodef, tris_p);
            pEndosVec.push_back(endo);

            pEndosMap[endoidx][zone->name()] = endo;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void PatchVesRaft::addTri(TriVesRaft* tri) {
    AssertLog(tri->patchdef() == def());

    index_t lidx_uint = pTris.size();
    triangle_local_id lidx(lidx_uint);
    triangle_global_id gidx(tri->idx());

    pTriidcs_L_to_G.emplace(lidx, gidx);
    pTriidcs_G_to_L.emplace(gidx, lidx);

    pTris.container().push_back(tri);
    pArea += tri->area();

    tri->setPatchVesRaft(this);
}

////////////////////////////////////////////////////////////////////////////////

TriVesRaft* PatchVesRaft::pickTriByArea(double rand01) const {
    if (countTris() == 0) {
        return nullptr;
    }
    if (countTris() == 1) {
        return pTris[triangle_local_id(0)];
    }
    double accum = 0.0;
    double selector = rand01 * area();

    for (const auto& t: tris()) {
        accum += t->area();
        if (selector <= accum) {
            return t;
        }
    }

    return tris().back();
}

////////////////////////////////////////////////////////////////////////////////

Endocytosis* PatchVesRaft::getEndocytosis(solver::endocytosis_local_id endoidx,
                                          const std::string& zone) {
    if (pZones.find(zone) == pZones.end()) {
        std::ostringstream os;
        os << "\nZone ID " << zone << " is unknown.\n";
        ProgErrLog(os.str());
    }

    auto endosIt = pEndosMap.find(endoidx);
    if (endosIt == pEndosMap.end()) {
        std::ostringstream os;
        os << "\nEndocytosis has not been added to zone.\n";
        ProgErrLog(os.str());
    }
    auto zoneIt = endosIt->second.find(zone);
    if (zoneIt == endosIt->second.end()) {
        std::ostringstream os;
        os << "\nEndocytosis has not been added to zone.\n";
        ProgErrLog(os.str());
    }
    return zoneIt->second;
}

////////////////////////////////////////////////////////////////////////////////

solver::raft_individual_id PatchVesRaft::addRaft(solver::Raftdef* raftdef, TriVesRaft* tri) {
    solver::raft_global_id rgidx = raftdef->gidx();
    if (pRafts.count(rgidx) == 0) {
        std::ostringstream os;
        os << "Raft addition error. Index " << rgidx << " is unknown in patch.\n";
        ProgErrLog(os.str());
    }

    triangle_global_id tri_gidx = tri->idx();

    math::position_abs rpos;

    bool avail = getPos(&rpos, raftdef, tri_gidx);
    if (avail == false) {
        return {};
    }
    solver::raft_individual_id unique_index = pVesRaft->getRaftNextIndex_();

    Raft* raft = new Raft(raftdef, this, tri, rpos, unique_index);

    raft->setupKProcs();

    pRafts[rgidx].emplace_back(raft);

    addTriRaftsRefs(tri_gidx, raft);

    return unique_index;
}

////////////////////////////////////////////////////////////////////////////////

solver::raft_individual_id PatchVesRaft::addRaft(solver::raft_global_id ridx, TriVesRaft* tri) {
    if (pRafts.count(ridx) == 0) {
        std::ostringstream os;
        os << "Raft index " << ridx << " is unknown in patch.\n";
        ProgErrLog(os.str());
    }

    solver::Raftdef* raftdef = def()->statedef()->raftdef(ridx);

    solver::raft_individual_id added_raft;

    uint attempts = 0;
    while (added_raft.unknown()) {
        attempts++;
        if (attempts > 1000) {
            return {};
        }
        added_raft = addRaft(raftdef, tri);
    }

    return added_raft;
}

////////////////////////////////////////////////////////////////////////////////

void PatchVesRaft::removeOneRaft(Raft* raft) {
    solver::raft_global_id raft_gidx = raft->def()->gidx();

    triangle_global_id raft_centraltri_gidx = raft->tri_central()->idx();

    removeTriRaftRefs(raft_centraltri_gidx, raft);

    // Finally remove the vesicle from list and delete it.
    pRafts[raft_gidx].remove(raft);

    delete (raft);
}

////////////////////////////////////////////////////////////////////////////////

void PatchVesRaft::setRaftCount(solver::raft_global_id ridx, uint count) {
    if (pRafts.count(ridx) == 0) {
        std::ostringstream os;
        os << "Raft index " << ridx << " is unknown in patch.\n";
        ProgErrLog(os.str());
    }

    // No need to call the remove function to ensure the references are removed
    // because raft destructor will do that
    for (auto r: pRafts[ridx]) {
        delete r;
    }
    pRafts[ridx].clear();

    solver::Raftdef* raftdef = def()->statedef()->raftdef(ridx);

    for (uint i = 0; i < count; ++i) {
        solver::raft_individual_id added_raft;

        uint attempts = 0;
        while (added_raft.unknown()) {
            attempts++;
            if (attempts > 10000) {
                ArgErrLog("Unable to set count of rafts: too many iterations.");
            }
            TriVesRaft* tri = *(util::random_element(tris().begin(), tris().end(), rng()));
            added_raft = addRaft(raftdef, tri);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

uint PatchVesRaft::getRaftCount(solver::raft_global_id ridx) const {
    auto map_it = pRafts.find(ridx);

    if (map_it == pRafts.end()) {
        std::ostringstream os;
        os << "Raft index " << ridx << " is unknown in patch.\n";
        ProgErrLog(os.str());
    }

    return map_it->second.size();
}

////////////////////////////////////////////////////////////////////////////////

void PatchVesRaft::setRaftCount(solver::raft_global_id ridx, TriVesRaft* tri, uint count) {
    if (pRafts.count(ridx) == 0) {
        std::ostringstream os;
        os << "Raft index " << ridx << " is unknown in patch.\n";
        ProgErrLog(os.str());
    }

    // Need to remove any existing Rafts and then replace them
    pRafts[ridx].remove_if([&tri](auto& raft) {
        if (raft->tri_central() == tri) {
            delete raft;
            return true;
        }
        return false;
    });

    solver::Raftdef* raftdef = def()->statedef()->raftdef(ridx);

    for (uint i = 0; i < count; ++i) {
        solver::raft_individual_id added_raft;
        uint attempts = 0;

        // addRaft will return RAFT_IIDX_UNDEFINED if unsuccessful
        while (added_raft.unknown()) {
            attempts++;
            if (attempts > 1000) {
                ArgErrLog("Unable to set count of rafts: too many iterations.");
            }
            added_raft = addRaft(raftdef, tri);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

uint PatchVesRaft::getRaftCount(solver::raft_global_id ridx, TriVesRaft* tri) const {
    auto map_it = pRafts.find(ridx);

    if (map_it == pRafts.end()) {
        std::ostringstream os;
        os << "Raft index " << ridx << " is unknown in patch.\n";
        ProgErrLog(os.str());
    }

    uint count = 0;
    for (auto const& raft: map_it->second) {
        if (raft->tri_central() == tri) {
            count += 1;
        }
    }

    return count;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::raft_individual_id> PatchVesRaft::getRaftIndices(
    solver::raft_global_id ridx) const {
    std::vector<solver::raft_individual_id> raft_indices;

    auto rit = pRafts.find(ridx);
    if (rit != pRafts.end()) {
        for (auto const& r: rit->second) {
            raft_indices.emplace_back(r->getUniqueIndex());
        }
    } else {
        std::ostringstream os;
        os << "Raft index " << ridx << " is unknown in patch.\n";
        ProgErrLog(os.str());
    }

    return raft_indices;
}

////////////////////////////////////////////////////////////////////////////////

std::map<solver::raft_individual_id, uint> PatchVesRaft::getRaftSpecCountMap(
    solver::raft_global_id ridx,
    solver::spec_global_id spec_gidx) const {
    auto map_it = pRafts.find(ridx);

    if (map_it == pRafts.end()) {
        std::ostringstream os;
        os << "Raft get species error. Index " << ridx << " is unknown in patch.\n";
        ProgErrLog(os.str());
    }

    std::map<solver::raft_individual_id, uint> spec_counts;

    for (auto const& raft: map_it->second) {
        spec_counts[raft->getUniqueIndex()] = raft->getSpecCountByGidx(spec_gidx);
    }

    return spec_counts;
}

////////////////////////////////////////////////////////////////////////////////

uint PatchVesRaft::getRaftSpecCount(solver::raft_global_id ridx,
                                    solver::spec_global_id spec_gidx) const {
    auto map_it = pRafts.find(ridx);

    if (map_it == pRafts.end()) {
        std::ostringstream os;
        os << "Raft get species error. Index " << ridx << " is unknown in patch.\n";
        ProgErrLog(os.str());
    }

    uint spec_counts = 0;

    for (auto const& raft: map_it->second) {
        spec_counts += raft->getSpecCountByGidx(spec_gidx);
    }

    return spec_counts;
}

////////////////////////////////////////////////////////////////////////////////

void PatchVesRaft::addTriRaftsRefs(triangle_global_id centraltri_gidx, Raft* raft) {
    solver::raft_global_id raft_gidx = raft->def()->gidx();

    // tri overlap holds global indices
    for (auto const& tri_gidx: pTriOverlap[raft_gidx][centraltri_gidx]) {
        triangle_local_id tri_lidx = triidx_G_to_L(tri_gidx);
        pTris[tri_lidx]->addRaftref(raft);
    }
}

////////////////////////////////////////////////////////////////////////////////

void PatchVesRaft::removeTriRaftRefs(triangle_global_id centraltri_gidx, Raft* raft) {
    solver::raft_global_id raft_gidx = raft->def()->gidx();

    // tri overlap holds global indices
    for (auto const& tri_gidx: pTriOverlap[raft_gidx][centraltri_gidx]) {
        triangle_local_id tri_lidx = triidx_G_to_L(tri_gidx);
        pTris[tri_lidx]->removeRaftref(raft);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool PatchVesRaft::getPos(math::position_abs* rpos,
                          solver::Raftdef* raftdef,
                          triangle_global_id central_tri_gidx,
                          Raft* raft) const {
    const auto& tri_verts = mesh()->getTri(central_tri_gidx);

    math::point3d v0 = mesh()->_getVertex(vertex_id_t(tri_verts[0]));
    math::point3d v1 = mesh()->_getVertex(vertex_id_t(tri_verts[1]));
    math::point3d v2 = mesh()->_getVertex(vertex_id_t(tri_verts[2]));

    math::position_abs pos;
    uint failures = 0;

    while (failures < 10)  // not really necessary, but to avoid while(1)
    {
        // Need 3 random numbers
        double s = rng()->getUnfII();
        double t = rng()->getUnfII();

        pos = math::tri_ranpnt(v0, v1, v2, s, t);

        double diam = raftdef->diameter();

        bool failed = false;

        for (auto const& rit: pRafts) {
            for (auto const& r: rit.second) {
                if (r == raft) {
                    continue;
                }
                math::position_abs r_pos = r->getPosition();
                double diam2 = r->getDiam();
                double distance2 = pow(pos[0] - r_pos[0], 2) + pow(pos[1] - r_pos[1], 2) +
                                   pow(pos[2] - r_pos[2], 2);
                if (sqrt(distance2) < (diam + diam2) / 2.0) {
                    failed = true;
                    break;
                }
            }
            if (failed) {
                break;
            }
        }
        if (failed) {
            failures += 1;
            if (failures == 10) {
                return false;
            }
        } else {
            break;
        }
    }

    // Copy the position
    (*rpos)[0] = pos[0];
    (*rpos)[1] = pos[1];
    (*rpos)[2] = pos[2];

    return true;
}

////////////////////////////////////////////////////////////////////////////////

void PatchVesRaft::runRaft(double dt) {
    std::vector<Raft*> raft_sched_del;

    // First we are going to give endo and dis a chance to fire
    for (auto const& rafts: pRafts) {
        for (auto const& raft: rafts.second) {
            bool applied = raft->applyEndoAndDis(dt);
            if (applied) {
                raft_sched_del.emplace_back(raft);
            }
        }
    }

    // raft_sched_del now holds list of rafts that underwent endocytosis
    // or dissolution, so they need to be deleted
    for (auto& raft: raft_sched_del) {
        removeOneRaft(raft);
    }

    // Can now go ahead with difusion of remaining rafts
    for (auto const& rit: pRafts) {
        for (auto const& r: rit.second) {
            if (r->getImmobility() == 0) {
                double cumulative_dt = 0.0;

                while (true) {
                    double next_dt = rng()->getExp(r->getScaledDcst());
                    if (cumulative_dt + next_dt > dt) {
                        break;
                    }
                    cumulative_dt += next_dt;

                    // Diffusion gonna happen if we can find a position
                    TriVesRaft* tri = r->selectDirectionTri(rng()->getUnfII());

                    math::position_abs rpos;

                    bool avail = getPos(&rpos, r->def(), tri->idx(), r);
                    if (avail == false) {
                        continue;  // Diffusion failed. Do nothing.
                    }

                    // Diffusion can happen. Remove the references, add the new ones
                    triangle_global_id raft_centraltri_gidx = r->tri_central()->idx();
                    removeTriRaftRefs(raft_centraltri_gidx, r);

                    r->setPosition(rpos, tri);

                    addTriRaftsRefs(tri->idx(), r);
                }
            }
        }
    }

    // A convenient place to give Endocytosis a chance, even though it's not
    // strictly part of the raft routines.

    for (auto& endo: pEndosVec) {
        double rate = endo->rate();
        if (rate > 0.0) {
            double endo_dt = rng()->getExp(rate);
            if (endo_dt < dt) {
                endo->apply(solverVesRaft());
                // it's possible endo fails if no space for vesicle.
            }
        }
    }

    // Here need to check raftgens and apply if necessary

    // Set up the overlap tables
    for (auto const& tri: pTris) {
        for (auto const& raftgen: tri->getAppliedRaftGens()) {
            // Prefetch some stuff
            solver::RaftGendef* raftgendef = def()->statedef()->raftgendef(raftgen.first);
            AssertLog(raftgendef != nullptr);
            solver::Raftdef* raftdef = raftgendef->raftdef();
            AssertLog(raftdef != nullptr);
            solver::raft_global_id rgidx = raftdef->gidx();

            // Can be more than one raft now
            for (uint count = 0; count < raftgen.second; ++count) {
                const auto& lhs_s_vec = raftgendef->lhs_S();  // OK- global species

                solver::raft_individual_id raft_unique_index = addRaft(raftdef, tri);
                if (raft_unique_index.unknown()) {
                    // We could not create the raft.
                    CLOG(WARNING, "general_log") << "RaftGen failed. No space for raft.\n";

                    // Have to roll back the removed species, make sure RDEF reapplies them.
                    // Note: RDEF checks clamp
                    for (auto sg: solver::spec_global_id::range(lhs_s_vec.size())) {
                        uint lhs = lhs_s_vec[sg];
                        if (lhs == 0) {
                            continue;
                        }
                        uint prev_count = solverVesRaft()->getTriSpecCount_(tri->idx(), sg);
                        solverVesRaft()->setTriSpecCount_(tri->idx(), sg, prev_count + lhs);
                    }
                    continue;
                }

                // if we are here then the Patch successfully created the Raft
                // Set the spcies on the raft- they have already been removed from the
                // triangle
                for (auto sg: solver::spec_global_id::range(lhs_s_vec.size())) {
                    uint lhs = lhs_s_vec[sg];
                    if (lhs == 0) {
                        continue;
                    }
                    // Set the raft count
                    solverVesRaft()->setSingleRaftSpecCount_(rgidx, raft_unique_index, sg, lhs);
                }
            }
        }  // End of applying raftgens over tris

        // We have applied any raftgens in the triangle. Lastly, clear them.
        tri->clearAppliedRaftGens();
    }
    // _updateLocal is always called after runRaft by solver, so
    // no need to do any update here.
}

////////////////////////////////////////////////////////////////////////////////

}  // namespace steps::mpi::tetvesicle
