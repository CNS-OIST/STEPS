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

#include "mpi/tetvesicle/comp_vesraft.hpp"

#include <unordered_set>

#include "math/sphere.hpp"
#include "math/tetrahedron.hpp"
#include "mpi/tetvesicle/tet_vesraft.hpp"
#include "mpi/tetvesicle/tetvesicle_vesraft.hpp"
#include "mpi/tetvesicle/vesicle.hpp"
#include "solver/compdef.hpp"
#include "solver/exocytosisdef.hpp"
#include "solver/specdef.hpp"
#include "util/checkpointing.hpp"

namespace steps::mpi::tetvesicle {

CompVesRaft::CompVesRaft(solver::Compdef* compdef,
                         tetmesh::Tetmesh* mesh,
                         TetVesicleVesRaft* vesraft)
    : pCompdef(compdef)
    , pVol(0.0)
    , pMesh(mesh)
    , pRNG(vesraft->rng())  // have to do this because it's const
    , pVesRaft(vesraft) {
    AssertLog(pCompdef != nullptr);
    AssertLog(pMesh != nullptr);
    AssertLog(pVesRaft != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

CompVesRaft::~CompVesRaft() {
    for (auto& vit: pVesicles) {
        for (auto const& v: vit.second) {
            delete v;
        }
        vit.second.clear();
    }
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pVol);
    std::map<solver::vesicle_global_id, std::vector<solver::vesicle_individual_id>> vesicles;
    for (auto const& ves: pVesicles) {
        for (auto const& ves_p: ves.second) {
            vesicles[ves.first].emplace_back(ves_p->getUniqueIndex());
        }
    }
    util::checkpoint(cp_file, vesicles);

    for (auto const& ves: pVesicles) {
        for (auto const& ves_p: ves.second) {
            ves_p->checkpoint(cp_file);
        }
    }

    std::map<solver::vesicle_global_id, std::set<solver::comp_global_id>> ves_perm_comps;
    for (auto const& ves_comp: pVesicles_permittedcomps) {
        for (auto const& comp: ves_comp.second) {
            ves_perm_comps[ves_comp.first].insert(comp->def()->gidx());
        }
    }
    util::checkpoint(cp_file, ves_perm_comps);
    util::checkpoint(cp_file, pVes_Tetskcst);
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::restore(std::fstream& cp_file) {
    util::compare(cp_file, pVol);
    std::map<solver::vesicle_global_id, std::vector<solver::vesicle_individual_id>> vesicles;
    util::restore(cp_file, vesicles);

    for (auto const& ves: vesicles) {
        auto& vesdef = def()->statedef().vesicledef(ves.first);
        for (auto const& ves_id: ves.second) {
            auto* ves_p = new Vesicle(&vesdef, this, ves_id, cp_file);
            const auto& tets_overlap = ves_p->getOverlap_gidx();
            pVesRaft->addOverlap_(tets_overlap, ves_p);
            pVesicles[ves.first].push_back(ves_p);
        }
    }

    std::map<solver::vesicle_global_id, std::set<solver::comp_global_id>> ves_perm_comps;
    util::restore(cp_file, ves_perm_comps);
    for (auto const& ves_comp: ves_perm_comps) {
        for (auto const& cidx: ves_comp.second) {
            CompVesRaft* comp = solverVesRaft()->getComp_(cidx);
            pVesicles_permittedcomps[ves_comp.first].insert(comp);
        }
    }
    util::restore(cp_file, pVes_Tetskcst);
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::reset() {
    def()->reset();

    for (auto& vit: pVesicles) {
        for (auto const& v: vit.second) {
            delete v;
        }
        vit.second.clear();
    }

    for (auto& vtdit: pVes_Tetskcst) {
        vtdit.second.clear();
    }

    // Tet objects now clear their own overlap
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::setupVesicles() {
    // Create the vesicles, empty lists initially.
    for (auto const& v: def()->statedef().vesicles()) {
        solver::vesicle_global_id vesicledef_gidx = v->gidx();

        if (pVesicles.find(vesicledef_gidx) != pVesicles.end()) {
            std::ostringstream os;
            os << "Vesicle initialisation error. Index " << vesicledef_gidx << " is not unique.\n";
            ProgErrLog(os.str());
        }

        std::list<Vesicle*> ves;

        pVesicles[vesicledef_gidx] = ves;

        std::map<tetrahedron_global_id, double> tetskcst;
        pVes_Tetskcst[vesicledef_gidx] = tetskcst;

        pVesicles_permittedcomps[vesicledef_gidx].insert(this);
    }
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::addTet(TetVesRaft* tet) {
    AssertLog(tet->compdef() == def());

    tet->setCompVesRaft(this);

    index_t lidx_uint = pTets.size();
    tetrahedron_local_id lidx(lidx_uint);
    tetrahedron_global_id gidx(tet->idx());

    pTetidcs_L_to_G.emplace(lidx, gidx);
    pTetidcs_G_to_L.emplace(gidx, lidx);

    pTets.push_back(tet);

    pVol += tet->vol();
}

////////////////////////////////////////////////////////////////////////////////

solver::vesicle_individual_id CompVesRaft::addVesicle(solver::Vesicledef* vesdef,
                                                      const math::position_abs& pos,
                                                      tetrahedron_global_id tet_gidx) {
    solver::vesicle_global_id vgidx = vesdef->gidx();
    if (pVesicles.count(vgidx) == 0) {
        std::ostringstream os;
        os << "Vesicle addition error. Index " << vgidx << " is unknown in compartment.\n";
        ProgErrLog(os.str());
    }

    // The position check needs to come before the creation of the vesicle. This
    // is due to the behaviour of endocytosis- if endocytosis fails then it tries
    // again at a different location.

    if (tet_gidx.unknown()) {
        tet_gidx = mesh()->findTetByPoint(pos);
        if (tet_gidx.unknown()) {
            return {};
        }
        auto const search = pTetidcs_G_to_L.find(tet_gidx);
        if (search == pTetidcs_G_to_L.end()) {
            std::ostringstream os;
            os << "Central Position of vesicle is outside compartment.\n";
            ProgErrLog(os.str());
        }
    }

    // For overlap library
    vector_t vpos(pos[0], pos[1], pos[2]);

    std::map<tetrahedron_global_id, double> tets_overlap_new;

    Vesicle* jama_vesicle = nullptr;

    // Check if position is available- if not return undefined individual ID
    // We know central tet is in this compartment, so we don't need to
    bool avail =
        checkPos(&vpos, vesdef->diameter(), vgidx, tets_overlap_new, jama_vesicle, tet_gidx);
    if (avail == false) {
        return {};
    }

    solver::vesicle_individual_id unique_index = pVesRaft->getVesicleNextIndex_();

    // Create the vesicle
    auto* ves = new Vesicle(vesdef, this, pos, unique_index, tets_overlap_new);

    // addOverlap takes care of the references for tets
    pVesRaft->addOverlap_(tets_overlap_new, ves);

    pVesicles[vgidx].push_back(ves);

    // Check if this vesicle is on a path and set up if so
    checkVesiclePath(ves, pos);

    return unique_index;
}

////////////////////////////////////////////////////////////////////////////////

// There are different ways that this function needs to be used, i.e the vesicle
// may or may not already be in existence
bool CompVesRaft::checkPos(vector_t* pos,
                           double diam,
                           solver::vesicle_global_id vesgidx,
                           std::map<tetrahedron_global_id, double>& tets_overlap_output,
                           Vesicle*& jama_vesicle_output,
                           tetrahedron_global_id central_tet_idx,
                           Vesicle* ves,
                           math::point3d move_vector,
                           bool check_permcomps) {
    // First check if any link specs on ves would go beyond their limits on this
    // move vector, and disallow the move if so.
    if (ves != nullptr) {
        if (not ves->linkSpecMoveAllowed(move_vector)) {
            return false;
        }
    }

    // Now need to check we're inside the boundaries
    // The way to do that is to check (and store) tetrahedron overlap for the test
    // position. If total overlap of the vesicle is less than 100%  then it would
    // be outside the boundary -> don't move it here (return false) If total
    // overlap is 100% (with some small leway) then the position is good. the
    // checkPos function can be called with the vesicle already in existence
    // (argument ves is not null pointer) for which we have the (current)
    // tetrahedron overlap already stored, which is a good starting point for
    // search in the case of diffusion (where to always call mesh()->findtetByPoint
    // would presuably be too expensive). From functions that add a vesicle or
    // change it's position arbitratily in spce, the best start is the encompassing
    // tetrahedron

    // Search tets are going to be stored by breadth layer. The initial search
    // tets (from previous overlap or just encompassing tet of pos) are layer 0,
    // their neighbours are layer 1, etc. This will all be done with indexing (not
    // e.g. a map). The reason is that we need to search the entire vector each
    // time to see if a tet has already been visited

    // So (because we can't start with an empty set) we need either the central
    // tet or a vesicle.
    if (central_tet_idx.unknown() && ves == nullptr) {
        std::ostringstream os;
        os << "Either a vesicle or starting tet required for CompVesRaft::checkPos search.\n";
        ProgErrLog(os.str());
    }

    static std::unordered_set<tetrahedron_global_id> visited;
    std::vector<tetrahedron_global_id> curr_layer, next_layer;
    visited.clear();

    Sphere vesicle(*pos, diam / 2.0);

    double nearly_100 = 99.99999999;
    double total_overlap = 0.0;  // percentage i.e. 100 (or over nearly_100) means
                                 // vesicle completely overlaps mesh

    double vol_sphere = (4.0 / 3) * math::PI * pow(diam / 2.0, 3);

    math::position_abs pos2((*pos)[0], (*pos)[1], (*pos)[2]);

    // If we have the central tet, just use that. If not, use the previous overlap.
    if (central_tet_idx.valid()) {
        auto const tet_p = pVesRaft->tet_(central_tet_idx);
        // Position no good if outside all defined compartments
        if (tet_p == nullptr) {
            return false;
        }
        // If diffusive move, central tet must be in diffusion group
        if (check_permcomps and !checkVesiclePermittedComp(vesgidx, tet_p->getCompVesRaft())) {
            return false;
        }

        scalar_t ovlp = overlap(vesicle, *pVesRaft->tet_ext_(central_tet_idx));
        if (ovlp != 0.0) {
            total_overlap += (100 * ovlp) / vol_sphere;
            tets_overlap_output.emplace(central_tet_idx, ovlp);

            curr_layer.emplace_back(central_tet_idx);
            visited.emplace(central_tet_idx);
        } else {
            ProgErrLog("Central tetrahedron does not have any overlap with vesicle.");
        }
    } else {
        // Try to use previous overlap is vesicle positions are close enough
        if (ves->getPosition().dist2(pos2) < diam * diam) {
            for (const auto& itet: ves->getOverlap_gidx()) {
                scalar_t ovlp = overlap(vesicle, *pVesRaft->tet_ext_(itet.first));
                if (ovlp != 0.0) {
                    total_overlap += (100 * ovlp) / vol_sphere;
                    tets_overlap_output.emplace(itet.first, ovlp);

                    visited.emplace(itet.first);
                    curr_layer.emplace_back(itet.first);

                    if (total_overlap >= nearly_100) {
                        break;
                    }
                }
            }
        }
        // If no overlap was found, walk to the point
        if (curr_layer.empty()) {
            math::point3d pos3((*pos)[0], (*pos)[1], (*pos)[2]);
            // Use predicate to restrict the walk to permitted compartments
            auto walkable = [=](const tetrahedron_global_id& tet) {
                auto const tet_p = pVesRaft->tet_(tet);
                if (tet_p != nullptr and not tet_p->isFullOverlap()) {
                    // If the tet is in another compartment, potentially check if the vesicle can
                    // diffuse there.
                    return pTetidcs_G_to_L.find(tet) == pTetidcs_G_to_L.end() or
                           not check_permcomps or
                           checkVesiclePermittedComp(vesgidx, tet_p->getCompVesRaft());
                }
                return false;
            };
            // Do not allow getting further than maxSqDist^(1/2) from the target during the walk
            double maxSqDist = pVesRaft->getMaxWalkDistSqFact_() *
                               pMesh->_getTetBarycenter(ves->getCentralTet()).dist2(pos3);

            auto start = mesh()->findTetByPointWalk(
                pos2, ves->getCentralTet(), walkable, maxSqDist, pVesRaft->getMinNbTetVisited_());

            if (start.unknown()) {
                return false;
            } else {
                scalar_t ovlp = overlap(vesicle, *pVesRaft->tet_ext_(start));
                if (ovlp != 0.0) {
                    total_overlap += (100 * ovlp) / vol_sphere;
                    tets_overlap_output.emplace(start, ovlp);

                    visited.emplace(start);
                    curr_layer.emplace_back(start);
                }
            }
        }
    }

    // Will stop if no overlap in a given layer
    bool overlap_in_layer = true;

    // Grow the overlapping tetrahedrons layer by layer
    while (not curr_layer.empty() and overlap_in_layer and total_overlap < nearly_100) {
        overlap_in_layer = false;
        for (const auto& tet: curr_layer) {
            // Check neighbours
            for (auto const& neighb: mesh()->_getTetTetNeighb(tet)) {
                if (neighb.unknown()) {
                    continue;
                }
                // Check that neighb was not already visited
                if (visited.emplace(neighb).second) {
                    // Check that neighb can contain the vesicle
                    if (pTetidcs_G_to_L.find(neighb) == pTetidcs_G_to_L.end()) {
                        auto const tet_p = pVesRaft->tet_(neighb);
                        if (tet_p == nullptr) {
                            continue;
                        }
                        if (check_permcomps and
                            !checkVesiclePermittedComp(vesgidx, tet_p->getCompVesRaft())) {
                            continue;
                        }
                    }
                    // Check whether neighb overlaps the vesicle
                    scalar_t ovlp = overlap(vesicle, *pVesRaft->tet_ext_(neighb));
                    if (ovlp != 0.0) {
                        overlap_in_layer = true;
                        total_overlap += (100 * ovlp) / vol_sphere;
                        tets_overlap_output.emplace(neighb, ovlp);

                        // Only add neighb if it overlaps
                        next_layer.emplace_back(neighb);
                    }
                }
            }
        }
        std::swap(curr_layer, next_layer);
        next_layer.clear();
    }

    if (total_overlap < nearly_100) {
        return false;
    }

    // finally checking distance of all vesicles that overlap the same tets as
    // this one. Will include vesicles centred in other compartments too
    for (auto tetov: tets_overlap_output) {
        // tetov.first is global index of tetrahedron
        TetVesRaft* tet = pVesRaft->tet_(tetov.first);

        for (auto const& v: tet->getVesrefs()) {
            if (v == ves) {  // don't check against self
                continue;
            }

            double diam2 = v->getDiam();
            double distance2 = v->getPosition().dist2(pos2);

            if (sqrt(distance2) < (diam + diam2) / 2.0) {
                // We have a vesicle in the way. Allow one for a potential forced
                // swap. But if we have two that's a fail. Trick if we do have two,
                // then set jama_vesicle to nullptr and the potential force will fail.
                if (jama_vesicle_output != nullptr) {
                    // It's possible to reference the same vesicle multiple times
                    // since they can overlap multiple tets
                    if (jama_vesicle_output != v) {
                        jama_vesicle_output = nullptr;
                        return false;
                    } else {
                        continue;
                    }
                } else {
                    jama_vesicle_output = v;
                }
            }
        }
    }

    // Since we allow the previous loop to check for multiple overlap vesicles, we
    // might get here with one vesicle in the way. In that case we have to return
    // false because that allows the calling function to know we have to do a swap
    if (jama_vesicle_output != nullptr) {
        return false;
    } else {
        // Everything good. Position fine and unoccupied.
        return true;
    }
}

////////////////////////////////////////////////////////////////////////////////

tetrahedron_global_id CompVesRaft::tetidx_L_to_G(tetrahedron_local_id lidx) const {
    auto map_it = pTetidcs_L_to_G.find(lidx);
    AssertLog(map_it != pTetidcs_L_to_G.end());

    return map_it->second;
}

////////////////////////////////////////////////////////////////////////////////

tetrahedron_local_id CompVesRaft::tetidx_G_to_L(tetrahedron_global_id gidx) const {
    auto map_it = pTetidcs_G_to_L.find(gidx);
    AssertLog(map_it != pTetidcs_G_to_L.end());

    return map_it->second;
}

////////////////////////////////////////////////////////////////////////////////

tetrahedron_global_id CompVesRaft::getRandPosByTetStaticVols(math::position_abs* new_pos) const {
    // Not totally uniform, just do it by picking a tet at random
    double rn = rng()->getUnfII();

    // Have to find the volume this way because of the excluded volume effect

    double accum = 0.0;
    double selector = rn * pVol;

    // Select by volume
    tetrahedron_local_id tet_lidx(0);
    for (auto t: tets()) {
        accum += t->staticVol();
        if (selector < accum) {
            break;
        }
        tet_lidx += 1;
    }

    AssertLog(tet_lidx < countTets());

    tetrahedron_global_id tet_gidx(tetidx_L_to_G(tet_lidx));

    const std::vector<index_t>& tet_verts = mesh()->getTet(tet_gidx);

    // Need 3 random numbers
    double s = rng()->getUnfII();
    double t = rng()->getUnfII();
    double u = rng()->getUnfII();

    math::position_abs position = math::tet_ranpnt(mesh()->_getVertex(vertex_id_t(tet_verts[0])),
                                                   mesh()->_getVertex(vertex_id_t(tet_verts[1])),
                                                   mesh()->_getVertex(vertex_id_t(tet_verts[2])),
                                                   mesh()->_getVertex(vertex_id_t(tet_verts[3])),
                                                   s,
                                                   t,
                                                   u);

    *new_pos = position;

    return tet_gidx;
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::setVesicleCount(solver::vesicle_global_id vidx, uint count) {
    if (pVesicles.count(vidx) == 0) {
        std::ostringstream os;
        os << "Vesicle index " << vidx << " is unknown in compartment.\n";
        ProgErrLog(os.str());
    }

    // delete v will romove references through destructor
    for (auto const& v: pVesicles[vidx]) {
        delete v;
    }
    pVesicles[vidx].clear();

    auto& vesdef = def()->statedef().vesicledef(vidx);

    for (uint i = 0; i < count; ++i) {
        math::position_abs pos;
        solver::vesicle_individual_id added_vesicle;
        uint attempts = 0;
        // added_vesicle will return std::nullopt if unsuccesful
        while (added_vesicle.unknown()) {
            attempts++;
            if (attempts > 10000) {
                ArgErrLog("Unable to set count of vesicles: too many iterations.");
            }
            tetrahedron_global_id tet_gidx = getRandPosByTetStaticVols(&pos);
            added_vesicle = addVesicle(&vesdef, pos, tet_gidx);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::setVesicleTetDcst(solver::vesicle_global_id vidx,
                                    tetrahedron_global_id tidx,
                                    double dcst) {
    if (pVesicles.count(vidx) == 0) {
        std::ostringstream os;
        os << "Vesicle index " << vidx << " is unknown in compartment.\n";
        ProgErrLog(os.str());
    }

    if (dcst < 0.0) {
        std::ostringstream os;
        os << "Negative diffusion constant not allowed.\n";
        ProgErrLog(os.str());
    }

    // No checks on the tet idx- if it's wrong it'll just get ignored
    pVes_Tetskcst[vidx][tidx] = dcst;
}

////////////////////////////////////////////////////////////////////////////////

TetVesRaft* CompVesRaft::pickTetByVol(double rand01) const {
    if (countTets() == 0) {
        return nullptr;
    }
    if (countTets() == 1) {
        return pTets[0];
    }
    double volume = 0.0;
    for (const auto& tet: pTets) {
        volume += tet->vol();
    }
    double accum = 0.0;
    double selector = rand01 * volume;

    for (const auto& tet: pTets) {
        accum += tet->vol();
        if (selector < accum) {
            return tet;
        }
    }
    AssertLog(false);
    return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

uint CompVesRaft::getVesicleCount(solver::vesicle_global_id vidx) const {
    auto map_it = pVesicles.find(vidx);

    if (map_it == pVesicles.end()) {
        std::ostringstream os;
        os << "Vesicle index " << vidx << " is unknown in compartment.\n";
        ProgErrLog(os.str());
    }

    return map_it->second.size();
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::runVesicle(double dt) {
    // First, check for any exocytosis events and give them priority

    std::map<Vesicle*, solver::exocytosis_global_id, util::DerefPtrLess<Vesicle>> exo_vesicles;

    for (auto const& vit: pVesicles) {
        for (auto const& v: vit.second) {
            solver::exocytosis_global_id exo_gidx = v->appliedExocytosis();
            if (exo_gidx.valid()) {
                exo_vesicles.insert({v, exo_gidx});
            }
        }
    }

    for (auto const& ev: exo_vesicles) {
        // despite our best efforts its still possible an exocytosis was applied in a vesproxy when
        // linkspecs were created on another vesproxy
        if (ev.first->containsLink() == false) {
            // in applying Exocytosis, the pVesicle strucure will be altered
            applyExocytosis(ev.first, ev.second);
        } else {
            ev.first->clearExocytosis();
        }
    }

    for (auto const& vit: pVesicles) {
        // Need to copy because vesicles can get removed from pVesicles now during
        // run if lost to another compartment
        auto const vesicles = vit.second;
        for (auto const& v: vesicles) {
            if (v->getImmobility() == 0) {
                solver::vesicle_global_id ves_gidx = v->def()->gidx();

                tetrahedron_global_id centraltet = v->getCentralTet();

                double dcst = 0.0;

                double x_new_pos, y_new_pos, z_new_pos;
                math::point3d move_vector;

                math::position_abs v_pos = v->getPosition();

                if (v->onPath()) {
                    math::position_abs next_position = v->getNextPosition(dt);

                    x_new_pos = next_position[0];
                    y_new_pos = next_position[1];
                    z_new_pos = next_position[2];

                    double dx = x_new_pos - v_pos[0];
                    double dy = y_new_pos - v_pos[1];
                    double dz = z_new_pos - v_pos[2];

                    move_vector = math::point3d(dx, dy, dz);
                } else {
                    if (pVes_Tetskcst[ves_gidx].count(centraltet) > 0) {
                        dcst = pVes_Tetskcst[ves_gidx][centraltet];
                    } else {
                        dcst = v->getDcst();
                    }
                    float scale = sqrt(2 * dcst * dt);  // because rng (next 3 lines) returns floats
                    auto dx = static_cast<double>(rng()->getStdNrm() * scale);
                    auto dy = static_cast<double>(rng()->getStdNrm() * scale);
                    auto dz = static_cast<double>(rng()->getStdNrm() * scale);

                    move_vector = math::point3d(dx, dy, dz);

                    x_new_pos = v_pos[0] + dx;
                    y_new_pos = v_pos[1] + dy;
                    z_new_pos = v_pos[2] + dz;
                }

                vector_t new_pos(x_new_pos, y_new_pos, z_new_pos);

                std::map<tetrahedron_global_id, double> tets_overlap_prev = v->getOverlap_gidx();
                std::map<tetrahedron_global_id, double> tets_overlap_new;

                Vesicle* jama_vesicle = nullptr;  // not using it but need a reference

                // here
                if (checkPos(&new_pos,
                             v->getDiam(),
                             ves_gidx,
                             tets_overlap_new,
                             jama_vesicle,
                             std::nullopt,
                             v,
                             move_vector)) {
                    // Update position and overlap.
                    // v->setposition will return false if it doesn't move due to spec
                    // overlap issue
                    bool moved = v->setPosition(new_pos, tets_overlap_new);

                    if (moved) {
                        v->setOverlap(tets_overlap_new);

                        pVesRaft->removeOverlap_(tets_overlap_prev, v);
                        pVesRaft->addOverlap_(tets_overlap_new, v);

                        if (v->onPath()) {
                            v->nextPositionOnPath();
                        } else {
                            // Check if this position starts a new path for vesicle
                            math::position_abs pos_ves{new_pos[0], new_pos[1], new_pos[2]};
                            checkVesiclePath(v, pos_ves);
                        }
                    }
                }
                // Do nothing if position is bad, enforcing the Silver boundary
                // condition
            }

            // Finally perform surface diffusion, even if immobile
            v->doSurfaceDiffusion();
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

std::map<solver::spec_global_id, std::vector<PointSpec*>> CompVesRaft::removeOneVesicle(
    Vesicle* ves,
    tetrahedron_global_id tet_gidx) {
    std::map<solver::spec_global_id, std::vector<PointSpec*>> ves_specs;

    if (ves->overlapsTet_gidx(tet_gidx)) {
        // Get species from vesicle to return that info if e.g exocytosis needs it
        ves_specs = ves->getSurfSpecs();

        solver::vesicle_global_id ves_gidx = ves->def()->gidx();

        // Finally remove the vesicle from list and delete it.
        pVesicles[ves_gidx].remove(ves);

        delete (ves);  // Vesicle destructor will sort out overlap and vesicle
                       // references in tetrahedrons

    } else {
        std::ostringstream os;
        os << "Failed to remove vesicle from tet with global index " << tet_gidx
           << ": no overlap.\n";
        ProgErrLog(os.str());
    }

    return ves_specs;
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::applyExocytosis(Vesicle* vesicle, solver::exocytosis_global_id exo_gidx) {
    AssertLog(exo_gidx < def()->statedef().countExocytosis());

    auto& exodef = def()->statedef().exocytosisdef(exo_gidx);
    auto vidx = vesicle->getUniqueIndex();

    // For chosing a triangle for application to
    std::vector<TriVesRaft*> poss_tris;
    // Need to store the tets too for the update
    std::vector<TetVesRaft*> conn_tets;

    const std::vector<tetrahedron_global_id>& tets_gidx = vesicle->getOverlapVec_gidx();

    for (auto const& tetgidx: tets_gidx) {
        TetVesRaft* tet = pVesRaft->tet_(tetgidx);

        for (auto const& tetnexttri: tet->nexttris()) {
            if (tetnexttri == nullptr) {
                continue;
            }
            poss_tris.push_back(tetnexttri);
            conn_tets.push_back(tet);
        }
    }

    uint tupdate_size = poss_tris.size();
    // Should have gotten at least one triangle if rate is non-zero, or only +ve
    // update
    AssertLog(tupdate_size > 0);

    // Here need to select a triangle by raft availablity. If none available then
    // Exocytosis fails.
    uint rIndex = std::numeric_limits<uint>::max();
    solver::Raftdef* raftdef = exodef.raftdef();
    solver::raft_individual_id raft_unique_index;

    if (raftdef != nullptr) {
        // Try all possible triangles for maximum chance exocytosis will be
        // successful
        for (uint tidx = 0; tidx < tupdate_size; ++tidx) {
            TriVesRaft* tri = poss_tris[tidx];
            PatchVesRaft* patch = tri->patchVesRaft();
            solver::raft_individual_id added_raft = patch->addRaft(raftdef, tri);
            if (added_raft.valid()) {
                rIndex = tidx;
                raft_unique_index = added_raft;
                break;
            }
        }

        if (rIndex == std::numeric_limits<uint>::max()) {
            // If we got here it means we could not create the raft.
            CLOG(WARNING, "general_log") << "Exocytosis failed. No space for raft.\n";
            // UpdVec is empty anyway so can just return that.
            return;
        }

        AssertLog(raft_unique_index.valid());
        AssertLog(rIndex != std::numeric_limits<uint>::max());
    }
    // Select randomly from vector and update
    else {
        rIndex = rng()->get() % tupdate_size;
    }

    TriVesRaft* update_tri = poss_tris[rIndex];
    TetVesRaft* update_tet = conn_tets[rIndex];

    tetrahedron_global_id tet_gidx = update_tet->idx();

    // Got to release the contents into outside compartment
    // First see if there is actually a compartment there or not
    TetVesRaft* itet = update_tri->iTet();
    TetVesRaft* otet = update_tri->oTet();
    TetVesRaft* extra_tet;
    if (itet != update_tet) {
        AssertLog(otet == update_tet);
        extra_tet = itet;
    } else {
        extra_tet = otet;
    }

    std::map<solver::spec_global_id, uint> kiss_and_run_inner_spec_loss;

    for (auto& inner_spec_count: vesicle->getInnerSpecCounts()) {
        solver::spec_global_id spec_gidx = inner_spec_count.first;
        solver::spec_local_id spec_lidx;
        if (extra_tet != nullptr) {
            spec_lidx = extra_tet->compdef()->specG2L(spec_gidx);
            if (spec_lidx.unknown()) {
                CLOG(WARNING, "general_log")
                    << "Applying Exocytosis '" << exodef.name() << "': Species '"
                    << solverVesRaft()->statedef().specdef(spec_gidx).name()
                    << "' undefined in Comp '" << extra_tet->compdef()->name()
                    << "'. Species will not be added to Comp.\n";
            }
        }

        uint inner_spec_count_loss = inner_spec_count.second;

        if (exodef.getKissAndRun()) {
            double n = inner_spec_count_loss * exodef.getKissAndRunPartRelease();
            double n_int = std::floor(n);
            double n_frc = n - n_int;
            inner_spec_count_loss = static_cast<uint>(n_int);
            if (n_frc > 0.0) {
                double rand01 = rng()->getUnfIE();
                if (rand01 < n_frc) {
                    inner_spec_count_loss++;
                }
            }
            kiss_and_run_inner_spec_loss[spec_gidx] = inner_spec_count_loss;
        }

        if (extra_tet != nullptr and spec_lidx.valid()) {
            pVesRaft->regTetPoolSync_(extra_tet->idx(), spec_gidx, inner_spec_count_loss);
        }
    }

    if (exodef.getKissAndRun()) {
        // In this version of kiss-and-run, some species (optionally) change ID.
        // Immobility of vesicle remains unaffected

        for (auto const& spec_gidxs: exodef.getKissAndRunSpecChanges()) {
            // Vesicle does the heavy lifting
            vesicle->changeSurfSpecGidx(spec_gidxs.first, spec_gidxs.second);
        }

        vesicle->clearExocytosis();
        vesicle->redInnerSpecCounts(kiss_and_run_inner_spec_loss);
    } else {
        std::map<solver::spec_global_id, std::vector<PointSpec*>> ves_specs =
            removeOneVesicle(vesicle, tet_gidx);

        solver::Patchdef* pdef = update_tri->patchdef();


        if (raftdef != nullptr) {
            solver::raft_global_id rgidx = raftdef->gidx();
            // This will have elements of the else part, but safest to keep them
            // separate I think
            for (auto const& s: ves_specs) {
                solver::spec_global_id spec_gidx = s.first;
                // We're creating a raft here so no need to check if this species already
                // has some count
                solverVesRaft()->setSingleRaftSpecCount_(rgidx,
                                                         raft_unique_index,
                                                         spec_gidx,
                                                         s.second.size());
            }
        } else {
            for (auto const& s: ves_specs) {
                solver::spec_global_id spec_gidx = s.first;
                solver::spec_local_id spec_lidx = pdef->specG2L(spec_gidx);

                if (spec_lidx.unknown()) {
                    std::ostringstream os;
                    CLOG(WARNING, "general_log")
                        << "Applying Exocytosis '" << exodef.name() << "': Species '"
                        << solverVesRaft()->statedef().specdef(spec_gidx).name()
                        << "' undefined in Patch '" << pdef->name()
                        << "'. Species will be not be added to Patch.\n";
                    continue;
                }

                // Note: clamped status doesn't matter because RDEF stores that and will check clanp
                // during sync
                uint prev_count = pVesRaft->getTriSpecCount_(update_tri->idx(), spec_gidx);
                pVesRaft->setTriSpecCount_(update_tri->idx(),
                                           spec_gidx,
                                           prev_count + s.second.size());
            }
        }
    }

    exodef.addEvent(def()->statedef().time(), vidx, update_tri->idx(), raft_unique_index);
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::loseVesicle(Vesicle* ves) {
    solver::vesicle_global_id ves_gidx = ves->def()->gidx();

    // Finally remove the vesicle from list and delete it.
    pVesicles[ves_gidx].remove(ves);
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::gainVesicle(Vesicle* ves) {
    solver::vesicle_global_id ves_gidx = ves->def()->gidx();

    // We're just gaining a vesicle from another compartment. All checks
    // on position should have passed already.
    pVesicles[ves_gidx].push_back(ves);
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::deleteSingleVesicle(Vesicle* ves) {
    if (ves->containsLink()) {
        CLOG(WARNING, "general_log") << "Vesicle unique index " << ves->getUniqueIndex()
                                     << " contains link species. Not deleted.\n";
        return;
    }

    pVesicles[ves->idx()].remove(ves);
    delete (ves);  // Vesicle destructor will sort out overlap and vesicle
                   // references in tetrahedrons
    return;
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::setVesiclePos(solver::vesicle_global_id vidx,
                                solver::vesicle_individual_id ves_unique_idx,
                                const std::vector<double>& pos,
                                bool force) {
    if (pVesicles.count(vidx) == 0) {
        std::ostringstream os;
        os << "Vesicle addition error. Index " << vidx << " is unknown in compartment.\n";
        ProgErrLog(os.str());
    }
    for (auto const& v: pVesicles[vidx]) {
        if (v->getUniqueIndex() == ves_unique_idx) {
            AssertLog(v->idx() == vidx);

            auto tet_gidx = mesh()->findTetByPoint(pos);
            if (tet_gidx.unknown()) {
                CLOG(WARNING, "general_log")
                    << "Did not change position: outside mesh boundaries.\n";
                return;
            }

            vector_t new_pos(pos[0], pos[1], pos[2]);

            math::position_abs old_pos = v->getPosition();

            vector_t old_pos_vect(old_pos[0], old_pos[1], old_pos[2]);

            math::point3d move_vector =
                math::point3d(pos[0] - old_pos[0], pos[1] - old_pos[1], pos[2] - old_pos[2]);

            std::map<tetrahedron_global_id, double> tets_overlap_prev = v->getOverlap_gidx();
            std::map<tetrahedron_global_id, double> tets_overlap_new;

            Vesicle* jama_vesicle = nullptr;

            if (checkPos(&new_pos,
                         v->getDiam(),
                         vidx,
                         tets_overlap_new,
                         jama_vesicle,
                         tet_gidx,
                         v,
                         move_vector,
                         false)) {
                bool moved = v->setPosition(new_pos, tets_overlap_new, true);
                if (moved) {
                    v->setOverlap(tets_overlap_new);
                    pVesRaft->removeOverlap_(tets_overlap_prev, v);
                    pVesRaft->addOverlap_(tets_overlap_new, v);

                    // Now need to remove from any paths it's currently on
                    v->removeFromPath();

                    // check if this position starts a new path for vesicle
                    math::position_abs pos_ves{new_pos[0], new_pos[1], new_pos[2]};

                    checkVesiclePath(v, pos_ves);
                } else {
                    CLOG(WARNING, "general_log") << "Did not change position: surface "
                                                 << "molecule has no overlap tet.\n";
                }
            } else {
                if ((!force) || (jama_vesicle == nullptr)) {
                    // Position is not available. If not force do nothing, or if force but
                    // we have no vesicle it means position is outside geometrical
                    // boundaries.
                    CLOG(WARNING, "general_log")
                        << "Did not change position: already occupied or outside"
                        << "permitted compartment boundaries.\n";

                } else {
                    // We want to force and we have a vesicle in the way. Swap.
                    bool moved1 = v->setPosition(new_pos, tets_overlap_new, true);
                    AssertLog(moved1);  // positions should be valid
                    v->setOverlap(tets_overlap_new);
                    pVesRaft->removeOverlap_(tets_overlap_prev, v);
                    pVesRaft->addOverlap_(tets_overlap_new, v);
                    v->removeFromPath();
                    math::position_abs pos_ves{new_pos[0], new_pos[1], new_pos[2]};

                    checkVesiclePath(v, pos_ves);

                    std::map<tetrahedron_global_id, double> tets_overlap_jv =
                        jama_vesicle->getOverlap_gidx();

                    bool moved2 = jama_vesicle->setPosition(old_pos_vect, tets_overlap_prev, true);
                    AssertLog(moved2);  // positions should be valid
                    jama_vesicle->setOverlap(tets_overlap_prev);
                    pVesRaft->removeOverlap_(tets_overlap_jv, jama_vesicle);
                    pVesRaft->addOverlap_(tets_overlap_prev, jama_vesicle);
                    jama_vesicle->removeFromPath();
                    math::position_abs pos_jamaves{old_pos[0], old_pos[1], old_pos[2]};
                    checkVesiclePath(jama_vesicle, pos_jamaves);
                }
            }
            // All roads lead here- safe to have only this return statement
            return;
        }
    }

    CLOG(WARNING, "general_log") << "vesicle unique index " << ves_unique_idx
                                 << " is not in use, nothing changed.";
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::checkVesiclePath(Vesicle* v, math::position_abs const& pos_ves) const {
    solver::vesicle_global_id ves_gidx = v->idx();

    double ves_radius = v->getDiam() / 2.0;

    auto ves_paths = pVesRaft->vesicleCrossedPaths_(pos_ves, ves_gidx, ves_radius);

    for (auto const& ves_path: ves_paths) {
        // One final check to see if the spec deps are there
        std::map<solver::spec_global_id, uint> const& path_ves_speceps =
            ves_path->getVesicleSpecDeps(ves_gidx);

        bool gotdeps = true;
        for (auto const& vsd: path_ves_speceps) {
            if (v->getSurfSpecCount(vsd.first) < vsd.second) {
                gotdeps = false;
                break;
            }
        }

        if (gotdeps) {
            auto route = ves_path->calculateRoute(pos_ves, ves_gidx, rng(), ves_radius);
            v->setPathPositions(route);
            return;  // vesicle can only be on one path at a time
        }
        // If the deps are not there we simply continue to the next path if there is one
    }
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::addVesicleSpecs(solver::vesicle_global_id vidx,
                                  solver::vesicle_individual_id ves_unique_idx,
                                  std::map<solver::spec_global_id, int> ves_specs) {
    if (pVesicles.count(vidx) == 0) {
        std::ostringstream os;
        os << "Vesicle addition error. Index " << vidx << " is unknown in compartment.\n";
        ProgErrLog(os.str());
    }
    for (auto const& ves: pVesicles[vidx]) {
        if (ves->getUniqueIndex() == ves_unique_idx) {
            ves->addSurfSpecs(ves_specs);
            return;
        }
    }

    std::ostringstream os;
    os << "Currently, vesicles of this type do not contain unique index " << ves_unique_idx
       << ".\n";
    ProgErrLog(os.str());

    return;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::vesicle_individual_id> CompVesRaft::getVesicleIndices(
    solver::vesicle_global_id vidx) const {
    std::vector<solver::vesicle_individual_id> ves_indices;

    auto vit = pVesicles.find(vidx);
    if (vit != pVesicles.end()) {
        for (auto const& v: vit->second) {
            ves_indices.emplace_back(v->getUniqueIndex());
        }
    } else {
        std::ostringstream os;
        os << "Vesicle index " << vidx << " is unknown in compartment.\n";
        ProgErrLog(os.str());
    }

    return ves_indices;
}

////////////////////////////////////////////////////////////////////////////////

std::map<solver::vesicle_individual_id, uint> CompVesRaft::getVesicleSurfaceSpecCountMap(
    solver::vesicle_global_id vidx,
    solver::spec_global_id spec_gidx) const {
    auto map_it = pVesicles.find(vidx);

    if (map_it == pVesicles.end()) {
        std::ostringstream os;
        os << "Vesicle addition error. Index " << vidx << " is unknown in compartment.\n";
        ProgErrLog(os.str());
    }

    std::map<solver::vesicle_individual_id, uint> spec_counts;

    for (auto const& ves: map_it->second) {
        spec_counts[ves->getUniqueIndex()] = ves->getSurfSpecCount(spec_gidx);
    }

    return spec_counts;
}

////////////////////////////////////////////////////////////////////////////////

uint CompVesRaft::getVesicleSurfaceSpecCount(solver::vesicle_global_id vidx,
                                             solver::spec_global_id spec_gidx) const {
    auto map_it = pVesicles.find(vidx);

    if (map_it == pVesicles.end()) {
        std::ostringstream os;
        os << "Vesicle addition error. Index " << vidx << " is unknown in compartment.\n";
        ProgErrLog(os.str());
    }

    uint spec_counts = 0;

    for (auto const& ves: map_it->second) {
        spec_counts += ves->getSurfSpecCount(spec_gidx);
    }

    return spec_counts;
}

////////////////////////////////////////////////////////////////////////////////

uint CompVesRaft::getVesicleInnerSpecCount(solver::vesicle_global_id vidx,
                                           solver::spec_global_id spec_gidx) const {
    auto map_it = pVesicles.find(vidx);

    if (map_it == pVesicles.end()) {
        std::ostringstream os;
        os << "Vesicle addition error. Index " << vidx << " is unknown in compartment.\n";
        ProgErrLog(os.str());
    }

    uint spec_counts = 0;

    for (auto const& ves: map_it->second) {
        spec_counts += ves->getInnerSpecCount(spec_gidx);
    }

    return spec_counts;
}

////////////////////////////////////////////////////////////////////////////////

uint CompVesRaft::getVesicleLinkSpecCount(solver::vesicle_global_id vidx,
                                          solver::linkspec_global_id linkspec_gidx) const {
    auto map_it = pVesicles.find(vidx);

    if (map_it == pVesicles.end()) {
        std::ostringstream os;
        os << "Vesicle addition error. Index " << vidx << " is unknown in compartment.\n";
        ProgErrLog(os.str());
    }

    uint lspec_counts = 0;

    for (auto const& ves: map_it->second) {
        lspec_counts += ves->getLinkSpecCount(linkspec_gidx);
    }

    return lspec_counts;
}

////////////////////////////////////////////////////////////////////////////////

std::map<solver::vesicle_individual_id, uint> CompVesRaft::getVesicleLinkSpecCountMap(
    solver::vesicle_global_id vidx,
    solver::linkspec_global_id linkspec_gidx) const {
    auto map_it = pVesicles.find(vidx);

    if (map_it == pVesicles.end()) {
        std::ostringstream os;
        os << "Vesicle addition error. Index " << vidx << " is unknown in compartment.\n";
        ProgErrLog(os.str());
    }

    std::map<solver::vesicle_individual_id, uint> lspec_counts;

    for (auto const& ves: map_it->second) {
        lspec_counts[ves->getUniqueIndex()] = ves->getLinkSpecCount(linkspec_gidx);
    }

    return lspec_counts;
}

////////////////////////////////////////////////////////////////////////////////

bool CompVesRaft::getVesicleTetOverlap(solver::vesicle_global_id vidx,
                                       tetrahedron_global_id tet_gidx) const {
    auto map_it = pVesicles.find(vidx);

    if (map_it == pVesicles.end()) {
        std::ostringstream os;
        os << "Vesicle index " << vidx << " is unknown in compartment.\n";
        ProgErrLog(os.str());
    }

    for (auto const& ves: map_it->second) {
        if (ves->overlapsTet_gidx(tet_gidx)) {
            return true;
        }
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////

void CompVesRaft::addVesiclePermittedComps(solver::vesicle_global_id vidx,
                                           const std::vector<CompVesRaft*>& comps) {
    auto compsIt = pVesicles_permittedcomps.find(vidx);
    if (compsIt == pVesicles_permittedcomps.end()) {
        std::ostringstream os;
        os << "Vesicle index " << vidx << " is unknown in compartment " << pCompdef->name();
        ProgErrLog(os.str());
    }

    compsIt->second.insert(comps.begin(), comps.end());
}

}  // namespace steps::mpi::tetvesicle
