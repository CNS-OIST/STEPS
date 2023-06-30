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

#include "mpi/tetvesicle/vesicle.hpp"

// STEPS headers.
#include "math/constants.hpp"
#include "math/point.hpp"
#include "math/sphere.hpp"
#include "mpi/mpi_common.hpp"
#include "mpi/tetvesicle/comp_vesraft.hpp"
#include "mpi/tetvesicle/linkspec.hpp"
#include "mpi/tetvesicle/tetvesicle_vesraft.hpp"
#include "solver/vesicledef.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"

namespace steps::mpi::tetvesicle {

Vesicle::Vesicle(solver::Vesicledef* vesdef,
                 CompVesRaft* comp,
                 const math::position_abs& pos,
                 solver::vesicle_individual_id unique_index,
                 const std::map<tetrahedron_global_id, double>& overlap)
    : pDef(vesdef)
    , pComp_central(comp)
    , pIndex(unique_index)
    , pImmobility(0)
    , pPathPosition_index(0)
    , pOnPath(false) {
    AssertLog(pDef != nullptr);
    AssertLog(comp != nullptr);

    pComp_central->solverVesRaft()->recordVesicle_(pIndex, this);

    pPos = pos;

    setOverlap(overlap);
}

////////////////////////////////////////////////////////////////////////////////

Vesicle::Vesicle(solver::Vesicledef* vesdef,
                 CompVesRaft* comp,
                 solver::vesicle_individual_id unique_index,
                 std::fstream& cp_file)
    : pDef(vesdef)
    , pComp_central(comp)
    , pIndex(unique_index) {
    AssertLog(pDef != nullptr);
    AssertLog(comp != nullptr);

    pComp_central->solverVesRaft()->recordVesicle_(pIndex, this);

    util::compare(cp_file, pIndex);
    util::restore(cp_file, pPos);
    util::restore(cp_file, pCentral_tet);
    util::restore(cp_file, pTets_overlap_gidx);
    util::restore(cp_file, pTets_overlap_vec_gidx);
    util::restore(cp_file, pInnerSpecCount);
    util::restore(cp_file, pImmobility);
    util::restore(cp_file, pPathPositions);
    util::restore(cp_file, pPathPosition_index);
    util::restore(cp_file, pPathPosition_index_next);
    util::restore(cp_file, pTime_accum);
    util::restore(cp_file, pTime_accum_next);
    util::restore(cp_file, pOnPath);

    std::map<solver::spec_global_id, uint> surfspecs;
    util::restore(cp_file, surfspecs);
    for (auto const& ss: surfspecs) {
        for (uint i = 0; i < ss.second; ++i) {
            // Constructor does the restore
            auto* ps = new PointSpec(ss.first, getDiam() / 2.0, cp_file);
            pSurfSpecs[ss.first].emplace_back(ps);
        }
    }

    std::vector<solver::linkspec_global_id> linkspecs;
    util::restore(cp_file, linkspecs);
    for (auto const& ls_id: linkspecs) {
        solver::LinkSpecdef* lspecdef = def()->statedef()->linkspecdef(ls_id);
        // Constructor does the restore
        auto ls = new LinkSpec(lspecdef, this, cp_file);
        pLinkSpecs[ls->getUniqueID()] = ls;
    }
}

////////////////////////////////////////////////////////////////////////////////

Vesicle::~Vesicle() {
    for (auto& sit: pSurfSpecs) {
        for (auto& s: sit.second) {
            delete s;
        }
        sit.second.clear();
    }

    // Need to remove overlap from compartment
    for (auto const& vov: getOverlap_gidx()) {
        solver()->tet_(vov.first)->changeOverlap(-(vov.second));
        solver()->tet_(vov.first)->removeVesref(this);
    }

    // During sim a vesicle with linkspecs should not be destroyed. HOWEVER during
    // solver destruction obviously a vesicle with linkspecs can be destroyed,
    // in which case we have to delete here to avoid memory leaks.
    for (auto ls: pLinkSpecs) {
        delete (ls.second);
    }
    pLinkSpecs.clear();

    comp()->solverVesRaft()->removeVesicle_(getUniqueIndex(), this);
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pIndex);
    util::checkpoint(cp_file, pPos);  // pPos will not be available at time of restore
    util::checkpoint(cp_file, pCentral_tet);
    util::checkpoint(cp_file, pTets_overlap_gidx);
    util::checkpoint(cp_file, pTets_overlap_vec_gidx);
    util::checkpoint(cp_file, pInnerSpecCount);
    util::checkpoint(cp_file, pImmobility);
    util::checkpoint(cp_file, pPathPositions);
    util::checkpoint(cp_file, pPathPosition_index);
    util::checkpoint(cp_file, pPathPosition_index_next);
    util::checkpoint(cp_file, pTime_accum);
    util::checkpoint(cp_file, pTime_accum_next);
    util::checkpoint(cp_file, pOnPath);
    // pAppliedExocytosis should be empty at point of checkpointing (cleared by clearExocytosis())
    AssertLog(pAppliedExocytosis.empty());

    std::map<solver::spec_global_id, uint> surfspecs;
    for (auto const& ss: pSurfSpecs) {
        surfspecs[ss.first] = ss.second.size();
    }
    util::checkpoint(cp_file, surfspecs);
    for (auto const& ss: pSurfSpecs) {
        for (auto const& p: ss.second) {
            p->checkpoint(cp_file);
        }
    }

    std::vector<solver::linkspec_global_id> linkspecs;
    for (auto const& ls: pLinkSpecs) {
        linkspecs.emplace_back(ls.second->getGidx());
    }
    util::checkpoint(cp_file, linkspecs);
    for (auto const& ls: pLinkSpecs) {
        ls.second->checkpoint(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

TetVesicleVesRaft* Vesicle::solver() const noexcept {
    return comp()->solverVesRaft();
}

////////////////////////////////////////////////////////////////////////////////

const rng::RNGptr Vesicle::rng() const noexcept {
    return comp()->rng();
}

////////////////////////////////////////////////////////////////////////////////

tetmesh::Tetmesh* Vesicle::mesh() const noexcept {
    return comp()->mesh();
}

////////////////////////////////////////////////////////////////////////////////

bool Vesicle::setPosition(const vector_t& new_pos,
                          std::map<tetrahedron_global_id, double> const& tets_overlap_temp) {
    // Had numerical issues with this. Setting new position, even with 100%
    // overlap, was sometimes unable to find tet for a species. Now do a check and
    // if that happens, just don't move

    // So currently does loop twice: once to check, once to set

    std::map<PointSpec*, tetrahedron_global_id> ps_updtets;

    for (auto const& sit: pSurfSpecs) {
        for (auto const& psit: sit.second) {
            math::position_rel_to_ves pos_rel = psit->getPosCartesian();
            math::position_abs pos_abs{
                pos_rel[0] + new_pos[0],
                pos_rel[1] + new_pos[1],
                pos_rel[2] + new_pos[2],
            };

            tetrahedron_global_id tetgidx = getTetSpecOverlap(pos_abs, tets_overlap_temp);
            if (tetgidx.unknown()) {
                return false;
            }
            ps_updtets[psit] = tetgidx;
        }
    }

    for (auto const& ps: ps_updtets) {
        ps.first->setOverlapTet_gidx(ps.second);
    }

    // All good in neighbourhood

    pPos[0] = new_pos[0];
    pPos[1] = new_pos[1];
    pPos[2] = new_pos[2];

    tetrahedron_global_id tetov_gidx = mesh()->findTetByPoint(pPos, tets_overlap_temp);

    AssertLog(tetov_gidx.valid());

    pCentral_tet = tetov_gidx;

    // Need to check central compartment. If that has changed then the comp needs
    // to know about it.
    CompVesRaft* central_comp = solver()->tet_(pCentral_tet)->getCompVesRaft();

    if (central_comp != pComp_central) {
        // Compartment has changed
        pComp_central->loseVesicle(this);
        central_comp->gainVesicle(this);

        pComp_central = central_comp;

        return true;
    }

    else {
        return true;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::updImmobility(int mob_upd) {
    if (pImmobility == 0 and mob_upd < 0) {
        std::ostringstream os;
        os << "Negative immobility is not possible for vesicle. Model error. ";
        ProgErrLog(os.str());
    }

    pImmobility += mob_upd;
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::setOverlap(const std::map<tetrahedron_global_id, double>& overlap) {
    pTets_overlap_gidx = overlap;

    pTets_overlap_vec_gidx.clear();

    for (auto const& tets: pTets_overlap_gidx) {
        pTets_overlap_vec_gidx.emplace_back(tets.first);
    }

    pCentral_tet = mesh()->findTetByPoint(pPos, overlap);
}

////////////////////////////////////////////////////////////////////////////////

tetrahedron_global_id Vesicle::getTetSpecOverlap(
    const math::position_abs& pos,
    std::map<tetrahedron_global_id, double> const& tets_overlap_temp) {
    tetrahedron_global_id tetov_gidx = mesh()->findTetByPoint(pos, tets_overlap_temp);

    // May be unknown, but will be checked
    return tetov_gidx;
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::doSurfaceDiffusion() {
    solver::vesicle_global_id ves_gidx = pDef->gidx();

    for (auto const& sit: pSurfSpecs) {
        for (auto const& psit: sit.second) {
            double phi = solver()->getQPhiSpec_(ves_gidx, sit.first);

            // phi may be zero if species not defined to have surface diffusion rate
            if (phi != 0.0) {
                uint attempts = 0;
                tetrahedron_global_id tetgidx;

                // record the original cartesians in case we need to return to them:
                math::position_rel_to_ves pos_rel_orig = psit->getPosCartesian();

                while (tetgidx.unknown()) {
                    attempts += 1;
                    if (attempts == 1000) {
                        std::ostringstream os;
                        os << "Failed to do surface diffusion; too many failed attempts. ";
                        ProgErrLog(os.str());
                    }

                    // Need to pick a random direction, theta, in range 0-2pi
                    double theta = rng()->getUnfII() * math::PI * 2.0;

                    psit->updatePos(theta, phi);

                    math::position_rel_to_ves pos_rel = psit->getPosCartesian();

                    // Those are relative positions- make them absolute
                    math::position_abs pos_abs{
                        pos_rel[0] + pPos[0],
                        pos_rel[1] + pPos[1],
                        pos_rel[2] + pPos[2],
                    };

                    tetgidx = getTetSpecOverlap(pos_abs, pTets_overlap_gidx);

                    // If we don't have tet we're going to try again so have to reset
                    // position
                    if (tetgidx.unknown()) {
                        psit->setPosCartesian(pos_rel_orig);
                    }
                }

                psit->setOverlapTet_gidx(tetgidx);
            }
        }
    }


    for (auto const& link_spec: pLinkSpecs) {
        auto ls = link_spec.second;
        // Here need to do above but check if length is still within bounds first,
        // if not reject
        double phi = solver()->getQPhiLinkspec_(ves_gidx, ls->getGidx());


        // phi may be zero if species not defined to have surface diffusion rate
        if (phi != 0.0) {
            // record the original cartesians in case we need to return to them:
            math::position_rel_to_ves pos_rel_orig = ls->getPosCartesian_rel();

            // Need to pick a random direction, theta, in range 0-2pi
            double theta = rng()->getUnfII() * math::PI * 2.0;

            ls->updatePos(theta, phi);

            if (not ls->withinBounds()) {
                ls->setPosCartesian_rel(pos_rel_orig);
                continue;
            }

            math::position_rel_to_ves pos_rel = ls->getPosCartesian_rel();

            // Those are relative positions- make them absolute
            math::position_abs pos_abs{
                pos_rel[0] + pPos[0],
                pos_rel[1] + pPos[1],
                pos_rel[2] + pPos[2],
            };
            tetrahedron_global_id tetgidx = getTetSpecOverlap(pos_abs, pTets_overlap_gidx);

            // If we don't have tet we're not going ahead so have to reset position
            if (tetgidx.unknown()) {
                ls->setPosCartesian_rel(pos_rel_orig);
                continue;
            }
            // we can go ahead with the new coordinates. Just need to update the
            // overlap tet
            ls->setOverlapTet_gidx(tetgidx);  // Still necessary??
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::setSurfSpecCount(solver::spec_global_id sgidx, uint count) {
    uint sscount = getSurfSpecCount(sgidx);

    // Do nothing if the count to set it to is the same as the current count
    if (sscount == count) {
        return;
    }

    uint starting_count = 0;
    if (pSurfSpecs.count(sgidx) != 0) {
        starting_count = pSurfSpecs[sgidx].size();
    } else {
        pSurfSpecs[sgidx] = std::vector<PointSpec*>();
    }

    if (count > starting_count) {
        // Add point specs one by one and randomise position
        for (uint i = starting_count; i < count; ++i) {
            auto* ps = new PointSpec(sgidx, getDiam() / 2.0, solver()->getPointSpecNextIndex_());
            math::position_rel_to_ves pos_rel;

            pSurfSpecs[sgidx].push_back(ps);

            // Positions are relative to the centre of the vesicle
            // Initial position is randomised

            tetrahedron_global_id tetgidx;
            while (tetgidx.unknown()) {
                pos_rel = (getDiam() / 2.0) * math::sphere_unit_randsurfpos(rng());

                ps->setPosCartesian(pos_rel);

                math::position_abs pos_abs{
                    pos_rel[0] + pPos[0],
                    pos_rel[1] + pPos[1],
                    pos_rel[2] + pPos[2],
                };
                // pointSpecs now have to be told in which tetrahedron they currently
                // reside
                tetgidx = getTetSpecOverlap(pos_abs, pTets_overlap_gidx);
            }

            ps->setOverlapTet_gidx(tetgidx);
        }
    } else if (starting_count > count) {
        // Removing from the back. TODO choose randomly?
        for (uint i = starting_count; i > count; --i) {
            delete pSurfSpecs[sgidx][i - 1];
            pSurfSpecs[sgidx].pop_back();
        }
    } else {
        // Shouldn't get here because we checked for equality
        AssertLog(false);
    }
}

////////////////////////////////////////////////////////////////////////////////

uint Vesicle::getSurfSpecCount(solver::spec_global_id spec_gidx, tetrahedron_global_id tet_gidx) {
    if (pSurfSpecs.count(spec_gidx) == 0) {
        return 0;
    } else {
        uint count = 0;
        for (auto const& ps: pSurfSpecs[spec_gidx]) {
            if (ps->getOverlapTet_gidx() == tet_gidx) {
                count += 1;
            }
        }

        return count;
    }
}

////////////////////////////////////////////////////////////////////////////////
/*
std::vector<PointSpec *>
Vesicle::getSurfSpecs(solver::spec_global_id spec_gidx,
                                    tetrahedron_global_id tet_gidx) {

  std::vector<PointSpec *> pointspecs;
  if (pSurfSpecs.count(spec_gidx) == 0)
    return pointspecs;
  else {
    for (auto const &ps : pSurfSpecs[spec_gidx]) {
      if (ps->getOverlapTet_gidx() == tet_gidx) {
        pointspecs.emplace_back(ps);
      }
    }
    return pointspecs;
  }
}
////////////////////////////////////////////////////////////////////////////////
*/
void Vesicle::addSurfSpecs(const std::map<solver::spec_global_id, int>& specs,
                           tetrahedron_global_id tet_gidx) {
    // First go through and simply record positions of anything that's getting
    // removed:
    std::vector<math::position_rel_to_ves> positions;

    for (auto const& vs: specs) {
        solver::spec_global_id spec_gidx = vs.first;  // just for clarity

        uint starting_count = 0;

        if (pSurfSpecs.count(spec_gidx) != 0) {
            starting_count = pSurfSpecs[spec_gidx].size();
        }

        int starting_count_int = static_cast<int>(starting_count);
        if (starting_count_int + vs.second < 0) {
            std::ostringstream os;
            os << "Trying to set count on vesicle idx " << pIndex
               << " to negative number. Starting count: " << starting_count
               << ", trying to 'add': " << vs.second;
            ProgErrLog(os.str());
        }

        // Safe to make this unsigned since we have checked for negativity above
        uint end_count = starting_count + vs.second;

        if (starting_count > end_count) {
            uint counted = 0;

            for (uint i = 0; i < pSurfSpecs[spec_gidx].size(); ++i) {
                PointSpec* ps = pSurfSpecs[spec_gidx][i];
                if (ps->getOverlapTet_gidx() == tet_gidx) {
                    math::position_rel_to_ves pos_rel = ps->getPosCartesian();
                    positions.push_back(pos_rel);
                    counted += 1;
                    if (counted == (starting_count - end_count)) {
                        break;
                    }
                }
            }
            AssertLog(counted == (starting_count - end_count));
        }
    }

    for (auto const& vs: specs) {
        solver::spec_global_id spec_gidx = vs.first;  // just for clarity

        uint starting_count = 0;

        if (pSurfSpecs.count(spec_gidx) != 0) {
            starting_count = pSurfSpecs[spec_gidx].size();
        } else {
            pSurfSpecs[spec_gidx] = std::vector<PointSpec*>();
        }

        int end_count_int = starting_count + vs.second;

        // We have already checked this in previous loop, but just for extra safety
        // for now
        AssertLog(end_count_int >= 0);

        uint end_count = static_cast<uint>(end_count_int);

        if (end_count > starting_count) {
            for (uint i = starting_count; i < end_count; ++i) {
                auto* ps =
                    new PointSpec(spec_gidx, getDiam() / 2.0, solver()->getPointSpecNextIndex_());
                math::position_rel_to_ves pos_carts;

                pSurfSpecs[spec_gidx].push_back(ps);

                // The first conditional should always be valid, so should be OK to
                // include them in the while loop- an infinite while should NOT be
                // possible. Do a counter just in case
                uint attempts = 0;

                tetrahedron_global_id tetgidx;
                while (tetgidx.unknown()) {
                    attempts += 1;
                    if (attempts == 1000) {
                        std::ostringstream os;
                        os << "Failed to add surface species to vesicle index: " << pIndex
                           << ", too many failed attempts. ";
                        ProgErrLog(os.str());
                    }

                    if (positions.size() != 0) {
                        uint randomIndex = rng()->get() % positions.size();
                        pos_carts = positions[randomIndex];
                    } else {
                        pos_carts = (getDiam() / 2.0) * math::sphere_unit_randsurfpos(rng());
                    }

                    ps->setPosCartesian(pos_carts);

                    math::position_abs pos_abs{
                        pos_carts[0] + pPos[0],
                        pos_carts[1] + pPos[1],
                        pos_carts[2] + pPos[2],
                    };

                    // pointSpecs now have to be told in which tetrahedron they
                    // currently reside
                    tetgidx = getTetSpecOverlap(pos_abs, pTets_overlap_gidx);
                }

                ps->setOverlapTet_gidx(tetgidx);
            }
        } else if (starting_count > end_count) {
            if (tet_gidx.valid()) {
                uint erased = 0;
                for (uint i = 0; i < pSurfSpecs[vs.first].size(); ++i) {
                    PointSpec* ps = pSurfSpecs[vs.first][i];

                    if (ps->getOverlapTet_gidx() == tet_gidx) {
                        delete pSurfSpecs[spec_gidx][i];
                        pSurfSpecs[spec_gidx].erase(pSurfSpecs[vs.first].begin() + i);
                        erased += 1;
                        if (erased == (starting_count - end_count)) {
                            break;
                        }
                    }
                }
                if (erased != starting_count - end_count) {
                    // If we get here we couldn't remove species from the tet we wanted
                    // to. ERROR!
                    std::ostringstream os;
                    os << "Not enough species available in tetrahedron for application "
                          "of vesicle surface reaction. ";
                    ProgErrLog(os.str());
                }

            } else {
                for (uint i = starting_count; i > end_count; --i) {
                    delete pSurfSpecs[spec_gidx][-1];
                    pSurfSpecs[spec_gidx].pop_back();
                }
            }
        } else {
            std::ostringstream os;
            os << "No change in vesicle spec update. This is currently an error. ";
            ProgErrLog(os.str());
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::addOneSurfSpec(solver::spec_global_id spec_gidx,
                             solver::pointspec_individual_id spec_idx,
                             tetrahedron_global_id /*tet_gidx*/,
                             const math::position_abs& pos_abs) {
    // Do a quick check the tet index fits
    tetrahedron_global_id tet_gidx_check = getTetSpecOverlap(pos_abs, pTets_overlap_gidx);

    if (tet_gidx_check.unknown()) {
        ProgErrLog("Overlap tetrahedron cannot be found for surface species.");
    }

    // Could do a warning here if tet_gidx_check != tet_gidx? But it seems unnecessary- tet_gidx
    // argument is more of a hint, right 999,999,999 times out of a billion, but it doesn't really
    // matter if there's a mismatch as long as ovrelap tet can be found

    math::position_rel_to_ves pos_rel{pos_abs[0] - pPos[0],
                                      pos_abs[1] - pPos[1],
                                      pos_abs[2] - pPos[2]};

    auto* ps = new PointSpec(spec_gidx, getDiam() / 2.0, spec_idx);
    ps->setPosCartesian(pos_rel);

    ps->setOverlapTet_gidx(tet_gidx_check);

    pSurfSpecs[spec_gidx].emplace_back(ps);
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::addOneSurfSpec(solver::spec_global_id spec_gidx,
                             tetrahedron_global_id tet_gidx,
                             const math::position_abs& pos_abs) {
    addOneSurfSpec(spec_gidx, solver()->getPointSpecNextIndex_(), tet_gidx, pos_abs);
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::addOneSurfSpec(solver::spec_global_id spec_gidx, math::position_rel_to_ves pos_rel) {
    math::position_abs pos_abs{pos_rel[0] + pPos[0], pos_rel[1] + pPos[1], pos_rel[2] + pPos[2]};

    // Chrck it fits
    tetrahedron_global_id tet_gidx = getTetSpecOverlap(pos_abs, pTets_overlap_gidx);

    if (tet_gidx.unknown()) {
        ProgErrLog("Overlap tetrahedron cannot be found for surface species.");
    }

    auto* ps = new PointSpec(spec_gidx, getDiam() / 2.0, solver()->getPointSpecNextIndex_());
    ps->setPosCartesian(pos_rel);

    ps->setOverlapTet_gidx(tet_gidx);

    pSurfSpecs[spec_gidx].emplace_back(ps);
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::changeSurfSpecGidx(solver::spec_global_id spec_src_gidx,
                                 solver::spec_global_id spec_dst_gidx) {
    for (auto const& ps_src: pSurfSpecs[spec_src_gidx]) {
        ps_src->setSpecIndex(spec_dst_gidx);
        pSurfSpecs[spec_dst_gidx].emplace_back(ps_src);
    }
    pSurfSpecs[spec_src_gidx].clear();
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::clearSurfSpecs() {
    for (auto& sit: pSurfSpecs) {
        for (auto& s: sit.second) {
            delete s;
        }
        sit.second.clear();
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> Vesicle::getSurfaceSpecPosSpherical(
    solver::spec_global_id spec_gidx) const {
    std::vector<std::vector<double>> positions;
    auto it = pSurfSpecs.find(spec_gidx);
    if (it != pSurfSpecs.end()) {
        for (auto const& ps: it->second) {
            auto pos = ps->getPosSpherical();
            positions.push_back({pos.getTheta(), pos.getPhi()});
        }
    }
    return positions;
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::setSurfaceSpecPosSpherical(solver::spec_global_id spec_gidx,
                                         const std::vector<std::vector<double>>& pos_spherical) {
    uint ind = 0;
    for (auto const& ps: pSurfSpecs[spec_gidx]) {
        if (pos_spherical[ind].size() != 2) {
            std::ostringstream os;
            os << "Require list of length 2 for spherical coordinates.";
            throw steps::ArgErr(os.str());
        }

        ps->setPosSpherical(pos_spherical[ind][0], pos_spherical[ind][1]);
        if (++ind >= pos_spherical.size()) {
            break;
        }
    }
    if (ind < pos_spherical.size()) {
        CLOG(WARNING, "general_log")
            << "Tried to set the spherical position of" << pos_spherical.size()
            << "species but only" << pSurfSpecs[spec_gidx].size() << "exist on the vesicle surface."
            << std::endl;
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> Vesicle::getSurfaceSpecPos(solver::spec_global_id spec_gidx,
                                                            tetrahedron_global_id tet_gidx) const {
    std::vector<std::vector<double>> positions;

    auto specs = pSurfSpecs.find(spec_gidx);
    if (specs != pSurfSpecs.end()) {
        for (const auto& ps: specs->second) {
            if (ps->getOverlapTet_gidx() == tet_gidx) {
                const math::position_rel_to_ves& pos_rel = ps->getPosCartesian();
                // pointspecs hold postion relative to vesicle centre
                const std::vector<double> carts{pos_rel[0] + pPos[0],
                                                pos_rel[1] + pPos[1],
                                                pos_rel[2] + pPos[2]};
                positions.emplace_back(std::move(carts));
            }
        }
    }

    return positions;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> Vesicle::getSurfaceSpecPos(
    solver::spec_global_id spec_gidx) const {
    std::vector<std::vector<double>> positions;

    auto specs = pSurfSpecs.find(spec_gidx);
    if (specs != pSurfSpecs.end()) {
        positions.reserve(specs->second.size());
        for (const auto& ps: specs->second) {
            const math::position_rel_to_ves& pos_rel = ps->getPosCartesian();
            // pointspecs hold position relative to vesicle centre
            std::vector<double> carts{pos_rel[0] + pPos[0],
                                      pos_rel[1] + pPos[1],
                                      pos_rel[2] + pPos[2]};
            positions.emplace_back(std::move(carts));
        }
    }

    return positions;
}


////////////////////////////////////////////////////////////////////////////////

std::vector<solver::pointspec_individual_id> Vesicle::getSurfaceSpecIndices(
    solver::spec_global_id spec_gidx) const {
    std::vector<solver::pointspec_individual_id> indices;
    for (auto const& ss: pSurfSpecs) {
        if (ss.first == spec_gidx) {
            for (auto const& ps: ss.second) {
                indices.emplace_back(ps->getUniqueIndex());
            }
        }
    }

    return indices;
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::setInnerSpecCount(solver::spec_global_id spec_gidx, uint count) {
    // Since this is called from API, simply overwrite anything already there
    pInnerSpecCount[spec_gidx] = count;
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::incInnerSpecCount(solver::spec_global_id spec_gidx, uint count) {
    if (pInnerSpecCount.count(spec_gidx) == 0) {
        pInnerSpecCount[spec_gidx] = count;
    } else {
        pInnerSpecCount[spec_gidx] += count;
    }
}

////////////////////////////////////////////////////////////////////////////////

uint Vesicle::getInnerSpecCount(solver::spec_global_id spec_gidx) const {
    const auto spec_it = pInnerSpecCount.find(spec_gidx);
    if (spec_it != pInnerSpecCount.end()) {
        return spec_it->second;
    } else {
        return 0;
    }
}

////////////////////////////////////////////////////////////////////////////////
/*
uint Vesicle::getLinkSpecCount(
    solver::linkspec_global_id spec_gidx, tetrahedron_global_id tet_gidx) {
  uint count = 0;

  for (auto const &ls : pLinkSpecs) {
    AssertLog(ls.second->getVesicle() == this); // TODO remove once tested

    if (ls.second->getGidx() == spec_gidx &&
        ls.second->getOverlapTet_gidx() == tet_gidx) {
      count += 1;
    }
  }

  return count;
}
*/
////////////////////////////////////////////////////////////////////////////////

uint Vesicle::getLinkSpecCount(solver::linkspec_global_id spec_gidx) const {
    uint count = 0;
    for (auto const& ls: pLinkSpecs) {
        AssertLog(ls.second->getVesicle() == this);  // TODO remove once tested

        if (ls.second->getGidx() == spec_gidx) {
            count += 1;
        }
    }

    return count;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::linkspec_individual_id> Vesicle::getLinkSpecIndices(
    solver::linkspec_global_id spec_gidx) const {
    std::vector<solver::linkspec_individual_id> indices;
    for (auto const& ls: pLinkSpecs) {
        if (ls.second->getGidx() == spec_gidx) {
            indices.emplace_back(ls.first);
        }
    }

    return indices;
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::addLinkSpec(solver::linkspec_individual_id linkspec_uniqueid, LinkSpec* link_spec) {
    if (pLinkSpecs.count(linkspec_uniqueid) != 0) {
        ProgErrLog("LinkSpec already added to Vesicle");
    }
    if (linkspec_uniqueid != link_spec->getUniqueID()) {
        ProgErrLog("LinkSpec mixup.");
    }
    pLinkSpecs[linkspec_uniqueid] = link_spec;
}

void Vesicle::updateLinkSpec(solver::linkspec_individual_id linkspec_uniqueid,
                             solver::LinkSpecdef* linkspec_def) {
    if (pLinkSpecs.count(linkspec_uniqueid) == 0) {
        ProgErrLog("LinkSpec not known in Vesicle");
    }

    pLinkSpecs[linkspec_uniqueid]->setDef(linkspec_def);
}

////////////////////////////////////////////////////////////////////////////////
/*
void Vesicle::remLinkSpec(
    solver::linkspec_individual_id linkspec_uniqueid) {
  if (pLinkSpecs.count(linkspec_uniqueid) == 0) {
    ProgErrLog("LinkSpec not known in Vesicle");
  }
  pLinkSpecs.erase(linkspec_uniqueid);
}
*/
////////////////////////////////////////////////////////////////////////////////

void Vesicle::remLinkSpec(solver::linkspec_individual_id linkspec_uniqueid, LinkSpec* linkspec) {
    if (pLinkSpecs.count(linkspec_uniqueid) == 0) {
        ProgErrLog("LinkSpec not known in Vesicle");
    }
    if (pLinkSpecs[linkspec_uniqueid] != linkspec) {
        ProgErrLog("LinkSpec mismatch");
    }
    pLinkSpecs.erase(linkspec_uniqueid);

    delete (linkspec);
}

////////////////////////////////////////////////////////////////////////////////

LinkSpec* Vesicle::getLinkSpec(solver::linkspec_individual_id linkspec_uniqueid) const {
    const auto specs = pLinkSpecs.find(linkspec_uniqueid);
    if (specs != pLinkSpecs.end()) {
        ProgErrLog("LinkSpec not known in Vesicle");
    }
    return specs->second;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> Vesicle::getLinkSpecPos(
    solver::linkspec_global_id linkspec_gidx) const {
    std::vector<std::vector<double>> positions;

    for (auto const& bs: pLinkSpecs) {
        if (bs.second->getGidx() == linkspec_gidx) {
            math::position_abs pos_3d = bs.second->getPosCartesian_abs();

            std::vector<double> pos_vec;
            pos_vec.emplace_back(pos_3d[0]);
            pos_vec.emplace_back(pos_3d[1]);
            pos_vec.emplace_back(pos_3d[2]);

            positions.emplace_back(pos_vec);
        }
    }

    return positions;
}

////////////////////////////////////////////////////////////////////////////////

bool Vesicle::linkSpecMoveAllowed(const math::point3d& move_vector) {
    // Need to check against all link species

    bool moveallowed;

    for (auto const& bs: pLinkSpecs) {
        moveallowed = bs.second->movePosAllowed(move_vector);
        if (not moveallowed) {
            return false;
        }
    }

    // means we are ok
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::setPathPositions(
    const std::vector<std::pair<double, math::position_abs>>& path_positions) {
    AssertLog(pPathPositions.empty());
    AssertLog(!path_positions.empty());
    AssertLog(!pOnPath);

    pPathPositions = path_positions;
    pPathPosition_index = 0;
    pPathPosition_index_next = 0;
    pTime_accum = 0.0;
    pTime_accum_next = 0.0;
    pOnPath = true;
}

////////////////////////////////////////////////////////////////////////////////

math::position_abs Vesicle::getNextPosition(double dt) noexcept {
    double deltas_summed = 0.0;  // have to sum any moves that we apply

    for (uint i = pPathPosition_index + 1; i < pPathPositions.size(); ++i) {
        if (deltas_summed + (pPathPositions[i].first - pTime_accum) > dt) {
            // We stop because next move is still in the future
            pTime_accum_next = pTime_accum + (dt - deltas_summed);
            return pPathPositions[pPathPosition_index_next].second;
        } else {
            // We have a hit, a move that is before the next dt
            // May get multiple hits. May be the last point
            deltas_summed += pPathPositions[i].first;
            pPathPosition_index_next = i;
        }
    }

    // May be at the end
    return pPathPositions[pPathPosition_index_next].second;
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::nextPositionOnPath() {
    pPathPosition_index = pPathPosition_index_next;
    pTime_accum = pTime_accum_next;

    if (pPathPosition_index == pPathPositions.size() - 1) {
        // Path has come to an end.
        removeFromPath();
    }
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::removeFromPath() {
    pPathPositions.clear();
    pPathPosition_index = 0;
    pPathPosition_index_next = 0;

    pTime_accum = 0.0;
    pTime_accum_next = 0.0;
    pOnPath = false;
}

////////////////////////////////////////////////////////////////////////////////

solver::exocytosis_global_id Vesicle::appliedExocytosis() {
    if (pAppliedExocytosis.empty()) {
        return {};
    } else {
        // First come, first served. TODO: random selection?
        return *pAppliedExocytosis.begin();
    }
}

}  // namespace steps::mpi::tetvesicle
