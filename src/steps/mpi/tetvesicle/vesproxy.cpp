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

#include "mpi/tetvesicle/vesproxy.hpp"

// STEPS headers.
#include "math/constants.hpp"
#include "math/point.hpp"
#include "math/sphere.hpp"
#include "mpi/tetvesicle/tetvesicle_rdef.hpp"
#include "solver/vesicledef.hpp"

namespace steps::mpi::tetvesicle {

VesProxy::VesProxy(solver::Vesicledef* vesdef,
                   TetRDEF* tet,
                   solver::vesicle_individual_id unique_index,
                   const math::position_abs& pos,
                   bool contains_link)
    : pDef(vesdef)
    , pUniqueIndex(unique_index)
    , pTet(tet)
    , pImmobilityUpdate(0)
    , pExoApplied{} {
    AssertLog(vesdef != nullptr);
    AssertLog(pTet != nullptr);

    pVesicleIndex = vesdef->gidx();

    pPos = pos;  // error-checking??

    pContainsLink = contains_link;
}

////////////////////////////////////////////////////////////////////////////////

VesProxy::~VesProxy() = default;

////////////////////////////////////////////////////////////////////////////////

void VesProxy::checkpoint(std::fstream& /*cp_file*/) {
    // Reserve. Nothing to do here because only created when vesicles are created
}

////////////////////////////////////////////////////////////////////////////////

void VesProxy::restore(std::fstream& /*cp_file*/) {
    // Reserve.
}

////////////////////////////////////////////////////////////////////////////////

void VesProxy::updImmobility(int mob_upd) {
    pImmobilityUpdate += mob_upd;
}

////////////////////////////////////////////////////////////////////////////////

inline const rng::RNGptr& VesProxy::rng() const noexcept {
    return getTet()->getCompRDEF()->rng();
}  // long-winded

////////////////////////////////////////////////////////////////////////////////

tetmesh::Tetmesh* VesProxy::mesh() const noexcept {
    return getTet()->getCompRDEF()->mesh();
}

////////////////////////////////////////////////////////////////////////////////

const std::vector<std::pair<solver::pointspec_individual_id, math::position_abs>>&
VesProxy::getSurfaceSpecPos_absolute(solver::spec_global_id spec_gidx) const {
    // it's possible this spec doesn't exist so have to allow for a non-entry in
    // the map

    auto specs = pSpecs_V.find(spec_gidx);
    if (specs != pSpecs_V.end()) {
        return specs->second;
    } else {
        return pEmptyPosVec;
    }
}

////////////////////////////////////////////////////////////////////////////////

void VesProxy::removeOneSurfSpec(solver::spec_global_id spec_gidx, uint array_index) {
    if (array_index >= pSpecs_V[spec_gidx].size()) {
        std::ostringstream os;
        os << "Cannot remove spec global index " << spec_gidx << " from vesicle idx "
           << pUniqueIndex << ". Internal error. Contact STEPS support!!!";
        ProgErrLog(os.str());
    }

    pSpecs_V[spec_gidx].erase(pSpecs_V[spec_gidx].begin() + array_index);
}

////////////////////////////////////////////////////////////////////////////////

void VesProxy::addSurfSpec(solver::spec_global_id spec_gidx,
                           solver::pointspec_individual_id idx,
                           const math::position_abs& pos_abs) {
    pSpecs_V[spec_gidx].emplace_back(idx, pos_abs);

    // Store in known positions
    pAllPositions.emplace_back(pos_abs);
}

////////////////////////////////////////////////////////////////////////////////

bool VesProxy::addSurfSpecs(const std::map<solver::spec_global_id, int>& specs) {
    // To allow this method to fail (instead of stopping with an error) it's necessary to check, if
    // we have a positive update, whether we need to generate a point or not
    bool pos_update = false;
    std::vector<math::position_abs> removed_pos;
    for (auto const& vs: specs) {
        if (vs.second > 0) {
            pos_update = true;
            continue;
        } else if (vs.second < 0) {
            auto it = pSpecs_V.find(vs.first);
            if (it != pSpecs_V.end() and static_cast<int>(it->second.size()) + vs.second >= 0) {
                for (int i = 0; i < -vs.second; ++i) {
                    removed_pos.emplace_back(it->second[i].second);
                }
            } else {
                int starting_count = (it != pSpecs_V.end()) ? it->second.size() : 0;
                std::ostringstream os;
                os << "Trying to set count on vesicle idx " << pVesicleIndex
                   << " to negative number. Starting count: " << starting_count
                   << ", trying to 'add': " << vs.second;
                ProgErrLog(os.str());
            }
        }
    }

    if (pos_update && pAllPositions.size() == 0) {
        // Need to generate a random point on the surface. If we can't, return a false and the
        // calling function can deal with that as it may
        uint attempts = 0;
        math::position_abs pos_carts_absolute;
        while (true) {
            if (attempts == 1000000) {
                // allow many attempts since we now check known positions so this
                // loop should only need to be done maximum one time per vesicle_dt
                CLOG(WARNING, "general_log")
                    << "Failed to add surface species to vesicle index: " << pVesicleIndex
                    << ", too many failed attempts.\n";
                return false;
            }
            attempts += 1;

            math::position_rel_to_ves pos_carts = (getDiam() / 2.0) *
                                                  math::sphere_unit_randsurfpos(rng());

            pos_carts_absolute[0] = pos_carts[0] + pPos[0];
            pos_carts_absolute[1] = pos_carts[1] + pPos[1];
            pos_carts_absolute[2] = pos_carts[2] + pPos[2];

            if (mesh()->isPointInTet(pos_carts_absolute, getTetGidx())) {
                pAllPositions.emplace_back(pos_carts_absolute);
                break;
            }
        }
    }

    uint all_pos_size = pAllPositions.size();
    if (pos_update) {
        AssertLog(all_pos_size > 0);
    }

    for (auto const& vs: specs) {
        solver::spec_global_id spec_gidx = vs.first;  // just for clarity

        uint starting_count = 0;

        if (pSpecs_V.count(spec_gidx) != 0) {
            starting_count = pSpecs_V[spec_gidx].size();
        } else {
            pSpecs_V[spec_gidx] =
                std::vector<std::pair<solver::pointspec_individual_id, math::position_abs>>();
        }

        int end_count_int = starting_count + vs.second;

        // We have already checked this in previous loop, but just for extra safety
        // for now
        AssertLog(end_count_int >= 0);

        uint end_count = static_cast<uint>(end_count_int);
        if (end_count > starting_count) {
            for (uint i = starting_count; i < end_count; ++i) {
                auto idx = pTet->solverRDEF()->getNextPointspecIndividualID_();
                if (not removed_pos.empty()) {
                    pSpecs_V[spec_gidx].emplace_back(idx, removed_pos.back());
                    removed_pos.pop_back();
                } else {
                    pSpecs_V[spec_gidx].emplace_back(idx,
                                                     pAllPositions[rng()->get() % all_pos_size]);
                }
            }
        } else if (starting_count > end_count) {
            uint num_to_erase = starting_count - end_count;
            pSpecs_V[spec_gidx].erase(pSpecs_V[spec_gidx].begin(),
                                      pSpecs_V[spec_gidx].begin() + num_to_erase);
        } else {
            std::ostringstream os;
            os << "No change in vesicle spec update. This is currently an error. ";
            ProgErrLog(os.str());
        }
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void VesProxy::addLinkSpec(solver::linkspec_individual_id ls_unique_idx,
                           solver::linkspec_global_id ls_gidx,
                           math::position_abs ls_pos_abs) {
    AssertLog(pLinkSpecs_V.count(ls_unique_idx) == 0);

    pLinkSpecs_V.emplace(ls_unique_idx, std::make_pair(ls_gidx, ls_pos_abs));
}

////////////////////////////////////////////////////////////////////////////////

uint VesProxy::getLinkSpecCount_V(solver::linkspec_global_id linkspec_gidx) const {
    AssertLog(linkspec_gidx < def()->countLinkSpecs_V());

    uint count = 0;

    for (auto const& ls_uidx: pLinkSpecs_V) {
        if (ls_uidx.second.first == linkspec_gidx) {
            count += 1;
        }
    }

    return count;
}

////////////////////////////////////////////////////////////////////////////////

void VesProxy::changeLinkSpecGidx(solver::linkspec_global_id linkspecbefore_gidx,
                                  solver::linkspec_global_id linkspecafter_gidx) {
    AssertLog(linkspecbefore_gidx < def()->countLinkSpecs_V());
    AssertLog(linkspecafter_gidx < def()->countLinkSpecs_V());

    for (auto& ls_uidx: pLinkSpecs_V) {
        // Simply, for now, doing it on the first match
        if (ls_uidx.second.first == linkspecbefore_gidx) {
            // We have found a match. Same position (ls_uidx.second.second)
            const auto& ls_newpair = std::make_pair(linkspecafter_gidx, ls_uidx.second.second);
            ls_uidx.second = ls_newpair;

            reqLinkSpecUpd(ls_uidx.first);

            return;
        }
    }

    // if we get here something went horribly wrong
    ProgErrLog("Unable to change LinkSpec ID in VesProxy. ");
}

////////////////////////////////////////////////////////////////////////////////
/*
void VesProxy::setLinkSpecCount_V(solver::linkspec_global_id
linkspec_gidx, uint count)
{
    // Using global indices but this should not be an issue. Check anyway for
safety. AssertLog(linkspec_gidx < def()->countLinkSpecs_V());
    pLinkSpecCount_V[linkspec_gidx.get()] = count;
}
*/

////////////////////////////////////////////////////////////////////////////////

void VesProxy::incSpecCount_I(solver::spec_global_id spec_gidx, uint count) {
    // Using global indices but this should not be an issue. Check anyway for
    // safety.
    AssertLog(spec_gidx < def()->countSpecs_I());

    if (pSpecCount_I.find(spec_gidx) == pSpecCount_I.end()) {
        pSpecCount_I[spec_gidx] += count;
    } else {
        pSpecCount_I[spec_gidx] = count;
    }
}

////////////////////////////////////////////////////////////////////////////////

void VesProxy::applyExo(solver::exocytosis_global_id exo_gidx) {
    if (pExoApplied.valid()) {
        ProgErrLog("Multiple exocytosis events not allowed for VesProxy. ");
    }
    pExoApplied = exo_gidx;
}

}  // namespace steps::mpi::tetvesicle
