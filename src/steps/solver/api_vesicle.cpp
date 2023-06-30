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
#include "solver/api.hpp"
#include "solver/compdef.hpp"
#include "solver/endocytosisdef.hpp"
#include "solver/exocytosisdef.hpp"
#include "solver/fwd.hpp"
#include "solver/patchdef.hpp"
#include "solver/raftendocytosisdef.hpp"
#include "solver/statedef.hpp"
#include "util/common.hpp"
#include "util/error.hpp"
#include "util/vocabulary.hpp"

// logging
#include "easylogging++.h"

namespace steps::solver {

uint API::getCompVesicleCount(std::string const& c, std::string const& v) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);

    return _getCompVesicleCount(cidx, vidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompVesicleCount(std::string const& c, std::string const& v, uint n) {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);

    _setCompVesicleCount(cidx, vidx, n);
}

////////////////////////////////////////////////////////////////////////////////

vesicle_individual_id API::addCompVesicle(std::string const& c, std::string const& v) {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);

    return _addCompVesicle(cidx, vidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::deleteSingleVesicle(std::string const& v, vesicle_individual_id ves_unique_index) {
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);

    _deleteSingleVesicle(vidx, ves_unique_index);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getSingleVesicleSurfaceLinkSpecCount(std::string const& v,
                                               vesicle_individual_id ves_unique_index,
                                               std::string const& ls) const {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    linkspec_global_id lsidx = pStatedef->getLinkSpecIdx(ls);

    return _getSingleVesicleSurfaceLinkSpecCount(vidx, ves_unique_index, lsidx);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<linkspec_individual_id> API::getSingleVesicleSurfaceLinkSpecIndices(
    std::string const& v,
    vesicle_individual_id ves_unique_index,
    std::string const& ls) const {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    linkspec_global_id lsidx = pStatedef->getLinkSpecIdx(ls);

    return _getSingleVesicleSurfaceLinkSpecIndices(vidx, ves_unique_index, lsidx);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<pointspec_individual_id> API::getSingleVesicleSurfaceSpecIndices(
    std::string const& v,
    vesicle_individual_id ves_unique_index,
    std::string const& s) const {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getSingleVesicleSurfaceSpecIndices(vidx, ves_unique_index, sidx);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<vesicle_individual_id> API::getCompVesicleIndices(std::string const& c,
                                                              std::string const& v) const {
    comp_global_id cidx = pStatedef->getCompIdx(c);
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);

    return _getCompVesicleIndices(cidx, vidx);
}

////////////////////////////////////////////////////////////////////////////////

std::string API::getSingleVesicleCompartment(std::string const& v,
                                             vesicle_individual_id ves_unique_index) const {
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);

    comp_global_id cidx = _getSingleVesicleCompartment(vidx, ves_unique_index);

    if (cidx.unknown()) {
        return {};
    } else {
        return pStatedef->compdef(cidx)->name();
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getSingleVesiclePos(std::string const& v,
                                             vesicle_individual_id ves_unique_index) const {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);

    return _getSingleVesiclePos(vidx, ves_unique_index);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompSingleVesiclePos(std::string const& c,
                                  std::string const& v,
                                  vesicle_individual_id ves_unique_index,
                                  const std::vector<double>& pos,
                                  bool force) {
    if (pos.size() != 3) {
        std::ostringstream os;
        os << "Position argument must be sequence of length 3.";
        throw steps::ArgErr(os.str());
    }

    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);

    _setCompSingleVesiclePos(cidx, vidx, ves_unique_index, pos, force);
}

////////////////////////////////////////////////////////////////////////////////

std::map<vesicle_individual_id, uint> API::getCompVesicleSurfaceSpecCountDict(
    std::string const& c,
    std::string const& v,
    std::string const& s) const

{
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getCompVesicleSurfaceSpecCountMap(cidx, vidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getCompVesicleSurfaceSpecCount(std::string const& c,
                                         std::string const& v,
                                         std::string const& s) const

{
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getCompVesicleSurfaceSpecCount(cidx, vidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getCompVesicleInnerSpecCount(std::string const& c,
                                       std::string const& v,
                                       std::string const& s) const

{
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getCompVesicleInnerSpecCount(cidx, vidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getSingleVesicleSurfaceSpecCount(std::string const& v,
                                           vesicle_individual_id ves_unique_index,
                                           std::string const& s) const {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getSingleVesicleSurfaceSpecCount(vidx, ves_unique_index, sidx);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getSingleVesicleInnerSpecCount(std::string const& v,
                                         vesicle_individual_id ves_unique_index,
                                         std::string const& s) const {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getSingleVesicleInnerSpecCount(vidx, ves_unique_index, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setSingleVesicleSurfaceSpecCount(std::string const& v,
                                           vesicle_individual_id ves_unique_index,
                                           std::string const& s,
                                           uint count) {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    _setSingleVesicleSurfaceSpecCount(vidx, ves_unique_index, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> API::getSingleVesicleSurfaceSpecPos(
    std::string const& v,
    vesicle_individual_id ves_unique_index,
    std::string const& s) {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getSingleVesicleSurfaceSpecPos(vidx, ves_unique_index, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setSingleVesicleInnerSpecCount(std::string const& v,
                                         vesicle_individual_id ves_unique_index,
                                         std::string const& s,
                                         uint count) {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    _setSingleVesicleInnerSpecCount(vidx, ves_unique_index, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> API::getSingleVesicleSurfaceSpecPosSpherical(
    std::string const& v,
    vesicle_individual_id ves_unique_index,
    std::string const& s) const {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getSingleVesicleSurfaceSpecPosSpherical(vidx, ves_unique_index, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setSingleVesicleSurfaceSpecPosSpherical(
    std::string const& v,
    vesicle_individual_id ves_unique_index,
    std::string const& s,
    const std::vector<std::vector<double>>& pos_spherical) {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    _setSingleVesicleSurfaceSpecPosSpherical(vidx, ves_unique_index, sidx, pos_spherical);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getSingleSpecPosSpherical(std::string const& s,
                                                   pointspec_individual_id ps_unique_id) const {
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getSingleSpecPosSpherical(sidx, ps_unique_id);
}

////////////////////////////////////////////////////////////////////////////////

std::map<vesicle_individual_id, uint> API::getCompVesicleSurfaceLinkSpecCountDict(
    std::string const& c,
    std::string const& v,
    std::string const& ls) const

{
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    linkspec_global_id lsidx = pStatedef->getLinkSpecIdx(ls);

    return _getCompVesicleSurfaceLinkSpecCountMap(cidx, vidx, lsidx);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getCompVesicleSurfaceLinkSpecCount(std::string const& c,
                                             std::string const& v,
                                             std::string const& ls) const

{
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    linkspec_global_id lsidx = pStatedef->getLinkSpecIdx(ls);

    return _getCompVesicleSurfaceLinkSpecCount(cidx, vidx, lsidx);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> API::getSingleVesicleSurfaceLinkSpecPos(
    std::string const& v,
    vesicle_individual_id ves_unique_index,
    std::string const& ls) {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    linkspec_global_id lsidx = pStatedef->getLinkSpecIdx(ls);

    return _getSingleVesicleSurfaceLinkSpecPos(vidx, ves_unique_index, lsidx);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getSingleLinkSpecPos(linkspec_individual_id ls_unique_id) const {
    return _getSingleLinkSpecPos(ls_unique_id);
}

////////////////////////////////////////////////////////////////////////////////

linkspec_individual_id API::getSingleLinkSpecLinkedTo(linkspec_individual_id ls_unique_id) const {
    return _getSingleLinkSpecLinkedTo(ls_unique_id);
}

////////////////////////////////////////////////////////////////////////////////

vesicle_individual_id API::getSingleLinkSpecVes(linkspec_individual_id ls_unique_id) const {
    return _getSingleLinkSpecVes(ls_unique_id);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getSingleVesicleImmobility(std::string const& v,
                                     vesicle_individual_id ves_unique_index) const {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);

    return _getSingleVesicleImmobility(vidx, ves_unique_index);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::tetrahedron_global_id> API::getSingleVesicleOverlapTets(
    std::string const& v,
    vesicle_individual_id ves_unique_index) const {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);

    return _getSingleVesicleOverlapTets(vidx, ves_unique_index);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetVesicleDcst(tetrahedron_global_id tidx, std::string const& v, double dcst) {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);

    _setTetVesicleDcst(tidx, vidx, dcst);
}

////////////////////////////////////////////////////////////////////////////////

void API::setVesicleSurfaceLinkSpecSDiffD(std::string const& v,
                                          std::string const& ls,
                                          double dcst) {
    // the following may throw exceptions if strings are unknown
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);
    linkspec_global_id lsidx = pStatedef->getLinkSpecIdx(ls);

    _setVesicleSurfaceLinkSpecSDiffD(vidx, lsidx, dcst);
}

////////////////////////////////////////////////////////////////////////////////
/*
void API::setCompVesicleSpecDiffD(std::string const & c, std::string const & v,
std::string const & s, double d) const
{
    // the following may throw exceptions if strings are unknown
    uint cidx = pStatedef->getCompIdx(c);
    uint vidx = pStatedef->getVesicleIdx(v);
    uint sidx = pStatedef->getSpecIdx(s);

    _setCompVesicleSpecDiffD(cidx, vidx, sidx, d);

}
*/

////////////////////////////////////////////////////////////////////////////////

void API::setVesSReacK(std::string const& vsr, double kf) {
    ArgErrLogIf(kf < 0.0, "Reaction constant cannot be negative.");

    vessreac_global_id vsridx = pStatedef->getVesSReacIdx(vsr);

    return _setVesSReacK(vsridx, kf);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getVesSReacExtent(std::string const& vsr) const {
    vessreac_global_id vsridx = pStatedef->getVesSReacIdx(vsr);

    return _getVesSReacExtent(vsridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setExocytosisK(std::string const& exo, double kf) {
    exocytosis_global_id exoidx = pStatedef->getExocytosisIdx(exo);

    return _setExocytosisK(exoidx, kf);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getExocytosisExtent(std::string const& exo) const {
    exocytosis_global_id exoidx = pStatedef->getExocytosisIdx(exo);

    return _getExocytosisExtent(exoidx);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<ExocytosisEvent> API::getExocytosisEvents(std::string const& exo) {
    exocytosis_global_id exoidx = pStatedef->getExocytosisIdx(exo);

    return _getExocytosisEvents(exoidx);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getRaftEndocytosisExtent(std::string const& rendo) const {
    raftendocytosis_global_id rendoidx = pStatedef->getRaftEndocytosisIdx(rendo);

    return _getRaftEndocytosisExtent(rendoidx);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<RaftEndocytosisEvent> API::getRaftEndocytosisEvents(std::string const& rendo) {
    raftendocytosis_global_id rendoidx = pStatedef->getRaftEndocytosisIdx(rendo);

    return _getRaftEndocytosisEvents(rendoidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setRaftEndocytosisK(std::string const& rendo, double kcst) {
    ArgErrLogIf(kcst < 0.0, "Reaction constant cannot be negative.");

    raftendocytosis_global_id rendoidx = pStatedef->getRaftEndocytosisIdx(rendo);

    _setRaftEndocytosisK(rendoidx, kcst);
}

////////////////////////////////////////////////////////////////////////////////

void API::addVesicleDiffusionGroup(std::string const& v, const std::vector<std::string>& comps) {
    vesicle_global_id vidx = pStatedef->getVesicleIdx(v);

    std::vector<comp_global_id> compindices;

    for (auto const& c: comps) {
        comp_global_id cidx = pStatedef->getCompIdx(c);
        compindices.emplace_back(cidx);
    }

    _addVesicleDiffusionGroup(vidx, compindices);
}

////////////////////////////////////////////////////////////////////////////////

void API::createPath(std::string const& /*id*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::addPathPoint(std::string const& /*path_name*/,
                       uint /*point_id*/,
                       const std::vector<double>& /*position*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::addPathBranch(std::string const& /*path_name*/,
                        uint /*point_id*/,
                        const std::map<uint, double>& /*dest_points*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::map<std::string, std::map<uint, std::pair<std::vector<double>, std::map<uint, double>>>>
API::getAllPaths() const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::getPatchRaftCount(std::string const& p, std::string const& r) const {
    // the following may throw exceptions if strings are unknown
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    raft_global_id ridx = pStatedef->getRaftIdx(r);

    return _getPatchRaftCount(pidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchRaftCount(std::string const& p, std::string const& r, uint n) {
    // the following may throw exceptions if strings are unknown
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    raft_global_id ridx = pStatedef->getRaftIdx(r);

    _setPatchRaftCount(pidx, ridx, n);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getSingleRaftPos(std::string const& r,
                                          raft_individual_id raft_unique_index) const {
    // the following may throw exceptions if strings are unknown
    raft_global_id ridx = pStatedef->getRaftIdx(r);

    return _getSingleRaftPos(ridx, raft_unique_index);
}

////////////////////////////////////////////////////////////////////////////////

std::map<raft_individual_id, uint> API::getPatchRaftSpecCountDict(std::string const& p,
                                                                  std::string const& r,
                                                                  std::string const& s) const

{
    // the following may throw exceptions if strings are unknown
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    raft_global_id ridx = pStatedef->getRaftIdx(r);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getPatchRaftSpecCountMap(pidx, ridx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getPatchRaftSpecCount(std::string const& p,
                                std::string const& r,
                                std::string const& s) const

{
    // the following may throw exceptions if strings are unknown
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    raft_global_id ridx = pStatedef->getRaftIdx(r);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getPatchRaftSpecCount(pidx, ridx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getSingleRaftSpecCount(std::string const& r,
                                 raft_individual_id raft_unique_index,
                                 std::string const& s) const {
    // the following may throw exceptions if strings are unknown
    raft_global_id ridx = pStatedef->getRaftIdx(r);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getSingleRaftSpecCount(ridx, raft_unique_index, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setSingleRaftSpecCount(std::string const& r,
                                 raft_individual_id raft_unique_index,
                                 std::string const& s,
                                 uint count) {
    // the following may throw exceptions if strings are unknown
    raft_global_id ridx = pStatedef->getRaftIdx(r);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    _setSingleRaftSpecCount(ridx, raft_unique_index, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getSingleRaftImmobility(std::string const& r,
                                  raft_individual_id raft_unique_index) const {
    // the following may throw exceptions if strings are unknown
    raft_global_id ridx = pStatedef->getRaftIdx(r);

    return _getSingleRaftImmobility(ridx, raft_unique_index);
}

////////////////////////////////////////////////////////////////////////////////

double API::getSingleRaftRaftEndocytosisK(std::string const& r,
                                          raft_individual_id raft_unique_index,
                                          std::string const& rendo) const {
    raft_global_id ridx = pStatedef->getRaftIdx(r);
    raftendocytosis_global_id rendoidx = pStatedef->getRaftEndocytosisIdx(rendo);
    return _getSingleRaftRaftEndocytosisK(ridx, raft_unique_index, rendoidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setSingleRaftRaftEndocytosisK(std::string const& r,
                                        raft_individual_id raft_unique_index,
                                        std::string const& rendo,
                                        double k) {
    raft_global_id ridx = pStatedef->getRaftIdx(r);
    raftendocytosis_global_id rendoidx = pStatedef->getRaftEndocytosisIdx(rendo);
    _setSingleRaftRaftEndocytosisK(ridx, raft_unique_index, rendoidx, k);
}

////////////////////////////////////////////////////////////////////////////////

void API::setSingleRaftSReacActive(std::string const& r,
                                   raft_individual_id raft_unique_index,
                                   std::string const& rsreac,
                                   bool active) {
    raft_global_id ridx = pStatedef->getRaftIdx(r);
    raftsreac_global_id rsreacidx = pStatedef->getRaftSReacIdx(rsreac);
    _setSingleRaftSReacActive(ridx, raft_unique_index, rsreacidx, active);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getSingleRaftSReacActive(std::string const& r,
                                   raft_individual_id raft_unique_index,
                                   std::string const& rsreac) const {
    raft_global_id ridx = pStatedef->getRaftIdx(r);
    raftsreac_global_id rsreacidx = pStatedef->getRaftSReacIdx(rsreac);
    return _getSingleRaftSReacActive(ridx, raft_unique_index, rsreacidx);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<raft_individual_id> API::getPatchRaftIndices(std::string const& p,
                                                         std::string const& r) const {
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    raft_global_id ridx = pStatedef->getRaftIdx(r);

    return _getPatchRaftIndices(pidx, ridx);
}


////////////////////////////////////////////////////////////////////////////////

std::string API::getSingleRaftPatch(std::string const& r,
                                    raft_individual_id raft_unique_index) const {
    raft_global_id ridx = pStatedef->getRaftIdx(r);

    patch_global_id pidx = _getSingleRaftPatch(ridx, raft_unique_index);
    if (pidx.unknown()) {
        return {};
    } else {
        return pStatedef->patchdef(pidx)->name();
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchEndocyticZoneEndocytosisActive(std::string const& patch,
                                                 std::string const& zone,
                                                 std::string const& endo,
                                                 bool active) {
    patch_global_id pidx = pStatedef->getPatchIdx(patch);
    endocytosis_global_id endoidx = pStatedef->getEndocytosisIdx(endo);
    _setPatchEndocyticZoneEndocytosisActive(pidx, zone, endoidx, active);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchEndocyticZoneEndocytosisK(std::string const& patch,
                                            std::string const& zone,
                                            std::string const& endo,
                                            double k) {
    patch_global_id pidx = pStatedef->getPatchIdx(patch);
    endocytosis_global_id endoidx = pStatedef->getEndocytosisIdx(endo);
    _setPatchEndocyticZoneEndocytosisK(pidx, zone, endoidx, k);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getPatchEndocyticZoneEndocytosisExtent(std::string const& patch,
                                                 std::string const& zone,
                                                 std::string const& endo) const {
    patch_global_id pidx = pStatedef->getPatchIdx(patch);
    endocytosis_global_id endoidx = pStatedef->getEndocytosisIdx(endo);
    return _getPatchEndocyticZoneEndocytosisExtent(pidx, zone, endoidx);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<EndocytosisEvent> API::getPatchEndocyticZoneEndocytosisEvents(
    std::string const& patch,
    std::string const& zone,
    std::string const& endo) const {
    patch_global_id pidx = pStatedef->getPatchIdx(patch);
    endocytosis_global_id endoidx = pStatedef->getEndocytosisIdx(endo);
    return _getPatchEndocyticZoneEndocytosisEvents(pidx, zone, endoidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::addPathVesicle(std::string const& path_name,
                         std::string const& ves,
                         double speed,
                         const std::map<std::string, uint>& spec_deps,
                         const std::vector<double>& stoch_stepsize) {
    if (stoch_stepsize.size() > 2) {
        std::ostringstream os;
        os << "Stochastic step size must be a single value or list of length 2 (if applying a "
              "double exponential).";
        throw steps::ArgErr(os.str());
    }
    if (stoch_stepsize.size() == 2) {
        if (!(stoch_stepsize[1] > 0.0 && stoch_stepsize[1] < 1.0)) {
            std::ostringstream os;
            os << "Stochastic step size double exponential factor must be between 0 and 1.";
            throw steps::ArgErr(os.str());
        }
    }

    vesicle_global_id vidx = pStatedef->getVesicleIdx(ves);

    std::map<spec_global_id, uint> spec_deps_gidx;

    for (const auto& [spec, gid]: spec_deps) {
        spec_deps_gidx[pStatedef->getSpecIdx(spec)] = gid;
    }

    _addPathVesicle(path_name, vidx, speed, spec_deps_gidx, stoch_stepsize);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getTriRaftCount(triangle_global_id tidx, std::string const& r) const {
    if (auto mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        if (tidx >= mesh->countTris()) {
            std::ostringstream os;
            os << "Triangle index out of range.";
            throw ArgErr(os.str());
        }
        // the following may raise exceptions if strings are unused
        raft_global_id ridx = pStatedef->getRaftIdx(r);

        return _getTriRaftCount(tidx, ridx);
    }

    else {
        std::ostringstream os;
        os << "Method not available for this solver.";
        throw NotImplErr();
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriRaftCount(triangle_global_id tidx, std::string const& r, uint n) {
    if (auto mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        if (tidx >= mesh->countTris()) {
            std::ostringstream os;
            os << "Triangle index out of range.";
            throw ArgErr(os.str());
        }
        // the following may raise exceptions if strings are unused
        raft_global_id ridx = pStatedef->getRaftIdx(r);

        _setTriRaftCount(tidx, ridx, n);
    }

    else {
        std::ostringstream os;
        os << "Method not available for this solver.";
        throw NotImplErr();
    }
}

////////////////////////////////////////////////////////////////////////////////

raft_individual_id API::addTriRaft(triangle_global_id tidx, std::string const& r) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        if (tidx >= mesh->countTris()) {
            std::ostringstream os;
            os << "Triangle index out of range.";
            throw steps::ArgErr(os.str());
        }
        // the following may raise exceptions if strings are unused
        raft_global_id ridx = pStatedef->getRaftIdx(r);

        return _addTriRaft(tidx, ridx);
    }

    else {
        std::ostringstream os;
        os << "Method not available for this solver.";
        throw NotImplErr();
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriExocytosisActive(triangle_global_id tidx, std::string const& e) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        if (tidx >= mesh->countTris()) {
            std::ostringstream os;
            os << "Triangle index out of range.";
            throw steps::ArgErr(os.str());
        }
        // the following may raise exception if string is unknown
        exocytosis_global_id eidx = pStatedef->getExocytosisIdx(e);

        return _getTriExocytosisActive(tidx, eidx);
    }

    else {
        std::ostringstream os;
        os << "Method not available for this solver.";
        throw NotImplErr();
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriExocytosisActive(triangle_global_id tidx, std::string const& e, bool act) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        if (tidx >= mesh->countTris()) {
            std::ostringstream os;
            os << "Triangle index out of range.";
            throw steps::ArgErr(os.str());
        }

        // the following may raise exception if string is unknown
        exocytosis_global_id eidx = pStatedef->getExocytosisIdx(e);

        _setTriExocytosisActive(tidx, eidx, act);
    }

    else {
        std::ostringstream os;
        os << "Method not available for this solver.";
        throw NotImplErr();
    }
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getCompVesicleCount(comp_global_id /*cidx*/, vesicle_global_id /*vidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setCompVesicleCount(comp_global_id /*cidx*/, vesicle_global_id /*vidx*/, uint /*n*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

vesicle_individual_id API::_addCompVesicle(comp_global_id /*cidx*/, vesicle_global_id /*vidx*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_deleteSingleVesicle(vesicle_global_id /*vidx*/,
                               vesicle_individual_id /*ves_unique_index*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getSingleVesicleSurfaceLinkSpecCount(vesicle_global_id /*vidx*/,
                                                vesicle_individual_id /*ves_unique_index*/,
                                                linkspec_global_id /*lsidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<linkspec_individual_id> API::_getSingleVesicleSurfaceLinkSpecIndices(
    vesicle_global_id /*vidx*/,
    vesicle_individual_id /*ves_unique_index*/,
    linkspec_global_id /*lsidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<pointspec_individual_id> API::_getSingleVesicleSurfaceSpecIndices(
    vesicle_global_id /*vidx*/,
    vesicle_individual_id /*ves_unique_index*/,
    spec_global_id /*sidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<vesicle_individual_id> API::_getCompVesicleIndices(comp_global_id /*cidx*/,
                                                               vesicle_global_id /*vidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

comp_global_id API::_getSingleVesicleCompartment(vesicle_global_id /*vidx*/,
                                                 vesicle_individual_id /*ves_unique_index*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::_getSingleVesiclePos(vesicle_global_id /*vidx*/,
                                              vesicle_individual_id /*ves_unique_index*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setCompSingleVesiclePos(comp_global_id /*cidx*/,
                                   vesicle_global_id /*vidx*/,
                                   vesicle_individual_id /*ves_unique_index*/,
                                   const std::vector<double>& /*pos*/,
                                   bool /*force*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::map<vesicle_individual_id, uint> API::_getCompVesicleSurfaceSpecCountMap(
    comp_global_id /*cidx*/,
    vesicle_global_id /*vidx*/,
    spec_global_id /*sidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getCompVesicleSurfaceSpecCount(comp_global_id /*cidx*/,
                                          vesicle_global_id /*vidx*/,
                                          spec_global_id /*sidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getCompVesicleInnerSpecCount(comp_global_id /*cidx*/,
                                        vesicle_global_id /*vidx*/,
                                        spec_global_id /*sidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getSingleVesicleSurfaceSpecCount(vesicle_global_id /*vidx*/,
                                            vesicle_individual_id /*ves_unique_index*/,
                                            spec_global_id /*sidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getSingleVesicleInnerSpecCount(vesicle_global_id /*vidx*/,
                                          vesicle_individual_id /*ves_unique_index*/,
                                          spec_global_id /*sidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setSingleVesicleSurfaceSpecCount(vesicle_global_id /*vidx*/,
                                            vesicle_individual_id /*ves_unique_index*/,
                                            spec_global_id /*sidx*/,
                                            uint /*c*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> API::_getSingleVesicleSurfaceSpecPos(
    vesicle_global_id /*vidx*/,
    vesicle_individual_id /*ves_unique_index*/,
    spec_global_id /*sidx*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> API::_getSingleVesicleSurfaceSpecPosSpherical(
    vesicle_global_id /*vidx*/,
    vesicle_individual_id /*ves_unique_index*/,
    spec_global_id /*sidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setSingleVesicleSurfaceSpecPosSpherical(
    vesicle_global_id /*vidx*/,
    vesicle_individual_id /*ves_unique_index*/,
    spec_global_id /*sidx*/,
    const std::vector<std::vector<double>>& /*pos_spherical*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setSingleVesicleInnerSpecCount(vesicle_global_id /*vidx*/,
                                          vesicle_individual_id /*ves_unique_index*/,
                                          spec_global_id /*sidx*/,
                                          uint /*c*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::_getSingleSpecPosSpherical(
    spec_global_id /*sidx*/,
    pointspec_individual_id /*ps_unique_id*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::map<vesicle_individual_id, uint> API::_getCompVesicleSurfaceLinkSpecCountMap(
    comp_global_id /*cidx*/,
    vesicle_global_id /*vidx*/,
    linkspec_global_id /*lsidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getCompVesicleSurfaceLinkSpecCount(comp_global_id /*cidx*/,
                                              vesicle_global_id /*vidx*/,
                                              linkspec_global_id /*lsidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> API::_getSingleVesicleSurfaceLinkSpecPos(
    vesicle_global_id /*vidx*/,
    vesicle_individual_id /*ves_unique_index*/,
    linkspec_global_id /*lsidx*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::_getSingleLinkSpecPos(linkspec_individual_id /*ls_unique_id*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

linkspec_individual_id API::_getSingleLinkSpecLinkedTo(
    linkspec_individual_id /*ls_unique_id*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

vesicle_individual_id API::_getSingleLinkSpecVes(linkspec_individual_id /*ls_unique_id*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getSingleVesicleImmobility(vesicle_global_id /*vidx*/,
                                      vesicle_individual_id /*ves_unique_index*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::tetrahedron_global_id> API::_getSingleVesicleOverlapTets(
    vesicle_global_id /*vidx*/,
    vesicle_individual_id /*ves_unique_index*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetVesicleDcst(tetrahedron_global_id /*tidx*/,
                             vesicle_global_id /*vidx*/,
                             double /*dcst*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setVesicleSurfaceLinkSpecSDiffD(vesicle_global_id /*vidx*/,
                                           linkspec_global_id /*lsidx*/,
                                           double /*dcst*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setVesSReacK(vessreac_global_id /*vsridx*/, double /*kf*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getVesSReacExtent(vessreac_global_id /*vsridx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setExocytosisK(exocytosis_global_id /*exoidx*/, double /*kf*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getExocytosisExtent(exocytosis_global_id /*exoidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<ExocytosisEvent> API::_getExocytosisEvents(solver::exocytosis_global_id /*exoidx*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getRaftEndocytosisExtent(raftendocytosis_global_id /*rendoidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<RaftEndocytosisEvent> API::_getRaftEndocytosisEvents(
    raftendocytosis_global_id /*rendoidx*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setRaftEndocytosisK(raftendocytosis_global_id /*rendoidx*/, double /*kcst*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_addVesicleDiffusionGroup(vesicle_global_id /*vidx*/,
                                    const std::vector<comp_global_id>& /*comp_indices*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_addPathVesicle(std::string const& /*path_name*/,
                          vesicle_global_id /*vidx*/,
                          double /*speed*/,
                          const std::map<spec_global_id, uint>& /*spec_deps*/,
                          const std::vector<double>& /*stoch_stepsize*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getPatchRaftCount(patch_global_id /*pidx*/, raft_global_id /*ridx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setPatchRaftCount(patch_global_id /*pidx*/, raft_global_id /*ridx*/, uint /*n*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::_getSingleRaftPos(raft_global_id /*ridx*/,
                                           raft_individual_id /*raft_unique_index*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::map<raft_individual_id, uint> API::_getPatchRaftSpecCountMap(patch_global_id /*pidx*/,
                                                                  raft_global_id /*ridx*/,
                                                                  spec_global_id /*sidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getPatchRaftSpecCount(patch_global_id /*pidx*/,
                                 raft_global_id /*ridx*/,
                                 spec_global_id /*sidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getSingleRaftSpecCount(raft_global_id /*ridx*/,
                                  raft_individual_id /*raft_unique_index*/,
                                  spec_global_id /*sidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setSingleRaftSpecCount(raft_global_id /*ridx*/,
                                  raft_individual_id /*raft_unique_index*/,
                                  spec_global_id /*sidx*/,
                                  uint /*c*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getSingleRaftImmobility(raft_global_id /*ridx*/,
                                   raft_individual_id /*raft_unique_index*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<raft_individual_id> API::_getPatchRaftIndices(patch_global_id /*pidx*/,
                                                          raft_global_id /*ridx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

patch_global_id API::_getSingleRaftPatch(raft_global_id /*ridx*/,
                                         raft_individual_id /*raft_unique_index*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setPatchEndocyticZoneEndocytosisActive(patch_global_id /*pidx*/,
                                                  std::string const& /*zone*/,
                                                  endocytosis_global_id /*endogidx*/,
                                                  bool /*active*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setPatchEndocyticZoneEndocytosisK(patch_global_id /*pidx*/,
                                             std::string const& /*zone*/,
                                             endocytosis_global_id /*endogidx*/,
                                             double /*k*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getPatchEndocyticZoneEndocytosisExtent(patch_global_id /*pidx*/,
                                                  std::string const& /*zone*/,
                                                  endocytosis_global_id /*endogidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

std::vector<EndocytosisEvent> API::_getPatchEndocyticZoneEndocytosisEvents(
    patch_global_id /*pidx*/,
    std::string const& /*zone*/,
    endocytosis_global_id /*endogidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getSingleRaftRaftEndocytosisK(raft_global_id /*ridx*/,
                                           raft_individual_id /*raft_unique_index*/,
                                           raftendocytosis_global_id /*rendoidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setSingleRaftRaftEndocytosisK(raft_global_id /*ridx*/,
                                         raft_individual_id /*raft_unique_index*/,
                                         raftendocytosis_global_id /*rendoidx*/,
                                         double /*k*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setSingleRaftSReacActive(raft_global_id /*ridx*/,
                                    raft_individual_id /*raft_unique_index*/,
                                    raftsreac_global_id /*rsreacidx*/,
                                    bool /*active*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getSingleRaftSReacActive(raft_global_id /*ridx*/,
                                    raft_individual_id /*raft_unique_index*/,
                                    raftsreac_global_id /*rsreacidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getTriRaftCount(triangle_global_id /*tidx*/, raft_global_id /*ridx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriRaftCount(triangle_global_id /*tidx*/, raft_global_id /*ridx*/, uint /*n*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

raft_individual_id API::_addTriRaft(triangle_global_id /*tidx*/, raft_global_id /*ridx*/) {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriExocytosisActive(triangle_global_id /*tidx*/,
                                  exocytosis_global_id /*eidx*/) const {
    throw NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriExocytosisActive(triangle_global_id /*tidx*/,
                                  exocytosis_global_id /*eidx*/,
                                  bool /*act*/) {
    throw NotImplErr();
}

}  // namespace steps::solver
