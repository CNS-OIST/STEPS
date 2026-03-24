/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2026 Okinawa Institute of Science and Technology, Japan.
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

#include "mpi/tetvesicle/path.hpp"

// STEPS headers.
#include "math/constants.hpp"  //PI
#include "math/point.hpp"
#include "solver/fwd.hpp"
#include "util/checkpointing.hpp"
#include "util/error.hpp"
#include <algorithm>
#include <limits>
#include <numeric>

namespace steps::mpi::tetvesicle {

////////////////////////////////////////////////////////////////////////////////

void PathVesicleInfo::checkpoint(std::ostream& cp_file) const {
    util::checkpoint(cp_file, speed);
    util::checkpoint(cp_file, dependencies);
    util::checkpoint(cp_file, stoch_steps);
    util::checkpoint(cp_file, binding_rate);
    util::checkpoint(cp_file, unbinding_rate);
    util::checkpoint(cp_file, allow_path_intersection);
    util::checkpoint(cp_file, _min_binding_radius);
    util::checkpoint(cp_file, _max_binding_radius);
}

////////////////////////////////////////////////////////////////////////////////

void PathVesicleInfo::restore(std::istream& cp_file) {
    util::restore(cp_file, speed);
    util::restore(cp_file, dependencies);
    util::restore(cp_file, stoch_steps);
    util::restore(cp_file, binding_rate);
    util::restore(cp_file, unbinding_rate);
    util::restore(cp_file, allow_path_intersection);
    util::restore(cp_file, _min_binding_radius);
    util::restore(cp_file, _max_binding_radius);
}

////////////////////////////////////////////////////////////////////////////////

double PathVesicleInfo::min_binding_radius(const Vesicle& ves) const {
    if (_min_binding_radius < 0) {
        return -_min_binding_radius * ves.getDiam() / 2.0;
    } else {
        return _min_binding_radius;
    }
}

////////////////////////////////////////////////////////////////////////////////

double PathVesicleInfo::max_binding_radius(const Vesicle& ves) const {
    if (_max_binding_radius < 0) {
        return -_max_binding_radius * ves.getDiam() / 2.0;
    } else {
        return _max_binding_radius;
    }
}

////////////////////////////////////////////////////////////////////////////////

void PathEdge::checkpoint(std::ostream& cp_file) const {
    util::checkpoint(cp_file, source);
    util::checkpoint(cp_file, destination);
    util::checkpoint(cp_file, allow_binding);
    util::checkpoint(cp_file, probability);
}

////////////////////////////////////////////////////////////////////////////////

void PathEdge::restore(std::istream& cp_file) {
    util::restore(cp_file, source);
    util::restore(cp_file, destination);
    util::restore(cp_file, allow_binding);
    util::restore(cp_file, probability);
}

////////////////////////////////////////////////////////////////////////////////

Path::Path(solver::path_global_id const& id, const std::string& name, bool bind_to_start)
    : pID(id)
    , pName(name)
    , pBind_to_start(bind_to_start) {}

////////////////////////////////////////////////////////////////////////////////

/**
 * \TODO(Weiliang) Question of whether Paths can be removed during a simulation to mimic
 * dynamic chenges in actin etc. And if so the cleanup that will be necessary.
 * Minimally any tethered vesicles will have to disperse and solver will need
 * to remove this Path from its memory banks.
 */
Path::~Path() = default;

////////////////////////////////////////////////////////////////////////////////

void Path::checkpoint(std::ostream& cp_file) const {
    util::checkpoint(cp_file, pID);
    util::checkpoint(cp_file, pName);
    util::checkpoint(cp_file, pBind_to_start);
    util::checkpoint(cp_file, pVesicles);
    util::checkpoint(cp_file, pPoints);
    util::checkpoint(cp_file, pEdges);
    util::checkpoint(cp_file, pEdgeMap);
}

////////////////////////////////////////////////////////////////////////////////

void Path::restore(std::istream& cp_file) {
    util::restore(cp_file, pID);
    util::restore(cp_file, pName);
    util::restore(cp_file, pBind_to_start);
    util::restore(cp_file, pVesicles);
    util::restore(cp_file, pPoints);
    util::restore(cp_file, pEdges);
    util::restore(cp_file, pEdgeMap);

    // Rebuild the spatial edge index
    edgeIndex.clear();
    for (auto eid: pEdges.range()) {
        const auto& edge = pEdges[eid];
        if (edge.allow_binding) {
            segment seg(pPoints.at(edge.source), pPoints.at(edge.destination));
            edgeIndex.insert(std::make_pair(seg, eid));
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void Path::addVesicle(solver::vesicle_global_id ves_idx,
                      double speed,
                      const std::map<solver::spec_global_id, uint>& spec_deps,
                      const std::vector<double>& stoch_stepsize,
                      double binding_rate,
                      double min_binding_radius,
                      double max_binding_radius,
                      double unbinding_rate,
                      bool allow_path_intersection) {
    auto [_, inserted] = pVesicles.emplace(ves_idx,
                                           PathVesicleInfo(speed,
                                                           spec_deps,
                                                           stoch_stepsize,
                                                           binding_rate,
                                                           min_binding_radius,
                                                           max_binding_radius,
                                                           unbinding_rate,
                                                           allow_path_intersection));
    ArgErrLogIf(not inserted, "Vesicle has already been added to Path.");
}

////////////////////////////////////////////////////////////////////////////////

void Path::addPoint(uint index, const std::array<double, 3>& position) {
    if (pPoints.find(index) == pPoints.end()) {
        math::position_abs pos{position[0], position[1], position[2]};
        pPoints[index] = pos;
    } else {
        ArgErrLog("Point has already been added to Path.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void Path::addEdge(uint source, uint dest, double weight, bool allow_binding) {
    if (pBind_to_start and pPoints.find(1) == pPoints.end()) {
        ArgErrLog(
            "Paths that use bind_to_start must contain root Point (index == 1) before edges can "
            "be added.");
    }

    auto sourceIt = pPoints.find(source);
    if (sourceIt == pPoints.end()) {
        std::ostringstream os;
        os << "Point " << source << " is unknown in Path.";
        ArgErrLog(os.str());
    }
    auto destIt = pPoints.find(dest);
    if (destIt == pPoints.end()) {
        std::ostringstream os;
        os << "Point " << dest << " is unknown in Path.";
        ArgErrLog(os.str());
    }
    if (weight < 0.0) {
        std::ostringstream os;
        os << "Negative weights are not allowed.";
        ArgErrLog(os.str());
    }

    auto& edgesIds = pEdgeMap[source];
    solver::path_edge_global_id edgeId{pEdges.size()};
    pEdges.container().emplace_back(source, dest, weight, allow_binding);
    edgesIds.emplace_back(edgeId);

    // Recompute probabilities based on weights
    double sum = std::transform_reduce(
        edgesIds.begin(), edgesIds.end(), 0.0, std::plus{}, [this](const auto& edgeId) {
            return pEdges[edgeId].weight;
        });
    for (auto& eid: edgesIds) {
        auto& edge = pEdges[eid];
        edge.probability = edge.weight / sum;
    }

    // Insert the corresponding segment in the spatial index
    segment seg(sourceIt->second, destIt->second);
    edgeIndex.insert(std::make_pair(seg, edgeId));
}

////////////////////////////////////////////////////////////////////////////////

std::map<uint, std::pair<std::vector<double>, std::map<uint, double>>> Path::getPathMap() const {
    std::map<uint, std::pair<std::vector<double>, std::map<uint, double>>> ret;
    for (const auto& p: pPoints) {
        auto ip = ret.emplace(p.first,
                              std::pair<std::vector<double>, std::map<uint, double>>(
                                  std::vector<double>(p.second.begin(), p.second.end()), {}));
        const auto edgesIdsIt = pEdgeMap.find(p.first);
        if (ip.second and edgesIdsIt != pEdgeMap.end()) {
            for (const auto& edgeId: edgesIdsIt->second) {
                auto& edge = pEdges[edgeId];
                ip.first->second.second.emplace(edge.destination, edge.probability);
            }
        }
    }
    return ret;
}

////////////////////////////////////////////////////////////////////////////////

solver::path_edge_global_id Path::selectEdge(uint source, const rng::RNGptr& rng) const {
    double accum = 0.0;
    double rn = rng->getUnfII();

    auto edgesIds = pEdgeMap.find(source);

    if (edgesIds != pEdgeMap.end()) {
        for (auto edgeId: edgesIds->second) {
            auto& edge = pEdges[edgeId];
            accum += edge.probability;
            if (accum >= rn) {
                return edgeId;
            }
        }
    }
    return {};
}

////////////////////////////////////////////////////////////////////////////////

std::pair<std::vector<std::pair<double, math::position_abs>>, math::position_abs>
Path::calculateRoute(const Vesicle& ves, const rng::RNGptr& rng) const {
    auto vesicle_index = ves.idx();

    auto it = pVesicles.find(vesicle_index);
    if (it == pVesicles.end()) {
        ProgErrLog("Vesicle has not been added to Path.");
    }
    const auto& info = it->second;

    if (pEdgeMap.size() == 0) {
        std::ostringstream os;
        os << "\nPath " << pID << " does not contain any Edges.";
        ProgErrLog(os.str());
    }

    solver::path_edge_global_id current_edge;
    double current_frac;
    math::position_abs pos_track{ves.getPosition()};

    // First attempt binding to the start, if possible
    if (pBind_to_start) {
        const auto dist = math::distance(ves.getPosition(), pPoints.at(1));
        if (dist <= info.max_binding_radius(ves) and dist >= info.min_binding_radius(ves)) {
            current_edge = selectEdge(1, rng);
            current_frac = 0;
        }
    }
    // If we did not bind to the start, attempt binding to path edges
    if (current_edge.unknown()) {
        const auto& vespos = ves.getPosition();
        math::position_abs shift{info.max_binding_radius(ves),
                                 info.max_binding_radius(ves),
                                 info.max_binding_radius(ves)};
        boost::geometry::model::box<math::position_abs> ves_box(vespos - shift, vespos + shift);
        double minDist = std::numeric_limits<double>::infinity();
        for (auto rit = edgeIndex.qbegin(boost::geometry::index::intersects(ves_box));
             rit != edgeIndex.qend();
             ++rit) {
            // Check whether the vesicle really intersects the path edge
            const auto& seg = rit->first;
            const auto& edge = pEdges[rit->second];
            double dist = boost::geometry::distance(vespos, seg);
            if (dist <= info.max_binding_radius(ves) and dist >= info.min_binding_radius(ves) and
                dist < minDist) {
                // Find closest point on segment
                const auto& source = pPoints.at(edge.source);
                const auto& dest = pPoints.at(edge.destination);
                minDist = dist;

                current_edge = rit->second;
                current_frac =
                    std::clamp((dest - source).dot(vespos - source) /
                                   (math::norm(dest - source) * math::norm(dest - source)),
                               0.0,
                               1.0);
            }
        }
    }
    ArgErrLogIf(current_edge.unknown(), "Could not bind to path.");

    std::vector<std::pair<double, math::position_abs>> positions;  // the delta time to next
                                                                   // position, and the position

    // Need to record the starting position because it's perfectly possible the vesicle
    // won't move on each vesicle_dt
    positions.emplace_back(0.0, pos_track);

    double speed = info.speed;
    const auto& stoch_stepsize = info.stoch_steps;
    double step_size = stoch_stepsize[0];
    double doubleexp_factor = stoch_stepsize.size() == 2 ? stoch_stepsize[1] : -1;
    std::optional<math::position_abs> starting_shift;

    do {
        const auto& edge = pEdges[current_edge];
        auto psrc = pPoints.at(edge.source);
        auto pdst = pPoints.at(edge.destination);
        double dist = math::distance(pdst, psrc);
        double scaled_dist = dist * (1.0 - current_frac);

        if (not starting_shift) {
            starting_shift = psrc + current_frac * (pdst - psrc) - ves.getPosition();
        }

        if (scaled_dist > 0) {
            auto unit_vect = (pdst - psrc) / dist;

            // Vesicle moves in pre-determined step lengths. Pre-calculate stochastic steps
            // depending on vesicle dt and speed. Do this by calculating the arrival times for each
            // step, and stopping to add to positions whenever we go past the next vesicle time
            // point

            uint n_steps_reg = static_cast<unsigned int>(
                floor(scaled_dist / step_size));  // the number of 'regular' steps of stoch_stepsize
            double dist_last = scaled_dist - (n_steps_reg * step_size);

            double delta_next;  // Time delta to the next step (of the stoch_stepsize or to end
                                // point of branch if that is shorter)
            double dist_next = step_size;  // Distance of the next step (usually equal to
                                           // stoch_stepsize)

            for (uint n = 0; n < n_steps_reg + 1; ++n) {
                if (n == n_steps_reg) {
                    if (dist_last <= std::numeric_limits<double>::epsilon()) {
                        break;
                    }
                    dist_next = dist_last;
                }

                if (doubleexp_factor > 0.0) {
                    delta_next = rng->getExp((1.0 / doubleexp_factor) * (speed / dist_next)) +
                                 rng->getExp((1.0 / (1.0 - doubleexp_factor)) *
                                             (speed / dist_next));
                } else {
                    delta_next = rng->getExp(1.0 / (dist_next / speed));
                }

                pos_track += unit_vect * dist_next;

                positions.emplace_back(delta_next, pos_track);
            }
        }
        AssertLog(math::distance(pos_track, pdst) <= info.max_binding_radius(ves) and
                  math::distance(pos_track, pdst) >= info.min_binding_radius(ves));

        current_edge = selectEdge(edge.destination, rng);
        current_frac = 0.0;
    } while (current_edge.valid());

    return {positions, *starting_shift};
}

////////////////////////////////////////////////////////////////////////////////

double Path::getBindingRate(const Vesicle& ves) const {
    // Check that the vesicle can bind to the path
    auto it = pVesicles.find(ves.idx());
    if (it == pVesicles.end()) {
        return 0;
    }
    const auto& pvInfo = it->second;

    // Check whether the vesicle has the species necessary to bind to the path
    for (auto const& vsd: pvInfo.dependencies) {
        if (ves.getSurfSpecCount(vsd.first) < vsd.second) {
            return 0;
        }
    }

    if (pBind_to_start) {
        const auto dist = math::distance(ves.getPosition(), pPoints.at(1));
        if (dist <= pvInfo.max_binding_radius(ves) and dist >= pvInfo.min_binding_radius(ves)) {
            return 1e10;  // High value for instantaneous binding
        }
    }

    // Check that the binding zone intersects with the path, and if the path does not allow
    // intersection of the vesicle, also check that the vesicle is not intersecting any part of the
    // path.
    if (intersectsSphere(ves.getPosition(),
                         pvInfo.min_binding_radius(ves),
                         pvInfo.max_binding_radius(ves),
                         true) and
        (pvInfo.allow_path_intersection or
         not intersectsSphere(ves.getPosition(), 0, ves.getDiam() / 2.0, false))) {
        return pvInfo.binding_rate;
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

bool Path::intersectsSphere(const math::position_abs& center,
                            const double& min_radius,
                            const double& max_radius,
                            bool binding_edges_only) const {
    // Find intersections with path edge bounding boxes
    math::position_abs shift{max_radius, max_radius, max_radius};
    boost::geometry::model::box<math::position_abs> ves_box(center - shift, center + shift);
    for (auto rit = edgeIndex.qbegin(boost::geometry::index::intersects(ves_box));
         rit != edgeIndex.qend();
         ++rit) {
        // Check whether the vesicle really intersects the path edge
        const auto& seg = rit->first;
        if (binding_edges_only and not pEdges[rit->second].allow_binding) {
            continue;
        }
        double dist = boost::geometry::distance(center, seg);
        if (dist <= max_radius and dist >= min_radius) {
            return true;
        }
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool Path::canIntersect(const Vesicle& ves) const {
    auto it = pVesicles.find(ves.idx());
    if (it == pVesicles.end()) {
        return true;
    }
    return it->second.allow_path_intersection;
}

////////////////////////////////////////////////////////////////////////////////

double Path::getUnBindingRate(const Vesicle& ves) const {
    // Check that the vesicle can bind to the path
    auto it = pVesicles.find(ves.idx());
    if (it == pVesicles.end()) {
        return 0;
    }
    return it->second.unbinding_rate;
}

}  // namespace steps::mpi::tetvesicle
