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

#include "mpi/tetvesicle/path.hpp"

// STEPS headers.
#include "math/constants.hpp"  //PI
#include "util/checkpointing.hpp"

namespace steps::mpi::tetvesicle {

////////////////////////////////////////////////////////////////////////////////

// Path::Path(void)
Path::Path(std::string const& id)
    : pID(id) {}

////////////////////////////////////////////////////////////////////////////////

/**
 * \TODO(Weiliang) Question of whether Paths can be removed during a simulation to mimic
 * dynamic chenges in actin etc. And if so the cleanup that will be necessary.
 * Minimally any tethered vesicles will have to disperse and solver will need
 * to remove this Path from its memory banks.
 */
Path::~Path() = default;

////////////////////////////////////////////////////////////////////////////////

void Path::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pVesicles);
    util::checkpoint(cp_file, pVesicle_specdeps);
    util::checkpoint(cp_file, pPoints);
    util::checkpoint(cp_file, pBranches);
    util::checkpoint(cp_file, pVesicle_stochsteps);
}

////////////////////////////////////////////////////////////////////////////////

void Path::restore(std::fstream& cp_file) {
    util::restore(cp_file, pVesicles);
    util::restore(cp_file, pVesicle_specdeps);
    util::restore(cp_file, pPoints);
    util::restore(cp_file, pBranches);
    util::restore(cp_file, pVesicle_stochsteps);
}

////////////////////////////////////////////////////////////////////////////////

void Path::addVesicle(solver::vesicle_global_id ves_idx,
                      double speed,
                      const std::map<solver::spec_global_id, uint>& spec_deps,
                      const std::vector<double>& stoch_stepsize) {
    if (pVesicles.find(ves_idx) == pVesicles.end()) {
        pVesicles[ves_idx] = speed;
        AssertLog(pVesicle_specdeps.find(ves_idx) == pVesicle_specdeps.end());
        pVesicle_specdeps[ves_idx] = spec_deps;
        AssertLog(pVesicle_stochsteps.find(ves_idx) == pVesicle_stochsteps.end());
        pVesicle_stochsteps[ves_idx] = stoch_stepsize;
    } else {
        ArgErrLog("Vesicle has already been added to Path.");
    }
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

void Path::addBranch(uint root_point, std::map<uint, double> term_points) {
    // Simplest check first
    if (pPoints.find(1) == pPoints.end()) {
        ArgErrLog(
            "Path must contain root Point (index == 1) before branches can "
            "be added.");
    }

    if (pPoints.find(root_point) == pPoints.end()) {
        std::ostringstream os;
        os << "Point " << root_point << " is unknown in Path.";
        ArgErrLog(os.str());
    }

    if (pBranches.find(root_point) != pBranches.end()) {
        std::ostringstream os;
        os << "Point " << root_point
           << " already has specified branching (this can only be done once per "
              "point).";
        ArgErrLog(os.str());
    }

    double total_weight = 0.0;
    for (auto t: term_points) {
        if (pPoints.find(t.first) == pPoints.end()) {
            std::ostringstream os;
            os << "Point " << t.first << " is unknown in Path.";
            ArgErrLog(os.str());
        }

        if (t.first == root_point) {
            std::ostringstream os;
            os << "Point cannot branch to itself.";
            ArgErrLog(os.str());
        }

        // Check the weights
        if (t.second < 0.0) {
            std::ostringstream os;
            os << "Negative weights are not allowed.";
            ArgErrLog(os.str());
        }

        total_weight += t.second;
    }

    if (total_weight == 0.0) {
        std::ostringstream os;
        os << "Total weight can't be zero.";
        ArgErrLog(os.str());
    }

    // Let's normalise
    for (auto t: term_points) {
        // doesn't work t.second/= total_weight;
        term_points[t.first] = t.second / total_weight;
    }

    // So, simple checks that the points are defined has passed. Now check the
    // branching. We start at 1 and allow divergence, not convergence. So e.g. if
    // we have a branch from 1 to 2 and 2 to 3, we cannot have a branch from 3
    // to 2. This means that Points must be a terminal in a branch before they can
    // become a root (with the exception of point1), and a point can only be a
    // termninal once.

    bool gotterm = false;

    if (root_point == 1) {
        gotterm = true;
    }

    for (auto b: pBranches) {
        // First check if root_point appears as a terminal somewhere
        for (auto r: b.second) {
            // r.first is the index of the point
            if (r.first == root_point) {
                // This should only happen once so let's check that with an AssertLog.
                AssertLog(gotterm == false);
                gotterm = true;
            }

            for (auto t: term_points) {
                if (t.first == r.first) {
                    std::ostringstream os;
                    os << "Point " << t.first
                       << " is already a terminal, so cannot be terminal in another "
                          "branch. ";
                    ArgErrLog(os.str());
                }
            }
        }

        // Now check the roots in other branches. A root cannot be a terminal here
        for (auto t: term_points) {
            if (t.first == b.first) {
                std::ostringstream os;
                os << "Point " << t.first
                   << " is a root in another branch, so cannot be a terminal here. ";
                ArgErrLog(os.str());
            }
        }
    }

    if (not gotterm) {
        std::ostringstream os;
        os << "Point " << root_point << " is not a terminal point yet, so cannot add this Branch.";
        ArgErrLog(os.str());
    }

    pBranches[root_point] = term_points;
}

////////////////////////////////////////////////////////////////////////////////

std::map<uint, std::pair<std::vector<double>, std::map<uint, double>>> Path::getPathMap() const {
    std::map<uint, std::pair<std::vector<double>, std::map<uint, double>>> ret;
    for (const auto& p: pPoints) {
        auto ip = ret.emplace(p.first,
                              std::pair<std::vector<double>, std::map<uint, double>>(
                                  std::vector<double>(p.second.begin(), p.second.end()), {}));
        const auto branches = pBranches.find(p.first);
        if (ip.second and branches != pBranches.end()) {
            ip.first->second.second.insert(branches->second.begin(), branches->second.end());
        }
    }
    return ret;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::pair<double, math::position_abs>> Path::calculateRoute(
    math::position_abs starting_pos,
    solver::vesicle_global_id vesicle_index,
    const rng::RNGptr& rng,
    double vesicle_radius) const {
    if (pVesicles.find(vesicle_index) == pVesicles.end()) {
        ProgErrLog("Vesicle has not been added to Path.");
    }

    if (math::distance(starting_pos, pPoints.at(1)) > vesicle_radius) {
        ProgErrLog(
            "Internal error: vesicle position does not overlap starting "
            "point of Path. Contact STEPS support!");
    }

    if (pBranches.size() == 0) {
        std::ostringstream os;
        os << "\nPath " << pID << " does not contain any Branches.";
        ProgErrLog(os.str());
    }

    std::vector<std::pair<double, math::position_abs>> positions;  // the delta time to next
                                                                   // position, and the position

    uint current_point = 1;
    math::position_abs pos_track = starting_pos;

    // Need to record the starting position because it's perfectly possible the vesicle
    // won't move on each vesicle_dt
    positions.emplace_back(std::pair<double, math::position_abs>(0.0, pos_track));

    // While we still have branching
    auto branch = pBranches.find(current_point);
    while (branch != pBranches.end()) {
        std::map<uint, double> destinations(branch->second);

        // Pick the destination
        double accum = 0.0;
        double rn = rng->getUnfII();

        uint destination = UINT_MAX;

        for (auto d: destinations) {
            accum += d.second;
            if (accum >= rn) {
                destination = d.first;
                break;
            }
        }

        AssertLog(destination != UINT_MAX);

        math::position_abs pos_root = pPoints.at(current_point);
        math::position_abs pos_dest = pPoints.at(destination);

        double speed = pVesicles.at(vesicle_index);
        double dist_total = math::distance(pos_root, pos_dest);

        const auto& stoch_stepsize = pVesicle_stochsteps.at(vesicle_index);

        // Vesicle moves in pre-determined step lengths. Pre-calculate stochastic steps
        // depending on vesicle dt and speed. Do this by calculating the arrival times for each
        // step, and stopping to add to positions whenever we go past the next vesicle time
        // point

        double step_size = stoch_stepsize[0];
        double doubleexp_factor = stoch_stepsize.size() == 2 ? stoch_stepsize[1] : -1;

        math::position_abs path_vector_unitary = (pos_dest - pos_root) /
                                                 math::distance(pos_root, pos_dest);

        uint n_steps_reg = static_cast<unsigned int>(
            floor(dist_total / step_size));  // the number of 'regular' steps of stoch_stepsize
        double dist_last = dist_total - (n_steps_reg * step_size);

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
                             rng->getExp((1.0 / (1.0 - doubleexp_factor)) * (speed / dist_next));
            } else {
                delta_next = rng->getExp(1.0 / (dist_next / speed));
            }

            pos_track += path_vector_unitary * dist_next;

            positions.emplace_back(delta_next, pos_track);
        }


        if (math::distance(pos_track, pos_dest) > vesicle_radius) {
            ProgErrLog(
                "Internal error: vesicle position does not overlap point of "
                "Path. Contact STEPS support!");
        }

        current_point = destination;
        branch = pBranches.find(current_point);
    }

    return positions;
}

////////////////////////////////////////////////////////////////////////////////

bool Path::crossedPath(math::position_abs ves_pos,
                       solver::vesicle_global_id ves_gidx,
                       double ves_rad) {
    if (pVesicles.find(ves_gidx) == pVesicles.end()) {
        return false;
    }

    if (math::distance(ves_pos, pPoints[1]) < ves_rad) {
        return true;
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////

std::map<solver::spec_global_id, uint> const& Path::getVesicleSpecDeps(
    solver::vesicle_global_id ves_idx) const {
    auto it = pVesicle_specdeps.find(ves_idx);
    if (it == pVesicle_specdeps.end()) {
        ProgErrLog(
            "Internal error: Vesicle has not been added to Path. Contact "
            "STEPS support!!");
    }

    return it->second;
}

}  // namespace steps::mpi::tetvesicle
