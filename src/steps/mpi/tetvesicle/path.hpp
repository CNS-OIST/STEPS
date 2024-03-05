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

#pragma once

// Standard library & STL headers.
//#include <string>
#include <math.h>
#include <random>
#include <vector>

// STEPS headers.
#include "math/point.hpp"
#include "mpi/tetvesicle/vesicle.hpp"

namespace steps::mpi::tetvesicle {

////////////////////////////////////////////////////////////////////////////////

class Path {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Path(std::string const& id);

    ~Path();

    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    void addVesicle(solver::vesicle_global_id ves_idx,
                    double speed,
                    const std::map<solver::spec_global_id, uint>& spec_deps,
                    const std::vector<double>& stoch_stepsize);

    void addPoint(uint index, const std::array<double, 3>& position);

    void addBranch(uint root_point, std::map<uint, double> dest_points);

    std::map<uint, std::pair<std::vector<double>, std::map<uint, double>>> getPathMap() const;

    bool crossedPath(math::position_abs ves_pos,
                     solver::vesicle_global_id ves_gidx,
                     double ves_rad);

    /// Return the patch id.
    ///
    /// \return ID of the patch.
    inline std::string getID() const noexcept {
        return pID;
    }

    // calculate a route through the path based on the vesicle speed, the vesicle
    // update dt and a random choice of branches
    // Keep this quite dynamic in order to allow for adding/removing/growing Paths
    // in the simulation in the future The vesicle radius is there really just for
    // some checks.
    std::vector<std::pair<double, math::position_abs>> calculateRoute(
        math::position_abs starting_pos,
        solver::vesicle_global_id vesicle_index,
        const rng::RNGptr& rng,
        double vesicle_radius) const;

    std::map<solver::spec_global_id, uint> const& getVesicleSpecDeps(
        solver::vesicle_global_id ves_idx) const;

    ////////////////////////////////////////////////////////////////////////

  private:
    std::string pID;

    // list of vesicles that belong to this path by index, along with their speeds
    std::map<solver::vesicle_global_id, double> pVesicles;
    std::map<solver::vesicle_global_id, std::map<solver::spec_global_id, uint>> pVesicle_specdeps;
    std::map<solver::vesicle_global_id, std::vector<double>> pVesicle_stochsteps;

    std::map<uint, math::position_abs> pPoints;

    std::map<uint, std::map<uint, double>> pBranches;
};

}  // namespace steps::mpi::tetvesicle
