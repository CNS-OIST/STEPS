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

#pragma once

// Standard library & STL headers.
#include <math.h>
#include <random>
#include <vector>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/index/rtree.hpp>

// STEPS headers.
#include "math/point.hpp"
#include "mpi/tetvesicle/vesicle.hpp"
#include "solver/fwd.hpp"

// Register math::absolute_position as a point type for boost::geometry
namespace boost::geometry::traits {

template <>
struct tag<steps::math::position_abs> {
    using type = point_tag;
};
template <>
struct dimension<steps::math::position_abs>: boost::mpl::int_<3> {};
template <>
struct coordinate_type<steps::math::position_abs> {
    using type = double;
};
template <>
struct coordinate_system<steps::math::position_abs> {
    using type = boost::geometry::cs::cartesian;
};
template <std::size_t Index>
struct access<steps::math::position_abs, Index> {
    static_assert(Index < 3, "Out of range");
    using Point = steps::math::position_abs;
    using CoordinateType = typename coordinate_type<Point>::type;
    static inline CoordinateType get(Point const& p) {
        return p[Index];
    }
    static inline void set(Point& p, CoordinateType const& value) {
        p[Index] = value;
    }
};

}  // namespace boost::geometry::traits

namespace steps::mpi::tetvesicle {

////////////////////////////////////////////////////////////////////////////////

struct PathVesicleInfo {
    PathVesicleInfo() = default;
    PathVesicleInfo(double _s,
                    const std::map<solver::spec_global_id, uint>& _d,
                    const std::vector<double>& _ss,
                    double _brat,
                    double _brad1,
                    double _brad2,
                    double _ubr,
                    bool _api)
        : speed(_s)
        , dependencies(_d)
        , stoch_steps(_ss)
        , binding_rate(_brat)
        , unbinding_rate(_ubr)
        , allow_path_intersection(_api)
        , _min_binding_radius(_brad1)
        , _max_binding_radius(_brad2) {}

    double speed{};
    std::map<solver::spec_global_id, uint> dependencies;
    std::vector<double> stoch_steps;
    double binding_rate{};
    double unbinding_rate{};
    bool allow_path_intersection{true};

    void checkpoint(std::ostream& cp_file) const;
    void restore(std::istream& cp_file);

    double min_binding_radius(const Vesicle& ves) const;
    double max_binding_radius(const Vesicle& ves) const;

  private:
    double _min_binding_radius{};
    double _max_binding_radius{};
};

////////////////////////////////////////////////////////////////////////////////

struct PathEdge {
    PathEdge() = default;
    PathEdge(uint _s, uint _d, double _w, bool _a)
        : source(_s)
        , destination(_d)
        , weight(_w)
        , allow_binding(_a)
        , probability(1) {}

    uint source;
    uint destination;
    double weight;
    bool allow_binding;
    double probability;

    void checkpoint(std::ostream& cp_file) const;
    void restore(std::istream& cp_file);
};

////////////////////////////////////////////////////////////////////////////////

class Path {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Path() = default;
    Path(solver::path_global_id const& id, const std::string& name, bool bind_to_start);

    ~Path();

    /// checkpoint data
    void checkpoint(std::ostream& cp_file) const;

    /// restore data
    void restore(std::istream& cp_file);

    void addVesicle(solver::vesicle_global_id ves_idx,
                    double speed,
                    const std::map<solver::spec_global_id, uint>& spec_deps,
                    const std::vector<double>& stoch_stepsize,
                    double binding_rate,
                    double min_binding_radius,
                    double max_binding_radius,
                    double unbinding_rate,
                    bool allow_path_intersection);

    void addPoint(uint index, const std::array<double, 3>& position);

    void addEdge(uint source, uint dest, double weight, bool allow_binding);

    std::map<uint, std::pair<std::vector<double>, std::map<uint, double>>> getPathMap() const;

    double getBindingRate(const Vesicle& ves) const;

    double getUnBindingRate(const Vesicle& ves) const;

    bool intersectsSphere(const math::position_abs& center,
                          const double& min_radius,
                          const double& max_radius,
                          bool binding_edges_only) const;

    bool canIntersect(const Vesicle& ves) const;

    solver::path_edge_global_id selectEdge(uint source, const rng::RNGptr& rng) const;

    /// Return the patch id.
    ///
    /// \return ID of the patch.
    inline solver::path_global_id getID() const noexcept {
        return pID;
    }

    inline const std::string& getName() const noexcept {
        return pName;
    }

    // calculate a route through the path based on the vesicle speed, the vesicle
    // update dt and a random choice of branches
    // Keep this quite dynamic in order to allow for adding/removing/growing Paths
    // in the simulation in the future The vesicle radius is there really just for
    // some checks.
    std::pair<std::vector<std::pair<double, math::position_abs>>, math::position_abs>
    calculateRoute(const Vesicle& ves, const rng::RNGptr& rng) const;

    ////////////////////////////////////////////////////////////////////////

  private:
    solver::path_global_id pID;
    std::string pName;

    // Whether vesicles should automatically bind to the path if they intersect point 1
    bool pBind_to_start{true};

    // List of vesicles that belong to this path by index, along with their information
    std::map<solver::vesicle_global_id, PathVesicleInfo> pVesicles;

    std::map<uint, math::position_abs> pPoints;
    util::strongid_vector<solver::path_edge_global_id, PathEdge> pEdges;

    // Map between points and all edges that start at that point
    std::map<uint, std::vector<solver::path_edge_global_id>> pEdgeMap;

    // Spatial index for path edges
    typedef boost::geometry::model::segment<math::position_abs> segment;
    typedef std::pair<segment, solver::path_edge_global_id> rtree_val;

    boost::geometry::index::rtree<rtree_val, boost::geometry::index::rstar<8>> edgeIndex;
};

}  // namespace steps::mpi::tetvesicle
