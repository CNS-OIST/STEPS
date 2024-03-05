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

// STL headers.
#include <fstream>
#include <iostream>
#include <memory>

// STEPS headers.
#include "efieldsolver.hpp"
#include "solver/fwd.hpp"
#include "util/vocabulary.hpp"

namespace steps::solver::efield {

// Forward declarations
class TetMesh;

// TODO: Using abstract base class for EField solver functionality
// but consider moving to template impl instead, at least until
// EField internal API is revamped.

////////////////////////////////////////////////////////////////////////////////

// NOTE: THis class provides access to the EFIield calculation from the solver,
// as such all values provided to and returned from EField object are in base
// s.i. units. The EField object converts these values to the required units
// for the matrix calculation.
//
// Two important parameters - the bulk resistivity and membrane capacitance -
// may be set with solver methods, or take default values.
// Default value for resistivity  = 1.0x10^-2 farad / m^2
// Default value for bulk resistivity is 1 ohm.m
class EField {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor for an EField solver object.
    ///
    /// \param impl Solver object for performing EField calculations.
    explicit EField(std::unique_ptr<EFieldSolver> impl);

    /// Initialise EField solver with mesh data.
    ///
    /// \param verts A 1D array of number of vertices in the mesh * 3 doubles in row-major order.
    /// \param tris A 1D array of number of triangle in the mesh * 3 uints in row-major order. These
    ///         indices point into the vertex array.
    /// \param tets  A 1D array of number of tetrahedral elements in the mesh * 4 uints in row-major
    /// order. These
    ///         indices point into the vertex array.
    ///
    void initMesh(const std::vector<double>& verts,
                  const std::vector<vertex_id_t>& tris,
                  const std::vector<vertex_id_t>& tets,
                  uint opt_method = 1,
                  std::string const& opt_file_name = "",
                  double search_percent = 100.0);

    /// Destructor
    ~EField();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    // Save optimal vertex configuration
    void saveOptimal(std::string const& opt_file_name);

    ////////////////////////////////////////////////////////////////////////

    /// Set the surface resistivity of the membrane.
    /// \param midx Membrane identifier
    /// \param res Surface resistivity (ohm.m^2)
    /// \param erev Reversal potential (volts)
    void setSurfaceResistivity(solver::membrane_global_id midx, double res, double erev);

    /// Get the surface resistivity of the membrane.
    /// \param midx membrane index
    /// \return A pair containing Surface resistivity (ohm.m^2) and Reversal potential (volts)
    std::pair<double, double> getSurfaceResistivity() const;

    /// Set the initial electric potential of all points in the conduction volume
    /// including surface nodes.
    /// \param v Potential (volts)
    void setMembPotential(solver::membrane_global_id midx, double v);

    /// Set the specific membrane capacitance of the membrane
    /// \param cm Specific membrane capacitance (farad / m^2)
    void setMembCapac(solver::membrane_global_id midx, double cm);

    /// Set the bulk electrical resistivity of the section of the mesh
    /// representing the volume conductor \param ro Electrical resistivity (ohm.m)
    void setMembVolRes(solver::membrane_global_id midx, double ro);

    ////////////////////////////////////////////////////////////////////////

    /// Return the electric potential of a vertex (volts)
    /// \param vidx Index of the vertex
    /// \return Electric potential of the vertex (volts)
    double getVertV(vertex_id_t vidx);

    /// Set the electric potential of a vertex (volts)
    /// \param vidx Index of the vertex
    /// \param v Electric potential (volts)
    void setVertV(vertex_id_t vidx, double v);

    /// Return whether the electric potential of a vertex is clamped.
    /// \param vidx Index of the vertex
    /// \Return True if electric potential is clamped, False otherwise
    bool getVertVClamped(vertex_id_t vidx);

    /// Set whether the electric potential of a vertex is clamped or not.
    /// \param vidx Index of the vertex
    /// \param cl Clamped (true) or unclamped (false)
    void setVertVClamped(vertex_id_t vidx, bool cl);

    /// Get the current clamp for a vertex element.
    /// \param vidx Index of the vertex
    double getVertIClamp(vertex_id_t vidx);

    /// Set the current clamp for a vertex element.
    /// \param vidx Index of the vertex
    /// \param cur Current clamp for the vertex
    void setVertIClamp(vertex_id_t vidx, double cur);

    ////////////////////////////////////////////////////////////////////////

    /// Return the electric potential of a triangle surface element (volts)
    /// \param tidx Index of the triangle surface element
    /// \return Electric potential of the triangle surface element (volts)
    double getTriV(triangle_local_id tidx);

    /// Set the electric potential of a triangle surface element (volts)
    /// \param tidx Index of the triangle surface element
    /// \param v Electric potential (volts)
    void setTriV(triangle_local_id tidx, double v);

    /// Return whether the electric potential of a triangle surface element is
    /// clamped or not. \param tidx Index of the triangle surface element \return
    /// True if electric potential clamped, False otherwise
    bool getTriVClamped(triangle_local_id tidx);

    /// Set whether the electric potential of a triangle surface element is
    /// clamped or not. \param tidx Index of the triangle \param cl Clamped (true)
    /// or unclamped (false)
    void setTriVClamped(triangle_local_id tidx, bool cl);

    /// Return the total current across a triangle surface element.
    /// \param tidx Index of the triangle surface element
    /// \return Current across the triangle (amps)
    double getTriI(triangle_local_id tidx);

    /// Set the current across a triangle surface element.
    /// \param tidx Index of the triangle surface element
    /// \param cur Current across the triangle (amps)
    void setTriI(triangle_local_id tidx, double cur);

    /// Get the current clamp for a triangle surface element.
    /// \param tidx Index of the triangle surface element
    double getTriIClamp(triangle_local_id tidx);

    /// Set the current clamp for a triangle surface element.
    /// \param tidx Index of the triangle surface element
    /// \param cur Current clamp for the triangle surface element
    void setTriIClamp(triangle_local_id tidx, double cur);

    /// Auxiliary function for setting current in all triangles at once.
    /// \param cur A 1D array, size = number of surface triangles,
    ///     of current across triangles
    // void    setTriI(double * cur);

    /// Set the specific capacitance of a triangle surface element.
    /// \param tidx Index of the triangle surface element
    /// \param cm Specific membrane capacitance (farad / m^2)
    void setTriCapac(triangle_local_id tidx, double cm);

    ////////////////////////////////////////////////////////////////////////

    /// Return electric potential for internal tetrahedral elements (volts).
    /// \param tidx Index of the tetrahedron
    /// \return Electric potential of tetrahedron (volts)
    double getTetV(tetrahedron_local_id tidx);

    /// Set the electric potential of a tetrahedral element (volts)
    /// \param tidx Index of the tetrahedral element
    /// \param v Electric potential (volts)
    void setTetV(tetrahedron_local_id tidx, double v);

    /// Return whether the electric potential of a tetrahedral element is clamped
    /// or not. \param tidx Index of the tetrahedral element \return True if
    /// electric potential clamped, False otherwise
    bool getTetVClamped(tetrahedron_local_id tidx);

    /// Set whether the electric potential of a tetrahedral element is clamped or
    /// not. \param tidx Index of the tetrahedral element \param cl Clamped (true)
    /// or unclamped (false)
    void setTetVClamped(tetrahedron_local_id tidx, bool cl);

    ////////////////////////////////////////////////////////////////////////
    // SIMULATION CONTROLS
    ////////////////////////////////////////////////////////////////////////

    /// Advance the EField simulation.
    /// \param sec The time to advance the EField simulation (seconds)
    void advance(double sec);

    ////////////////////////////////////////////////////////////////////////

  private:
    TetMesh* pMesh{nullptr};
    std::unique_ptr<EFieldSolver> pVProp;
    std::vector<vertex_id_t> pCPerm;

    index_t pNVerts{0};
    index_t pNTris{0};
    index_t pNTets{0};

    std::vector<vertex_id_t> pTritoVert;
};

////////////////////////////////////////////////////////////////////////////////

/** Convenience function for constructing EFieldSolver object and
 * correspondingly initialised EField object.
 */

template <typename S, typename... Args>
std::unique_ptr<EField> make_EField(Args&&... args) {
    std::unique_ptr<EFieldSolver> ef_solver(new S(std::forward<Args>(args)...));
    return std::unique_ptr<EField>(new EField(std::move(ef_solver)));
}

}  // namespace steps::solver::efield
