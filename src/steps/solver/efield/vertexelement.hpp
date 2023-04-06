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


#ifndef STEPS_SOLVER_EFIELD_VERTEXELEMENT_HPP
#define STEPS_SOLVER_EFIELD_VERTEXELEMENT_HPP 1

// STL headers.
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

// STEPS headers.
#include "util/common.h"

////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace solver {
namespace efield {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class VertexElement;
class VertexConnection;
class Mesh;

// Auxiliary declarations.
typedef VertexElement *                     VertexElementP;
typedef std::vector<VertexElementP>         VertexElementPVec;
typedef VertexElementPVec::iterator         VertexElementPVecI;
typedef VertexElementPVec::const_iterator   VertexElementPVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Stores information for a mesh vertex (= mesh point, mesh node, ...).
/// This information includes:
/// <UL>
/// <LI> Its index.
/// <LI> Its coordinates (expressed in micrometer).
/// <LI> The volume surrounding this vertex.
/// <LI> The surface around this vertex.
/// <LI> A list of neighbouring vertices, together with their coupling
///      constants.
/// </UL>
///
class VertexElement
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor.
    ///
    /// \param idx
    ///         The integer index of this vertex.
    /// \param vpos
    ///         A 1D array of size 3, giving the vertex's coordinates.
    ///
    VertexElement(uint idx, double * vpos);

    /// Destructor.
    ///
    ~VertexElement();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////

    /// This function sets the index of this vertex to a new value.
    /// It's not clear why this should happen... Maybe delete in the
    /// future?
    ///
    inline void setIDX(uint i) noexcept
    { pIDX = i; }

    /// Adds an amount of surface area to the surface associated with
    /// this vertex.
    ///
    inline void incrementSurfaceArea(double sa) noexcept
    { pSurface += sa; }

    /// Called by VertexConnection::attachToVertices().
    ///
    inline void addConnection(VertexConnection* vc) noexcept
    { pConnections.push_back(vc); }

    /// Called by TetMesh::extractConnections(). This basically sets up
    /// some additional data structures locally in VertexElement based
    /// on connections to other vertices.
    ///
    void fix();

    inline void setVolume(double d) noexcept
    { pVolume = d; }

    ////////////////////////////////////////////////////////////////////////

    inline void applySurfaceCapacitance(double c) noexcept
    { pCapacitance = c * pSurface; }

    inline void updateCapacitance(double c) noexcept {
        pCapacitance += c;
    }

    void applyConductance(double);

    ////////////////////////////////////////////////////////////////////////
    // GENERAL INFORMATION
    ////////////////////////////////////////////////////////////////////////

    inline uint getIDX() const noexcept
    { return pIDX; }

    inline double getX() const noexcept
    { return pXPos; }

    inline double getY() const noexcept
    { return pYPos; }

    inline double getZ() const noexcept
    { return pZPos; }

    inline double getSurfaceArea() const noexcept
    { return pSurface; }

    inline double getCapacitance() const noexcept
    { return pCapacitance; }

    ////////////////////////////////////////////////////////////////////////
    // CONNECTIVITY INFORMATION
    ////////////////////////////////////////////////////////////////////////

    inline VertexElement* getNeighbor(uint i) const noexcept
    { return pNbrs[i]; }

    inline VertexElement** getNeighbours() const noexcept
    { return pNbrs; }

    inline uint nbrIdx(uint i) const noexcept
    { return pNbrs[i]->getIDX(); }

    inline uint getNCon() const noexcept
    { return pNCon; }

    inline double getCC(uint i) const noexcept
    { return pCcs[i]; }

    ////////////////////////////////////////////////////////////////////////
    // OUTPUT
    ////////////////////////////////////////////////////////////////////////

    /// Print a summary of the vertex to an output stream.
    ///
    // friend std::ostream & operator<< (std::ostream & os, VertexElement const&);

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////
    // GENERAL DATA FIELDS
    ////////////////////////////////////////////////////////////////////////

    uint                        pIDX;

    double                      pXPos;
    double                      pYPos;
    double                      pZPos;

    /// During initialization, the surface area is computed for each
    /// triangle specified in the mesh. One third of this value is added
    /// to the surface area associated with each vertex of the triangle.
    /// The total value stored here thus contains contributions from multiple
    /// triangles. An internal vertex (a vertex that is not on the corner
    /// of any triangle will have a surface area of 0).
    ///
    double                      pSurface;
    double                      pVolume;
    double                      pCapacitance;

    ////////////////////////////////////////////////////////////////////////
    // CONNECTIVITY DATA
    ////////////////////////////////////////////////////////////////////////

    std::vector<VertexConnection*>      pConnections;

    /// Set during VertexElement::fix().
    ///
    uint                                pNCon;

    /// Set during VertexElement::fix().
    ///
    VertexElement                    ** pNbrs;

    /// Set to during VertexElement::fix().
    ///
    /// TODO: do we really need a local copy of the coupling constants?
    /// Couldn't we just get them from the VertexConnection objects?
    /// Check this later (depends on whether the CC's are read out
    /// only during setup, or during simulation).
    ///
    double                            * pCcs;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN
std::ostream & operator<< (std::ostream & os, steps::solver::efield::VertexElement const&);

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SOLVER_EFIELD_VERTEXELEMENT_HPP

// END
