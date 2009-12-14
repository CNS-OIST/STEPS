////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
//
// This file is part of STEPS.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
//
//
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_TETMESH_TMCOMP_HPP
#define STEPS_TETMESH_TMCOMP_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STEPS headers.
#include <steps/common.h>
#include <steps/geom/comp.hpp>
#include <steps/geom/tetmesh.hpp>
#include <steps/model/volsys.hpp>


// STL headers
#include <vector>
#include <string>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(tetmesh)

////////////////////////////////////////////////////////////////////////////////

// Forward & auxiliary declarations.
class Tetmesh;

////////////////////////////////////////////////////////////////////////////////

/// Provides annotation for a group of tetrahedron in a Tetmesh.
///
///
/// \warning Methods start with an underscore are not exposed to Python.
class TmComp : public steps::wm::Comp
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor.
    ///
    /// \param id ID of the TmComp.
	/// \param container Temesh container for the tetrahedrons.
    /// \param tets A sequence of tetrahedron (by index) as a vector
	///             of unsigned integers which is represented as a
    ///             sequence of positive integer values) in Python.
    /// \param volsys Pointer to the volume system associated.
	///
    TmComp(std::string const & id, Tetmesh * container,
    	 std::vector<uint> const & tets);

    /// Destructor.
    ~TmComp(void);

    ////////////////////////////////////////////////////////////////////////
    // BASE CLASS METHODS
    ////////////////////////////////////////////////////////////////////////

	void setVol(double vol);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

	/// Return a list of all tetrahedron by indices.
    ///
    /// \return List of indices of the tetrahedrons.
    inline std::vector<uint> getAllTetIndices(void) const
    { return pTet_indices; }

    // Return whether tetrahedrons (specified by index) are inside this compartment.
    ///
    /// \param tet List of indices of tetrahedrons.
    /// \return List of results of the tetrahedrons are inside the compartment.
    std::vector<bool> isTetInside(std::vector<uint> tet) const;

    /// Get the minimal coordinate of the rectangular bounding box.
    ///
    /// \return Minimal coordinate of the rectangular bounding box.
    std::vector<double> getBoundMin(void) const;
    /// Get the maximal coordinate of the rectangular bounding box.
    ///
    /// \return Maximal coordinate of the rectangular bounding box.
    std::vector<double> getBoundMax(void) const;

	////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO C++)
    ////////////////////////////////////////////////////////////////////////

    /// Return all tetrahedrons (by index) in the compartment.
    ///
    /// \return List of indices of tetrahedrons.
    inline std::vector<uint> const & _getAllTetIndices(void) const
    { return pTet_indices; }

private:

	////////////////////////////////////////////////////////////////////////

	Tetmesh                           * pTetmesh;
	std::vector<uint>                   pTet_indices;
	uint								pTetsN;

	/// Information about the minimal and maximal boundary values
	///
	double                              pXmin;
    double                              pXmax;
    double                              pYmin;
    double                              pYmax;
    double                              pZmin;
    double                              pZmax;

	////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetmesh)
END_NAMESPACE(steps)

#endif
// STEPS_TETMESH_TMCOMP_HPP

// END
