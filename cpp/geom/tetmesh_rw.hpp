////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */

#ifndef STEPS_TETMESH_TETMESH_RW_HPP
#define STEPS_TETMESH_TETMESH_RW_HPP 1


// STEPS headers.
#include "../common.h"

// STL headers
#include <cassert>
#include <string>

START_NAMESPACE(steps)
START_NAMESPACE(tetmesh)

////////////////////////////////////////////////////////////////////////////////

// Forward & auxiliary declarations.
class Tetmesh;

////////////////////////////////////////////////////////////////////////////////

//@{
/// loadASCII() and saveASCII() are functions that read and write a
/// tetmesh to a very simple STEPS-centered ASCII (i.e. text based)
/// format. The format is straightforward and consists of five parts,
/// following the five parts of the tetmesh itself.
///
/// <OL>
/// <LI>A list of vertices: first the number of vertices is given. Then
///     for each vertex there is a line with its x, y and z coordinates
///     in scientific notation.
/// <LI>A list of triangles. First the total number of triangles is given.
///     After this, we have 1 line per triangle giving its three corner
///     points as indices into the vertex list.
/// <LI>A list of tetrahedrons. On a first line, the number of tets is
///     specified. After this, for each tet, follows a line with its four
///     corner points, given as indices into the vertex list.
/// <LI>The compartments in the mesh. First the number of compartments is
///     given on a single line. Then each of these compartments is
///     specified using the following structure:
///     <OL>
///     <LI>Its name (ID string) on a single line.
///     <LI>The number of volume systems that have been added to the
///         compartment.
///     <LI>For each volume system, its name (ID string) is given on a
///         separate line.
///     <LI>The number of tetrahedrons associated with the compartment.
///     <LI>The indices of these tetrahedrons (8 per line).
///     </OL>
/// <LI>A similar approach is taken to list all patches in the meshes.
///     first, we specify the number of patches attached to the mesh.
///     Then, for each patch, we describe like this:
///     <OL>
///     <LI>The patch's name on a single line.
///     <LI>If the patch has NO inner compartment, we print a "0"
///         on a single line. Otherwise the line will start with a "1",
///         followed by the name (ID string) of this inner compartment.
///     <LI>The same for the outer compartment.
///     <LI>The number of surface systems that have been added to the
///         patch.
///     <LI>For each of these surface systems, we specify its name on
///         a separate line.
///     <LI>The number of triangles in the patch.
///     <LI>The indices of these triangles (8 per line).
///     </OL>
/// <LI>
/// </OL>
///
/// All indices start from zero in this format, as they do in
/// STEPS (and C++) internally.
///
/// \todo More care could be taken in handling wrongly specified
/// ASCII files. A binary format, or better even a format based a
/// library like NetCDF (or HDF), would be more efficient in terms
/// occupying disk space and could make it easier to read meshes from
/// other software such as Matlab. (But first see whether Matlab's
/// support for these formats is actually as easy as they say it is.)
///
Tetmesh * loadASCII(std::string pathname);
void saveASCII(std::string pathname, Tetmesh * m);
//@}

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetmesh)
END_NAMESPACE(steps)

#endif
// STEPS_TETMESH_TETMESH_RW_HPP

// END
