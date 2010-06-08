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

/// \namespace steps::solver
///
/// This namespace contains all the code that is shared by each solver
/// implemented in STEPS. This currently includes the code for the
/// solver API and the code that stores information about dependencies
/// and indices for the kinetic model.
///
/// Classes CompDef, ReacDef, SpecDef and StateDef together define the
/// layout of a state. Unlike the Python classes in package steps.model,
/// this classes are designed for . This means they have little or no
/// flexibility to allow for changes to the state during simulation (i.e.
/// adding/deleting new species, reactions, compartments etc). They
/// work largely on the basis of integer indices.
///
/// When studying these classes, the key method to examine is
/// StateDef::
///
/// The order in which setupFinal() methods are called is as follows:
/// <OL>
/// <LI> SpecDef::setupFinal() </LI>
/// <LI> ReacDef::setupFinal() </LI>
/// <LI> DiffDef::setupFinal() </LI>
/// <LI> CompDef::setupLocalIndices() </LI>
/// <LI> PatchDef::setupLocalIndices() </LI>
/// <LI> CompDef::setupDependencies() </LI>
/// <LI> PatchDef::setupDependencies() </LI>
/// </OL>
///
/// Generally speaking, this group of classes (or rather, this whole
/// directory) can use some review to deal with error checking/handling
/// a bit better. Do this when it has matured, and then do it from the
/// perspective of the using methods. Kinda vague, I know...
///

// END
