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

/// \namespace steps::console
///
/// The purpose of steps::console is to provide a mechanism for STEPS
/// to report its current status even in non-interactive modes of use.
/// Currently, this namespace provides two 'channels' on which STEPS
/// code can report:
///
/// <UL>
/// <LI>The 'info' channel reports information at fixed points in the
///     program's life. Examples: a welcome message upon succesful
///     initialization, start of a new simulation, etc.
///
///     Also reported on the info channel are warnings; these are
///     issued when STEPS detects that something might not be entirely
///     in order (usually something of a numerical nature), however in
///     a way that does not warrant immediate interruption of the current
///     activity. (In contrast, real errors are reported by throwing
///     exceptions that do interrupt the current activity -- as is the
///     custom in Python-based software systems).
///
/// <LI>The 'debug' channel reports the current activities of the software
///     in a (much) more detailed manner.
/// </UL>
///
/// Both channels are implemented as output streams.
///

// END
