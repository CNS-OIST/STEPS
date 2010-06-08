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

#ifndef STEPS_ERROR_HPP
#define STEPS_ERROR_HPP 1


// Standard library & STL headers.
#include <cassert>
#include <iostream>
#include <string>

// STEPS headers.
#include "common.h"

START_NAMESPACE(steps)

////////////////////////////////////////////////////////////////////////////////

/// Base STEPS exception class. All 'real' exceptions are derived from this.
///
struct Err
{
    Err(std::string const & msg = "")
    : pMessage(msg) { }
    const char * getMsg(void);
private:
    std::string                 pMessage;
};

////////////////////////////////////////////////////////////////////////////////

/// This exception gets thrown when some STEPS interface function is not
/// implemented, usually in a solver module. This can be because the
/// function does not make sense for that particular solver, or because
/// the programmer took a shortcut...
///
struct NotImplErr
: public Err
{
    NotImplErr(std::string const & msg = "")
    : Err(msg) { }
};

////////////////////////////////////////////////////////////////////////////////

/// This exception usually gets thrown when some caller-provided argument
/// doesn't make sense immediately. This 'caller' should be restricted
/// to the actual user as much as possible.
///
struct ArgErr
: public Err
{
    ArgErr(std::string const & msg = "")
    : Err(msg) { }
};

////////////////////////////////////////////////////////////////////////////////

/// This exception gets thrown whenever the program knows it has failed,
/// in other words not because of some mistake by the user, but e.g.
/// because of some hard coded limit (running out of index numbers,...)
///
/// In principle, it should be only be thrown in absurdly rare cases.
///
struct ProgErr
: public Err
{
    ProgErr(std::string const & msg = "")
    : Err(msg) { }
};

////////////////////////////////////////////////////////////////////////////////

/// Generic base class for any system-induced error. These exceptions signal
/// truly unexpected situations, such as out-of-memory or I/O problems.
///
struct SysErr
: public Err
{
    SysErr(std::string const & msg = "")
    : Err(msg) { }
};

////////////////////////////////////////////////////////////////////////////////

/// Gets thrown whenever there is a system failure involving I/O, such
/// as dealing with files.
///
struct IOErr
: public SysErr
{
    IOErr(std::string const & msg = "")
    : SysErr(msg) { }
};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(steps)

#endif
// STEPS_ERROR_HPP

// END
