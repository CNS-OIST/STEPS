//////////////////////////////////////////////////////
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
// $Id$
//////////////////////////////////////////////////////

#ifndef STEPS_CONSOLE_CONSOLE_HPP
#define STEPS_CONSOLE_CONSOLE_HPP 1
/*
// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STEPS headers.
#include <steps/common.h>
#include <steps/console/channel.hpp>

// Standard library & STL headers.
#include <cassert>
#include <string>

START_NAMESPACE(steps)
START_NAMESPACE(console)

//////////////////////////////////////////////////////

// Initializes the steps::console namespace. This method is called by
// steps::init().
//
STEPS_EXTERN
void init(void);

// Clean up the streams before STEPS exits/unloads. It's called
// automatically by steps::finish().
//
STEPS_EXTERN
void finish(void);

//////////////////////////////////////////////////////

// Returns a reference to the 'info' channel. Implemented using the
// singleton pattern.
//
// This should be made available to STEPS internals.
//
STEPS_EXTERN
Channel & info(void);

// Tie the info stream to nothing. This will discard any info messages.
//
STEPS_EXTERN
void info_to_null(void);

// Tie the info stream to standard output (this is the default after
// calling steps::console::init()). If the stream was tied to a file
// before, this file will be closed first.
//
STEPS_EXTERN
void info_to_stdout(void);

// Tie the info stream to standard error. If the stream was tied to
// a file before, this file will be closed first.
//
STEPS_EXTERN
void info_to_stderr(void);

// Tie the info stream to some file, specified by its filename. If
// the stream was tied to another file before, this file will be
// closed first.
//
// \param filename
// The name of the file to which the channel should be streamed. If
// this parameter is empty, a filename is chosen automatically
// (steps.info.xxxxxxxx) with xxxxxxxx starting from 00000000.
//
// \param app
// If true, the channel output is appended to an existing file.
// If false and the file already exists, it is truncated. The default
// value is false.
//
// \throw steps::ProgErr
// When it's impossible to automatically generate a filename because
// all filenames have been taken. (Shouldn't be a very common
// exception :-)
//
STEPS_EXTERN
void info_to_file(std::string const & filename = "", bool app = false);

//////////////////////////////////////////////////////

// Returns a reference to the debug ('dbg' in short) channel. Implemented
// using the singleton pattern.
//
// This should be made available only to STEPS internals.
//
STEPS_EXTERN
Channel & dbg(void);

// Tie the debug stream to nothing. This will discard any debug messages.
// This is the default after calling steps::console::init().
//
STEPS_EXTERN
void dbg_to_null(void);

// Tie the debug stream to standard output. If the stream was tied to a
// file before, this file will be closed first.
//
STEPS_EXTERN
void dbg_to_stdout(void);

// Tie the debug stream to standard error. If the stream was tied to a file
// before, this file will be closed first.
//
STEPS_EXTERN
void dbg_to_stderr(void);

// Tie the debug stream to some file, specified by its filename. If
// the stream was tied to another file before, this file will be
// closed first.
//
// \param filename
// The name of the file to which the channel should be streamed. If
// this parameter is empty, a filename is chosen automatically
// (steps.dbg.xxxxxxxx) with xxxxxxxx starting from 00000000.
//
// \param app
// If true, the channel output is appended to an existing file.
// If false and the file already exists, it is truncated. The default
// value is false.
//
// \throw steps::ProgErr
// When it's impossible to automatically generate a filename because
// all filenames have been taken. (Shouldn't be a very common
// exception :-)
//
STEPS_EXTERN
void dbg_to_file(std::string const & filename = "", bool app = false);

//////////////////////////////////////////////////////

END_NAMESPACE(console)
END_NAMESPACE(steps)

#endif
// STEPS_CONSOLE_CONSOLE_HPP

// END
*/
