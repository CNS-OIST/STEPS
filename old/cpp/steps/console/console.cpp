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
// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STEPS headers.
#include <steps/common.h>
#include <steps/error.hpp>
#include <steps/console/channel.hpp>
#include <steps/console/console.hpp>

// Standard library & STL headers.
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// Boost libraries.
#include <boost/filesystem.hpp>

NAMESPACE_ALIAS(boost::filesystem, bfs);
NAMESPACE_ALIAS(steps::console, scons);
USING_NAMESPACE(std);

////////////////////////////////////////////////////////////////////////////////

static ofstream * sInfoFile;
static ofstream * sDbgFile;

static void close_info_file(void)
{
    if (sInfoFile == 0) return;
    sInfoFile->close();
    delete sInfoFile;
    sInfoFile = 0;
}

static void close_dbg_file(void)
{
    if (sDbgFile == 0) return;
    sDbgFile->close();
    delete sDbgFile;
    sDbgFile = 0;
}

static void make_info_filename(string & filename)
{
    ostringstream os;
    for (uint i = 0; i < 99999999; ++i)
    {
        os.str("");
        os << "steps.info.";
        os.fill('0');
        os.width(8);
        os << i;
        filename = os.str();
        if (bfs::exists(filename) == false) return;
    }
    throw steps::ProgErr("Cannot generate filename for 'info' channel output.");
}

static void make_dbg_filename(string & filename)
{
    ostringstream os;
    for (uint i = 0; i < 99999999; ++i)
    {
        os.str("");
        os << "steps.dbg.";
        os.fill('0');
        os.width(8);
        os << i;
        filename = os.str();
        if (bfs::exists(filename) == false) return;
    }
    throw steps::ProgErr("Cannot generate filename for 'dbg' channel output.");
}

////////////////////////////////////////////////////////////////////////////////

static bool & initialized(void)
{
    static bool i = false;
    return i;
}

void scons::init(void)
{
    if (initialized() == true) return;
    sInfoFile = 0;
    sDbgFile = 0;
    scons::info().setStream(&cout);
    scons::dbg().setStream(0);
    initialized() = true;
}

////////////////////////////////////////////////////////////////////////////////

void scons::finish(void)
{
    close_info_file();
    close_dbg_file();
}

////////////////////////////////////////////////////////////////////////////////

scons::Channel & scons::info(void)
{
    static scons::Channel info_channel;
    return info_channel;
}

////////////////////////////////////////////////////////////////////////////////

void scons::info_to_null(void)
{
    close_info_file();
    scons::info().setStream(0);
}

////////////////////////////////////////////////////////////////////////////////

void scons::info_to_stdout(void)
{
    close_info_file();
    scons::info().setStream(&cout);
}

////////////////////////////////////////////////////////////////////////////////

void scons::info_to_stderr(void)
{
    close_info_file();
    scons::info().setStream(&cerr);
}

////////////////////////////////////////////////////////////////////////////////

void scons::info_to_file(string const & filename, bool app)
{
    string filename2 = filename;
    if (filename.is_empty() == true)
    {
        make_info_filename(filename2);
    }
    ofstream * ofs;
    if (app == true) ofs = new ofstream(filename2, ios_base::app);
    else ofs = new ofstream(filename2, ios_base::trunc);
    assert(ofs != 0);
    if (!(*ofs))
    {
        delete ofs;
        ostringstream os;
        os << "Opening info channel file \"" << filename2 << "\" failed.";
        throw steps::IOErr(os.str());
    }
    close_info_file();
    sInfoFile = ofs;
    scons::info().setStream(sInfoFile);
}

////////////////////////////////////////////////////////////////////////////////

Channel & scons::dbg(void)
{
    static scons::Channel dbg_channel;
    return dbg_channel;
}

////////////////////////////////////////////////////////////////////////////////

void scons::dbg_to_null(void)
{
    close_dbg_file();
    scons::dbg().setStream(0);
}

////////////////////////////////////////////////////////////////////////////////

void scons::dbg_to_stdout(void)
{
    close_dbg_file();
    scons::dbg().setStream(&cout);
}

////////////////////////////////////////////////////////////////////////////////

void scons::dbg_to_stderr(void)
{
    close_dbg_file();
    scons::dbg().setStream(&cerr);
}

////////////////////////////////////////////////////////////////////////////////

void scons::dbg_to_file(string const & filename, bool app)
{
    string filename2 = filename;
    if (filename.is_empty() == true)
    {
        make_dbg_filename(filename2);
    }
    ofstream * ofs;
    if (app == true) ofs = new ofstream(filename2, ios_base::app);
    else ofs = new ofstream(filename2, ios_base::trunc);
    assert(ofs != 0);
    if (!(*ofs))
    {
        delete ofs;
        ostringstream os;
        os << "Opening debug channel file \"" << filename2 << "\" failed.";
        throw steps::IOErr(os.str());
    }
    close_dbg_file();
    sDbgFile = ofs;
    scons::dbg().setStream(sDbgFile);
}

////////////////////////////////////////////////////////////////////////////////

// END
*/
