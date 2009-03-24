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
// $Id$
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>
#include <sstream>

// STEPS headers.
#include <steps/common.h>
#include <steps/error.hpp>
#include <steps/solver/api.hpp>
#include <steps/solver/statedef.hpp>

////////////////////////////////////////////////////////////////////////////////

USING(std, string);
USING_NAMESPACE(steps::solver);

////////////////////////////////////////////////////////////////////////////////

double API::getPatchArea(string const & p) const
{
	// the following may raise an exception if string is unused
	uint pidx = pStatedef->getPatchIdx(p);

	return _getPatchArea(pidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchArea(string const & p, double area)
{
	if (area <= 0.0)
	{
		std::ostringstream os;
		os << "Area cannot be negative or zero.";
		throw steps::ArgErr(os.str());
	}
	// the following may raise an exception if string is unused
	uint pidx = pStatedef->getPatchIdx(p);

	_setPatchArea(pidx, area);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchCount(string const & p, string const & s) const
{
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sidx = pStatedef->getSpecIdx(s);

	return _getPatchCount(pidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchCount(string const & p, string const & s, double n)
{
	if (n < 0.0)
	{
		std::ostringstream os;
		os << "Number of molecules cannot be negative.";
		throw steps::ArgErr(os.str());
	}
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sidx = pStatedef->getSpecIdx(s);

	_setPatchCount(pidx, sidx, n);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchAmount(string const & p, string const & s) const
{
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sidx = pStatedef->getSpecIdx(s);

	return _getPatchAmount(pidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchAmount(string const & p, string const & s, double a)
{
	if (a < 0.0)
	{
		std::ostringstream os;
		os << "Amount of mols cannot be negative.";
	    throw steps::ArgErr(os.str());
	}
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sidx = pStatedef->getSpecIdx(s);

	_setPatchAmount(pidx, sidx, a);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getPatchClamped(string const & p, string const & s) const
{
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sidx = pStatedef->getSpecIdx(s);

	return _getPatchClamped(pidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchClamped(string const & p, string const & s, bool buf)
{
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sidx = pStatedef->getSpecIdx(s);

	_setPatchClamped(pidx, sidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchSReacK(string const & p, string const & sr) const
{
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sridx = pStatedef->getSReacIdx(sr);

	return _getPatchSReacK(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchSReacK(string const & p, string const & sr, double kf)
{
	if (kf < 0.0)
	{
		std::ostringstream os;
		os << "Reaction constant cannot be negative.";
	    throw steps::ArgErr(os.str());
	}
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sridx = pStatedef->getSReacIdx(sr);

	_setPatchSReacK(pidx, sridx, kf);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getPatchSReacActive(string const & p, string const & sr) const
{
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sridx = pStatedef->getSReacIdx(sr);

	return _getPatchSReacActive(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchSReacActive(string const & p, string const & sr, bool a)
{
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sridx = pStatedef->getSReacIdx(sr);

	_setPatchSReacActive(pidx, sridx, a);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchSReacH(string const & p, string const & sr) const
{
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sridx = pStatedef->getSReacIdx(sr);

	return _getPatchSReacH(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchSReacC(string const & p, string const & sr) const
{
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sridx = pStatedef->getSReacIdx(sr);

	return _getPatchSReacC(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchSReacA(string const & p, string const & sr) const
{
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sridx = pStatedef->getSReacIdx(sr);

	return _getPatchSReacA(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getPatchSReacExtent(string const & p, string const & sr) const
{
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sridx = pStatedef->getSReacIdx(sr);

	return _getPatchSReacExtent(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::resetPatchSReacExtent(string const & p, string const & sr)
{
	// the following may raise exceptions if strings are unused
	uint pidx = pStatedef->getPatchIdx(p);
	uint sridx = pStatedef->getSReacIdx(sr);

	 _resetPatchSReacExtent(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::_setPatchArea(uint pidx, double area)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getPatchSReacH(uint pidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getPatchSReacC(uint pidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getPatchSReacA(uint pidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getPatchSReacExtent(uint pidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_resetPatchSReacExtent(uint pidx, uint ridx)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

// END
