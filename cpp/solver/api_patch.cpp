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

// STL headers.
#include <string>
#include <sstream>

// STEPS headers.
#include "../common.h"
#include "../error.hpp"
#include "api.hpp"
#include "statedef.hpp"

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
