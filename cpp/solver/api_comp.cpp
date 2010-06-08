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

double API::getCompVol(string const & c) const
{
	// the following may throw an exception if string is unknown
	uint cidx = pStatedef->getCompIdx(c);

	return _getCompVol(cidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompVol(string const & c, double vol)
{
	if (vol <= 0.0)
	{
		std::ostringstream os;
		os << "Volume cannot be negative or zero.";
	    throw steps::ArgErr(os.str());
	}
	// the following may throw an exception if string is unknown
	uint cidx = pStatedef->getCompIdx(c);

	_setCompVol(cidx, vol);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompCount(string const & c, string const & s) const
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint sidx = pStatedef->getSpecIdx(s);

	return _getCompCount(cidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompCount(string const & c, string const & s, double n)
{
	if (n < 0.0)
	{
		std::ostringstream os;
		os << "Number of molecules cannot be negative.";
	    throw steps::ArgErr(os.str());
	}
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint sidx = pStatedef->getSpecIdx(s);

	_setCompCount(cidx, sidx, n);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompAmount(string const & c, string const & s) const
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint sidx = pStatedef->getSpecIdx(s);

	return _getCompAmount(cidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompAmount(string const & c, string const & s, double a)
{
	if (a < 0.0)
	{
		std::ostringstream os;
		os << "Amount of mols cannot be negative.";
	    throw steps::ArgErr(os.str());
	}
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint sidx = pStatedef->getSpecIdx(s);

	_setCompAmount(cidx, sidx, a);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompConc(string const & c, string const & s) const
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint sidx = pStatedef->getSpecIdx(s);

	return _getCompConc(cidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompConc(string const & c, string const & s, double conc)
{
	if (conc < 0.0)
	{
		std::ostringstream os;
		os << "Concentration cannot be negative.";
	    throw steps::ArgErr(os.str());
	}
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint sidx = pStatedef->getSpecIdx(s);

	_setCompConc(cidx, sidx, conc);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getCompClamped(string const & c, string const & s) const
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint sidx = pStatedef->getSpecIdx(s);

	return _getCompClamped(cidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompClamped(string const & c, string const & s, bool b)
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint sidx = pStatedef->getSpecIdx(s);

	_setCompClamped(cidx, sidx, b);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompReacK(string const & c, string const & r) const
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint ridx = pStatedef->getReacIdx(r);

	return _getCompReacK(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompReacK(string const & c, string const & r, double kf)
{
	if (kf < 0.0)
	{
		std::ostringstream os;
		os << "Reaction constant cannot be negative.";
	    throw steps::ArgErr(os.str());
	}
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint ridx = pStatedef->getReacIdx(r);

	_setCompReacK(cidx, ridx, kf);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getCompReacActive(string const & c, string const & r) const
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint ridx = pStatedef->getReacIdx(r);

	return _getCompReacActive(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompReacActive(string const & c, string const & r, bool a)
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint ridx = pStatedef->getReacIdx(r);

	_setCompReacActive(cidx, ridx, a);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompDiffD(string const & c, string const & d) const
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint didx = pStatedef->getDiffIdx(d);

	return _getCompDiffD(cidx, didx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompDiffD(string const & c, string const & d, double dcst)
{
	if (dcst < 0.0)
	{
		std::ostringstream os;
		os << "Diffusion constant cannot be negative.";
	    throw steps::ArgErr(os.str());
	}
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint didx = pStatedef->getDiffIdx(d);

	_setCompDiffD(cidx, didx, dcst);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getCompDiffActive(string const & c, string const & d) const
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint didx = pStatedef->getDiffIdx(d);

	return _getCompDiffActive(cidx, didx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompDiffActive(string const & c, string const & d, bool act)
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint didx = pStatedef->getDiffIdx(d);

	_setCompDiffActive(cidx, didx, act);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompReacH(string const & c, string const & r) const
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint ridx = pStatedef->getReacIdx(r);

	return _getCompReacH(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompReacC(string const & c, string const & r) const
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint ridx = pStatedef->getReacIdx(r);

	return _getCompReacC(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompReacA(string const & c, string const & r) const
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint ridx = pStatedef->getReacIdx(r);

	return _getCompReacA(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

uint API::getCompReacExtent(string const & c, string const & r) const
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint ridx = pStatedef->getReacIdx(r);

	return _getCompReacExtent(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::resetCompReacExtent(string const & c, string const & r)
{
	// the following may throw exceptions if strings are unknown
	uint cidx = pStatedef->getCompIdx(c);
	uint ridx = pStatedef->getReacIdx(r);

	_resetCompReacExtent(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::_setCompVol(uint cidx, double vol)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getCompDiffD(uint cidx, uint didx) const
{
	throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setCompDiffD(uint cidx, uint didx, double dcst)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getCompDiffActive(uint cidx, uint didx) const
{
	throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setCompDiffActive(uint cidx, uint didx, bool act)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getCompReacH(uint cidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getCompReacC(uint cidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getCompReacA(uint cidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::_getCompReacExtent(uint cidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_resetCompReacExtent(uint cidx, uint ridx)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

// END

