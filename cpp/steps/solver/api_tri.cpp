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
#include <steps/geom/tetmesh.hpp>

////////////////////////////////////////////////////////////////////////////////

USING(std, string);
USING_NAMESPACE(steps::solver);

////////////////////////////////////////////////////////////////////////////////

double API::getTriArea(uint tidx) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

	if (tidx >= mesh->countTris())
	{
		std::ostringstream os;
		os << "Triangle index out of range.";
	    throw steps::ArgErr(os.str());
	}
	return _getTriArea(tidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriArea(uint tidx, double area)
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

	if (tidx >= mesh->countTris())
	{
		std::ostringstream os;
		os << "Triangle index out of range.";
	    throw steps::ArgErr(os.str());
	}
	// NOTE: the following method may never be implemented
	_setTriArea(tidx, area);
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriCount(uint tidx, string const & s) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

	if (tidx >= mesh->countTris())
	{
		std::ostringstream os;
		os << "Triangle index out of range.";
	    throw steps::ArgErr(os.str());
	}
	// the following may raise exceptions if strings are unused
	uint sidx = pStatedef->getSpecIdx(s);

	return _getTriCount(tidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriCount(uint tidx, string const & s, double n)
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

	if (tidx >= mesh->countTris())
	{
		std::ostringstream os;
		os << "Triangle index out of range.";
	    throw steps::ArgErr(os.str());
	}
	if (n < 0.0)
	{
		std::ostringstream os;
		os << "Number of molecules cannot be negative.";
	    throw steps::ArgErr(os.str());
	}
	// the following may raise exception if string is uknown
	uint sidx = pStatedef->getSpecIdx(s);

	_setTriCount(tidx, sidx, n);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriClamped(uint tidx, string const & s) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

	if (tidx >= mesh->countTris())
	{
		std::ostringstream os;
		os << "Triangle index out of range.";
	    throw steps::ArgErr(os.str());
	}
	// the following may raise exception if string is unknown
	uint sidx = pStatedef->getSpecIdx(s);

	return _getTriClamped(tidx, sidx);

}

////////////////////////////////////////////////////////////////////////////////

void API::setTriClamped(uint tidx, string const & s, bool buf)
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

	if (tidx >= mesh->countTris())
	{
		std::ostringstream os;
		os << "Triangle index out of range.";
	    throw steps::ArgErr(os.str());
	}
	// the following may raise exception if string is unknown
	uint sidx = pStatedef->getSpecIdx(s);

	_setTriClamped(tidx, sidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacK(uint tidx, string const & r) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

	if (tidx >= mesh->countTris())
	{
		std::ostringstream os;
		os << "Triangle index out of range.";
	    throw steps::ArgErr(os.str());
	}
	// the following may raise exception if string is unknown
	uint sridx = pStatedef->getSReacIdx(r);

	return _getTriSReacK(tidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriSReacK(uint tidx, string const & r, double kf)
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

	if (tidx >= mesh->countTris())
	{
		std::ostringstream os;
		os << "Triangle index out of range.";
	    throw steps::ArgErr(os.str());
	}
	if (kf < 0.0)
	{
		std::ostringstream os;
		os << "Reaction constant cannot be negative.";
	    throw steps::ArgErr(os.str());
	}
	// the following may raise exception if string is unknown
	uint sridx = pStatedef->getSReacIdx(r);

	_setTriSReacK(tidx, sridx, kf);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriSReacActive(uint tidx, string const & r) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

	if (tidx >= mesh->countTris())
	{
		std::ostringstream os;
		os << "Triangle index out of range.";
	    throw steps::ArgErr(os.str());
	}
	// the following may raise exception if string is unknown
	uint sridx = pStatedef->getSReacIdx(r);

	return _getTriSReacActive(tidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriSReacActive(uint tidx, string const & r, bool act)
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

	if (tidx >= mesh->countTris())
	{
		std::ostringstream os;
		os << "Triangle index out of range.";
	    throw steps::ArgErr(os.str());
	}

	// the following may raise exception if string is unknown
	uint sridx = pStatedef->getSReacIdx(r);

	_setTriSReacActive(tidx, sridx, act);
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacH(uint tidx, string const & r) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

	if (tidx >= mesh->countTris())
	{
		std::ostringstream os;
		os << "Triangle index out of range.";
	    throw steps::ArgErr(os.str());
	}
	// the following may raise exception if string is unknown
	uint sridx = pStatedef->getSReacIdx(r);

	return _getTriSReacH(tidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacC(uint tidx, string const & r) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

	if (tidx >= mesh->countTris())
	{
		std::ostringstream os;
		os << "Triangle index out of range.";
	    throw steps::ArgErr(os.str());
	}
	// the following may raise exception if string is unknown
	uint sridx = pStatedef->getSReacIdx(r);

	return _getTriSReacC(tidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacA(uint tidx, string const & r) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

	if (tidx >= mesh->countTris())
	{
		std::ostringstream os;
		os << "Triangle index out of range.";
	    throw steps::ArgErr(os.str());
	}
	// the following may raise exception if string is unknown
	uint sridx = pStatedef->getSReacIdx(r);

	return _getTriSReacA(tidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriArea(uint tidx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriArea(uint tidx, double area)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriCount(uint tidx, uint sidx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriCount(uint tidx, uint sidx, double n)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriClamped(uint tidx, uint sidx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriClamped(uint tidx, uint sidx, bool buf)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacK(uint tidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriSReacK(uint tidx, uint ridx, double kf)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriSReacActive(uint tidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriSReacActive(uint tidx, uint ridx, bool act)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacH(uint tidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacC(uint tidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacA(uint tidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

// END

