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
#include "../geom/tetmesh.hpp"

////////////////////////////////////////////////////////////////////////////////

USING(std, string);
USING_NAMESPACE(steps::solver);

////////////////////////////////////////////////////////////////////////////////

double API::getTriArea(uint tidx) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTris())
		{
			std::ostringstream os;
			os << "Triangle index out of range.";
			throw steps::ArgErr(os.str());
		}
		return _getTriArea(tidx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriArea(uint tidx, double area)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTris())
		{
			std::ostringstream os;
			os << "Triangle index out of range.";
			throw steps::ArgErr(os.str());
		}
		// NOTE: the following method may never be implemented
		_setTriArea(tidx, area);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriCount(uint tidx, string const & s) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
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

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriCount(uint tidx, string const & s, double n)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
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

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriAmount(uint tidx, string const & s) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTris())
		{
			std::ostringstream os;
			os << "Triangle index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint sidx = pStatedef->getSpecIdx(s);

		return _getTriAmount(tidx, sidx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriAmount(uint tidx, string const & s, double m)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTris())
		{
			std::ostringstream os;
			os << "Triangle index out of range.";
			throw steps::ArgErr(os.str());
		}
		if (m < 0.0)
		{
			std::ostringstream os;
			os << "Amount of mols cannot be negative.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint sidx = pStatedef->getSpecIdx(s);

		_setTriAmount(tidx, sidx, m);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriClamped(uint tidx, string const & s) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
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

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriClamped(uint tidx, string const & s, bool buf)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
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

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacK(uint tidx, string const & r) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
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

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriSReacK(uint tidx, string const & r, double kf)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
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

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriSReacActive(uint tidx, string const & r) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
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

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriSReacActive(uint tidx, string const & r, bool act)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
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

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacH(uint tidx, string const & r) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
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

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacC(uint tidx, string const & r) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
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

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacA(uint tidx, string const & r) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
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

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
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

double API::_getTriAmount(uint tidx, uint sidx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriAmount(uint tidx, uint sidx, double m)
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

