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

double API::getTetVol(uint tidx) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}

		return _getTetVol(tidx);
	}
	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetVol(uint tidx, double vol)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// NOTE: the following method may never be implemented
		_setTetVol(tidx, vol);
	}
	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetCount(uint tidx, string const & s) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint sidx = pStatedef->getSpecIdx(s);

		return _getTetCount(tidx, sidx);
	}
	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetCount(uint tidx, string const & s, double n)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		if (n < 0.0)
		{
			std::ostringstream os;
			os << "Number of molecules cannot be negative.";
			throw steps::ArgErr(os.str());
		}

	// the following may throw exception if string is unknown
	uint sidx = pStatedef->getSpecIdx(s);

	_setTetCount(tidx, sidx, n);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetAmount(uint tidx, string const & s) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint sidx = pStatedef->getSpecIdx(s);

		return _getTetAmount(tidx, sidx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetAmount(uint tidx, string const & s, double m)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
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

		_setTetAmount(tidx, sidx, m);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetConc(uint tidx, string const & s) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint sidx = pStatedef->getSpecIdx(s);

		return _getTetConc(tidx, sidx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetConc(uint tidx, string const & s, double c)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		if (c < 0.0)
		{
			std::ostringstream os;
			os << "Concentration cannot be negative.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint sidx = pStatedef->getSpecIdx(s);

		_setTetConc(tidx, sidx, c);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTetClamped(uint tidx, string const & s) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint sidx = pStatedef->getSpecIdx(s);

		return _getTetClamped(tidx, sidx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetClamped(uint tidx, string const & s, bool buf)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint sidx = pStatedef->getSpecIdx(s);

		_setTetClamped(tidx, sidx, buf);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacK(uint tidx, string const & r) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint ridx = pStatedef->getReacIdx(r);

		return _getTetReacK(tidx, ridx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetReacK(uint tidx, string const & r, double kf)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		if (kf < 0.0)
		{
			std::ostringstream os;
			os << "Reaction constant cannot be negative.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint ridx = pStatedef->getReacIdx(r);

		_setTetReacK(tidx, ridx, kf);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTetReacActive(uint tidx, string const & r) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint ridx = pStatedef->getReacIdx(r);

		return _getTetReacActive(tidx, ridx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetReacActive(uint tidx, string const & r, bool act)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint ridx = pStatedef->getReacIdx(r);

		_setTetReacActive(tidx, ridx, act);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetDiffD(uint tidx, string const & d) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint didx = pStatedef->getDiffIdx(d);

		return _getTetDiffD(tidx, didx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetDiffD(uint tidx, string const & d, double dk)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		if (dk < 0.0)
		{
			std::ostringstream os;
			os << "Diffusion constant cannot be negative.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint didx = pStatedef->getDiffIdx(d);

		_setTetDiffD(tidx, didx, dk);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTetDiffActive(uint tidx, string const & d) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint didx = pStatedef->getDiffIdx(d);

		return _getTetDiffActive(tidx, didx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetDiffActive(uint tidx, string const & d, bool act)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{

		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint didx = pStatedef->getDiffIdx(d);

		_setTetDiffActive(tidx, didx, act);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacH(uint tidx, string const & r) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{

		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint ridx = pStatedef->getReacIdx(r);

		return _getTetReacH(tidx, ridx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacC(uint tidx, string const & r) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint ridx = pStatedef->getReacIdx(r);

		return _getTetReacC(tidx, ridx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacA(uint tidx, string const & r) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{

		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint ridx = pStatedef->getReacIdx(r);

		return _getTetReacA(tidx, ridx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetDiffA(uint tidx, string const & d) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (tidx >= mesh->countTets())
		{
			std::ostringstream os;
			os << "Tetrahedron index out of range.";
			throw steps::ArgErr(os.str());
		}
		// the following may throw exception if string is unknown
		uint didx = pStatedef->getDiffIdx(d);

		return _getTetDiffA(tidx, didx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetVol(uint tidx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetVol(uint tidx, double vol)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetCount(uint tidx, uint sidx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetCount(uint tidx, uint sidx, double n)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetAmount(uint tidx, uint sidx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetAmount(uint tidx, uint sidx, double m)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetConc(uint tidx, uint sidx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetConc(uint tidx, uint sidx, double c)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetClamped(uint tidx, uint sidx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetClamped(uint tidx, uint sidx, bool buf)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacK(uint tidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetReacK(uint tidx, uint ridx, double kf)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetReacActive(uint tidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetReacActive(uint tidx, uint ridx, bool act)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetDiffD(uint tidx, uint didx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetDiffD(uint tidx, uint didx, double dk)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetDiffActive(uint tidx, uint didx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetDiffActive(uint tidx, uint didx, bool act)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacH(uint tidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacC(uint tidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacA(uint tidx, uint ridx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetDiffA(uint tidx, uint didx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

// END

