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

double API::getTetVol(uint tidx) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

	if (tidx >= mesh->countTets())
	{
		std::ostringstream os;
		os << "Tetrahedron index out of range.";
	    throw steps::ArgErr(os.str());
	}
	return _getTetVol(tidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetVol(uint tidx, double vol)
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());
	try{
	if (tidx >= mesh->countTets())
	{
		std::ostringstream os;
		os << "Tetrahedron index out of range.";
	    throw steps::ArgErr(os.str());
	}}
	catch(...)
	{}
	// NOTE: the following method may never be implemented
	_setTetVol(tidx, vol);
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetCount(uint tidx, string const & s) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

void API::setTetCount(uint tidx, string const & s, double n)
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

double API::getTetAmount(uint tidx, string const & s) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

void API::setTetAmount(uint tidx, string const & s, double m)
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

double API::getTetConc(uint tidx, string const & s) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

void API::setTetConc(uint tidx, string const & s, double c)
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

bool API::getTetClamped(uint tidx, string const & s) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

void API::setTetClamped(uint tidx, string const & s, bool buf)
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacK(uint tidx, string const & r) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

void API::setTetReacK(uint tidx, string const & r, double kf)
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

bool API::getTetReacActive(uint tidx, string const & r) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

void API::setTetReacActive(uint tidx, string const & r, bool act)
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

double API::getTetDiffD(uint tidx, string const & d) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

void API::setTetDiffD(uint tidx, string const & d, double dk)
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

bool API::getTetDiffActive(uint tidx, string const & d) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

void API::setTetDiffActive(uint tidx, string const & d, bool act)
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacH(uint tidx, string const & r) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacC(uint tidx, string const & r) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacA(uint tidx, string const & r) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

////////////////////////////////////////////////////////////////////////////////

double API::getTetDiffA(uint tidx, string const & d) const
{
	steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());

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

