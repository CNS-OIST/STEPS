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
#include <cassert>
#include <sstream>
#include <string>

// STEPS headers.
#include "../common.h"
#include "../error.hpp"
#include "model.hpp"
#include "spec.hpp"
#include "surfsys.hpp"
#include "volsys.hpp"

////////////////////////////////////////////////////////////////////////////////

USING_NAMESPACE(std);
USING_NAMESPACE(steps::model);

////////////////////////////////////////////////////////////////////////////////

bool steps::model::isValidID(string const & id)
{
    int idlen = id.length();
    if (idlen == 0) return false;
    const char * s = id.c_str();
    char s0 = s[0];
    if (((s0 >= 'a') && (s0 <= 'z')) ||
        ((s0 >= 'A') && (s0 <= 'Z')) ||
        (s0 == '_')) goto l0;
    return false;
l0:
    int i = 0;
l1:
    if (++i == idlen) return true;
    char si = s[i];
    if (((si >= 'a') && (si <= 'z')) ||
        ((si >= 'A') && (si <= 'Z')) ||
        ((si >= '0') && (si <= '9')) ||
        (si == '_')) goto l1;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void steps::model::checkID(string const & id)
{
    bool r = steps::model::isValidID(id);
    if (r == false)
    {
        ostringstream os;
        os << "'" << id << "' is not a valid id";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

Model::Model(void)
: pSpecs()
, pVolsys()
, pSurfsys()
{
}

////////////////////////////////////////////////////////////////////////////////

Model::~Model(void)
{
    while (pSpecs.empty() == false)
    {
        SpecPMapCI spec = pSpecs.begin();
        delete(spec->second);
    }

    while (pVolsys.empty() == false)
    {
        VolsysPMapCI vsys = pVolsys.begin();
        delete(vsys->second);
    }

    while (pSurfsys.empty() == false)
    {
        SurfsysPMapCI ssys = pSurfsys.begin();
        delete(ssys->second);
    }
}

////////////////////////////////////////////////////////////////////////////////

Spec * Model::getSpec(string const & id) const
{
    SpecPMapCI spec = pSpecs.find(id);
    if (spec == pSpecs.end())
    {
        ostringstream os;
		os << "Model does not contain species with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
    assert(spec->second != 0);
    return spec->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delSpec(string const & id)
{
    Spec * spec = getSpec(id);
    delete(spec);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec *> Model::getAllSpecs(void) const
{
	SpecPVec specs = SpecPVec();
	SpecPMapCI spec_end = pSpecs.end();
	for (SpecPMapCI s = pSpecs.begin(); s != spec_end; ++s)
	{
		specs.push_back(s->second);
	}
	return specs;
}

////////////////////////////////////////////////////////////////////////////////

Volsys * Model::getVolsys(string const & id) const
{
    VolsysPMapCI volsys = pVolsys.find(id);
    if (volsys == pVolsys.end())
    {
        ostringstream os;
        os << "Model does not contain volume system with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
	assert(volsys->second != 0);
	return volsys->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delVolsys(string const & id)
{
    Volsys * volsys = getVolsys(id);
    delete(volsys);
}


////////////////////////////////////////////////////////////////////////////////

Surfsys * Model::getSurfsys(string const & id) const
{
	SurfsysPMapCI surfsys = pSurfsys.find(id);
    if (surfsys == pSurfsys.end())
    {
		ostringstream os;
		os << "Model does not contain surface system with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
	assert(surfsys->second != 0);
    return surfsys->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delSurfsys(string const & id)
{
    Surfsys * surfsys = getSurfsys(id);
    delete(surfsys);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkSpecID(string const & id) const
{
    checkID(id);
    if (pSpecs.find(id) != pSpecs.end())
    {
		ostringstream os;
        os << "'" << id << "' is already in use";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSpecIDChange(string const & o, string const & n)
{
    SpecPMapCI s_old = pSpecs.find(o);
    assert(s_old != pSpecs.end());

    if (o == n) return;
    _checkSpecID(n);

    Spec * s = s_old->second;
    assert(s != 0);
    pSpecs.erase(s->getID());						// or s_old->first
    pSpecs.insert(SpecPMap::value_type(n, s));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSpecAdd(Spec * spec)
{
    assert(spec->getModel() == this);
    _checkSpecID(spec->getID());
    pSpecs.insert(SpecPMap::value_type(spec->getID(), spec));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSpecDel(Spec * spec)
{
    VolsysPMapCI vsys_end = pVolsys.end();
    for (VolsysPMapCI vsys = pVolsys.begin(); vsys != vsys_end; ++vsys)
    {
        vsys->second->_handleSpecDelete(spec);
    }

    SurfsysPMapCI ssys_end = pSurfsys.end();
    for (SurfsysPMapCI ssys = pSurfsys.begin(); ssys != ssys_end; ++ssys)
    {
        ssys->second->_handleSpecDelete(spec);
    }

    pSpecs.erase(spec->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkVolsysID(string const & id) const
{
    checkID(id);
    if (pVolsys.find(id) != pVolsys.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVolsysIDChange(string const & o, string const & n)
{
    VolsysPMapCI v_old = pVolsys.find(o);
    assert(v_old != pVolsys.end());

    if (o == n) return;
    _checkVolsysID(n);

    Volsys * v = v_old->second;
    assert(v != 0);
    pVolsys.erase(v->getID());						// or v_old->first
    pVolsys.insert(VolsysPMap::value_type(n, v));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVolsysAdd(Volsys * volsys)
{
    assert(volsys->getModel() == this);
    _checkVolsysID(volsys->getID());
	pVolsys.insert(VolsysPMap::value_type(volsys->getID(), volsys));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVolsysDel(Volsys * volsys)
{
    assert (volsys->getModel() == this);
    pVolsys.erase(volsys->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkSurfsysID(string const & id) const
{
    checkID(id);
    if (pSurfsys.find(id) != pSurfsys.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSurfsysIDChange(string const & o, string const & n)
{
    SurfsysPMapCI s_old = pSurfsys.find(o);
    assert(s_old != pSurfsys.end());

    if (o==n) return;
    _checkSurfsysID(n);

    Surfsys * s = s_old->second;
    assert(s != 0);
    pSurfsys.erase(s->getID());
    pSurfsys.insert(SurfsysPMap::value_type(n, s));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSurfsysAdd(Surfsys * surfsys)
{
    assert(surfsys->getModel() == this);
    _checkSurfsysID(surfsys->getID());
    pSurfsys.insert(SurfsysPMap::value_type(surfsys->getID(), surfsys));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSurfsysDel(Surfsys * surfsys)
{
	assert (surfsys->getModel() == this);
    pSurfsys.erase(surfsys->getID());
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countReacs(void) const
{
	uint nreacs = 0;

	VolsysPMapCI vs_end = pVolsys.end();
	for (VolsysPMapCI vs = pVolsys.begin(); vs != vs_end; ++vs)
	{
		nreacs += vs->second->_countReacs();
	}
	return nreacs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countSReacs(void) const
{
	uint nsreacs = 0;

	SurfsysPMapCI ss_end = pSurfsys.end();
	for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ss_end; ++ss)
	{
		nsreacs += ss->second->_countSReacs();
	}
	return nsreacs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countDiffs(void) const
{
	uint ndiffs = 0;

	VolsysPMapCI vs_end = pVolsys.end();
	for (VolsysPMapCI vs = pVolsys.begin(); vs != vs_end; ++vs)
	{
		ndiffs += vs->second->_countDiffs();
	}
	return ndiffs;
}

////////////////////////////////////////////////////////////////////////////////

Spec * Model::_getSpec(uint gidx) const
{
	assert (gidx < pSpecs.size());
	std::map<std::string, Spec *>::const_iterator sp_it = pSpecs.begin();
	for (uint i=0; i< gidx; ++i) ++sp_it;
	return sp_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Reac * Model::_getReac(uint gidx) const
{
	// first find which volsys this reac (by global index) belongs to
	uint lidx = gidx;
	VolsysPMapCI vs_end = pVolsys.end();
	for (VolsysPMapCI vs = pVolsys.begin(); vs != vs_end; ++vs)
	{
		uint reacs_tot = vs->second->_countReacs();
		if (reacs_tot > lidx) return vs->second->_getReac(lidx);
		lidx -=reacs_tot;
	}

	// we shouldn't have gotten here
	assert(false);
}

////////////////////////////////////////////////////////////////////////////////

SReac * Model::_getSReac(uint gidx) const
{
	// first find which surfsys this sreac (by global index) belongs to
	uint lidx = gidx;
	SurfsysPMapCI ss_end = pSurfsys.end();
	for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ss_end; ++ss)
	{
		uint sreacs_tot = ss->second->_countSReacs();
		if (sreacs_tot > lidx) return ss->second->_getSReac(lidx);
		lidx -=sreacs_tot;
	}

	// we shouldn't have gotten here
	assert(false);
}

////////////////////////////////////////////////////////////////////////////////

Diff * Model::_getDiff(uint gidx) const
{
	// first find which volsys this diff (by global index) belongs to
	uint lidx = gidx;
	VolsysPMapCI vs_end = pVolsys.end();
	for (VolsysPMapCI vs = pVolsys.begin(); vs != vs_end; ++vs)
	{
		uint diffs_tot = vs->second->_countDiffs();
		if (diffs_tot > lidx) return vs->second->_getDiff(lidx);
		lidx -=diffs_tot;
	}
	// we shouldn't have gotten to the end
	assert(false);
}

////////////////////////////////////////////////////////////////////////////////

// END
