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

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <sstream>
#include <string>
#include <map>

// STEPS headers.
#include "../common.h"
#include "../error.hpp"
#include "model.hpp"
#include "spec.hpp"
#include "volsys.hpp"
#include "reac.hpp"
#include "diff.hpp"

////////////////////////////////////////////////////////////////////////////////

USING_NAMESPACE(std);
USING_NAMESPACE(steps::model);

////////////////////////////////////////////////////////////////////////////////

Volsys::Volsys(string const & id, Model * model)
: pID(id)
, pModel(model)
, pReacs()
, pDiffs()
{
    if (pModel == 0)
    {
        ostringstream os;
        os << "No model provided to Volsys initializer function";
        throw steps::ArgErr(os.str());
    }
    pModel->_handleVolsysAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Volsys::~Volsys(void)
{
    if (pModel == 0) return;
	_handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleSelfDelete(void)
{
	std::vector<steps::model::Reac *> allreacs = getAllReacs();
	ReacPVecCI reac_end = allreacs.end();
	for(ReacPVecCI reac = allreacs.begin(); reac != reac_end; ++reac)
	{
		delete(*reac);
	}

	std::vector<steps::model::Diff *> alldiffs = getAllDiffs();
	DiffPVecCI diff_end = alldiffs.end();
	for (DiffPVecCI diff = alldiffs.begin(); diff != diff_end; ++diff)
	{
		delete(*diff);
	}
	pModel->_handleVolsysDel(this);
	pReacs.clear();
	pDiffs.clear();
	pModel = 0;
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::setID(string const & id)
{
    assert(pModel != 0);
    if (id == pID) return;
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pModel->_handleVolsysIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

Reac * Volsys::getReac(string const & id) const
{
    ReacPMapCI reac = pReacs.find(id);
    if (reac == pReacs.end())
    {
        ostringstream os;
        os << "Model does not contain reaction with name '" << id << "'";
		throw steps::ArgErr(os.str());
	}
	assert(reac->second != 0);
    return reac->second;
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::delReac(string const & id)
{
	Reac * reac = getReac(id);
	// Delete reac object since it is owned by c++, not python
	delete(reac);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Reac *> Volsys::getAllReacs(void) const
{
	ReacPVec reacs = ReacPVec();
	ReacPMapCI r_end = pReacs.end();
	for (ReacPMapCI r = pReacs.begin(); r != r_end; ++r)
	{
		reacs.push_back(r->second);
	}
	return reacs;
}

////////////////////////////////////////////////////////////////////////////////

Diff * Volsys::getDiff(string const & id) const
{
    DiffPMapCI diff = pDiffs.find(id);
    if (diff == pDiffs.end())
    {
        ostringstream os;
        os << "Model does not contain diffusion with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
    assert(diff->second != 0);
    return diff->second;
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::delDiff(string const & id)
{
    Diff * diff = getDiff(id);
	// delete diff object since it is owned by c++, not python
	delete(diff);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Diff *> Volsys::getAllDiffs(void) const
{
	DiffPVec diffs = DiffPVec();
	DiffPMapCI d_end = pDiffs.end();
	for (DiffPMapCI d = pDiffs.begin(); d != d_end; ++d)
	{
		diffs.push_back(d->second);
	}
	return diffs;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec *> Volsys::getAllSpecs(void) const
{
	SpecPVec specs = SpecPVec();
	bool first_occ = true;

	ReacPVec reacs = getAllReacs();
	ReacPVecCI reac_end = reacs.end();
	for(ReacPVecCI reac = reacs.begin(); reac != reac_end; ++reac)
	{
		SpecPVec r_specs = (*reac)->getAllSpecs();
		SpecPVecCI r_spec_end = r_specs.end();
		for (SpecPVecCI r_spec = r_specs.begin();
		     r_spec != r_spec_end; ++r_spec)
		{
			first_occ = true;
			SpecPVecCI allspecs_end = specs.end();
			for (SpecPVecCI allspecs = specs.begin();
				allspecs != allspecs_end; ++allspecs)
			{
				if ((*r_spec) == (*allspecs))
				{
					first_occ = false;
					break;
				}
			}
			if (first_occ == true) specs.push_back(*r_spec);
		}
	}

	DiffPVec diffs = getAllDiffs();
	DiffPVecCI diff_end = diffs.end();
	for(DiffPVecCI diff = diffs.begin();diff != diff_end; ++diff)
	{
		SpecPVec d_specs = (*diff)->getAllSpecs();
		SpecPVecCI d_spec_end = d_specs.end();
		for (SpecPVecCI d_spec = d_specs.begin();
			d_spec != d_spec_end; ++d_spec)
		{
			first_occ = true;
			SpecPVecCI allspecs_end = specs.end();
			for (SpecPVecCI allspecs = specs.begin();
				allspecs != allspecs_end; ++allspecs)
			{
				if ((*d_spec) == (*allspecs))
				{
					first_occ = false;
					break;
				}
			}
			if (first_occ == true) specs.push_back(*d_spec);
		}
	}

	return specs;
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_checkReacID(string const & id) const
{
    checkID(id);
    if (pReacs.find(id) != pReacs.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleReacIDChange(string const & o, string const & n)
{
    ReacPMapCI r_old = pReacs.find(o);
    assert(r_old != pReacs.end());

    if(o==n) return;
    _checkReacID(n);

    Reac * r = r_old->second;
    assert(r != 0);
    pReacs.erase(r->getID());
    pReacs.insert(ReacPMap::value_type(n,r));
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleReacAdd(Reac * reac)
{
    assert(reac->getVolsys() == this);
    _checkReacID(reac->getID());
    pReacs.insert(ReacPMap::value_type(reac->getID(), reac));
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleReacDel(Reac * reac)
{
	assert(reac->getVolsys() == this);
    pReacs.erase(reac->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_checkDiffID(string const & id) const
{
    checkID(id);
    if (pDiffs.find(id) != pDiffs.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleDiffIDChange(string const & o, string const & n)
{
    DiffPMapCI d_old = pDiffs.find(o);
    assert(d_old != pDiffs.end());

    if (o==n) return;
    _checkDiffID(n);

    Diff * d = d_old->second;
    assert(d != 0);
    pDiffs.erase(d->getID());
    pDiffs.insert(DiffPMap::value_type(n,d));
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleDiffAdd(Diff * diff)
{
    assert(diff->getVolsys() == this);
    _checkDiffID(diff->getID());
    pDiffs.insert(DiffPMap::value_type(diff->getID(), diff));
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleDiffDel(Diff * diff)
{
	assert (diff->getVolsys() == this);
    pDiffs.erase(diff->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleSpecDelete(Spec * spec)
{
	std::vector<std::string> reacs_del = std::vector<std::string>();
    ReacPMapCI reac_end = pReacs.end();
    for (ReacPMapCI reac = pReacs.begin(); reac != reac_end; ++reac)
    {
		SpecPVec specs = (reac->second->getAllSpecs());
		SpecPVecCI r_spec_end = specs.end();
		for (SpecPVecCI r_spec = specs.begin();
			 r_spec != r_spec_end; ++r_spec)
		{
			if ((*r_spec) == spec)
			{
				reacs_del.push_back(reac->second->getID());
				break;
			}
		}
	}

	std::vector<std::string> diffs_del = std::vector<std::string>();
	DiffPMapCI diff_end = pDiffs.end();
    for (DiffPMapCI diff = pDiffs.begin(); diff != diff_end; ++diff)
    {
		SpecPVec specs = (diff->second->getAllSpecs());
		SpecPVecCI d_spec_end = specs.end();
		for (SpecPVecCI d_spec = specs.begin();
			 d_spec != d_spec_end; ++d_spec)
		{
			if ((*d_spec) == spec)
			{
				diffs_del.push_back(diff->second->getID());
				break;
			}
		}
	}

	std::vector<std::string>::const_iterator r_del_end = reacs_del.end();
	for (std::vector<std::string>::const_iterator r_del = reacs_del.begin();
		 r_del != r_del_end; ++r_del)
	{
		delReac(*r_del);
	}

	std::vector<std::string>::const_iterator d_del_end = diffs_del.end();
	for (std::vector<std::string>::const_iterator d_del = diffs_del.begin();
		 d_del != d_del_end; ++d_del)
	{
		delDiff(*d_del);
	}
}

////////////////////////////////////////////////////////////////////////////////

Reac * Volsys::_getReac(uint lidx) const
{
	assert (lidx < pReacs.size());
	std::map<std::string, Reac *>::const_iterator rc_it = pReacs.begin();
	for (uint i=0; i< lidx; ++i) ++rc_it;
	return rc_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Diff * Volsys::_getDiff(uint lidx) const
{
	assert (lidx < pDiffs.size());
	std::map<std::string, Diff *>::const_iterator df_it = pDiffs.begin();
	for (uint i=0; i< lidx; ++i) ++df_it;
	return df_it->second;
}

////////////////////////////////////////////////////////////////////////////////

// END
