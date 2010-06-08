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
#include "../geom/geom.hpp"
#include "../model/model.hpp"
#include "../rng/rng.hpp"
#include "api.hpp"
#include "statedef.hpp"

////////////////////////////////////////////////////////////////////////////////

USING(std, string);
USING_NAMESPACE(steps::solver);

////////////////////////////////////////////////////////////////////////////////

API::API(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r)
: pModel(m)
, pGeom(g)
, pRNG(r)
, pStatedef(0)
{
    if (pModel == 0)
    {
        std::ostringstream os;
        os << "No model provided to solver initializer function";
        throw steps::ArgErr(os.str());
    }
    if (pGeom == 0)
    {
        std::ostringstream os;
        os << "No geometry provided to solver initializer function";
        throw steps::ArgErr(os.str());
    }
    if (pRNG == 0)
    {
        std::ostringstream os;
        os << "No RNG provided to solver initializer function";
        throw steps::ArgErr(os.str());
    }

    if (m->_countSpecs() == 0)
    {
    	std::ostringstream os;
    	os << "Cannot create solver object with this ";
    	os << "steps.model.Model description object. ";
    	os << "Model must contain at least one chemical Species.";
    	throw steps::ArgErr(os.str());
    }

    if (g->_countComps() == 0)
    {
    	std::ostringstream os;
    	os << "Cannot create solver object with this ";
    	os << "steps.geom.Geom geometry description object. ";
    	os << "Geometry must contain at least one Compartment.";
    	throw steps::ArgErr(os.str());

    }

    std::vector<steps::wm::Comp *> comps = g->getAllComps();
    std::vector<steps::wm::Comp *>::const_iterator c_end = comps.end();
    for (std::vector<steps::wm::Comp *>::const_iterator c = comps.begin(); c != c_end; ++c)
    {
    	if ((*c)->getVol() == 0.0)
    	{
        	std::ostringstream os;
        	os << "Cannot create solver object with this ";
        	os << "steps.geom.Geom geometry description object. ";
        	os << "All Compartments must have non-zero volume.";
        	throw steps::ArgErr(os.str());
    	}
    }

    // create state object, which will in turn create compdef, specdef etc
    //objects and initialise
    pStatedef = new Statedef(m, g, r);
    assert (pStatedef != 0);
}

////////////////////////////////////////////////////////////////////////////////

API::~API(void)
{
	delete pStatedef;
}


////////////////////////////////////////////////////////////////////////////////

void API::step(void)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::advance(double adv)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void  API::setDT(double dt)
{
	throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::getDT(void) const
{
	throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double API::getA0(void) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

uint API::getNSteps(void) const
{
	throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::setTime(double time)
{
	throw steps::NotImplErr();
}
////////////////////////////////////////////////////////////////////////////////

void API::setNSteps(uint nsteps)
{
	throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////
// END
