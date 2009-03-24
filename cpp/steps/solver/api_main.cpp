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
#include <steps/geom/geom.hpp>
#include <steps/model/model.hpp>
#include <steps/rng/rng.hpp>
#include <steps/solver/api.hpp>
#include <steps/solver/statedef.hpp>

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

void API::saveState(string const & filename)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::step(void)
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

// END
