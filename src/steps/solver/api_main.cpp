/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */


// STL headers.
#include <sstream>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/geom/geom.hpp"
#include "steps/model/model.hpp"
#include "steps/rng/rng.hpp"
#include "steps/solver/api.hpp"
#include "steps/solver/statedef.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

USING(std, string);
using namespace steps::solver;

////////////////////////////////////////////////////////////////////////////////

API::API(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r)
: pModel(m)
, pGeom(g)
, pRNG(r)
, pStatedef(nullptr)
{
    if (pModel == nullptr)
    {
        std::ostringstream os;
        os << "No model provided to solver initializer function";
        ArgErrLog(os.str());
    }
    if (pGeom == nullptr)
    {
        std::ostringstream os;
        os << "No geometry provided to solver initializer function";
        ArgErrLog(os.str());
    }

    // r is allowed to be null pointer for deterministic solvers

    if (m->_countSpecs() == 0)
    {
        std::ostringstream os;
        os << "Cannot create solver object with this ";
        os << "steps.model.Model description object. ";
        os << "Model must contain at least one chemical Species.";
        ArgErrLog(os.str());
    }

    if (g->_countComps() == 0)
    {
        std::ostringstream os;
        os << "Cannot create solver object with this ";
        os << "steps.geom.Geom geometry description object. ";
        os << "Geometry must contain at least one Compartment.";
        ArgErrLog(os.str());

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
            ArgErrLog(os.str());
        }
    }

    // create state object, which will in turn create compdef, specdef etc
    //objects and initialise
    pStatedef = new Statedef(m, g, r);
    AssertLog(pStatedef != 0);
}

////////////////////////////////////////////////////////////////////////////////

API::~API()
{
    delete pStatedef;
}


////////////////////////////////////////////////////////////////////////////////

void API::step()
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::advance(double adv)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void  API::setRk4DT(double dt)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getRk4DT() const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void  API::setDT(double dt)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getDT() const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setEfieldDT(double efdt)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getEfieldDT() const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setTemp(double temp)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getTemp() const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::getA0() const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

uint API::getNSteps() const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::setTime(double time)
{
    NotImplErrLog("");
}
////////////////////////////////////////////////////////////////////////////////

void API::setNSteps(uint nsteps)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////
// END
