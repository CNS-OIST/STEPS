/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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
#include "api.hpp"
#include "statedef.hpp"
// util
#include "util/error.hpp"
// logging
#include <easylogging++.h>
////////////////////////////////////////////////////////////////////////////////

USING(std, string);
using namespace steps::solver;

////////////////////////////////////////////////////////////////////////////////

API::API(steps::model::Model *m, steps::wm::Geom *g, const rng::RNGptr &r)
    : pModel(m), pGeom(g), pRNG(r), pStatedef(nullptr) {
  ArgErrLogIf(pModel == nullptr,
              "No model provided to solver initializer function");
  ArgErrLogIf(pGeom == nullptr,
              "No geometry provided to solver initializer function");
  ArgErrLogIf(
      m->_countSpecs() == 0,
      "Cannot create solver object with this steps.model.Model description "
      "object. Model must contain at least one chemical Species.");

  ArgErrLogIf(
      g->_countComps() == 0,
      "Cannot create solver object with this steps.geom.Geom geometry "
      "description object. Geometry must contain at least one Compartment.");

  std::vector<steps::wm::Comp *> comps = g->getAllComps();
  std::vector<steps::wm::Comp *>::const_iterator c_end = comps.end();
  for (std::vector<steps::wm::Comp *>::const_iterator c = comps.begin();
       c != c_end; ++c) {
    ArgErrLogIf(
        (*c)->getVol() == 0.0,
        "Cannot create solver object with this steps.geom.Geom geometry "
        "description object. All Compartments must have non-zero volume.");
  }

  // create state object, which will in turn create compdef, specdef etc
  // objects and initialise
  pStatedef = new Statedef(m, g, r);
}

////////////////////////////////////////////////////////////////////////////////

API::~API() { delete pStatedef; }

////////////////////////////////////////////////////////////////////////////////

void API::step() { NotImplErrLog(""); }

////////////////////////////////////////////////////////////////////////////////

void API::advance(double /*adv*/) { NotImplErrLog(""); }

////////////////////////////////////////////////////////////////////////////////

void API::setRk4DT(double /*dt*/) { NotImplErrLog(""); }

////////////////////////////////////////////////////////////////////////////////

double API::getRk4DT() const { NotImplErrLog(""); }

////////////////////////////////////////////////////////////////////////////////

void API::setDT(double /*dt*/) { NotImplErrLog(""); }

////////////////////////////////////////////////////////////////////////////////

double API::getDT() const { NotImplErrLog(""); }

////////////////////////////////////////////////////////////////////////////////

void API::setEfieldDT(double /*efdt*/) { NotImplErrLog(""); }

////////////////////////////////////////////////////////////////////////////////

double API::getEfieldDT() const { NotImplErrLog(""); }

////////////////////////////////////////////////////////////////////////////////

void API::setTemp(double /*temp*/) { NotImplErrLog(""); }

////////////////////////////////////////////////////////////////////////////////

double API::getTemp() const { NotImplErrLog(""); }

////////////////////////////////////////////////////////////////////////////////

double API::getA0() const { NotImplErrLog(""); }

////////////////////////////////////////////////////////////////////////////////

uint API::getNSteps() const { NotImplErrLog(""); }

////////////////////////////////////////////////////////////////////////////////

void API::setTime(double /*time*/) { NotImplErrLog(""); }
////////////////////////////////////////////////////////////////////////////////

void API::setNSteps(uint /*nsteps*/) { NotImplErrLog(""); }

////////////////////////////////////////////////////////////////////////////////
// END
