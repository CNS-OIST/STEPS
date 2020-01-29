/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
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
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/model/chanstate.hpp"
#include "steps/model/model.hpp"
#include "steps/model/ohmiccurr.hpp"
#include "steps/model/surfsys.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

////////////////////////////////////////////////////////////////////////////////

OhmicCurr::OhmicCurr(string const & id, Surfsys * surfsys,
          ChanState * chanstate, double erev, double g)
: pID(id)
, pModel(nullptr)
, pSurfsys(surfsys)
, pChanState(chanstate)
, pERev(erev)
, pG(g)
{
    if (pSurfsys == nullptr)
    {
        ostringstream os;
        os << "No surfsys provided to OhmicCurr initializer function";
        ArgErrLog(os.str());
    }
    if (pChanState == nullptr)
    {
        ostringstream os;
        os << "No channel state provided to OhmicCurr initializer function";
        ArgErrLog(os.str());
    }
    if (pG < 0.0)
    {
        ostringstream os;
        os << "Channel conductance can't be negative";
        ArgErrLog(os.str());
    }

    pModel = pSurfsys->getModel();
    AssertLog(pModel != nullptr);

    pSurfsys->_handleOhmicCurrAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurr::~OhmicCurr()
{
    if (pSurfsys == nullptr) { return;
}
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void OhmicCurr::setID(string const & id)
{
    AssertLog(pSurfsys != nullptr);
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pSurfsys->_handleOhmicCurrIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void OhmicCurr::setChanState(ChanState * chanstate)
{
    AssertLog(chanstate != nullptr);
    pChanState = chanstate;
}

////////////////////////////////////////////////////////////////////////////////

void OhmicCurr::setERev(double erev)
{
    AssertLog(pSurfsys != nullptr);
    pERev = erev;
}

////////////////////////////////////////////////////////////////////////////////

void OhmicCurr::setG(double g)
{
    AssertLog(pSurfsys != nullptr);
    if(g < 0.0)
    {
        ostringstream os;
        os << "Conductance provided to OhmicCurr::setG function can't be negative";
        ArgErrLog(os.str());
    }
    pG = g;
}

////////////////////////////////////////////////////////////////////////////////

void OhmicCurr::_handleSelfDelete()
{
    pSurfsys->_handleOhmicCurrDel(this);
    pG = 0.0;
    pERev = 0;
    pSurfsys = nullptr;
    pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

// END

