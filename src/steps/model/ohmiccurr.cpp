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
#include <cassert>
#include <sstream>
#include <string>
#include <iostream>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/model/model.hpp"
#include "steps/model/surfsys.hpp"
#include "steps/model/ohmiccurr.hpp"
#include "steps/model/chanstate.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

////////////////////////////////////////////////////////////////////////////////

OhmicCurr::OhmicCurr(string const & id, Surfsys * surfsys,
          ChanState * chanstate, double erev, double g)
: pID(id)
, pModel(0)
, pSurfsys(surfsys)
, pChanState(chanstate)
, pERev(erev)
, pG(g)
{
    if (pSurfsys == 0)
    {
        ostringstream os;
        os << "No surfsys provided to OhmicCurr initializer function";
        ArgErrLog(os.str());
    }
    if (pChanState == 0)
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
    AssertLog(pModel != 0);

    pSurfsys->_handleOhmicCurrAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurr::~OhmicCurr(void)
{
    if (pSurfsys == 0) return;
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void OhmicCurr::setID(string const & id)
{
    AssertLog(pSurfsys != 0);
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
    AssertLog(chanstate != 0);
    pChanState = chanstate;
}

////////////////////////////////////////////////////////////////////////////////

void OhmicCurr::setERev(double erev)
{
    AssertLog(pSurfsys != 0);
    pERev = erev;
}

////////////////////////////////////////////////////////////////////////////////

void OhmicCurr::setG(double g)
{
    AssertLog(pSurfsys != 0);
    if(g < 0.0)
    {
        ostringstream os;
        os << "Conductance provided to OhmicCurr::setG function can't be negative";
        ArgErrLog(os.str());
    }
    pG = g;
}

////////////////////////////////////////////////////////////////////////////////

void OhmicCurr::_handleSelfDelete(void)
{
    pSurfsys->_handleOhmicCurrDel(this);
    pG = 0.0;
    pERev = 0;
    pSurfsys = 0;
    pModel = 0;
}

////////////////////////////////////////////////////////////////////////////////

// END

