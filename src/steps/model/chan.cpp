/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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
#include <map>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/model/model.hpp"
#include "steps/model/spec.hpp"
#include "steps/model/chan.hpp"
#include "steps/model/chanstate.hpp"
#include "steps/util/checkid.hpp"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

using steps::util::checkID;

////////////////////////////////////////////////////////////////////////////////

Chan::Chan(string const & id, Model * model)
: pID(id)
, pModel(model)
, pChanStates()
{
    if (pModel == 0)
    {
        ostringstream os;
        os << "No model provided to Channel initializer function.";
        throw steps::ArgErr(os.str());
    }
    pModel->_handleChanAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Chan::~Chan(void)
{
    if (pModel == 0) return;
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Chan::_handleSelfDelete(void)
{
    std::vector<steps::model::ChanState *> allstates = getAllChanStates();
    ChanStatePVecCI cstate_end = allstates.end();
    for(ChanStatePVecCI cstate = allstates.begin(); cstate != cstate_end; ++cstate)
    {
        delete(*cstate);
    }

    pModel->_handleChanDel(this);
    pChanStates.clear();
    pModel = 0;
}

////////////////////////////////////////////////////////////////////////////////

void Chan::setID(string const & id)
{
    assert(pModel != 0);
    if (id == pID) return;
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pModel->_handleChanIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

ChanState * Chan::getChanState(string const & id) const
{
    ChanStatePMapCI cstate = pChanStates.find(id);
    if (cstate == pChanStates.end())
    {
        ostringstream os;
        os << "Model does not contain channel state with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
    assert(cstate->second != 0);
    return cstate->second;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<ChanState *> Chan::getAllChanStates(void) const
{
    ChanStatePVec cstates = ChanStatePVec();
    ChanStatePMapCI cs_end = pChanStates.end();
    for (ChanStatePMapCI cs = pChanStates.begin(); cs != cs_end; ++cs)
    {
        cstates.push_back(cs->second);
    }
    return cstates;
}

////////////////////////////////////////////////////////////////////////////////

void Chan::_checkChanStateID(string const & id) const
{
    checkID(id);
    if (pChanStates.find(id) != pChanStates.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Chan::_handleChanStateIDChange(string const & o, string const & n)
{
    ChanStatePMapCI cs_old = pChanStates.find(o);
    assert(cs_old != pChanStates.end());

    if(o==n) return;
    _checkChanStateID(n);

    ChanState * cs = cs_old->second;
    assert(cs != 0);
    pChanStates.erase(cs->getID());
    pChanStates.insert(ChanStatePMap::value_type(n,cs));
}

////////////////////////////////////////////////////////////////////////////////

void Chan::_handleChanStateAdd(ChanState * cstate)
{
    assert(cstate->getChan() == this);
    _checkChanStateID(cstate->getID());
    pChanStates.insert(ChanStatePMap::value_type(cstate->getID(), cstate));
}

////////////////////////////////////////////////////////////////////////////////

void Chan::_handleChanStateDel(ChanState * cstate)
{
    assert(cstate->getChan() == this);
    pChanStates.erase(cstate->getID());
}

////////////////////////////////////////////////////////////////////////////////

// END
