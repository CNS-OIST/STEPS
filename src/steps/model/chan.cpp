/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

#include "chan.hpp"

#include "chanstate.hpp"
#include "model.hpp"

#include "util/checkid.hpp"
#include "util/error.hpp"

namespace steps::model {

using util::checkID;

////////////////////////////////////////////////////////////////////////////////

Chan::Chan(std::string const& id, Model& model)
    : pID(id)
    , pModel(model) {
    pModel._handleChanAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

Chan::~Chan() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Chan::_handleSelfDelete() {
    for (const auto& state: getAllChanStates()) {
        delete state;
    }
    pModel._handleChanDel(*this);
}

////////////////////////////////////////////////////////////////////////////////

void Chan::setID(std::string const& id) {
    if (id == pID) {
        return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pModel._handleChanIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

ChanState& Chan::getChanState(std::string const& id) const {
    auto cstate = pChanStates.find(id);

    ArgErrLogIf(cstate == pChanStates.end(),
                "Model does not contain channel state with name '" + id + "'");

    return *cstate->second;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<ChanState*> Chan::getAllChanStates() const {
    std::vector<ChanState*> cstates;
    cstates.reserve(pChanStates.size());
    for (auto const& cs: pChanStates) {
        cstates.push_back(cs.second);
    }
    return cstates;
}

////////////////////////////////////////////////////////////////////////////////

void Chan::_checkChanStateID(std::string const& id) const {
    checkID(id);

    ArgErrLogIf(pChanStates.find(id) != pChanStates.end(), "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Chan::_handleChanStateIDChange(std::string const& o, std::string const& n) {
    auto cs_old = pChanStates.find(o);
    AssertLog(cs_old != pChanStates.end());

    if (o == n) {
        return;
    }
    _checkChanStateID(n);

    ChanState* cs = cs_old->second;
    AssertLog(cs != nullptr);
    pChanStates.erase(cs_old);
    pChanStates.emplace(n, cs);
}

////////////////////////////////////////////////////////////////////////////////

void Chan::_handleChanStateAdd(ChanState& cstate) {
    AssertLog(&cstate.getChan() == this);
    _checkChanStateID(cstate.getID());
    pChanStates.emplace(cstate.getID(), &cstate);
}

////////////////////////////////////////////////////////////////////////////////

void Chan::_handleChanStateDel(ChanState& cstate) {
    AssertLog(&cstate.getChan() == this);
    pChanStates.erase(cstate.getID());
}

}  // namespace steps::model
