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

#include "model/raftendocytosis.hpp"

#include "raftsys.hpp"
#include "spec.hpp"

#include "util/error.hpp"

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

RaftEndocytosis::RaftEndocytosis(std::string const& id,
                                 Raftsys& raftsys,
                                 Vesicle* irhs,
                                 Vesicle* orhs,
                                 std::vector<Spec*> const& dep_surfacespec,
                                 double kcst)
    : pID(id)
    , pModel(raftsys.getModel())
    , pRaftsys(raftsys)
    , pKcst(kcst) {
    if (irhs == nullptr and orhs == nullptr) {
        std::ostringstream os;
        os << "No associated Vesicle provided to RaftEndocytosis initializer "
              "function";
        ArgErrLog(os.str());
    }

    if (irhs != nullptr and orhs != nullptr) {
        std::ostringstream os;
        os << "Both inner and outer compartment Vesicles defined in "
              "RaftEndocytosis initializer function";
        ArgErrLog(os.str());
    }

    if (irhs == nullptr) {
        pVesicle = orhs;
        pInner = false;
    } else {
        pVesicle = irhs;
        pInner = true;
    }

    if (pKcst < 0.0) {
        std::ostringstream os;
        os << "RaftEndocytosis reaction constant can't be negative";
        ArgErrLog(os.str());
    }

    pDepSurface.reserve(dep_surfacespec.size());
    for (auto const& deps: dep_surfacespec) {
        AssertLog(&deps->getModel() == &pModel);
        pDepSurface.push_back(deps);
    }

    pRaftsys._handleRaftEndocytosisAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

RaftEndocytosis::~RaftEndocytosis() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosis::_handleSelfDelete() {
    pRaftsys._handleRaftEndocytosisDel(*this);
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosis::setID(std::string const& id) {
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pRaftsys._handleRaftEndocytosisIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosis::setKcst(double kcst) {
    if (kcst < 0.0) {
        std::ostringstream os;
        os << "RaftEndocytosis reaction constant can't be negative";
        ArgErrLog(os.str());
    }

    pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

util::flat_set<Spec*> RaftEndocytosis::getAllSpecs() const {
    util::flat_set<Spec*> specs;
    const auto& specdeps = getSpecDeps();
    specs.insert(specdeps.begin(), specdeps.end());
    return specs;
}

}  // namespace steps::model
