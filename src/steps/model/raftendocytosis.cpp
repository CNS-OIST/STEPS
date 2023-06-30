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

// STL headers.
#include <iostream>
#include <sstream>
#include <string>

// STEPS headers.
#include "model/model.hpp"
#include "model/raftendocytosis.hpp"
#include "model/raftsys.hpp"
#include "model/spec.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

RaftEndocytosis::RaftEndocytosis(std::string const& id,
                                 Raftsys* raftsys,
                                 Vesicle* irhs,
                                 Vesicle* orhs,
                                 std::vector<Spec*> const& dep_surfacespec,
                                 double kcst)
    : pID(id)
    , pModel(nullptr)
    , pRaftsys(raftsys)
    , pVesicle()
    , pKcst(kcst) {
    if (pRaftsys == nullptr) {
        std::ostringstream os;
        os << "No surfsys provided to RaftEndocytosis initializer function";
        ArgErrLog(os.str());
    }

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

    pModel = pRaftsys->getModel();
    AssertLog(pModel != nullptr);

    pDepSurface.reserve(dep_surfacespec.size());
    for (auto const& deps: dep_surfacespec) {
        AssertLog(deps->getModel() == pModel);
        pDepSurface.push_back(deps);
    }

    pRaftsys->_handleRaftEndocytosisAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

RaftEndocytosis::~RaftEndocytosis() {
    if (pRaftsys == nullptr) {
        return;
    }
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosis::_handleSelfDelete() {
    pRaftsys->_handleRaftEndocytosisDel(this);
    pKcst = 0.0;
    pRaftsys = nullptr;
    pModel = nullptr;
    pVesicle = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosis::setID(std::string const& id) {
    AssertLog(pRaftsys != nullptr);
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pRaftsys->_handleRaftEndocytosisIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosis::setKcst(double kcst) {
    AssertLog(pRaftsys != nullptr);

    if (kcst < 0.0) {
        std::ostringstream os;
        os << "RaftEndocytosis reaction constant can't be negative";
        ArgErrLog(os.str());
    }

    pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec*> RaftEndocytosis::getAllSpecs() const {
    SpecPVec specs = SpecPVec();
    bool first_occ = true;

    SpecPVec specdeps = getSpecDeps();

    for (auto const& sd: specdeps) {
        first_occ = true;

        for (auto const& s: specs) {
            if (s == sd) {
                first_occ = false;
                break;
            }
        }
        if (first_occ == true) {
            specs.emplace_back(sd);
        }
    }

    return specs;
}

}  // namespace steps::model
