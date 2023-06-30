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
#include <set>
#include <sstream>
#include <string>

// STEPS headers.
#include "model/exocytosis.hpp"
#include "model/model.hpp"
#include "model/spec.hpp"
#include "model/vessurfsys.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

Exocytosis::Exocytosis(std::string const& id,
                       VesSurfsys* vessurfsys,
                       std::vector<Spec*> const& dep_surfacespec,
                       Raft* raft,
                       double kcst,
                       bool kiss_and_run,
                       std::map<Spec*, Spec*> const& kiss_and_run_spec_changes)
    : pID(id)
    , pModel(nullptr)
    , pVesSurfsys(vessurfsys)
    , pRaft(raft)
    , pKissAndRun(kiss_and_run)
    , pKissAndRunSpecChanges(kiss_and_run_spec_changes)
    , pKcst(kcst) {
    if (pVesSurfsys == nullptr) {
        std::ostringstream os;
        os << "No vessurfsys provided to Exocytosis initializer function";
        ArgErrLog(os.str());
    }

    if (pRaft != nullptr && pKissAndRun == true) {
        std::ostringstream os;
        os << "Cannot create Raft from kiss-n-run type Exocytosis";
        ArgErrLog(os.str());
    }

    if (pKcst < 0.0) {
        std::ostringstream os;
        os << "Exocytosis reaction constant can't be negative";
        ArgErrLog(os.str());
    }

    if (!pKissAndRun && !pKissAndRunSpecChanges.empty()) {
        CLOG(WARNING, "general_log")
            << "Exocytosis " << pID << " is not kiss-and-run so kiss-and-run "
            << "spec changes will be ignored.\n";
    }

    if (pKissAndRun) {
        std::set<Spec*> knr_specs;
        for (auto const& spec: kiss_and_run_spec_changes) {
            if (spec.first == spec.second) {
                std::ostringstream os;
                os << "In Exocytosis kiss-and-run spec changes, supplied information does not "
                      "model a change.";
                ArgErrLog(os.str());
            }
            if (!knr_specs.emplace(spec.first).second || !knr_specs.emplace(spec.second).second) {
                std::ostringstream os;
                os << "Any given species can only appear once in Exocytosis kiss-and-run spec "
                      "changes.";
                ArgErrLog(os.str());
            }
        }
    }


    pModel = pVesSurfsys->getModel();
    AssertLog(pModel != nullptr);

    for (auto const& deps: dep_surfacespec) {
        AssertLog(deps->getModel() == pModel);
        pDepSurface.emplace_back(deps);
    }

    pVesSurfsys->_handleExocytosisAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Exocytosis::~Exocytosis() {
    if (pVesSurfsys == nullptr) {
        return;
    }
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Exocytosis::_handleSelfDelete() {
    pVesSurfsys->_handleExocytosisDel(this);
    pKcst = 0.0;
    pVesSurfsys = nullptr;
    pModel = nullptr;
    pRaft = nullptr;
    pDepSurface.clear();
}

////////////////////////////////////////////////////////////////////////////////

void Exocytosis::setID(std::string const& id) {
    AssertLog(pVesSurfsys != nullptr);
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pVesSurfsys->_handleExocytosisIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Exocytosis::setKcst(double kcst) {
    assert(pVesSurfsys != nullptr);
    if (kcst < 0.0) {
        std::ostringstream os;
        os << "Exocytosis reaction constant can't be negative";
        ArgErrLog(os.str());
    }
    pKcst = kcst;
}

}  // namespace steps::model
