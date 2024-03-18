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

#include "model/endocytosis.hpp"

#include "model.hpp"
#include "spec.hpp"
#include "surfsys.hpp"

#include "util/error.hpp"

namespace steps::model {

Endocytosis::Endocytosis(std::string const& id,
                         Surfsys& surfsys,
                         Vesicle* irhs,
                         Vesicle* orhs,
                         std::vector<Spec*> const& dep_surfacespec,
                         double kcst)
    : pID(id)
    , pModel(surfsys.getModel())
    , pSurfsys(surfsys)
    , pVesicle()
    , pKcst(kcst) {
    if (irhs == nullptr and orhs == nullptr) {
        std::ostringstream os;
        os << "No associated Vesicle provided to Endocytosis initializer function";
        ArgErrLog(os.str());
    }

    if (irhs != nullptr and orhs != nullptr) {
        std::ostringstream os;
        os << "Both inner and outer compartment Vesicles defined in Endocytosis "
              "initializer function";
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
        os << "Endocytosis reaction constant can't be negative";
        ArgErrLog(os.str());
    }

    for (auto const& deps: dep_surfacespec) {
        AssertLog(&deps->getModel() == &pModel);
        pDepSurface.emplace_back(deps);
    }

    pSurfsys._handleEndocytosisAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

Endocytosis::~Endocytosis() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Endocytosis::_handleSelfDelete() {
    pSurfsys._handleEndocytosisDel(*this);
    pKcst = 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void Endocytosis::setID(std::string const& id) {
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pSurfsys._handleEndocytosisIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Endocytosis::setKcst(double kcst) {
    if (kcst < 0.0) {
        std::ostringstream os;
        os << "Endocytosis reaction constant can't be negative";
        ArgErrLog(os.str());
    }
    pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

util::flat_set<Spec*> Endocytosis::getAllSpecs() const {
    util::flat_set<Spec*> specs;
    specs.insert(getSpecDeps().begin(), getSpecDeps().end());
    return specs;
}

}  // namespace steps::model
