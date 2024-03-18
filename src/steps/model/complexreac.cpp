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

#include "model/complexreac.hpp"

#include "model/model.hpp"
#include "model/spec.hpp"
#include "model/volsys.hpp"
#include "util/error.hpp"

namespace steps::model {

ComplexReac::ComplexReac(std::string const& id,
                         Volsys& volsys,
                         std::vector<Spec*> const& lhs,
                         std::vector<Spec*> const& rhs,
                         std::vector<ComplexEvent*> const& compEvs,
                         double kcst)
    : pID(id)
    , pModel(volsys.getModel())
    , pVolsys(volsys)
    , pLHS(lhs)
    , pRHS(rhs) {
    setKcst(kcst);

    for (auto const& l: pLHS) {
        AssertLog(&l->getModel() == &pModel);
    }
    for (auto const& r: pRHS) {
        AssertLog(&r->getModel() == &pModel);
    }

    pOrder = pLHS.size();

    for (auto* ev: compEvs) {
        if (dynamic_cast<ComplexUpdateEvent*>(ev)) {
            pCompUPD.push_back(dynamic_cast<ComplexUpdateEvent*>(ev));
            ++pOrder;
        } else if (dynamic_cast<ComplexDeleteEvent*>(ev)) {
            pCompDEL.push_back(dynamic_cast<ComplexDeleteEvent*>(ev));
            ++pOrder;
        } else if (dynamic_cast<ComplexCreateEvent*>(ev)) {
            pCompCRE.push_back(dynamic_cast<ComplexCreateEvent*>(ev));
        }
    }

    pVolsys._handleComplexReacAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

void ComplexReac::setKcst(double kcst) {
    ArgErrLogIf(kcst < 0.0, "Reaction constant can't be negative");
    pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

util::flat_set<Spec*> ComplexReac::getAllSpecs() const {
    util::flat_set<Spec*> specSet;
    specSet.insert(pLHS.begin(), pLHS.end());
    specSet.insert(pRHS.begin(), pRHS.end());
    return specSet;
}

}  // namespace steps::model
