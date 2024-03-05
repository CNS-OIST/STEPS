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

#include "vesbind.hpp"

#include "linkspec.hpp"
#include "model.hpp"
#include "spec.hpp"
#include "vesicle.hpp"
#include "volsys.hpp"

#include "util/error.hpp"

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

VesBind::VesBind(std::string const& id,
                 Volsys& volsys,
                 Vesicle& vesicle1,
                 Spec& reactant1,
                 Vesicle& vesicle2,
                 Spec& reactant2,
                 LinkSpec& product1,
                 LinkSpec& product2,
                 double length_max,
                 double length_min,
                 std::vector<Spec*> const& vdeps1,
                 std::vector<Spec*> const& vdeps2,
                 std::vector<LinkSpec*> const& ldeps1,
                 std::vector<LinkSpec*> const& ldeps2,
                 double kcst,
                 Immobilization immobilization)
    : pID(id)
    , pModel(volsys.getModel())
    , pVolsys(volsys)
    , pKcst(kcst)
    , pImmobilization(immobilization) {
    if (pKcst < 0.0) {
        std::ostringstream os;
        os << "Vesicle binding constant can't be negative";
        ArgErrLog(os.str());
    }

    if (pImmobilization == MOBILIZING) {
        std::ostringstream os;
        os << "Unsupported immobilization flag. A VesBind event cannot mobilize vesicles.";
        ArgErrLog(os.str());
    }

    if (length_min >= length_max || length_min < 0.0) {
        std::ostringstream os;
        os << "Maximum length must be greater than minimum length, and neither "
              "length can be negative.";
        ArgErrLog(os.str());
    }

    AssertLog(&reactant1.getModel() == &pModel);
    AssertLog(&reactant2.getModel() == &pModel);

    pReactants1 = {&vesicle1, &reactant1};
    pReactants2 = {&vesicle2, &reactant2};

    AssertLog(&product1.getModel() == &pModel);
    AssertLog(&product2.getModel() == &pModel);

    pProducts1 = {&vesicle1, &product1};
    pProducts2 = {&vesicle2, &product2};

    pLength_max = length_max;
    pLength_min = length_min;

    setVDeps1(vdeps1);
    setVDeps2(vdeps2);
    setLDeps1(ldeps1);
    setLDeps2(ldeps2);

    pVolsys._handleVesBindAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

VesBind::~VesBind() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void VesBind::_handleSelfDelete() {
    pVolsys._handleVesBindDel(*this);
}

////////////////////////////////////////////////////////////////////////////////

void VesBind::setID(std::string const& id) {
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pVolsys._handleVesBindIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void VesBind::setKcst(double kcst) {
    if (kcst < 0.0) {
        std::ostringstream os;
        os << "Vesicle binding rate constant can't be negative";
        ArgErrLog(os.str());
    }
    pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

void VesBind::setVDeps1(std::vector<Spec*> const& vdeps) {
    for (auto const& vd: vdeps) {
        AssertLog(&vd->getModel() == &pModel);
    }
    pVDeps1 = vdeps;
}

////////////////////////////////////////////////////////////////////////////////

void VesBind::setVDeps2(std::vector<Spec*> const& vdeps) {
    for (auto const& vd: vdeps) {
        AssertLog(&vd->getModel() == &pModel);
    }
    pVDeps2 = vdeps;
}

////////////////////////////////////////////////////////////////////////////////

void VesBind::setLDeps1(std::vector<LinkSpec*> const& ldeps) {
    for (auto const& ld: ldeps) {
        AssertLog(&ld->getModel() == &pModel);
    }
    pLDeps1 = ldeps;
}

////////////////////////////////////////////////////////////////////////////////

void VesBind::setLDeps2(std::vector<LinkSpec*> const& ldeps) {
    for (auto const& ld: ldeps) {
        AssertLog(&ld->getModel() == &pModel);
    }
    pLDeps2 = ldeps;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec*> VesBind::getAllSpecs() const {
    std::vector<Spec*> specs;

    specs.emplace_back(std::get<1>(pReactants1));

    Spec* spec2 = std::get<1>(pReactants2);

    if (specs[0] != spec2) {
        specs.emplace_back(spec2);
    }
    return specs;
}

}  // namespace steps::model
