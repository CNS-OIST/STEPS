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

#include "vessreac.hpp"

#include "linkspec.hpp"
#include "model.hpp"
#include "spec.hpp"
#include "vessurfsys.hpp"

#include "util/error.hpp"

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

VesSReac::VesSReac(std::string const& id,
                   VesSurfsys& vessurfsys,
                   std::vector<Spec*> const& olhs,
                   std::vector<Spec*> const& slhs,
                   std::vector<Spec*> const& vlhs,
                   std::vector<LinkSpec*> const& llhs,
                   std::vector<LinkSpec*> const& lrhs,
                   std::vector<Spec*> const& vrhs,
                   std::vector<Spec*> const& srhs,
                   std::vector<Spec*> const& orhs,
                   std::vector<Spec*> const& irhs,
                   std::vector<Spec*> const& vdeps,
                   double kcst,
                   Immobilization immobilization,
                   double max_distance)
    : pID(id)
    , pModel(vessurfsys.getModel())
    , pVesSurfsys(vessurfsys)
    , pOrder(0)
    , pKcst(kcst)
    , pImmobilization(immobilization)
    , pMax_distance(max_distance) {
    if (pKcst < 0.0) {
        std::ostringstream os;
        os << "Vesicle Surface reaction constant can't be negative";
        ArgErrLog(os.str());
    }

    if (slhs.size() == 0 && pMax_distance > 0.0) {
        CLOG(WARNING, "general_log") << "WARNING: Minimum distance will be ignored "
                                        "for vesicle surface reaction "
                                     << id << "\n";
        pMax_distance = -1.0;
    }

    if (pMax_distance <= 0.0 && pMax_distance != -1.0) {
        CLOG(WARNING, "general_log") << "WARNING: Negative or zero minimum "
                                        "distance for vesicle surface reaction "
                                     << id << " will be ignored."
                                     << "\n";
    }

    if (llhs.size() != lrhs.size()) {
        std::ostringstream os;
        os << "Currently vesicle surface reactions must not change the number of "
              "link species.";
        ArgErrLog(os.str());
    }

    if (llhs.size() > 1 or lrhs.size() > 1) {
        std::ostringstream os;
        os << "Currently only one link species on lhs or rhs supported.";
        ArgErrLog(os.str());
    }

    if (olhs.size() > 0) {
        setOLHS(olhs);
    }
    if (slhs.size() > 0) {
        setSLHS(slhs);
    }
    setVLHS(vlhs);
    setLLHS(llhs);
    setLRHS(lrhs);
    setVRHS(vrhs);
    setSRHS(srhs);
    setORHS(orhs);
    setIRHS(irhs);
    setVDeps(vdeps);

    pVesSurfsys._handleVesSReacAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

VesSReac::~VesSReac() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::_handleSelfDelete() {
    pVesSurfsys._handleVesSReacDel(*this);
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::_setOrder() {
    pOrder = pOLHS.size() + pSLHS.size() + pVLHS.size() + pLLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setID(std::string const& id) {
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pVesSurfsys._handleVesSReacIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setOLHS(std::vector<Spec*> const& olhs) {
    for (auto const& ol: olhs) {
        AssertLog(&ol->getModel() == &pModel);
    }
    pOLHS = olhs;

    _setOrder();
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setVLHS(std::vector<Spec*> const& vlhs) {
    for (auto const& vl: vlhs) {
        AssertLog(&vl->getModel() == &pModel);
    }
    pVLHS = vlhs;

    _setOrder();
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setLLHS(std::vector<LinkSpec*> const& llhs) {
    for (auto const& ll: llhs) {
        AssertLog(&ll->getModel() == &pModel);
    }
    pLLHS = llhs;

    _setOrder();
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setSLHS(std::vector<Spec*> const& slhs) {
    for (auto const& sl: slhs) {
        AssertLog(&sl->getModel() == &pModel);
    }
    pSLHS = slhs;

    _setOrder();
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setVRHS(std::vector<Spec*> const& vrhs) {
    for (auto const& vr: vrhs) {
        AssertLog(&vr->getModel() == &pModel);
    }
    pVRHS = vrhs;
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setLRHS(std::vector<LinkSpec*> const& lrhs) {
    for (auto const& lr: lrhs) {
        AssertLog(&lr->getModel() == &pModel);
    }
    pLRHS = lrhs;
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setSRHS(std::vector<Spec*> const& srhs) {
    for (auto const& sr: srhs) {
        AssertLog(&sr->getModel() == &pModel);
    }
    pSRHS = srhs;
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setORHS(std::vector<Spec*> const& orhs) {
    for (auto const& ors: orhs) {
        AssertLog(&ors->getModel() == &pModel);
    }
    pORHS = orhs;
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setIRHS(std::vector<Spec*> const& irhs) {
    for (auto const& irs: irhs) {
        AssertLog(&irs->getModel() == &pModel);
    }
    pIRHS = irhs;
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setVDeps(std::vector<Spec*> const& vdeps) {
    for (auto const& vd: vdeps) {
        AssertLog(&vd->getModel() == &pModel);
    }
    pVDeps = vdeps;
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setKcst(double kcst) {
    if (kcst < 0.0) {
        std::ostringstream os;
        os << "Vesicle surface reaction constant can't be negative";
        ArgErrLog(os.str());
    }
    pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

util::flat_set<Spec*> VesSReac::getAllSpecs() const {
    util::flat_set<Spec*> specs;
    specs.insert(getOLHS().begin(), getOLHS().end());
    specs.insert(getSLHS().begin(), getSLHS().end());
    specs.insert(getVLHS().begin(), getVLHS().end());
    specs.insert(getIRHS().begin(), getIRHS().end());
    specs.insert(getSRHS().begin(), getSRHS().end());
    specs.insert(getORHS().begin(), getORHS().end());
    specs.insert(getVRHS().begin(), getVRHS().end());
    return specs;
}

////////////////////////////////////////////////////////////////////////////////

util::flat_set<LinkSpec*> VesSReac::getAllLinkSpecs() const {
    util::flat_set<LinkSpec*> linkspecs;
    linkspecs.insert(getLLHS().begin(), getLLHS().end());
    linkspecs.insert(getLRHS().begin(), getLRHS().end());
    return linkspecs;
}

}  // namespace steps::model
