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
#include "model/linkspec.hpp"
#include "model/model.hpp"
#include "model/spec.hpp"
#include "model/vessreac.hpp"
#include "model/vessurfsys.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

VesSReac::VesSReac(std::string const& id,
                   VesSurfsys* vessurfsys,
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
    , pModel(nullptr)
    , pVesSurfsys(vessurfsys)
    , pOrder(0)
    , pKcst(kcst)
    , pImmobilization(immobilization)
    , pMax_distance(max_distance) {
    if (pVesSurfsys == nullptr) {
        std::ostringstream os;
        os << "No vessurfsys provided to VesSReac initializer function";
        ArgErrLog(os.str());
    }

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

    pModel = pVesSurfsys->getModel();
    AssertLog(pModel != nullptr);

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

    pVesSurfsys->_handleVesSReacAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

VesSReac::~VesSReac() {
    if (pVesSurfsys == nullptr) {
        return;
    }
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::_handleSelfDelete() {
    pVesSurfsys->_handleVesSReacDel(this);
    pKcst = 0.0;
    pOrder = 0;
    pIRHS.clear();
    pORHS.clear();
    pSRHS.clear();
    pVRHS.clear();
    pVLHS.clear();
    pSLHS.clear();
    pOLHS.clear();
    pVDeps.clear();
    pVesSurfsys = nullptr;
    pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::_setOrder() {
    pOrder = pOLHS.size() + pSLHS.size() + pVLHS.size() + pLLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setID(std::string const& id) {
    AssertLog(pVesSurfsys != nullptr);

    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pVesSurfsys->_handleVesSReacIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setOLHS(std::vector<Spec*> const& olhs) {
    AssertLog(pVesSurfsys != nullptr);

    pOLHS.clear();
    for (auto const& ol: olhs) {
        AssertLog(ol->getModel() == pModel);
        pOLHS.emplace_back(ol);
    }

    _setOrder();
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setVLHS(std::vector<Spec*> const& vlhs) {
    AssertLog(pVesSurfsys != nullptr);

    pVLHS.clear();
    for (auto const& vl: vlhs) {
        AssertLog(vl->getModel() == pModel);
        pVLHS.emplace_back(vl);
    }

    _setOrder();
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setLLHS(std::vector<LinkSpec*> const& llhs) {
    AssertLog(pVesSurfsys != nullptr);

    pLLHS.clear();
    for (auto const& ll: llhs) {
        AssertLog(ll->getModel() == pModel);
        pLLHS.emplace_back(ll);
    }

    _setOrder();
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setSLHS(std::vector<Spec*> const& slhs) {
    AssertLog(pVesSurfsys != nullptr);

    pSLHS.clear();
    for (auto const& sl: slhs) {
        AssertLog(sl->getModel() == pModel);
        pSLHS.emplace_back(sl);
    }

    _setOrder();
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setVRHS(std::vector<Spec*> const& vrhs) {
    AssertLog(pVesSurfsys != nullptr);
    pVRHS.clear();
    for (auto const& vr: vrhs) {
        AssertLog(vr->getModel() == pModel);
        pVRHS.emplace_back(vr);
    }
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setLRHS(std::vector<LinkSpec*> const& lrhs) {
    AssertLog(pVesSurfsys != nullptr);
    pLRHS.clear();
    for (auto const& lr: lrhs) {
        AssertLog(lr->getModel() == pModel);
        pLRHS.emplace_back(lr);
    }
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setSRHS(std::vector<Spec*> const& srhs) {
    AssertLog(pVesSurfsys != nullptr);
    pSRHS.clear();
    for (auto const& sr: srhs) {
        AssertLog(sr->getModel() == pModel);
        pSRHS.emplace_back(sr);
    }
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setORHS(std::vector<Spec*> const& orhs) {
    AssertLog(pVesSurfsys != nullptr);
    pORHS.clear();
    for (auto const& ors: orhs) {
        AssertLog(ors->getModel() == pModel);
        pORHS.emplace_back(ors);
    }
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setIRHS(std::vector<Spec*> const& irhs) {
    AssertLog(pVesSurfsys != nullptr);
    pIRHS.clear();
    for (auto const& irs: irhs) {
        AssertLog(irs->getModel() == pModel);
        pIRHS.emplace_back(irs);
    }
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setVDeps(std::vector<Spec*> const& vdeps) {
    AssertLog(pVesSurfsys != nullptr);
    pVDeps.clear();
    for (auto const& vd: vdeps) {
        AssertLog(vd->getModel() == pModel);
        pVDeps.emplace_back(vd);
    }
}

////////////////////////////////////////////////////////////////////////////////

void VesSReac::setKcst(double kcst) {
    AssertLog(pVesSurfsys != nullptr);

    if (kcst < 0.0) {
        std::ostringstream os;
        os << "Vesicle surface reaction constant can't be negative";
        ArgErrLog(os.str());
    }
    pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec*> VesSReac::getAllSpecs() const {
    SpecPVec specs = SpecPVec();
    bool first_occ = true;

    for (auto const& ol: getOLHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == ol) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) {
            specs.emplace_back(ol);
        }
    }

    for (auto const& sl: getSLHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == sl) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) {
            specs.emplace_back(sl);
        }
    }

    for (auto const& vl: getVLHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == vl) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) {
            specs.emplace_back(vl);
        }
    }

    for (auto const& ir: getIRHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == ir) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) {
            specs.emplace_back(ir);
        }
    }

    for (auto const& sr: getSRHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == sr) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) {
            specs.emplace_back(sr);
        }
    }

    for (auto const& ors: getORHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == ors) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) {
            specs.emplace_back(ors);
        }
    }

    for (auto const& vr: getVRHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == vr) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) {
            specs.emplace_back(vr);
        }
    }

    return specs;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<LinkSpec*> VesSReac::getAllLinkSpecs() const {
    LinkSpecPVec specs = LinkSpecPVec();
    bool first_occ = true;

    for (auto const& ll: getLLHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == ll) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) {
            specs.emplace_back(ll);
        }
    }

    for (auto const& lr: getLRHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == lr) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) {
            specs.emplace_back(lr);
        }
    }

    return specs;
}

}  // namespace steps::model
