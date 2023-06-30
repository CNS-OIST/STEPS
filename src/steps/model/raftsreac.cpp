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
#include "model/raftsreac.hpp"
#include "model/raftsys.hpp"
#include "model/spec.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

RaftSReac::RaftSReac(std::string const& id,
                     Raftsys* raftsys,
                     std::vector<Spec*> const& ilhs,
                     std::vector<Spec*> const& olhs,
                     std::vector<Spec*> const& slhs,
                     std::vector<Spec*> const& rslhs,
                     std::vector<Spec*> const& rsrhs,
                     std::vector<Spec*> const& srhs,
                     std::vector<Spec*> const& orhs,
                     std::vector<Spec*> const& irhs,
                     std::vector<Spec*> const& rsdeps,
                     std::vector<Spec*> const& anti_rsdeps,
                     double kcst,
                     Immobilization immobilization)
    : pID(id)
    , pModel(nullptr)
    , pRaftsys(raftsys)
    , pImmobilization(immobilization)
    , pOuter(true)
    , pOrder(0)
    , pKcst(kcst) {
    if (pRaftsys == nullptr) {
        std::ostringstream os;
        os << "No raftsys provided to RaftSReac initializer function";
        ArgErrLog(os.str());
    }

    if (pKcst < 0.0) {
        std::ostringstream os;
        os << "Raft surface reaction constant can't be negative";
        ArgErrLog(os.str());
    }

    // Can't have species on the lhs in the inner and outer compartment
    if (olhs.size() != 0 && ilhs.size() != 0) {
        std::ostringstream os;
        os << "Reactant species cannot consist of both ilhs and olhs species.";
        ArgErrLog(os.str());
    }

    pModel = pRaftsys->getModel();
    AssertLog(pModel != nullptr);

    if (!olhs.empty()) {
        setOLHS(olhs);
    }
    if (!ilhs.empty()) {
        setILHS(ilhs);
    }
    setSLHS(slhs);
    setRsLHS(rslhs);
    setRsRHS(rsrhs);
    setSRHS(srhs);
    setORHS(orhs);
    setIRHS(irhs);
    setRsDeps(rsdeps);
    setAntiRsDeps(anti_rsdeps);

    pRaftsys->_handleRaftSReacAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

RaftSReac::~RaftSReac() {
    if (pRaftsys == nullptr) {
        return;
    }
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::_handleSelfDelete() {
    pRaftsys->_handleRaftSReacDel(this);
    pKcst = 0.0;
    pOrder = 0;
    pIRHS.clear();
    pORHS.clear();
    pSRHS.clear();
    pRsRHS.clear();
    pRsLHS.clear();
    pSLHS.clear();
    pOLHS.clear();
    pILHS.clear();
    pRsDeps.clear();
    pRaftsys = nullptr;
    pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::setID(std::string const& id) {
    AssertLog(pRaftsys != nullptr);
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pRaftsys->_handleRaftSReacIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::setILHS(std::vector<Spec*> const& ilhs) {
    AssertLog(pRaftsys != nullptr);

    if (pOLHS.size() != 0) {
        CLOG(WARNING, "general_log") << "Removing outer compartment species from "
                                        "lhs stoichiometry for RaftSReac "
                                     << getID() << ".\n";
        pOLHS.clear();
    }

    pILHS.clear();

    for (auto const& il: ilhs) {
        AssertLog(il->getModel() == pModel);
        pILHS.push_back(il);
    }
    pOuter = false;
    pOrder = pILHS.size() + pRsLHS.size() + pSLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::setOLHS(std::vector<Spec*> const& olhs) {
    AssertLog(pRaftsys != nullptr);

    if (pILHS.size() != 0) {
        CLOG(WARNING, "general_log") << "Removing inner compartment species from "
                                        "lhs stoichiometry for RaftSReac "
                                     << getID() << ".\n";
        pILHS.clear();
    }

    pOLHS.clear();

    for (auto const& ol: olhs) {
        AssertLog(ol->getModel() == pModel);
        pOLHS.push_back(ol);
    }
    pOuter = true;
    pOrder = pOLHS.size() + pSLHS.size() + pRsLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::setRsLHS(std::vector<Spec*> const& rslhs) {
    AssertLog(pRaftsys != nullptr);

    pRsLHS.clear();

    for (auto const& rsl: rslhs) {
        AssertLog(rsl->getModel() == pModel);
        pRsLHS.push_back(rsl);
    }

    if (pOuter) {
        pOrder = pOLHS.size() + pSLHS.size() + pRsLHS.size();
    } else {
        pOrder = pILHS.size() + pRsLHS.size() + pSLHS.size();
    }
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::setSLHS(std::vector<Spec*> const& slhs) {
    AssertLog(pRaftsys != nullptr);

    pSLHS.clear();

    for (auto const& sl: slhs) {
        AssertLog(sl->getModel() == pModel);
        pSLHS.push_back(sl);
    }

    if (pOuter) {
        pOrder = pOLHS.size() + pSLHS.size() + pRsLHS.size();
    } else {
        pOrder = pILHS.size() + pRsLHS.size() + pSLHS.size();
    }
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::setRsRHS(std::vector<Spec*> const& rsrhs) {
    AssertLog(pRaftsys != nullptr);

    pRsRHS.clear();

    for (auto const& rsr: rsrhs) {
        AssertLog(rsr->getModel() == pModel);
        pRsRHS.push_back(rsr);
    }
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::setSRHS(std::vector<Spec*> const& srhs) {
    AssertLog(pRaftsys != nullptr);

    pSRHS.clear();

    for (auto const& sr: srhs) {
        AssertLog(sr->getModel() == pModel);
        pSRHS.push_back(sr);
    }
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::setORHS(std::vector<Spec*> const& orhs) {
    AssertLog(pRaftsys != nullptr);

    pORHS.clear();

    for (auto const& ors: orhs) {
        AssertLog(ors->getModel() == pModel);
        pORHS.push_back(ors);
    }
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::setIRHS(std::vector<Spec*> const& irhs) {
    AssertLog(pRaftsys != nullptr);

    pIRHS.clear();

    for (auto const& irs: irhs) {
        AssertLog(irs->getModel() == pModel);
        pIRHS.push_back(irs);
    }
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::setRsDeps(std::vector<Spec*> const& rsdeps) {
    AssertLog(pRaftsys != nullptr);

    pRsDeps.clear();

    for (auto const& rsd: rsdeps) {
        AssertLog(rsd->getModel() == pModel);
        pRsDeps.push_back(rsd);
    }
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::setAntiRsDeps(std::vector<Spec*> const& anti_rsdeps) {
    AssertLog(pRaftsys != nullptr);

    pAntiRsDeps.clear();

    for (auto const& arsd: anti_rsdeps) {
        AssertLog(arsd->getModel() == pModel);
        pAntiRsDeps.push_back(arsd);
    }
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::setKcst(double kcst) {
    AssertLog(pRaftsys != nullptr);

    if (kcst < 0.0) {
        std::ostringstream os;
        os << "Surface reaction constant can't be negative";
        ArgErrLog(os.str());
    }
    pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec*> RaftSReac::getAllSpecs() const {
    SpecPVec specs = SpecPVec();
    bool first_occ = true;

    SpecPVec ilhs = getILHS();

    for (auto const& il: ilhs) {
        first_occ = true;

        for (auto const& s: specs) {
            if (s == il) {
                first_occ = false;
                break;
            }
        }
        if (first_occ == true) {
            specs.push_back(il);
        }
    }

    SpecPVec olhs = getOLHS();

    for (auto const& ol: olhs) {
        first_occ = true;

        for (auto const& s: specs) {
            if (s == ol) {
                first_occ = false;
                break;
            }
        }
        if (first_occ == true) {
            specs.push_back(ol);
        }
    }

    SpecPVec rslhs = getRsLHS();

    for (auto const& rsl: rslhs) {
        first_occ = true;

        for (auto const& s: specs) {
            if (s == rsl) {
                first_occ = false;
                break;
            }
        }
        if (first_occ == true) {
            specs.push_back(rsl);
        }
    }

    SpecPVec slhs = getSLHS();

    for (auto const& sl: slhs) {
        first_occ = true;

        for (auto const& s: specs) {
            if (s == sl) {
                first_occ = false;
                break;
            }
        }
        if (first_occ == true) {
            specs.push_back(sl);
        }
    }

    SpecPVec rsrhs = getRsRHS();

    for (auto const& rsr: rsrhs) {
        first_occ = true;

        for (auto const& s: specs) {
            if (s == rsr) {
                first_occ = false;
                break;
            }
        }
        if (first_occ == true) {
            specs.push_back(rsr);
        }
    }

    SpecPVec srhs = getSRHS();

    for (auto const& sr: srhs) {
        first_occ = true;

        for (auto const& s: specs) {
            if (s == sr) {
                first_occ = false;
                break;
            }
        }
        if (first_occ == true) {
            specs.push_back(sr);
        }
    }

    SpecPVec orhs = getORHS();

    for (auto const& ors: orhs) {
        first_occ = true;

        for (auto const& s: specs) {
            if (s == ors) {
                first_occ = false;
                break;
            }
        }
        if (first_occ == true) {
            specs.push_back(ors);
        }
    }

    SpecPVec irhs = getIRHS();

    for (auto const& irs: irhs) {
        first_occ = true;

        for (auto const& s: specs) {
            if (s == irs) {
                first_occ = false;
                break;
            }
        }
        if (first_occ == true) {
            specs.push_back(irs);
        }
    }

    return specs;
}

}  // namespace steps::model
