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

#include "model/vessurfsys.hpp"

// STL headers.
#include <set>
#include <sstream>

// STEPS headers.
#include "model/chan.hpp"
#include "model/chanstate.hpp"
#include "model/exocytosis.hpp"
#include "model/linkspec.hpp"
#include "model/model.hpp"
#include "model/spec.hpp"
#include "model/vessdiff.hpp"
#include "model/vessreac.hpp"
#include "util/checkid.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace steps::model {

using util::checkID;

////////////////////////////////////////////////////////////////////////////////

VesSurfsys::VesSurfsys(std::string const& id, Model* model)
    : pID(id)
    , pModel(model) {
    if (pModel == nullptr) {
        std::ostringstream os;
        os << "No model provided to VesSurfsys initializer function";
        ArgErrLog(os.str());
    }

    pModel->_handleVesSurfsysAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

VesSurfsys::~VesSurfsys() {
    if (pModel == nullptr) {
        return;
    }
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::setID(std::string const& id) {
    AssertLog(pModel != nullptr);

    if (id == pID) {
        return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pModel->_handleVesSurfsysIDChange(pID, id);

    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_handleSelfDelete() {
    for (auto vessreac: getAllVesSReacs()) {
        delete vessreac;
    }

    for (auto vessdiff: getAllVesSDiffs()) {
        delete vessdiff;
    }

    for (auto exo: getAllExocytosis()) {
        delete exo;
    }

    pModel->_handleVesSurfsysDel(this);

    pVesSReacs.clear();
    pVesSDiffs.clear();
    pExocytosis.clear();
    pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

VesSReac* VesSurfsys::getVesSReac(std::string const& id) const {
    auto vessreac = pVesSReacs.find(id);
    if (vessreac == pVesSReacs.end()) {
        std::ostringstream os;
        os << "Model does not contain vesicle surface "
              "reaction with name '"
           << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(vessreac->second != nullptr);

    return vessreac->second;
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::delVesSReac(std::string const& id) const {
    VesSReac* vessreac = getVesSReac(id);
    delete (vessreac);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<VesSReac*> VesSurfsys::getAllVesSReacs() const {
    VesSReacPVec vessreacs = VesSReacPVec();

    for (auto const& vsr: pVesSReacs) {
        vessreacs.emplace_back(vsr.second);
    }
    return vessreacs;
}

////////////////////////////////////////////////////////////////////////////////

VesSDiff* VesSurfsys::getVesSDiff(std::string const& id) const {
    auto vessdiff = pVesSDiffs.find(id);
    if (vessdiff == pVesSDiffs.end()) {
        std::ostringstream os;
        os << "Model does not contain vesicle surface diffusion with name '" << id << "'";
        ArgErrLog(os.str());
    }

    AssertLog(vessdiff->second != nullptr);
    return vessdiff->second;
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::delVesSDiff(std::string const& id) const {
    VesSDiff* vsdiff = getVesSDiff(id);
    // Previous function has checked it's not a nullptr

    // delete VesSDiff object since it is owned by c++, not python
    delete (vsdiff);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<VesSDiff*> VesSurfsys::getAllVesSDiffs() const {
    VesSDiffPVec diffs = VesSDiffPVec();

    for (auto const& d: pVesSDiffs) {
        diffs.emplace_back(d.second);
    }
    return diffs;
}

////////////////////////////////////////////////////////////////////////////////

Exocytosis* VesSurfsys::getExocytosis(std::string const& id) const {
    auto exo = pExocytosis.find(id);
    if (exo == pExocytosis.end()) {
        std::ostringstream os;
        os << "Model does not contain exocytotic "
              "reaction with name '"
           << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(exo->second != nullptr);
    return exo->second;
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::delExocytosis(std::string const& id) const {
    Exocytosis* exo = getExocytosis(id);
    // Previous function has checked it's not a nullptr

    delete (exo);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Exocytosis*> VesSurfsys::getAllExocytosis() const {
    ExocytosisPVec exos = ExocytosisPVec();

    for (auto const& exo: pExocytosis) {
        exos.emplace_back(exo.second);
    }
    return exos;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec*> VesSurfsys::getAllSpecs() const {
    std::set<Spec*> specs_set;

    for (auto const& sreac: getAllVesSReacs()) {
        const auto& specs = sreac->getAllSpecs();
        specs_set.insert(specs.begin(), specs.end());
    }

    for (auto const& diff: getAllVesSDiffs()) {
        const auto& specs = diff->getAllSpecs();
        specs_set.insert(specs.begin(), specs.end());
    }
    return {specs_set.begin(), specs_set.end()};
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_checkVesSReacID(std::string const& id) const {
    checkID(id);
    if (pVesSReacs.find(id) != pVesSReacs.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_handleVesSReacIDChange(std::string const& o, std::string const& n) {
    auto sr_old = pVesSReacs.find(o);
    AssertLog(sr_old != pVesSReacs.end());

    if (o == n) {
        return;
    }
    _checkVesSReacID(n);

    VesSReac* sr = sr_old->second;
    AssertLog(sr != nullptr);
    pVesSReacs.erase(sr->getID());
    pVesSReacs.emplace(n, sr);
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_handleVesSReacAdd(VesSReac* sreac) {
    AssertLog(sreac->getVesSurfsys() == this);
    _checkVesSReacID(sreac->getID());
    pVesSReacs.emplace(sreac->getID(), sreac);
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_handleVesSReacDel(VesSReac* sreac) {
    AssertLog(sreac->getVesSurfsys() == this);
    pVesSReacs.erase(sreac->getID());
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_checkVesSDiffID(std::string const& id) const {
    checkID(id);
    if (pVesSDiffs.find(id) != pVesSDiffs.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_handleVesSDiffIDChange(std::string const& o, std::string const& n) {
    auto d_old = pVesSDiffs.find(o);
    AssertLog(d_old != pVesSDiffs.end());

    if (o == n) {
        return;
    }
    _checkVesSDiffID(n);

    VesSDiff* d = d_old->second;
    AssertLog(d != nullptr);
    pVesSDiffs.erase(d->getID());
    pVesSDiffs.emplace(n, d);
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_handleVesSDiffAdd(VesSDiff* diff) {
    AssertLog(diff->getVesSurfsys() == this);
    _checkVesSDiffID(diff->getID());
    pVesSDiffs.emplace(diff->getID(), diff);
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_handleVesSDiffDel(VesSDiff* diff) {
    AssertLog(diff->getVesSurfsys() == this);
    pVesSDiffs.erase(diff->getID());
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_checkExocytosisID(std::string const& id) const {
    checkID(id);
    if (pExocytosis.find(id) != pExocytosis.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_handleExocytosisIDChange(std::string const& o, std::string const& n) {
    auto e_old = pExocytosis.find(o);
    AssertLog(e_old != pExocytosis.end());

    if (o == n) {
        return;
    }
    _checkExocytosisID(n);

    Exocytosis* e = e_old->second;
    AssertLog(e != nullptr);
    pExocytosis.erase(e->getID());
    pExocytosis.emplace(n, e);
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_handleExocytosisAdd(Exocytosis* exocyt) {
    AssertLog(exocyt->getVesSurfsys() == this);
    _checkExocytosisID(exocyt->getID());
    pExocytosis.emplace(exocyt->getID(), exocyt);
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_handleExocytosisDel(Exocytosis* exocyt) {
    AssertLog(exocyt->getVesSurfsys() == this);
    pExocytosis.erase(exocyt->getID());
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_handleSpecDelete(Spec* spec) {
    std::vector<std::string> sreacs_del;

    for (auto const& sreac: pVesSReacs) {
        for (auto const& sr_spec: sreac.second->getAllSpecs()) {
            if (sr_spec == spec) {
                sreacs_del.emplace_back(sreac.second->getID());
                break;
            }
        }
    }

    for (const auto& sr_del: sreacs_del) {
        delVesSReac(sr_del);
    }

    std::vector<std::string> diffs_del;

    for (auto const& diff: pVesSDiffs) {
        for (auto const& d_spec: diff.second->getAllSpecs()) {
            if (d_spec == spec) {
                diffs_del.emplace_back(diff.second->getID());
                break;
            }
        }
    }

    for (const auto& d_del: diffs_del) {
        delVesSDiff(d_del);
    }
}

////////////////////////////////////////////////////////////////////////////////

void VesSurfsys::_handleLinkSpecDelete(LinkSpec* spec) {
    std::vector<std::string> sreacs_del;

    for (auto const& sreac: pVesSReacs) {
        for (auto const& sr_spec: sreac.second->getAllLinkSpecs()) {
            if (sr_spec == spec) {
                sreacs_del.emplace_back(sreac.second->getID());
                break;
            }
        }
    }

    for (const auto& sr_del: sreacs_del) {
        delVesSReac(sr_del);
    }
}

////////////////////////////////////////////////////////////////////////////////

VesSReac* VesSurfsys::_getVesSReac(uint lidx) const {
    AssertLog(lidx < pVesSReacs.size());
    auto sr_it = pVesSReacs.begin();
    std::advance(sr_it, lidx);
    return sr_it->second;
}

////////////////////////////////////////////////////////////////////////////////

VesSDiff* VesSurfsys::_getVesSDiff(uint lidx) const {
    AssertLog(lidx < pVesSDiffs.size());
    auto df_it = pVesSDiffs.begin();
    std::advance(df_it, lidx);
    return df_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Exocytosis* VesSurfsys::_getExocytosis(uint lidx) const {
    AssertLog(lidx < pExocytosis.size());
    auto exo_it = pExocytosis.begin();
    std::advance(exo_it, lidx);
    return exo_it->second;
}

}  // namespace steps::model
