/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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

/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <map>
#include <sstream>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/model/diff.hpp"
#include "steps/model/model.hpp"
#include "steps/model/reac.hpp"
#include "steps/model/spec.hpp"
#include "steps/model/volsys.hpp"
#include "steps/util/checkid.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

using steps::util::checkID;

////////////////////////////////////////////////////////////////////////////////

Volsys::Volsys(string const & id, Model * model)
: pID(id)
, pModel(model)
{
    if (pModel == nullptr)
    {
        ostringstream os;
        os << "No model provided to Volsys initializer function";
        ArgErrLog(os.str());
    }
    pModel->_handleVolsysAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Volsys::~Volsys()
{
    if (pModel == nullptr) { return;
}
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleSelfDelete()
{
    for (auto const& reac: getAllReacs()) {
        delete reac;
    }

    for (auto const& diff: getAllDiffs()) {
        delete diff;
    }
    pModel->_handleVolsysDel(this);
    pReacs.clear();
    pDiffs.clear();
    pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::setID(string const & id)
{
    AssertLog(pModel != nullptr);
    if (id == pID) return;
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pModel->_handleVolsysIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

Reac * Volsys::getReac(string const & id) const
{
    auto reac = pReacs.find(id);
    if (reac == pReacs.end())
    {
        ostringstream os;
        os << "Model does not contain reaction with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(reac->second != nullptr);
    return reac->second;
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::delReac(string const & id)
{
    auto reac = getReac(id);
    // Delete reac object since it is owned by c++, not python
    delete reac;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Reac *> Volsys::getAllReacs() const
{
    ReacPVec reacs;
    reacs.reserve(pReacs.size());
    for (auto const& r: pReacs) {
      reacs.push_back(r.second);
    }
    return reacs;
}

////////////////////////////////////////////////////////////////////////////////

Diff * Volsys::getDiff(string const & id) const
{
    auto diff = pDiffs.find(id);
    if (diff == pDiffs.end())
    {
        ostringstream os;
        os << "Model does not contain diffusion with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(diff->second != nullptr);
    return diff->second;
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::delDiff(string const & id)
{
    auto diff = getDiff(id);
    // delete diff object since it is owned by c++, not python
    delete diff;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Diff *> Volsys::getAllDiffs() const
{
    DiffPVec diffs;
    diffs.reserve(pDiffs.size());
    for (auto const& d: pDiffs) {
      diffs.push_back(d.second);
    }
    return diffs;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec *> Volsys::getAllSpecs() const
{
    SpecPVec specs;

    for (auto const& reac: getAllReacs()) {
        for (auto const& r_spec: reac->getAllSpecs()) {
            bool first_occ = true;
            for (auto const& allspecs: specs) {
                if (r_spec == allspecs) {
                    first_occ = false;
                    break;
                }
            }
            if (first_occ) specs.push_back(r_spec);
        }
    }

    for (auto const& diff: getAllDiffs()) {
        for (auto const & d_spec: diff->getAllSpecs()) {
            bool first_occ = true;
            for (auto const& allspecs: specs) {
                if (d_spec == allspecs) {
                    first_occ = false;
                    break;
                }
            }
            if (first_occ) specs.push_back(d_spec);
        }
    }

    return specs;
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_checkReacID(string const & id) const
{
    checkID(id);
    if (pReacs.find(id) != pReacs.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleReacIDChange(string const & o, string const & n)
{
    auto r_old = pReacs.find(o);
    AssertLog(r_old != pReacs.end());

    if(o==n) return;
    _checkReacID(n);

    Reac * r = r_old->second;
    AssertLog(r != nullptr);
    pReacs.erase(r->getID());
    pReacs.insert(ReacPMap::value_type(n,r));
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleReacAdd(Reac * reac)
{
    AssertLog(reac->getVolsys() == this);
    _checkReacID(reac->getID());
    pReacs.insert(ReacPMap::value_type(reac->getID(), reac));
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleReacDel(Reac * reac)
{
    AssertLog(reac->getVolsys() == this);
    pReacs.erase(reac->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_checkDiffID(string const & id) const
{
    checkID(id);
    if (pDiffs.find(id) != pDiffs.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleDiffIDChange(string const & o, string const & n)
{
    auto d_old = pDiffs.find(o);
    AssertLog(d_old != pDiffs.end());

    if (o==n) return;
    _checkDiffID(n);

    Diff * d = d_old->second;
    AssertLog(d != nullptr);
    pDiffs.erase(d->getID());
    pDiffs.insert(DiffPMap::value_type(n,d));
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleDiffAdd(Diff * diff)
{
    AssertLog(diff->getVolsys() == this);
    _checkDiffID(diff->getID());
    pDiffs.insert(DiffPMap::value_type(diff->getID(), diff));
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleDiffDel(Diff * diff)
{
    AssertLog(diff->getVolsys() == this);
    pDiffs.erase(diff->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleSpecDelete(Spec * spec)
{
  {
    std::vector<std::string> reacs_del;
    for (auto const& reac: pReacs) {
      for (auto const& r_spec: reac.second->getAllSpecs()) {
        if (r_spec == spec) {
          reacs_del.push_back(reac.second->getID());
          break;
        }
      }
    }
    for (auto const &r_del: reacs_del) {
      delReac(r_del);
    }
  }
  {
    std::vector<std::string> diffs_del;
    for (auto const& diff: pDiffs) {
        for (auto const& d_spec: diff.second->getAllSpecs()) {
            if (d_spec == spec) {
                diffs_del.push_back(diff.second->getID());
                break;
            }
        }
    }
    for (auto const& d_del: diffs_del) {
        delDiff(d_del);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

Reac * Volsys::_getReac(uint lidx) const
{
    AssertLog(lidx < pReacs.size());
    auto rc_it = pReacs.begin();
    std::advance(rc_it, lidx);
    return rc_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Diff * Volsys::_getDiff(uint lidx) const
{
    AssertLog(lidx < pDiffs.size());
    auto df_it = pDiffs.begin();
    std::advance(df_it, lidx);
    return df_it->second;
}

////////////////////////////////////////////////////////////////////////////////

// END
