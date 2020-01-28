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


// STL headers.
#include <cassert>
#include <map>
#include <sstream>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/model/chan.hpp"
#include "steps/model/chanstate.hpp"
#include "steps/model/diff.hpp"
#include "steps/model/ghkcurr.hpp"
#include "steps/model/model.hpp"
#include "steps/model/ohmiccurr.hpp"
#include "steps/model/spec.hpp"
#include "steps/model/sreac.hpp"
#include "steps/model/surfsys.hpp"
#include "steps/model/vdepsreac.hpp"
#include "steps/model/vdeptrans.hpp"
#include "steps/util/checkid.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

using steps::util::checkID;

////////////////////////////////////////////////////////////////////////////////

Surfsys::Surfsys(string const & id, Model * model)
: pID(id)
, pModel(model)
{
    if (pModel == nullptr)
    {
        ostringstream os;
        os << "No model provided to Surfsys initializer function";
        ArgErrLog(os.str());
    }
    pModel->_handleSurfsysAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Surfsys::~Surfsys()
{
    if (pModel == nullptr) { return;
}
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::setID(string const & id)
{
    AssertLog(pModel != nullptr);
    if (id == pID) return;
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pModel->_handleSurfsysIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSelfDelete()
{
    for (auto const& sreac: getAllSReacs()) {
        delete sreac;
    }
    for (auto const& vdtrans: getAllVDepTrans()) {
        delete vdtrans;
    }
    for (auto const& vdsreac: getAllVDepSReacs()) {
        delete vdsreac;
    }
    for (auto const& oc: getAllOhmicCurrs()) {
        delete oc;
    }
    for (auto const& ghk: getAllGHKcurrs()) {
        delete ghk;
    }

    for (auto const& diff: getAllDiffs()) {
        delete diff;
    }

    pModel->_handleSurfsysDel(this);

    pSReacs.clear();
    pVDepTrans.clear();
    pVDepSReacs.clear();
    pOhmicCurrs.clear();
    pGHKcurrs.clear();

    pDiffs.clear();

    pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

SReac * Surfsys::getSReac(string const & id) const
{
    auto sreac = pSReacs.find(id);
    if (sreac == pSReacs.end())
    {
        ostringstream os;
        os << "Model does not contain surface "
        "reaction with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(sreac->second != nullptr);
    return sreac->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delSReac(string const & id)
{
    SReac * sreac = getSReac(id);
    delete sreac;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<SReac *> Surfsys::getAllSReacs() const
{
    SReacPVec sreacs;
    sreacs.reserve(pSReacs.size());
    for (auto const& sr: pSReacs) {
        sreacs.push_back(sr.second);
    }
    return sreacs;
}

////////////////////////////////////////////////////////////////////////////////

Diff * Surfsys::getDiff(string const & id) const
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

void Surfsys::delDiff(string const & id)
{
    Diff * diff = getDiff(id);
    // delete diff object since it is owned by c++, not python
    delete diff;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Diff *> Surfsys::getAllDiffs() const
{
    DiffPVec diffs;
    diffs.reserve(pDiffs.size());
    for (auto const& d: pDiffs) {
        diffs.push_back(d.second);
    }
    return diffs;
}

////////////////////////////////////////////////////////////////////////////////

SpecPVec Surfsys::getAllSpecs() const
{
    SpecPVec specs = SpecPVec();
    bool first_occ = true;
    for (auto const& sreac: getAllSReacs()) {
        for (auto const& sr_spec: sreac->getAllSpecs()) {
            first_occ = true;
            for (auto const& allspecs: specs) {
                if (sr_spec == allspecs) {
                    first_occ = false;
                    break;
                }
            }
            if (first_occ) specs.push_back(sr_spec);
        }
    }

    for (auto const& vdepsreac : getAllVDepSReacs()) {
        for (auto const& sr_spec: vdepsreac->getAllSpecs()) {
            for (auto const& allspecs: specs) {
                if (sr_spec == allspecs) {
                    first_occ = false;
                    break;
                }
            }
            if (first_occ) specs.push_back(sr_spec);
        }
    }

    for (auto const& vdeptrans: getAllVDepTrans()) {
        SpecP dst = vdeptrans->getDst();
        SpecP src = vdeptrans->getSrc();
        bool first_occ_dst = true;
        bool first_occ_src = true;
        for (auto const& allspecs: specs) {
            if (dst == allspecs) {
                first_occ_dst = false;
                if (!first_occ_src) break;
            } else if (src == allspecs) {
                first_occ_src = false;
                if (!first_occ_dst) break;
            }
        }
        if (first_occ_dst) specs.push_back(dst);
        if (first_occ_src) specs.push_back(src);
    }

    for (auto const& ghk: getAllGHKcurrs())
    {
        SpecP ghk_spec = ghk->getIon();

        first_occ = true;
        for (auto const& allspecs: specs) {
            if (ghk_spec == allspecs) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) specs.push_back(ghk_spec);
    }

    for (auto const& diff: getAllDiffs()) {
        for (auto const& d_spec: diff->getAllSpecs()) {
            first_occ = true;
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

    //CLOG(INFO, "general_log") << "\nSurfsys::getAllSpecs() called. Need to add stuff for channel states??";
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkSReacID(string const & id) const
{
    checkID(id);
    if (pSReacs.find(id) != pSReacs.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacIDChange(string const & o, string const & n)
{
    SReacPMapCI sr_old = pSReacs.find(o);
    AssertLog(sr_old != pSReacs.end());

    if(o==n) return;
    _checkSReacID(n);

    SReac * sr = sr_old->second;
    AssertLog(sr != nullptr);
    pSReacs.erase(sr->getID());
    pSReacs.emplace(n, sr);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacAdd(SReac * sreac)
{
    AssertLog(sreac->getSurfsys() == this);
    _checkSReacID(sreac->getID());
    pSReacs.insert(SReacPMap::value_type(sreac->getID(), sreac));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacDel(SReac * sreac)
{
    AssertLog(sreac->getSurfsys() == this);
    pSReacs.erase(sreac->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkDiffID(string const & id) const
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

void Surfsys::_handleDiffIDChange(string const & o, string const & n)
{
    DiffPMapCI d_old = pDiffs.find(o);
    AssertLog(d_old != pDiffs.end());

    if (o==n) return;
    _checkDiffID(n);

    Diff * d = d_old->second;
    AssertLog(d != nullptr);
    AssertLog(pDiffs.erase(d->getID()) == 1);
    pDiffs.insert(DiffPMap::value_type(n,d));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleDiffAdd(Diff * diff)
{
    AssertLog(diff->getSurfsys() == this);
    _checkDiffID(diff->getID());
    pDiffs.insert(DiffPMap::value_type(diff->getID(), diff));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleDiffDel(Diff * diff)
{
    AssertLog(diff->getSurfsys() == this);
    pDiffs.erase(diff->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSpecDelete(Spec * spec)
{
    {
        std::vector<std::string> sreacs_del;
        for (auto const& sreac: pSReacs) {
            SpecPVec specs = sreac.second->getAllSpecs();
            for (auto const& sr_spec: specs) {
                if (sr_spec == spec) {
                    sreacs_del.push_back(sreac.second->getID());
                    break;
                }
            }
        }
        for (auto const& sr_del: sreacs_del) {
            delSReac(sr_del);
        }
    }
    {
        std::vector<std::string> ghks_del;
        for (auto const& ghk: pGHKcurrs) {
            SpecP ion = ghk.second->getIon();
            if (ion == spec) {
                ghks_del.push_back(ghk.second->getID());
            }
            // spec may be a channel state
            SpecP cstate = ghk.second->getChanState();
            if (cstate == spec) {
                ghks_del.push_back(ghk.second->getID());
            }
        }
        for (auto const& ghk: ghks_del) {
            delGHKcurr(ghk);
        }
    }
    {
        // spec may also be a derived ChanState object -> need to delete any
        // vdeptrans and ohmic currents that include this channel state
        std::vector<std::string> oc_del;
        for (auto const& oc: pOhmicCurrs) {
            SpecP cstate = oc.second->getChanState();
            if (cstate == spec) {
                oc_del.push_back(cstate->getID());
            }
        }
        for (auto const& occurr_del: oc_del) {
            delOhmicCurr(occurr_del);
        }
    }
    {
        std::vector<std::string> vdt_del;
        for (auto const& vdt: pVDepTrans) {
            SpecP dst = vdt.second->getDst();
            SpecP src = vdt.second->getSrc();
            if (dst == spec or src == spec) {
                vdt_del.push_back(vdt.second->getID());
            }
        }
        for (auto const& vdept_del: vdt_del)  {
            delVDepTrans(vdept_del);
        }
    }
    {
        std::vector<std::string> vdepsreacs_del;
        for (auto const& vdepsreac: pVDepSReacs) {
            SpecPVec specs = vdepsreac.second->getAllSpecs();
            for (auto const& sr_spec: specs) {
                if (sr_spec == spec) {
                    vdepsreacs_del.push_back(vdepsreac.second->getID());
                    break;
                }
            }
        }
        for (auto const& vdsr_del: vdepsreacs_del) {
            delVDepSReac(vdsr_del);
        }
    }
    {
        std::vector<std::string> diffs_del;
        for (auto const& diff: pDiffs) {
            SpecPVec specs = diff.second->getAllSpecs();
            for (auto const& d_spec: specs) {
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

void Surfsys::_handleChanDelete(Chan * chan)
{
    {
        std::vector<std::string> vdtrans_del;
        for (auto const &vdtrans: pVDepTrans) {
            ChanP chans = vdtrans.second->getChan();
            if (chans == chan) {
                vdtrans_del.push_back(vdtrans.second->getID());
            }
        }
        for (auto const &vd_del: vdtrans_del) {
            delVDepTrans(vd_del);
        }
    }
    {
        std::vector<std::string> ohmcurr_del;
        for (auto const &ohmcurr: pOhmicCurrs) {
            ChanP chans = ohmcurr.second->getChanState()->getChan();
            if (chans == chan) {
                ohmcurr_del.push_back(ohmcurr.second->getID());
            }
        }
        for (auto const& oc_del: ohmcurr_del) {
            delOhmicCurr(oc_del);
        }
    }
    {
        std::vector<std::string> ghkcurr_del;
        for (auto const& ghkcurr: pGHKcurrs) {
            ChanP chans = ghkcurr.second->getChanState()->getChan();
            if (chans == chan) {
                ghkcurr_del.push_back(ghkcurr.second->getID());
            }
        }
        for (auto const& ghk_del: ghkcurr_del) {
            delGHKcurr(ghk_del);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

VDepTrans * Surfsys::getVDepTrans(std::string const & id) const
{
    auto vdeptrans = pVDepTrans.find(id);
    if (vdeptrans == pVDepTrans.end())
    {
        ostringstream os;
        os << "Model does not contain voltage-dependent "
        "transition with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(vdeptrans->second != nullptr);
    return vdeptrans->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delVDepTrans(std::string const & id)
{
    VDepTrans * vdeptrans = getVDepTrans(id);
    delete vdeptrans;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<VDepTrans *> Surfsys::getAllVDepTrans() const
{
    VDepTransPVec vdeptrans;
    vdeptrans.reserve(pVDepTrans.size());
    for (auto const& vd: pVDepTrans) {
        vdeptrans.push_back(vd.second);
    }
    return vdeptrans;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepTransIDChange(string const & o, string const & n)
{
    VDepTransPMapCI vd_old = pVDepTrans.find(o);
    AssertLog(vd_old != pVDepTrans.end());

    if (o==n) return;
    _checkVDepTransID(n);

    VDepTrans * vd = vd_old->second;
    AssertLog(vd != nullptr);
    pVDepTrans.erase(vd->getID());
    pVDepTrans.insert(VDepTransPMap::value_type(n,vd));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkVDepTransID(string const & id) const
{
    checkID(id);
    if (pVDepTrans.find(id) != pVDepTrans.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepTransAdd(VDepTrans * vdeptrans)
{
    AssertLog(vdeptrans->getSurfsys() == this);
    _checkVDepTransID(vdeptrans->getID());
    pVDepTrans.insert(VDepTransPMap::value_type(vdeptrans->getID(), vdeptrans));

}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepTransDel(VDepTrans * vdeptrans)
{
    AssertLog(vdeptrans->getSurfsys() == this);
    pVDepTrans.erase(vdeptrans->getID());
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac * Surfsys::getVDepSReac(std::string const & id) const
{
    auto vdepsreac = pVDepSReacs.find(id);
    if (vdepsreac == pVDepSReacs.end())
    {
        ostringstream os;
        os << "Model does not contain voltage-dependent "
        "surface reaction with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(vdepsreac->second != nullptr);
    return vdepsreac->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delVDepSReac(std::string const & id)
{
    VDepSReac * vdepsreac = getVDepSReac(id);
    delete vdepsreac;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<VDepSReac *> Surfsys::getAllVDepSReacs() const
{
    VDepSReacPVec vdepsreac;
    vdepsreac.reserve(pVDepSReacs.size());
    for (auto const& vd: pVDepSReacs) {
        vdepsreac.push_back(vd.second);
    }
    return vdepsreac;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepSReacIDChange(string const & o, string const & n)
{
    auto vd_old = pVDepSReacs.find(o);
    AssertLog(vd_old != pVDepSReacs.end());

    if (o==n) return;
    _checkVDepSReacID(n);

    VDepSReac * vd = vd_old->second;
    AssertLog(vd != nullptr);
    pVDepSReacs.erase(vd->getID());
    pVDepSReacs.insert(VDepSReacPMap::value_type(n,vd));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkVDepSReacID(string const & id) const
{
    checkID(id);
    if (pVDepSReacs.find(id) != pVDepSReacs.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepSReacAdd(VDepSReac * vdepsreac)
{
    AssertLog(vdepsreac->getSurfsys() == this);
    _checkVDepSReacID(vdepsreac->getID());
    pVDepSReacs.insert(VDepSReacPMap::value_type(vdepsreac->getID(), vdepsreac));

}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepSReacDel(VDepSReac * vdepsreac)
{
    AssertLog(vdepsreac->getSurfsys() == this);
    pVDepSReacs.erase(vdepsreac->getID());
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurr * Surfsys::getOhmicCurr(std::string const & id) const
{
    auto ohmiccurr = pOhmicCurrs.find(id);
    if (ohmiccurr == pOhmicCurrs.end())
    {
        ostringstream os;
        os << "Model does not contain ohmic current with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(ohmiccurr->second != nullptr);
    return ohmiccurr->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delOhmicCurr(std::string const & id)
{
    OhmicCurr * ohmiccurr = getOhmicCurr(id);
    delete ohmiccurr;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<OhmicCurr *> Surfsys::getAllOhmicCurrs() const
{
    OhmicCurrPVec ohmiccurr;
    ohmiccurr.reserve(pOhmicCurrs.size());
    for (auto const& oc: pOhmicCurrs) {
        ohmiccurr.push_back(oc.second);
    }
    return ohmiccurr;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleOhmicCurrIDChange(string const & o, string const & n)
{
    OhmicCurrPMapCI oc_old = pOhmicCurrs.find(o);
    AssertLog(oc_old != pOhmicCurrs.end());

    if (o==n) return;
    _checkOhmicCurrID(n);

    OhmicCurr * oc = oc_old->second;
    AssertLog(oc != nullptr);
    pOhmicCurrs.erase(oc->getID());
    pOhmicCurrs.insert(OhmicCurrPMap::value_type(n,oc));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkOhmicCurrID(string const & id) const
{
    checkID(id);
    if (pOhmicCurrs.find(id) != pOhmicCurrs.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleOhmicCurrAdd(OhmicCurr * ohmiccurr)
{
    AssertLog(ohmiccurr->getSurfsys() == this);
    _checkOhmicCurrID(ohmiccurr->getID());
    pOhmicCurrs.insert(OhmicCurrPMap::value_type(ohmiccurr->getID(), ohmiccurr));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleOhmicCurrDel(OhmicCurr * ohmiccurr)
{
    AssertLog(ohmiccurr->getSurfsys() == this);
    pOhmicCurrs.erase(ohmiccurr->getID());
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr * Surfsys::getGHKcurr(std::string const & id) const
{
    auto ghkcurr = pGHKcurrs.find(id);
    if (ghkcurr == pGHKcurrs.end())
    {
        ostringstream os;
        os << "Model does not contain ghk current with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(ghkcurr->second != nullptr);
    return ghkcurr->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delGHKcurr(std::string const & id)
{
    GHKcurr * ghkcurr = getGHKcurr(id);
    delete ghkcurr;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<GHKcurr *> Surfsys::getAllGHKcurrs() const
{
    GHKcurrPVec ghkcurr;
    ghkcurr.reserve(pGHKcurrs.size());
    for (auto const& ghk: pGHKcurrs) {
        ghkcurr.push_back(ghk.second);
    }
    return ghkcurr;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleGHKcurrIDChange(string const & o, string const & n)
{
    auto ghk_old = pGHKcurrs.find(o);
    AssertLog(ghk_old != pGHKcurrs.end());

    if (o==n) return;
    _checkGHKcurrID(n);

    GHKcurr * ghk = ghk_old->second;
    AssertLog(ghk != nullptr);
    pGHKcurrs.erase(ghk->getID());
    pGHKcurrs.insert(GHKcurrPMap::value_type(n,ghk));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkGHKcurrID(string const & id) const
{
    checkID(id);
    if (pGHKcurrs.find(id) != pGHKcurrs.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleGHKcurrAdd(GHKcurr * ghkcurr)
{
    AssertLog(ghkcurr->getSurfsys() == this);
    _checkGHKcurrID(ghkcurr->getID());
    pGHKcurrs.insert(GHKcurrPMap::value_type(ghkcurr->getID(), ghkcurr));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleGHKcurrDel(GHKcurr * ghkcurr)
{
    AssertLog(ghkcurr->getSurfsys() == this);
    pGHKcurrs.erase(ghkcurr->getID());
}

////////////////////////////////////////////////////////////////////////////////

SReac * Surfsys::_getSReac(uint lidx) const
{
    AssertLog(lidx < pSReacs.size());
    auto sr_it = pSReacs.begin();
    std::advance(sr_it, lidx);
    return sr_it->second;
}

////////////////////////////////////////////////////////////////////////////////

VDepTrans * Surfsys::_getVDepTrans(uint lidx) const
{
    AssertLog(lidx < pVDepTrans.size());
    auto vd_it = pVDepTrans.begin();
    std::advance(vd_it, lidx);
    return vd_it->second;
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac * Surfsys::_getVDepSReac(uint lidx) const
{
    AssertLog(lidx < pVDepSReacs.size());
    auto vd_it = pVDepSReacs.begin();
    std::advance(vd_it, lidx);
    return vd_it->second;
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurr * Surfsys::_getOhmicCurr(uint lidx) const
{
    AssertLog(lidx < pOhmicCurrs.size());
    auto oc_it = pOhmicCurrs.begin();
    std::advance(oc_it, lidx);
    return oc_it->second;
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr * Surfsys::_getGHKcurr(uint lidx) const
{
    AssertLog(lidx < pGHKcurrs.size());
    auto ghk_it = pGHKcurrs.begin();
    std::advance(ghk_it, lidx);
    return ghk_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Diff * Surfsys::_getDiff(uint lidx) const
{
    AssertLog(lidx < pDiffs.size());
    auto df_it = pDiffs.begin();
    std::advance(df_it, lidx);
    return df_it->second;
}

////////////////////////////////////////////////////////////////////////////////

// END
