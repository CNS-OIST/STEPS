/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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
#include <sstream>
#include <string>
#include <map>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/model/model.hpp"
#include "steps/model/spec.hpp"
#include "steps/model/chan.hpp"
#include "steps/model/chanstate.hpp"
#include "steps/model/surfsys.hpp"
#include "steps/model/sreac.hpp"
#include "steps/model/vdeptrans.hpp"
#include "steps/model/ohmiccurr.hpp"
#include "steps/model/ghkcurr.hpp"
#include "steps/model/vdepsreac.hpp"
#include "steps/model/diff.hpp"
#include "steps/util/checkid.hpp"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

using steps::util::checkID;

////////////////////////////////////////////////////////////////////////////////

Surfsys::Surfsys(string const & id, Model * model)
: pID(id)
, pModel(model)
, pSReacs()
, pVDepTrans()
, pVDepSReacs()
, pOhmicCurrs()
, pGHKcurrs()
{
    if (pModel == 0)
    {
        ostringstream os;
        os << "No model provided to Surfsys initializer function";
        throw steps::ArgErr(os.str());
    }
    pModel->_handleSurfsysAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Surfsys::~Surfsys(void)
{
    if (pModel == 0) return;
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::setID(string const & id)
{
    assert(pModel != 0);
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

void Surfsys::_handleSelfDelete(void)
{
    std::vector<steps::model::SReac *> allsreacs = getAllSReacs();
    SReacPVecCI sreac_end = allsreacs.end();
    for(SReacPVecCI sreac = allsreacs.begin(); sreac != sreac_end; ++sreac)
    {
        delete(*sreac);
    }

    std::vector<steps::model::VDepTrans *> allvdeptrans = getAllVDepTrans();
    VDepTransPVecCI vdtrans_end = allvdeptrans.end();
    for(VDepTransPVecCI vdtrans = allvdeptrans.begin(); vdtrans != vdtrans_end; ++vdtrans)
    {
        delete(*vdtrans);
    }

    std::vector<steps::model::VDepSReac *> allvdepsreacs = getAllVDepSReacs();
    VDepSReacPVecCI vdsreac_end = allvdepsreacs.end();
    for(VDepSReacPVecCI vdsreac = allvdepsreacs.begin(); vdsreac != vdsreac_end; ++vdsreac)
    {
        delete(*vdsreac);
    }

    std::vector<steps::model::OhmicCurr *> allocurrs = getAllOhmicCurrs();
    OhmicCurrPVecCI oc_end = allocurrs.end();
    for(OhmicCurrPVecCI oc = allocurrs.begin(); oc != oc_end; ++oc)
    {
        delete(*oc);
    }

    std::vector<steps::model::GHKcurr *> allghks = getAllGHKcurrs();
    GHKcurrPVecCI ghk_end = allghks.end();
    for(GHKcurrPVecCI ghk = allghks.begin(); ghk != ghk_end; ++ghk)
    {
        delete(*ghk);
    }

    std::vector<steps::model::Diff *> alldiffs = getAllDiffs();
    DiffPVecCI diff_end = alldiffs.end();
    for (DiffPVecCI diff = alldiffs.begin(); diff != diff_end; ++diff)
    {
        delete(*diff);
    }

    pModel->_handleSurfsysDel(this);

    pSReacs.clear();
    pVDepTrans.clear();
    pVDepSReacs.clear();
    pOhmicCurrs.clear();
    pGHKcurrs.clear();

    pDiffs.clear();

    pModel = 0;
}

////////////////////////////////////////////////////////////////////////////////

SReac * Surfsys::getSReac(string const & id) const
{
    SReacPMapCI sreac = pSReacs.find(id);
    if (sreac == pSReacs.end())
    {
        ostringstream os;
        os << "Model does not contain surface "
        "reaction with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
    assert(sreac->second != 0);
    return sreac->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delSReac(string const & id)
{
    SReac * sreac = getSReac(id);
    delete(sreac);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<SReac *> Surfsys::getAllSReacs(void) const
{
    SReacPVec sreacs = SReacPVec();
    SReacPMapCI sr_end = pSReacs.end();
    for (SReacPMapCI sr = pSReacs.begin(); sr != sr_end; ++sr)
    {
        sreacs.push_back(sr->second);
    }
    return sreacs;
}

////////////////////////////////////////////////////////////////////////////////

Diff * Surfsys::getDiff(string const & id) const
{
    DiffPMapCI diff = pDiffs.find(id);
    if (diff == pDiffs.end())
    {
        ostringstream os;
        os << "Model does not contain diffusion with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
    assert(diff->second != 0);
    return diff->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delDiff(string const & id)
{
    Diff * diff = getDiff(id);
    // delete diff object since it is owned by c++, not python
    delete(diff);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Diff *> Surfsys::getAllDiffs(void) const
{
    DiffPVec diffs = DiffPVec();
    DiffPMapCI d_end = pDiffs.end();
    for (DiffPMapCI d = pDiffs.begin(); d != d_end; ++d)
    {
        diffs.push_back(d->second);
    }
    return diffs;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec *> Surfsys::getAllSpecs(void) const
{
    SpecPVec specs = SpecPVec();
    bool first_occ = true;
    bool first_occ_dst = true;
    bool first_occ_src = true;

    SReacPVec sreacs = getAllSReacs();
    SReacPVecCI sreac_end = sreacs.end();
    for (SReacPVecCI sreac = sreacs.begin(); sreac != sreac_end; ++sreac)
    {
        SpecPVec sr_specs = (*sreac)->getAllSpecs();
        SpecPVecCI sr_spec_end = sr_specs.end();
        for (SpecPVecCI sr_spec = sr_specs.begin();
             sr_spec != sr_spec_end; ++sr_spec)
        {
            first_occ = true;
            SpecPVecCI allspecs_end = specs.end();
            for (SpecPVecCI allspecs = specs.begin();
                allspecs != allspecs_end; ++allspecs)
            {
                if ((*sr_spec) == (*allspecs))
                {
                    first_occ = false;
                    break;
                }
            }
            if (first_occ == true) specs.push_back(*sr_spec);
        }
    }

    VDepSReacPVec vdepsreacs = getAllVDepSReacs();
    VDepSReacPVecCI vdepsreac_end = vdepsreacs.end();
    for (VDepSReacPVecCI vdepsreac = vdepsreacs.begin(); vdepsreac != vdepsreac_end; ++vdepsreac)
    {
        SpecPVec sr_specs = (*vdepsreac)->getAllSpecs();
        SpecPVecCI sr_spec_end = sr_specs.end();
        for (SpecPVecCI sr_spec = sr_specs.begin();
             sr_spec != sr_spec_end; ++sr_spec)
        {
            first_occ = true;
            SpecPVecCI allspecs_end = specs.end();
            for (SpecPVecCI allspecs = specs.begin();
                allspecs != allspecs_end; ++allspecs)
            {
                if ((*sr_spec) == (*allspecs))
                {
                    first_occ = false;
                    break;
                }
            }
            if (first_occ == true) specs.push_back(*sr_spec);
        }
    }

    VDepTransPVec vdeptranss = getAllVDepTrans();
    VDepTransPVecCI vdeptrans_end = vdeptranss.end();
    for (VDepTransPVecCI vdeptrans = vdeptranss.begin(); vdeptrans != vdeptrans_end; ++vdeptrans)
    {
        SpecP dst = (*vdeptrans)->getDst();
        SpecP src = (*vdeptrans)->getSrc();
        first_occ_dst = true;
        first_occ_src = true;
        SpecPVecCI allspecs_end = specs.end();
        for (SpecPVecCI allspecs = specs.begin();
                allspecs != allspecs_end; ++allspecs)
        {
            if (dst == (*allspecs))
            {
                first_occ_dst = false;
                if (first_occ_src == false) break;
            }
            else if (src == (*allspecs))
            {
                first_occ_src = false;
                if (first_occ_dst == false) break;
            }
        }
        if (first_occ_dst == true) specs.push_back(dst);
        if (first_occ_src == true) specs.push_back(src);
    }

    GHKcurrPVec ghks = getAllGHKcurrs();
    GHKcurrPVecCI ghk_end = ghks.end();
    for(GHKcurrPVecCI ghk = ghks.begin(); ghk != ghk_end; ++ghk)
    {
        SpecP ghk_spec = (*ghk)->getIon();

        first_occ = true;
        SpecPVecCI allspecs_end = specs.end();
        for (SpecPVecCI allspecs = specs.begin();
            allspecs != allspecs_end; ++allspecs)
        {
            if (ghk_spec == (*allspecs))
            {
                first_occ = false;
                break;
            }
        }
        if (first_occ == true) specs.push_back(ghk_spec);
    }

    DiffPVec diffs = getAllDiffs();
    DiffPVecCI diff_end = diffs.end();
    for(DiffPVecCI diff = diffs.begin();diff != diff_end; ++diff)
    {
        SpecPVec d_specs = (*diff)->getAllSpecs();
        SpecPVecCI d_spec_end = d_specs.end();
        for (SpecPVecCI d_spec = d_specs.begin();
            d_spec != d_spec_end; ++d_spec)
        {
            first_occ = true;
            SpecPVecCI allspecs_end = specs.end();
            for (SpecPVecCI allspecs = specs.begin();
                allspecs != allspecs_end; ++allspecs)
            {
                if ((*d_spec) == (*allspecs))
                {
                    first_occ = false;
                    break;
                }
            }
            if (first_occ == true) specs.push_back(*d_spec);
        }
    }

    return specs;

    //std::cout << "\nSurfsys::getAllSpecs() called. Need to add stuff for channel states??";
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkSReacID(string const & id) const
{
    checkID(id);
    if (pSReacs.find(id) != pSReacs.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacIDChange(string const & o, string const & n)
{
    SReacPMapCI sr_old = pSReacs.find(o);
    assert(sr_old != pSReacs.end());

    if(o==n) return;
    _checkSReacID(n);

    SReac * sr = sr_old->second;
    assert(sr != 0);
    pSReacs.erase(sr->getID());
    pSReacs.insert(SReacPMap::value_type(n,sr));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacAdd(SReac * sreac)
{
    assert(sreac->getSurfsys() == this);
    _checkSReacID(sreac->getID());
    pSReacs.insert(SReacPMap::value_type(sreac->getID(), sreac));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacDel(SReac * sreac)
{
    assert (sreac->getSurfsys() == this);
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
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleDiffIDChange(string const & o, string const & n)
{
    DiffPMapCI d_old = pDiffs.find(o);
    assert(d_old != pDiffs.end());

    if (o==n) return;
    _checkDiffID(n);

    Diff * d = d_old->second;
    assert(d != 0);
    pDiffs.erase(d->getID());
    pDiffs.insert(DiffPMap::value_type(n,d));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleDiffAdd(Diff * diff)
{
    assert(diff->getSurfsys() == this);
    _checkDiffID(diff->getID());
    pDiffs.insert(DiffPMap::value_type(diff->getID(), diff));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleDiffDel(Diff * diff)
{
    assert (diff->getSurfsys() == this);
    pDiffs.erase(diff->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSpecDelete(Spec * spec)
{
    SReacPMapCI sreac_end = pSReacs.end();
    std::vector<std::string> sreacs_del = std::vector<std::string>();
    for (SReacPMapCI sreac = pSReacs.begin(); sreac != sreac_end; ++sreac)
    {
        SpecPVec specs = (sreac->second->getAllSpecs());
        SpecPVecCI sr_spec_end = specs.end();
        for (SpecPVecCI sr_spec = specs.begin();
             sr_spec != sr_spec_end; ++sr_spec)
        {
            if ((*sr_spec)== spec)
            {
                sreacs_del.push_back(sreac->second->getID());
                break;
            }
        }
    }
    std::vector<std::string>::const_iterator sr_del_end = sreacs_del.end();
    for (std::vector<std::string>::const_iterator sr_del = sreacs_del.begin();
         sr_del != sr_del_end; ++sr_del)
    {
        delSReac(*sr_del);
    }

    GHKcurrPMapCI ghk_end = pGHKcurrs.end();
    std::vector<std::string> ghks_del = std::vector<std::string>();
    for (GHKcurrPMapCI ghk = pGHKcurrs.begin(); ghk != ghk_end; ++ghk)
    {
        SpecP ion = ghk->second->getIon();
        if (ion == spec)
        {
            ghks_del.push_back(ghk->second->getID());
        }
        // spec may be a channel state
        SpecP cstate = ghk->second->getChanState();
        if (cstate == spec)
        {
            ghks_del.push_back(ghk->second->getID());
        }
    }
    std::vector<std::string>::const_iterator ghk_del_end = ghks_del.end();
    for (std::vector<std::string>::const_iterator ghkcurr_del = ghks_del.begin();
         ghkcurr_del != ghk_del_end; ++ghkcurr_del)
    {
        delGHKcurr(*ghkcurr_del);
    }

    // spec may also be a derived ChanState object -> need to delete any
    // vdeptrans and ohmic currents that include this channel state
    OhmicCurrPMapCI oc_end = pOhmicCurrs.end();
    std::vector<std::string> oc_del = std::vector<std::string>();
    for (OhmicCurrPMapCI oc = pOhmicCurrs.begin(); oc != oc_end; ++oc)
    {
        SpecP cstate = oc->second->getChanState();
        if (cstate == spec)
        {
            oc_del.push_back(oc->second->getID());
        }
    }
    std::vector<std::string>::const_iterator oc_del_end = oc_del.end();
    for (std::vector<std::string>::const_iterator occurr_del = oc_del.begin();
         occurr_del != oc_del_end; ++occurr_del)
    {
        delOhmicCurr(*occurr_del);
    }

    VDepTransPMapCI vdt_end = pVDepTrans.end();
    std::vector<std::string> vdt_del = std::vector<std::string>();
    for(VDepTransPMapCI vdt = pVDepTrans.begin(); vdt != vdt_end; ++vdt)
    {
        SpecP dst = vdt->second->getDst();
        SpecP src = vdt->second->getSrc();
        if(dst == spec or src == spec)
        {
            vdt_del.push_back(vdt->second->getID());
        }
    }
    std::vector<std::string>::const_iterator vdt_del_end = vdt_del.end();
    for(std::vector<std::string>::const_iterator vdept_del = vdt_del.begin();
        vdept_del != vdt_del_end; ++vdept_del)
    {
        delVDepTrans(*vdept_del);
    }

    VDepSReacPMapCI vdepsreac_end = pVDepSReacs.end();
    std::vector<std::string> vdepsreacs_del = std::vector<std::string>();
    for (VDepSReacPMapCI vdepsreac = pVDepSReacs.begin(); vdepsreac != vdepsreac_end; ++vdepsreac)
    {
        SpecPVec specs = (vdepsreac->second->getAllSpecs());
        SpecPVecCI sr_spec_end = specs.end();
        for (SpecPVecCI sr_spec = specs.begin();
             sr_spec != sr_spec_end; ++sr_spec)
        {
            if ((*sr_spec) == spec)
            {
                vdepsreacs_del.push_back(vdepsreac->second->getID());
                break;
            }
        }
    }
    std::vector<std::string>::const_iterator vdsr_del_end = vdepsreacs_del.end();
    for (std::vector<std::string>::const_iterator vdsr_del = vdepsreacs_del.begin();
         vdsr_del != vdsr_del_end; ++vdsr_del)
    {
        delVDepSReac(*vdsr_del);
    }

    std::vector<std::string> diffs_del = std::vector<std::string>();
    DiffPMapCI diff_end = pDiffs.end();
    for (DiffPMapCI diff = pDiffs.begin(); diff != diff_end; ++diff)
    {
        SpecPVec specs = (diff->second->getAllSpecs());
        SpecPVecCI d_spec_end = specs.end();
        for (SpecPVecCI d_spec = specs.begin();
             d_spec != d_spec_end; ++d_spec)
        {
            if ((*d_spec) == spec)
            {
                diffs_del.push_back(diff->second->getID());
                break;
            }
        }
    }
    std::vector<std::string>::const_iterator d_del_end = diffs_del.end();
    for (std::vector<std::string>::const_iterator d_del = diffs_del.begin();
         d_del != d_del_end; ++d_del)
    {
        delDiff(*d_del);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleChanDelete(Chan * chan)
{
    VDepTransPMapCI vdtrans_end = pVDepTrans.end();
    std::vector<std::string> vdtrans_del = std::vector<std::string>();
    for (VDepTransPMapCI vdtrans = pVDepTrans.begin(); vdtrans != vdtrans_end; ++vdtrans)
    {
        ChanP chans = (vdtrans->second->getChan());
        if (chans == chan)
        {
            vdtrans_del.push_back(vdtrans->second->getID());
        }
    }
    std::vector<std::string>::const_iterator vd_del_end = vdtrans_del.end();
    for (std::vector<std::string>::const_iterator vd_del = vdtrans_del.begin();
         vd_del != vd_del_end; ++vd_del)
    {
        delVDepTrans(*vd_del);
    }

    OhmicCurrPMapCI ohmcurr_end = pOhmicCurrs.end();
    std::vector<std::string> ohmcurr_del = std::vector<std::string>();
    for (OhmicCurrPMapCI ohmcurr = pOhmicCurrs.begin(); ohmcurr != ohmcurr_end; ++ohmcurr)
    {
        ChanP chans = (ohmcurr->second->getChanState()->getChan());
        if (chans == chan)
        {
            ohmcurr_del.push_back(ohmcurr->second->getID());
        }
    }
    std::vector<std::string>::const_iterator oc_del_end = ohmcurr_del.end();
    for (std::vector<std::string>::const_iterator oc_del = ohmcurr_del.begin();
         oc_del != oc_del_end; ++oc_del)
    {
        delOhmicCurr(*oc_del);
    }

    GHKcurrPMapCI ghkcurr_end = pGHKcurrs.end();
    std::vector<std::string> ghkcurr_del = std::vector<std::string>();
    for (GHKcurrPMapCI ghkcurr = pGHKcurrs.begin(); ghkcurr != ghkcurr_end; ++ghkcurr)
    {
        ChanP chans = (ghkcurr->second->getChanState()->getChan());
        if (chans == chan)
        {
            ghkcurr_del.push_back(ghkcurr->second->getID());
        }
    }
    std::vector<std::string>::const_iterator ghk_del_end = ghkcurr_del.end();
    for (std::vector<std::string>::const_iterator ghk_del = ghkcurr_del.begin();
         ghk_del != ghk_del_end; ++ghk_del)
    {
        delGHKcurr(*ghk_del);
    }
}

////////////////////////////////////////////////////////////////////////////////

VDepTrans * Surfsys::getVDepTrans(std::string const & id) const
{
    VDepTransPMapCI vdeptrans = pVDepTrans.find(id);
    if (vdeptrans == pVDepTrans.end())
    {
        ostringstream os;
        os << "Model does not contain voltage-dependent "
        "transition with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
    assert(vdeptrans->second != 0);
    return vdeptrans->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delVDepTrans(std::string const & id)
{
    VDepTrans * vdeptrans = getVDepTrans(id);
    delete(vdeptrans);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<VDepTrans *> Surfsys::getAllVDepTrans(void) const
{
    VDepTransPVec vdeptrans = VDepTransPVec();
    VDepTransPMapCI vd_end = pVDepTrans.end();
    for (VDepTransPMapCI vd = pVDepTrans.begin(); vd != vd_end; ++vd)
    {
        vdeptrans.push_back(vd->second);
    }
    return vdeptrans;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepTransIDChange(string const & o, string const & n)
{
    VDepTransPMapCI vd_old = pVDepTrans.find(o);
    assert(vd_old != pVDepTrans.end());

    if (o==n) return;
    _checkVDepTransID(n);

    VDepTrans * vd = vd_old->second;
    assert(vd != 0);
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
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepTransAdd(VDepTrans * vdeptrans)
{
    assert(vdeptrans->getSurfsys() == this);
    _checkVDepTransID(vdeptrans->getID());
    pVDepTrans.insert(VDepTransPMap::value_type(vdeptrans->getID(), vdeptrans));

}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepTransDel(VDepTrans * vdeptrans)
{
    assert(vdeptrans->getSurfsys() == this);
    pVDepTrans.erase(vdeptrans->getID());
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac * Surfsys::getVDepSReac(std::string const & id) const
{
    VDepSReacPMapCI vdepsreac = pVDepSReacs.find(id);
    if (vdepsreac == pVDepSReacs.end())
    {
        ostringstream os;
        os << "Model does not contain voltage-dependent "
        "surface reaction with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
    assert(vdepsreac->second != 0);
    return vdepsreac->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delVDepSReac(std::string const & id)
{
    VDepSReac * vdepsreac = getVDepSReac(id);
    delete(vdepsreac);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<VDepSReac *> Surfsys::getAllVDepSReacs(void) const
{
    VDepSReacPVec vdepsreac = VDepSReacPVec();
    VDepSReacPMapCI vd_end = pVDepSReacs.end();
    for (VDepSReacPMapCI vd = pVDepSReacs.begin(); vd != vd_end; ++vd)
    {
        vdepsreac.push_back(vd->second);
    }
    return vdepsreac;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepSReacIDChange(string const & o, string const & n)
{
    VDepSReacPMapCI vd_old = pVDepSReacs.find(o);
    assert(vd_old != pVDepSReacs.end());

    if (o==n) return;
    _checkVDepSReacID(n);

    VDepSReac * vd = vd_old->second;
    assert(vd != 0);
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
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepSReacAdd(VDepSReac * vdepsreac)
{
    assert(vdepsreac->getSurfsys() == this);
    _checkVDepSReacID(vdepsreac->getID());
    pVDepSReacs.insert(VDepSReacPMap::value_type(vdepsreac->getID(), vdepsreac));

}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepSReacDel(VDepSReac * vdepsreac)
{
    assert(vdepsreac->getSurfsys() == this);
    pVDepSReacs.erase(vdepsreac->getID());
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurr * Surfsys::getOhmicCurr(std::string const & id) const
{
    OhmicCurrPMapCI ohmiccurr = pOhmicCurrs.find(id);
    if (ohmiccurr == pOhmicCurrs.end())
    {
        ostringstream os;
        os << "Model does not contain ohmic current with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
    assert(ohmiccurr->second != 0);
    return ohmiccurr->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delOhmicCurr(std::string const & id)
{
    OhmicCurr * ohmiccurr = getOhmicCurr(id);
    delete(ohmiccurr);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<OhmicCurr *> Surfsys::getAllOhmicCurrs(void) const
{
    OhmicCurrPVec ohmiccurr = OhmicCurrPVec();
    OhmicCurrPMapCI oc_end = pOhmicCurrs.end();
    for (OhmicCurrPMapCI oc = pOhmicCurrs.begin(); oc != oc_end; ++oc)
    {
        ohmiccurr.push_back(oc->second);
    }
    return ohmiccurr;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleOhmicCurrIDChange(string const & o, string const & n)
{
    OhmicCurrPMapCI oc_old = pOhmicCurrs.find(o);
    assert(oc_old != pOhmicCurrs.end());

    if (o==n) return;
    _checkOhmicCurrID(n);

    OhmicCurr * oc = oc_old->second;
    assert(oc != 0);
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
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleOhmicCurrAdd(OhmicCurr * ohmiccurr)
{
    assert(ohmiccurr->getSurfsys() == this);
    _checkOhmicCurrID(ohmiccurr->getID());
    pOhmicCurrs.insert(OhmicCurrPMap::value_type(ohmiccurr->getID(), ohmiccurr));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleOhmicCurrDel(OhmicCurr * ohmiccurr)
{
    assert(ohmiccurr->getSurfsys() == this);
    pOhmicCurrs.erase(ohmiccurr->getID());
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr * Surfsys::getGHKcurr(std::string const & id) const
{
    GHKcurrPMapCI ghkcurr = pGHKcurrs.find(id);
    if (ghkcurr == pGHKcurrs.end())
    {
        ostringstream os;
        os << "Model does not contain ghk current with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
    assert(ghkcurr->second != 0);
    return ghkcurr->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delGHKcurr(std::string const & id)
{
    GHKcurr * ghkcurr = getGHKcurr(id);
    delete(ghkcurr);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<GHKcurr *> Surfsys::getAllGHKcurrs(void) const
{
    GHKcurrPVec ghkcurr = GHKcurrPVec();
    GHKcurrPMapCI ghk_end = pGHKcurrs.end();
    for (GHKcurrPMapCI ghk = pGHKcurrs.begin(); ghk != ghk_end; ++ghk)
    {
        ghkcurr.push_back(ghk->second);
    }
    return ghkcurr;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleGHKcurrIDChange(string const & o, string const & n)
{
    GHKcurrPMapCI ghk_old = pGHKcurrs.find(o);
    assert(ghk_old != pGHKcurrs.end());

    if (o==n) return;
    _checkGHKcurrID(n);

    GHKcurr * ghk = ghk_old->second;
    assert(ghk != 0);
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
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleGHKcurrAdd(GHKcurr * ghkcurr)
{
    assert(ghkcurr->getSurfsys() == this);
    _checkGHKcurrID(ghkcurr->getID());
    pGHKcurrs.insert(GHKcurrPMap::value_type(ghkcurr->getID(), ghkcurr));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleGHKcurrDel(GHKcurr * ghkcurr)
{
    assert(ghkcurr->getSurfsys() == this);
    pGHKcurrs.erase(ghkcurr->getID());
}

////////////////////////////////////////////////////////////////////////////////

SReac * Surfsys::_getSReac(uint lidx) const
{
    assert (lidx < pSReacs.size());
    std::map<std::string, SReac *>::const_iterator sr_it = pSReacs.begin();
    for (uint i=0; i< lidx; ++i) ++sr_it;
    return sr_it->second;
}

////////////////////////////////////////////////////////////////////////////////

VDepTrans * Surfsys::_getVDepTrans(uint lidx) const
{
    assert (lidx < pVDepTrans.size());
    std::map<std::string, VDepTrans *>::const_iterator vd_it = pVDepTrans.begin();
    for (uint i=0; i< lidx; ++i) ++vd_it;
    return vd_it->second;
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac * Surfsys::_getVDepSReac(uint lidx) const
{
    assert (lidx < pVDepSReacs.size());
    std::map<std::string, VDepSReac *>::const_iterator vd_it = pVDepSReacs.begin();
    for (uint i=0; i< lidx; ++i) ++vd_it;
    return vd_it->second;
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurr * Surfsys::_getOhmicCurr(uint lidx) const
{
    assert (lidx < pOhmicCurrs.size());
    std::map<std::string, OhmicCurr *>::const_iterator oc_it = pOhmicCurrs.begin();
    for (uint i=0; i< lidx; ++i) ++oc_it;
    return oc_it->second;
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr * Surfsys::_getGHKcurr(uint lidx) const
{
    assert (lidx < pGHKcurrs.size());
    std::map<std::string, GHKcurr *>::const_iterator ghk_it = pGHKcurrs.begin();
    for (uint i=0; i< lidx; ++i) ++ghk_it;
    return ghk_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Diff * Surfsys::_getDiff(uint lidx) const
{
    assert (lidx < pDiffs.size());
    std::map<std::string, Diff *>::const_iterator df_it = pDiffs.begin();
    for (uint i=0; i< lidx; ++i) ++df_it;
    return df_it->second;
}

////////////////////////////////////////////////////////////////////////////////

// END
