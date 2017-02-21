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

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/model/model.hpp"
#include "steps/model/spec.hpp"
#include "steps/model/chan.hpp"
#include "steps/model/vdeptrans.hpp"
#include "steps/model/vdepsreac.hpp"
#include "steps/model/ghkcurr.hpp"
#include "steps/model/ohmiccurr.hpp"
#include "steps/model/surfsys.hpp"
#include "steps/model/volsys.hpp"

#include "steps/util/checkid.hpp"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

using steps::util::checkID;
//
////////////////////////////////////////////////////////////////////////////////

Model::Model(void)
: pSpecs()
, pChans()
, pVolsys()
, pSurfsys()
{
}

////////////////////////////////////////////////////////////////////////////////

Model::~Model(void)
{
    while (pSpecs.empty() == false)
    {
        SpecPMapCI spec = pSpecs.begin();
        delete(spec->second);
    }

    while (pChans.empty() == false)
    {
        ChanPMapCI chan = pChans.begin();
        delete(chan->second);
    }

    while (pVolsys.empty() == false)
    {
        VolsysPMapCI vsys = pVolsys.begin();
        delete(vsys->second);
    }

    while (pSurfsys.empty() == false)
    {
        SurfsysPMapCI ssys = pSurfsys.begin();
        delete(ssys->second);
    }
}

////////////////////////////////////////////////////////////////////////////////

Spec * Model::getSpec(string const & id) const
{
    SpecPMapCI spec = pSpecs.find(id);
    if (spec == pSpecs.end())
    {
        ostringstream os;
        os << "Model does not contain species with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
    assert(spec->second != 0);
    return spec->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delSpec(string const & id)
{
    Spec * spec = getSpec(id);
    delete(spec);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec *> Model::getAllSpecs(void) const
{
    SpecPVec specs = SpecPVec();
    SpecPMapCI spec_end = pSpecs.end();
    for (SpecPMapCI s = pSpecs.begin(); s != spec_end; ++s)
    {
        specs.push_back(s->second);
    }
    return specs;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Volsys *> Model::getAllVolsyss(void) const
{
    VolsysPVec volsyss = VolsysPVec();
    VolsysPMapCI vsys_end = pVolsys.end();
    for (VolsysPMapCI vs = pVolsys.begin(); vs != vsys_end; ++vs)
    {
        volsyss.push_back(vs->second);
    }
    return volsyss;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Surfsys *> Model::getAllSurfsyss(void) const
{
    SurfsysPVec surfsyss = SurfsysPVec();
    SurfsysPMapCI ssys_end = pSurfsys.end();
    for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ssys_end; ++ss)
    {
        surfsyss.push_back(ss->second);
    }
    return surfsyss;
}

////////////////////////////////////////////////////////////////////////////////

Chan * Model::getChan(string const & id) const
{
    ChanPMapCI chan = pChans.find(id);
    if (chan == pChans.end())
    {
        ostringstream os;
        os << "Model does not contain channel with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
    assert(chan->second != 0);
    return chan->second;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Chan *> Model::getAllChans(void) const
{
    ChanPVec chans = ChanPVec();
    ChanPMapCI chan_end = pChans.end();
    for (ChanPMapCI c = pChans.begin(); c != chan_end; ++c)
    {
        chans.push_back(c->second);
    }
    return chans;
}

////////////////////////////////////////////////////////////////////////////////

Volsys * Model::getVolsys(string const & id) const
{
    VolsysPMapCI volsys = pVolsys.find(id);
    if (volsys == pVolsys.end())
    {
        ostringstream os;
        os << "Model does not contain volume system with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
    assert(volsys->second != 0);
    return volsys->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delVolsys(string const & id)
{
    Volsys * volsys = getVolsys(id);
    delete(volsys);
}


////////////////////////////////////////////////////////////////////////////////

Surfsys * Model::getSurfsys(string const & id) const
{
    SurfsysPMapCI surfsys = pSurfsys.find(id);
    if (surfsys == pSurfsys.end())
    {
        ostringstream os;
        os << "Model does not contain surface system with name '" << id << "'";
        throw steps::ArgErr(os.str());
    }
    assert(surfsys->second != 0);
    return surfsys->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delSurfsys(string const & id)
{
    Surfsys * surfsys = getSurfsys(id);
    delete(surfsys);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkSpecID(string const & id) const
{
    checkID(id);
    if (pSpecs.find(id) != pSpecs.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSpecIDChange(string const & o, string const & n)
{
    SpecPMapCI s_old = pSpecs.find(o);
    assert(s_old != pSpecs.end());

    if (o == n) return;
    _checkSpecID(n);

    Spec * s = s_old->second;
    assert(s != 0);
    pSpecs.erase(s->getID());                        // or s_old->first
    pSpecs.insert(SpecPMap::value_type(n, s));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSpecAdd(Spec * spec)
{
    assert(spec->getModel() == this);
    _checkSpecID(spec->getID());
    pSpecs.insert(SpecPMap::value_type(spec->getID(), spec));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSpecDel(Spec * spec)
{
    VolsysPMapCI vsys_end = pVolsys.end();
    for (VolsysPMapCI vsys = pVolsys.begin(); vsys != vsys_end; ++vsys)
    {
        vsys->second->_handleSpecDelete(spec);
    }

    SurfsysPMapCI ssys_end = pSurfsys.end();
    for (SurfsysPMapCI ssys = pSurfsys.begin(); ssys != ssys_end; ++ssys)
    {
        ssys->second->_handleSpecDelete(spec);
    }

    pSpecs.erase(spec->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleChanDel(Chan * chan)
{

    SurfsysPMapCI ssys_end = pSurfsys.end();
    for (SurfsysPMapCI ssys = pSurfsys.begin(); ssys != ssys_end; ++ssys)
    {
        ssys->second->_handleChanDelete(chan);
    }

    pChans.erase(chan->getID());
}
////////////////////////////////////////////////////////////////////////////////

void Model::_checkChanID(string const & id) const
{
    checkID(id);
    if (pChans.find(id) != pChans.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleChanIDChange(string const & o, string const & n)
{
    ChanPMapCI c_old = pChans.find(o);
    assert(c_old != pChans.end());

    if (o == n) return;
    _checkChanID(n);

    Chan * c = c_old->second;
    assert(c != 0);
    pChans.erase(c->getID());                        // or c_old->first
    pChans.insert(ChanPMap::value_type(n, c));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleChanAdd(Chan * chan)
{
    assert(chan->getModel() == this);
    _checkChanID(chan->getID());
    pChans.insert(ChanPMap::value_type(chan->getID(), chan));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkVolsysID(string const & id) const
{
    checkID(id);
    if (pVolsys.find(id) != pVolsys.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVolsysIDChange(string const & o, string const & n)
{
    VolsysPMapCI v_old = pVolsys.find(o);
    assert(v_old != pVolsys.end());

    if (o == n) return;
    _checkVolsysID(n);

    Volsys * v = v_old->second;
    assert(v != 0);
    pVolsys.erase(v->getID());                        // or v_old->first
    pVolsys.insert(VolsysPMap::value_type(n, v));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVolsysAdd(Volsys * volsys)
{
    assert(volsys->getModel() == this);
    _checkVolsysID(volsys->getID());
    pVolsys.insert(VolsysPMap::value_type(volsys->getID(), volsys));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVolsysDel(Volsys * volsys)
{
    assert (volsys->getModel() == this);
    pVolsys.erase(volsys->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkSurfsysID(string const & id) const
{
    checkID(id);
    if (pSurfsys.find(id) != pSurfsys.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSurfsysIDChange(string const & o, string const & n)
{
    SurfsysPMapCI s_old = pSurfsys.find(o);
    assert(s_old != pSurfsys.end());

    if (o==n) return;
    _checkSurfsysID(n);

    Surfsys * s = s_old->second;
    assert(s != 0);
    pSurfsys.erase(s->getID());
    pSurfsys.insert(SurfsysPMap::value_type(n, s));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSurfsysAdd(Surfsys * surfsys)
{
    assert(surfsys->getModel() == this);
    _checkSurfsysID(surfsys->getID());
    pSurfsys.insert(SurfsysPMap::value_type(surfsys->getID(), surfsys));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSurfsysDel(Surfsys * surfsys)
{
    assert (surfsys->getModel() == this);
    pSurfsys.erase(surfsys->getID());
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countReacs(void) const
{
    uint nreacs = 0;

    VolsysPMapCI vs_end = pVolsys.end();
    for (VolsysPMapCI vs = pVolsys.begin(); vs != vs_end; ++vs)
    {
        nreacs += vs->second->_countReacs();
    }
    return nreacs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countSReacs(void) const
{
    uint nsreacs = 0;

    SurfsysPMapCI ss_end = pSurfsys.end();
    for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ss_end; ++ss)
    {
        nsreacs += ss->second->_countSReacs();
    }
    return nsreacs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countVDiffs(void) const
{
    uint ndiffs = 0;

    VolsysPMapCI vs_end = pVolsys.end();
    for (VolsysPMapCI vs = pVolsys.begin(); vs != vs_end; ++vs)
    {
        ndiffs += vs->second->_countDiffs();
    }
    return ndiffs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countSDiffs(void) const
{
    uint ndiffs = 0;

    SurfsysPMapCI ss_end = pSurfsys.end();
    for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ss_end; ++ss)
    {
        ndiffs += ss->second->_countDiffs();
    }
    return ndiffs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countVDepTrans(void) const
{
    uint nvdts = 0;

    SurfsysPMapCI ss_end = pSurfsys.end();
    for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ss_end; ++ss)
    {
        nvdts += ss->second->_countVDepTrans();
    }
    return nvdts;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countVDepSReacs(void) const
{
    uint nvdsrs = 0;

    SurfsysPMapCI ss_end = pSurfsys.end();
    for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ss_end; ++ss)
    {
        nvdsrs += ss->second->_countVDepSReacs();
    }
    return nvdsrs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countOhmicCurrs(void) const
{
    uint nocs = 0;

    SurfsysPMapCI ss_end = pSurfsys.end();
    for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ss_end; ++ss)
    {
        nocs += ss->second->_countOhmicCurrs();
    }
    return nocs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countGHKcurrs(void) const
{
    uint nghks = 0;

    SurfsysPMapCI ss_end = pSurfsys.end();
    for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ss_end; ++ss)
    {
        nghks += ss->second->_countGHKcurrs();
    }
    return nghks;
}

////////////////////////////////////////////////////////////////////////////////

Spec * Model::_getSpec(uint gidx) const
{
    assert (gidx < pSpecs.size());
    std::map<std::string, Spec *>::const_iterator sp_it = pSpecs.begin();
    for (uint i=0; i< gidx; ++i) ++sp_it;
    return sp_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Chan * Model::_getChan(uint gidx) const
{
    assert (gidx < pChans.size());
    std::map<std::string, Chan *>::const_iterator ch_it = pChans.begin();
    for (uint i=0; i< gidx; ++i) ++ch_it;
    return ch_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Reac * Model::_getReac(uint gidx) const
{
    // first find which volsys this reac (by global index) belongs to
    uint lidx = gidx;
    VolsysPMapCI vs_end = pVolsys.end();
    for (VolsysPMapCI vs = pVolsys.begin(); vs != vs_end; ++vs)
    {
        uint reacs_tot = vs->second->_countReacs();
        if (reacs_tot > lidx) return vs->second->_getReac(lidx);
        lidx -=reacs_tot;
    }

    // we shouldn't have gotten here
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////

SReac * Model::_getSReac(uint gidx) const
{
    // first find which surfsys this sreac (by global index) belongs to
    uint lidx = gidx;
    SurfsysPMapCI ss_end = pSurfsys.end();
    for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ss_end; ++ss)
    {
        uint sreacs_tot = ss->second->_countSReacs();
        if (sreacs_tot > lidx) return ss->second->_getSReac(lidx);
        lidx -=sreacs_tot;
    }

    // we shouldn't have gotten here
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////

Diff * Model::_getVDiff(uint gidx) const
{
    // first find which volsys this diff (by global index) belongs to
    uint lidx = gidx;
    VolsysPMapCI vs_end = pVolsys.end();
    for (VolsysPMapCI vs = pVolsys.begin(); vs != vs_end; ++vs)
    {
        uint diffs_tot = vs->second->_countDiffs();
        if (diffs_tot > lidx) return vs->second->_getDiff(lidx);
        lidx -=diffs_tot;
    }
    // we shouldn't have gotten to the end
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////

Diff * Model::_getSDiff(uint gidx) const
{
    // first find which surfsys this diff (by global index) belongs to
    uint lidx = gidx;
    SurfsysPMapCI ss_end = pSurfsys.end();
    for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ss_end; ++ss)
    {
        uint diffs_tot = ss->second->_countDiffs();
        if (diffs_tot > lidx) return ss->second->_getDiff(lidx);
        lidx -=diffs_tot;
    }
    // we shouldn't have gotten to the end
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////

VDepTrans * Model::_getVDepTrans(uint gidx) const
{
    // first find which surfsys this v-dep-trans (by global index) belongs to
    uint lidx = gidx;
    SurfsysPMapCI ss_end = pSurfsys.end();
    for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ss_end; ++ss)
    {
        uint vdts_tot = ss->second->_countVDepTrans();
        if (vdts_tot > lidx) return ss->second->_getVDepTrans(lidx);
        lidx -=vdts_tot;
    }

    // we shouldn't have gotten here
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac * Model::_getVDepSReac(uint gidx) const
{
    // first find which surfsys this v-dep-sreac (by global index) belongs to
    uint lidx = gidx;
    SurfsysPMapCI ss_end = pSurfsys.end();
    for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ss_end; ++ss)
    {
        uint vdsrs_tot = ss->second->_countVDepSReacs();
        if (vdsrs_tot > lidx) return ss->second->_getVDepSReac(lidx);
        lidx -=vdsrs_tot;
    }

    // we shouldn't have gotten here
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurr * Model::_getOhmicCurr(uint gidx) const
{
    // first find which surfsys this ohmic current (by global index) belongs to
    uint lidx = gidx;
    SurfsysPMapCI ss_end = pSurfsys.end();
    for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ss_end; ++ss)
    {
        uint ocs_tot = ss->second->_countOhmicCurrs();
        if (ocs_tot > lidx) return ss->second->_getOhmicCurr(lidx);
        lidx -=ocs_tot;
    }

    // we shouldn't have gotten here
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr * Model::_getGHKcurr(uint gidx) const
{
    // first find which surfsys this ghk current (by global index) belongs to
    uint lidx = gidx;
    SurfsysPMapCI ss_end = pSurfsys.end();
    for (SurfsysPMapCI ss = pSurfsys.begin(); ss != ss_end; ++ss)
    {
        uint ghks_tot = ss->second->_countGHKcurrs();
        if (ghks_tot > lidx) return ss->second->_getGHKcurr(lidx);
        lidx -=ghks_tot;
    }

    // we shouldn't have gotten here
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////

// END
