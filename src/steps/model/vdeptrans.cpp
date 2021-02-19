/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

// STL headers.
#include <cassert>
#include <sstream>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/model/chan.hpp"
#include "steps/model/chanstate.hpp"
#include "steps/model/model.hpp"
#include "steps/model/surfsys.hpp"
#include "steps/model/vdeptrans.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

////////////////////////////////////////////////////////////////////////////////

VDepTrans::VDepTrans(std::string const & id, Surfsys * surfsys,
                     ChanState * src, ChanState * dst,
                     std::vector<double> rate, double vmin, double vmax,
                     double dv, uint tablesize)
: pID(id)
, pModel(nullptr)
, pSurfsys(surfsys)
, pChan(nullptr)
, pSrc(src)
, pDst(dst)
, pRate()
, pVMin(vmin)
, pVMax(vmax)
, pDV(dv)
, pTablesize(tablesize)
{
    if (pSurfsys == nullptr)
    {
        ostringstream os;
        os << "No surfsys provided to VDepTrans initializer function";
        ArgErrLog(os.str());
    }
    if (pSrc->getChan() != pDst->getChan())
    {
        ostringstream os;
        os << "Source channel state and destination channel state do not ";
        os << "belong to the same channel";
        ArgErrLog(os.str());
    }

    if (rate.size() != pTablesize)
    {
        ostringstream os;
        os << "Table of transition rates is not of expected size";
        ArgErrLog(os.str());
    }
    pModel = pSurfsys->getModel();
    AssertLog(pModel != nullptr);

    pChan = pSrc->getChan();

    AssertLog(pDV > 0.0);

    // Copy the rate information to local array
    pRate = new double[pTablesize];
    std::memcpy(pRate, rate.data(), pTablesize * sizeof(double));

    pSurfsys->_handleVDepTransAdd(this);

}

////////////////////////////////////////////////////////////////////////////////

VDepTrans::~VDepTrans()
{
    if (pSurfsys == nullptr) { return;
}
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void VDepTrans::_handleSelfDelete()
{
    pSurfsys->_handleVDepTransDel(this);
    delete[] pRate;
    pSrc = nullptr;
    pDst = nullptr;
    pSurfsys = nullptr;
    pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void VDepTrans::setID(string const & id)
{
    AssertLog(pSurfsys != nullptr);
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pSurfsys->_handleVDepTransIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void VDepTrans::setSrc(ChanState * src)
{
    AssertLog(src != nullptr);

    if (src->getChan() != pDst->getChan())
    {
        ostringstream os;
        os << "Source channel state and destination channel state do not ";
        os << "belong to the same channel";
        ArgErrLog(os.str());
    }

    pSrc = src;
}

////////////////////////////////////////////////////////////////////////////////

void VDepTrans::setDst(ChanState * dst)
{
    AssertLog(dst != nullptr);

    if (dst->getChan() != pSrc->getChan())
    {
        ostringstream os;
        os << "Source channel state and destination channel state do not ";
        os << "belong to the same channel";
        ArgErrLog(os.str());
    }

    pDst = dst;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> VDepTrans::getRate() const
{
    std::vector<double> rate(pRate, pRate + pTablesize);
    return rate;
}

////////////////////////////////////////////////////////////////////////////////

// END
