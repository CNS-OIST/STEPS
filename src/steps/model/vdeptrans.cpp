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

#include "vdeptrans.hpp"

#include <cassert>
#include <sstream>
#include <string>

#include <easylogging++.h>

#include "chan.hpp"
#include "chanstate.hpp"
#include "model.hpp"
#include "surfsys.hpp"

#include "util/error.hpp"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

////////////////////////////////////////////////////////////////////////////////

VDepTrans::VDepTrans(std::string const &id, Surfsys *surfsys, ChanState *src,
                     ChanState *dst, std::vector<double> rate, double vmin,
                     double vmax, double dv, uint tablesize)
    : pID(id), pModel(nullptr), pSurfsys(surfsys), pChan(nullptr), pSrc(src),
      pDst(dst), pRate(), pVMin(vmin), pVMax(vmax), pDV(dv),
      pTablesize(tablesize) {
  ArgErrLogIf(pSurfsys == nullptr,
              "No surfsys provided to VDepTrans initializer function");

  ArgErrLogIf(pSrc->getChan() != pDst->getChan(),
              "Source channel state and destination channel state do not "
              "belong to the same channel");
  ArgErrLogIf(rate.size() != pTablesize,
              "Table of transition rates is not of expected size");

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

VDepTrans::~VDepTrans() {
  if (pSurfsys == nullptr) {
    return;
  }
  _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void VDepTrans::_handleSelfDelete() {
  pSurfsys->_handleVDepTransDel(this);
  delete[] pRate;
  pSrc = nullptr;
  pDst = nullptr;
  pSurfsys = nullptr;
  pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void VDepTrans::setID(string const &id) {
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

void VDepTrans::setSrc(ChanState *src) {
  AssertLog(src != nullptr);

  ArgErrLogIf(src->getChan() != pDst->getChan(),
              "Source channel state and destination channel state do not "
              "belong to the same channel");

  pSrc = src;
}

////////////////////////////////////////////////////////////////////////////////

void VDepTrans::setDst(ChanState *dst) {
  AssertLog(dst != nullptr);

  ArgErrLogIf(dst->getChan() != pSrc->getChan(),
              "Source channel state and destination channel state do not "
              "belong to the same channel");

  pDst = dst;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> VDepTrans::getRate() const {
  std::vector<double> rate(pRate, pRate + pTablesize);
  return rate;
}

////////////////////////////////////////////////////////////////////////////////
