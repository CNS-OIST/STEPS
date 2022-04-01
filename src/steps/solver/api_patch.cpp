/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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
#include <sstream>
#include <string>

// STEPS headers.
#include "api.hpp"
#include "statedef.hpp"
// util
#include "util/error.hpp"
// logging
#include <easylogging++.h>
////////////////////////////////////////////////////////////////////////////////

USING(std, string);
using namespace steps::solver;

////////////////////////////////////////////////////////////////////////////////

double API::getPatchArea(string const &p) const {
  // the following may raise an exception if string is unused
  uint pidx = pStatedef->getPatchIdx(p);

  return _getPatchArea(pidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchArea(string const &p, double area) {
  ArgErrLogIf(area <= 0.0, "Area cannot be negative or zero.");

  // the following may raise an exception if string is unused
  uint pidx = pStatedef->getPatchIdx(p);

  _setPatchArea(pidx, area);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchCount(string const &p, string const &s) const {
  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sidx = pStatedef->getSpecIdx(s);

  return _getPatchCount(pidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchCount(string const &p, string const &s, double n) {
  ArgErrLogIf(n < 0.0, "Number of molecules cannot be negative.");

  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sidx = pStatedef->getSpecIdx(s);

  _setPatchCount(pidx, sidx, n);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchAmount(string const &p, string const &s) const {
  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sidx = pStatedef->getSpecIdx(s);

  return _getPatchAmount(pidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchAmount(string const &p, string const &s, double a) {
  ArgErrLogIf(a < 0.0, "Amount of mols cannot be negative.");

  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sidx = pStatedef->getSpecIdx(s);

  _setPatchAmount(pidx, sidx, a);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getPatchClamped(string const &p, string const &s) const {
  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sidx = pStatedef->getSpecIdx(s);

  return _getPatchClamped(pidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchClamped(string const &p, string const &s, bool buf) {
  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sidx = pStatedef->getSpecIdx(s);

  _setPatchClamped(pidx, sidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchSReacK(string const &p, string const &sr) const {
  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sridx = pStatedef->getSReacIdx(sr);

  return _getPatchSReacK(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchSReacK(string const &p, string const &sr, double kf) {
  ArgErrLogIf(kf < 0.0, "Reaction constant cannot be negative.");

  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sridx = pStatedef->getSReacIdx(sr);

  _setPatchSReacK(pidx, sridx, kf);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getPatchSReacActive(string const &p, string const &sr) const {
  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sridx = pStatedef->getSReacIdx(sr);

  return _getPatchSReacActive(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchSReacActive(string const &p, string const &sr, bool a) {
  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sridx = pStatedef->getSReacIdx(sr);

  _setPatchSReacActive(pidx, sridx, a);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getPatchVDepSReacActive(string const &p, string const &vsr) const {
  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint vsridx = pStatedef->getVDepSReacIdx(vsr);

  return _getPatchVDepSReacActive(pidx, vsridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchVDepSReacActive(string const &p, string const &vsr, bool a) {
  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint vsridx = pStatedef->getVDepSReacIdx(vsr);

  _setPatchVDepSReacActive(pidx, vsridx, a);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchSReacH(string const &p, string const &sr) const {
  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sridx = pStatedef->getSReacIdx(sr);

  return _getPatchSReacH(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchSReacC(string const &p, string const &sr) const {
  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sridx = pStatedef->getSReacIdx(sr);

  return _getPatchSReacC(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchSReacA(string const &p, string const &sr) const {
  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sridx = pStatedef->getSReacIdx(sr);

  return _getPatchSReacA(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long API::getPatchSReacExtent(string const &p,
                                            string const &sr) const {
  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sridx = pStatedef->getSReacIdx(sr);

  return _getPatchSReacExtent(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::resetPatchSReacExtent(string const &p, string const &sr) {
  // the following may raise exceptions if strings are unused
  uint pidx = pStatedef->getPatchIdx(p);
  uint sridx = pStatedef->getSReacIdx(sr);

  _resetPatchSReacExtent(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::_setPatchArea(uint /*pidx*/, double /*area*/) { NotImplErrLog(""); }

////////////////////////////////////////////////////////////////////////////////

double API::_getPatchSReacH(uint /*pidx*/, uint /*ridx*/) const {
  NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getPatchSReacC(uint /*pidx*/, uint /*ridx*/) const {
  NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getPatchSReacA(uint /*pidx*/, uint /*ridx*/) const {
  NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long API::_getPatchSReacExtent(uint /*pidx*/,
                                             uint /*ridx*/) const {
  NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_resetPatchSReacExtent(uint /*pidx*/, uint /*ridx*/) {
  NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getPatchVDepSReacActive(uint /*pidx*/, uint /*vsridx*/) const {
  NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setPatchVDepSReacActive(uint /*pidx*/, uint /*vsridx*/, bool /*a*/) {
  NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

// END
