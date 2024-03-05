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

#include "ghkcurr.hpp"

#include "chanstate.hpp"
#include "model.hpp"
#include "spec.hpp"
#include "surfsys.hpp"

#include "math/ghk.hpp"
#include "util/error.hpp"

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

GHKcurr::GHKcurr(std::string const& id,
                 Surfsys& surfsys,
                 ChanState& chanstate,
                 Spec& ion,
                 bool computeflux,
                 double virtual_oconc,
                 double vshift)
    : pID(id)
    , pModel(surfsys.getModel())
    , pSurfsys(surfsys)
    , pChanState(&chanstate)
    , pIon(&ion)
    , pRealFlux(computeflux)
    , pG(0.0)
    , pValence(0)
    , pV(0.0)
    , pTemp(0.0)
    , pInnerConc(0.0)
    , pOuterConc(0.0)
    , pP(0.0)
    , pInfoSupplied(false)
    , pVirtual_conc(virtual_oconc)
    , pVshift(vshift) {
    pValence = pIon->getValence();

    ArgErrLogIf(pValence == 0, "Ion provided to GHKcurr initializer function has valence zero");

    pSurfsys._handleGHKcurrAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr::~GHKcurr() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::setID(std::string const& id) {
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pSurfsys._handleGHKcurrIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::setChanState(ChanState& chanstate) {
    pChanState = &chanstate;
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::setIon(Spec& ion) {
    ArgErrLogIf(ion.getValence() == 0, "Ion provided to GHK::setIon function has valence zero");

    pValence = ion.getValence();  /// TODO Tristan fixed bug?
    pIon = &ion;

    if (pG != 0.0) {
        pP = math::permeability(pG, pV, pValence, pTemp, pInnerConc, pOuterConc);
    }
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::setP(double p) {
    ArgErrLogIf(p <= 0.0,
                "Permeability provided to GHKcurr::setP function can't "
                "be negative or zero");

    if (pG != 0.0) {
        CLOG(WARNING, "general_log") << "Permeability information previously defined for GHKcurr "
                                        "object will be overwritten.";
        pG = 0.0;
        pV = 0.0;
        pTemp = 0.0;
        pOuterConc = 0.0;
        pInnerConc = 0.0;
    }

    pP = p;

    pInfoSupplied = true;
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::setPInfo(double g, double V, double T, double oconc, double iconc) {
    if (pP != 0.0) {
        CLOG(WARNING, "general_log") << "Permeability previously defined for "
                                        "GHKcurr object will be overwritten.";
        pP = 0.0;
    }

    ArgErrLogIf(g <= 0.0,
                "Conductance provided to GHKcurr::setPInfo function "
                "can't be negative or zero");

    pG = g;

    ArgErrLogIf(V == 0.0, "Potential provided to GHKcurr::setPInfo function can't be zero.");

    pV = V;

    // Must be absurdly rare, but lets check for negative temperature.

    ArgErrLogIf(T < 0.0,
                "Temperature provided to GHKcurr::setPInfo function can't be "
                "negative. \nTemperature is required in Kelvin.");

    pTemp = T;

    ArgErrLogIf(oconc < 0.0,
                "Outer concentration provided to GHKcurr::setPInfo "
                "function can't be negative");

    pOuterConc = oconc;

    ArgErrLogIf(iconc < 0.0,
                "Inner concentration provided to GHKcurr::setPInfo "
                "function can't be negative");

    pInnerConc = iconc;

    pP = math::permeability(pG, pV, pValence, pTemp, pInnerConc, pOuterConc);
    pInfoSupplied = true;
}

////////////////////////////////////////////////////////////////////////////////

double GHKcurr::_P() const {
    AssertLog(_infosupplied() == true);
    return pP;
}

////////////////////////////////////////////////////////////////////////////////

int GHKcurr::_valence() const {
    AssertLog(_infosupplied() == true);
    return pValence;
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::_handleSelfDelete() {
    pSurfsys._handleGHKcurrDel(*this);
    pG = 0.0;
    pValence = 0;
    pV = 0.0;
    pTemp = 0.0;
    pInnerConc = 0.0;
    pOuterConc = 0.0;
    pP = 0.0;
    pInfoSupplied = false;
}

}  // namespace steps::model
