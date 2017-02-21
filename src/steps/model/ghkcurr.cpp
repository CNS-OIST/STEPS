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

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <sstream>
#include <string>
#include <iostream>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/model/model.hpp"
#include "steps/model/surfsys.hpp"
#include "steps/model/ghkcurr.hpp"
#include "steps/model/spec.hpp"
#include "steps/model/chanstate.hpp"


////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

////////////////////////////////////////////////////////////////////////////////

GHKcurr::GHKcurr(string const & id, Surfsys * surfsys, ChanState * chanstate,
                 Spec * ion, bool computeflux, double virtual_oconc, double vshift)
: pID(id)
, pModel(0)
, pSurfsys(surfsys)
, pChanState(chanstate)
, pIon(ion)
, pG(0.0)
, pValence(0)
, pV(0.0)
, pTemp(0.0)
, pInnerConc(0.0)
, pOuterConc(0.0)
, pP(0.0)
, pRealFlux(computeflux)
, pVshift(vshift)
, pInfoSupplied(false)
, pVirtual_conc(virtual_oconc)
{
    if (pSurfsys == 0)
    {
        ostringstream os;
        os << "No surfsys provided to GHKcurr initializer function";
        throw steps::ArgErr(os.str());
    }
    if (pChanState == 0)
    {
        ostringstream os;
        os << "No channel state provided to GHKcurr initializer function";
        throw steps::ArgErr(os.str());
    }
    if (pIon == 0)
    {
        ostringstream os;
        os << "No ion provided to GHKcurr initializer function";
        throw steps::ArgErr(os.str());
    }

    pValence = pIon->getValence();
    if (pValence == 0)
    {
        ostringstream os;
        os << "Ion provided to GHKcurr initializer function has valence zero";
        throw steps::ArgErr(os.str());
    }

    /*
    pGInfo.insert(pair<string,double>("valence", static_cast<double>(pIon->getValence())));

    if (pG < 0.0)
    {
        ostringstream os;
        os << "Channel conductance provided to GHKcurr initializer";
        os << " function can't be negative";
        throw steps::ArgErr(os.str());
    }

    // Scan the conductance info and if there are no errors, fill pGInfo

    map<string, double>::const_iterator gi_end = ginfo.end();
    for (map<string, double>::const_iterator gi = ginfo.begin(); gi != gi_end; ++gi)
    {
        if(gi->first == "V" or gi->first == "temp"
            or gi->first == "iconc" or gi->first == "oconc")
        {
            pGInfo.insert(pair<string,double>(gi->first,gi->second));
        }
        else
        {
            ostringstream os;
            os << "Unknown key: '" << gi->first << "' in conductance";
            os << " measurement information in GHKcurr initialiser function.";
            os << "\nAccepted keys are: 'temp', 'V', 'iconc' and 'oconc'.";
            throw steps::ArgErr(os.str());
        }
    }
    */

    pModel = pSurfsys->getModel();
    assert (pModel != 0);

    pSurfsys->_handleGHKcurrAdd(this);

}

////////////////////////////////////////////////////////////////////////////////

GHKcurr::~GHKcurr(void)
{
    if (pSurfsys == 0) return;
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::setID(string const & id)
{
    assert(pSurfsys != 0);
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pSurfsys->_handleGHKcurrIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::setChanState(ChanState * chanstate)
{
    assert(chanstate != 0);
    pChanState = chanstate;
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::setIon(Spec * ion)
{
    assert (pSurfsys != 0);

    if (ion->getValence() == 0)
    {
        ostringstream os;
        os << "Ion provided to GHK::setIon function has valence zero";
        throw steps::ArgErr(os.str());
    }

    pValence = pIon->getValence();
    pIon = ion;
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::setP(double p)
{
    assert(pSurfsys != 0);
    if (p <= 0.0)
    {
        ostringstream os;
        os << "Permeability provided to GHKcurr::setP function can't be negative or zero";
        throw steps::ArgErr(os.str());
    }

    if (pG != 0.0)
    {
        ostringstream os;
        os << "\nWARNING: Permability information previously defined for GHKcurr object will be overwritten.\n";
        pG = 0.0;
        pV=0.0;
        pTemp=0.0;
        pOuterConc=0.0;
        pInnerConc=0.0;
    }

    pP = p;

    pInfoSupplied = true;
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::setPInfo(double g, double V, double T, double oconc, double iconc)
{
    assert(pSurfsys != 0);

    if (pP != 0.0)
    {
        ostringstream os;
        os << "\nWARNING: Permeability previously defined for GHKcurr object will be overwritten.\n";
        pP = 0.0;
    }

    if (g <= 0.0)
    {
        ostringstream os;
        os << "Conductance provided to GHKcurr::setPInfo function can't be negative or zero";
        throw steps::ArgErr(os.str());
    }
    pG = g;

    if (V == 0.0)
    {
        ostringstream os;
        os << "Potential provided to GHKcurr::setPInfo function can't be zero.";
        throw steps::ArgErr(os.str());
    }
    pV = V;

    // Must be absurdly rare, but lets check for negative temperature.
    if (T < 0.0)
    {
        ostringstream os;
        os << "Temperature provided to GHKcurr::setPInfo function can't be negative. ";
        os << "Temperature is required in Kelvin.";
        throw steps::ArgErr(os.str());
    }
    pTemp = T;

    if (oconc < 0.0)
    {
        ostringstream os;
        os << "Outer concentration provided to GHKcurr::setPInfo function can't be negative";
        throw steps::ArgErr(os.str());
    }

    pOuterConc = oconc;

    if (iconc < 0.0)
    {
        ostringstream os;
        os << "Inner concentration provided to GHKcurr::setPInfo function can't be negative";
        throw steps::ArgErr(os.str());
    }
    pInnerConc = iconc;

    pInfoSupplied = true;
}

////////////////////////////////////////////////////////////////////////////////

double GHKcurr::_G(void) const
{
    assert (_infosupplied() == true);
    return pG;
}

////////////////////////////////////////////////////////////////////////////////

double GHKcurr::_P(void) const
{
    assert (_infosupplied() == true);
    return pP;
}

////////////////////////////////////////////////////////////////////////////////

int GHKcurr::_valence(void) const
{
    assert (_infosupplied() == true);
    return pValence;
}

////////////////////////////////////////////////////////////////////////////////

double GHKcurr::_V(void) const
{
    assert (_infosupplied() == true);
    return pV;;
}

////////////////////////////////////////////////////////////////////////////////

double GHKcurr::_temp(void) const
{
    assert (_infosupplied() == true);
    return pTemp;
}

////////////////////////////////////////////////////////////////////////////////

double GHKcurr::_oconc(void) const
{
    assert (_infosupplied() == true);
    return pOuterConc;
}

////////////////////////////////////////////////////////////////////////////////

double GHKcurr::_iconc(void) const
{
    assert (_infosupplied() == true);
    return pInnerConc;
}

////////////////////////////////////////////////////////////////////////////////
/*
void GHKcurr::setGInfo(double g)
{
    assert (pSurfsys != 0);

    if(g < 0.0)
    {
        ostringstream os;
        os << "Conductance provided to GHKcurr::setG function can't be negative";
        throw steps::ArgErr(os.str());
    }
    pG = g;
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::setGMeasInfo(std::map<std::string, double> const & ginfo)
{
    assert (pSurfsys != 0);

    // Keeping any values that are not overwritten, allowing
    // partial information to be provided with this function.
    map<string, double>::const_iterator gi_end = ginfo.end();
    for (map<string, double>::const_iterator gi = ginfo.begin(); gi != gi_end; ++gi)
    {
        if(gi->first == "V" or gi->first == "temp"
            or gi->first == "iconc" or gi->first == "oconc")
        {
            pGInfo.erase(gi->first);
            pGInfo.insert(pair<string,double>(gi->first,gi->second));
        }
        else
        {
            ostringstream os;
            os << "Unknown key: '" << gi->first << "' in conductance";
            os << " measurement information in GHKcurr::setGMeasInfo function.";
            os << "\nAccepted keys are: 'temp', 'V', 'iconc' and 'oconc'.";
            throw steps::ArgErr(os.str());
        }
    }
}
*/
////////////////////////////////////////////////////////////////////////////////

void GHKcurr::_handleSelfDelete(void)
{
    pSurfsys->_handleGHKcurrDel(this);
    pG = 0.0;
    pValence = 0;
    pV = 0.0;
    pTemp = 0.0;
    pInnerConc = 0.0;
    pOuterConc = 0.0;
    pP = 0.0;
    pInfoSupplied = false;
    pIon = 0;
    pSurfsys = 0;
    pModel = 0;
}

////////////////////////////////////////////////////////////////////////////////

// END
