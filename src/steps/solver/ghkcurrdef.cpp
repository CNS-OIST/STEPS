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


// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/math/ghk.hpp"
#include "steps/model/chanstate.hpp"
#include "steps/model/ghkcurr.hpp"
#include "steps/model/spec.hpp"
#include "steps/solver/ghkcurrdef.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/types.hpp"

// logging
#include "easylogging++.h"

namespace ssolver = steps::solver;
namespace smod = steps::model;

////////////////////////////////////////////////////////////////////////////////

ssolver::GHKcurrdef::GHKcurrdef(Statedef * sd, uint gidx, smod::GHKcurr * ghk)
: pStatedef(sd)
, pIdx(gidx)
{
    AssertLog(pStatedef != nullptr);
    AssertLog(ghk != nullptr);

    pName = ghk->getID();
    pChanState = ghk->getChanState()->getID();
    pIon = ghk->getIon()->getID();
    pRealFlux = ghk->_realflux();
    pVirtual_oconc = ghk->_voconc();
    pVshift = ghk->_vshift();

    if (! ghk->_infosupplied())
    {
        std::ostringstream os;
        os << "\nPermeability not defined for GHK current object.";
        ArgErrLog(os.str());
    }

    pValence = ghk->_valence();
    AssertLog(pValence != 0);

    // Fetch the conductance measurement information from the GHKcurr object
    // in order to calculate the permeability.
    double G = ghk->_G();
    if (G > 0.0)
    {
        double V = ghk->_V();
        double temp = ghk->_temp();
        AssertLog(temp >=0.0);
        double oconc = ghk->_oconc();
        AssertLog(oconc >= 0.0);
        double iconc = ghk->_iconc();
        AssertLog(iconc >= 0.0);

        // Convert concentrations in Molar input to units for GHK equation, namely mol/m^3 or mmol/l
        oconc *= 1000.0;
        iconc *= 1000.0;

        // Calculate the permeability from the GHK flux equation.
        pPerm = steps::math::permeability(G, V, pValence, temp, iconc, oconc);
        if (pPerm != pPerm) // Implies pPerm is nan
        {
            std::ostringstream os;
            os << "\nFailed to find permeability for GHK current object, check parameters.";
            ArgErrLog(os.str());
        }
        // TODO: some check here that the permeability isn't ridiculous.
    }
    else
    {
        double perm = ghk->_P();
        AssertLog(perm>0.0);
        pPerm = perm;
    }

    uint nspecs = pStatedef->countSpecs();
    if (nspecs == 0) { return; // Would be weird, but okay.
}
    pSpec_DEP = new int[nspecs];
    std::fill_n(pSpec_DEP, nspecs, DEP_NONE);

    pSpec_VOL_DEP = new int[nspecs];
    std::fill_n(pSpec_VOL_DEP, nspecs, DEP_NONE);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::GHKcurrdef::~GHKcurrdef()
{
    if (pStatedef->countSpecs() > 0)
    {
        delete[] pSpec_DEP;
        delete[] pSpec_VOL_DEP;
    }
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::GHKcurrdef::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(&pRealFlux), sizeof(bool));
    cp_file.write(reinterpret_cast<char*>(&pVirtual_oconc), sizeof(double));
    cp_file.write(reinterpret_cast<char*>(&pPerm), sizeof(double));
    cp_file.write(reinterpret_cast<char*>(&pValence), sizeof(int));
    cp_file.write(reinterpret_cast<char*>(&pVshift), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::GHKcurrdef::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(&pRealFlux), sizeof(bool));
    cp_file.read(reinterpret_cast<char*>(&pVirtual_oconc), sizeof(double));
    cp_file.read(reinterpret_cast<char*>(&pPerm), sizeof(double));
    cp_file.read(reinterpret_cast<char*>(&pValence), sizeof(int));
    cp_file.read(reinterpret_cast<char*>(&pVshift), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::GHKcurrdef::setup()
{
    AssertLog(pSetupdone == false);

    uint chidx = pStatedef->getSpecIdx(pChanState);
    uint ionidx = pStatedef->getSpecIdx(pIon);

    pSpec_CHANSTATE = chidx;
    pSpec_ION = ionidx;

    pSpec_DEP[chidx] |= DEP_STOICH;    // TODO: come back to this and check.
    // Note: The dependency here is indirect. The flux is not modelled as
    // a 2nd order reaction ion + channel -> movement of ion
    // instead the flux comes from the GHK flux equation and is then
    // modelled as a 1st order reaction, so the channel is the only dependency
    // Update is a difference matter though, each event will change the number of
    // molecules in bordering tetrahedrons.
    // And, the concentration of ion in the inner and outer compartment
    // has an affect on the RATE
    // pSpec_DEP[ionidx] |= DEP_RATE;

    // And, the concentration of ion in the inner and outer compartment
    // has an affect on the RATE
    pSpec_VOL_DEP[ionidx] |= DEP_RATE;

    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::GHKcurrdef::chanstate() const
{
    AssertLog(pSetupdone == true);
    return pSpec_CHANSTATE;
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::GHKcurrdef::ion() const
{
    AssertLog(pSetupdone == true);
    return pSpec_ION;
}
////////////////////////////////////////////////////////////////////////////////

int ssolver::GHKcurrdef::dep(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_DEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::GHKcurrdef::dep_v(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_VOL_DEP[gidx];
}
////////////////////////////////////////////////////////////////////////////////

bool ssolver::GHKcurrdef::req(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    if (pSpec_DEP[gidx] != DEP_NONE) { return true;
}
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool ssolver::GHKcurrdef::req_v(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    if (pSpec_VOL_DEP[gidx] != DEP_NONE) { return true;
}
    return false;
}

////////////////////////////////////////////////////////////////////////////////

// END
