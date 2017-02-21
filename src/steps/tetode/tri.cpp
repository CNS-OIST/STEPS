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



// Standard library & STL headers.
#include <cassert>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/patchdef.hpp"
#include "steps/solver/ohmiccurrdef.hpp"
#include "steps/solver/ghkcurrdef.hpp"

#include "steps/tetode/tetode.hpp"

#include "steps/tetode/tet.hpp"
#include "steps/tetode/tri.hpp"
#include "steps/math/constants.hpp"
#include "steps/math/ghk.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace stode = steps::tetode;
namespace ssolver = steps::solver;
namespace sm = steps::math;

////////////////////////////////////////////////////////////////////////////////

stode::Tri::Tri(uint idx, steps::solver::Patchdef * patchdef, double area,
            double l0, double l1, double l2, double d0, double d1, double d2,
            int tetinner, int tetouter, int tri0, int tri1, int tri2)
: pIdx(idx)
, pPatchdef(patchdef)
, pArea(area)
, pLengths()
, pDist()
, pInnerTet(0)
, pOuterTet(0)
, pTets()
, pNextTri()
{
    assert(pPatchdef != 0);
    assert(pArea > 0.0);

    assert (l0 > 0.0 && l1 > 0.0 && l2 > 0.0);
    assert (d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0);

    pTets[0] = tetinner;
    pTets[1] = tetouter;

    pTris[0] = tri0;
    pTris[1] = tri1;
    pTris[2] = tri2;

    pNextTri[0] = 0;
    pNextTri[1] = 0;
    pNextTri[2] = 0;

    pLengths[0] = l0;
    pLengths[1] = l1;
    pLengths[2] = l2;

    pDist[0] = d0;
    pDist[1] = d1;
    pDist[2] = d2;

}

////////////////////////////////////////////////////////////////////////////////

stode::Tri::~Tri(void)
{

}

////////////////////////////////////////////////////////////////////////////////

void stode::Tri::setInnerTet(stode::Tet * t)
{
    pInnerTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void stode::Tri::setOuterTet(stode::Tet * t)
{
    pOuterTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void stode::Tri::setNextTri(uint i, stode::Tri * t)
{
    assert (i <= 2);

    pNextTri[i]= t;
}

////////////////////////////////////////////////////////////////////////////////

void stode::Tri::checkpoint(std::fstream & cp_file)
{
}

////////////////////////////////////////////////////////////////////////////////

void stode::Tri::restore(std::fstream & cp_file)
{
}

////////////////////////////////////////////////////////////////////////////////

double stode::Tri::getOhmicI(double v, steps::tetode::TetODE * solver) const
{
    double current = 0.0;
    uint nocs = patchdef()->countOhmicCurrs();
    for (uint i = 0; i < nocs; ++i)
    {
        ssolver::OhmicCurrdef * ocdef = patchdef()->ohmiccurrdef(i);
        // The next is ok because Patchdef returns local index
        //if (idx() %1000 == 0) std::cout << "\nOhmic current: " << patchdef()->ohmiccurrdef(i)->name();

        // Now need to get the states from TetODE object, and remember to convert local indices to the global ones it needs
        uint spec_gidx = patchdef()->specL2G(patchdef()->ohmiccurr_chanstate(i));
        uint patch_gidx = patchdef()->gidx();

        double n = solver->_getTriCount(pIdx, spec_gidx);
        //if (idx() %1000 == 0) std::cout << "\nN# " << i << ": " << n;

        current += (n*ocdef->getG())*(v-ocdef->getERev());
        //if (idx() %1000 == 0) std::cout << "\nCurrent# " << i << ": " << (n*ocdef->getG())*(v-ocdef->getERev());
    }
    //if (idx() %1000 == 0) std::cout << "\n";

    return current;
}


double stode::Tri::getGHKI(double v,double dt, steps::tetode::TetODE * solver) const
{
    double current = 0.0;

    uint nghkcurrs = patchdef()->countGHKcurrs();
    for (uint i =0; i < nghkcurrs; ++i)
    {
        ssolver::GHKcurrdef * ghkdef = patchdef()->ghkcurrdef(i);

        // The rate comes from the GHK flux equation. The flux equation
        // returns a single-channel current which must be converted to a rate,
        // remembering that flux can be positive or negative (bi-directional)
        const uint gidxion = ghkdef->ion();
        double voconc = ghkdef->voconc();


        // HOW TO GET THE CONCENTRATIONS?

        // Get concentrations in Molar units: convert to Mol/m^3
        double iconc = solver->_getTetConc(iTet()->idx(), gidxion)*1.0e3;
        double oconc = 0.0;

        if (voconc < 0.0)  oconc = solver->_getTetConc(oTet()->idx(), gidxion)*1.0e3;
        else  oconc = voconc*1.0e3;

        //double v = solver->getTriV(idx()); // check indices are global or local
        double T = solver->getTemp();

        double flux_per_channel = sm::GHKcurrent(ghkdef->perm(), v+ghkdef->vshift(), ghkdef->valence(),
                                     T, iconc, oconc);

        // Fetch global index of channel state
        uint cs_gidx = ghkdef->chanstate();

        double flux = solver->_getTriCount(idx(), cs_gidx) * flux_per_channel;

        current+=flux;

        //Now got to move ions
        if (ghkdef->realflux())
        {
            // Note: For a positive flux, this could be an efflux of +ve cations,
            // or an influx of -ve anions. Need to check the valence.

            // Get the rate of ion flux, remembering valence may by other than 1.
            double rt = flux/(sm::E_CHARGE * static_cast<double>(ghkdef->valence()));
            // Now a positive rate is always an efflux and a negative rate is an influx

            // rt is number of ions per second; positive is an efflux and a negative is an influx
            double count = rt*dt;

            if (voconc < 0.0) solver->_setTetCount(oTet()->idx(), gidxion, solver->_getTetCount(oTet()->idx(), gidxion)+count);
            solver->_setTetCount(iTet()->idx(), gidxion, solver->_getTetCount(iTet()->idx(), gidxion)-count);


        }

    }

    return current;

}

////////////////////////////////////////////////////////////////////////////////
//END
