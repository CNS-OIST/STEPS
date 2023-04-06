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



/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

// STL headers.
#include <cassert>
#include <string>

// STEPS headers.
#include "reacdef.hpp"
#include "compdef.hpp"
#include "geom/comp.hpp"
#include "model/reac.hpp"

#include "util/error.hpp"
// logging
#include <easylogging++.h>

////////////////////////////////////////////////////////////////////////////////

namespace ssolver = steps::solver;
namespace smod = steps::model;

////////////////////////////////////////////////////////////////////////////////

ssolver::Reacdef::Reacdef(Statedef * sd, uint idx, steps::model::Reac * r)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pOrder()
, pKcst()
, pLhs()
, pRhs()
, pSetupdone(false)
, pSpec_DEP(nullptr)
, pSpec_LHS(nullptr)
, pSpec_RHS(nullptr)
, pSpec_UPD(nullptr)
, pSpec_UPD_Coll()
{
    AssertLog(pStatedef != nullptr);
    AssertLog(r != nullptr);

    pName = r->getID();
    pOrder = r->getOrder();
    pKcst = r->getKcst();
    pLhs = r->getLHS();
    pRhs = r->getRHS();

    uint nspecs = pStatedef->countSpecs();
    if (nspecs == 0) { return; // Would be weird, but okay.
}
    pSpec_DEP = new int[nspecs];
    std::fill_n(pSpec_DEP, nspecs, DEP_NONE);
    pSpec_LHS = new uint[nspecs];
    std::fill_n(pSpec_LHS, nspecs, 0);
    pSpec_RHS = new uint[nspecs];
    std::fill_n(pSpec_RHS, nspecs, 0);
    pSpec_UPD = new int[nspecs];
    std::fill_n(pSpec_UPD, nspecs, 0);

}

////////////////////////////////////////////////////////////////////////////////

ssolver::Reacdef::~Reacdef()
{
    if (pStatedef->countSpecs() > 0)
    {
        delete[] pSpec_DEP;
        delete[] pSpec_LHS;
        delete[] pSpec_RHS;
        delete[] pSpec_UPD;
    }
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Reacdef::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(&pKcst), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Reacdef::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(&pKcst), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Reacdef::setup()
{
    AssertLog(pSetupdone == false);

    // first copy the information about the reaction stoichiometry from Reac object
    for (auto const& l: pLhs) {
        uint sidx = pStatedef->getSpecIdx(l);
        pSpec_LHS[sidx] += 1;
    }
    for (auto const& r: pRhs) {
        uint sidx = pStatedef->getSpecIdx(r);
        pSpec_RHS[sidx] += 1;
    }

    // Now set up the update vector
    const uint nspecs = pStatedef->countSpecs();
    for (uint i = 0; i < nspecs; ++i)
    {
        auto lhs = static_cast<int>(pSpec_LHS[i]);
        auto rhs = static_cast<int>(pSpec_RHS[i]);
        int aux = pSpec_UPD[i] = (rhs - lhs);
        if (lhs != 0) { pSpec_DEP[i] |= DEP_STOICH;
}
        if (aux != 0) pSpec_UPD_Coll.push_back(i);
    }

    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

std::string const ssolver::Reacdef::name() const
{
    return pName;
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Reacdef::order() const
{
    return pOrder;
}

////////////////////////////////////////////////////////////////////////////////

double ssolver::Reacdef::kcst() const
{
    return pKcst;
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Reacdef::lhs(uint gidx) const
{
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_LHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Reacdef::dep(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_DEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Reacdef::rhs(uint gidx) const
{
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_RHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Reacdef::upd(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_UPD[gidx];
}

////////////////////////////////////////////////////////////////////////////////

bool ssolver::Reacdef::reqspec(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    if (pSpec_DEP[gidx] != DEP_NONE) { return true;
}
    if (pSpec_RHS[gidx] != 0) { return true;
}
    return false;
}

////////////////////////////////////////////////////////////////////////////////

// END
