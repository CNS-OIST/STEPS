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
#include <iostream>
#include <sstream>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/model/spec.hpp"
#include "steps/model/vdepsreac.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/types.hpp"
#include "steps/solver/vdepsreacdef.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace ssolver = steps::solver;
namespace smod = steps::model;

////////////////////////////////////////////////////////////////////////////////

ssolver::VDepSReacdef::VDepSReacdef(Statedef * sd, uint idx, smod::VDepSReac * vdsr)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pIlhs()
, pOlhs()
, pSlhs()
, pIrhs()
, pOrhs()
, pSrhs()
, pSurface_surface(true)

// pOrient?

, pSpec_I_DEP(nullptr)
, pSpec_S_DEP(nullptr)
, pSpec_O_DEP(nullptr)
, pSpec_I_LHS(nullptr)
, pSpec_S_LHS(nullptr)
, pSpec_O_LHS(nullptr)
, pSpec_I_RHS(nullptr)
, pSpec_S_RHS(nullptr)
, pSpec_O_RHS(nullptr)
, pSpec_I_UPD(nullptr)
, pSpec_S_UPD(nullptr)
, pSpec_O_UPD(nullptr)
, pSpec_I_UPD_Coll()
, pSpec_S_UPD_Coll()
, pSpec_O_UPD_Coll()
, pVMin(0.0)
, pVMax(0.0)
, pDV(0.0)
, pVKTab(nullptr)
, pReqInside(false)
, pReqOutside(false)
{
    AssertLog(pStatedef != 0);
    AssertLog(vdsr != 0);

    pName = vdsr->getID();
    pOrder = vdsr->getOrder();

    if (pOrder == 0)
    {
        std::ostringstream os;
        os << "Model contains zero-order voltage-dependent surface reaction, which are not permitted. ";
        ArgErrLog(os.str());
    }

    // Copy rate information from model object
    pVMin = vdsr->_getVMin();
    pVMax = vdsr->_getVMax();
    pDV = vdsr->_getDV();
    uint tablesize = vdsr->_getTablesize();
    AssertLog(tablesize == static_cast<uint>(std::floor((pVMax - pVMin) / pDV)) + 1);

    pVKTab = new double[tablesize];
    // Just temporarily store the pointer:
    //double * k = vdsr->_getK();

    for (uint i = 0; i < tablesize; ++i)
    {
        pVKTab[i] = vdsr->_getK()[i];
    }

    pIlhs = vdsr->getILHS();
    pOlhs = vdsr->getOLHS();
    pSlhs = vdsr->getSLHS();

    pIrhs = vdsr->getIRHS();
    pOrhs = vdsr->getORHS();
    pSrhs = vdsr->getSRHS();

    if (vdsr->getInner() == true) { pOrient = VDepSReacdef::INSIDE;
    } else { pOrient = VDepSReacdef::OUTSIDE;
}

    uint nspecs = pStatedef->countSpecs();
    if (nspecs == 0) { return; // Would be weird, but okay.
}
    pSpec_S_DEP = new int[nspecs];
    std::fill_n(pSpec_S_DEP, nspecs, DEP_NONE);
    pSpec_S_LHS = new uint[nspecs];
    std::fill_n(pSpec_S_LHS, nspecs, 0);
    if (pOrient ==  VDepSReacdef::INSIDE)
    {
        pSpec_I_DEP = new int[nspecs];
        std::fill_n(pSpec_I_DEP, nspecs, DEP_NONE);
        pSpec_I_LHS = new uint[nspecs];
        std::fill_n(pSpec_I_LHS, nspecs, 0);
    }
    else
    {
        pSpec_O_DEP = new depT[nspecs];
        std::fill_n(pSpec_O_DEP, nspecs, DEP_NONE);
        pSpec_O_LHS = new uint[nspecs];
        std::fill_n(pSpec_O_LHS, nspecs, 0);
    }
    pSpec_I_RHS = new uint[nspecs];
    std::fill_n(pSpec_I_RHS, nspecs, 0);
    pSpec_S_RHS = new uint[nspecs];
    std::fill_n(pSpec_S_RHS, nspecs, 0);
    pSpec_O_RHS = new uint[nspecs];
    std::fill_n(pSpec_O_RHS, nspecs, 0);
    pSpec_I_UPD = new int[nspecs];
    std::fill_n(pSpec_I_UPD, nspecs, 0);
    pSpec_S_UPD = new int[nspecs];
    std::fill_n(pSpec_S_UPD, nspecs, 0);
    pSpec_O_UPD = new int[nspecs];
    std::fill_n(pSpec_O_UPD, nspecs, 0);

}

////////////////////////////////////////////////////////////////////////////////

ssolver::VDepSReacdef::~VDepSReacdef()
{
    if (pStatedef->countSpecs() > 0)
    {
        if (pOrient == VDepSReacdef::INSIDE)
        {
            delete[] pSpec_I_DEP;
            delete[] pSpec_I_LHS;
        }
        else
        {
            delete[] pSpec_O_DEP;
            delete[] pSpec_O_LHS;

        }
        delete[] pSpec_S_DEP;
        delete[] pSpec_S_LHS;
        delete[] pSpec_I_RHS;
        delete[] pSpec_S_RHS;
        delete[] pSpec_O_RHS;
        delete[] pSpec_I_UPD;
        delete[] pSpec_S_UPD;
        delete[] pSpec_O_UPD;

    }

    delete[] pVKTab;

}

////////////////////////////////////////////////////////////////////////////////

void ssolver::VDepSReacdef::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(&pVMin), sizeof(double));
    cp_file.write(reinterpret_cast<char*>(&pVMax), sizeof(double));
    cp_file.write(reinterpret_cast<char*>(&pDV), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::VDepSReacdef::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(&pVMin), sizeof(double));
    cp_file.read(reinterpret_cast<char*>(&pVMax), sizeof(double));
    cp_file.read(reinterpret_cast<char*>(&pDV), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::VDepSReacdef::setup()
{
    AssertLog(pSetupdone == false);

    if (outside()) { AssertLog(pIlhs.empty()); }
    else if (inside()) { AssertLog(pOlhs.empty()); }
    else AssertLog(false);

    for (auto const& ol: pOlhs) {
        pSurface_surface = false;
        uint sidx = pStatedef->getSpecIdx(ol);
        pSpec_O_LHS[sidx] += 1;
    }

    for (auto const& il: pIlhs) {
        pSurface_surface = false;
        uint sidx = pStatedef->getSpecIdx(il);
        pSpec_I_LHS[sidx] += 1;
    }

    for (auto const& sl: pSlhs) {
        uint sidx = pStatedef->getSpecIdx(sl);
        pSpec_S_LHS[sidx] += 1;
    }

    for (auto const& ir: pIrhs) {
        uint sidx = pStatedef->getSpecIdx(ir);
        pSpec_I_RHS[sidx] += 1;
    }

    for (auto const& sr: pSrhs) {
        uint sidx = pStatedef->getSpecIdx(sr);
        pSpec_S_RHS[sidx] += 1;
    }

    for (auto const& orh: pOrhs) {
        uint sidx = pStatedef->getSpecIdx(orh);
        pSpec_O_RHS[sidx] += 1;
    }

    // Now set up the update vector
    uint nspecs = pStatedef->countSpecs();
    // Deal with surface.
    for (uint i = 0; i < nspecs; ++i)
    {
        auto lhs = static_cast<int>(pSpec_S_LHS[i]);
        auto rhs = static_cast<int>(pSpec_S_RHS[i]);
        int aux = pSpec_S_UPD[i] = (rhs - lhs);
        if (lhs != 0) { pSpec_S_DEP[i] |= DEP_STOICH;
}
        if (aux != 0) pSpec_S_UPD_Coll.push_back(i);
    }

    // Deal with inside.
    for (uint i = 0; i < nspecs; ++i)
    {
        int lhs = (inside() ? static_cast<int>(pSpec_I_LHS[i]) : 0);
        auto rhs = static_cast<int>(pSpec_I_RHS[i]);
        int aux = pSpec_I_UPD[i] = (rhs - lhs);
        if (lhs != 0) { pSpec_I_DEP[i] |= DEP_STOICH;
}
        if (aux != 0) pSpec_I_UPD_Coll.push_back(i);
    }

    // Deal with outside.
    for (uint i = 0; i < nspecs; ++i)
    {
        int lhs = (outside() ? static_cast<int>(pSpec_O_LHS[i]) : 0);
        auto rhs = static_cast<int>(pSpec_O_RHS[i]);
        int aux = pSpec_O_UPD[i] = (rhs - lhs);
        if (lhs != 0) { pSpec_O_DEP[i] |= DEP_STOICH;
}
        if (aux != 0) pSpec_O_UPD_Coll.push_back(i);
    }

    // This has to come before the final loop
    pSetupdone = true;

    for (uint i = 0; i < nspecs; ++i)
    {
        if (reqspec_I(i) == true) { pReqInside  = true;
}
        if (reqspec_O(i) == true) { pReqOutside  = true;
}
    }

}

////////////////////////////////////////////////////////////////////////////////

double ssolver::VDepSReacdef::getVDepK(double v) const
{
    AssertLog(pSetupdone == true);
    AssertLog(pVKTab != nullptr);
    if (v > pVMax)
    {
        std::ostringstream os;
        os << "Voltage to VDepSReac::getVDepRate higher than maximum: ";
        os << v << " > " << pVMax;
        ProgErrLog(os.str());
    }
    if (v < pVMin)
    {
        std::ostringstream os;
        os << "Voltage to VDepSReac::getVDepRate lower than minimum: ";
        os << v << " < " << pVMin;
        ProgErrLog(os.str());
    }

    double v2 = ((v - pVMin) / pDV);
    double lv = floor(v2);
    auto lvidx = static_cast<uint>(lv);
    uint uvidx = static_cast<uint>(ceil(v2));
    double r = v2-lv;

    return (((1.0 - r) * pVKTab[lvidx]) + (r * pVKTab[uvidx]));

}

////////////////////////////////////////////////////////////////////////////////
/*
bool ssolver::VDepSReacdef::reqInside() const
{
    AssertLog(pSetupdone == true);

    // This can be checked by seeing if DEP_I or RHS_I is non-zero
    // for any species.
    uint nspecs = pStatedef->countSpecs();
    for (uint i = 0; i < nspecs; ++i) if (reqspec_I(i) == true) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool ssolver::VDepSReacdef::reqOutside() const
{
    AssertLog(pSetupdone == true);

    // This can be checked by seeing if DEP_O or RHS_O is non-zero
    // for any species.
    uint nspecs = pStatedef->countSpecs();
    for (uint i = 0; i < nspecs; ++i)
        if (reqspec_O(i) == true) return true;
    return false;
}
*/
////////////////////////////////////////////////////////////////////////////////

uint ssolver::VDepSReacdef::lhs_I(uint gidx) const
{
    if (outside()) { return 0;
}
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_I_LHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::VDepSReacdef::lhs_S(uint gidx) const
{
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_S_LHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::VDepSReacdef::lhs_O(uint gidx) const
{
    if (inside()) { return 0;
}
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_O_LHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::VDepSReacdef::dep_I(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    if (outside()) { return DEP_NONE;
}
    return pSpec_I_DEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::VDepSReacdef::dep_S(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_S_DEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::VDepSReacdef::dep_O(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    if (inside()) { return DEP_NONE;
}
    return pSpec_O_DEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::VDepSReacdef::rhs_I(uint gidx) const
{
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_I_RHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::VDepSReacdef::rhs_S(uint gidx) const
{
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_S_RHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::VDepSReacdef::rhs_O(uint gidx) const
{
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_O_RHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::VDepSReacdef::upd_I(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_I_UPD[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::VDepSReacdef::upd_S(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_S_UPD[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::VDepSReacdef::upd_O(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_O_UPD[gidx];
}

////////////////////////////////////////////////////////////////////////////////

bool ssolver::VDepSReacdef::reqspec_I(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    if (inside()) {
        if (pSpec_I_DEP[gidx] != DEP_NONE) { return true;
}
}
    if (pSpec_I_RHS[gidx] != 0) { return true;
}
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool ssolver::VDepSReacdef::reqspec_S(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    if (pSpec_S_DEP[gidx] != DEP_NONE) { return true;
}
    return pSpec_S_RHS[gidx] != 0;
}

////////////////////////////////////////////////////////////////////////////////

bool ssolver::VDepSReacdef::reqspec_O(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    if (outside()) {
        if (pSpec_O_DEP[gidx] != DEP_NONE) { return true;
}
}
    return pSpec_O_RHS[gidx] != 0;
}

////////////////////////////////////////////////////////////////////////////////
// END
