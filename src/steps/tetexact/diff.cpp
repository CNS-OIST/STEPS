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


// Standard library & STL headers.
#include <iostream>
#include <vector>
// STEPS headers.
#include "diff.hpp"
#include "tet.hpp"
#include "tri.hpp"
#include "solver/compdef.hpp"

// logging
#include <easylogging++.h>
#include "util/common.h"
////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;
namespace ssolver = steps::solver;
namespace smath = steps::math;

////////////////////////////////////////////////////////////////////////////////

stex::Diff::Diff(ssolver::Diffdef * ddef, stex::Tet * tet)
:
pDiffdef(ddef)
, pTet(tet)
{
    AssertLog(pDiffdef != nullptr);
    AssertLog(pTet != nullptr);
    std::array<stex::Tet*, 4> next{pTet->nextTet(0),
                                    pTet->nextTet(1),
                                    pTet->nextTet(2),
                                    pTet->nextTet(3)
                                    };

    ssolver::Compdef * cdef = pTet->compdef();
    lidxTet = cdef->specG2L(pDiffdef->lig());

    // Precalculate part of the scaled diffusion constant.
    uint ldidx = pTet->compdef()->diffG2L(pDiffdef->gidx());
    double dcst = pTet->compdef()->dcst(ldidx);
    pDcst = dcst;

    std::array<double, 4> d{ 0.0, 0.0, 0.0, 0.0 };

    for (uint i = 0; i < 4; ++i) {
        pDiffBndDirection[i] = pTet->getDiffBndDirection(i);
        if (next[i] != nullptr) {
            pNeighbCompLidx[i] = next[i]->compdef()->specG2L(pDiffdef->lig());
        }

        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTet->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr))
        {
            // d[i] only need to set if
            // 1) not towards a boundary, and
            // 2) next[i] in the same compartment as pTet
            // d[i] changes when setDiffBndActive() is called
            // and pDiffBndActive[i] becomes active
            if (!pDiffBndDirection[i] && next[i]->compdef() == cdef)
            {
                d[i] = (pTet->area(i) * dcst) / (pTet->vol() * dist);
                pScaledDcst += d[i];
            }
        }
    }

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst > 0.0) {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + d[1] / pScaledDcst;
        pCDFSelector[2] = pCDFSelector[1] + d[2] / pScaledDcst;
    }
}

////////////////////////////////////////////////////////////////////////////////

stex::Diff::~Diff()
= default;

////////////////////////////////////////////////////////////////////////////////

void stex::Diff::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(&rExtent), sizeof(unsigned long long));
    cp_file.write(reinterpret_cast<char*>(&pFlags), sizeof(uint));

    auto n_direct_dcsts = directionalDcsts.size();
    cp_file.write(reinterpret_cast<char*>(&n_direct_dcsts), sizeof(uint));
    for (auto& item : directionalDcsts) {
        cp_file.write(reinterpret_cast<const char*>(&item.first), sizeof(uint));
        cp_file.write(reinterpret_cast<char*>(&item.second), sizeof(double));
    }


    cp_file.write(reinterpret_cast<char*>(&pScaledDcst), sizeof(double));
    cp_file.write(reinterpret_cast<char*>(&pDcst), sizeof(double));
    cp_file.write(reinterpret_cast<char*>(pDiffBndActive.data()), sizeof(bool) * 4);
    cp_file.write(reinterpret_cast<char*>(pDiffBndDirection.data()), sizeof(bool) * 4);
    cp_file.write(reinterpret_cast<char*>(pNeighbCompLidx.data()), sizeof(ssolver::lidxT) * 4);
    cp_file.write(reinterpret_cast<char*>(pCDFSelector.data()), sizeof(double) * 3);
    cp_file.write(reinterpret_cast<char*>(&crData.recorded), sizeof(bool));
    cp_file.write(reinterpret_cast<char*>(&crData.pow), sizeof(int));
    cp_file.write(reinterpret_cast<char*>(&crData.pos), sizeof(unsigned));
    cp_file.write(reinterpret_cast<char*>(&crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Diff::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(&rExtent), sizeof(unsigned long long));
    cp_file.read(reinterpret_cast<char*>(&pFlags), sizeof(uint));

    uint n_direct_dcsts = 0;
    cp_file.read(reinterpret_cast<char*>(&n_direct_dcsts), sizeof(uint));
    for (uint i = 0; i < n_direct_dcsts; i++) {
        uint id = 0;
        double value = 0.0;
        cp_file.read(reinterpret_cast<char*>(&id), sizeof(uint));
        cp_file.read(reinterpret_cast<char*>(&value), sizeof(double));
        directionalDcsts[id] = value;
    }

    cp_file.read(reinterpret_cast<char*>(&pScaledDcst), sizeof(double));
    cp_file.read(reinterpret_cast<char*>(&pDcst), sizeof(double));
    cp_file.read(reinterpret_cast<char*>(pDiffBndActive.data()), sizeof(bool) * 4);
    cp_file.read(reinterpret_cast<char*>(pDiffBndDirection.data()), sizeof(bool) * 4);
    cp_file.read(reinterpret_cast<char*>(pNeighbCompLidx.data()), sizeof(ssolver::lidxT) * 4);
    cp_file.read(reinterpret_cast<char*>(pCDFSelector.data()), sizeof(double) * 3);
    cp_file.read(reinterpret_cast<char*>(&crData.recorded), sizeof(bool));
    cp_file.read(reinterpret_cast<char*>(&crData.pow), sizeof(int));
    cp_file.read(reinterpret_cast<char*>(&crData.pos), sizeof(unsigned));
    cp_file.read(reinterpret_cast<char*>(&crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Diff::setupDeps()
{
    // We will check all KProcs of the following simulation elements:
    //   * the 'source' tetrahedron
    //   * any neighbouring triangles
    //
    // But also in the possible 'destination' tetrahedrons (leading to
    // four different dependency lists, each containing a copy of the
    // dependencies in the 'source' tet):
    //   * any neighbouring tetrahedrons
    //   * any neighbouring triangles of these neighbouring tets
    //
    // Since there can be no diffusion between tetrahedrons blocked by
    // a triangle, there is no need to filter out duplicate dependent
    // kprocs.

    // Search for dependencies in the 'source' tetrahedron.
    std::set<stex::KProc*> local;

    for (auto const& k: pTet->kprocs()) {
        // Check locally.
        if (k->depSpecTet(pDiffdef->lig(), pTet) == true) {
            local.insert(k);
        }
    }
    // Check the neighbouring triangles.
    for (uint i = 0; i < 4; ++i)
    {
        stex::Tri * next = pTet->nextTri(i);
        if (next == nullptr) {
            continue;
        }
        for (auto const& k: next->kprocs()) {
            if (k->depSpecTet(pDiffdef->lig(), pTet) == true) {
                local.insert(k);
            }
        }
    }

    // Search for dependencies in neighbouring tetrahedra.
    for (uint i = 0; i < 4; ++i)
    {
        // Fetch next tetrahedron, if it exists.
        stex::Tet * next = pTet->nextTet(i);
        if (next == nullptr) {
            continue;
        }
        if (pTet->nextTri(i) != nullptr) {
            continue;
        }

        // Copy local dependencies.
        std::set<stex::KProc*> local2(local.begin(), local.end());

        // Find the ones 'locally' in the next tet.
        for (auto const& k: next->kprocs()) {
            if (k->depSpecTet(pDiffdef->lig(), next) == true) {
                local2.insert(k);
            }
        }

        // Find deps in neighbouring triangles in the next tet.
        // As said before, this cannot logically include the shared
        // triangle.
        for (uint j = 0; j < 4; ++j)
        {
            // Fetch next triangle, if it exists.
            stex::Tri * next2 = next->nextTri(j);
            if (next2 == nullptr) { continue;
        }

            // Find deps.
            for (auto const& k: next2->kprocs()) {
                if (k->depSpecTet(pDiffdef->lig(), next) == true) {
                    local2.insert(k);
                }
            }
        }

        // Copy the set to the update vector.
        pUpdVec[i].assign(local2.begin(), local2.end());
    }
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Diff::depSpecTet(uint gidx, stex::WmVol * tet)
{
    if (pTet != tet) return false;
    if (gidx != pDiffdef->lig()) { return false;
}
    return true;
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Diff::depSpecTri(uint /*gidx*/, stex::Tri * /*tri*/)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Diff::reset()
{
    resetExtent();

    // NOTE: These must become the dcst calculation for obvious reasons
    pDiffBndActive = {false, false, false, false};

    uint ldidx = pTet->compdef()->diffG2L(pDiffdef->gidx());
    double dcst = pTet->compdef()->dcst(ldidx);

    // directional dcst will also be clear by setDcst
    setDcst(dcst);

    setActive(true);

    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;

}

////////////////////////////////////////////////////////////////////////////////

void stex::Diff::setDiffBndActive(uint i, bool active)
{
    AssertLog(i < 4);
    AssertLog(pDiffBndDirection[i]);

    // Only need to update if the flags are changing
    if (pDiffBndActive[i] != active)
    {
        pDiffBndActive[i] = active;
        setDcst(pDcst);
    }

}

////////////////////////////////////////////////////////////////////////////////

bool stex::Diff::getDiffBndActive(uint i) const
{
    AssertLog(i < 4);
    AssertLog(pDiffBndDirection[i]);

    return pDiffBndActive[i];
}

////////////////////////////////////////////////////////////////////////////////

double stex::Diff::dcst(int direction)
{
    auto search_result = directionalDcsts.find(direction);
    if (search_result != directionalDcsts.end()) {
        return search_result->second;
    }
    else {
        return pDcst;
    }
}

////////////////////////////////////////////////////////////////////////////////

void stex::Diff::setDcst(double dcst)
{
    AssertLog(dcst >= 0.0);
    pDcst = dcst;
    directionalDcsts.clear();

    std::array<stex::Tet *, 4> next{pTet->nextTet(0),
                                        pTet->nextTet(1),
                                        pTet->nextTet(2),
                                        pTet->nextTet(3)
                                        };

    std::array<double, 4> d{ 0.0, 0.0, 0.0, 0.0 };
    pScaledDcst = 0.0;

    for (uint i = 0; i < 4; ++i)
    {
        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTet->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr))
        {
            if ((pDiffBndDirection[i] && pDiffBndActive[i]) || (!pDiffBndDirection[i] && next[i]->compdef() == pTet->compdef()))
            {
                d[i] = (pTet->area(i) * dcst) / (pTet->vol() * dist);
            }
        }
        pScaledDcst += d[i];
    }

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pCDFSelector = {0.0, 0.0, 0.0};
    }
    else
    {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + d[1] / pScaledDcst;
        pCDFSelector[2] = pCDFSelector[1] + d[2] / pScaledDcst;
    }
}

////////////////////////////////////////////////////////////////////////////////

void stex::Diff::setDirectionDcst(int direction, double dcst)
{
    AssertLog(direction < 4);
    AssertLog(direction >= 0);
    AssertLog(dcst >= 0.0);
    directionalDcsts[direction] = dcst;

    // Automatically activate boundary diffusion if necessary
    if (pDiffBndDirection[direction] == true) {
        pDiffBndActive[direction] = true;
    }

    std::array<stex::Tet *, 4> next{pTet->nextTet(0),
                                        pTet->nextTet(1),
                                        pTet->nextTet(2),
                                        pTet->nextTet(3)
                                        };

    std::array<double, 4> d{ 0.0, 0.0, 0.0, 0.0 };
    pScaledDcst = 0.0;

    for (uint i = 0; i < 4; ++i)
    {
        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTet->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr))
        {
            if ((pDiffBndDirection[i] && pDiffBndActive[i]) || (!pDiffBndDirection[i] && next[i]->compdef() == pTet->compdef()))
            {
                auto search_result = directionalDcsts.find(i);
                if (search_result != directionalDcsts.end())
                    d[i] = (pTet->area(i) * search_result->second) / (pTet->vol() * dist);
                else
                    d[i] = (pTet->area(i) * pDcst) / (pTet->vol() * dist);
            }
        }
        pScaledDcst += d[i];
    }

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pCDFSelector = {0.0, 0.0, 0.0};
    }
    else
    {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + d[1] / pScaledDcst;
        pCDFSelector[2] = pCDFSelector[1] + d[2] / pScaledDcst;
    }
}

////////////////////////////////////////////////////////////////////////////////

double stex::Diff::rate(steps::tetexact::Tetexact * /*solver*/)
{
    if (inactive()) return 0.0;

    // Compute the rate.
    double rate = (pScaledDcst) * static_cast<double>(pTet->pools()[lidxTet]);
    AssertLog(std::isnan(rate) == false);

    // Return.
    return rate;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<stex::KProc*> const & stex::Diff::apply(const rng::RNGptr &rng, double /*dt*/, double /*simtime*/)
{
    //uint lidxTet = this->lidxTet;
    // Pre-fetch some general info.


    // Apply local change.
    auto * local = pTet->pools().data() + lidxTet;
    bool clamped = pTet->clamped(lidxTet);

    if (clamped == false)
    {
        AssertLog(*local > 0);
    }

    // Apply change in next voxel: select a direction.
    double sel = rng->getUnfEE();

    uint iSel = 0;
    for (; iSel < 3; ++iSel)
        if(sel < pCDFSelector[iSel])
            break;

    // Direction iSel.
    stex::Tet * nexttet = pTet->nextTet(iSel);
    // If there is no next tet 0, pCDFSelector[0] should be zero
    // So we can assert that nextet 0 does indeed exist
    AssertLog(nexttet != nullptr);
    AssertLog(pNeighbCompLidx[iSel] != solver::LIDX_UNDEFINED);

    if (nexttet->clamped(pNeighbCompLidx[iSel]) == false)
        nexttet->incCount(pNeighbCompLidx[iSel],1);

    if (clamped == false)
        pTet->incCount(lidxTet, -1);

    rExtent++;

    return pUpdVec[iSel];
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Diff::updVecSize() const
{
    auto maxsize = pUpdVec[0].size();
    for (uint i=1; i <= 3; ++i)
    {
        if (pUpdVec[i].size() > maxsize)
            maxsize = pUpdVec[i].size();
    }
    return maxsize;
}

////////////////////////////////////////////////////////////////////////////////

// END
