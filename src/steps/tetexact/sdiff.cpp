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
#include <array>
#include <iostream>
#include <vector>

// STEPS headers.
#include "sdiff.hpp"
#include "tri.hpp"
#include "solver/patchdef.hpp"


// logging
#include <easylogging++.h>
#include "util/error.hpp"


////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;
namespace ssolver = steps::solver;
namespace smath = steps::math;

////////////////////////////////////////////////////////////////////////////////

stex::SDiff::SDiff(ssolver::Diffdef * sdef, stex::Tri * tri)
:
 pSDiffdef(sdef)
, pTri(tri)
{
    AssertLog(pSDiffdef != nullptr);
    AssertLog(pTri != nullptr);
    std::array<stex::Tri *, 3> next{pTri->nextTri(0),
                                    pTri->nextTri(1),
                                    pTri->nextTri(2)
                                    };

    ssolver::Patchdef * pdef = pTri->patchdef();
    lidxTri = pdef->specG2L(pSDiffdef->lig());

    for (uint i = 0; i < 3; ++i)
    {
        pSDiffBndDirection[i] = pTri->getSDiffBndDirection(i);
        if (next[i] != nullptr) {
            pNeighbPatchLidx[i] = next[i]->patchdef()->specG2L(pSDiffdef->lig());
        }
    }

    // Precalculate part of the scaled diffusion constant.
    uint ldidx = pTri->patchdef()->surfdiffG2L(pSDiffdef->gidx());
    double dcst = pTri->patchdef()->dcst(ldidx);
    pDcst = dcst;

    std::array<double, 3> d{0.0, 0.0, 0.0};
    for (uint i = 0; i < 3; ++i)
    {
        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTri->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr))
        {
            if (!pSDiffBndDirection[i] && next[i]->patchdef() == pTri->patchdef())
            {
                d[i] = (pTri->length(i) * dcst) / (pTri->area() * dist);
                pScaledDcst += d[i];
            }
        }
    }

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst > 0.0)
    {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + d[1] / pScaledDcst;
    }
}

////////////////////////////////////////////////////////////////////////////////

stex::SDiff::~SDiff()
= default;

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(&rExtent), sizeof(unsigned long long));
    cp_file.write(reinterpret_cast<char*>(&pFlags), sizeof(uint));

    auto n_direct_dcsts = static_cast<int>(directionalDcsts.size());
    cp_file.write(reinterpret_cast<char*>(&n_direct_dcsts), sizeof(uint));
    for (auto& item : directionalDcsts) {
        cp_file.write(reinterpret_cast<const char*>(&item.first), sizeof(uint));
        cp_file.write(reinterpret_cast<char*>(&item.second), sizeof(double));
    }


    cp_file.write(reinterpret_cast<char*>(&pScaledDcst), sizeof(double));
    cp_file.write(reinterpret_cast<char*>(&pDcst), sizeof(double));
    cp_file.write(reinterpret_cast<char*>(pCDFSelector.data()), sizeof(double) * 2);
    cp_file.write(reinterpret_cast<char*>(pSDiffBndActive.data()), sizeof(bool) * 3);
    cp_file.write(reinterpret_cast<char*>(pSDiffBndDirection.data()), sizeof(bool) * 3);
    cp_file.write(reinterpret_cast<char*>(pNeighbPatchLidx.data()), sizeof(ssolver::lidxT) * 3);

    cp_file.write(reinterpret_cast<char*>(&crData.recorded), sizeof(bool));
    cp_file.write(reinterpret_cast<char*>(&crData.pow), sizeof(int));
    cp_file.write(reinterpret_cast<char*>(&crData.pos), sizeof(unsigned));
    cp_file.write(reinterpret_cast<char*>(&crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::restore(std::fstream & cp_file)
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
    cp_file.read(reinterpret_cast<char*>(pCDFSelector.data()), sizeof(double) * 2);
    cp_file.read(reinterpret_cast<char*>(pSDiffBndActive.data()), sizeof(bool) * 3);
    cp_file.read(reinterpret_cast<char*>(pSDiffBndDirection.data()), sizeof(bool) * 3);
    cp_file.read(reinterpret_cast<char*>(pNeighbPatchLidx.data()), sizeof(ssolver::lidxT) * 3);

    cp_file.read(reinterpret_cast<char*>(&crData.recorded), sizeof(bool));
    cp_file.read(reinterpret_cast<char*>(&crData.pow), sizeof(int));
    cp_file.read(reinterpret_cast<char*>(&crData.pos), sizeof(unsigned));
    cp_file.read(reinterpret_cast<char*>(&crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::setupDeps()
{
    // We will check all KProcs of the following simulation elements:
    //   * the 'source' triangle
    //   * any neighbouring tetrahedrons
    //
    // But also in the possible 'destination' triangles (leading to
    // four different dependency lists, each containing a copy of the
    // dependencies in the 'source' tet):
    //   * any neighbouring triangles
    //   * any neighbouring tetrahedrons of these neighbouring tris
    //


    // Search for dependencies in the 'source' triangle.
    std::set<stex::KProc*> local;

    for (auto const& k :pTri->kprocs()) {
        // Check locally.
        if (k->depSpecTri(pSDiffdef->lig(), pTri)) {
            local.insert(k);
        }
    }

    // Check the neighbouring tetrahedrons.
    {
        stex::WmVol *itet[2] = {pTri->iTet(), pTri->oTet()};
        for (auto &i : itet) {
            if (i == nullptr) {
                continue;
            }
            for (auto const& k : i->kprocs()) {
                if (k->depSpecTri(pSDiffdef->lig(), pTri)) {
                    local.insert(k);
                }
            }
        }
    }

    // Search for dependencies in neighbouring triangles.
    for (uint i = 0; i < 3; ++i)
    {
        // Fetch next triangle, if it exists.
        stex::Tri * next = pTri->nextTri(i);
        if (next == nullptr) {
          continue;
        }

        // Copy local dependencies.
        std::set<stex::KProc*> local2(local.begin(), local.end());

        // Find the ones 'locally' in the next tri.
        for (auto const& k : next->kprocs()) {
            if (k->depSpecTri(pSDiffdef->lig(), next)) {
                local2.insert(k);
            }
        }

        // Fetch inner tetrahedron, if it exists.
        stex::WmVol * itet[2] = {pTri->iTet(), pTri->oTet()};
        for (auto &j : itet) {
            if (j == nullptr) {
                continue;
            }

            // Find deps.
            for (auto const & k : j->kprocs()) {
                if (k->depSpecTri(pSDiffdef->lig(), next)) {
                    local2.insert(k);
                }
            }
        }

        // Copy the set to the update vector.
        pUpdVec[i].assign(local2.begin(), local2.end());
    }

}

////////////////////////////////////////////////////////////////////////////////

bool stex::SDiff::depSpecTet(uint /*gidx*/, stex::WmVol * /*tet*/)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool stex::SDiff::depSpecTri(uint gidx, stex::Tri * tri)
{
    if (pTri != tri) { return false;
}
    if (gidx != pSDiffdef->lig()) { return false;
}
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::reset()
{
    resetExtent();

    pSDiffBndActive = {false, false, false};

    uint ldidx = pTri->patchdef()->surfdiffG2L(pSDiffdef->gidx());
    double dcst = pTri->patchdef()->dcst(ldidx);

    setDcst(dcst);

    setActive(true);

    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;

}

////////////////////////////////////////////////////////////////////////////////

double stex::SDiff::dcst(int direction)
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

void stex::SDiff::setDcst(double dcst)
{
    AssertLog(dcst >= 0.0);
    pDcst = dcst;
    directionalDcsts.clear();

    std::array<stex::Tri *, 3> next{pTri->nextTri(0),
                                    pTri->nextTri(1),
                                    pTri->nextTri(2)
                                    };

    std::array<double, 3> d{0.0, 0.0, 0.0};
    pScaledDcst = 0.0;
    for (uint i = 0; i < 3; ++i)
    {
        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTri->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr))
        {
            if ((pSDiffBndDirection[i] && pSDiffBndActive[i]) || (!pSDiffBndDirection[i] && next[i]->patchdef() == pTri->patchdef()))
            {
                d[i] = (pTri->length(i) * dcst) / (pTri->area() * dist);
            }
        }
        pScaledDcst += d[i];
    }

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pCDFSelector = {0.0, 0.0};
    }
    else
    {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + d[1] / pScaledDcst;
    }

}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::setDirectionDcst(int direction, double dcst)
{
    AssertLog(direction < 3);
    AssertLog(direction >= 0);
    AssertLog(dcst >= 0.0);
    directionalDcsts[direction] = dcst;

    // Automatically activate boundary diffusion if necessary
    if (pSDiffBndDirection[direction]) {
        pSDiffBndActive[direction] = true;
    }

    std::array<stex::Tri *, 3> next{pTri->nextTri(0),
                                    pTri->nextTri(1),
                                    pTri->nextTri(2)
                                    };

    std::array<double, 3> d{0.0, 0.0, 0.0};
    pScaledDcst = 0.0;

    for (uint i = 0; i < 3; ++i)
    {
        // if directional diffusion dcst exists use directional dcst, else use the standard one
        double dist = pTri->dist(i);

        if ((dist > 0.0) && (next[i] != nullptr))
        {
            if ((pSDiffBndDirection[i] && pSDiffBndActive[i]) || (!pSDiffBndDirection[i] && next[i]->patchdef() == pTri->patchdef()))
            {
                auto search_result = directionalDcsts.find(i);
                if (search_result != directionalDcsts.end()) {

                    d[i] = (pTri->length(i) * search_result->second) / (pTri->area() * dist);
                }
                else {
                    // This part must use the default pDcst, not the function argument dcst
                    d[i] = (pTri->length(i) * pDcst) / (pTri->area() * dist);
                }
            }
        }
        pScaledDcst += d[i];
    }

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pCDFSelector = {0.0, 0.0};
    }
    else
    {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + d[1] / pScaledDcst;
    }
}
////////////////////////////////////////////////////////////////////////////////

double stex::SDiff::rate(steps::tetexact::Tetexact * /*solver*/)
{
    if (inactive()) return 0.0;

    // Compute the rate.
    double rate = pScaledDcst * static_cast<double>(pTri->pools()[lidxTri]);
    AssertLog(std::isnan(rate) == false);

    // Return.
    return rate;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<stex::KProc*> const & stex::SDiff::apply(const rng::RNGptr &rng, double /*dt*/, double /*simtime*/)
{
    //uint lidxTet = this->lidxTet;
    // Pre-fetch some general info.


    // Apply local change.
    uint * local = pTri->pools() + lidxTri;
    bool clamped = pTri->clamped(lidxTri);

    if (clamped == false)
    {
        AssertLog(*local > 0);
    }

    // Apply change in next voxel: select a direction.
    double sel = rng->getUnfEE();

    int iSel = 0;
    for (; iSel < 2; ++iSel)
        if(sel < pCDFSelector[iSel])
            break;

    // Direction iSel.
    stex::Tri * nexttri = pTri->nextTri(iSel);
    // If there is no next tet iSel, pCDFSelector[iSel] should be zero
    // So we can assert that nextet 0 does indeed exist
    AssertLog(nexttri != nullptr);
    AssertLog(pNeighbPatchLidx[iSel] != ssolver::LIDX_UNDEFINED);

    if (nexttri->clamped(pNeighbPatchLidx[iSel]) == false)
        nexttri->incCount(pNeighbPatchLidx[iSel],1);

    if (clamped == false)
        pTri->incCount(lidxTri, -1);

    rExtent++;

    return pUpdVec[iSel];
}

////////////////////////////////////////////////////////////////////////////////

uint stex::SDiff::updVecSize() const
{
    auto maxsize = pUpdVec[0].size();
    for (uint i=1; i <= 2; ++i)
    {
        if (pUpdVec[i].size() > maxsize)
            maxsize = pUpdVec[i].size();
    }
    return maxsize;
}

////////////////////////////////////////////////////////////////////////////////


void stex::SDiff::setSDiffBndActive(uint i, bool active)
{
    AssertLog(i < 3);
    AssertLog(pSDiffBndDirection[i] == true);

    // Only need to update if the flags are changing
    if (pSDiffBndActive[i] != active)
    {

        pSDiffBndActive[i] = active;
        setDcst(pDcst);
    }

}

////////////////////////////////////////////////////////////////////////////////

bool stex::SDiff::getSDiffBndActive(uint i) const
{
    AssertLog(i < 3);
    AssertLog(pSDiffBndDirection[i] == true);

    return pSDiffBndActive[i];
}

////////////////////////////////////////////////////////////////////////////////

// END
