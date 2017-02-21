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
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <memory>
#include <queue>
#include <fstream>
#include <iomanip>

// STEPS headers.
#include "steps/common.h"
#include "steps/tetexact/tetexact.hpp"
#include "steps/tetexact/kproc.hpp"
#include "steps/tetexact/tet.hpp"
#include "steps/tetexact/tri.hpp"
#include "steps/tetexact/reac.hpp"
#include "steps/tetexact/sreac.hpp"
#include "steps/tetexact/diff.hpp"
#include "steps/tetexact/sdiff.hpp"
#include "steps/tetexact/comp.hpp"
#include "steps/tetexact/patch.hpp"
#include "steps/tetexact/wmvol.hpp"
#include "steps/tetexact/ghkcurr.hpp"
#include "steps/tetexact/vdeptrans.hpp"
#include "steps/tetexact/vdepsreac.hpp"
#include "steps/tetexact/diffboundary.hpp"
#include "steps/math/constants.hpp"
#include "steps/math/point.hpp"
#include "steps/error.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/compdef.hpp"
#include "steps/solver/patchdef.hpp"
#include "steps/solver/reacdef.hpp"
#include "steps/solver/sreacdef.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/solver/chandef.hpp"
#include "steps/solver/ghkcurrdef.hpp"
#include "steps/solver/ohmiccurrdef.hpp"
#include "steps/solver/vdeptransdef.hpp"
#include "steps/solver/vdepsreacdef.hpp"
#include "steps/solver/types.hpp"
#include "steps/geom/tetmesh.hpp"
#include "steps/util/distribute.hpp"

#include "steps/solver/efield/efield.hpp"
#include "steps/solver/efield/dVsolver.hpp"

#include "third_party/easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace ssolver = steps::solver;
namespace smath = steps::math;
namespace stex = steps::tetexact;

using steps::math::point3d;

////////////////////////////////////////////////////////////////////////////////

void stex::schedIDXSet_To_Vec(stex::SchedIDXSet const & s, stex::SchedIDXVec & v)
{
    v.resize(s.size());
    std::copy(s.begin(), s.end(), v.begin());
}

////////////////////////////////////////////////////////////////////////////////

stex::Tetexact::Tetexact(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r,
                         int calcMembPot)
: API(m, g, r)
, pMesh(0)
, pKProcs()
, pComps()
, pCompMap()
, pPatches()
, pDiffBoundaries()
, pTets()
, pTris()
, pWmVols()
, pA0(0.0)
//, pBuilt(false)
, pEFoption(static_cast<EF_solver>(calcMembPot))
, pTemp(0.0)
, pEFDT(1.0e-5)
, pEFNVerts(0)
, pEFVerts(0)
, pEFNTris(0)
, pEFTris(0)
, pEFTris_vec(0)
, pEFNTets(0)
, pEFTets(0)
, pEFVert_GtoL()
, pEFTri_GtoL()
, pEFTet_GtoL()
, pEFTri_LtoG()
{
    if (rng() == 0)
    {
        std::ostringstream os;
        os << "No RNG provided to solver initializer function";
        throw steps::ArgErr(os.str());
    }
    
    // All initialization code now in _setup() to allow EField solver to be
    // derived and create EField local objects within the constructor
    _setup();
}

////////////////////////////////////////////////////////////////////////////////

stex::Tetexact::~Tetexact(void)
{
    for (auto c: pComps) delete c;
    for (auto p: pPatches) delete p;
    for (auto db: pDiffBoundaries) delete db;
    for (auto wvol: pWmVols) delete wvol;
    for (auto t: pTets) delete t;
    for (auto t: pTris) delete t;
    for (auto g: nGroups) {
        g->free_indices();
        delete g;
    }
    for (auto g: pGroups) {
        g->free_indices();
        delete g;
    }

    if (efflag())
    {
        delete[] pEFVerts;
        delete[] pEFTris;
        delete[] pEFTets;
        delete[] pEFVert_GtoL;
        delete[] pEFTri_GtoL;
        delete[] pEFTet_GtoL;
        delete[] pEFTri_LtoG;
    }
}

///////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::checkpoint(std::string const & file_name)
{
    std::cout << "Checkpoint to " << file_name  << "...";
    std::fstream cp_file;

    cp_file.open(file_name.c_str(),
                std::fstream::out | std::fstream::binary | std::fstream::trunc);

    statedef()->checkpoint(cp_file);

    CompPVecCI comp_e = pComps.end();
    for (CompPVecCI c = pComps.begin(); c != comp_e; ++c) (*c)->checkpoint(cp_file);

    PatchPVecCI patch_e = pPatches.end();
    for (PatchPVecCI p = pPatches.begin(); p != patch_e; ++p) (*p)->checkpoint(cp_file);

    DiffBoundaryPVecCI db_e = pDiffBoundaries.end();
    for (DiffBoundaryPVecCI db = pDiffBoundaries.begin(); db != db_e; ++db) {
        (*db)->checkpoint(cp_file);
    }

    WmVolPVecCI wmv_e = pWmVols.end();
    for (WmVolPVecCI wmv = pWmVols.begin(); wmv != wmv_e; ++wmv)
    {
        if ((*wmv) != 0) {
            (*wmv)->checkpoint(cp_file);
        }
    }


    TetPVecCI tet_e = pTets.end();
    for (TetPVecCI t = pTets.begin(); t != tet_e; ++t)
    {
        if ((*t) != 0) {
        (*t)->checkpoint(cp_file);
        }
    }

    TriPVecCI tri_e = pTris.end();
    for (TriPVecCI t = pTris.begin(); t != tri_e; ++t)
    {
        if ((*t) != 0) {
            (*t)->checkpoint(cp_file);
        }
    }

    KProcPVecCI e = pKProcs.end();
    for (KProcPVecCI i = pKProcs.begin(); i != e; ++i) (*i)->checkpoint(cp_file);

    if (efflag()) {
        cp_file.write((char*)&pTemp, sizeof(double));
        cp_file.write((char*)&pEFDT, sizeof(double));
        pEField->checkpoint(cp_file);
    }

    cp_file.write((char*)&nEntries, sizeof(uint));

    // checkpoint CR SSA

    cp_file.write((char*)&pSum, sizeof(double));
    cp_file.write((char*)&nSum, sizeof(double));
    cp_file.write((char*)&pA0, sizeof(double));

    uint n_ngroups = nGroups.size();
    uint n_pgroups = pGroups.size();

    cp_file.write((char*)&n_ngroups, sizeof(uint));
    cp_file.write((char*)&n_pgroups, sizeof(uint));

    for (uint i = 0; i < n_ngroups; i++) {
        CRGroup* group = nGroups[i];
        cp_file.write((char*)&(group->capacity), sizeof(unsigned));
        cp_file.write((char*)&(group->size), sizeof(unsigned));
        cp_file.write((char*)&(group->max), sizeof(double));
        cp_file.write((char*)&(group->sum), sizeof(double));

        for (uint j = 0; j < group->size; j++) {
            uint idx = group->indices[j]->schedIDX();
            cp_file.write((char*)&idx, sizeof(uint));
        }
    }

    for (uint i = 0; i < n_pgroups; i++) {
        CRGroup* group = pGroups[i];
        cp_file.write((char*)&(group->capacity), sizeof(unsigned));
        cp_file.write((char*)&(group->size), sizeof(unsigned));
        cp_file.write((char*)&(group->max), sizeof(double));
        cp_file.write((char*)&(group->sum), sizeof(double));

        for (uint j = 0; j < group->size; j++) {
            uint idx = group->indices[j]->schedIDX();
            cp_file.write((char*)&idx, sizeof(uint));
        }
    }

    cp_file.close();
    std::cout << "complete.\n";
}

///////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::restore(std::string const & file_name)
{
    std::cout << "Restore from " << file_name << "...";
    std::fstream cp_file;

    cp_file.open(file_name.c_str(),
                std::fstream::in | std::fstream::binary);

    cp_file.seekg(0);

    statedef()->restore(cp_file);

    CompPVecCI comp_e = pComps.end();
    for (CompPVecCI c = pComps.begin(); c != comp_e; ++c) (*c)->restore(cp_file);

    PatchPVecCI patch_e = pPatches.end();
    for (PatchPVecCI p = pPatches.begin(); p != patch_e; ++p) (*p)->restore(cp_file);

    DiffBoundaryPVecCI db_e = pDiffBoundaries.end();
    for (DiffBoundaryPVecCI db = pDiffBoundaries.begin(); db != db_e; ++db) {
        (*db)->restore(cp_file);
    }

    WmVolPVecCI wmv_e = pWmVols.end();
    for (WmVolPVecCI wmv = pWmVols.begin(); wmv != wmv_e; ++wmv)
    {
        if ((*wmv) != 0) {
            (*wmv)->restore(cp_file);
        }
    }

    TetPVecCI tet_e = pTets.end();
    for (TetPVecCI t = pTets.begin(); t != tet_e; ++t)
    {
        if ((*t) != 0) {
            (*t)->restore(cp_file);
        }
    }
    TriPVecCI tri_e = pTris.end();
    for (TriPVecCI t = pTris.begin(); t != tri_e; ++t)
    {
        if ((*t) != 0) {
            (*t)->restore(cp_file);
        }
    }

    KProcPVecCI e = pKProcs.end();
    for (KProcPVecCI i = pKProcs.begin(); i != e; ++i) (*i)->restore(cp_file);



    if (efflag()) {
        cp_file.read((char*)&pTemp, sizeof(double));
        cp_file.read((char*)&pEFDT, sizeof(double));
        pEField->restore(cp_file);
    }

    uint stored_entries = 0;
    cp_file.read((char*)&stored_entries, sizeof(uint));

    if (stored_entries != nEntries) {
        std::ostringstream os;
        os << "Unknown Restore Error!";
        throw steps::ArgErr(os.str());
    }

    // restore CR SSA
    cp_file.read((char*)&pSum, sizeof(double));
    cp_file.read((char*)&nSum, sizeof(double));
    cp_file.read((char*)&pA0, sizeof(double));

    uint n_ngroups;
    uint n_pgroups;

    cp_file.read((char*)&n_ngroups, sizeof(uint));
    cp_file.read((char*)&n_pgroups, sizeof(uint));

    nGroups.resize(n_ngroups);
    pGroups.resize(n_pgroups);

    for (uint i = 0; i < n_ngroups; i++) {
        unsigned capacity;
        unsigned size;
        double max;
        double sum;

        cp_file.read((char*)&capacity, sizeof(unsigned));
        cp_file.read((char*)&size, sizeof(unsigned));
        cp_file.read((char*)&max, sizeof(double));
        cp_file.read((char*)&sum, sizeof(double));

        nGroups[i] = new CRGroup(0, capacity);
        nGroups[i]->size = size;
        nGroups[i]->max = max;
        nGroups[i]->sum = sum;

        for (uint j = 0; j < size; j++) {
            uint idx;
            cp_file.read((char*)&idx, sizeof(uint));
            nGroups[i]->indices[j] = pKProcs[idx];
        }
    }

    for (uint i = 0; i < n_pgroups; i++) {
        unsigned capacity;
        unsigned size;
        double max;
        double sum;

        cp_file.read((char*)&capacity, sizeof(unsigned));
        cp_file.read((char*)&size, sizeof(unsigned));
        cp_file.read((char*)&max, sizeof(double));
        cp_file.read((char*)&sum, sizeof(double));

        pGroups[i] = new CRGroup(0, capacity);
        pGroups[i]->size = size;
        pGroups[i]->max = max;
        pGroups[i]->sum = sum;

        for (uint j = 0; j < size; j++) {
            uint idx;
            cp_file.read((char*)&idx, sizeof(uint));
            pGroups[i]->indices[j] = pKProcs[idx];
        }
    }

    cp_file.close();

    std::cout << "complete.\n";
}

////////////////////////////////////////////////////////////////////////////////


std::string stex::Tetexact::getSolverName(void) const
{
    return "tetexact";
}

////////////////////////////////////////////////////////////////////////////////

std::string stex::Tetexact::getSolverDesc(void) const
{
    return "SSA Composition and Rejection Exact Method in tetrahedral mesh";
}

////////////////////////////////////////////////////////////////////////////////

std::string stex::Tetexact::getSolverAuthors(void) const
{
    return "Stefan Wils, Iain Hepburn, Weiliang Chen";
}

////////////////////////////////////////////////////////////////////////////////

std::string stex::Tetexact::getSolverEmail(void) const
{
    return "steps.dev@gmail.com";
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setup(void)
{
    // Perform upcast.
    pMesh = dynamic_cast<steps::tetmesh::Tetmesh *>(geom());
    if (!pMesh)
        throw steps::ArgErr("Geometry description to steps::solver::Tetexact solver "
                "constructor is not a valid steps::tetmesh::Tetmesh object.");

    // First initialise the pTets, pTris vector, because
    // want tets and tris to maintain indexing from Geometry
    uint ntets = pMesh->countTets();
    uint ntris = pMesh->countTris();
    uint ncomps = pMesh->_countComps();

    pTets.assign(ntets, NULL);
    pTris.assign(ntris, NULL);
    pWmVols.assign(ncomps, NULL);

    // Now create the actual compartments.
    ssolver::CompDefPVecCI c_end = statedef()->endComp();
    for (ssolver::CompDefPVecCI c = statedef()->bgnComp(); c != c_end; ++c)
    {
        uint compdef_gidx = (*c)->gidx();
        uint comp_idx = _addComp(*c);
        assert(compdef_gidx == comp_idx);
    }
    // Create the actual patches.
    ssolver::PatchDefPVecCI p_end = statedef()->endPatch();
    for (ssolver::PatchDefPVecCI p = statedef()->bgnPatch(); p != p_end; ++p)
    {
        uint patchdef_gidx = (*p)->gidx();
        uint patch_idx = _addPatch(*p);
        assert(patchdef_gidx == patch_idx);
    }

    // Create the diffusion boundaries
    ssolver::DiffBoundarydefPVecCI db_end = statedef()->endDiffBoundary();
    for (ssolver::DiffBoundaryDefPVecCI db = statedef()->bgnDiffBoundary(); db != db_end; ++db)
    {
        uint diffboundary_gidx = (*db)->gidx();
        uint diffb_idx = _addDiffBoundary(*db);
        assert(diffboundary_gidx == diffb_idx);
    }

    uint npatches = pPatches.size();
    assert (pMesh->_countPatches() == npatches);
    for (uint p = 0; p < npatches; ++p)
    {
        // Add the tris for this patch
        // We have checked the indexing - p is the global index
        steps::wm::Patch * wmpatch = pMesh->_getPatch(p);

        // Perform upcast
        steps::tetmesh::TmPatch *tmpatch = dynamic_cast<steps::tetmesh::TmPatch*>(wmpatch);
        if (!tmpatch)
            throw steps::ArgErr("Well-mixed patches not supported in steps::solver::Tetexact solver.");
        steps::tetexact::Patch *localpatch = pPatches[p];

        for (uint tri: tmpatch->_getAllTriIndices()) 
        {
            assert (pMesh->getTriPatch(tri) == tmpatch);

            double area = pMesh->getTriArea(tri);

            // NB: Tri vertices may not be in consistent order, so use bar interface.
            const uint *tri_bars = pMesh->_getTriBars(tri);
            double l[3] = {0, 0, 0};

            for (uint j=0; j<3; ++j) {
                const uint *v = pMesh->_getBar(tri_bars[j]);
                l[j] = distance(pMesh->_getVertex(v[0]), pMesh->_getVertex(v[1]));
            }

            std::vector<int> tris = pMesh->getTriTriNeighb(tri, tmpatch);
            point3d baryc = pMesh->_getTriBarycenter(tri);

            double d[3] = {0, 0, 0};
            for (uint j = 0; j < 3; ++j) {
                if (tris[j]==-1) continue;
                d[j] = distance(baryc, pMesh->_getTriBarycenter(tris[j]));
            }

            const int *tri_tets = pMesh->_getTriTetNeighb(tri);
            _addTri(tri, localpatch, area, l[0], l[1], l[2], d[0], d[1], d[2], tri_tets[0], tri_tets[1], tris[0], tris[1], tris[2]);
        }
    }

    ncomps = pComps.size();
    assert (pMesh->_countComps() == ncomps);

    for (uint c = 0; c < ncomps; ++c)
    {
        // Now add the tets for this comp
         // We have checked the indexing- c is the global index
        steps::wm::Comp * wmcomp = pMesh->_getComp(c);

        // Perform upcast
        steps::tetmesh::TmComp *tmcomp = dynamic_cast<steps::tetmesh::TmComp*>(wmcomp);
        if (tmcomp) {
             steps::tetexact::Comp * localcomp = pComps[c];

             for (uint tet: tmcomp->_getAllTetIndices())
             {
                 assert (pMesh->getTetComp(tet) == tmcomp);

                 double vol = pMesh->getTetVol(tet);

                 const uint *tris = pMesh->_getTetTriNeighb(tet);

                 double a[4] = {0, 0, 0, 0};
                 for (uint j = 0; j < 4; ++j) {
                     a[j] = pMesh->getTriArea(tris[j]);
                 }

                 const int *tets = pMesh->_getTetTetNeighb(tet);
                 point3d baryc = pMesh->_getTetBarycenter(tet);

                 double d[4] = {0, 0, 0, 0};
                 for (uint j = 0; j < 4; ++j) {
                     if (tets[j]==-1) continue;
                     d[j] = distance(baryc, pMesh->_getTetBarycenter(tets[j]));
                 }

                 _addTet(tet, localcomp, vol, a[0], a[1], a[2], a[3], d[0], d[1], d[2], d[3],
                         tets[0], tets[1], tets[2], tets[3]);
             }
        }
        else
        {
            // This means that this compartment is a well-mixed compartment
            // It will behave like a tetrahedral-based compartment, but
            // contain only one 'voxel' that is connected to all surface
            // triangles and has the same volume as the whole compartment
            steps::tetexact::Comp * localcomp = pComps[c];
            uint cidx = c;
            _addWmVol(cidx, localcomp, localcomp->def()->vol());
            assert(pWmVols[c] != 0);

            // Now find all the triangles that reference this well-mixed volume
            // and set the inner or outer tetrahedron index accordingly.

            uint nopatches = wmcomp->_countOPatches();
            for (uint i = 0; i < nopatches; ++i)
            {
                steps::wm::Patch * op = wmcomp->_getOPatch(i);
                //     Comp may have no outer patch
                if (op == nullptr) continue;

                steps::tetmesh::TmPatch *comp_opatch = dynamic_cast<steps::tetmesh::TmPatch*>(op);
                if (!comp_opatch)
                    throw steps::ProgErr("Compartment outer patch is not a TmPatch.");

                for (uint tri: comp_opatch->_getAllTriIndices())
                {
                    pTris[tri]->setInnerTet(pWmVols[c]);
                    // Add triangle to WmVols' table of neighbouring triangles.
                    pWmVols[c]->setNextTri(pTris[tri]);
                }
            }

            uint nipatches = wmcomp->_countIPatches();
            for (uint i = 0; i < nipatches; ++i)
            {
                steps::wm::Patch * ip = wmcomp->_getIPatch(i);
                // Comp may not have an inner patch
                if (ip == nullptr) continue;

                steps::tetmesh::TmPatch *comp_ipatch = dynamic_cast<steps::tetmesh::TmPatch*>(ip);
                if (!comp_ipatch)
                    throw steps::ProgErr("Compartment inner patch is not a TmPatch.");

                for (uint tri: comp_ipatch->_getAllTriIndices())
                {
                    pTris[tri]->setOuterTet(pWmVols[c]);
                    // Add triangle to WmVols' table of neighbouring triangles.
                    pWmVols[c]->setNextTri(pTris[tri]);
                }
            }
        }
    }


    // All tets and tris that belong to some comp or patch have been created
    // locally- now we can connect them locally
    // NOTE: currently if a tetrahedron's neighbour belongs to a different
    // comp they do not talk to each other (see stex::Tet::setNextTet())
    //

    assert (ntets == pTets.size());
    // pTets member size of all tets in geometry, but may not be filled with
    // local tets if they have not been added to a compartment
    for (uint t = 0; t < ntets; ++t)
    {
        if (pTets[t] == 0) continue;

        for (uint j = 0; j < 4; ++j) {
            int tet = pTets[t]->tet(j);
            if (tet >= 0 && pTets[tet] != 0) pTets[t]->setNextTet(j, pTets[tet]);
        }
        // Not setting Tet triangles at this point- only want to set
        // for surface triangles
    }
    assert (ntris == pTris.size());

    for (uint t = 0; t < ntris; ++t)
    {
        // Looping over all possible tris, but only some have been added to a patch
        if (pTris[t] == 0) continue;

        for (uint j = 0; j < 3; ++j) {
            int tri = pTris[t]->tri(j);
            if (tri >= 0 && pTris[tri] != 0) pTris[t]->setNextTri(j, pTris[tri]);
        }

        // By convention, triangles in a patch should have an inner tetrahedron defined
        // (neighbouring tets 'flipped' if necessary in Tetmesh)
        // but not necessarily an outer tet
        // 17/3/10- actually this is not the case any more with well-mixed compartments
        //
        int tetinner = pTris[t]->tet(0);
        int tetouter = pTris[t]->tet(1);

        // Now inside and outside tetrahedrons may be normal tetrahedrons, which
        // means compartment is a mesh compartment, or wmvols describing a
        // well-mixed compartment with multiple triangle connections.


        if (tetinner >= 0)
        {
            // NEW FOR THIS VERSION: Tris store index of inner and outer tet (outer may not exist if on
            // surface) but tets may not belong to a compartment, even inner tets now
            // since they may be well-mixed compartments
            //
            if (pTets[tetinner] != 0)
            {
                // A triangle may already have an inner tet defined as a well-mixed
                // volume, but that should not be the case here:
                assert (pTris[t]->iTet() == 0);

                pTris[t]->setInnerTet(pTets[tetinner]);
                // Now add this triangle to inner tet's list of neighbours
                for (uint i=0; i <= 4; ++i)
                {
                    // include assert for debugging purposes and remove
                    // once this is tested
                    assert (i < 4);                                                        //////////
                    // check if there is already a neighbouring tet or tri
                    // In theory if there is a tri to add, the tet should
                    // have less than 4 neighbouring tets added because
                    // a neighbouring tet(s) is in a different compartment

                    // THIS IS NOT THE CASE ANYMORE: tets in different compartments can be neighbours
                    // so as to allow for diffusion boundaries

                    //     Also check tris because in some cases a surface tet
                    // may have more than 1 neighbouring tri
                    // NOTE: The order here will end up being different to the
                    // neighbour order at the Tetmesh level

                    // Now with diffusion boundaries, meaning tets can have neighbours that
                    // are in different comps, we must check the compartment
                    steps::tetexact::Tet * tet_in = pTets[tetinner];
                    if (tet_in->nextTet(i) != 0 && tet_in->compdef() == tet_in->nextTet(i)->compdef()) continue;

                    if (tet_in->nextTri(i) != 0) continue;
                    tet_in->setNextTri(i, pTris[t]);
                    break;
                }
            }
        }

        // DEBUG 18/03/09:
        // Now correct check, previously didn't allow for tet index == 0
        if (tetouter >= 0)
        {
            if (pTets[tetouter] != 0)
            {
                // A triangle may already have an inner tet defined as a well-mixed
                // volume, but that should not be the case here:
                assert (pTris[t]->oTet() == 0);

                pTris[t]->setOuterTet(pTets[tetouter]);
                // Add this triangle to outer tet's list of neighbours
                for (uint i=0; i <= 4; ++i)
                {
                    assert (i < 4);

                    // See above in that tets now store tets from different comps
                    steps::tetexact::Tet * tet_out = pTets[tetouter];

                    if (tet_out->nextTet(i) != 0 && tet_out->compdef() == tet_out->nextTet(i)->compdef()) continue;

                    if (tet_out->nextTri(i) != 0) continue;
                    tet_out->setNextTri(i, pTris[t]);
                    break;
                }
            }
        }
    }

    // Now loop over the diffusion boundaries:
    // 1) get all the triangles and get the two tetrahedrons
    // 2) figure out which direction is the direction for a tetrahedron
    // 3) add the tetrahedron and the direction to local object

    // This is here because we need all tets to have been assigned correctly
    // to compartments. Check every one and set the compA and compB for the db
    uint ndiffbnds = pDiffBoundaries.size();
    assert(ndiffbnds ==    pMesh->_countDiffBoundaries());

    for (uint db = 0; db < ndiffbnds; ++db)
    {
        steps::tetexact::DiffBoundary * localdiffb = pDiffBoundaries[db];

        uint compAidx = localdiffb->def()->compa();
        uint compBidx = localdiffb->def()->compb();
        steps::solver::Compdef * compAdef = statedef()->compdef(compAidx);
        steps::solver::Compdef * compBdef = statedef()->compdef(compBidx);

        for (uint dbtri: localdiffb->def()->tris())
        {
            const int *tri_tets = pMesh->_getTriTetNeighb(dbtri);

            int tetAidx = tri_tets[0];
            int tetBidx = tri_tets[1];
            assert(tetAidx >= 0 && tetBidx >= 0);

            steps::tetexact::Tet * tetA = _tet(tetAidx);
            steps::tetexact::Tet * tetB = _tet(tetBidx);
            assert(tetA != 0 && tetB != 0);

            steps::solver::Compdef *tetA_cdef = tetA->compdef();
            steps::solver::Compdef *tetB_cdef = tetB->compdef();
            assert(tetA_cdef != 0);
            assert(tetB_cdef != 0);

            if (tetA_cdef != compAdef)
            {
                assert(tetB_cdef == compAdef);
                assert(tetA_cdef == compBdef);
            }
            else
            {
                assert(tetB_cdef == compBdef);
                assert(tetA_cdef == compAdef);
            }

            // Ok, checks over, lets get down to business
            int direction_idx_a = -1;
            int direction_idx_b = -1;

            const uint *tetA_tris = pMesh->_getTetTriNeighb(tetAidx);
            const uint *tetB_tris = pMesh->_getTetTriNeighb(tetBidx);

            for (uint i = 0; i < 4; ++i)
            {
                if (tetA_tris[i] == dbtri)
                {
                    assert(direction_idx_a == -1);
                    direction_idx_a = i;
                }
                if (tetB_tris[i] == dbtri)
                {
                    assert(direction_idx_b == -1);
                    direction_idx_b = i;
                }
            }
            assert (direction_idx_a != -1);
            assert (direction_idx_b != -1);

            // Set the tetrahedron and direction to the Diff Boundary object
            localdiffb->setTetDirection(tetAidx, direction_idx_a);
            localdiffb->setTetDirection(tetBidx, direction_idx_b);
        }

        localdiffb->setComps(_comp(compAidx), _comp(compBidx));

        // Before the kprocs are set up ( in _setup) the tetrahedrons need to know the diffusion
        // boundary direction, so let's do it here  - the diff bounday has had all
        // tetrahedrons added

        // Might as well copy the vectors because we need to index through
        std::vector<uint> tets = localdiffb->getTets();
        std::vector<uint> tets_direction = localdiffb->getTetDirection();

        ntets = tets.size();
        assert (ntets <= pTets.size());
        assert (tets_direction.size() == ntets);

        for (uint t = 0; t < ntets; ++t)
            _tet(tets[t])->setDiffBndDirection(tets_direction[t]);
    }

    for (auto t: pTets)
        if (t) t->setupKProcs(this);

    for (auto wmv: pWmVols)
        if (wmv) wmv->setupKProcs(this);

    for (auto t: pTris)
        if (t) t->setupKProcs(this, efflag());

    // Resolve all dependencies
    for (auto t: pTets) {
        // DEBUG: vector holds all possible tetrahedrons,
        // but they have not necessarily been added to a compartment.
        if (!t) continue;
        for (auto k: t->kprocs()) k->setupDeps();
    }

    for (auto wmv: pWmVols) {
        // Vector allows for all compartments to be well-mixed, so
        // hold null-pointer for mesh compartments
        if (!wmv) continue;
        for (auto k: wmv->kprocs()) k->setupDeps();
    }

    for (auto t: pTris) {
        // DEBUG: vector holds all possible triangles, but
        // only patch triangles are filled
        if (!t) continue;
        for (auto k: t->kprocs()) k->setupDeps();
    }

    // Create EField structures if EField is to be calculated
    if (efflag() == true) _setupEField();

    nEntries = pKProcs.size();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setupEField(void)
{
    using steps::math::point3d;
    using namespace steps::solver::efield;

    //// Note to self: for now roughly following flow from original code in sim/controller.py.
    //// code for setting up a mesh was in func_tetmesh constructor and called functions.

    assert(efflag() == true);

    switch (pEFoption) {
    case EF_DEFAULT:
    case EF_DV_BDSYS:
        pEField = make_EField<dVSolverBanded>();
        break;
    default:
        throw steps::ArgErr("Unsupported E-Field solver.");
    }

    // Give temperature a default value of 20c
    pTemp = 293.15;

    uint nmembs = mesh()->_countMembs();

    if  (nmembs != 1)
    {
        std::ostringstream os;
        os << "Membrane potential solver currently supports only one ";
        os << "membrane description object.";
        throw steps::ArgErr(os.str());
    }

    steps::tetmesh::Memb * memb = mesh()->_getMemb(0);
    assert(memb != 0);

    // TODO: Decide what checks are needed for the membrane and implement them here

    pEFNTets = memb->countVolTets();
    pEFNTris = memb->countTris();
    pEFNVerts = memb->countVerts();

    pEFTets = new uint[neftets() * 4];
    assert(pEFTets != 0);

    // All the triangles we will count here are the real membrane triangles,
    // virtual triangles will not require a capacitance.
    pEFTris = new uint[neftris() * 3];
    assert(pEFTris != 0);

    pEFVerts = new double[nefverts() * 3];
    assert(pEFVerts != 0);

    uint nverts = mesh()->countVertices();
    uint ntris = mesh()->countTris();
    uint ntets= mesh()->countTets();

    pEFVert_GtoL = new int[nverts];
    for (uint i=0; i < nverts; ++i) pEFVert_GtoL[i] = -1;
    pEFTri_GtoL = new int[ntris];
    for (uint i=0; i< ntris; ++i) pEFTri_GtoL[i] = -1;
    pEFTet_GtoL = new int[ntets];
    for (uint i=0; i < ntets; ++i) pEFTet_GtoL[i] = -1;

    pEFTri_LtoG = new uint[neftris()];

    // Copy the data to local structures.

    std::vector<uint> membverts = memb->_getAllVertIndices();
    assert(membverts.size() == nefverts());
    for (uint efv = 0; efv < nefverts(); ++efv)
    {
        uint vertidx = membverts[efv];
        point3d verttemp = mesh()->_getVertex(vertidx);
        uint efv2 = efv*3;

        // CONVERTING TO MICRONS HERE. EFIELD OBJECT WILL NOT PERFORM THIS CONVERSION
        verttemp *= 1.0e6;
        pEFVerts[efv2] = verttemp[0];
        pEFVerts[efv2+1] = verttemp[1];
        pEFVerts[efv2+2] = verttemp[2];

        pEFVert_GtoL[vertidx] = efv;
    }

    std::vector<uint> membtets = memb->_getAllVolTetIndices();
    assert(membtets.size() == neftets());
    for (uint eft=0; eft < neftets(); ++eft)
    {
        uint tetidx = membtets[eft];
        const uint* tettemp = mesh()->_getTet(tetidx);
        uint eft2 = eft*4;

        // Convert to indices used by EField object
        int tv0 =  pEFVert_GtoL[tettemp[0]];
        int tv1 = pEFVert_GtoL[tettemp[1]];
        int tv2 = pEFVert_GtoL[tettemp[2]];
        int tv3 = pEFVert_GtoL[tettemp[3]];
        if  (tv0 ==-1 || tv1 == -1 || tv2 == -1 || tv3 == -1)
        {
            std::ostringstream os;
            os << "Failed to create EField structures.";
            throw steps::ProgErr(os.str());
        }

        pEFTets[eft2] = tv0;
        pEFTets[eft2+1] = tv1;
        pEFTets[eft2+2] = tv2;
        pEFTets[eft2+3] = tv3;

        pEFTet_GtoL[tetidx] = eft;
    }

    std::vector<uint> membtris = memb->_getAllTriIndices();
    assert(membtris.size() == neftris());

    pEFTris_vec.resize(neftris());

    for (uint eft = 0; eft < neftris(); ++eft)
    {
        uint triidx = membtris[eft];
        const uint* tritemp = mesh()->_getTri(triidx);
        uint eft2 = eft*3;

        // Convert to indices used by EField object
        int tv0 =  pEFVert_GtoL[tritemp[0]];
        int tv1 = pEFVert_GtoL[tritemp[1]];
        int tv2 = pEFVert_GtoL[tritemp[2]];
        if  (tv0 ==-1 || tv1 == -1 || tv2 == -1)
        {
            std::ostringstream os;
            os << "Failed to create EField structures.";
            throw steps::ProgErr(os.str());
        }

        pEFTris[eft2] = tv0;
        pEFTris[eft2+1] = tv1;
        pEFTris[eft2+2] = tv2;

        pEFTri_GtoL[triidx] = eft;
        pEFTri_LtoG[eft] = triidx;

        // This is added now for quicker iteration during run()
        // Extremely important for larger meshes, orders of magnitude times faster
        pEFTris_vec[eft] = pTris[triidx];
    }

    pEField->initMesh(nefverts(), pEFVerts, neftris(), pEFTris, neftets(), pEFTets, memb->_getOpt_method(), memb->_getOpt_file_name(), memb->_getSearch_percent());
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::saveMembOpt(std::string const & opt_file_name)
{
    if  (efflag() != true)
    {
        std::ostringstream os;
        os << "saveMembOpt method only available if running EField ";
        throw steps::ArgErr(os.str());
    }

    pEField->saveOptimal(opt_file_name);

}

////////////////////////////////////////////////////////////////////////////////

// 'Safe' global to local index translation methods that throw on error.

inline uint stex::Tetexact::specG2L_or_throw(Comp *comp, uint gidx) const {
    assert(gidx < statedef()->countSpecs());
    uint lidx = comp->def()->specG2L(gidx);

    if (lidx == steps::solver::LIDX_UNDEFINED)
        throw steps::ArgErr("species undefined in compartment");
    return lidx;
}

inline uint stex::Tetexact::specG2L_or_throw(Patch *patch, uint gidx) const {
    assert(gidx < statedef()->countSpecs());
    uint lidx = patch->def()->specG2L(gidx);

    if (lidx == steps::solver::LIDX_UNDEFINED)
        throw steps::ArgErr("species undefined in patch");
    return lidx;
}

#if 0
inline uint stex::Tetexact::specG2L_or_throw(Tet *tet, uint gidx) const {
    assert(gidx < statedef()->countSpecs());
    uint lidx = tet->compdef()->specG2L(gidx);

    if (lidx == steps::solver::LIDX_UNDEFINED)
        throw steps::ArgErr("species undefined in tetrahedron");
    return lidx;
}

inline uint stex::Tetexact::specG2L_or_throw(Tri *tri, uint gidx) const {
    assert(gidx < statedef()->countSpecs());
    uint lidx = tri->patchdef()->specG2L(gidx);

    if (lidx == steps::solver::LIDX_UNDEFINED)
        throw steps::ArgErr("species undefined in triangle");
    return lidx;
}
#endif

inline uint stex::Tetexact::reacG2L_or_throw(Comp *comp, uint gidx) const {
    assert(gidx < statedef()->countReacs());
    uint lidx = comp->def()->reacG2L(gidx);

    if (lidx == steps::solver::LIDX_UNDEFINED)
        throw steps::ArgErr("reaction undefined in compartment");
    return lidx;
}

inline uint stex::Tetexact::sreacG2L_or_throw(Patch *patch, uint gidx) const {
    assert(gidx < statedef()->countSReacs());
    uint lidx = patch->def()->sreacG2L(gidx);

    if (lidx == solver::LIDX_UNDEFINED)
        throw steps::ArgErr("surface reaction undefined in patch");
    return lidx;
}

inline uint stex::Tetexact::diffG2L_or_throw(Comp *comp, uint gidx) const {
    assert(gidx < statedef()->countDiffs());
    uint lidx = comp->def()->diffG2L(gidx);

    if (lidx == steps::solver::LIDX_UNDEFINED)
        throw steps::ArgErr("diffusion rule undefined in compartment");
    return lidx;
}

inline uint stex::Tetexact::sdiffG2L_or_throw(Patch *patch, uint gidx) const {
    assert(gidx < statedef()->countSurfDiffs());
    uint lidx = patch->def()->surfdiffG2L(gidx);

    if (lidx == steps::solver::LIDX_UNDEFINED)
        throw steps::ArgErr("diffusion rule undefined in patch");
    return lidx;
}

inline uint stex::Tetexact::vdepsreacG2L_or_throw(Patch *patch, uint gidx) const {
    assert(gidx < statedef()->countVDepSReacs());
    uint lidx = patch->def()->vdepsreacG2L(gidx);

    if (lidx == steps::solver::LIDX_UNDEFINED)
        throw steps::ArgErr("voltage-dependent surface reation undefined in patch");
    return lidx;
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::_addComp(steps::solver::Compdef * cdef)
{
    stex::Comp * comp = new Comp(cdef);
    assert(comp != 0);
    uint compidx = pComps.size();
    pComps.push_back(comp);
    pCompMap[cdef] = comp;
    return compidx;
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::_addPatch(steps::solver::Patchdef * pdef)
{
    stex::Patch * patch = new Patch(pdef);
    assert(patch != 0);
    uint patchidx = pPatches.size();
    pPatches.push_back(patch);
    return patchidx;
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::_addDiffBoundary(steps::solver::DiffBoundarydef * dbdef)
{
    stex::DiffBoundary * diffb = new DiffBoundary(dbdef);
    assert(diffb != 0);
    uint dbidx = pDiffBoundaries.size();
    pDiffBoundaries.push_back(diffb);
    return dbidx;
}


////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_addTet(uint tetidx,
                             steps::tetexact::Comp * comp, double vol,
                             double a1, double a2, double a3, double a4,
                             double d1, double d2, double d3, double d4,
                             int tet0, int tet1, int tet2, int tet3)
{
    steps::solver::Compdef * compdef  = comp->def();
    stex::Tet * localtet = new stex::Tet(tetidx, compdef, vol, a1, a2, a3, a4, d1, d2, d3, d4,
                                         tet0, tet1, tet2, tet3);
    assert(localtet != 0);
    assert(tetidx < pTets.size());
    assert(pTets[tetidx] == 0);
    pTets[tetidx] = localtet;
    comp->addTet(localtet);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_addWmVol(uint cidx, steps::tetexact::Comp * comp, double vol)
{
    steps::solver::Compdef * compdef  = comp->def();
    stex::WmVol * localtet = new stex::WmVol(cidx, compdef, vol);
    assert(localtet != 0);
    assert(cidx < pWmVols.size());
    pWmVols[cidx] = localtet;
    comp->addTet(localtet);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_addTri(uint triidx, steps::tetexact::Patch * patch, double area,
                             double l0, double l1, double l2, double d0, double d1, double d2,  int tinner, int touter, int tri0, int tri1, int tri2)
{
    steps::solver::Patchdef * patchdef = patch->def();
    stex::Tri * tri = new stex::Tri(triidx, patchdef, area, l0, l1, l2, d0, d1, d2,  tinner, touter, tri0, tri1, tri2);
    assert(tri != 0);
    assert (triidx < pTris.size());
    assert (pTris[triidx] == 0);
    pTris[triidx] = tri;
    patch->addTri(tri);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::reset(void)
{
    std::for_each(pComps.begin(), pComps.end(), std::mem_fun(&Comp::reset));
    std::for_each(pPatches.begin(), pPatches.end(), std::mem_fun(&Patch::reset));

    TetPVecCI tet_end = pTets.end();
    for (TetPVecCI tet = pTets.begin(); tet != tet_end; ++tet)
    {
        if (*tet == 0) continue;
        (*tet)->reset();
    }

    WmVolPVecCI wmvol_end = pWmVols.end();
    for (WmVolPVecCI wmvol = pWmVols.begin(); wmvol != wmvol_end; ++wmvol)
    {
        if ((*wmvol) == 0) continue;
        (*wmvol)->reset();
    }

    TriPVecCI tri_end = pTris.end();
    for (TriPVecCI t = pTris.begin(); t != tri_end; ++t)
    {
        if ((*t) == 0) continue;
        (*t)->reset();
    }

    uint ngroups = nGroups.size();
    for (uint i = 0; i < ngroups; i++) {
        free(nGroups[i]->indices);
        delete nGroups[i];
    }
    nGroups.clear();

    ngroups = pGroups.size();
    for (uint i = 0; i < ngroups; i++) {
        free(pGroups[i]->indices);
        delete pGroups[i];
    }
    pGroups.clear();

    pSum = 0.0;
    nSum = 0.0;
    pA0 = 0.0;

    statedef()->resetTime();
    statedef()->resetNSteps();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::run(double endtime)
{
    if (efflag() == false)
    {
        if (endtime < statedef()->time())
        {
            std::ostringstream os;
            os << "Endtime is before current simulation time";
            throw steps::ArgErr(os.str());
        }
        while (statedef()->time() < endtime)
        {
            stex::KProc * kp = _getNext();
            if (kp == 0) break;
            double a0 = getA0();
            if (a0 == 0.0) break;
            double dt = rng()->getExp(a0);
            if ((statedef()->time() + dt) > endtime) break;
            _executeStep(kp, dt);
        }
        statedef()->setTime(endtime);
    }
    else if (efflag() == true)
    {
        // Run the simulation, including the EField calculation.
        // This loop will assume that the SSA dt is sufficiently small so
        // that a number of SSA events execute between every EField calculation.
        // A warning message will be printed if the SSA dt is larger than the EField dt.
        // The EField dt is actually a MAXIMUM dt- the actual time for the EField
        // calculation will be exact with respect to the last event time in the
        // SSA before reaching the EField dt.
        while (statedef()->time() < endtime)
        {
            // The zero propensity
            double a0 = getA0();
            // We need a bool to check if the SSA contains no possible events. In
            // this rare case (continue to) execute the EField calculation to the endtime.
            bool ssa_on = true;
            double ssa_dt = 0.0;
            if (a0 != 0.0) ssa_dt = rng()->getExp(a0);
            else (ssa_on = false);
            // Set the actual efield dt. This value will take a maximum pEFDT.
            double ef_dt = 0.0;

            while (ssa_on && (ef_dt + ssa_dt) < pEFDT )
            {
                stex::KProc * kp = _getNext();
                if (kp == 0) break;
                _executeStep(kp, ssa_dt);
                ef_dt += ssa_dt;

                a0 = getA0();
                if (a0 != 0.0) ssa_dt = rng()->getExp(a0);
                else (ssa_on = false);

            }
            assert(ef_dt < pEFDT);

            // It's possible that ef_dt is zero here: ssa_dt is large, or has become large.
            // In that case print a warning but continue, running the EField simulation for EFDT
            if (ssa_on == false || ef_dt == 0.0)
            {
                std::ostringstream os;
                //os << "\nWARNING: SSA tau is larger than EField dt.";
                //std::cout << os << std::endl;
            }

            if (ef_dt == 0.0)
            {
                // This means that tau is larger than EField dt. We have no choice but to
                // increase the state time by pEFDT.
                ef_dt = pEFDT;
                statedef()->incTime(pEFDT);
            }

            // Now to perform the EField calculation. This means finding ohmic and GHK
            // currents from triangles during the ef_dt and applying these to the EField
            // object.

            TriPVecCI eftri_end = pEFTris_vec.end();
            uint tlidx = 0;
            double sttime = statedef()->time();
            #ifdef SERIAL_EFIELD_DEBUG
            CLOG(DEBUG, "steps_debug") << "Received currents:\n";
            #endif
            for (TriPVecCI eft = pEFTris_vec.begin(); eft != eftri_end; ++eft)
            {
                double v = pEField->getTriV(tlidx);
                double cur = (*eft)->computeI(v, ef_dt, sttime);
                #ifdef SERIAL_EFIELD_DEBUG
                if (cur != 0.0) {
                    CLOG(DEBUG, "steps_debug") << "lid: " << tlidx << " cur: " << cur;
                }
                #endif
                pEField->setTriI(tlidx, cur);
                tlidx++;
            }

            pEField->advance(ef_dt);
            
            #ifdef SERIAL_EFIELD_DEBUG
            std::vector<double> EFTrisV;
            
            tlidx = 0;
            for (uint t = 0; t < pEFNTris; t++)
            {
                EFTrisV.push_back(pEField->getTriV(tlidx));
                ++tlidx;
            }
            CLOG(DEBUG, "steps_debug") << "computed voltages: " << EFTrisV;
            #endif
            // TODO: Replace this with something that only resets voltage-dependent things
            _update();
        }
    }

    else assert(false);
}

////////////////////////////////////////////////////////////////////////

void stex::Tetexact::advance(double adv)
{
    if (adv < 0.0)
    {
        std::ostringstream os;
        os << "Time to advance cannot be negative";
        throw steps::ArgErr(os.str());
    }

    double endtime = statedef()->time() + adv;
    run(endtime);
}

////////////////////////////////////////////////////////////////////////

void stex::Tetexact::step(void)
{
    if (efflag() == true)
    {
        std::ostringstream os;
        os << "Method not available with EField calculation.";
        throw steps::ArgErr(os.str());
    }

    stex::KProc * kp = _getNext();
    if (kp == 0) return;
    double a0 = getA0();
    if (a0 == 0.0) return;
    double dt = rng()->getExp(a0);
    _executeStep(kp, dt);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::getTime(void) const
{
    return statedef()->time();
}

////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::getNSteps(void) const
{
    return statedef()->nsteps();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setTime(double time)
{
    statedef()->setTime(time);
}

////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setNSteps(uint nsteps)
{
    statedef()->setNSteps(nsteps);
}

////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setTemp(double t)
{
    if (efflag() == false)
    {
        std::ostringstream os;
        os << "\nWARNING: Temperature set in simulation without membrane ";
        os << "potential calculation will be ignored.\n";
        std::cout << os.str() << std::endl;
    }
    assert(t >= 0.0);
    pTemp = t;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompVol(uint cidx) const
{
    return _comp(cidx)->vol();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompCount(uint cidx, uint sidx) const
{
    stex::Comp *comp = _comp(cidx);
    uint slidx = specG2L_or_throw(comp, sidx);

    uint count = 0;
    for (auto &tet: comp->tets()) count += tet->pools()[slidx];

    return count;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompCount(uint cidx, uint sidx, double n)
{
    stex::Comp * comp = _comp(cidx);
    uint slidx = specG2L_or_throw(comp, sidx);

    // functions for distribution:
    auto set_count = [slidx](WmVol *tet, uint c) { tet->setCount(slidx, c); };
    auto inc_count = [slidx](WmVol *tet, int c) { tet->incCount(slidx, c); };
    auto weight = [](WmVol *tet) { return tet->vol(); };

    steps::util::distribute_quantity(n, comp->bgnTet(), comp->endTet(),
        weight, set_count, inc_count, *rng(), comp->def()->vol());

    for (auto &tet: comp->tets()) _updateSpec(tet, slidx);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompAmount(uint cidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getCompCount(cidx, sidx);
    return (count / smath::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompAmount(uint cidx, uint sidx, double a)
{
    // convert amount in mols to number of molecules
    double a2 = a * steps::math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setCompCount(cidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompConc(uint cidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    return _getCompCount(cidx, sidx) / (1.0e3 * _comp(cidx)->vol() * steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompConc(uint cidx, uint sidx, double c)
{
    assert(c >= 0.0);
    double count = c * (1.0e3 * _comp(cidx)->vol() * steps::math::AVOGADRO);
    // the following method does all the necessary argument checking
    _setCompCount(cidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getCompClamped(uint cidx, uint sidx) const
{
    stex::Comp * comp = _comp(cidx);
    uint lsidx = specG2L_or_throw(comp, sidx);

    for (auto &tet: comp->tets())
        if (!tet->clamped(lsidx)) return false;

    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompClamped(uint cidx, uint sidx, bool b)
{
    stex::Comp *comp = _comp(cidx);
    uint lsidx = specG2L_or_throw(comp, sidx);

    // Set the flag in def object, though this may not be necessary
    comp->def()->setClamped(lsidx, b);
    for (auto &tet: comp->tets()) tet->setClamped(lsidx, b);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompReacK(uint cidx, uint ridx) const
{
    stex::Comp *comp = _comp(cidx);
    uint lridx = reacG2L_or_throw(comp, ridx);

    // We're just returning the default value for this comp, individual
    // tets may have different Kcsts set individually
    return comp->def()->kcst(lridx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompReacK(uint cidx, uint ridx, double kf)
{
    assert(kf >= 0.0);
    stex::Comp *comp = _comp(cidx);
    uint lridx = reacG2L_or_throw(comp, ridx);

    // First set the default value for the comp
    comp->def()->setKcst(lridx, kf);

    // Now update all tetrahedra in this comp
    for (auto &tet: comp->tets()) tet->reac(lridx)->setKcst(kf);

    // Rates have changed
    _update();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getCompReacActive(uint cidx, uint ridx) const
{
    stex::Comp *comp = _comp(cidx);
    uint lridx = reacG2L_or_throw(comp, ridx);

    for (auto &tet: comp->tets())
        if (tet->reac(lridx)->inactive()) return false;

    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompReacActive(uint cidx, uint ridx, bool a)
{
    stex::Comp *comp = _comp(cidx);
    uint lridx = reacG2L_or_throw(comp, ridx);

    // Set the default value for the comp, though this is not entirely
    // necessary
    comp->def()->setActive(lridx, a);

    for (auto &tet: comp->tets()) tet->reac(lridx)->setActive(a);

    // It's cheaper to just recompute everything.
    _update();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompDiffD(uint cidx, uint didx) const
{
    stex::Comp *comp = _comp(cidx);
    uint ldidx = diffG2L_or_throw(comp, didx);

    // We're just returning the default value for this comp, individual
    // tets may have different Dcsts set individually
    return comp->def()->dcst(ldidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompDiffD(uint cidx, uint didx, double dk)
{
    assert(dk >= 0.0);
    stex::Comp *comp = _comp(cidx);
    uint ldidx = diffG2L_or_throw(comp, didx);

    // First set the default value for the comp
    comp->def()->setDcst(ldidx, dk);

    // Now update all tets in this comp
    for (auto &wmvol: comp->tets()) {
        stex::Tet *tet = dynamic_cast<stex::Tet *>(wmvol);
        if (!tet)
            throw steps::ArgErr("cannot change diffusion constant in well-mixed compartment");

        tet->diff(ldidx)->setDcst(dk);
    }

    // Rates have changed
    _update();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getCompDiffActive(uint cidx, uint didx) const
{
    stex::Comp *comp = _comp(cidx);
    uint ldidx = diffG2L_or_throw(comp, didx);

    for (auto &wmvol: comp->tets()) {
        stex::Tet *tet = dynamic_cast<stex::Tet *>(wmvol);
        if (!tet)
            throw steps::ArgErr("diffusion activation not defined in well-mixed compartment");

        if (tet->diff(ldidx)->inactive()) return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompDiffActive(uint cidx, uint didx, bool act)
{
    stex::Comp *comp = _comp(cidx);
    uint ldidx = diffG2L_or_throw(comp, didx);

    for (auto &wmvol: comp->tets()) {
        stex::Tet *tet = dynamic_cast<stex::Tet *>(wmvol);
        if (!tet)
            throw steps::ArgErr("diffusion activation not defined in well-mixed compartment");

        tet->diff(ldidx)->setActive(act);
    }

    // It's cheaper to just recompute everything.
    _update();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchArea(uint pidx) const
{
    return _patch(pidx)->area();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchCount(uint pidx, uint sidx) const
{
    stex::Patch *patch = _patch(pidx);
    uint slidx = specG2L_or_throw(patch, sidx);

    uint count = 0;
    for (auto &tri: patch->tris()) count += tri->pools()[slidx];
    return count;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchCount(uint pidx, uint sidx, double n)
{
////// 
    stex::Patch *patch = _patch(pidx);
    uint slidx = specG2L_or_throw(patch, sidx);

    // functions for distribution:
    auto set_count = [slidx](Tri *tri, uint c) { tri->setCount(slidx, c); };
    auto inc_count = [slidx](Tri *tri, int c) { tri->incCount(slidx, c); };
    auto weight = [](Tri *tri) { return tri->area(); };

    steps::util::distribute_quantity(n, patch->bgnTri(), patch->endTri(),
        weight, set_count, inc_count, *rng(), patch->def()->area());

    for (auto &tri: patch->tris()) _updateSpec(tri, slidx);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchAmount(uint pidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getPatchCount(pidx, sidx);
    return (count / steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchAmount(uint pidx, uint sidx, double a)
{
    assert(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * steps::math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setPatchCount(pidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getPatchClamped(uint pidx, uint sidx) const
{
    stex::Patch * patch = _patch(pidx);
    uint lsidx = specG2L_or_throw(patch, sidx);

    for (auto &tri: patch->tris()) {
        if (!tri->clamped(lsidx)) return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchClamped(uint pidx, uint sidx, bool buf)
{
    stex::Patch *patch = _patch(pidx);
    uint lsidx = specG2L_or_throw(patch, sidx);

    // Set the flag in def object for consistency, though this is not
    // entirely necessary
    patch->def()->setClamped(lsidx, buf);

    for (auto &tri: patch->tris()) tri->setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchSReacK(uint pidx, uint ridx) const
{
    stex::Patch *patch = _patch(pidx);
    uint lsridx = sreacG2L_or_throw(patch, ridx);

    // We're just returning the default value for this patch, individual
    // triangles may have different Kcsts set
    return patch->def()->kcst(lsridx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchSReacK(uint pidx, uint ridx, double kf)
{
    assert(kf >= 0.0);
    stex::Patch *patch = _patch(pidx);
    uint lsridx = sreacG2L_or_throw(patch, ridx);

    // First set the default values for this patch
    patch->def()->setKcst(lsridx, kf);

    // Now update all triangles in this patch
    for (auto &tri: patch->tris()) tri->sreac(lsridx)->setKcst(kf);

    // Rates have changed
    _update();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getPatchSReacActive(uint pidx, uint ridx) const
{
    stex::Patch *patch = _patch(pidx);
    uint lsridx = sreacG2L_or_throw(patch, ridx);

    for (auto &tri: patch->tris()) {
        if (tri->sreac(lsridx)->inactive()) return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setDiffBoundaryDiffusionActive(uint dbidx, uint sidx, bool act)
{
    // Need to do two things:
    // 1) check if the species is defined in both compartments conencted
    // by the diffusion boundary
    // 2) loop over all tetrahedrons around the diff boundary and then the
    // diffusion rules and activate diffusion if the diffusion rule
    // relates to this species

    stex::DiffBoundary * diffb = _diffboundary(dbidx);
    stex::Comp * compA = diffb->compA();
    stex::Comp * compB = diffb->compB();

    uint lsidxA = specG2L_or_throw(compA, sidx);
    uint lsidxB = specG2L_or_throw(compB, sidx);

    std::vector<uint> bdtets = diffb->getTets();
    std::vector<uint> bdtetsdir = diffb->getTetDirection();

    // Have to use indices rather than iterator because need access to the
    // tet direction
    uint ntets = bdtets.size();

    for (uint bdt = 0; bdt != ntets; ++bdt)
    {
        stex::Tet * tet = _tet(bdtets[bdt]);
        uint direction = bdtetsdir[bdt];
        assert(direction >= 0 and direction < 4);

        // Each diff kproc then has access to the species through it's defined parent
        uint ndiffs = tet->compdef()->countDiffs();
        for (uint d = 0; d != ndiffs; ++d)
        {
            stex::Diff * diff = tet->diff(d);
            // sidx is the global species index; so is the lig() return from diffdef
            uint specgidx = diff->def()->lig();
            if (specgidx == sidx)
            {
                diff->setDiffBndActive(direction, act);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getDiffBoundaryDiffusionActive(uint dbidx, uint sidx) const
{
    stex::DiffBoundary * diffb = _diffboundary(dbidx);
    stex::Comp * compA = diffb->compA();
    stex::Comp * compB = diffb->compB();

    uint lsidxA = specG2L_or_throw(compA,sidx);
    uint lsidxB = specG2L_or_throw(compB,sidx);

    std::vector<uint> bdtets = diffb->getTets();
    std::vector<uint> bdtetsdir = diffb->getTetDirection();

    // Have to use indices rather than iterator because need access to the
    // tet direction

    uint ntets = bdtets.size();

    for (uint bdt = 0; bdt != ntets; ++bdt)
    {
        stex::Tet * tet = _tet(bdtets[bdt]);
        uint direction = bdtetsdir[bdt];
        assert(direction >= 0 and direction < 4);

        // Each diff kproc then has access to the species through it's defined parent
        uint ndiffs = tet->compdef()->countDiffs();
        for (uint d = 0; d != ndiffs; ++d)
        {
            stex::Diff * diff = tet->diff(d);
            // sidx is the global species index; so is the lig() return from diffdef
            uint specgidx = diff->def()->lig();
            if (specgidx == sidx)
            {
                // Just need to check the first one
                if (diff->getDiffBndActive(direction)) return true;
                else return false;
            }
        }
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setDiffBoundaryDcst(uint dbidx, uint sidx, double dcst, uint direction_comp)
{
    stex::DiffBoundary * diffb = _diffboundary(dbidx);
    stex::Comp * compA = diffb->compA();
    stex::Comp * compB = diffb->compB();
    
    uint lsidxA = specG2L_or_throw(compA,sidx);
    uint lsidxB = specG2L_or_throw(compB,sidx);

    steps::solver::Compdef * dirc_compdef = NULL;
    if (direction_comp != std::numeric_limits<uint>::max()) {
        dirc_compdef = _comp(direction_comp)->def();
    }
    
    std::vector<uint> bdtets = diffb->getTets();
    std::vector<uint> bdtetsdir = diffb->getTetDirection();
    
    uint ntets = bdtets.size();
    
    for (uint bdt = 0; bdt != ntets; ++bdt)
    {
        stex::Tet * tet = _tet(bdtets[bdt]);
        // if tet compdef equals to dirc_compdef,
        //it is the desination tet so diff should not be changed
        // NULL (bidirection) and source tet are both different
        // fromdirc_compdef
        if (dirc_compdef == tet->compdef()) {
            continue;
        }
        uint direction = bdtetsdir[bdt];
        assert(direction >= 0 and direction < 4);
        
        // Each diff kproc then has access to the species through it's defined parent
        uint ndiffs = tet->compdef()->countDiffs();
        for (uint d = 0; d != ndiffs; ++d)
        {
            stex::Diff * diff = tet->diff(d);
            // sidx is the global species index; so is the lig() return from diffdef
            uint specgidx = diff->def()->lig();
            if (specgidx == sidx)
            {
                #ifdef DIRECTIONAL_DCST_DEBUG
                CLOG(DEBUG, "steps_debug") << "direction: " << direction << " dcst: " << dcst << "\n";
                #endif
                
                diff->setDiffBndActive(direction, true);
                diff->setDirectionDcst(direction, dcst);
                _updateElement(diff);
            }
        }
    }
    
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchSReacActive(uint pidx, uint ridx, bool a)
{
    stex::Patch * patch = _patch(pidx);
    uint lsridx = sreacG2L_or_throw(patch, ridx);

    // First set the flags in def object for consistency, though this is
    // not entirely necessary for this solver
    patch->def()->setActive(lsridx, a);

    for (auto &tri: patch->tris()) tri->sreac(lsridx)->setActive(a);

    // It's cheaper to just recompute everything.
    _update();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getPatchVDepSReacActive(uint pidx, uint vsridx) const
{
    stex::Patch *patch = _patch(pidx);
    uint lvsridx = vdepsreacG2L_or_throw(patch, vsridx);

    for (auto &tri: patch->tris()) {
        if (tri->vdepsreac(lvsridx)->inactive()) return false;
    }

    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchVDepSReacActive(uint pidx, uint vsridx, bool a)
{
    stex::Patch *patch = _patch(pidx);
    assert(patch != 0);
    uint lvsridx = vdepsreacG2L_or_throw(patch, vsridx);


    // Not necessary and not possible to set the flags in def object
    for (auto &tri: patch->tris()) tri->vdepsreac(lvsridx)->setActive(a);

    // It's cheaper to just recompute everything.
    _update();
}
////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::addKProc(steps::tetexact::KProc * kp)
{
    assert (kp != 0);

    SchedIDX nidx = pKProcs.size();
    pKProcs.push_back(kp);
    kp->setSchedIDX(nidx);
}

////////////////////////////////////////////////////////////////////////////////
/*
void stex::Tetexact::_build(void)
{
    assert (pBuilt == false);

    pBuilt = true;
}
*/
////////////////////////////////////////////////////////////////////////////////

steps::tetexact::KProc * stex::Tetexact::_getNext(void) const
{

    assert(pA0 >= 0.0);
    // Quick check to see whether nothing is there.
    if (pA0 == 0.0) return NULL;

    double selector = pA0 * rng()->getUnfII();

    double partial_sum = 0.0;

    uint n_neg_groups = nGroups.size();
    uint n_pos_groups = pGroups.size();

    for (uint i = 0; i < n_neg_groups; i++) {
        CRGroup* group = nGroups[i];
        if (group->size == 0) continue;

        if (selector > partial_sum + group->sum) {
            partial_sum += group->sum;

            continue;
        }

        double g_max = group->max;
        double random_rate = g_max * rng()->getUnfII();;
        uint group_size = group->size;
        uint random_pos = rng()->get() % group_size;
        KProc* random_kp = group->indices[random_pos];

        while (random_kp->crData.rate <= random_rate) {
            random_rate = g_max * rng()->getUnfII();
            random_pos = rng()->get() % group_size;
            random_kp = group->indices[random_pos];
        }

        return random_kp;
    }


    for (uint i = 0; i < n_pos_groups; i++) {
        CRGroup* group = pGroups[i];
        if (group->size == 0) continue;

        if (selector > partial_sum + group->sum) {
            partial_sum += group->sum;
            continue;
        }

        double g_max = group->max;
        double random_rate = g_max * rng()->getUnfII();;
        uint group_size = group->size;
        uint random_pos = rng()->get() % group_size;
        KProc* random_kp = group->indices[random_pos];


        while (random_kp->crData.rate <= random_rate) {
            random_rate = g_max * rng()->getUnfII();
            random_pos = rng()->get() % group_size;
            random_kp = group->indices[random_pos];
        }

        return random_kp;
    }

    // Precision rounding error force clean up
    // Force the search in the last non-empty group
    for (int i = n_pos_groups - 1; i >= 0; i--) {
        CRGroup* group = pGroups[i];
        if (group->size == 0) continue;

        double g_max = group->max;
        double random_rate = g_max * rng()->getUnfII();;
        uint group_size = group->size;
        uint random_pos = rng()->get() % group_size;
        KProc* random_kp = group->indices[random_pos];


        while (random_kp->crData.rate <= random_rate) {
            random_rate = g_max * rng()->getUnfII();
            random_pos = rng()->get() % group_size;
            random_kp = group->indices[random_pos];

        }

        return random_kp;
    }

    for (int i = n_neg_groups - 1; i >= 0; i--) {
        CRGroup* group = nGroups[i];
        if (group->size == 0) continue;

        double g_max = group->max;
        double random_rate = g_max * rng()->getUnfII();;
        uint group_size = group->size;
        uint random_pos = rng()->get() % group_size;
        KProc* random_kp = group->indices[random_pos];



        while (random_kp->crData.rate <= random_rate) {
            random_rate = g_max * rng()->getUnfII();
            random_pos = rng()->get() % group_size;
            random_kp = group->indices[random_pos];

        }

        return random_kp;
    }

    // Precision rounding error force clean up - Complete

    std::cerr << "Cannot find any suitable entry.\n";
    std::cerr << "A0: " << std::setprecision (15) << pA0 << "\n";
    std::cerr << "Selector: " << std::setprecision (15) << selector << "\n";
    std::cerr << "Current Partial Sum: " << std::setprecision (15) << partial_sum << "\n";

    std::cerr << "Distribution of group sums\n";
    std::cerr << "Negative groups\n";

    for (uint i = 0; i < n_neg_groups; i++) {
        std::cerr << i << ": " << std::setprecision (15) << nGroups[i]->sum << "\n";
    }
    std::cerr << "Positive groups\n";
    for (uint i = 0; i < n_pos_groups; i++) {
        std::cerr << i << ": " << std::setprecision (15) << pGroups[i]->sum << "\n";
    }

    throw;
}

////////////////////////////////////////////////////////////////////////////////
/*
void stex::Tetexact::_reset(void)
{

}
*/
////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_executeStep(steps::tetexact::KProc * kp, double dt)
{
    std::vector<KProc*> const & upd = kp->apply(rng(), dt, statedef()->time());
    _update(upd.begin(), upd.end());
    statedef()->incTime(dt);
    statedef()->incNSteps(1);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_updateSpec(steps::tetexact::WmVol * tet, uint spec_lidx)
{
    std::set<KProc*> updset;

    // Loop over tet.
    for (auto &kproc: tet->kprocs()) updset.insert(kproc);

    for (auto &tri: tet->nexttris()) {
        if (!tri) continue;
        for (auto &kproc: tri->kprocs()) updset.insert(kproc);
    }

    // Send the list of kprocs that need to be updated to the schedule.
    _update(updset.begin(), updset.end());
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_updateSpec(steps::tetexact::Tri * tri, uint spec_lidx)
{
    _update(tri->kprocBegin(), tri->kprocEnd());
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompReacH(uint cidx, uint ridx) const
{
    Comp *comp = _comp(cidx);
    uint lridx = reacG2L_or_throw(comp, ridx);

    double h = 0.0;
    for (auto &tet: comp->tets()) h += tet->reac(lridx)->h();
    return h;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompReacC(uint cidx, uint ridx) const
{
    Comp *comp = _comp(cidx);
    uint lridx = reacG2L_or_throw(comp, ridx);

    if (comp->tets().empty()) return 0.0;

    double c = 0.0;
    double v = 0.0;
    for (auto &tet: comp->tets()) {
        double v2 = tet->vol();
        c += tet->reac(lridx)->c() * v2;
        v += v2;
    }
    assert(v > 0.0);
    return c/v;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompReacA(uint cidx, uint ridx) const
{
    Comp *comp = _comp(cidx);
    uint lridx = reacG2L_or_throw(comp, ridx);

    double a = 0.0;
    for (auto &tet: comp->tets()) a += tet->reac(lridx)->rate();
    return a;
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::_getCompReacExtent(uint cidx, uint ridx) const
{
    stex::Comp *comp = _comp(cidx);
    uint lridx = reacG2L_or_throw(comp, ridx);

    uint x = 0.0;
    for (auto &tet: comp->tets()) x += tet->reac(lridx)->getExtent();
    return x;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_resetCompReacExtent(uint cidx, uint ridx)
{
    Comp *comp = _comp(cidx);
    uint lridx = reacG2L_or_throw(comp, ridx);

    // The 'local' Comp object has same index as solver::Compdef object
    for (auto &tet: comp->tets()) tet->reac(lridx)->resetExtent();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchSReacH(uint pidx, uint ridx) const
{
    Patch *patch = _patch(pidx);
    uint lsridx = sreacG2L_or_throw(patch, ridx);

    double h = 0.0;
    for (auto &tri: patch->tris()) h += tri->sreac(lsridx)->h();
    return h;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchSReacC(uint pidx, uint ridx) const
{
    Patch *patch = _patch(pidx);
    uint lsridx = sreacG2L_or_throw(patch, ridx);

    if (patch->tris().empty()) return 0.0;

    double c = 0.0;
    double a = 0.0;
    for (auto &tri: patch->tris()) {
        double a2 = tri->area();
        c += tri->sreac(lsridx)->c() * a2;
        a += a2;
    }
    assert(a > 0.0);
    return c/a;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchSReacA(uint pidx, uint ridx) const
{
    Patch *patch = _patch(pidx);
    uint lsridx = sreacG2L_or_throw(patch, ridx);

    double a = 0.0;
    for (auto &tri: patch->tris()) a += tri->sreac(lsridx)->rate();
    return a;
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::_getPatchSReacExtent(uint pidx, uint ridx) const
{
    Patch *patch = _patch(pidx);
    uint lsridx = sreacG2L_or_throw(patch, ridx);

    double x = 0.0;
    for (auto &tri: patch->tris()) x += tri->sreac(lsridx)->getExtent();
    return x;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_resetPatchSReacExtent(uint pidx, uint ridx)
{
    Patch *patch = _patch(pidx);
    uint lsridx = sreacG2L_or_throw(patch, ridx);

    for (auto &tri: patch->tris()) tri->sreac(lsridx)->resetExtent();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetVol(uint tidx) const
{
    assert (tidx < pTets.size());
    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
        throw steps::ArgErr(os.str());
    }
    return pTets[tidx]->vol();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetVol(uint tidx, double vol)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTetSpecDefined(uint tidx, uint sidx) const
{
    assert (tidx < pTets.size());
    assert (sidx < statedef()->countSpecs());

    if (pTets[tidx] == 0) return false;

    stex::Tet * tet = pTets[tidx];
    uint lsidx = tet->compdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return false;
    else return true;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetCount(uint tidx, uint sidx) const
{
    assert (tidx < pTets.size());
    assert (sidx < statedef()->countSpecs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];
    uint lsidx = tet->compdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    return tet->pools()[lsidx];
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetCount(uint tidx, uint sidx, double n)
{
    assert (tidx < pTets.size());
    assert (sidx < statedef()->countSpecs());
    assert (n >= 0.0);

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }
    if (n > std::numeric_limits<unsigned int>::max( ))
    {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << std::numeric_limits<unsigned int>::max( ) << ").\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];

    uint lsidx = tet->compdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    double n_int = std::floor(n);
    double n_frc = n - n_int;
    uint c = static_cast<uint>(n_int);
    if (n_frc > 0.0)
    {
        double rand01 = rng()->getUnfIE();
        if (rand01 < n_frc) c++;
    }

    // Tet object updates def level Comp object counts
    tet->setCount(lsidx, c);
    _updateSpec(tet, lsidx);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetAmount(uint tidx, uint sidx) const
{
    // following method does all necessary argument checking
    double count = _getTetCount(tidx, sidx);
    return count/steps::math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetAmount(uint tidx, uint sidx, double m)
{
    // convert amount in mols to number of molecules
    double m2 = m * steps::math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setTetCount(tidx, sidx, m2);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetConc(uint tidx, uint sidx) const
{
    // following method does all necessary argument checking
    double count = _getTetCount(tidx, sidx);
    stex::Tet * tet = pTets[tidx];
    double vol = tet->vol();
    return (count/(1.0e3 * vol * steps::math::AVOGADRO));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetConc(uint tidx, uint sidx, double c)
{
    assert (c >= 0.0);
    assert (tidx < pTets.size());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];
    double count = c * (1.0e3 * tet->vol() * steps::math::AVOGADRO);
    // the following method does all the necessary argument checking
    _setTetCount(tidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTetClamped(uint tidx, uint sidx) const
{
    assert (tidx < pTets.size());
    assert (sidx < statedef()->countSpecs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];

    uint lsidx = tet->compdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    return tet->clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetClamped(uint tidx, uint sidx, bool buf)
{
    assert (tidx < pTets.size());
    assert (sidx < statedef()->countSpecs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];

    uint lsidx = tet->compdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    tet->setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetReacK(uint tidx, uint ridx) const
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];

    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    return (tet->reac(lridx)->kcst());
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetReacK(uint tidx, uint ridx, double kf)
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());
    assert (kf >= 0.0);

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];

    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "\nReaction undefined in tetrahedron.";
        throw steps::ArgErr(os.str());
    }

    tet->reac(lridx)->setKcst(kf);

    _updateElement(tet->reac(lridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTetReacActive(uint tidx, uint ridx) const
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];

    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    if (tet->reac(lridx)->inactive() == true) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetReacActive(uint tidx, uint ridx, bool act)
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];

    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    tet->reac(lridx)->setActive(act);

    _updateElement(tet->reac(lridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetDiffD(uint tidx, uint didx, uint direction_tet) const
{
    assert (tidx < pTets.size());
    assert (didx < statedef()->countDiffs());
    
    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }
    
    stex::Tet * tet = pTets[tidx];
    
    uint ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }
    
    if (direction_tet == std::numeric_limits<uint>::max()) {
        return tet->diff(ldidx)->dcst();
    }
    else {
        int direction = tet->getTetDirection(direction_tet);
        if (direction == -1) {
            std::ostringstream os;
            os << "Tetrahedron " << direction_tet << " is not a neighbor of tetrahedron " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }

        return tet->diff(ldidx)->dcst(direction);
    }
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetDiffD(uint tidx, uint didx, double dk, uint direction_tet)
{
    #ifdef DIRECTIONAL_DCST_DEBUG
    CLOG_IF(direction_tet !=  std::numeric_limits<uint>::max(), DEBUG, "steps_debug") << "tidx: " << tidx << " didx: " << didx << " dk: " << dk << " direction tet: " << direction_tet << "\n";
    
    CLOG_IF(direction_tet ==  std::numeric_limits<uint>::max(), DEBUG, "steps_debug") << "tidx: " << tidx << " didx: " << didx << " dk: " << dk << " for all directions.\n";
    #endif
    
    assert (tidx < pTets.size());
    assert (didx < statedef()->countDiffs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];

    uint ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    if (direction_tet == std::numeric_limits<uint>::max()) {
        tet->diff(ldidx)->setDcst(dk);
    }
    else {
        int direction = tet->getTetDirection(direction_tet);
        if (direction == -1) {
            std::ostringstream os;
            os << "Tetrahedron " << direction_tet << " is not a neighbor of tetrahedron " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        #ifdef DIRECTIONAL_DCST_DEBUG
        CLOG(DEBUG, "steps_debug") << "use tet " << direction_tet << " to set direction " << direction << ".\n";
        #endif
        
        tet->diff(ldidx)->setDirectionDcst(direction, dk);
    }
    _updateElement(tet->diff(ldidx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTetDiffActive(uint tidx, uint didx) const
{
    assert (tidx < pTets.size());
    assert (didx < statedef()->countDiffs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];

    uint ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    if (tet->diff(ldidx)->inactive() == true) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetDiffActive(uint tidx, uint didx, bool act)
{
    assert (tidx < pTets.size());
    assert (didx < statedef()->countDiffs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];

    uint ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    tet->diff(ldidx)->setActive(act);

    _updateElement(tet->diff(ldidx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetReacH(uint tidx, uint ridx) const
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];

    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    return tet->reac(lridx)->h();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetReacC(uint tidx, uint ridx) const
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];

    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    return tet->reac(lridx)->c();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetReacA(uint tidx, uint ridx) const
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];

    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    return tet->reac(lridx)->rate();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetDiffA(uint tidx, uint didx) const
{
    assert (tidx < pTets.size());
    assert (didx < statedef()->countDiffs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tet * tet = pTets[tidx];

    uint ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    return tet->diff(ldidx)->rate();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriArea(uint tidx) const
{
    assert (tidx < pTris.size());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.";
        throw steps::ArgErr(os.str());
    }

    return pTris[tidx]->area();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriArea(uint tidx, double area)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTriSpecDefined(uint tidx, uint sidx) const
{
    assert (tidx < pTris.size());
    assert (sidx < statedef()->countSpecs());

    if (pTris[tidx] == 0) return false;

    stex::Tri * tri = pTris[tidx];
    uint lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return false;
    else return true;
}

////////////////////////////////////////////////////////////////////////////////


double stex::Tetexact::_getTriCount(uint tidx, uint sidx) const
{
    assert (tidx < pTris.size());
    assert (sidx < statedef()->countSpecs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];
    uint lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    return tri->pools()[lsidx];
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriCount(uint tidx, uint sidx, double n)
{
    assert (tidx < pTris.size());
    assert (sidx < statedef()->countSpecs());
    assert (n >= 0.0);

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }
    if (n > std::numeric_limits<unsigned int>::max( ))
    {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << std::numeric_limits<unsigned int>::max( ) << ").\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];
    uint lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    double n_int = std::floor(n);
    double n_frc = n - n_int;
    uint c = static_cast<uint>(n_int);
    if (n_frc > 0.0)
    {
        double rand01 = rng()->getUnfIE();
        if (rand01 < n_frc) c++;
    }

    // Tri object updates counts in def level Comp object
    tri->setCount(lsidx, c);
    _updateSpec(tri, lsidx);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriAmount(uint tidx, uint sidx) const
{
    // following method does all necessary argument checking
    return _getTriCount(tidx, sidx) / steps::math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriAmount(uint tidx, uint sidx, double m)
{
    // the following method does all the necessary argument checking
    _setTriCount(tidx, sidx, m * steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////


bool stex::Tetexact::_getTriClamped(uint tidx, uint sidx) const
{
    assert (tidx < pTris.size());
    assert (sidx < statedef()->countSpecs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    uint lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    return tri->clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriClamped(uint tidx, uint sidx, bool buf)
{
    assert (tidx < pTris.size());
    assert (sidx < statedef()->countSpecs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    uint lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    tri->setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriSReacK(uint tidx, uint ridx) const
{
    assert (tidx < pTris.size());
    assert (ridx < statedef()->countSReacs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    uint lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    return (tri->sreac(lsridx)->kcst());
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriSReacK(uint tidx, uint ridx, double kf)
{
    assert (tidx < pTris.size());
    assert (ridx < statedef()->countSReacs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    uint lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    tri->sreac(lsridx)->setKcst(kf);

    _updateElement(tri->sreac(lsridx));
    _updateSum();

}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTriSReacActive(uint tidx, uint ridx) const
{
    assert (tidx < pTris.size());
    assert (ridx < statedef()->countSReacs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    uint lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    if (tri->sreac(lsridx)->inactive() == true) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriSReacActive(uint tidx, uint ridx, bool act)
{
    assert (tidx < pTris.size());
    assert (ridx < statedef()->countSReacs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    uint lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    tri->sreac(lsridx)->setActive(act);

    _updateElement(tri->sreac(lsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriDiffD(uint tidx, uint didx, uint direction_tri) const
{

    assert (tidx < pTris.size());
    assert (didx < statedef()->countSurfDiffs());
    
    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }
    
    stex::Tri * tri = pTris[tidx];
    
    uint ldidx = tri->patchdef()->surfdiffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }
    
    if (direction_tri == std::numeric_limits<uint>::max()) {
        return tri->sdiff(ldidx)->dcst();
        
    }
    else {
        int direction = tri->getTriDirection(direction_tri);
        if (direction == -1) {
            std::ostringstream os;
            os << "Triangle " << direction_tri << " is not a neighbor of triangle " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        return tri->sdiff(ldidx)->dcst(direction);
    }
    
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriDiffD(uint tidx, uint didx, double dk, uint direction_tri)
{
    #ifdef DIRECTIONAL_DCST_DEBUG
    CLOG_IF(direction_tri !=  std::numeric_limits<uint>::max(), DEBUG, "steps_debug") << "tidx: " << tidx << " didx: " << didx << " dk: " << dk << " direction tri: " << direction_tri << "\n";
    
    CLOG_IF(direction_tri ==  std::numeric_limits<uint>::max(), DEBUG, "steps_debug") << "tidx: " << tidx << " didx: " << didx << " dk: " << dk << " for all directions.\n";
    #endif
    
    assert (tidx < pTris.size());
    assert (didx < statedef()->countSurfDiffs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    uint ldidx = tri->patchdef()->surfdiffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }
    
    if (direction_tri == std::numeric_limits<uint>::max()) {
        tri->sdiff(ldidx)->setDcst(dk);

    }
    else {
        int direction = tri->getTriDirection(direction_tri);
        if (direction == -1) {
            std::ostringstream os;
            os << "Triangle " << direction_tri << " is not a neighbor of triangle " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        #ifdef DIRECTIONAL_DCST_DEBUG
        CLOG(DEBUG, "steps_debug") << "use tri " << direction_tri << " to set direction " << direction << ".\n";
        #endif
        
        tri->sdiff(ldidx)->setDirectionDcst(direction, dk);
    }
    _updateElement(tri->sdiff(ldidx));
    _updateSum();
    
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTriVDepSReacActive(uint tidx, uint vsridx) const
{
    assert (tidx < pTris.size());
    assert (vsridx < statedef()->countVDepSReacs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    uint lvsridx = tri->patchdef()->vdepsreacG2L(vsridx);
    if (lvsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Voltage-dependent surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    if (tri->vdepsreac(lvsridx)->inactive() == true) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriVDepSReacActive(uint tidx, uint vsridx, bool act)
{
    assert (tidx < pTris.size());
    assert (vsridx < statedef()->countVDepSReacs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    uint lvsridx = tri->patchdef()->vdepsreacG2L(vsridx);
    if (lvsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Voltage-dependent surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    tri->vdepsreac(lvsridx)->setActive(act);

    _updateElement(tri->vdepsreac(lvsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriSReacH(uint tidx, uint ridx) const
{
    assert (tidx < pTris.size());
    assert (ridx < statedef()->countSReacs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    uint lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    return tri->sreac(lsridx)->h();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriSReacC(uint tidx, uint ridx) const
{
    assert (tidx < pTris.size());
    assert (ridx < statedef()->countSReacs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    uint lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    return tri->sreac(lsridx)->c();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriSReacA(uint tidx, uint ridx) const
{
    assert (tidx < pTris.size());
    assert (ridx < statedef()->countSReacs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    uint lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    return tri->sreac(lsridx)->rate();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setEfieldDT(double efdt)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    if (efdt <= 0.0)
    {
        std::ostringstream os;
        os << "EField dt must be graeter than zero.";
        throw steps::ArgErr(os.str());
    }
    pEFDT = efdt;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetV(uint tidx) const
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    int loctidx = pEFTet_GtoL[tidx];
    if (loctidx == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        throw steps::ArgErr(os.str());
    }

    // EField object should convert value to base s.i. units
    return pEField->getTetV(loctidx);

}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetV(uint tidx, double v)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    int loctidx = pEFTet_GtoL[tidx];
    if (loctidx == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        throw steps::ArgErr(os.str());
    }

    // EField object should convert to millivolts
    pEField->setTetV(loctidx, v);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTetVClamped(uint tidx) const
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    int loctidx = pEFTet_GtoL[tidx];
    if (loctidx == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        throw steps::ArgErr(os.str());
    }

    return pEField->getTetVClamped(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetVClamped(uint tidx, bool cl)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    int loctidx = pEFTet_GtoL[tidx];
    if (loctidx == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        throw steps::ArgErr(os.str());
    }

    pEField->setTetVClamped(loctidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriV(uint tidx) const
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    int loctidx = pEFTri_GtoL[tidx];
    if (loctidx == -1)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        throw steps::ArgErr(os.str());
    }
    // EField object should convert value to base s.i. units
    return pEField->getTriV(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriV(uint tidx, double v)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    int loctidx = pEFTri_GtoL[tidx];
    if (loctidx == -1)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        throw steps::ArgErr(os.str());
    }

    // EField object should convert to millivolts
    pEField->setTriV(loctidx, v);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTriVClamped(uint tidx) const
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    int loctidx = pEFTri_GtoL[tidx];
    if (loctidx == -1)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        throw steps::ArgErr(os.str());
    }

    return pEField->getTriVClamped(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriVClamped(uint tidx, bool cl)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    int loctidx = pEFTri_GtoL[tidx];
    if (loctidx == -1)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        throw steps::ArgErr(os.str());
    }

    pEField->setTriVClamped(loctidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriOhmicI(uint tidx)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    int loctidx = pEFTri_GtoL[tidx];
    if (loctidx == -1)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        throw steps::ArgErr(os.str());
    }

    return tri->getOhmicI(pEField->getTriV(loctidx), efdt());
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriOhmicI(uint tidx, uint ocidx)
{
    assert (tidx < pTris.size());
    assert (ocidx < statedef()->countOhmicCurrs());

    int loctidx = pEFTri_GtoL[tidx];
    if (loctidx == -1)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    uint locidx = tri->patchdef()->ohmiccurrG2L(ocidx);
    if (locidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Ohmic current undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    return tri->getOhmicI(locidx, pEField->getTriV(loctidx), efdt());
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriGHKI(uint tidx)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    return tri->getGHKI(efdt());
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriGHKI(uint tidx, uint ghkidx)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }

    stex::Tri * tri = pTris[tidx];

    uint locidx = tri->patchdef()->ghkcurrG2L(ghkidx);
    if (locidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "GHK current undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    return tri->getGHKI(locidx, efdt());
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriI(uint tidx) const
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }

    int loctidx = pEFTri_GtoL[tidx];
    if (loctidx == -1)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        throw steps::ArgErr(os.str());
    }

    // EField object should convert to required units
    return pEField->getTriI(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setVertIClamp(uint vidx, double cur)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    int locvidx = pEFVert_GtoL[vidx];
    if (locvidx == -1)
    {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        throw steps::ArgErr(os.str());
    }

    // EField object should convert to required units
    pEField->setVertIClamp(locvidx, cur);
}


////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriIClamp(uint tidx, double cur)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    int loctidx = pEFTri_GtoL[tidx];
    if (loctidx == -1)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        throw steps::ArgErr(os.str());
    }

    // EField object should convert to required units
    pEField->setTriIClamp(loctidx, cur);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getVertV(uint vidx) const
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    int locvidx = pEFVert_GtoL[vidx];
    if (locvidx == -1)
    {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        throw steps::ArgErr(os.str());
    }
    // EField object should convert value to base s.i. units
    return pEField->getVertV(locvidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setVertV(uint vidx, double v)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    int locvidx = pEFVert_GtoL[vidx];
    if (locvidx == -1)
    {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        throw steps::ArgErr(os.str());
    }
    // EField object should convert to millivolts
    pEField->setVertV(locvidx, v);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getVertVClamped(uint vidx) const
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    int locvidx = pEFVert_GtoL[vidx];
    if (locvidx == -1)
    {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        throw steps::ArgErr(os.str());
    }

    return pEField->getVertVClamped(locvidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setVertVClamped(uint vidx, bool cl)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    int locvidx = pEFVert_GtoL[vidx];
    if (locvidx == -1)
    {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        throw steps::ArgErr(os.str());
    }
    // EField object should convert to millivolts
    pEField->setVertVClamped(locvidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setMembRes(uint midx, double ro, double vrev)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    if (ro <= 0.0)
    {
        std::ostringstream os;
        os << "Resistivity must be greater than zero.";
        throw steps::ArgErr(os.str());
    }
    // EField object should convert to required units
    assert (midx == 0);
    pEField->setSurfaceResistivity(midx, ro, vrev);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setMembPotential(uint midx, double v)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    // EField object should convert to millivolts
    assert (midx == 0);
    pEField->setMembPotential(midx, v);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setMembCapac(uint midx, double cm)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    if (cm < 0.0)
    {
        std::ostringstream os;
        os << "Capacitance must be greater than or equal to zero.";
        throw steps::ArgErr(os.str());
    }


    // EField object should convert to required units
    assert (midx == 0);
    pEField->setMembCapac(midx, cm);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setMembVolRes(uint midx, double ro)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    if (ro < 0.0)
    {
        std::ostringstream os;
        os << "Resistivity must be greater than or equal to zero.";
        throw steps::ArgErr(os.str());
    }
    // EField object should convert to required units
    assert (midx == 0);
    pEField->setMembVolRes(midx, ro);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_updateElement(KProc* kp)
{

    double new_rate = kp->rate(this);

    CRKProcData & data = kp->crData;
    double old_rate = data.rate;

    data.rate = new_rate;

    if (old_rate == new_rate)  return;


    // new rate in positive groups
    if (new_rate >= 0.5) {

        // pow is the same
        int old_pow = data.pow;
        int new_pow;
        double temp = frexp(new_rate, &new_pow);

        if (old_pow == new_pow && data.recorded) {

            CRGroup* old_group = _getGroup(old_pow);

            old_group->sum += (new_rate - old_rate);
        }
        // pow is not the same
        else {
            data.pow = new_pow;


            if (data.recorded) {

                // remove old
                CRGroup* old_group = _getGroup(old_pow);
                (old_group->size) --;

                if (old_group->size == 0) old_group->sum = 0.0;
                else {
                    old_group->sum -= old_rate;

                    KProc* last = old_group->indices[old_group->size];
                    old_group->indices[data.pos] = last;
                    last->crData.pos = data.pos;
                }
            }

            // add new
            if (pGroups.size() <= new_pow) _extendPGroups(new_pow + 1);

            CRGroup* new_group = pGroups[new_pow];

            assert(new_group != NULL);
            if (new_group->size == new_group->capacity) _extendGroup(new_group);
            uint pos = new_group->size;
            new_group->indices[pos] = kp;
            new_group->size++;
            new_group->sum += new_rate;
            data.pos = pos;

        }
        data.recorded = true;

    }
    // new rate in negative group
    else if (new_rate < 0.5 && new_rate > 1e-20) {
        int old_pow = data.pow;
        int new_pow;
        double temp = frexp(new_rate, &new_pow);

        if (old_pow == new_pow && data.recorded) {

            CRGroup* old_group = _getGroup(old_pow);

            old_group->sum += (new_rate - old_rate);
        }
        // pow is not the same
        else {
            data.pow = new_pow;

            if (data.recorded) {
                CRGroup* old_group = _getGroup(old_pow);
                (old_group->size) --;

                if (old_group->size == 0) old_group->sum = 0.0;
                else {
                    old_group->sum -= old_rate;

                    KProc* last = old_group->indices[old_group->size];
                    old_group->indices[data.pos] = last;
                    last->crData.pos = data.pos;
                }
            }

            // add new

            if (nGroups.size() <= -new_pow) _extendNGroups(-new_pow + 1);

            CRGroup* new_group = nGroups[-new_pow];

            if (new_group->size == new_group->capacity) _extendGroup(new_group);
            uint pos = new_group->size;
            new_group->indices[pos] = kp;
            new_group->size++;
            new_group->sum += new_rate;
            data.pos = pos;

        }
        data.recorded = true;
    }

    else {

        if (data.recorded) {

            CRGroup* old_group = _getGroup(data.pow);

            // remove old
            old_group->size --;

            if (old_group->size == 0) old_group->sum = 0.0;
            else {
                old_group->sum -= old_rate;

                KProc* last = old_group->indices[old_group->size];
                old_group->indices[data.pos] = last;
                last->crData.pos = data.pos;
            }
        }
        data.recorded = false;
    }

}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// ROI Data Access
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stex::Tetexact::getBatchTetCounts(std::vector<uint> const & tets, std::string const & s) const
{
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;

    uint ntets = tets.size();
    uint sgidx = statedef()->getSpecIdx(s);
    std::vector<double> data(ntets, 0.0);

    for (uint t = 0; t < ntets; t++) {
        uint tidx = tets[t];

        if (tidx >= pTets.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }

        if (pTets[tidx] == 0)
        {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        stex::Tet * tet = pTets[tidx];
        uint slidx = tet->compdef()->specG2L(sgidx);
        if (slidx == ssolver::LIDX_UNDEFINED)
        {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }

        data[t] = tet->pools()[slidx];
    }

    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, fill in zeros at target positions:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        std::cerr << "Warning: Species " << s << " has not been defined in the following tetrahedrons, fill in zeros at target positions:\n";
        std::cerr << spec_undefined.str() << "\n";
    }
    return data;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stex::Tetexact::getBatchTriCounts(std::vector<uint> const & tris, std::string const & s) const
{
    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;


    uint ntris = tris.size();
    uint sgidx = statedef()->getSpecIdx(s);
    std::vector<double> data(ntris, 0.0);

    for (uint t = 0; t < ntris; t++) {
        uint tidx = tris[t];

        if (tidx >= pTris.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }

        if (pTris[tidx] == 0)
        {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        stex::Tri * tri = pTris[tidx];
        uint slidx = tri->patchdef()->specG2L(sgidx);
        if (slidx == ssolver::LIDX_UNDEFINED)
        {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }

        data[t] = tri->pools()[slidx];
    }

    if (has_tri_warning) {
        std::cerr << "Warning: The following triangles have not been assigned to a patch, fill in zeros at target positions:\n";
        std::cerr << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        std::cerr << "Warning: Species " << s << " has not been defined in the following triangles, fill in zeros at target positions:\n";
        std::cerr << spec_undefined.str() << "\n";
    }
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::getBatchTetCountsNP(unsigned int* indices, int input_size, std::string const & s, double* counts, int output_size) const
{
    if (input_size != output_size)
    {
        std::ostringstream os;
        os << "Error: output array (counts) size should be the same as input array (indices) size.\n";
        throw steps::ArgErr(os.str());
    }

    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;

    uint sgidx = statedef()->getSpecIdx(s);

    for (int t = 0; t < input_size; t++) {
        uint tidx = indices[t];

        if (tidx >= pTets.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }

        if (pTets[tidx] == 0)
        {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }

        stex::Tet * tet = pTets[tidx];
        uint slidx = tet->compdef()->specG2L(sgidx);
        if (slidx == ssolver::LIDX_UNDEFINED)
        {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }

        counts[t] = tet->pools()[slidx];
    }

    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, fill in zeros at target positions:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        std::cerr << "Warning: Species " << s << " has not been defined in the following tetrahedrons, fill in zeros at target positions:\n";
        std::cerr << spec_undefined.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::getBatchTriCountsNP(unsigned int* indices, int input_size, std::string const & s, double* counts, int output_size) const
{
    if (input_size != output_size)
    {
        std::ostringstream os;
        os << "Error: output array (counts) size should be the same as input array (indices) size.\n";
        throw steps::ArgErr(os.str());
    }

    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;


    uint sgidx = statedef()->getSpecIdx(s);

    for (int t = 0; t < input_size; t++) {
        uint tidx = indices[t];

        if (tidx >= pTris.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }

        if (pTris[tidx] == 0)
        {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }

        stex::Tri * tri = pTris[tidx];
        uint slidx = tri->patchdef()->specG2L(sgidx);
        if (slidx == ssolver::LIDX_UNDEFINED)
        {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }

        counts[t] = tri->pools()[slidx];
    }

    if (has_tri_warning) {
        std::cerr << "Warning: The following triangles have not been assigned to a patch, fill in zeros at target positions:\n";
        std::cerr << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        std::cerr << "Warning: Species " << s << " has not been defined in the following triangles, fill in zeros at target positions:\n";
        std::cerr << spec_undefined.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// ROI Data Access
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stex::Tetexact::getROITetCounts(std::string ROI_id, std::string const & s) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int inputsize = mesh()->getROIDataSize(ROI_id);
    std::vector<double> data(inputsize);
    
    getBatchTetCountsNP(indices, inputsize, s, const_cast<double*>(&data.front()), data.size());
    return data;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stex::Tetexact::getROITriCounts(std::string ROI_id, std::string const & s) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TRI)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int inputsize = mesh()->getROIDataSize(ROI_id);
    std::vector<double> data(inputsize);
    
    getBatchTriCountsNP(indices, inputsize, s, const_cast<double*>(&data.front()), data.size());
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::getROITetCountsNP(std::string ROI_id, std::string const & s, double* counts, int output_size) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int inputsize = mesh()->getROIDataSize(ROI_id);
    
    getBatchTetCountsNP(indices, inputsize, s, counts, output_size);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::getROITriCountsNP(std::string ROI_id, std::string const & s, double* counts, int output_size) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TRI)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int inputsize = mesh()->getROIDataSize(ROI_id);
    
    getBatchTriCountsNP(indices, inputsize, s, counts, output_size);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::getROIVol(std::string ROI_id) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    double sum = 0.0;
    for (uint t = 0; t < datasize; t++) {
        sum += pTets[indices[t]]->vol();
    }
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::getROIArea(std::string ROI_id) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TRI)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    double sum = 0.0;
    for (uint t = 0; t < datasize; t++) {
        sum += pTris[indices[t]]->area();
    }
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::getROICount(std::string ROI_id, std::string const & s) const
{
    steps::tetmesh::ElementType type = mesh()->getROIType(ROI_id);
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    double sum = 0.0;
    
    if (type == steps::tetmesh::ELEM_TRI) {
        bool has_tri_warning = false;
        bool has_spec_warning = false;
        std::ostringstream tri_not_assign;
        std::ostringstream spec_undefined;
        
        uint sgidx = statedef()->getSpecIdx(s);
        
        for (uint t = 0; t < datasize; t++) {
            uint tidx = indices[t];
            
            if (tidx >= pTris.size())
            {
                std::ostringstream os;
                os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
                throw steps::ArgErr(os.str());
            }
            
            if (pTris[tidx] == 0)
            {
                tri_not_assign << tidx << " ";
                has_tri_warning = true;
                continue;
            }
            
            stex::Tri * tri = pTris[tidx];
            uint slidx = tri->patchdef()->specG2L(sgidx);
            if (slidx == ssolver::LIDX_UNDEFINED)
            {
                spec_undefined << tidx << " ";
                has_spec_warning = true;
                continue;
            }
            
            sum += tri->pools()[slidx];
        }
        
        if (has_tri_warning) {
            std::cerr << "Warning: The following triangles have not been assigned to a patch, fill in zeros at target positions:\n";
            std::cerr << tri_not_assign.str() << "\n";
        }
        
        if (has_spec_warning) {
            std::cerr << "Warning: Species " << s << " has not been defined in the following triangles, fill in zeros at target positions:\n";
            std::cerr << spec_undefined.str() << "\n";
        }
    }
    else if (type == steps::tetmesh::ELEM_TET) {
        bool has_tet_warning = false;
        bool has_spec_warning = false;
        std::ostringstream tet_not_assign;
        std::ostringstream spec_undefined;
        
        uint sgidx = statedef()->getSpecIdx(s);
        
        for (uint t = 0; t < datasize; t++) {
            uint tidx = indices[t];
            
            if (tidx >= pTets.size())
            {
                std::ostringstream os;
                os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
                throw steps::ArgErr(os.str());
            }
            
            if (pTets[tidx] == 0)
            {
                tet_not_assign << tidx << " ";
                has_tet_warning = true;
                continue;
            }
            
            stex::Tet * tet = pTets[tidx];
            uint slidx = tet->compdef()->specG2L(sgidx);
            if (slidx == ssolver::LIDX_UNDEFINED)
            {
                spec_undefined << tidx << " ";
                has_spec_warning = true;
                continue;
            }
            
            sum += tet->pools()[slidx];
        }
        
        if (has_tet_warning) {
            std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, fill in zeros at target positions:\n";
            std::cerr << tet_not_assign.str() << "\n";
        }
        
        if (has_spec_warning) {
            std::cerr << "Warning: Species " << s << " has not been defined in the following tetrahedrons, fill in zeros at target positions:\n";
            std::cerr << spec_undefined.str() << "\n";
        }
    }
    else {
        std::ostringstream os;
        os << "Error: Cannot find suitable ROI for the function call getROICount.\n";
        throw steps::ArgErr(os.str());
    }
    
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setROICount(std::string ROI_id, std::string const & s, double count)
{
    steps::tetmesh::ElementType roi_type = mesh()->getROIType(ROI_id);
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    double total_weight = 0.0;
    uint sgidx = statedef()->getSpecIdx(s);
        
    switch (roi_type) {
    case steps::tetmesh::ELEM_TRI:
        {
            std::vector<Tri *> apply;
            for (uint t = 0; t < datasize; t++) {
                uint tidx = indices[t];
                if (tidx >= pTris.size())
                    throw steps::ArgErr("ROI refers to nonexistent triangle "+std::to_string(tidx));
                
                stex::Tri *tri = pTris[tidx];
                if (!tri || tri->patchdef()->specG2L(sgidx) == ssolver::LIDX_UNDEFINED)
                    continue;
                
                apply.push_back(tri);
                total_weight += tri->area();
            }

            steps::util::distribute_quantity(count, apply.begin(), apply.end(),
                [](Tri *tri) { return tri->area(); },
                [sgidx](Tri *tri, uint c) { tri->setCount(tri->patchdef()->specG2L(sgidx), c); },
                [sgidx](Tri *tri, int c)  { tri->incCount(tri->patchdef()->specG2L(sgidx), c); },
                *rng(),
                total_weight);

            for (auto &tri: apply) _updateSpec(tri, tri->patchdef()->specG2L(sgidx));
        }
        break;

    case steps::tetmesh::ELEM_TET:
        {
            std::vector<Tet *> apply;
            for (uint t = 0; t < datasize; t++) {
                uint tidx = indices[t];
                if (tidx >= pTets.size())
                    throw steps::ArgErr("ROI refers to nonexistent tetrahedron "+std::to_string(tidx));

                stex::Tet *tet = pTets[tidx];
                if (!tet || tet->compdef()->specG2L(sgidx) == ssolver::LIDX_UNDEFINED)
                    continue;
                
                apply.push_back(tet);
                total_weight += tet->vol();
            }

            steps::util::distribute_quantity(count, apply.begin(), apply.end(),
                [](Tet *tet) { return tet->vol(); },
                [sgidx](Tet *tet, uint c) { tet->setCount(tet->compdef()->specG2L(sgidx), c); },
                [sgidx](Tet *tet, int c)  { tet->incCount(tet->compdef()->specG2L(sgidx), c); },
                *rng(),
                total_weight);

            for (auto &tet: apply) _updateSpec(tet, tet->compdef()->specG2L(sgidx));
        }
        break;

    default:
        throw steps::ArgErr("can only set counts in tetrahedra or triangle ROIs");
    }
}

////////////////////////////////////////////////////////////////////////////////


double stex::Tetexact::getROIAmount(std::string ROI_id, std::string const & s) const
{
    double count = getROICount(ROI_id, s);
    return (count / smath::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::getROIConc(std::string ROI_id, std::string const & s) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    double count = getROICount(ROI_id, s);
    double vol = getROIVol(ROI_id);
    return count/ (1.0e3 * vol * steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setROIConc(std::string ROI_id, std::string const & s, double conc)
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET))
        throw steps::ArgErr("can only set concentrations in tetrahedra ROIs");

    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    double total_weight = 0.0;
    uint sgidx = statedef()->getSpecIdx(s);

    std::vector<Tet *> apply;
    for (uint t = 0; t < datasize; t++) {
        uint tidx = indices[t];
        if (tidx >= pTets.size())
            throw steps::ArgErr("ROI refers to nonexistent tetrahedron "+std::to_string(tidx));

        stex::Tet *tet = pTets[tidx];
        if (!tet || tet->compdef()->specG2L(sgidx) == ssolver::LIDX_UNDEFINED)
            continue;
        
        apply.push_back(tet);
        total_weight += tet->vol();
    }

    double count = conc * (1.0e3 * total_weight * steps::math::AVOGADRO);

    steps::util::distribute_quantity(count, apply.begin(), apply.end(),
        [](Tet *tet) { return tet->vol(); },
        [sgidx](Tet *tet, uint c) { tet->setCount(tet->compdef()->specG2L(sgidx), c); },
        [sgidx](Tet *tet, int c)  { tet->incCount(tet->compdef()->specG2L(sgidx), c); },
        *rng(),
        total_weight);

    for (auto &tet: apply) _updateSpec(tet, tet->compdef()->specG2L(sgidx));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setROIClamped(std::string ROI_id, std::string const & s, bool b)
{
    steps::tetmesh::ElementType type = mesh()->getROIType(ROI_id);
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    if (type == steps::tetmesh::ELEM_TRI) {
        bool has_tri_warning = false;
        bool has_spec_warning = false;
        std::ostringstream tri_not_assign;
        std::ostringstream spec_undefined;
        
        uint sgidx = statedef()->getSpecIdx(s);
        
        for (uint t = 0; t < datasize; t++) {
            uint tidx = indices[t];
            
            if (tidx >= pTris.size())
            {
                std::ostringstream os;
                os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
                throw steps::ArgErr(os.str());
            }
            
            if (pTris[tidx] == 0)
            {
                tri_not_assign << tidx << " ";
                has_tri_warning = true;
                continue;
            }
            
            stex::Tri * tri = pTris[tidx];
            uint slidx = tri->patchdef()->specG2L(sgidx);
            if (slidx == ssolver::LIDX_UNDEFINED)
            {
                spec_undefined << tidx << " ";
                has_spec_warning = true;
                continue;
            }
            
            tri->setClamped(slidx, b);
        }
        
        if (has_tri_warning) {
            std::cerr << "Warning: The following triangles have not been assigned to a patch, fill in zeros at target positions:\n";
            std::cerr << tri_not_assign.str() << "\n";
        }
        
        if (has_spec_warning) {
            std::cerr << "Warning: Species " << s << " has not been defined in the following triangles, fill in zeros at target positions:\n";
            std::cerr << spec_undefined.str() << "\n";
        }
    }
    else if (type == steps::tetmesh::ELEM_TET) {
        bool has_tet_warning = false;
        bool has_spec_warning = false;
        std::ostringstream tet_not_assign;
        std::ostringstream spec_undefined;
        
        uint sgidx = statedef()->getSpecIdx(s);
        
        for (uint t = 0; t < datasize; t++) {
            uint tidx = indices[t];
            
            if (tidx >= pTets.size())
            {
                std::ostringstream os;
                os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
                throw steps::ArgErr(os.str());
            }
            
            if (pTets[tidx] == 0)
            {
                tet_not_assign << tidx << " ";
                has_tet_warning = true;
                continue;
            }
            
            stex::Tet * tet = pTets[tidx];
            uint slidx = tet->compdef()->specG2L(sgidx);
            if (slidx == ssolver::LIDX_UNDEFINED)
            {
                spec_undefined << tidx << " ";
                has_spec_warning = true;
                continue;
            }
            
            tet->setClamped(slidx, b);
        }
        
        if (has_tet_warning) {
            std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, fill in zeros at target positions:\n";
            std::cerr << tet_not_assign.str() << "\n";
        }
        
        if (has_spec_warning) {
            std::cerr << "Warning: Species " << s << " has not been defined in the following tetrahedrons, fill in zeros at target positions:\n";
            std::cerr << spec_undefined.str() << "\n";
        }
    }
    else {
        std::ostringstream os;
        os << "Error: Cannot find suitable ROI for the function call getROICount.\n";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setROIReacK(std::string ROI_id, std::string const & r, double kf)
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    bool has_tet_warning = false;
    bool has_reac_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream reac_undefined;
    
    uint rgidx = statedef()->getReacIdx(r);
    
    for (uint t = 0; t < datasize; t++) {
        uint tidx = indices[t];
        
        if (tidx >= pTets.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        if (pTets[tidx] == 0)
        {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }
        
        stex::Tet * tet = pTets[tidx];
        uint rlidx = tet->compdef()->reacG2L(rgidx);
        if (rlidx == ssolver::LIDX_UNDEFINED)
        {
            reac_undefined << tidx << " ";
            has_reac_warning = true;
            continue;
        }
        
        tet->reac(rlidx)->setKcst(kf);
    }
    
    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, no change is applied to them:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }
    
    if (has_reac_warning) {
        std::cerr << "Warning: Reac " << r << " has not been defined in the following tetrahedrons, no change is applied to them:\n";
        std::cerr << reac_undefined.str() << "\n";
    }
    
    _update();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setROISReacK(std::string ROI_id, std::string const & sr, double kf)
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TRI)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    bool has_tri_warning = false;
    bool has_sreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream sreac_undefined;
    
    uint srgidx = statedef()->getSReacIdx(sr);
    
    for (uint t = 0; t < datasize; t++) {
        uint tidx = indices[t];
        
        if (tidx >= pTris.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        if (pTris[tidx] == 0)
        {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }
        
        stex::Tri * tri = pTris[tidx];
        uint srlidx = tri->patchdef()->sreacG2L(srgidx);
        if (srlidx == ssolver::LIDX_UNDEFINED)
        {
            sreac_undefined << tidx << " ";
            has_sreac_warning = true;
            continue;
        }
        
        tri->sreac(srlidx)->setKcst(kf);
    }
    
    if (has_tri_warning) {
        std::cerr << "Warning: The following triangles have not been assigned to a patch, no change is applied to them:\n";
        std::cerr << tri_not_assign.str() << "\n";
    }
    
    if (has_sreac_warning) {
        std::cerr << "Warning: SReac " << sr << " has not been defined in the following patch, no change is applied to them:\n";
        std::cerr << sreac_undefined.str() << "\n";
    }
    
    _update();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setROIDiffD(std::string ROI_id, std::string const & d, double dk)
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    bool has_tet_warning = false;
    bool has_diff_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream diff_undefined;
    
    uint dgidx = statedef()->getDiffIdx(d);
    
    for (uint t = 0; t < datasize; t++) {
        uint tidx = indices[t];
        
        if (tidx >= pTets.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        if (pTets[tidx] == 0)
        {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }
        
        stex::Tet * tet = pTets[tidx];
        uint dlidx = tet->compdef()->diffG2L(dgidx);
        if (dlidx == ssolver::LIDX_UNDEFINED)
        {
            diff_undefined << tidx << " ";
            has_diff_warning = true;
            continue;
        }
        
        tet->diff(dlidx)->setDcst(dk);
    }
    
    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, no change is applied to them:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }
    
    if (has_diff_warning) {
        std::cerr << "Warning: Diff " << d << " has not been defined in the following tetrahedrons, no change is applied to them:\n";
        std::cerr << diff_undefined.str() << "\n";
    }
    
    _update();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setROIReacActive(std::string ROI_id, std::string const & r, bool a)
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    bool has_tet_warning = false;
    bool has_reac_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream reac_undefined;
    
    uint rgidx = statedef()->getReacIdx(r);
    
    for (uint t = 0; t < datasize; t++) {
        uint tidx = indices[t];
        
        if (tidx >= pTets.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        if (pTets[tidx] == 0)
        {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }
        
        stex::Tet * tet = pTets[tidx];
        uint rlidx = tet->compdef()->reacG2L(rgidx);
        if (rlidx == ssolver::LIDX_UNDEFINED)
        {
            reac_undefined << tidx << " ";
            has_reac_warning = true;
            continue;
        }
        
        tet->reac(rlidx)->setActive(a);
    }
    
    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, no change is applied to them:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }
    
    if (has_reac_warning) {
        std::cerr << "Warning: Reac " << r << " has not been defined in the following tetrahedrons, no change is applied to them:\n";
        std::cerr << reac_undefined.str() << "\n";
    }
    
    _update();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setROISReacActive(std::string ROI_id, std::string const & sr, bool a)
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TRI)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    bool has_tri_warning = false;
    bool has_sreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream sreac_undefined;
    
    uint srgidx = statedef()->getSReacIdx(sr);
    
    for (uint t = 0; t < datasize; t++) {
        uint tidx = indices[t];
        
        if (tidx >= pTris.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        if (pTris[tidx] == 0)
        {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }
        
        stex::Tri * tri = pTris[tidx];
        uint srlidx = tri->patchdef()->sreacG2L(srgidx);
        if (srlidx == ssolver::LIDX_UNDEFINED)
        {
            sreac_undefined << tidx << " ";
            has_sreac_warning = true;
            continue;
        }
        
        tri->sreac(srlidx)->setActive(a);
    }
    
    if (has_tri_warning) {
        std::cerr << "Warning: The following triangles have not been assigned to a patch, no change is applied to them:\n";
        std::cerr << tri_not_assign.str() << "\n";
    }
    
    if (has_sreac_warning) {
        std::cerr << "Warning: SReac " << sr << " has not been defined in the following patch, no change is applied to them:\n";
        std::cerr << sreac_undefined.str() << "\n";
    }
    
    _update();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setROIDiffActive(std::string ROI_id, std::string const & d, bool a)
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    bool has_tet_warning = false;
    bool has_diff_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream diff_undefined;
    
    uint dgidx = statedef()->getDiffIdx(d);
    
    for (uint t = 0; t < datasize; t++) {
        uint tidx = indices[t];
        
        if (tidx >= pTets.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        if (pTets[tidx] == 0)
        {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }
        
        stex::Tet * tet = pTets[tidx];
        uint dlidx = tet->compdef()->diffG2L(dgidx);
        if (dlidx == ssolver::LIDX_UNDEFINED)
        {
            diff_undefined << tidx << " ";
            has_diff_warning = true;
            continue;
        }
        
        tet->diff(dlidx)->setActive(a);
    }
    
    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, no change is applied to them:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }
    
    if (has_diff_warning) {
        std::cerr << "Warning: Diff " << d << " has not been defined in the following tetrahedrons, no change is applied to them:\n";
        std::cerr << diff_undefined.str() << "\n";
    }
    
    _update();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setROIVDepSReacActive(std::string ROI_id, std::string const & vsr, bool a)
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TRI)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    bool has_tri_warning = false;
    bool has_vsreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream vsreac_undefined;
    
    uint vsrgidx = statedef()->getVDepSReacIdx(vsr);
    
    for (uint t = 0; t < datasize; t++) {
        uint tidx = indices[t];
        
        if (tidx >= pTris.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        if (pTris[tidx] == 0)
        {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }
        
        stex::Tri * tri = pTris[tidx];
        uint vsrlidx = tri->patchdef()->vdepsreacG2L(vsrgidx);
        if (vsrlidx == ssolver::LIDX_UNDEFINED)
        {
            vsreac_undefined << tidx << " ";
            has_vsreac_warning = true;
            continue;
        }
        
        tri->vdepsreac(vsrlidx)->setActive(a);
    }
    
    if (has_tri_warning) {
        std::cerr << "Warning: The following triangles have not been assigned to a patch, no change is applied to them:\n";
        std::cerr << tri_not_assign.str() << "\n";
    }
    
    if (has_vsreac_warning) {
        std::cerr << "Warning: VDepSReac " << vsr << " has not been defined in the following patch, no change is applied to them:\n";
        std::cerr << vsreac_undefined.str() << "\n";
    }
    _update();
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::getROIReacExtent(std::string ROI_id, std::string const & r) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    bool has_tet_warning = false;
    bool has_reac_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream reac_undefined;
    
    uint rgidx = statedef()->getReacIdx(r);
    
    uint sum = 0;
    
    for (uint t = 0; t < datasize; t++) {
        uint tidx = indices[t];
        
        if (tidx >= pTets.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        if (pTets[tidx] == 0)
        {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }
        
        stex::Tet * tet = pTets[tidx];
        uint rlidx = tet->compdef()->reacG2L(rgidx);
        if (rlidx == ssolver::LIDX_UNDEFINED)
        {
            reac_undefined << tidx << " ";
            has_reac_warning = true;
            continue;
        }
        
        sum += tet->reac(rlidx)->getExtent();
    }
    
    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, no change is applied to them:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }
    
    if (has_reac_warning) {
        std::cerr << "Warning: Reac " << r << " has not been defined in the following tetrahedrons, no change is applied to them:\n";
        std::cerr << reac_undefined.str() << "\n";
    }
    
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::resetROIReacExtent(std::string ROI_id, std::string const & r)
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    bool has_tet_warning = false;
    bool has_reac_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream reac_undefined;
    
    uint rgidx = statedef()->getReacIdx(r);
    
    for (uint t = 0; t < datasize; t++) {
        uint tidx = indices[t];
        
        if (tidx >= pTets.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        if (pTets[tidx] == 0)
        {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }
        
        stex::Tet * tet = pTets[tidx];
        uint rlidx = tet->compdef()->reacG2L(rgidx);
        if (rlidx == ssolver::LIDX_UNDEFINED)
        {
            reac_undefined << tidx << " ";
            has_reac_warning = true;
            continue;
        }
        
        tet->reac(rlidx)->resetExtent();
    }
    
    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, no change is applied to them:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }
    
    if (has_reac_warning) {
        std::cerr << "Warning: Reac " << r << " has not been defined in the following tetrahedrons, no change is applied to them:\n";
        std::cerr << reac_undefined.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::getROISReacExtent(std::string ROI_id, std::string const & sr) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TRI)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    bool has_tri_warning = false;
    bool has_sreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream sreac_undefined;
    
    uint srgidx = statedef()->getSReacIdx(sr);
    
    uint sum = 0;
    
    for (uint t = 0; t < datasize; t++) {
        uint tidx = indices[t];
        
        if (tidx >= pTris.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        if (pTris[tidx] == 0)
        {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }
        
        stex::Tri * tri = pTris[tidx];
        uint srlidx = tri->patchdef()->sreacG2L(srgidx);
        if (srlidx == ssolver::LIDX_UNDEFINED)
        {
            sreac_undefined << tidx << " ";
            has_sreac_warning = true;
            continue;
        }
        
        sum += tri->sreac(srlidx)->getExtent();
    }
    
    if (has_tri_warning) {
        std::cerr << "Warning: The following triangles have not been assigned to a patch, no change is applied to them:\n";
        std::cerr << tri_not_assign.str() << "\n";
    }
    
    if (has_sreac_warning) {
        std::cerr << "Warning: SReac " << sr << " has not been defined in the following patch, no change is applied to them:\n";
        std::cerr << sreac_undefined.str() << "\n";
    }
    
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::resetROISReacExtent(std::string ROI_id, std::string const & sr)
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TRI)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    bool has_tri_warning = false;
    bool has_sreac_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream sreac_undefined;
    
    uint srgidx = statedef()->getSReacIdx(sr);
    
    for (uint t = 0; t < datasize; t++) {
        uint tidx = indices[t];
        
        if (tidx >= pTris.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no triangle with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        if (pTris[tidx] == 0)
        {
            tri_not_assign << tidx << " ";
            has_tri_warning = true;
            continue;
        }
        
        stex::Tri * tri = pTris[tidx];
        uint srlidx = tri->patchdef()->sreacG2L(srgidx);
        if (srlidx == ssolver::LIDX_UNDEFINED)
        {
            sreac_undefined << tidx << " ";
            has_sreac_warning = true;
            continue;
        }
        
        tri->sreac(srlidx)->resetExtent();
    }
    
    if (has_tri_warning) {
        std::cerr << "Warning: The following triangles have not been assigned to a patch, no change is applied to them:\n";
        std::cerr << tri_not_assign.str() << "\n";
    }
    
    if (has_sreac_warning) {
        std::cerr << "Warning: SReac " << sr << " has not been defined in the following patch, no change is applied to them:\n";
        std::cerr << sreac_undefined.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::getROIDiffExtent(std::string ROI_id, std::string const & d) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    bool has_tet_warning = false;
    bool has_diff_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream diff_undefined;
    
    uint dgidx = statedef()->getDiffIdx(d);
    
    uint sum = 0;
    
    for (uint t = 0; t < datasize; t++) {
        uint tidx = indices[t];
        
        if (tidx >= pTets.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        if (pTets[tidx] == 0)
        {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }
        
        stex::Tet * tet = pTets[tidx];
        uint dlidx = tet->compdef()->diffG2L(dgidx);
        if (dlidx == ssolver::LIDX_UNDEFINED)
        {
            diff_undefined << tidx << " ";
            has_diff_warning = true;
            continue;
        }
        
        sum += tet->diff(dlidx)->getExtent();
    }
    
    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, no change is applied to them:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }
    
    if (has_diff_warning) {
        std::cerr << "Warning: Diff " << d << " has not been defined in the following tetrahedrons, no change is applied to them:\n";
        std::cerr << diff_undefined.str() << "\n";
    }
    
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::resetROIDiffExtent(std::string ROI_id, std::string const & d)
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    bool has_tet_warning = false;
    bool has_diff_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream diff_undefined;
    
    uint dgidx = statedef()->getDiffIdx(d);
    
    for (uint t = 0; t < datasize; t++) {
        uint tidx = indices[t];
        
        if (tidx >= pTets.size())
        {
            std::ostringstream os;
            os << "Error (Index Overbound): There is no tetrahedron with index " << tidx << ".\n";
            throw steps::ArgErr(os.str());
        }
        
        if (pTets[tidx] == 0)
        {
            tet_not_assign << tidx << " ";
            has_tet_warning = true;
            continue;
        }
        
        stex::Tet * tet = pTets[tidx];
        uint dlidx = tet->compdef()->diffG2L(dgidx);
        if (dlidx == ssolver::LIDX_UNDEFINED)
        {
            diff_undefined << tidx << " ";
            has_diff_warning = true;
            continue;
        }
        
        tet->diff(dlidx)->resetExtent();
    }
    
    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, no change is applied to them:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }
    
    if (has_diff_warning) {
        std::cerr << "Warning: Diff " << d << " has not been defined in the following tetrahedrons, no change is applied to them:\n";
        std::cerr << diff_undefined.str() << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

// END
