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


// Standard library headers.
#include <cmath>
#include <vector>
#include <map>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <queue>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <numeric>

#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// STEPS headers.
#include "steps/common.h"
#include "steps/mpi/mpi_common.hpp"
#include "steps/mpi/tetopsplit/tetopsplit.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
#include "steps/mpi/tetopsplit/tet.hpp"
#include "steps/mpi/tetopsplit/tri.hpp"
#include "steps/mpi/tetopsplit/reac.hpp"
#include "steps/mpi/tetopsplit/sreac.hpp"
#include "steps/mpi/tetopsplit/diff.hpp"
#include "steps/mpi/tetopsplit/sdiff.hpp"
#include "steps/mpi/tetopsplit/comp.hpp"
#include "steps/mpi/tetopsplit/patch.hpp"
#include "steps/mpi/tetopsplit/wmvol.hpp"
#include "steps/mpi/tetopsplit/ghkcurr.hpp"
#include "steps/mpi/tetopsplit/vdeptrans.hpp"
#include "steps/mpi/tetopsplit/vdepsreac.hpp"
#include "steps/mpi/tetopsplit/diffboundary.hpp"
#include "steps/mpi/tetopsplit/sdiffboundary.hpp"
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
#include "steps/solver/diffboundarydef.hpp"
#include "steps/solver/sdiffboundarydef.hpp"
#include "steps/geom/tetmesh.hpp"

#include "steps/solver/efield/efield.hpp"
#include "steps/solver/efield/dVsolver.hpp"
#include "steps/solver/efield/dVsolver_slu.hpp"
#ifdef USE_PETSC
#include "steps/solver/efield/dVsolver_petsc.hpp"
#endif
#include "steps/util/distribute.hpp"

#include "third_party/easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;
namespace ssolver = steps::solver;
namespace smath = steps::math;

using steps::math::point3d;

////////////////////////////////////////////////////////////////////////////////

void smtos::schedIDXSet_To_Vec(smtos::SchedIDXSet const & s, smtos::SchedIDXVec & v)
{
    v.resize(s.size());
    std::copy(s.begin(), s.end(), v.begin());
}

////////////////////////////////////////////////////////////////////////////////

smtos::TetOpSplitP::TetOpSplitP(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r,
        int calcMembPot, std::vector<uint> const &tet_hosts, std::map<uint, uint> const &tri_hosts,
        std::vector<uint> const &wm_hosts)
: API(m, g, r)
, pMesh(0)
, pKProcs()
, pComps()
, pCompMap()
, pPatches()
, pDiffBoundaries()
, pSDiffBoundaries()
, pTets()
, pTris()
, pWmVols()
, pA0(0.0)
, pEFoption(static_cast<EF_solver>(calcMembPot))
, pTemp(0.0)
, pEFDT(1.0e-5)
, pEFNVerts(0)
, pEFNTris(0)
, pEFTris_vec(0)
, pEFNTets(0)
, pEFVert_GtoL()
, pEFTri_GtoL()
, pEFTet_GtoL()
, pEFTri_LtoG()
, pEFTrisVStale(true)
, tetHosts(tet_hosts)
, triHosts(tri_hosts)
, wmHosts(wm_hosts)
, diffApplyThreshold(10)
, diffSep(0)
, sdiffSep(0)
, updPeriod(0.0)
, recomputeUpdPeriod(true)
, reacExtent(0.0)
, diffExtent(0.0)
, nIteration(0.0)
, rd()
, gen(rd())
, compTime(0.0)
, syncTime(0.0)
, idleTime(0.0)
, efieldTime(0.0)
, rdTime(0.0)
, dataExchangeTime(0.0)
{
    if (rng() == 0)
    {
        std::ostringstream os;
        os << "No RNG provided to solver initializer function";
        throw steps::ArgErr(os.str());
    }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    
    MPI_Comm_size(MPI_COMM_WORLD, &nHosts);
    
    
    // All initialization code now in _setup() to allow EField solver to be
    // derived and create EField local objects within the constructor
    _setup();
    
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

smtos::TetOpSplitP::~TetOpSplitP(void)
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
        delete[] pEFVert_GtoL;
        delete[] pEFTri_GtoL;
        delete[] pEFTet_GtoL;
        delete[] pEFTri_LtoG;
    }
}

///////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::checkpoint(std::string const & file_name)
{
    std::ostringstream os;
    os << "This function has not been implemented!";
    throw steps::NotImplErr(os.str());

    /*
    
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

    SDiffBoundaryPVecCI sdb_e = pSDiffBoundaries.end();
    for (SDiffBoundaryPVecCI sdb = pSDiffBoundaries.begin(); sdb != sdb_e; ++sdb) {
        (*sdb)->checkpoint(cp_file);
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
    //std::cout << "complete.\n";

    */
}

///////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::restore(std::string const & file_name)
{
    std::ostringstream os;
    os << "This function has not been implemented!";
    throw steps::NotImplErr(os.str());
    
    /*
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

    SDiffBoundaryPVecCI sdb_e = pSDiffBoundaries.end();
    for (SDiffBoundaryPVecCI sdb = pSDiffBoundaries.begin(); sdb != sdb_e; ++sdb) {
        (*sdb)->restore(cp_file);
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

    //std::cout << "complete.\n";

    */
}

////////////////////////////////////////////////////////////////////////////////


std::string smtos::TetOpSplitP::getSolverName(void) const
{
    return "Parallel TetOpSplit";
}

////////////////////////////////////////////////////////////////////////////////

std::string smtos::TetOpSplitP::getSolverDesc(void) const
{
    return "Parallel approximate stochastic method in tetrahedral mesh";
}

////////////////////////////////////////////////////////////////////////////////

std::string smtos::TetOpSplitP::getSolverAuthors(void) const
{
    return "Iain Hepburn, Weiliang Chen, Stefan Wils, Sam Yates";
}

////////////////////////////////////////////////////////////////////////////////

std::string smtos::TetOpSplitP::getSolverEmail(void) const
{
    return "steps.dev@gmail.com";
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setup(void)
{
    // Perform upcast.
    pMesh = dynamic_cast<steps::tetmesh::Tetmesh *>(geom());
    if (!pMesh)
        throw steps::ArgErr("Geometry description to steps::solver::Tetexact solver "
                "constructor is not a valid steps::tetmesh::Tetmesh object.");

    // First initialise the pTets, pTris vector, because
    // want tets and tris to maintain indexing from Geometry
    uint ntets = mesh()->countTets();
    uint ntris = mesh()->countTris();
    uint ncomps = mesh()->_countComps();

    pTets.assign(ntets, NULL);
    pTris.assign(ntris, NULL);
    pWmVols.assign(ncomps, NULL);
    diffSep = 0;
    sdiffSep = 0;
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

    // Create the surface diffusion boundaries
    ssolver::SDiffBoundarydefPVecCI sdb_end = statedef()->endSDiffBoundary();
    for (ssolver::SDiffBoundaryDefPVecCI sdb = statedef()->bgnSDiffBoundary(); sdb != sdb_end; ++sdb)
    {
        uint sdiffboundary_gidx = (*sdb)->gidx();
        uint sdiffb_idx = _addSDiffBoundary(*sdb);
        assert(sdiffboundary_gidx == sdiffb_idx);
    }

    uint npatches = pPatches.size();
    assert (mesh()->_countPatches() == npatches);
    int rk; MPI_Comm_rank(MPI_COMM_WORLD, &rk);

    for (uint p = 0; p < npatches; ++p)
    {
        // Add the tris for this patch
        // We have checked the indexing - p is the global index
        steps::wm::Patch *wmpatch = mesh()->_getPatch(p);

        // Perform upcast
        steps::tetmesh::TmPatch *tmpatch = dynamic_cast<steps::tetmesh::TmPatch*>(wmpatch);
        if (!tmpatch)
            throw steps::ArgErr("Well-mixed patches not supported in steps::solver::TetOpSplitP solver.");

        steps::mpi::tetopsplit::Patch *localpatch = pPatches[p];

        /*** Previous impl
        // Create a map between edges and adjacent tris in this patch
        std::map<uint, std::vector<uint> > bar2tri;
        for (uint tri: tmpatch->_getAllTriIndices()) 
        {
            const uint *bars = pMesh->_getTriBars(tri);
            for (int i = 0; i < 3; ++i)
                bar2tri[bars[i]].push_back(tri);
        }

        for (uint tri: tmpatch->_getAllTriIndices()) 
        ***/

        // Create a map between edges and adjacent tris in this patch
        std::map<uint, std::vector<uint> > bar2tri;

        // We need to go through all patches to record bar2tri mapping
        // for all connected triangle neighbors even they are in different
        // patches, because their information is needed for surface diffusion boundary

        for (uint bar_p = 0; bar_p < npatches; ++bar_p) {

            steps::tetmesh::TmPatch *bar_patch = dynamic_cast<steps::tetmesh::TmPatch*>(pMesh->_getPatch(bar_p));

            for (uint tri: bar_patch->_getAllTriIndices())
            {
                const uint *bars = pMesh->_getTriBars(tri);
                for (int i = 0; i < 3; ++i)
                    bar2tri[bars[i]].push_back(tri);
            }
        }

        auto tri_idxs = tmpatch->_getAllTriIndices();

#pragma omp parallel for
        for (int i = 0; i< tri_idxs.size(); ++i) 
        //*** end of new impl
        {
            auto tri = tri_idxs[i];
            assert (pMesh->getTriPatch(tri) == tmpatch);

            double area = pMesh->getTriArea(tri);

            // NB: Tri vertices may not be in consistent order, so use bar interface.
            const uint *tri_bars = pMesh->_getTriBars(tri);
            double l[3] = {0, 0, 0};

            for (uint j=0; j<3; ++j) {
                const uint *v = pMesh->_getBar(tri_bars[j]);
                l[j] = distance(pMesh->_getVertex(v[0]), pMesh->_getVertex(v[1]));
            }

            // Get neighboring tris
            std::vector<int> tris(3, -1);
            for (int j = 0; j < 3; ++j)
            {
                std::vector<uint> neighb_tris = bar2tri[tri_bars[j]];
                for (int k = 0; k < neighb_tris.size(); ++k)
                {
                    if (neighb_tris[k] == tri || pMesh->getTriPatch(neighb_tris[k])  == nullptr)
                        continue;

                    tris[j] = neighb_tris[k];
                    break;
                }
            }

            point3d baryc = pMesh->_getTriBarycenter(tri);

            double d[3] = {0, 0, 0};
            for (uint j = 0; j < 3; ++j) {
                if (tris[j]==-1) continue;
                d[j] = distance(baryc, pMesh->_getTriBarycenter(tris[j]));
            }

            const int *tri_tets = pMesh->_getTriTetNeighb(tri);
#pragma omp critical
                _addTri(tri, localpatch, area, l[0], l[1], l[2], d[0], d[1], d[2], tri_tets[0], tri_tets[1], tris[0], tris[1], tris[2]);
        }
    }

    ncomps = pComps.size();
    assert (mesh()->_countComps() == ncomps);

    for (uint c = 0; c < ncomps; ++c)
    {
        // Add the tets for this comp.
        // We have checked the indexing - c is the global index.
        steps::wm::Comp * wmcomp = mesh()->_getComp(c);
        // Perform upcast
        steps::tetmesh::TmComp *tmcomp = dynamic_cast<steps::tetmesh::TmComp*>(wmcomp);
        if (tmcomp) {
             steps::mpi::tetopsplit::Comp * localcomp = pComps[c];

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
            steps::mpi::tetopsplit::Comp * localcomp = pComps[c];
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
    // comp they do not talk to each other (see smtos::Tet::setNextTet())
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
                    steps::mpi::tetopsplit::Tet * tet_in = pTets[tetinner];
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
                    steps::mpi::tetopsplit::Tet * tet_out = pTets[tetouter];

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
    assert(ndiffbnds ==    mesh()->_countDiffBoundaries());

    for (uint db = 0; db < ndiffbnds; ++db)
    {
        steps::mpi::tetopsplit::DiffBoundary * localdiffb = pDiffBoundaries[db];

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

            steps::mpi::tetopsplit::Tet * tetA = _tet(tetAidx);
            steps::mpi::tetopsplit::Tet * tetB = _tet(tetBidx);
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


    // Now loop over the surface diffusion boundaries:
    // 1) get all the bars and get the two triangles
    // 2) figure out which direction is the direction for a triangle
    // 3) add the triangle and the direction to local object

    // This is here because we need all tris to have been assigned correctly
    // to patches. Check every one and set the patchA and patchB for the db
    uint nsdiffbnds = pSDiffBoundaries.size();
    assert(nsdiffbnds == mesh()->_countSDiffBoundaries());

    for (uint sdb = 0; sdb < nsdiffbnds; ++sdb)
    {
        steps::mpi::tetopsplit::SDiffBoundary * localsdiffb = pSDiffBoundaries[sdb];

        uint patchAidx = localsdiffb->def()->patcha();
        uint patchBidx = localsdiffb->def()->patchb();
        steps::solver::Patchdef * patchAdef = statedef()->patchdef(patchAidx);
        steps::solver::Patchdef * patchBdef = statedef()->patchdef(patchBidx);

        for (uint sdbbar: localsdiffb->def()->bars())
        {
            const int *bar_tris = pMesh->_getBarTriNeighb(sdbbar);

            int triAidx = bar_tris[0];
            int triBidx = bar_tris[1];
            assert(triAidx >= 0 && triBidx >= 0);

            steps::mpi::tetopsplit::Tri * triA = _tri(triAidx);
            steps::mpi::tetopsplit::Tri * triB = _tri(triBidx);
            assert(triA != 0 && triB != 0);

            steps::solver::Patchdef *triA_pdef = triA->patchdef();
            steps::solver::Patchdef *triB_pdef = triB->patchdef();
            assert(triA_pdef != 0);
            assert(triB_pdef != 0);

            if (triA_pdef != patchAdef)
            {
                assert(triB_pdef == patchAdef);
                assert(triA_pdef == patchBdef);
            }
            else
            {
                assert(triB_pdef == patchBdef);
                assert(triA_pdef == patchAdef);
            }

            // Ok, checks over, lets get down to business
            int direction_idx_a = -1;
            int direction_idx_b = -1;

            const uint *triA_bars = pMesh->_getTriBars(triAidx);
            const uint *triB_bars = pMesh->_getTriBars(triBidx);

            for (uint i = 0; i < 3; ++i)
            {
                if (triA_bars[i] == sdbbar)
                {
                    assert(direction_idx_a == -1);
                    direction_idx_a = i;
                }
                if (triB_bars[i] == sdbbar)
                {
                    assert(direction_idx_b == -1);
                    direction_idx_b = i;
                }
            }
            assert (direction_idx_a != -1);
            assert (direction_idx_b != -1);

            // Set the tetrahedron and direction to the Diff Boundary object
            localsdiffb->setTriDirection(triAidx, direction_idx_a);
            localsdiffb->setTriDirection(triBidx, direction_idx_b);
        }
        localsdiffb->setPatches(_patch(patchAidx), _patch(patchBidx));


        // Might as well copy the vectors because we need to index through
        std::vector<uint> tris = localsdiffb->getTris();
        std::vector<uint> tris_direction = localsdiffb->getTriDirection();

        ntris = tris.size();
        assert (ntris <= pTris.size());
        assert (tris_direction.size() == ntris);

        for (uint t = 0; t < ntris; ++t)
            _tri(tris[t])->setSDiffBndDirection(tris_direction[t]);
    }

    for (auto t: pTets)
        if (t) t->setupKProcs(this);

    for (auto wmv: pWmVols)
        if (wmv) wmv->setupKProcs(this);

    for (auto t: pTris)
        if (t) t->setupKProcs(this, efflag());

    // Resolve all dependencies

    // DEBUG: vector holds all possible tetrahedrons,
    // but they have not necessarily been added to a compartment.
    for (auto t: pTets)
        if (t && t->getInHost()) t->setupDeps();

    // Vector allows for all compartments to be well-mixed, so
    // hold null-pointer for mesh compartments
    for (auto wmv: pWmVols)
        if (wmv && wmv->getInHost()) wmv->setupDeps();

    // DEBUG: vector holds all possible triangles, but
    // only patch triangles are filled
    for (auto t: pTris)
        if (t && t->getInHost()) t->setupDeps();

    // Create EField structures if EField is to be calculated
    if (efflag() == true) _setupEField();

    for (auto tet : boundaryTets) {
        tet->setupBufferLocations();
    }
    for (auto tri : boundaryTris) {
        tri->setupBufferLocations();
    }
    // just in case
    neighbHosts.erase(myRank);
    nNeighbHosts = neighbHosts.size();
    
    // construct remote molecule change buffers
    remoteChanges.clear();
    for (auto neighbor : neighbHosts) {
        remoteChanges[neighbor] = std::vector<uint> ();
    }
    
    nEntries = pKProcs.size();
    diffSep=pDiffs.size();
    sdiffSep=pSDiffs.size();
    _updateLocal();
    
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setupEField(void)
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
    case EF_DV_SLUSYS:
        pEField = make_EField<dVSolverSLU>(MPI_COMM_WORLD);
        break;
#ifdef USE_PETSC 
    case EF_DV_PETSC:
        pEField = make_EField<dVSolverPETSC>();
        break;
#endif
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
    
    std::vector<uint> pEFTets(pEFNTets * 4);

    // All the triangles we will count here are the real membrane triangles,
    // virtual triangles will not require a capacitance.
    std::vector<uint> pEFTris(pEFNTris * 3);

    std::vector<double> pEFVerts(pEFNVerts * 3);

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

    pEFTris_vec.resize(pEFNTris);
    EFTrisV.resize(pEFNTris);

    EFTrisI_permuted.resize(pEFNTris);
    EFTrisI_idx.resize(pEFNTris);

    EFTrisI_offset.assign(nHosts,0);
    EFTrisI_count.assign(nHosts,0);

    std::vector<int> local_eftri_indices;
    for (uint eft = 0; eft < pEFNTris; ++eft)
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
        smtos::Tri *tri_p = pTris[triidx];
        pEFTris_vec[eft] = tri_p;

        int tri_host = tri_p->getHost();
        ++EFTrisI_count[tri_host];
        if (myRank == tri_host) local_eftri_indices.push_back(eft);
    }

    const int *count_begin=&EFTrisI_count[0];
    std::partial_sum(count_begin, count_begin+(nHosts-1), 1+&EFTrisI_offset[0]);

    assert(local_eftri_indices.size()==EFTrisI_count[myRank]);

    MPI_Allgatherv(&local_eftri_indices[0], (int)local_eftri_indices.size(), MPI_INT,
            &EFTrisI_idx[0], &EFTrisI_count[0], &EFTrisI_offset[0], MPI_INT, MPI_COMM_WORLD);
    
    pEField->initMesh(pEFNVerts, &(pEFVerts.front()), pEFNTris, &(pEFTris.front()), pEFNTets, &(pEFTets.front()), memb->_getOpt_method(), memb->_getOpt_file_name(), memb->_getSearch_percent());
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::saveMembOpt(std::string const & opt_file_name)
{
    if (myRank != 0) return;
    
    if  (efflag() != true)
    {
        std::ostringstream os;
        os << "saveMembOpt method only available if running EField ";
        throw steps::ArgErr(os.str());
    }

    pEField->saveOptimal(opt_file_name);

}

////////////////////////////////////////////////////////////////////////////////

uint smtos::TetOpSplitP::_addComp(steps::solver::Compdef * cdef)
{
    smtos::Comp * comp = new Comp(cdef);
    assert(comp != 0);
    uint compidx = pComps.size();
    pComps.push_back(comp);
    pCompMap[cdef] = comp;
    return compidx;
}

////////////////////////////////////////////////////////////////////////////////

uint smtos::TetOpSplitP::_addPatch(steps::solver::Patchdef * pdef)
{
    /* Comp * icomp = 0;
     Comp * ocomp = 0;
     if (pdef->icompdef()) icomp = pCompMap[pdef->icompdef()];
     if (pdef->ocompdef()) ocomp = pCompMap[pdef->ocompdef()];
     */
    smtos::Patch * patch = new Patch(pdef);
    assert(patch != 0);
    uint patchidx = pPatches.size();
    pPatches.push_back(patch);
    return patchidx;
}

////////////////////////////////////////////////////////////////////////////////

uint smtos::TetOpSplitP::_addDiffBoundary(steps::solver::DiffBoundarydef * dbdef)
{
    smtos::DiffBoundary * diffb = new DiffBoundary(dbdef);
    assert(diffb != 0);
    uint dbidx = pDiffBoundaries.size();
    pDiffBoundaries.push_back(diffb);
    return dbidx;
}

////////////////////////////////////////////////////////////////////////////////

uint smtos::TetOpSplitP::_addSDiffBoundary(steps::solver::SDiffBoundarydef * sdbdef)
{
	smtos::SDiffBoundary * sdiffb = new SDiffBoundary(sdbdef);
    assert(sdiffb != 0);
    uint sdbidx = pSDiffBoundaries.size();
    pSDiffBoundaries.push_back(sdiffb);
    return sdbidx;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_addTet(uint tetidx,
                             steps::mpi::tetopsplit::Comp * comp, double vol,
                             double a1, double a2, double a3, double a4,
                             double d1, double d2, double d3, double d4,
                             int tet0, int tet1, int tet2, int tet3)
{
    steps::solver::Compdef * compdef  = comp->def();
    smtos::Tet * localtet = new smtos::Tet(tetidx, compdef, vol, a1, a2, a3, a4, d1, d2, d3, d4,
                                         tet0, tet1, tet2, tet3, myRank, tetHosts[tetidx]);
    assert(localtet != 0);
    assert(tetidx < pTets.size());
    assert(pTets[tetidx] == 0);
    pTets[tetidx] = localtet;
    comp->addTet(localtet);

    // MPISTEPS
    localtet->setSolver(this);
    // MPISTEPS
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_addWmVol(uint cidx, steps::mpi::tetopsplit::Comp * comp, double vol)
{
    steps::solver::Compdef * compdef  = comp->def();
    smtos::WmVol * localtet = new smtos::WmVol(cidx, compdef, vol, myRank, wmHosts[cidx]);
    assert(localtet != 0);
    assert(cidx < pWmVols.size());
    pWmVols[cidx] = localtet;
    comp->addTet(localtet);
    
    // MPISTEPS
    localtet->setSolver(this);
    // MPISTEPS
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_addTri(uint triidx, steps::mpi::tetopsplit::Patch * patch, double area,
                             double l0, double l1, double l2, double d0, double d1, double d2,  int tinner, int touter, int tri0, int tri1, int tri2)
{
    steps::solver::Patchdef * patchdef = patch->def();
    smtos::Tri * tri = new smtos::Tri(triidx, patchdef, area, l0, l1, l2, d0, d1, d2,  tinner, touter, tri0, tri1, tri2, myRank, triHosts[triidx]);
    assert(tri != 0);
    assert (triidx < pTris.size());
    assert (pTris[triidx] == 0);
    pTris[triidx] = tri;
    patch->addTri(tri);

    // MPISTEPS
    tri->setSolver(this);
    // MPISTEPS
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::reset(void)
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
    
    reacExtent = 0.0;
    diffExtent = 0.0;
    nIteration = 0.0;
    

    
    statedef()->resetTime();
    statedef()->resetNSteps();
	_updateLocal();
    
    compTime = 0.0;
    syncTime = 0.0;
    idleTime = 0.0;

    efieldTime = 0.0;
    rdTime = 0.0;
    dataExchangeTime = 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::run(double endtime)
{
    if (endtime < statedef()->time())
    {
        std::ostringstream os;
        os << "Endtime is before current simulation time";
        throw steps::ArgErr(os.str());
    }
    
    if (recomputeUpdPeriod) _computeUpdPeriod();
    if (efflag()) _runWithEField(endtime);
    else _runWithoutEField(endtime);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_runWithoutEField(double sim_endtime)
{
    MPI_Barrier(MPI_COMM_WORLD);
    
    // This is the time (in seconds) to the next diffusion update. The upper limit,
    // so as to avoid systematic slowing of diffusion, is the inverse of the highest
    // local single-molecule diffusion rate of occupied tets.
    
    // This bool tracks if the update time has been aligned to the endtime. Trying to
    // do this by instead comparing two doubles can cause infinite loops.
    
    // decide globally
    bool aligned = false;
    
    double update_period = updPeriod;
    
    
    MPI_Request* requests = NULL;
    
    // here we assume that all molecule counts have been updated so the rates are accurate
    while (statedef()->time() < sim_endtime and not aligned) {
        #ifdef MPI_DEBUG
        CLOG(DEBUG, "mpi_debug") << "Current Sim Time: " << statedef()->time() << "\n";
        #endif
        
        #ifdef MPI_PROFILING
        double starttime = MPI_Wtime();
        #endif
        
        double pre_ssa_time = statedef()->time();
        // Update period may take us past the endtime- adjust if so
        if (pre_ssa_time + updPeriod > sim_endtime) {
            update_period = sim_endtime - pre_ssa_time;
            aligned=true;
        }
        
        #ifdef MPI_DEBUG
        CLOG(DEBUG, "mpi_debug") << "Global update period: " << update_period << "\n";
        #endif


        // *********************** Operator Split: SSA *********************************
        
        // Run SSA for the update period
        
        double cumulative_dt=0.0;
        
        // Store a sequence of kprocs actually applied
        std::set<smtos::KProc*> applied_ssa_kprocs;
        while(1)
        {
            smtos::KProc * kp = _getNext();
            if (kp == 0) break;
                          
            double a0 = getA0();
            if (a0 == 0.0) break;
            
            double dt=rng()->getExp(a0);
            if (cumulative_dt +dt > update_period) break;
            cumulative_dt += dt;
            

            _executeStep(kp, dt, cumulative_dt);                
            reacExtent +=1;

            
            applied_ssa_kprocs.insert(kp);
            
        }

        // Now for each process advance to jump time and get jump randomly
        // ****************************************************************************
            
        // *********************** Operator Split: Diffusion ***************************

        // Apply diffusion after the update period

        #ifdef MPI_DEBUG
        CLOG(DEBUG, "mpi_debug") << "OP Split: Diffusion\n";
        #endif
        
        #ifdef MPI_PROFILING
        double endtime = MPI_Wtime();
        compTime += (endtime - starttime);
        starttime = MPI_Wtime();
        #endif

        // wait until previous loop finishes sending diffusion data
        if (requests != NULL) {
            MPI_Waitall(nNeighbHosts, requests, MPI_STATUSES_IGNORE);
            delete[] requests;
        }
        
        // create new requests for this loop
        requests = new MPI_Request[nNeighbHosts];
#ifdef MPI_DEBUG
        for (uint i = 0; i < nNeighbHosts; i++) {
            CLOG(DEBUG, "mpi_debug") << "request address" << &(requests[i]) << "\n";
        }
        
#endif
        
        for (auto neighbor : neighbHosts) {
            remoteChanges[neighbor].clear();
        }
        
        #ifdef MPI_PROFILING
        endtime = MPI_Wtime();
        idleTime += (endtime - starttime);
        starttime = MPI_Wtime();
        #endif
        
        // Track how many diffusion 'steps' we do, simply for bookkeeping
        uint nsteps=0;

        
        // to reduce memory cost we use directions to retrieve the update list in upd process
        std::vector<KProc*> applied_diffs;
        std::vector<int> directions;
        #ifdef MPI_DEBUG
        CLOG(DEBUG, "mpi_debug") << "Diffusion for period: " << update_period << "\n";
        #endif
        
        for (uint pos = 0; pos < diffSep; pos++)
        {
            Diff* d = pDiffs[pos];
            double rate = d->crData.rate;
            if (rate == 0) continue;
            // rate is the rate (scaled_dcst * population)
            double scaleddcst = d->getScaledDcst();

            // The number of molecules available for diffusion for this diffusion rule
            double population = rate/scaleddcst;

            // t1, AKA 'X', is a fractional number between 0 and 1: the update period divided
            // by the local mean single-molecule dwellperiod. This fraction gives the mean
            // proportion of molecules to diffuse.
            double t1 = update_period * scaleddcst;
            
            
            if (t1>=1.0) {
                t1=1.0;
            }
            
            // Calculate the occupancy, that is the integrated molecules over the period (units s)
            double occupancy = d->getTet()->getPoolOccupancy(d->getLigLidx()) + population* (update_period - d->getTet()->getLastUpdate(d->getLigLidx()) );
            
            // n is, correctly, a binomial, but the binomial function requires rounding to
            // an integer.
            
            // occupancy/update_period gives the mean number of molecules during the period
            double n_double = occupancy/update_period;

            // could be higher than those available - a source of error
            if (n_double > population) n_double = population;

            double n_int = std::floor(n_double);
            double n_frc = n_double - n_int;
            uint mean_n = static_cast<uint>(n_int);

            // deal linearly with the fraction
            if (n_frc > 0.0)
            {
                double rand01 = rng()->getUnfIE();
                if (rand01 < n_frc) mean_n++;
            }
            
            // Find the binomial n
            uint nmolcs = rng()->getBinom(mean_n, t1);
            
            if (nmolcs == 0) continue;
            
            // we apply here
            if (nmolcs > diffApplyThreshold)
            {
                int direction = d->apply(rng(), nmolcs);
                if (applied_diffs.empty() or applied_diffs.back() != d or directions.back() != direction) {
                    applied_diffs.push_back(d);
                    directions.push_back(direction);
                }
            }
            else
            {
                for (uint ai = 0; ai < nmolcs; ++ai)
                {
                    int direction = d->apply(rng());
                    if (applied_diffs.empty() or applied_diffs.back() != d or directions.back() != direction) {
                        applied_diffs.push_back(d);
                        directions.push_back(direction);
                    }
                }
                
            }
            nsteps += nmolcs;
            diffExtent += nmolcs;
        }
        
        // surface diffusion
        
        for (uint pos = 0; pos < sdiffSep; pos++)
        {
            SDiff* d = pSDiffs[pos];
            double rate = d->crData.rate;

            if (rate == 0) continue;
            // rate is the rate (scaled_dcst * population)
            double scaleddcst = d->getScaledDcst();

            // The number of molecules available for diffusion for this diffusion rule
            double population = rate/scaleddcst;

            // t1, AKA 'X', is a fractional number between 0 and 1: the update period divided
            // by the local mean single-molecule dwellperiod. This fraction gives the mean
            // proportion of molecules to diffuse.
            double t1 = update_period * scaleddcst;
            
            if (t1>=1.0)
            {
                t1=1.0;
            }

            double occupancy = d->getTri()->getPoolOccupancy(d->getLigLidx()) + population* (update_period-  d->getTri()->getLastUpdate(d->getLigLidx()) );

            // n is, correctly, a binomial, but the binomial function requires rounding to
            // an integer.
            
            // occupancy/update_period gives the mean number of molecules during the period
            double n_double = occupancy/update_period;

            // could be higher than those available - a source of error
            if (n_double > population) n_double = population;

            double n_int = std::floor(n_double);
            double n_frc = n_double - n_int;
            uint mean_n = static_cast<uint>(n_int);

            // deal linearly with the fraction
            if (n_frc > 0.0)
            {
                double rand01 = rng()->getUnfIE();
                if (rand01 < n_frc) mean_n++;
            }
            
            // Find the binomial n
            uint nmolcs = rng()->getBinom(mean_n, t1);
            
            if (nmolcs == 0) continue;
            
            // we apply here
            if (nmolcs > diffApplyThreshold)
            {
                int direction = d->apply(rng(), nmolcs);
                if (applied_diffs.empty() or applied_diffs.back() != d or directions.back() != direction) {
                    applied_diffs.push_back(d);
                    directions.push_back(direction);
                }
            }
            else
            {
                for (uint ai = 0; ai < nmolcs; ++ai)
                {
                    int direction = d->apply(rng());
                    if (applied_diffs.empty() or applied_diffs.back() != d or directions.back() != direction) {
                        applied_diffs.push_back(d);
                        directions.push_back(direction);
                    }
                }
            }
            nsteps += nmolcs;
            diffExtent += nmolcs;
        }
        
        #ifdef MPI_PROFILING
        endtime = MPI_Wtime();
        compTime += (endtime - starttime);
        #endif
        
        _remoteSyncAndUpdate(requests, applied_diffs, directions);
        
        // *********************** Operator Split: SSA *********************************
        #ifdef MPI_PROFILING
        starttime = MPI_Wtime();
        #endif
        
        KProcPSetCI akp_end = applied_ssa_kprocs.end();
        for (KProcPSetCI akp = applied_ssa_kprocs.begin(); akp != akp_end; ++akp)
        {
            (*akp)->resetOccupancies();
        }
        
        statedef()->setTime(pre_ssa_time + update_period);
        if (nsteps > 0) statedef()->incNSteps(nsteps);
        
        #ifdef MPI_DEBUG
        CLOG(DEBUG, "mpi_debug") << "OP Split: End of Iteration" << nIteration << "\n";
        #endif
        nIteration += 1;
        
        #ifdef MPI_PROFILING
        endtime = MPI_Wtime();
        compTime += (endtime - starttime);
        #endif
    }
    if (requests != NULL) {
        MPI_Waitall(nNeighbHosts, requests, MPI_STATUSES_IGNORE);
        delete[] requests;
    }
    MPI_Barrier(MPI_COMM_WORLD);

}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_refreshEFTrisV() {
    for (uint tlidx = 0; tlidx < pEFNTris; tlidx++) EFTrisV[tlidx] = pEField->getTriV(tlidx);
    pEFTrisVStale = false;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_runWithEField(double endtime)
{
    #ifdef MPI_DEBUG
    CLOG(DEBUG, "mpi_debug") << "EField computation\n";
    #endif

    #ifdef MPI_PROFILING
    double timing_start = MPI_Wtime();
    #endif

    if (pEFTrisVStale) _refreshEFTrisV();

    #ifdef MPI_PROFILING
    double timing_end = MPI_Wtime();
    efieldTime += (timing_end - timing_start);
    #endif

    while (statedef()->time() < endtime) {

        #ifdef MPI_PROFILING
        timing_start = MPI_Wtime();
        #endif

        double t0 = statedef()->time();
        _runWithoutEField( std::min(t0+pEFDT, endtime));

        #ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        rdTime += (timing_end - timing_start);
        #endif

        #ifdef MPI_PROFILING
        timing_start = MPI_Wtime();
        #endif
        // update host-local currents
        int i_begin = EFTrisI_offset[myRank];
        int i_end = i_begin + EFTrisI_count[myRank];

        double sttime = statedef()->time();
        for (int i = i_begin; i < i_end; ++i) {
            int tlidx = EFTrisI_idx[i];
            EFTrisI_permuted[i] = pEFTris_vec[tlidx]->computeI(EFTrisV[tlidx], sttime-t0, sttime);
        }

        #ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        efieldTime += (timing_end - timing_start);
        #endif

        #ifdef MPI_PROFILING
        timing_start = MPI_Wtime();
        #endif
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                &EFTrisI_permuted[0], &EFTrisI_count[0], &EFTrisI_offset[0], MPI_DOUBLE, MPI_COMM_WORLD);

        #ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        dataExchangeTime += (timing_end - timing_start);
        #endif

        #ifdef MPI_PROFILING
        timing_start = MPI_Wtime();
        #endif

        for (uint i = 0; i < pEFNTris; i++)
                pEField->setTriI(EFTrisI_idx[i], EFTrisI_permuted[i]);

        pEField->advance(sttime-t0);
        _refreshEFTrisV();

        #ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        efieldTime += (timing_end - timing_start);
        #endif
        // TODO: Replace this with something that only resets voltage-dependent things

        #ifdef MPI_PROFILING
        timing_start = MPI_Wtime();
        #endif

        _updateLocal();

        #ifdef MPI_PROFILING
        timing_end = MPI_Wtime();
        rdTime += (timing_end - timing_start);
        #endif
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::advance(double adv)
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

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::step(void)
{
    std::ostringstream os;
    os << "This function is not available for this solver!";
    throw steps::NotImplErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::getTime(void) const
{
    return statedef()->time();
}

////////////////////////////////////////////////////////////////////////

uint smtos::TetOpSplitP::getNSteps(void) const
{
    return statedef()->nsteps();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::setTime(double time)
{
    statedef()->setTime(time);
}

////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::setNSteps(uint nsteps)
{
    statedef()->setNSteps(nsteps);
}

////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::setTemp(double t)
{

    if (efflag() == false)
    {
        std::ostringstream os;
        os << "\nWARNING: Temperature set in simulation without membrane ";
        os << "potential calculation will be ignored.\n";
        //std::cout << os << std::endl;
    }
    assert(t >= 0.0);
    pTemp = t;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getCompVol(uint cidx) const
{
    assert(cidx < statedef()->countComps());
    assert (statedef()->countComps() == pComps.size());
    smtos::Comp * comp = _comp(cidx);
    assert(comp != 0);
    return comp->vol();
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getCompCount(uint cidx, uint sidx) const
{
    assert(cidx < statedef()->countComps());
    assert(sidx < statedef()->countSpecs());
    assert (statedef()->countComps() == pComps.size());
    smtos::Comp * comp = _comp(cidx);
    assert(comp != 0);
    uint slidx = comp->def()->specG2L(sidx);
    if (slidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    uint total_count = 0;
    WmVolPVecCI t_end = comp->endTet();
    for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
    {
        uint svol_count = (*t)->pools()[slidx];
        MPI_Bcast(&svol_count, 1, MPI_UNSIGNED, (*t)->getHost(), MPI_COMM_WORLD);
        total_count += svol_count;
    }

    return total_count;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setCompCount(uint cidx, uint sidx, double n)
{
    MPI_Barrier(MPI_COMM_WORLD);
    assert(cidx < statedef()->countComps());
    assert(sidx < statedef()->countSpecs());
    assert (n >= 0.0);
    assert (statedef()->countComps() == pComps.size());
    smtos::Comp * comp = _comp(cidx);
    uint slidx = comp->def()->specG2L(sidx);

    if (slidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }
    if (n > std::numeric_limits<unsigned int>::max( ))
    {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << std::numeric_limits<unsigned int>::max( ) << ").\n";
        throw steps::ArgErr(os.str());
    }

    // only do the distribution in rank 0
    // then bcast to other ranks
    
    if (myRank == 0) {
        // functions for distribution:
        auto set_count = [slidx](WmVol *tet, uint c) { tet->setCount(slidx, c); };
        auto inc_count = [slidx](WmVol *tet, int c) { tet->incCount(slidx, c, 0.0, true); };
        auto weight = [](WmVol *tet) { return tet->vol(); };

        steps::util::distribute_quantity(n, comp->bgnTet(), comp->endTet(), weight, set_count, inc_count, *rng(), comp->def()->vol());

    }
    
    uint ncomptets = comp->countTets();
    std::vector<uint> counts(ncomptets);
    
    WmVolPVecCI t_end = comp->endTet();
    if (myRank == 0) {
        uint curr_pos = 0;
        for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
        {
            counts[curr_pos] =(*t)->pools()[slidx];
            curr_pos++;
        }
    }

    MPI_Bcast(&(counts.front()), ncomptets, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    
    #ifdef MPI_DEBUG
    CLOG(DEBUG, "mpi_debug") << "distribute comp counts" << counts << "\n";
    #endif

    uint curr_pos = 0;
    for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
    {
        if (myRank != 0) {
            (*t)->setCount(slidx, counts[curr_pos]);
        }
        _updateSpec((*t), sidx);
        curr_pos++;
    }
    _updateSum();
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getCompAmount(uint cidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getCompCount(cidx, sidx);
    return (count / smath::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setCompAmount(uint cidx, uint sidx, double a)
{
    // convert amount in mols to number of molecules
    double a2 = a * steps::math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setCompCount(cidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getCompConc(uint cidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getCompCount(cidx, sidx);
    ssolver::Compdef * comp = statedef()->compdef(cidx);
    assert(comp != 0);
    double vol = comp->vol();
    return count/ (1.0e3 * vol * steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setCompConc(uint cidx, uint sidx, double c)
{
    assert(c >= 0.0);
    assert (cidx < statedef()->countComps());
    ssolver::Compdef * comp = statedef()->compdef(cidx);
    assert(comp != 0);
    double count = c * (1.0e3 * comp->vol() * steps::math::AVOGADRO);
    // the following method does all the necessary argument checking
    _setCompCount(cidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::TetOpSplitP::_getCompClamped(uint cidx, uint sidx) const
{
    assert(cidx < statedef()->countComps());
    assert(sidx < statedef()->countSpecs());
    assert(statedef()->countComps() == pComps.size());
    smtos::Comp * comp = _comp(cidx);
    assert(comp != 0);
    uint lsidx = comp->def()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    WmVolPVecCI t_end = comp->endTet();
    
    bool local_clamped = true;
    for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        if ((*t)->clamped(lsidx) == false) local_clamped = false;
    }
    
    bool global_clamped = false;
    
    MPI_Allreduce(&local_clamped, &global_clamped, 1, MPI_UNSIGNED_SHORT, MPI_LAND, MPI_COMM_WORLD);
    
    return global_clamped;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setCompClamped(uint cidx, uint sidx, bool b)
{
    assert(cidx < statedef()->countComps());
    assert(sidx < statedef()->countSpecs());
    assert(statedef()->countComps() == pComps.size());
    smtos::Comp * comp = _comp(cidx);
    assert(comp != 0);
    uint lsidx = comp->def()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    // Set the flag in def object, though this may not be necessary
    comp->def()->setClamped(lsidx, b);

    WmVolPVecCI t_end = comp->endTet();
    for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        (*t)->setClamped(lsidx, b);
    }
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getCompReacK(uint cidx, uint ridx) const
{
    assert(cidx < statedef()->countComps());
    assert(ridx < statedef()->countReacs());
    assert(statedef()->countComps() == pComps.size());
    smtos::Comp * comp = _comp(cidx);
    assert(comp != 0);
    uint lridx = comp->def()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    // We're just returning the default value for this comp, individual
    // tets may have different Kcsts set individually
    return (comp->def()->kcst(lridx));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setCompReacK(uint cidx, uint ridx, double kf)
{
    assert(cidx < statedef()->countComps());
    assert(ridx < statedef()->countReacs());
    assert(statedef()->countComps() == pComps.size());
    assert(kf >= 0.0);
    smtos::Comp * comp = _comp(cidx);
    assert(comp != 0);
    uint lridx = comp->def()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    // First set the default value for the comp
    comp->def()->setKcst(lridx, kf);

    // Now update all tetrahedra in this comp
    WmVolPVecCI t_end = comp->endTet();
    for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        (*t)->reac(lridx)->setKcst(kf);
    }

    // Rates have changed
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::TetOpSplitP::_getCompReacActive(uint cidx, uint ridx) const
{
    assert(cidx < statedef()->countComps());
    assert(ridx < statedef()->countReacs());
    assert(statedef()->countComps() == pComps.size());
    smtos::Comp * comp = _comp(cidx);
    assert(comp != 0);
    uint lridx = comp->def()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    bool local_active = true;
    
    WmVolPVecCI t_end = comp->endTet();
    for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        if ((*t)->reac(lridx)->inactive() == true) local_active = false;
    }
    
    bool global_active = false;
    
    MPI_Allreduce(&local_active, &global_active, 1, MPI_UNSIGNED_SHORT, MPI_LAND, MPI_COMM_WORLD);
    return global_active;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setCompReacActive(uint cidx, uint ridx, bool a)
{
    assert(cidx < statedef()->countComps());
    assert(ridx < statedef()->countReacs());
    assert(statedef()->countComps() == pComps.size());
    smtos::Comp * comp = _comp(cidx);
    assert(comp != 0);
    uint lridx = comp->def()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    // Set the default value for the comp, though this is not entirely
    // necessary
    comp->def()->setActive(lridx, a);

    WmVolPVecCI t_end = comp->endTet();
    for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        (*t)->reac(lridx)->setActive(a);
    }
    // It's cheaper to just recompute everything.
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getCompDiffD(uint cidx, uint didx) const
{
    assert (cidx < statedef()->countComps());
    assert (didx < statedef()->countDiffs());
    assert (statedef()->countComps() == pComps.size());
    smtos::Comp * comp = _comp(cidx);
    assert (comp != 0);
    uint ldidx = comp->def()->diffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    // We're just returning the default value for this comp, individual
    // tets may have different Dcsts set individually
    return (comp->def()->dcst(ldidx));

}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setCompDiffD(uint cidx, uint didx, double dk)
{
    assert (cidx < statedef()->countComps());
    assert (didx < statedef()->countDiffs());
    assert (statedef()->countComps() == pComps.size());
    assert(dk >= 0.0);
    smtos::Comp * comp = _comp(cidx);
    assert (comp != 0);
    uint ldidx = comp->def()->diffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }
    recomputeUpdPeriod = true;
    // First set the default value for the comp
    comp->def()->setDcst(ldidx, dk);

    // Now update all tets in this comp
    WmVolPVecCI t_end = comp->endTet();
    for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        if (smtos::Tet * tet = dynamic_cast<smtos::Tet *>(*t))
        {
            tet->diff(ldidx)->setDcst(dk);
        }
        else
        {
            std::ostringstream os;
            os << "Cannot change diffusion constant in well-mixed compartment.";
            throw steps::ArgErr(os.str());
        }
    }

    // Rates have changed
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::TetOpSplitP::_getCompDiffActive(uint cidx, uint didx) const
{
	assert (cidx < statedef()->countComps());
	assert (didx < statedef()->countDiffs());
	assert (statedef()->countComps() == pComps.size());
	smtos::Comp * comp = _comp(cidx);
	assert (comp != 0);
	uint ldidx = comp->def()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Diffusion rule undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

    bool local_active = true;
    
	WmVolPVecCI t_end = comp->endTet();
	for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
		//smtos::WmVol * wmcomp = (*t);
		if (smtos::Tet * tet = dynamic_cast<smtos::Tet *>((*t)))
		{
			if (tet->diff(ldidx)->inactive() == true) local_active = false;
		}
		else
		{
			std::ostringstream os;
			os << "Diffusion activation not defined in well-mixed compartment.\n";
			throw steps::ArgErr(os.str());
		}
	}
    bool global_active = false;
    MPI_Allreduce(&local_active, &global_active, 1, MPI_UNSIGNED_SHORT, MPI_LAND, MPI_COMM_WORLD);
	return global_active;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setCompDiffActive(uint cidx, uint didx, bool act)
{
    assert (cidx < statedef()->countComps());
    assert (didx < statedef()->countDiffs());
    assert (statedef()->countComps() == pComps.size());
    smtos::Comp * comp = _comp(cidx);
    assert (comp != 0);
    uint ldidx = comp->def()->diffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    WmVolPVecCI t_end = comp->endTet();
    for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        //smtos::WmVol * wmcomp = (*t);
        if (smtos::Tet * tet = dynamic_cast<smtos::Tet *>((*t)))
        {
            tet->diff(ldidx)->setActive(act);
        }
        else
        {
            std::ostringstream os;
            os << "Cannot change diffusion constant in well-mixed compartment.\n";
            throw steps::ArgErr(os.str());
        }
    }
    
    recomputeUpdPeriod = true;

    // It's cheaper to just recompute everything.
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getPatchArea(uint pidx) const
{
    assert (pidx < statedef()->countPatches());
    assert (statedef()->countPatches() == pPatches.size());
    smtos::Patch * patch = _patch(pidx);
    assert (patch != 0);
    return patch->area();
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getPatchCount(uint pidx, uint sidx) const
{
    assert (pidx < statedef()->countPatches());
    assert (sidx < statedef()->countSpecs());
    assert (statedef()->countPatches() == pPatches.size());
    smtos::Patch * patch = _patch(pidx);
    assert (patch != 0);
    uint slidx = patch->def()->specG2L(sidx);
    if (slidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }

    uint total_count = 0;
    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        uint svol_count = (*t)->pools()[slidx];
        MPI_Bcast(&svol_count, 1, MPI_UNSIGNED, (*t)->getHost(), MPI_COMM_WORLD);
        total_count += svol_count;
    }
    
    return total_count;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setPatchCount(uint pidx, uint sidx, double n)
{
    MPI_Barrier(MPI_COMM_WORLD);
    assert (pidx < statedef()->countPatches());
	assert (sidx < statedef()->countSpecs());
	assert (statedef()->countPatches() == pPatches.size());
	assert (n >= 0.0);
	smtos::Patch * patch = _patch(pidx);
	assert (patch != 0);
	uint slidx = patch->def()->specG2L(sidx);
    
    if (slidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}
	if (n > std::numeric_limits<unsigned int>::max( ))
	{
		std::ostringstream os;
		os << "Can't set count greater than maximum unsigned integer (";
		os << std::numeric_limits<unsigned int>::max( ) << ").\n";
		throw steps::ArgErr(os.str());
	}

    // only do the distribution in rank 0
    // then bcast to other ranks
    
    if (myRank == 0) {
        // functions for distribution:
        auto set_count = [slidx](Tri *tri, uint c) { tri->setCount(slidx, c); };
        auto inc_count = [slidx](Tri *tri, int c) { tri->incCount(slidx, c, 0.0, true); };
        auto weight = [](Tri *tri) { return tri->area(); };

        steps::util::distribute_quantity(n, patch->bgnTri(), patch->endTri(), weight, set_count, inc_count, *rng(), patch->def()->area());
    }
    uint npatchtris = patch->countTris();
    std::vector<uint> counts(npatchtris);
    TriPVecCI t_end = patch->endTri();
    if (myRank == 0) {
        uint curr_pos = 0;

        for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
        {
            counts[curr_pos] =(*t)->pools()[slidx];
            curr_pos++;
        }

    }
    
    MPI_Bcast(&(counts.front()), npatchtris, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    
    #ifdef MPI_DEBUG
    CLOG(DEBUG, "mpi_debug") << "distribute patch counts " << counts << "\n";
    #endif
    
    uint curr_pos = 0;
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        if (myRank != 0) {
            // set count only don't need sync
            (*t)->setCount(slidx, counts[curr_pos]);
        }
        _updateSpec((*t), sidx);
        curr_pos++;
    }
    _updateSum();
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getPatchAmount(uint pidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getPatchCount(pidx, sidx);
    return (count / steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setPatchAmount(uint pidx, uint sidx, double a)
{
    assert(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * steps::math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setPatchCount(pidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::TetOpSplitP::_getPatchClamped(uint pidx, uint sidx) const
{
    assert (pidx < statedef()->countPatches());
    assert (sidx < statedef()->countSpecs());
    assert(statedef()->countPatches() == pPatches.size());
    smtos::Patch * patch = _patch(pidx);
    assert(patch != 0);
    uint lsidx = patch->def()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }

    bool local_clamped = true;
    
    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        if ((*t)->clamped(lsidx) == false) local_clamped = false;
    }
    bool global_clamped = false;
    
    MPI_Allreduce(&local_clamped, &global_clamped, 1, MPI_UNSIGNED_SHORT, MPI_LAND, MPI_COMM_WORLD);
    
    return global_clamped;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setPatchClamped(uint pidx, uint sidx, bool buf)
{
    assert (pidx < statedef()->countPatches());
    assert (sidx < statedef()->countSpecs());
    assert(statedef()->countPatches() == pPatches.size());
    smtos::Patch * patch = _patch(pidx);
    assert(patch != 0);
    uint lsidx = patch->def()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }

    // Set the flag in def object for consistency, though this is not
    // entirely necessary
    patch->def()->setClamped(lsidx, buf);

    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        (*t)->setClamped(lsidx, buf);
    }
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getPatchSReacK(uint pidx, uint ridx) const
{
    assert (pidx < statedef()->countPatches());
    assert (ridx < statedef()->countSReacs());
    assert(statedef()->countPatches() == pPatches.size());
    smtos::Patch * patch = _patch(pidx);
    assert(patch != 0);
    uint lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }

    // We're just returning the default value for this patch, individual
    // triangles may have different Kcsts set
    return (patch->def()->kcst(lsridx));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setPatchSReacK(uint pidx, uint ridx, double kf)
{
    assert (pidx < statedef()->countPatches());
    assert (ridx < statedef()->countSReacs());
    assert(statedef()->countPatches() == pPatches.size());
    assert(kf >= 0.0);
    smtos::Patch * patch = _patch(pidx);
    assert(patch != 0);
    uint lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }

    // First set the default values for this patch
    patch->def()->setKcst(lsridx, kf);

    // Now update all triangles in this patch
    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        (*t)->sreac(lsridx)->setKcst(kf);
    }

    // Rates have changed
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::TetOpSplitP::_getPatchSReacActive(uint pidx, uint ridx) const
{
    assert (pidx < statedef()->countPatches());
    assert (ridx < statedef()->countSReacs());
    assert(statedef()->countPatches() == pPatches.size());
    smtos::Patch * patch = _patch(pidx);
    assert(patch != 0);
    uint lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }
    bool local_active = true;
    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        if ((*t)->sreac(lsridx)->inactive() == true) local_active = false;
    }
    bool global_active = false;
    MPI_Allreduce(&local_active, &global_active, 1, MPI_UNSIGNED_SHORT, MPI_LAND, MPI_COMM_WORLD);
    return global_active;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setDiffBoundaryDiffusionActive(uint dbidx, uint sidx, bool act)
{
    assert (dbidx < statedef()->countDiffBoundaries());
    assert (sidx < statedef()->countSpecs());

    // Need to do two things:
    // 1) check if the species is defined in both compartments conencted
    // by the diffusion boundary
    // 2) loop over all tetrahedrons around the diff boundary and then the
    // diffusion rules and activate diffusion if the diffusion rule
    // relates to this species

    smtos::DiffBoundary * diffb = _diffboundary(dbidx);
    smtos::Comp * compA = diffb->compA();
    smtos::Comp * compB = diffb->compB();

    /*
       ssolver::Diffdef * diffdef = statedef()->diffdef(didx);
       uint specgidx = diffdef->lig();
       */
    uint lsidxA = compA->def()->specG2L(sidx);
    uint lsidxB = compB->def()->specG2L(sidx);


    if (lsidxA == ssolver::LIDX_UNDEFINED or lsidxB == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartments connected by diffusion boundary.\n";
        throw steps::ArgErr(os.str());
    }

    std::vector<uint> bdtets = diffb->getTets();
    std::vector<uint> bdtetsdir = diffb->getTetDirection();

    // Have to use indices rather than iterator because need access to the
    // tet direction
    uint ntets = bdtets.size();

    for (uint bdt = 0; bdt != ntets; ++bdt)
    {
        smtos::Tet * tet = _tet(bdtets[bdt]);
        if (!tet->getInHost()) continue;
        uint direction = bdtetsdir[bdt];
        assert(direction >= 0 and direction < 4);

        // Each diff kproc then has access to the species through it's defined parent
        uint ndiffs = tet->compdef()->countDiffs();
        for (uint d = 0; d != ndiffs; ++d)
        {
            smtos::Diff * diff = tet->diff(d);
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

bool smtos::TetOpSplitP::_getDiffBoundaryDiffusionActive(uint dbidx, uint sidx) const
{
    assert (dbidx < statedef()->countDiffBoundaries());
    assert (sidx < statedef()->countSpecs());

    smtos::DiffBoundary * diffb = _diffboundary(dbidx);
    smtos::Comp * compA = diffb->compA();
    smtos::Comp * compB = diffb->compB();

    uint lsidxA = compA->def()->specG2L(sidx);
    uint lsidxB = compB->def()->specG2L(sidx);

    if (lsidxA == ssolver::LIDX_UNDEFINED or lsidxB == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartments connected by diffusion boundary.\n";
        throw steps::ArgErr(os.str());
    }

    std::vector<uint> bdtets = diffb->getTets();
    std::vector<uint> bdtetsdir = diffb->getTetDirection();

    // Have to use indices rather than iterator because need access to the
    // tet direction
    
    bool local_active = true;
    
    uint ntets = bdtets.size();
    
    for (uint bdt = 0; bdt != ntets; ++bdt)
    {
        smtos::Tet * tet = _tet(bdtets[bdt]);
        if (!tet->getInHost()) continue;
        uint direction = bdtetsdir[bdt];
        assert(direction >= 0 and direction < 4);

        // Each diff kproc then has access to the species through it's defined parent
        uint ndiffs = tet->compdef()->countDiffs();
        for (uint d = 0; d != ndiffs; ++d)
        {
            smtos::Diff * diff = tet->diff(d);
            // sidx is the global species index; so is the lig() return from diffdef
            uint specgidx = diff->def()->lig();
            if (specgidx == sidx)
            {
                // Just need to check the first one
                if (diff->getDiffBndActive(direction)) local_active = true;
                else local_active = false;
                break;
            }
        }
    }
    bool global_active = false;
    MPI_Allreduce(&local_active, &global_active, 1, MPI_UNSIGNED_SHORT, MPI_LAND, MPI_COMM_WORLD);
    return global_active;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setDiffBoundaryDcst(uint dbidx, uint sidx, double dcst, uint direction_comp)
{
    assert (dbidx < statedef()->countDiffBoundaries());
    assert (sidx < statedef()->countSpecs());

    smtos::DiffBoundary * diffb = _diffboundary(dbidx);
    smtos::Comp * compA = diffb->compA();
    smtos::Comp * compB = diffb->compB();
    
    uint lsidxA = compA->def()->specG2L(sidx);
    uint lsidxB = compB->def()->specG2L(sidx);
    
    if (lsidxA == ssolver::LIDX_UNDEFINED or lsidxB == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartments connected by diffusion boundary.\n";
        throw steps::ArgErr(os.str());
    }
    recomputeUpdPeriod = true;
    
    steps::solver::Compdef * dirc_compdef = NULL;
    if (direction_comp != std::numeric_limits<uint>::max()) {
        dirc_compdef = _comp(direction_comp)->def();
    }
    
    std::vector<uint> bdtets = diffb->getTets();
    std::vector<uint> bdtetsdir = diffb->getTetDirection();
    
    uint ntets = bdtets.size();
    
    for (uint bdt = 0; bdt != ntets; ++bdt)
    {
        smtos::Tet * tet = _tet(bdtets[bdt]);
        if (!tet->getInHost()) continue;
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
            smtos::Diff * diff = tet->diff(d);
            // sidx is the global species index; so is the lig() return from diffdef
            uint specgidx = diff->def()->lig();
            if (specgidx == sidx)
            {
                #ifdef DIRECTIONAL_DCST_DEBUG
                CLOG(DEBUG, "steps_debug") << "direction: " << direction << " dcst: " << dcst << "\n";
                #endif
                
                // The following function will automatically activate diffusion
                // in this direction if necessary
                diff->setDirectionDcst(direction, dcst);
                _updateElement(diff);
            }
        }
    }
    
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setPatchSReacActive(uint pidx, uint ridx, bool a)
{
    assert (pidx < statedef()->countPatches());
    assert (ridx < statedef()->countSReacs());
    assert(statedef()->countPatches() == pPatches.size());
    smtos::Patch * patch = _patch(pidx);
    assert(patch != 0);
    uint lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }

    // First set the flags in def object for consistency, though this is
    // not entirely necessary for this solver
    patch->def()->setActive(lsridx, a);

    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        (*t)->sreac(lsridx)->setActive(a);
    }
    // It's cheaper to just recompute everything.
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::TetOpSplitP::_getPatchVDepSReacActive(uint pidx, uint vsridx) const
{
    assert (pidx < statedef()->countPatches());
    assert (vsridx < statedef()->countVDepSReacs());
    assert(statedef()->countPatches() == pPatches.size());
    smtos::Patch * patch = _patch(pidx);
    assert(patch != 0);
    uint lvsridx = patch->def()->vdepsreacG2L(vsridx);
    if (lvsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Voltage-dependent surface reaction undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }
    
    bool local_active = true;

    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        if ((*t)->vdepsreac(lvsridx)->inactive() == true)  local_active = false;
    }
    bool global_active = false;
    MPI_Allreduce(&local_active, &global_active, 1, MPI_UNSIGNED_SHORT, MPI_LAND, MPI_COMM_WORLD);
    return global_active;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setPatchVDepSReacActive(uint pidx, uint vsridx, bool a)
{
    assert (pidx < statedef()->countPatches());
    assert (vsridx < statedef()->countVDepSReacs());
    assert(statedef()->countPatches() == pPatches.size());
    smtos::Patch * patch = _patch(pidx);
    assert(patch != 0);
    uint lvsridx = patch->def()->vdepsreacG2L(vsridx);
    if (lvsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Voltage-dependent surface reaction undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }

    // Not necessary and not possible to set the flags in def object

    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        (*t)->vdepsreac(lvsridx)->setActive(a);
    }
    // It's cheaper to just recompute everything.
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

// WEILIANG: check the following 3 surface diffusion functions please!
void smtos::TetOpSplitP::_setSDiffBoundaryDiffusionActive(uint sdbidx, uint sidx, bool act)
{
    // Need to do two things:
    // 1) check if the species is defined in both patches connected
    // by the surface diffusion boundary
    // 2) loop over all triangles around the sdiff boundary and then the
    // diffusion rules and activate diffusion if the diffusion rule
    // relates to this species

	smtos::SDiffBoundary * sdiffb = _sdiffboundary(sdbidx);
	smtos::Patch * patchA = sdiffb->patchA();
	smtos::Patch * patchB = sdiffb->patchB();

    uint lsidxA = patchA->def()->specG2L(sidx);
    uint lsidxB = patchB->def()->specG2L(sidx);

    if (lsidxA == ssolver::LIDX_UNDEFINED or lsidxB == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in patches connected by surface diffusion boundary.\n";
        throw steps::ArgErr(os.str());
    }

    std::vector<uint> sbdtris = sdiffb->getTris();
    std::vector<uint> sbdtrisdir = sdiffb->getTriDirection();

    // Have to use indices rather than iterator because need access to the
    // tri direction
    uint ntris = sbdtris.size();

    for (uint sbdt = 0; sbdt != ntris; ++sbdt)
    {
    	smtos::Tri * tri = _tri(sbdtris[sbdt]);
        if (!tri->getInHost()) continue;
        uint direction = sbdtrisdir[sbdt];
        assert(direction >= 0 and direction < 3);

        // Each sdiff kproc then has access to the species through its defined parent
        uint nsdiffs = tri->patchdef()->countSurfDiffs();
        for (uint sd = 0; sd != nsdiffs; ++sd)
        {
        	smtos::SDiff * sdiff = tri->sdiff(sd);
            // sidx is the global species index; so is the lig() return from diffdef
            uint specgidx = sdiff->def()->lig();
            if (specgidx == sidx)
            {
                sdiff->setSDiffBndActive(direction, act);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::TetOpSplitP::_getSDiffBoundaryDiffusionActive(uint sdbidx, uint sidx) const
{
	smtos::SDiffBoundary * sdiffb = _sdiffboundary(sdbidx);
	smtos::Patch * patchA = sdiffb->patchA();
	smtos::Patch * patchB = sdiffb->patchB();

    uint lsidxA = patchA->def()->specG2L(sidx);
    uint lsidxB = patchB->def()->specG2L(sidx);

    if (lsidxA == ssolver::LIDX_UNDEFINED or lsidxB == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in patches connected by surface diffusion boundary.\n";
        throw steps::ArgErr(os.str());
    }

    std::vector<uint> sbdtris = sdiffb->getTris();
    std::vector<uint> sbdtrisdir = sdiffb->getTriDirection();

    // Have to use indices rather than iterator because need access to the
    // tet direction

    short local_active = 1; // true

    uint ntris = sbdtris.size();

    for (uint sbdt = 0; sbdt != ntris; ++sbdt)
    {
    	smtos::Tri * tri = _tri(sbdtris[sbdt]);
        if (!tri->getInHost()) continue;
        uint direction = sbdtrisdir[sbdt];
        assert(direction >= 0 and direction < 3);

        // Each sdiff kproc then has access to the species through its defined parent
        uint nsdiffs = tri->patchdef()->countSurfDiffs();
        for (uint sd = 0; sd != nsdiffs; ++sd)
        {
        	smtos::SDiff * sdiff = tri->sdiff(sd);
            // sidx is the global species index; so is the lig() return from diffdef
            uint specgidx = sdiff->def()->lig();
            if (specgidx == sidx)
            {
                // Just need to check the first one
                if (sdiff->getSDiffBndActive(direction)) local_active = 1;
                else local_active = 0;
                break;
            }
        }
    }
    short global_active = 0;
    MPI_Allreduce(&local_active, &global_active, 1, MPI_SHORT, MPI_LAND, MPI_COMM_WORLD);
    return global_active;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setSDiffBoundaryDcst(uint sdbidx, uint sidx, double dcst, uint direction_patch)
{
	smtos::SDiffBoundary * sdiffb = _sdiffboundary(sdbidx);
	smtos::Patch * patchA = sdiffb->patchA();
	smtos::Patch * patchB = sdiffb->patchB();

    uint lsidxA = patchA->def()->specG2L(sidx);
    uint lsidxB = patchB->def()->specG2L(sidx);

    if (lsidxA == ssolver::LIDX_UNDEFINED or lsidxB == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in patches connected by surface diffusion boundary.\n";
        throw steps::ArgErr(os.str());
    }

    recomputeUpdPeriod = true;

    steps::solver::Patchdef * dirp_patchdef = NULL;
    if (direction_patch != std::numeric_limits<uint>::max()) {
    	dirp_patchdef = _patch(direction_patch)->def();
    }

    std::vector<uint> sbdtris = sdiffb->getTris();
    std::vector<uint> sbdtrisdir = sdiffb->getTriDirection();

    uint ntris = sbdtris.size();

    for (uint sbdt = 0; sbdt != ntris; ++sbdt)
    {
    	smtos::Tri * tri = _tri(sbdtris[sbdt]);

        if (!tri->getInHost()) continue;

        if (dirp_patchdef == tri->patchdef()) {
            continue;
        }
        uint direction = sbdtrisdir[sbdt];
        assert(direction >= 0 and direction < 3);

        // Each diff kproc then has access to the species through its defined parent
        uint nsdiffs = tri->patchdef()->countSurfDiffs();
        for (uint sd = 0; sd != nsdiffs; ++sd)
        {
        	smtos::SDiff * sdiff = tri->sdiff(sd);
            // sidx is the global species index; so is the lig() return from diffdef
            uint specgidx = sdiff->def()->lig();
            if (specgidx == sidx)
            {
                #ifdef DIRECTIONAL_DCST_DEBUG
                CLOG(DEBUG, "steps_debug") << "direction: " << direction << " dcst: " << dcst << "\n";
                #endif

                // The following function will automatically activate diffusion
                // in this direction if necessary
                sdiff->setDirectionDcst(direction, dcst);
                _updateElement(sdiff);
            }
        }
    }

    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

uint smtos::TetOpSplitP::addKProc(steps::mpi::tetopsplit::KProc * kp, int host)
{
    SchedIDX nidx = pKProcs.size();
    pKProcs.push_back(kp);
    return nidx;
}

////////////////////////////////////////////////////////////////////////////////

steps::mpi::tetopsplit::KProc * smtos::TetOpSplitP::_getNext(void) const
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

void smtos::TetOpSplitP::_executeStep(steps::mpi::tetopsplit::KProc * kp, double dt, double period)
{
    kp->apply(rng(), dt, statedef()->time(), period);
    statedef()->incTime(dt);
    
    // as in 0.6.1 reaction and surface reaction only require updates of local
    // KProcs, it may change if VDepSurface reaction is added in the future
    std::vector<smtos::KProc*> upd = kp->getLocalUpdVec();
    _updateLocal(upd);
    statedef()->incNSteps(1);

}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_updateSpec(steps::mpi::tetopsplit::WmVol * tet, uint spec_gidx)
{
    // NOTE: this function does not update the Sum of popensity, _updateSum() is required after calling it.
    
    #ifdef MPI_DEBUG
    //CLOG(DEBUG, "mpi_debug") << "Update kps related to species " << spec_gidx << " of subvolume " << tet->idx() <<".\n";
    #endif
    if (!tet->getInHost()) return;
    
    std::set<smtos::KProc*> updset;

    // Loop over tet.
    uint nkprocs = tet->countKProcs();
    uint startKProcIdx = tet->getStartKProcIdx();
    
    for (uint k = 0; k < nkprocs; k++)
    {
        if (tet->KProcDepSpecTet(k, tet, spec_gidx)) updset.insert(tet->getKProc(k));
    }
    
    std::vector<smtos::Tri *>::const_iterator tri_end = tet->nexttriEnd();
    for (std::vector<smtos::Tri *>::const_iterator tri = tet->nexttriBegin();
         tri != tri_end; ++tri)
    {
        if ((*tri) == 0) continue;
        nkprocs = (*tri)->countKProcs();
        startKProcIdx = (*tri)->getStartKProcIdx();
        for (uint sk = 0; sk < nkprocs; sk++) {
            if ((*tri)->KProcDepSpecTet(sk, tet, spec_gidx)) updset.insert((*tri)->getKProc(sk));
        }
    }

    #ifdef MPI_DEBUG
    //CLOG(DEBUG, "mpi_debug") << "length of update kp " << updvec.size() <<".\n";
    #endif
    for (auto & kp : updset) {
        _updateElement(kp);
    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_updateSpec(steps::mpi::tetopsplit::Tri * tri, uint spec_gidx)
{
    // NOTE: this function does not update the Sum of popensity, _updateSum() is required after calling it.
    
    #ifdef MPI_DEBUG
    //CLOG(DEBUG, "mpi_debug") << "Update kps related to species " << spec_gidx << " of tri " << tri->idx() << ".\n";
    #endif
    if (!tri->getInHost()) return;
    std::set<smtos::KProc*> updset;

    uint nkprocs = tri->countKProcs();
    uint startKProcIdx = tri->getStartKProcIdx();
    
    for (uint sk = 0; sk < nkprocs; sk++)
    {
        if (tri->KProcDepSpecTri(sk, tri, spec_gidx)) updset.insert(tri->getKProc(sk));
    }
    for (auto & kp : updset) {
        _updateElement(kp);
    }

}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getCompReacH(uint cidx, uint ridx) const
{
    assert(cidx < statedef()->countComps());
    assert(ridx < statedef()->countReacs());
    ssolver::Compdef * comp = statedef()->compdef(cidx);
    assert(comp != 0);
    uint lridx = comp->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    smtos::Comp * lcomp = pComps[cidx];
    assert (lcomp->def() == comp);

    WmVolPVecCI t_bgn = lcomp->bgnTet();
    WmVolPVecCI t_end = lcomp->endTet();
    if (t_bgn == t_end) return 0.0;

    double local_h = 0.0;
    for (WmVolPVecCI t = t_bgn; t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        smtos::Reac * reac = (*t)->reac(lridx);
        local_h += reac->h();
    }
    double global_h = 0.0;
    MPI_Allreduce(&local_h, &global_h, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_h;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getCompReacC(uint cidx, uint ridx) const
{
    assert(cidx < statedef()->countComps());
    assert(ridx < statedef()->countReacs());
    ssolver::Compdef * comp = statedef()->compdef(cidx);
    assert(comp != 0);
    uint lridx = comp->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    smtos::Comp * lcomp = pComps[cidx];
    assert (lcomp->def() == comp);

    WmVolPVecCI t_bgn = lcomp->bgnTet();
    WmVolPVecCI t_end = lcomp->endTet();
    if (t_bgn == t_end) return 0.0;
    double local_c = 0.0;
    double local_v = 0.0;
    for (WmVolPVecCI t = t_bgn; t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        double v = (*t)->vol();
        smtos::Reac * reac = (*t)->reac(lridx);
        local_c += (reac->c() * v);
        local_v += v;
    }
    double global_c = 0.0;
    double global_v = 0.0;
    MPI_Allreduce(&local_c, &global_c, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_v, &global_v, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_c/global_v;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getCompReacA(uint cidx, uint ridx) const
{
    assert(cidx < statedef()->countComps());
    assert(ridx < statedef()->countReacs());
    ssolver::Compdef * comp = statedef()->compdef(cidx);
    assert(comp != 0);
    uint lridx = comp->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    smtos::Comp * lcomp = pComps[cidx];
    assert (lcomp->def() == comp);

    WmVolPVecCI t_bgn = lcomp->bgnTet();
    WmVolPVecCI t_end = lcomp->endTet();
    if (t_bgn == t_end) return 0.0;

    double local_a = 0.0;
    for (WmVolPVecCI t = t_bgn; t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        smtos::Reac * reac = (*t)->reac(lridx);
        local_a += reac->rate();
    }
    double global_a = 0.0;
    MPI_Allreduce(&local_a, &global_a, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_a;
}

////////////////////////////////////////////////////////////////////////////////

uint smtos::TetOpSplitP::_getCompReacExtent(uint cidx, uint ridx) const
{
    assert(cidx < statedef()->countComps());
    assert(ridx < statedef()->countReacs());
    ssolver::Compdef * comp = statedef()->compdef(cidx);
    assert(comp != 0);
    uint lridx = comp->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    smtos::Comp * lcomp = pComps[cidx];
    assert (lcomp->def() == comp);

    WmVolPVecCI t_bgn = lcomp->bgnTet();
    WmVolPVecCI t_end = lcomp->endTet();
    if (t_bgn == t_end) return 0.0;

    uint local_x = 0.0;
    for (WmVolPVecCI t = t_bgn; t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        smtos::Reac * reac = (*t)->reac(lridx);
        local_x += reac->getExtent();
    }

    double global_x = 0.0;
    MPI_Allreduce(&local_x, &global_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_x;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_resetCompReacExtent(uint cidx, uint ridx)
{
    assert(cidx < statedef()->countComps());
    assert(ridx < statedef()->countReacs());
    ssolver::Compdef * comp = statedef()->compdef(cidx);
    assert(comp != 0);
    uint lridx = comp->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    // The 'local' Comp object has same index as solver::Compdef object
    smtos::Comp * lcomp = pComps[cidx];
    assert (lcomp->def() == comp);

    WmVolPVecCI t_bgn = lcomp->bgnTet();
    WmVolPVecCI t_end = lcomp->endTet();
    if (t_bgn == t_end) return;

    for (WmVolPVecCI t = t_bgn; t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        smtos::Reac * reac = (*t)->reac(lridx);
        reac->resetExtent();
    }
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getPatchSReacH(uint pidx, uint ridx) const
{
    assert (pidx < statedef()->countPatches());
    assert (ridx < statedef()->countSReacs());
    ssolver::Patchdef * patch = statedef()->patchdef(pidx);
    uint lsridx = patch->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    smtos::Patch * lpatch = pPatches[pidx];
    assert (lpatch->def() == patch);

    TriPVecCI t_bgn = lpatch->bgnTri();
    TriPVecCI t_end = lpatch->endTri();
    if (t_bgn == t_end) return 0.0;

    double local_h = 0.0;
    for (TriPVecCI t = t_bgn; t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        smtos::SReac * sreac = (*t)->sreac(lsridx);
        local_h += sreac->h();
    }

    double global_h = 0.0;
    MPI_Allreduce(&local_h, &global_h, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_h;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getPatchSReacC(uint pidx, uint ridx) const
{
    assert (pidx < statedef()->countPatches());
    assert (ridx < statedef()->countSReacs());
    ssolver::Patchdef * patch = statedef()->patchdef(pidx);
    uint lsridx = patch->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    smtos::Patch * lpatch = pPatches[pidx];
    assert (lpatch->def() == patch);

    TriPVecCI t_bgn = lpatch->bgnTri();
    TriPVecCI t_end = lpatch->endTri();
    if (t_bgn == t_end) return 0.0;

    double local_c = 0.0;
    double local_a = 0.0;
    for (TriPVecCI t = t_bgn; t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        double a = (*t)->area();
        smtos::SReac * sreac = (*t)->sreac(lsridx);
        local_c += (sreac->c() * a);
        local_a += a;
    }
    double global_c = 0.0;
    double global_a = 0.0;
    MPI_Allreduce(&local_c, &global_c, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_a, &global_a, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_c/global_a;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getPatchSReacA(uint pidx, uint ridx) const
{
    assert (pidx < statedef()->countPatches());
    assert (ridx < statedef()->countSReacs());
    ssolver::Patchdef * patch = statedef()->patchdef(pidx);
    uint lsridx = patch->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    smtos::Patch * lpatch = pPatches[pidx];
    assert (lpatch->def() == patch);

    TriPVecCI t_bgn = lpatch->bgnTri();
    TriPVecCI t_end = lpatch->endTri();
    if (t_bgn == t_end) return 0.0;

    double local_a = 0.0;
    for (TriPVecCI t = t_bgn; t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        smtos::SReac * sreac = (*t)->sreac(lsridx);
        local_a += sreac->rate();
    }
    double global_a = 0.0;
    MPI_Allreduce(&local_a, &global_a, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_a;
}

////////////////////////////////////////////////////////////////////////////////

uint smtos::TetOpSplitP::_getPatchSReacExtent(uint pidx, uint ridx) const
{
    assert (pidx < statedef()->countPatches());
    assert (ridx < statedef()->countSReacs());
    ssolver::Patchdef * patch = statedef()->patchdef(pidx);
    uint lsridx = patch->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    smtos::Patch * lpatch = pPatches[pidx];
    assert (lpatch->def() == patch);

    TriPVecCI t_bgn = lpatch->bgnTri();
    TriPVecCI t_end = lpatch->endTri();
    if (t_bgn == t_end) return 0;

    double local_x = 0.0;
    for (TriPVecCI t = t_bgn; t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        smtos::SReac * sreac = (*t)->sreac(lsridx);
        local_x += sreac->getExtent();
    }
    double global_x = 0.0;
    MPI_Allreduce(&local_x, &global_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_x;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_resetPatchSReacExtent(uint pidx, uint ridx)
{
    assert (pidx < statedef()->countPatches());
    assert (ridx < statedef()->countSReacs());
    ssolver::Patchdef * patch = statedef()->patchdef(pidx);
    uint lsridx = patch->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }

    // The 'local' Patch object has same index as solver::Patchdef object
    smtos::Patch * lpatch = pPatches[pidx];
    assert (lpatch->def() == patch);

    TriPVecCI t_bgn = lpatch->bgnTri();
    TriPVecCI t_end = lpatch->endTri();
    if (t_bgn == t_end) return;

    for (TriPVecCI t = t_bgn; t != t_end; ++t)
    {
        if (!(*t)->getInHost()) continue;
        smtos::SReac * sreac = (*t)->sreac(lsridx);
        sreac->resetExtent();
    }
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTetVol(uint tidx) const
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

void smtos::TetOpSplitP::_setTetVol(uint tidx, double vol)
{
    std::ostringstream os;
    os << "Can not change tetrahedron volume in a mesh based solver.\n";
    throw steps::NotImplErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::TetOpSplitP::_getTetSpecDefined(uint tidx, uint sidx) const
{
    assert (tidx < pTets.size());
    assert (sidx < statedef()->countSpecs());

    if (pTets[tidx] == 0) return false;

    smtos::Tet * tet = pTets[tidx];
    uint lsidx = tet->compdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return false;
    else return true;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTetCount(uint tidx, uint sidx) const
{
    MPI_Barrier(MPI_COMM_WORLD);
    #ifdef MPI_DEBUG
    //CLOG(DEBUG, "mpi_debug") << "Get Tet " << tidx << " count of spec " << sidx << "\n";
    #endif
    assert (tidx < pTets.size());
    assert (sidx < statedef()->countSpecs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    smtos::Tet * tet = pTets[tidx];
    uint lsidx = tet->compdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    uint count = tet->pools()[lsidx];
    MPI_Bcast(&count, 1, MPI_UNSIGNED, tetHosts[tidx], MPI_COMM_WORLD);
    #ifdef MPI_DEBUG
    //CLOG(DEBUG, "mpi_debug") << "count: " << count << "\n";
    #endif
	MPI_Barrier(MPI_COMM_WORLD);
    return count;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setTetCount(uint tidx, uint sidx, double n)
{
    #ifdef MPI_DEBUG
    //CLOG(DEBUG, "mpi_debug") << "Set Tet " << tidx << " count of spec " << sidx << ": " << n << "\n";
    #endif
    
    MPI_Barrier(MPI_COMM_WORLD);
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
    
    smtos::Tet * tet = pTets[tidx];
    
    uint lsidx = tet->compdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    uint count = 0;
    
    if (myRank == 0) {

        double n_int = std::floor(n);
        double n_frc = n - n_int;
        count = static_cast<uint>(n_int);
        if (n_frc > 0.0)
        {
            double rand01 = rng()->getUnfIE();
            if (rand01 < n_frc) count++;
        }
    }
    
    // Tet object updates def level Comp object counts
    MPI_Bcast(&count, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    
    // don't need sync
    tet->setCount(lsidx, count);
    _updateSpec(tet, sidx);
    _updateSum();
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTetAmount(uint tidx, uint sidx) const
{
    // following method does all necessary argument checking
    double count = _getTetCount(tidx, sidx);
    return count/steps::math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setTetAmount(uint tidx, uint sidx, double m)
{
    // convert amount in mols to number of molecules
    double m2 = m * steps::math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setTetCount(tidx, sidx, m2);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTetConc(uint tidx, uint sidx) const
{
    // following method does all necessary argument checking
    double count = _getTetCount(tidx, sidx);
    smtos::Tet * tet = pTets[tidx];
    double vol = tet->vol();
    return (count/(1.0e3 * vol * steps::math::AVOGADRO));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setTetConc(uint tidx, uint sidx, double c)
{
    assert (c >= 0.0);
    assert (tidx < pTets.size());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
        throw steps::ArgErr(os.str());
    }

    smtos::Tet * tet = pTets[tidx];
    double count = c * (1.0e3 * tet->vol() * steps::math::AVOGADRO);
    // the following method does all the necessary argument checking
    _setTetCount(tidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::TetOpSplitP::_getTetClamped(uint tidx, uint sidx) const
{
    assert (tidx < pTets.size());
    assert (sidx < statedef()->countSpecs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    smtos::Tet * tet = pTets[tidx];

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

void smtos::TetOpSplitP::_setTetClamped(uint tidx, uint sidx, bool buf)
{
    assert (tidx < pTets.size());
    assert (sidx < statedef()->countSpecs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    smtos::Tet * tet = pTets[tidx];

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

double smtos::TetOpSplitP::_getTetReacK(uint tidx, uint ridx) const
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());

    if (pTets[tidx] == 0 && tetHosts[tidx] == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }
    
    int host = tetHosts[tidx];
    
    smtos::Tet * tet = pTets[tidx];

    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }
    double kcst = 0;
    if (tet->getInHost()) {
        kcst = tet->reac(lridx)->kcst();
    }
    MPI_Bcast(&kcst, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return kcst;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setTetReacK(uint tidx, uint ridx, double kf)
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());
    assert (kf >= 0.0);

    if (pTets[tidx] == 0  && tetHosts[tidx] == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }
    
    smtos::Tet * tet = pTets[tidx];

    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "\nReaction undefined in tetrahedron.";
        throw steps::ArgErr(os.str());
    }

    if (!tet->getInHost()) return;
    tet->reac(lridx)->setKcst(kf);
    _updateElement(tet->reac(lridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::TetOpSplitP::_getTetReacActive(uint tidx, uint ridx) const
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());

    if (pTets[tidx] == 0 && tetHosts[tidx] == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }
    int host = tetHosts[tidx];
    smtos::Tet * tet = pTets[tidx];

    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }
    
    bool active = false;
    if (tet->getInHost()) {
        if (tet->reac(lridx)->inactive() == true) active = false;
        else active = true;
    }
    MPI_Bcast(&active, 1, MPI_UNSIGNED_SHORT, host, MPI_COMM_WORLD);
    return active;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setTetReacActive(uint tidx, uint ridx, bool act)
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());

    if (pTets[tidx] == 0 && tetHosts[tidx] == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    smtos::Tet * tet = pTets[tidx];

    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }
    if (!tet->getInHost()) return;
    tet->reac(lridx)->setActive(act);
    _updateElement(tet->reac(lridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTetDiffD(uint tidx, uint didx, uint direction_tet) const
{
    assert (tidx < pTets.size());
    assert (didx < statedef()->countDiffs());
    
    if (pTets[tidx] == 0 && tetHosts[tidx] == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }
    int host = tetHosts[tidx];
    smtos::Tet * tet = pTets[tidx];
    
    uint ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }
    double dcst = 0.0;
    if (tet->getInHost()) {
        if (direction_tet == std::numeric_limits<uint>::max()) {
            dcst = tet->diff(ldidx)->dcst();
        }
        else {
            int direction = tet->getTetDirection(direction_tet);
            if (direction == -1) {
                std::ostringstream os;
                os << "Tetrahedron " << direction_tet << " is not a neighbor of tetrahedron " << tidx << ".\n";
                throw steps::ArgErr(os.str());
            }

            dcst = tet->diff(ldidx)->dcst(direction);
        }
    }
    MPI_Bcast(&dcst, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return dcst;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setTetDiffD(uint tidx, uint didx, double dk, uint direction_tet)
{
    #ifdef DIRECTIONAL_DCST_DEBUG
    CLOG_IF(direction_tet !=  std::numeric_limits<uint>::max(), DEBUG, "steps_debug") << "tidx: " << tidx << " didx: " << didx << " dk: " << dk << " direction tet: " << direction_tet << "\n";
    
    CLOG_IF(direction_tet ==  std::numeric_limits<uint>::max(), DEBUG, "steps_debug") << "tidx: " << tidx << " didx: " << didx << " dk: " << dk << " for all directions.\n";
    #endif
    
    assert (tidx < pTets.size());
    assert (didx < statedef()->countDiffs());

    if (pTets[tidx] == 0 && tetHosts[tidx] == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }
    recomputeUpdPeriod = true;
    smtos::Tet * tet = pTets[tidx];

    uint ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }
    
    if (!tet->getInHost()) return;
    
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

bool smtos::TetOpSplitP::_getTetDiffActive(uint tidx, uint didx) const
{
    assert (tidx < pTets.size());
    assert (didx < statedef()->countDiffs());

    if (pTets[tidx] == 0 && tetHosts[tidx] == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }
    int host = tetHosts[tidx];
    smtos::Tet * tet = pTets[tidx];

    uint ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }
    bool active = false;
    if (tet->getInHost()) {
        if (tet->diff(ldidx)->inactive() == true) active = false;
        else active = true;
    }
    MPI_Bcast(&active, 1, MPI_UNSIGNED_SHORT, host, MPI_COMM_WORLD);
    return active;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setTetDiffActive(uint tidx, uint didx, bool act)
{
    assert (tidx < pTets.size());
    assert (didx < statedef()->countDiffs());

    if (pTets[tidx] == 0 && tetHosts[tidx] == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    smtos::Tet * tet = pTets[tidx];

    uint ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }
    if (!tet->getInHost()) return;
    tet->diff(ldidx)->setActive(act);

    recomputeUpdPeriod = true;

	_updateElement(tet->diff(ldidx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTetReacH(uint tidx, uint ridx) const
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());

    if (pTets[tidx] == 0 && tetHosts[tidx] == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }
    int host = tetHosts[tidx];
    smtos::Tet * tet = pTets[tidx];

    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }
    double h = 0;
    if (tet->getInHost()) {
        h = tet->reac(lridx)->h();
    }
    MPI_Bcast(&h, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return h;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTetReacC(uint tidx, uint ridx) const
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());

    if (pTets[tidx] == 0 && tetHosts[tidx] == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }
    int host = tetHosts[tidx];
    smtos::Tet * tet = pTets[tidx];

    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    double c = 0;
    if (tet->getInHost()) {
        c = tet->reac(lridx)->c();
    }
    MPI_Bcast(&c, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return c;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTetReacA(uint tidx, uint ridx) const
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());

    if (pTets[tidx] == 0 && tetHosts[tidx] == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }
    int host = tetHosts[tidx];
    smtos::Tet * tet = pTets[tidx];

    uint lridx = tet->compdef()->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }
    double a = 0;
    if (tet->getInHost()) {
        a = tet->reac(lridx)->rate();
    }
    MPI_Bcast(&a, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return a;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTetDiffA(uint tidx, uint didx) const
{
    assert (tidx < pTets.size());
    assert (didx < statedef()->countDiffs());

    if (pTets[tidx] == 0 && tetHosts[tidx] == -1)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }
    int host = tetHosts[tidx];
    smtos::Tet * tet = pTets[tidx];

    uint ldidx = tet->compdef()->diffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }
    double a = 0;
    if (tet->getInHost()) {
        a = tet->diff(ldidx)->rate();
    }
    MPI_Bcast(&a, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return a;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTriArea(uint tidx) const
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

void smtos::TetOpSplitP::_setTriArea(uint tidx, double area)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::TetOpSplitP::_getTriSpecDefined(uint tidx, uint sidx) const
{
    assert (tidx < pTris.size());
    assert (sidx < statedef()->countSpecs());

    if (pTris[tidx] == 0) return false;

    smtos::Tri * tri = pTris[tidx];
    uint lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return false;
    else return true;
}

////////////////////////////////////////////////////////////////////////////////


double smtos::TetOpSplitP::_getTriCount(uint tidx, uint sidx) const
{
    MPI_Barrier(MPI_COMM_WORLD);
    assert (tidx < pTris.size());
    assert (sidx < statedef()->countSpecs());

    if (pTris[tidx] <= 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    smtos::Tri * tri = pTris[tidx];
    uint lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    uint count = tri->pools()[lsidx];
    std::map<uint,uint>::const_iterator it = triHosts.find(tidx);
    if (it == triHosts.end()) {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a host.\n";
        throw steps::ArgErr(os.str());
    }
    
    MPI_Bcast(&count, 1, MPI_UNSIGNED, it->second, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    return count;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setTriCount(uint tidx, uint sidx, double n)
{
    #ifdef MPI_DEBUG
    //CLOG(DEBUG, "mpi_debug") << "Set Tri " << tidx << " count of spec " << sidx << ": " << n << "\n";
    #endif
    MPI_Barrier(MPI_COMM_WORLD);
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

    smtos::Tri * tri = pTris[tidx];
    uint lsidx = tri->patchdef()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    uint count = 0;
    
    if (myRank == 0) {
        double n_int = std::floor(n);
        double n_frc = n - n_int;
        count = static_cast<uint>(n_int);
        if (n_frc > 0.0)
        {
            double rand01 = rng()->getUnfIE();
            if (rand01 < n_frc) count++;
        }
    }

    MPI_Bcast(&count, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    
    tri->setCount(lsidx, count);
    _updateSpec(tri, sidx);
    _updateSum();
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTriAmount(uint tidx, uint sidx) const
{

    // following method does all necessary argument checking
    double count = _getTriCount(tidx, sidx);
    return count/steps::math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setTriAmount(uint tidx, uint sidx, double m)
{
    // convert amount in mols to number of molecules
    double m2 = m * steps::math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setTriCount(tidx, sidx, m2);
}

////////////////////////////////////////////////////////////////////////////////


bool smtos::TetOpSplitP::_getTriClamped(uint tidx, uint sidx) const
{
    assert (tidx < pTris.size());
    assert (sidx < statedef()->countSpecs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    smtos::Tri * tri = pTris[tidx];

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

void smtos::TetOpSplitP::_setTriClamped(uint tidx, uint sidx, bool buf)
{
    assert (tidx < pTris.size());
    assert (sidx < statedef()->countSpecs());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    smtos::Tri * tri = pTris[tidx];

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

double smtos::TetOpSplitP::_getTriSReacK(uint tidx, uint ridx)
{
    assert (tidx < pTris.size());
    assert (ridx < statedef()->countSReacs());

    if (pTris[tidx] == 0 && triHosts.find(tidx) == triHosts.end())
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }
    uint host = triHosts[tidx];
    smtos::Tri * tri = pTris[tidx];

    uint lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }
    double kcst = 0;
    if (tri->getInHost()) {
        kcst = tri->sreac(lsridx)->kcst();
    }
    MPI_Bcast(&kcst, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return kcst;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setTriSReacK(uint tidx, uint ridx, double kf)
{
    assert (tidx < pTris.size());
    assert (ridx < statedef()->countSReacs());

    if (pTris[tidx] == 0 && triHosts.find(tidx) == triHosts.end())
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }
    smtos::Tri * tri = pTris[tidx];

    uint lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }
    if (!tri->getInHost()) return;
    
    tri->sreac(lsridx)->setKcst(kf);
    _updateElement(tri->sreac(lsridx));
    _updateSum();

}

////////////////////////////////////////////////////////////////////////////////

bool smtos::TetOpSplitP::_getTriSReacActive(uint tidx, uint ridx)
{
    assert (tidx < pTris.size());
    assert (ridx < statedef()->countSReacs());

    if (pTris[tidx] == 0 && triHosts.find(tidx) == triHosts.end())
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }
    uint host = triHosts[tidx];
    smtos::Tri * tri = pTris[tidx];

    uint lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }
    bool active = false;
    if (tri->getInHost()) {
        if (tri->sreac(lsridx)->inactive() == true)  active = false;
        else  active = true;
    }
    MPI_Bcast(&active, 1, MPI_UNSIGNED_SHORT, host, MPI_COMM_WORLD);
    return active;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setTriSReacActive(uint tidx, uint ridx, bool act)
{
    assert (tidx < pTris.size());
    assert (ridx < statedef()->countSReacs());

    if (pTris[tidx] == 0 && tetHosts[tidx] == -1)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    smtos::Tri * tri = pTris[tidx];

    uint lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }
    if (!tri->getInHost()) return;
    tri->sreac(lsridx)->setActive(act);
    _updateElement(tri->sreac(lsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTriSDiffD(uint tidx, uint didx, uint direction_tri)
{
    assert (tidx < pTris.size());
    assert (didx < statedef()->countSurfDiffs());
    
    if (pTris[tidx] == 0 && triHosts.find(tidx) == triHosts.end())
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }
    uint host = triHosts[tidx];
    smtos::Tri * tri = pTris[tidx];
    
    uint ldidx = tri->patchdef()->surfdiffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }
    double dcst = 0.0;
    if (tri->getInHost()) {
        if (direction_tri == std::numeric_limits<uint>::max()) {
            dcst = tri->sdiff(ldidx)->dcst();
            
        }
        else {
            int direction = tri->getTriDirection(direction_tri);
            if (direction == -1) {
                std::ostringstream os;
                os << "Triangle " << direction_tri << " is not a neighbor of triangle " << tidx << ".\n";
                throw steps::ArgErr(os.str());
            }
            
            dcst = tri->sdiff(ldidx)->dcst(direction);
        }
    }
    MPI_Bcast(&dcst, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return dcst;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setTriSDiffD(uint tidx, uint didx, double dk, uint direction_tri)
{
    #ifdef DIRECTIONAL_DCST_DEBUG
    CLOG_IF(direction_tri !=  std::numeric_limits<uint>::max(), DEBUG, "steps_debug") << "tidx: " << tidx << " didx: " << didx << " dk: " << dk << " direction tri: " << direction_tri << "\n";
    
    CLOG_IF(direction_tri ==  std::numeric_limits<uint>::max(), DEBUG, "steps_debug") << "tidx: " << tidx << " didx: " << didx << " dk: " << dk << " for all directions.\n";
    #endif
    
    assert (tidx < pTris.size());
    assert (didx < statedef()->countSurfDiffs());

    if (pTris[tidx] == 0 && triHosts.find(tidx) == triHosts.end())
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    smtos::Tri * tri = pTris[tidx];

    uint ldidx = tri->patchdef()->surfdiffG2L(didx);
    if (ldidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Diffusion rule undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }
    recomputeUpdPeriod = true;
    if (!tri->getInHost()) return;
    
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

bool smtos::TetOpSplitP::_getTriVDepSReacActive(uint tidx, uint vsridx)
{
    assert (tidx < pTris.size());
    assert (vsridx < statedef()->countVDepSReacs());

    if (pTris[tidx] == 0 && triHosts.find(tidx) == triHosts.end())
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }
    uint host = triHosts[tidx];
    smtos::Tri * tri = pTris[tidx];

    uint lvsridx = tri->patchdef()->vdepsreacG2L(vsridx);
    if (lvsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Voltage-dependent surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }
    bool active = false;
    if (tri->getInHost()) {
        if (tri->vdepsreac(lvsridx)->inactive() == true)  active = false;
        else  active = true;
    }
    MPI_Bcast(&active, 1, MPI_UNSIGNED_SHORT, host, MPI_COMM_WORLD);
    return active;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setTriVDepSReacActive(uint tidx, uint vsridx, bool act)
{
    assert (tidx < pTris.size());
    assert (vsridx < statedef()->countVDepSReacs());

    if (pTris[tidx] == 0 && triHosts.find(tidx) == triHosts.end())
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }
    uint host = triHosts[tidx];
    smtos::Tri * tri = pTris[tidx];

    uint lvsridx = tri->patchdef()->vdepsreacG2L(vsridx);
    if (lvsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Voltage-dependent surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }
    if (!tri->getInHost()) return;
    tri->vdepsreac(lvsridx)->setActive(act);
    _updateElement(tri->vdepsreac(lvsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTriSReacH(uint tidx, uint ridx)
{
    assert (tidx < pTris.size());
    assert (ridx < statedef()->countSReacs());

    if (pTris[tidx] == 0 && triHosts.find(tidx) == triHosts.end())
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }
    uint host = triHosts[tidx];
    smtos::Tri * tri = pTris[tidx];

    uint lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }
    double h = 0;
    if (tri->getInHost()) h = tri->sreac(lsridx)->h();
    MPI_Bcast(&h, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return h;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTriSReacC(uint tidx, uint ridx)
{
    assert (tidx < pTris.size());
    assert (ridx < statedef()->countSReacs());

    if (pTris[tidx] == 0 && triHosts.find(tidx) == triHosts.end())
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }
    uint host = triHosts[tidx];
    smtos::Tri * tri = pTris[tidx];

    uint lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    double c = 0;
    if (tri->getInHost()) c = tri->sreac(lsridx)->c();
    MPI_Bcast(&c, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return c;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTriSReacA(uint tidx, uint ridx)
{
    assert (tidx < pTris.size());
    assert (ridx < statedef()->countSReacs());

    if (pTris[tidx] == 0 && triHosts.find(tidx) == triHosts.end())
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }
    uint host = triHosts[tidx];
    smtos::Tri * tri = pTris[tidx];

    uint lsridx = tri->patchdef()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    double a = 0;
    if (tri->getInHost()) a =  tri->sreac(lsridx)->rate();
    MPI_Bcast(&a, 1, MPI_DOUBLE, host, MPI_COMM_WORLD);
    return a;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::setEfieldDT(double efdt)
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

double smtos::TetOpSplitP::_getTetV(uint tidx) const
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

void smtos::TetOpSplitP::_setTetV(uint tidx, double v)
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

bool smtos::TetOpSplitP::_getTetVClamped(uint tidx) const
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

void smtos::TetOpSplitP::_setTetVClamped(uint tidx, bool cl)
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

double smtos::TetOpSplitP::_getTriV(uint tidx) const
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
    return EFTrisV[loctidx];
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setTriV(uint tidx, double v)
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
    EFTrisV[loctidx] = v;
    pEField->setTriV(loctidx, v);
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::TetOpSplitP::_getTriVClamped(uint tidx) const
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

void smtos::TetOpSplitP::_setTriVClamped(uint tidx, bool cl)
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

double smtos::TetOpSplitP::_getTriOhmicI(uint tidx)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }

    smtos::Tri * tri = pTris[tidx];

    int loctidx = pEFTri_GtoL[tidx];
    if (loctidx == -1)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        throw steps::ArgErr(os.str());
    }
    int tri_host = triHosts[tidx];
    double cur = 0.0;
    if (tri->getInHost()) {
        cur = tri->getOhmicI(EFTrisV[loctidx], efdt());
    }
    MPI_Bcast(&cur, 1, MPI_DOUBLE, tri_host, MPI_COMM_WORLD);
    return cur;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTriOhmicI(uint tidx, uint ocidx)
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

    smtos::Tri * tri = pTris[tidx];

    uint locidx = tri->patchdef()->ohmiccurrG2L(ocidx);
    if (locidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Ohmic current undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }
    int tri_host = triHosts[tidx];
    double cur = 0.0;
    if (tri->getInHost()) {
        cur = tri->getOhmicI(locidx, EFTrisV[loctidx], efdt());
    }
    MPI_Bcast(&cur, 1, MPI_DOUBLE, tri_host, MPI_COMM_WORLD);
    return cur;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTriGHKI(uint tidx)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }

    smtos::Tri * tri = pTris[tidx];
    
    int tri_host = triHosts[tidx];
    double cur = 0.0;
    if (tri->getInHost()) {
        cur = tri->getGHKI(efdt());
    }
    MPI_Bcast(&cur, 1, MPI_DOUBLE, tri_host, MPI_COMM_WORLD);
    return cur;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTriGHKI(uint tidx, uint ghkidx)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }

    smtos::Tri * tri = pTris[tidx];

    uint locidx = tri->patchdef()->ghkcurrG2L(ghkidx);
    if (locidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "GHK current undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    int tri_host = triHosts[tidx];
    double cur = 0.0;
    if (tri->getInHost()) {
        cur = tri->getGHKI(locidx, efdt());
    }
    MPI_Bcast(&cur, 1, MPI_DOUBLE, tri_host, MPI_COMM_WORLD);
    return cur;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getTriI(uint tidx) const
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

void smtos::TetOpSplitP::_setVertIClamp(uint vidx, double cur)
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

void smtos::TetOpSplitP::_setTriIClamp(uint tidx, double cur)
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

void smtos::TetOpSplitP::_setTriCapac(uint tidx, double cap)
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
    pEField->setTriCapac(loctidx, cap);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::_getVertV(uint vidx) const
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
    return pEField->getVertV(locvidx);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_setVertV(uint vidx, double v)
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

bool smtos::TetOpSplitP::_getVertVClamped(uint vidx) const
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

void smtos::TetOpSplitP::_setVertVClamped(uint vidx, bool cl)
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

void smtos::TetOpSplitP::_setMembRes(uint midx, double ro, double vrev)
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

void smtos::TetOpSplitP::_setMembPotential(uint midx, double v)
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

void smtos::TetOpSplitP::_setMembCapac(uint midx, double cm)
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

void smtos::TetOpSplitP::_setMembVolRes(uint midx, double ro)
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

void smtos::TetOpSplitP::_computeUpdPeriod(void)
{
    double local_max_rate = 0.0;
    
    for (uint pos = 0; pos < diffSep; pos++) {
        Diff* d = pDiffs[pos];
        // Now ignoring inactive diffusion
        double scaleddcst = 0.0;
        if(d->active()) scaleddcst = d->getScaledDcst();
        if (scaleddcst > local_max_rate) local_max_rate = scaleddcst;
    }
    
    for (uint pos = 0; pos < sdiffSep; pos++) {
        SDiff* d = pSDiffs[pos];
        // Now ignoring inactive diffusion
        double scaleddcst = 0.0;
        if(d->active()) scaleddcst = d->getScaledDcst();
        if (scaleddcst > local_max_rate) local_max_rate = scaleddcst;
    }
    // get global max rate
    double global_max_rate = 0;
    MPI_Allreduce(&local_max_rate, &global_max_rate, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    if (global_max_rate < 0.0)
    {
        std::ostringstream os;
        os << "Maximum scaled diffusion constant is " << global_max_rate << ". This should not happen in this solver.\n";
        throw steps::ArgErr(os.str());
    }
    updPeriod = 1.0 / global_max_rate;
    recomputeUpdPeriod = false;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_updateLocal(std::set<KProc*> const & upd_entries) {
    for (auto kp : upd_entries) {
        assert(kp != NULL);
        _updateElement(kp);
    }
    
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_updateLocal(std::vector<KProc*> const & upd_entries) {
    for (auto kp : upd_entries) {
        assert(kp != NULL);
        _updateElement(kp);
    }
    
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_updateLocal(std::vector<uint> const & upd_entries) {
    for (auto upd_idx : upd_entries) {
        KProc* kp = pKProcs[upd_idx];
        if (kp != NULL)
        {
            _updateElement(kp);
        }
    }
    
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_updateLocal(void) {
    for (uint i = 0; i < nEntries; i++) {
        KProc* kp = pKProcs[i];
        if (kp != NULL)
        {
            _updateElement(kp);
        }
    }
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

smtos::CRGroup* smtos::TetOpSplitP::_getGroup(int pow) {
    if (pow >= 0) {
        return pGroups[pow];
    }
    else {
        return nGroups[-pow];
    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_extendPGroups(uint new_size) {
    uint curr_size = pGroups.size();

    while (curr_size < new_size) {
        pGroups.push_back(new CRGroup(curr_size));
        curr_size ++;
    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_extendNGroups(uint new_size) {

    uint curr_size = nGroups.size();

    while (curr_size < new_size) {
        nGroups.push_back(new CRGroup(-curr_size));
        curr_size ++;
    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_extendGroup(CRGroup* group, uint size) {
    group->capacity += size;
    KProc** old_ptr = group->indices;
    group->indices = (KProc**)realloc(group->indices,
                                      sizeof(KProc*) * group->capacity);
    if (group->indices == NULL) {
        std::cerr << "DirectCR: unable to allocate memory for SSA group.\n";
        throw;
    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_updateSum(void) {
    pA0 = 0.0;

    uint n_neg_groups = nGroups.size();
    uint n_pos_groups = pGroups.size();

    for (uint i = 0; i < n_neg_groups; i++) {
        pA0 += nGroups[i]->sum;
    }

    for (uint i = 0; i < n_pos_groups; i++) {
        pA0 += pGroups[i]->sum;
    }
    #ifdef MPI_DEBUG
    //CLOG(DEBUG, "mpi_debug") << "New A0 " << pA0 << ".\n";
    #endif
}


////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_updateElement(KProc* kp)
{
    
    if (kp->getType() == KP_DIFF || kp->getType() == KP_SDIFF) {
        kp->crData.rate = kp->rate(this);;
        return;
    }
  
    #ifdef MPI_DEBUG
    //Reac* reac_upcast = dynamic_cast<Reac*>(kp);
    //if (reac_upcast != NULL) {
    //    CLOG(DEBUG, "mpi_debug") << "Update Reac " << kp->schedIDX() << ".\n";
    //}
    //SReac* sreac_upcast = dynamic_cast<SReac*>(kp);
    //if (sreac_upcast != NULL) {
    //    CLOG(DEBUG, "mpi_debug") << "Update SReac " << kp->schedIDX() << ".\n";
    //}
    #endif

    double new_rate = kp->rate(this);

    CRKProcData & data = kp->crData;
    double old_rate = data.rate;

    data.rate = new_rate;

    #ifdef MPI_DEBUG
    //CLOG(DEBUG, "mpi_debug") << "Old rate: " << old_rate << "new rate: " << new_rate <<".\n";
    #endif

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
// Batch data recording
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

std::vector<double> smtos::TetOpSplitP::getBatchTetCounts(std::vector<uint> const & tets, std::string const & s) const
{
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;

    uint ntets = tets.size();
    uint sgidx = statedef()->getSpecIdx(s);
    std::vector<double> local_counts(ntets, 0.0);
    std::vector<double> global_counts(ntets, 0.0);

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

        smtos::Tet * tet = pTets[tidx];
        uint slidx = tet->compdef()->specG2L(sgidx);
        if (slidx == ssolver::LIDX_UNDEFINED)
        {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }
        if (tet->getInHost()) {
            local_counts[t] = tet->pools()[slidx];
        }
        
    }

    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, fill in zeros at target positions:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        std::cerr << "Warning: Species " << s << " has not been defined in the following tetrahedrons, fill in zeros at target positions:\n";
        std::cerr << spec_undefined.str() << "\n";
    }
    MPI_Allreduce(&local_counts.front(), &global_counts.front(), ntets, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global_counts;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> smtos::TetOpSplitP::getBatchTriCounts(std::vector<uint> const & tris, std::string const & s) const
{
    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;


    uint ntris = tris.size();
    uint sgidx = statedef()->getSpecIdx(s);
    std::vector<double> local_counts(ntris, 0.0);
    std::vector<double> global_counts(ntris, 0.0);
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

        smtos::Tri * tri = pTris[tidx];
        uint slidx = tri->patchdef()->specG2L(sgidx);
        if (slidx == ssolver::LIDX_UNDEFINED)
        {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }
        if (tri->getInHost()) {
            local_counts[t] = tri->pools()[slidx];
        }
        
    }

    if (has_tri_warning) {
        std::cerr << "Warning: The following triangles have not been assigned to a patch, fill in zeros at target positions:\n";
        std::cerr << tri_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        std::cerr << "Warning: Species " << s << " has not been defined in the following triangles, fill in zeros at target positions:\n";
        std::cerr << spec_undefined.str() << "\n";
    }
    
    MPI_Allreduce(&local_counts.front(), &global_counts.front(), ntris, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global_counts;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::getBatchTetCountsNP(unsigned int* indices, int input_size, std::string const & s, double* counts, int output_size) const
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
    std::vector<double> local_counts(input_size, 0.0);
    
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

        smtos::Tet * tet = pTets[tidx];
        uint slidx = tet->compdef()->specG2L(sgidx);
        if (slidx == ssolver::LIDX_UNDEFINED)
        {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }
        
        if (tet->getInHost()) {
            local_counts[t] = tet->pools()[slidx];
        }
        
    }

    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, fill in zeros at target positions:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }

    if (has_spec_warning) {
        std::cerr << "Warning: Species " << s << " has not been defined in the following tetrahedrons, fill in zeros at target positions:\n";
        std::cerr << spec_undefined.str() << "\n";
    }
    MPI_Allreduce(&local_counts.front(), counts, input_size, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::getBatchTriCountsNP(unsigned int* indices, int input_size, std::string const & s, double* counts, int output_size) const
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
    std::vector<double> local_counts(input_size, 0.0);
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

        smtos::Tri * tri = pTris[tidx];
        uint slidx = tri->patchdef()->specG2L(sgidx);
        if (slidx == ssolver::LIDX_UNDEFINED)
        {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }
        if (tri->getInHost()) {
            local_counts[t] = tri->pools()[slidx];
        }
        
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

double smtos::TetOpSplitP::sumBatchTetCountsNP(unsigned int* indices, int input_size, std::string const & s)
{
    bool has_tet_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tet_not_assign;
    std::ostringstream spec_undefined;
    
    uint sgidx = statedef()->getSpecIdx(s);
    double partial_sum = 0.0;
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
        
        smtos::Tet * tet = pTets[tidx];
        uint slidx = tet->compdef()->specG2L(sgidx);
        if (slidx == ssolver::LIDX_UNDEFINED)
        {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }
        
        if (tet->getInHost()) {
            partial_sum += tet->pools()[slidx];
        }
    }
    
    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, fill in zeros at target positions:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }
    
    if (has_spec_warning) {
        std::cerr << "Warning: Species " << s << " has not been defined in the following tetrahedrons, fill in zeros at target positions:\n";
        std::cerr << spec_undefined.str() << "\n";
    }
    double global_sum = 0.0;
    MPI_Allreduce(&partial_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_sum;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::sumBatchTriCountsNP(unsigned int* indices, int input_size, std::string const & s)
{
    bool has_tri_warning = false;
    bool has_spec_warning = false;
    std::ostringstream tri_not_assign;
    std::ostringstream spec_undefined;
    
    double partial_sum = 0.0;
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
        
        smtos::Tri * tri = pTris[tidx];
        uint slidx = tri->patchdef()->specG2L(sgidx);
        if (slidx == ssolver::LIDX_UNDEFINED)
        {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }
        if (tri->getInHost()) {
            partial_sum += tri->pools()[slidx];
        }
    }
    
    if (has_tri_warning) {
        std::cerr << "Warning: The following triangles have not been assigned to a patch, fill in zeros at target positions:\n";
        std::cerr << tri_not_assign.str() << "\n";
    }
    
    if (has_spec_warning) {
        std::cerr << "Warning: Species " << s << " has not been defined in the following triangles, fill in zeros at target positions:\n";
        std::cerr << spec_undefined.str() << "\n";
    }
    
    double global_sum = 0.0;
    MPI_Allreduce(&partial_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_sum;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::sumBatchTriGHKIsNP(unsigned int* indices, int input_size, std::string const & ghk)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    // the following may raise exception if string is unknown
    uint ghkidx = statedef()->getGHKcurrIdx(ghk);
    
    
    double partial_sum = 0.0;
    
    for (uint t = 0; t < input_size; t++) {
        uint tidx = indices[t];
        
        if (tidx >= mesh()->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            throw steps::ArgErr(os.str());
        }
        
        smtos::Tri * tri = pTris[tidx];
        
        uint locidx = tri->patchdef()->ghkcurrG2L(ghkidx);
        if (locidx == ssolver::LIDX_UNDEFINED)
        {
            std::ostringstream os;
            os << "GHK current undefined in triangle.\n";
            throw steps::ArgErr(os.str());
        }
        
        if (tri->getInHost()) {
            partial_sum += tri->getGHKI(locidx, efdt());
        }
    }
    double global_sum = 0.0;
    MPI_Allreduce(&partial_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_sum;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::sumBatchTriOhmicIsNP(unsigned int* indices, int input_size, std::string const & oc)
{
    if (efflag() != true)
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        throw steps::ArgErr(os.str());
    }
    
    // the following may raise exception if string is unknown
    uint ocidx = statedef()->getOhmicCurrIdx(oc);
    
    
    double partial_sum = 0.0;
    
    for (uint t = 0; t < input_size; t++) {
        uint tidx = indices[t];
        int loctidx = pEFTri_GtoL[tidx];
        if (tidx >= mesh()->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            throw steps::ArgErr(os.str());
        }
        
        smtos::Tri * tri = pTris[tidx];
        
        uint locidx = tri->patchdef()->ohmiccurrG2L(ocidx);
        if (locidx == ssolver::LIDX_UNDEFINED)
        {
            std::ostringstream os;
            os << "Ohmic current undefined in triangle.\n";
            throw steps::ArgErr(os.str());
        }
        
        if (tri->getInHost()) {
            partial_sum += tri->getOhmicI(locidx, EFTrisV[loctidx], efdt());
        }
    }
    double global_sum = 0.0;
    MPI_Allreduce(&partial_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_sum;
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// ROI Data Access
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

std::vector<double> smtos::TetOpSplitP::getROITetCounts(std::string ROI_id, std::string const & s) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int inputsize = mesh()->getROIDataSize(ROI_id);
    std::vector<double> data(inputsize);
    
    getBatchTetCountsNP(indices, inputsize, s, const_cast<double*>(&data.front()), data.size());
    return data;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> smtos::TetOpSplitP::getROITriCounts(std::string ROI_id, std::string const & s) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TRI)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int inputsize = mesh()->getROIDataSize(ROI_id);
    std::vector<double> data(inputsize);
    
    getBatchTriCountsNP(indices, inputsize, s, const_cast<double*>(&data.front()), data.size());
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::getROITetCountsNP(std::string ROI_id, std::string const & s, double* counts, int output_size) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int inputsize = mesh()->getROIDataSize(ROI_id);
    
    getBatchTetCountsNP(indices, inputsize, s, counts, output_size);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::getROITriCountsNP(std::string ROI_id, std::string const & s, double* counts, int output_size) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TRI)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int inputsize = mesh()->getROIDataSize(ROI_id);
    
    getBatchTriCountsNP(indices, inputsize, s, counts, output_size);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::getROIVol(std::string ROI_id) const
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

double smtos::TetOpSplitP::getROIArea(std::string ROI_id) const
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

double smtos::TetOpSplitP::getROICount(std::string ROI_id, std::string const & s) const
{
    steps::tetmesh::ElementType type = mesh()->getROIType(ROI_id);
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    double local_sum = 0.0;
    double global_sum = 0.0;
    
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
            
            smtos::Tri * tri = pTris[tidx];
            uint slidx = tri->patchdef()->specG2L(sgidx);
            if (slidx == ssolver::LIDX_UNDEFINED)
            {
                spec_undefined << tidx << " ";
                has_spec_warning = true;
                continue;
            }
            
            // compute local sum for each process
            if (tri->getInHost()) local_sum += tri->pools()[slidx];
        }
        
        // gather global sum
        MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        
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
            
            smtos::Tet * tet = pTets[tidx];
            uint slidx = tet->compdef()->specG2L(sgidx);
            if (slidx == ssolver::LIDX_UNDEFINED)
            {
                spec_undefined << tidx << " ";
                has_spec_warning = true;
                continue;
            }
            
            // compute local sum for each process
            if (tet->getInHost()) local_sum += tet->pools()[slidx];
        }
        
        // gather global sum
        MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        
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
    
    return global_sum;
}

////////////////////////////////////////////////////////////////////////////////

// WEILIANG: Can we apply Sam's set count method (as in e.g. setCompCount) here?
void smtos::TetOpSplitP::setROICount(std::string ROI_id, std::string const & s, double count)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (count > std::numeric_limits<unsigned int>::max( ))
    {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << std::numeric_limits<unsigned int>::max( ) << ").\n";
        throw steps::ArgErr(os.str());
    }

    steps::tetmesh::ElementType type = mesh()->getROIType(ROI_id);
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    std::vector<uint> apply_indices;
    
    if (type == steps::tetmesh::ELEM_TRI) {
        double totalarea = 0.0;
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
            
            smtos::Tri * tri = pTris[tidx];
            uint slidx = tri->patchdef()->specG2L(sgidx);
            if (slidx == ssolver::LIDX_UNDEFINED)
            {
                spec_undefined << tidx << " ";
                has_spec_warning = true;
                continue;
            }
            
            apply_indices.push_back(tidx);
            totalarea += tri->area();
        }
        
        if (has_tri_warning) {
            std::cerr << "Warning: The following triangles have not been assigned to a patch, fill in zeros at target positions:\n";
            std::cerr << tri_not_assign.str() << "\n";
        }
        
        if (has_spec_warning) {
            std::cerr << "Warning: Species " << s << " has not been defined in the following triangles, fill in zeros at target positions:\n";
            std::cerr << spec_undefined.str() << "\n";
        }

        uint ind_size = apply_indices.size();
        
        std::vector<double> apply_count(ind_size, 0.0);
        
        // Create distribution at node 0
        if (myRank == 0) {
            uint c = static_cast<uint>(count);
            uint nremoved = 0;
            
            for (uint t = 0; t < ind_size; t++)
            {
                uint tidx = apply_indices[t];
                smtos::Tri * tri = pTris[tidx];
                
                if ((count == 0.0) || (nremoved == c)) break;
                
                double fract = static_cast<double>(c) * (tri->area() / totalarea);
                uint n3 = static_cast<uint>(std::floor(fract));
                
                double n3_frac = fract - static_cast<double>(n3);
                if (n3_frac > 0.0)
                {
                    double rand01 = rng()->getUnfIE();
                    if (rand01 < n3_frac) n3++;
                }
                
                nremoved += n3;
                
                if (nremoved >= c)
                {
                    n3 -= (nremoved-c);
                    nremoved = c;
                }
                
                apply_count[t] = n3;
            }
            assert(nremoved <= c);
            c -= nremoved;
            while (c != 0)
            {
                double accum = 0.0;
                double selector = rng()->getUnfIE() * totalarea;
                for (uint t = 0; t < ind_size; t++)
                {
                    uint tidx = apply_indices[t];
                    smtos::Tri * tri = pTris[tidx];
                    accum += tri->area();
                    if (selector < accum) {
                        apply_count[t] += 1.0;
                        break;
                    }
                }
                c--;
            }
        }
        
        MPI_Bcast(&(apply_count.front()), ind_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        // counts need to be set globally so that sync can be avoided
        for (uint t = 0; t < ind_size; t++)
        {
            uint tidx = apply_indices[t];
            smtos::Tri * tri = pTris[tidx];
            
            uint slidx = tri->patchdef()->specG2L(sgidx);
            tri->setCount(slidx, apply_count[t]);
            _updateSpec(tri, slidx);
        }
        _updateSum();
    }
    else if (type == steps::tetmesh::ELEM_TET) {
        double totalvol = 0.0;
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
            
            smtos::Tet * tet = pTets[tidx];
            uint slidx = tet->compdef()->specG2L(sgidx);
            if (slidx == ssolver::LIDX_UNDEFINED)
            {
                spec_undefined << tidx << " ";
                has_spec_warning = true;
                continue;
            }
            
            apply_indices.push_back(tidx);
            totalvol += tet->vol();
        }
        
        if (has_tet_warning) {
            std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, fill in zeros at target positions:\n";
            std::cerr << tet_not_assign.str() << "\n";
        }
        
        if (has_spec_warning) {
            std::cerr << "Warning: Species " << s << " has not been defined in the following tetrahedrons, fill in zeros at target positions:\n";
            std::cerr << spec_undefined.str() << "\n";
        }

        uint ind_size = apply_indices.size();
        
        std::vector<double> apply_count(ind_size, 0.0);
        
        // Create distribution at node 0
        if (myRank == 0) {
            uint c = static_cast<uint>(count);
            uint nremoved = 0;
            
            for (uint t = 0; t < ind_size; t++)
            {
                uint tidx = apply_indices[t];
                smtos::Tet * tet = pTets[tidx];
                
                if ((count == 0.0) || (nremoved == c)) break;
                
                double fract = static_cast<double>(c) * (tet->vol() / totalvol);
                uint n3 = static_cast<uint>(std::floor(fract));
                
                double n3_frac = fract - static_cast<double>(n3);
                if (n3_frac > 0.0)
                {
                    double rand01 = rng()->getUnfIE();
                    if (rand01 < n3_frac) n3++;
                }
                
                nremoved += n3;
                
                if (nremoved >= c)
                {
                    n3 -= (nremoved-c);
                    nremoved = c;
                }
                
                apply_count[t] = n3;
            }
            assert(nremoved <= c);
            c -= nremoved;
            while (c != 0)
            {
                double accum = 0.0;
                double selector = rng()->getUnfIE() * totalvol;
                for (uint t = 0; t < ind_size; t++)
                {
                    uint tidx = apply_indices[t];
                    smtos::Tet * tet = pTets[tidx];
                    accum += tet->vol();
                    if (selector < accum) {
                        apply_count[t] += 1.0;
                        break;
                    }
                }
                c--;
            }
        }
        
        MPI_Bcast(&(apply_count.front()), ind_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        // set the counts golbally and update local KProcs
        for (uint t = 0; t < ind_size; t++)
        {
            uint tidx = apply_indices[t];
            smtos::Tet * tet = pTets[tidx];
            uint slidx = tet->compdef()->specG2L(sgidx);
            tet->setCount(slidx, apply_count[t]);
            _updateSpec(tet, slidx);
        }
        
        _updateSum();
    }
    else {
        std::ostringstream os;
        os << "Error: Cannot find suitable ROI for the function call getROICount.\n";
        throw steps::ArgErr(os.str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::getROIAmount(std::string ROI_id, std::string const & s) const
{
    double count = getROICount(ROI_id, s);
    return (count / smath::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////
// WEILIANG: Can we apply Sam's set count method (as in e.g. setCompCount) here?
void smtos::TetOpSplitP::setROIConc(std::string ROI_id, std::string const & s, double conc)
{
    if (conc > std::numeric_limits<unsigned int>::max( ))
    {
        std::ostringstream os;
        os << "Can't set count greater than maximum unsigned integer (";
        os << std::numeric_limits<unsigned int>::max( ) << ").\n";
        throw steps::ArgErr(os.str());
    }
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = mesh()->_getROIData(ROI_id);
    int datasize = mesh()->getROIDataSize(ROI_id);
    
    std::vector<uint> apply_indices;
    
    double totalvol = 0.0;
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
        
        smtos::Tet * tet = pTets[tidx];
        uint slidx = tet->compdef()->specG2L(sgidx);
        if (slidx == ssolver::LIDX_UNDEFINED)
        {
            spec_undefined << tidx << " ";
            has_spec_warning = true;
            continue;
        }
        
        apply_indices.push_back(tidx);
        totalvol += tet->vol();
    }
    
    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, fill in zeros at target positions:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }
    
    if (has_spec_warning) {
        std::cerr << "Warning: Species " << s << " has not been defined in the following tetrahedrons, fill in zeros at target positions:\n";
        std::cerr << spec_undefined.str() << "\n";
    }
    
    int count = conc * (1.0e3 * totalvol * steps::math::AVOGADRO);
    
    setROICount(ROI_id, s, count);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::getROIConc(std::string ROI_id, std::string const & s) const
{
    if (!mesh()->checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    double count = getROICount(ROI_id, s);
    double vol = getROIVol(ROI_id);
    return count/ (1.0e3 * vol * steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::setROIClamped(std::string ROI_id, std::string const & s, bool b)
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
            
            smtos::Tri * tri = pTris[tidx];
            uint slidx = tri->patchdef()->specG2L(sgidx);
            if (slidx == ssolver::LIDX_UNDEFINED)
            {
                spec_undefined << tidx << " ";
                has_spec_warning = true;
                continue;
            }
            if (tri->getInHost()) tri->setClamped(slidx, b);
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
            
            smtos::Tet * tet = pTets[tidx];
            uint slidx = tet->compdef()->specG2L(sgidx);
            if (slidx == ssolver::LIDX_UNDEFINED)
            {
                spec_undefined << tidx << " ";
                has_spec_warning = true;
                continue;
            }
            
            if (tet->getInHost()) tet->setClamped(slidx, b);
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

void smtos::TetOpSplitP::setROIReacK(std::string ROI_id, std::string const & r, double kf)
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
        
        smtos::Tet * tet = pTets[tidx];
        uint rlidx = tet->compdef()->reacG2L(rgidx);
        if (rlidx == ssolver::LIDX_UNDEFINED)
        {
            reac_undefined << tidx << " ";
            has_reac_warning = true;
            continue;
        }
        
        if (tet->getInHost()) tet->reac(rlidx)->setKcst(kf);
    }
    
    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, no change is applied to them:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }
    
    if (has_reac_warning) {
        std::cerr << "Warning: Reac " << r << " has not been defined in the following tetrahedrons, no change is applied to them:\n";
        std::cerr << reac_undefined.str() << "\n";
    }
    
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::setROISReacK(std::string ROI_id, std::string const & sr, double kf)
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
        
        smtos::Tri * tri = pTris[tidx];
        uint srlidx = tri->patchdef()->sreacG2L(srgidx);
        if (srlidx == ssolver::LIDX_UNDEFINED)
        {
            sreac_undefined << tidx << " ";
            has_sreac_warning = true;
            continue;
        }
        
        if (tri->getInHost()) tri->sreac(srlidx)->setKcst(kf);
    }
    
    if (has_tri_warning) {
        std::cerr << "Warning: The following triangles have not been assigned to a patch, no change is applied to them:\n";
        std::cerr << tri_not_assign.str() << "\n";
    }
    
    if (has_sreac_warning) {
        std::cerr << "Warning: SReac " << sr << " has not been defined in the following patch, no change is applied to them:\n";
        std::cerr << sreac_undefined.str() << "\n";
    }
    
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::setROIDiffD(std::string ROI_id, std::string const & d, double dk)
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
        
        smtos::Tet * tet = pTets[tidx];
        uint dlidx = tet->compdef()->diffG2L(dgidx);
        if (dlidx == ssolver::LIDX_UNDEFINED)
        {
            diff_undefined << tidx << " ";
            has_diff_warning = true;
            continue;
        }
        
        if (tet->getInHost()) tet->diff(dlidx)->setDcst(dk);
    }
    
    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, no change is applied to them:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }
    
    if (has_diff_warning) {
        std::cerr << "Warning: Diff " << d << " has not been defined in the following tetrahedrons, no change is applied to them:\n";
        std::cerr << diff_undefined.str() << "\n";
    }

    recomputeUpdPeriod = true;

    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::setROIReacActive(std::string ROI_id, std::string const & r, bool a)
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
        
        smtos::Tet * tet = pTets[tidx];
        uint rlidx = tet->compdef()->reacG2L(rgidx);
        if (rlidx == ssolver::LIDX_UNDEFINED)
        {
            reac_undefined << tidx << " ";
            has_reac_warning = true;
            continue;
        }
        
        if (tet->getInHost()) tet->reac(rlidx)->setActive(a);
    }
    
    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, no change is applied to them:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }
    
    if (has_reac_warning) {
        std::cerr << "Warning: Reac " << r << " has not been defined in the following tetrahedrons, no change is applied to them:\n";
        std::cerr << reac_undefined.str() << "\n";
    }
    
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::setROISReacActive(std::string ROI_id, std::string const & sr, bool a)
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
        
        smtos::Tri * tri = pTris[tidx];
        uint srlidx = tri->patchdef()->sreacG2L(srgidx);
        if (srlidx == ssolver::LIDX_UNDEFINED)
        {
            sreac_undefined << tidx << " ";
            has_sreac_warning = true;
            continue;
        }
        
        if (tri->getInHost()) tri->sreac(srlidx)->setActive(a);
    }
    
    if (has_tri_warning) {
        std::cerr << "Warning: The following triangles have not been assigned to a patch, no change is applied to them:\n";
        std::cerr << tri_not_assign.str() << "\n";
    }
    
    if (has_sreac_warning) {
        std::cerr << "Warning: SReac " << sr << " has not been defined in the following patch, no change is applied to them:\n";
        std::cerr << sreac_undefined.str() << "\n";
    }
    
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::setROIDiffActive(std::string ROI_id, std::string const & d, bool a)
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
        
        smtos::Tet * tet = pTets[tidx];
        uint dlidx = tet->compdef()->diffG2L(dgidx);
        if (dlidx == ssolver::LIDX_UNDEFINED)
        {
            diff_undefined << tidx << " ";
            has_diff_warning = true;
            continue;
        }
        
        if (tet->getInHost()) tet->diff(dlidx)->setActive(a);
    }
    
    if (has_tet_warning) {
        std::cerr << "Warning: The following tetrahedrons have not been assigned to a compartment, no change is applied to them:\n";
        std::cerr << tet_not_assign.str() << "\n";
    }
    
    if (has_diff_warning) {
        std::cerr << "Warning: Diff " << d << " has not been defined in the following tetrahedrons, no change is applied to them:\n";
        std::cerr << diff_undefined.str() << "\n";
    }
    
    recomputeUpdPeriod = true;

    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::setROIVDepSReacActive(std::string ROI_id, std::string const & vsr, bool a)
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
        
        smtos::Tri * tri = pTris[tidx];
        uint vsrlidx = tri->patchdef()->vdepsreacG2L(vsrgidx);
        if (vsrlidx == ssolver::LIDX_UNDEFINED)
        {
            vsreac_undefined << tidx << " ";
            has_vsreac_warning = true;
            continue;
        }
        
        if (tri->getInHost()) tri->vdepsreac(vsrlidx)->setActive(a);
    }
    
    if (has_tri_warning) {
        std::cerr << "Warning: The following triangles have not been assigned to a patch, no change is applied to them:\n";
        std::cerr << tri_not_assign.str() << "\n";
    }
    
    if (has_vsreac_warning) {
        std::cerr << "Warning: VDepSReac " << vsr << " has not been defined in the following patch, no change is applied to them:\n";
        std::cerr << vsreac_undefined.str() << "\n";
    }
    _updateLocal();
}

////////////////////////////////////////////////////////////////////////////////

uint smtos::TetOpSplitP::getROIReacExtent(std::string ROI_id, std::string const & r) const
{
    std::ostringstream os;
    os << "This function has not been implemented!";
    throw steps::NotImplErr(os.str());
    
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
        
        smtos::Tet * tet = pTets[tidx];
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

void smtos::TetOpSplitP::resetROIReacExtent(std::string ROI_id, std::string const & r)
{
    std::ostringstream os;
    os << "This function has not been implemented!";
    throw steps::NotImplErr(os.str());
    
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
        
        smtos::Tet * tet = pTets[tidx];
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

uint smtos::TetOpSplitP::getROISReacExtent(std::string ROI_id, std::string const & sr) const
{
    std::ostringstream os;
    os << "This function has not been implemented!";
    throw steps::NotImplErr(os.str());
    
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
        
        smtos::Tri * tri = pTris[tidx];
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

void smtos::TetOpSplitP::resetROISReacExtent(std::string ROI_id, std::string const & sr)
{
    std::ostringstream os;
    os << "This function has not been implemented!";
    throw steps::NotImplErr(os.str());
    
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
        
        smtos::Tri * tri = pTris[tidx];
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

uint smtos::TetOpSplitP::getROIDiffExtent(std::string ROI_id, std::string const & d) const
{
    std::ostringstream os;
    os << "This function has not been implemented!";
    throw steps::NotImplErr(os.str());
    
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
        
        smtos::Tet * tet = pTets[tidx];
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

void smtos::TetOpSplitP::resetROIDiffExtent(std::string ROI_id, std::string const & d)
{
    std::ostringstream os;
    os << "This function has not been implemented!";
    throw steps::NotImplErr(os.str());
    
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
        
        smtos::Tet * tet = pTets[tidx];
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

uint smtos::TetOpSplitP::getTetHostRank(uint tidx)
{
    return tetHosts[tidx];
}

////////////////////////////////////////////////////////////////////////////////

uint smtos::TetOpSplitP::getTriHostRank(uint tidx)
{
    return triHosts[tidx];
}

////////////////////////////////////////////////////////////////////////////////

uint smtos::TetOpSplitP::getWMVolHostRank(uint idx)
{
    return wmHosts[idx];
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_updateLocal(uint* upd_entries, uint buffer_size) {
    for (uint i = 0; i < buffer_size; i++) {
        if (pKProcs[upd_entries[i]] != NULL)
        {
            _updateElement(pKProcs[upd_entries[i]]);
        }
    }
    
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::addDiff(Diff* diff)
{
    diff->crData.pos = pDiffs.size();
    pDiffs.push_back(diff);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_updateDiff(Diff* diff)
{
    double new_rate = diff->rate(this);

    diff->crData.rate = new_rate;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::addSDiff(SDiff* sdiff)
{
    sdiff->crData.pos = pSDiffs.size();
    pSDiffs.push_back(sdiff);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_updateSDiff(SDiff* sdiff)
{
    double new_rate = sdiff->rate(this);

    sdiff->crData.rate = new_rate;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::addNeighHost(int host)
{
    neighbHosts.insert(host);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::registerBoundaryTet(smtos::Tet* tet)
{
    boundaryTets.insert(tet);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::registerBoundaryTri(smtos::Tri* tri)
{
    boundaryTris.insert(tri);
}

////////////////////////////////////////////////////////////////////////////////

uint smtos::TetOpSplitP::registerRemoteMoleculeChange(int svol_host, uint loc, SubVolType svol_type, uint idx, uint slidx, uint change)
{
    uint new_loc = remoteChanges[svol_host].size();

    if (new_loc == 0 || new_loc - 4 < loc) {
        remoteChanges[svol_host].push_back(svol_type);
        remoteChanges[svol_host].push_back(idx);
        remoteChanges[svol_host].push_back(slidx);
        remoteChanges[svol_host].push_back(change);
    }
    else {
        uint stored_type = remoteChanges[svol_host][loc];
        uint stored_idx = remoteChanges[svol_host][loc + 1];
        uint stored_slidx = remoteChanges[svol_host][loc + 2];

        if (stored_type == svol_type && stored_idx ==idx && stored_slidx == slidx) {
            remoteChanges[svol_host][loc + 3] += change;
            new_loc = loc;
        }
        else {
            remoteChanges[svol_host].push_back(svol_type);
            remoteChanges[svol_host].push_back(idx);
            remoteChanges[svol_host].push_back(slidx);
            remoteChanges[svol_host].push_back(change);
        }
    }
    return new_loc;
}

////////////////////////////////////////////////////////////////////////////////
/*
void smtos::TetOpSplitP::registerSyncWmVol(smtos::WmVol * wmvol)
{
    #ifdef MPI_DEBUG
    //CLOG(DEBUG, "mpi_debug") << "register sync tet " << wmvol->idx()<< "\n";
    #endif
    syncWmVols.insert(wmvol);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::registerSyncTet(smtos::Tet * tet)
{
    #ifdef MPI_DEBUG
    //CLOG(DEBUG, "mpi_debug") << "register sync tet " << tet->idx()<< "\n";
    #endif
    syncTets.insert(tet);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::registerSyncTri(smtos::Tri * tri)
{
    #ifdef MPI_DEBUG
    //CLOG(DEBUG, "mpi_debug") << "register sync tri " << tri->idx()<< "\n";
    #endif
    syncTris.insert(tri);
}
*/


////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP:: _remoteSyncAndUpdate(void* requests, std::vector<KProc*> & applied_diffs, std::vector<int> & directions)
{
    #ifdef MPI_DEBUG
    CLOG(DEBUG, "mpi_debug") << "Start applying molecule changes.\n";
    #endif
    
    #ifdef MPI_PROFILING
    double starttime = MPI_Wtime();
    #endif
    
    #ifdef MPI_DEBUG
    CLOG(DEBUG, "mpi_debug") << "change buffer:" << remoteChanges;
    #endif
    
    MPI_Request* requestsPtr = static_cast<MPI_Request*>(requests);
    
    uint request_count = 0;
    for (auto dest : neighbHosts) {
#ifdef MPI_DEBUG
        CLOG(DEBUG, "mpi_debug") << "request pointer:" << &(requestsPtr[request_count]);
#endif
        MPI_Isend(&remoteChanges[dest].front(), remoteChanges[dest].size(), MPI_UNSIGNED, dest, OPSPLIT_MOLECULE_CHANGE, MPI_COMM_WORLD, &(requestsPtr[request_count]));
        request_count ++;
    }
    
    uint remain_comfirm = nNeighbHosts;
    MPI_Status status;
    std::set<KProc*> upd_kprocs;
    
    std::set<int> await_neighbors(neighbHosts);
    
    #ifdef MPI_PROFILING
    double endtime = MPI_Wtime();
    syncTime += (endtime - starttime);
    #endif
    
    while (!await_neighbors.empty()) {
        #ifdef MPI_PROFILING
        starttime = MPI_Wtime();
        #endif
        int flag = 0;
        int data_source = 0;
        for (auto neighbor : await_neighbors) {
            MPI_Iprobe(neighbor, OPSPLIT_MOLECULE_CHANGE, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                data_source = neighbor;
                break;
            }
        }
        #ifdef MPI_PROFILING
        endtime = MPI_Wtime();
        idleTime += (endtime - starttime);
        #endif
        if (!flag) continue;

        #ifdef MPI_PROFILING
        starttime = MPI_Wtime();
        #endif
        // receive data
        int change_size = 0;
        MPI_Get_count(&status, MPI_UNSIGNED, &change_size);
        std::vector<uint> changes(change_size);
        MPI_Recv(&changes.front(), change_size, MPI_UNSIGNED, status.MPI_SOURCE, OPSPLIT_MOLECULE_CHANGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        #ifdef MPI_DEBUG
        CLOG(DEBUG, "mpi_debug") << "Recving change "<< changes << " from " << status.MPI_SOURCE << ".\n";
        #endif
        
        // apply changes
        uint nchanges = change_size / 4;
        for (uint c = 0; c < nchanges; c++) {
            uint type = changes[c * 4];
            uint idx = changes[c * 4 + 1];
            uint slidx = changes[c * 4 + 2];
            uint value = changes[c * 4 + 3];
            
            #ifdef MPI_DEBUG
            CLOG(DEBUG, "mpi_debug") << "unpack: idx "<< idx << " slidx " << slidx << " value " << value << ".\n";
            #endif
            
            if (type == SUB_WM) {
                pWmVols[idx]->incCount(slidx, value);
            }
            if (type == SUB_TET) {
                pTets[idx]->incCount(slidx, value);
                std::vector<smtos::KProc*> const & remote_upd = pTets[idx]->getSpecUpdKProcs(slidx);
                upd_kprocs.insert(remote_upd.begin(), remote_upd.end());
            }
            if (type == SUB_TRI) {
                pTris[idx]->incCount(slidx, value);
                std::vector<smtos::KProc*> const & remote_upd = pTris[idx]->getSpecUpdKProcs(slidx);
                upd_kprocs.insert(remote_upd.begin(), remote_upd.end());
            }
        }
        
        await_neighbors.erase(data_source);
        #ifdef MPI_PROFILING
        endtime = MPI_Wtime();
        syncTime += (endtime - starttime);
        #endif
    }
    
    #ifdef MPI_PROFILING
    starttime = MPI_Wtime();
    #endif
    
    #ifdef MPI_DEBUG
    CLOG(DEBUG, "mpi_debug") << "Molecule changes have been applied.\n";
    #endif
    
    uint napply = applied_diffs.size();
    for (uint i = 0; i < napply; i++) {
        KProc* kp = applied_diffs[i];
        int direction = directions[i];

        std::vector<smtos::KProc*> const & local_upd = kp->getLocalUpdVec(direction);

        for (auto & upd_kp : local_upd) {
            _updateElement(upd_kp);
        }
    }
    
    // update kprocs caused by remote molecule changes
    for (auto & upd_kp : upd_kprocs) {
        _updateElement(upd_kp);
    }
    _updateSum();
    
    #ifdef MPI_DEBUG
    CLOG(DEBUG, "mpi_debug") << "Remote update completed.\n";
    #endif

    #ifdef MPI_PROFILING
    endtime = MPI_Wtime();
    compTime += (endtime - starttime);
    #endif
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::repartitionAndReset(std::vector<uint> const &tet_hosts,
                 std::map<uint, uint> const &tri_hosts,
                 std::vector<uint> const &wm_hosts)
{
    pKProcs.clear();
    pDiffs.clear();
    pSDiffs.clear();
    neighbHosts.clear();
    boundaryTets.clear();
    boundaryTris.clear();
    
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
    
    tetHosts.assign(tet_hosts.begin(), tet_hosts.end());
    triHosts.clear();
    triHosts.insert(tri_hosts.begin(), tri_hosts.end());
    wmHosts.assign(wm_hosts.begin(), wm_hosts.end());
    
    uint ntets = pTets.size();
    for (uint t = 0; t < ntets; ++t) {
        if (pTets[t] == 0) continue;
        pTets[t]->repartition(this, myRank, tetHosts[t]);
    }
    
    uint nwms = pWmVols.size();
    for (uint wm = 0; wm < nwms; ++wm) {
        if (pWmVols[wm] == 0) continue;
        pWmVols[wm]->repartition(this, myRank, wmHosts[wm]);
    }
    
    for (auto t: pTris) {
        if (t == 0) continue;
        t->repartition(this, myRank, triHosts[t->idx()]);
    }
    
    for (auto t: pTets)
    if (t && t->getInHost()) t->setupDeps();
    
    // Vector allows for all compartments to be well-mixed, so
    // hold null-pointer for mesh compartments
    for (auto wmv: pWmVols)
    if (wmv && wmv->getInHost()) wmv->setupDeps();
    
    // DEBUG: vector holds all possible triangles, but
    // only patch triangles are filled
    for (auto t: pTris)
    if (t && t->getInHost()) t->setupDeps();

    for (auto tet : boundaryTets) {
        tet->setupBufferLocations();
    }
    for (auto tri : boundaryTris) {
        tri->setupBufferLocations();
    }
    
    if (efflag() == true) {
        std::ostringstream os;
        os << "Repartition of EField is not implemented:\n";
        throw steps::ArgErr(os.str());
    }
    
    neighbHosts.erase(myRank);
    nNeighbHosts = neighbHosts.size();
    
    // construct remote molecule change buffers
    remoteChanges.clear();
    for (auto neighbor : neighbHosts) {
        remoteChanges[neighbor] = std::vector<uint> ();
    }
    
    nEntries = pKProcs.size();
    diffSep=pDiffs.size();
    sdiffSep=pSDiffs.size();
    reset();
    MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////
/*
void smtos::TetOpSplitP::_updateKProcRates(std::vector<KProc*> & applylist, std::vector<int> & directions, std::vector<MPI_Request> & requests)
{
    remoteUpdateKProcs.clear();
    uint napply = applylist.size();
    for (uint i = 0; i < napply; i++) {
        KProc* kp = applylist[i];
        int direction = directions[i];
        #ifdef MPI_DEBUG
        CLOG(DEBUG, "mpi_debug") << "direction: " << direction;
        #endif
        std::vector<smtos::KProc*> const & local_upd = kp->getLocalUpdVec(direction);
        #ifdef MPI_DEBUG
        CLOG(DEBUG, "mpi_debug") << "local upd: " << local_upd;
        #endif
        std::vector<uint> const & remote_upd = kp->getRemoteUpdVec(direction);
        #ifdef MPI_DEBUG
        CLOG(DEBUG, "mpi_debug") << "remote upd: " << remote_upd;
        #endif
        for (auto & upd_kp : local_upd) {
#ifdef MPI_DEBUG
            CLOG(DEBUG, "mpi_debug") << "local upd: " << upd_kp->schedIDX();
#endif
            _updateElement(upd_kp);
        }
        remoteUpdateKProcs.insert(remote_upd.begin(), remote_upd.end());
    }
    
    _updateRemoteKProcRates(requests);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::_updateRemoteKProcRates(std::vector<MPI_Request> & requests)
{
    #ifdef MPI_DEBUG
    CLOG(DEBUG, "mpi_debug") << "Recvive remote update\n";
    #endif

    std::vector<uint> remote_upd_kps(remoteUpdateKProcs.begin(), remoteUpdateKProcs.end());
    
    uint request_count = nNeighbHosts * 2;
    
    for (int dest : neighbHosts) {
        MPI_Isend(&remote_upd_kps.front(), remote_upd_kps.size(), MPI_UNSIGNED, dest, OPSPLIT_KPROC_UPD, MPI_COMM_WORLD, &(requests[request_count]));
        request_count++;
    }
    
    // wait for either remote changes or a complete confirm
    uint remain_comfirm = nNeighbHosts;
    MPI_Status status;
    while (remain_comfirm != 0) {
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == OPSPLIT_KPROC_UPD) {
            // receive data
            int data_size = 0;
            MPI_Get_count(&status, MPI_UNSIGNED, &data_size);
            std::vector<uint> data(data_size);
            MPI_Recv(&data.front(), data_size, MPI_UNSIGNED, status.MPI_SOURCE, OPSPLIT_KPROC_UPD, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // apply changes
            for (auto & kpidx : data) {
                smtos::KProc* kp = pKProcs[kpidx];
                if (kp != NULL) _updateElement(kp);
            }
            remain_comfirm --;
        }
        else {
        #ifdef MPI_DEBUG
            CLOG(DEBUG, "mpi_debug") << "Recvive Unknown signal: " << status.MPI_TAG << "\n";
        #endif
            throw;
        }
    }
    remoteUpdateKProcs.clear();
}
*/
////////////////////////////////////////////////////////////////////////////////

void smtos::TetOpSplitP::setDiffApplyThreshold(int threshold)
{
    diffApplyThreshold = threshold;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::getReacExtent(bool local)
{
    if (local) {
        return reacExtent;
    }
    
    double sum;
    MPI_Reduce(&reacExtent, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::getDiffExtent(bool local)
{
    if (local) {
        return diffExtent;
    }
    double sum;
    MPI_Reduce(&diffExtent, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::getNIteration(void)
{
    return nIteration;
}


////////////////////////////////////////////////////////////////////////////////


double smtos::TetOpSplitP::getCompTime(void)
{
    return compTime;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::getSyncTime(void)
{
    return syncTime;
}

////////////////////////////////////////////////////////////////////////////////
double smtos::TetOpSplitP::getIdleTime(void)
{
    return idleTime;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::getEFieldTime(void)
{
    return efieldTime;
}
////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::getRDTime(void)
{
    return rdTime;
}
////////////////////////////////////////////////////////////////////////////////

double smtos::TetOpSplitP::getDataExchangeTime(void)
{
    return dataExchangeTime;
}
////////////////////////////////////////////////////////////////////////////////
// END

