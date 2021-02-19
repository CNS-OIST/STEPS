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


#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>

#include <cvode/cvode.h>                 /* prototypes for CVODE fcts., consts. */
#include <cvode/cvode_dense.h>          /* prototype for CVDense */
#include <nvector/nvector_serial.h>      /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_dense.h>     /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>     /* definition of type realtype */

#include "steps/tetode/tetode.hpp"

#include "steps/tetode/comp.hpp"
#include "steps/tetode/patch.hpp"
#include "steps/tetode/tet.hpp"
#include "steps/tetode/tri.hpp"

#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/math/constants.hpp"
#include "steps/math/point.hpp"

#include "steps/solver/compdef.hpp"
#include "steps/solver/patchdef.hpp"
#include "steps/solver/reacdef.hpp"
#include "steps/solver/sreacdef.hpp"
#include "steps/solver/vdepsreacdef.hpp"

#include "steps/geom/tetmesh.hpp"

#include "steps/error.hpp"

#include "steps/solver/efield/dVsolver.hpp"
#include "steps/solver/efield/efield.hpp"

// logging
#include <easylogging++.h>

// CVODE definitions
#define Ith(v,i)    NV_Ith_S(v,i)
#define IJth(A,i,j) DENSE_ELEM(A,i,j) // IJth numbers rows,cols 1..NEQ

////////////////////////////////////////////////////////////////////////////////

// A vector all the reaction information, hopefully ingeniously
// removing the need for a sparse matrix at all
// This is stored in no-man's land because CVode function (f_cvode) needs access
// to it and it wasn't possible to use as a class member, even as a friend
// function (function needs access to Tetode pointer, which means re-writing the
// function definition from the c side- messy)
std::vector<std::vector<steps::tetode::structA> >  pSpec_matrixsub;

////////////////////////////////////////////////////////////////////////////////

using steps::math::point3d;

using steps::solver::LIDX_UNDEFINED;
using steps::math::AVOGADRO;

////////////////////////////////////////////////////////////////////////////////

static int f_cvode(realtype /*t*/, N_Vector y, N_Vector ydot, void */*user_data*/)
{
  uint i = 0;
  //for (uint i = 0; i < tetode->pSpecs_tot; ++i)
  for (auto const& sp: pSpec_matrixsub) {
    double dydt = 0.0;
    for (auto const& r: sp) {
      double dydt_r = r.upd * r.ccst;
      for (auto const& p: r.players) {
        for (auto const& q: p.info) {
          double val = Ith(y,q.spec_idx);
          uint order = q.order;
          if (order == 1) dydt_r*=val;
          else dydt_r*=pow(val, order);
        }
      }
      dydt+=dydt_r;
    }
    // Update the ydot vector with the calculated dydt
    Ith(ydot, i) = dydt;
    ++i;
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace tetode {

 // CVODE stuff
 struct CVodeState {
   // count (size of y_cvode, abstol_cvode vectors)
   uint N;
   // max CVode steps
   uint Nmax_cvode;
   // relative tolerance
   realtype reltol_cvode;
   // Absolute tolerance vector
   N_Vector abstol_cvode;
   // Vector of values, for us species counts
   N_Vector y_cvode;
   // Memory block for CVODE
   void     * cvode_mem_cvode;

   CVodeState(uint N_, uint maxn, double atol, double rtol);
   ~CVodeState();

   void setTolerances(double atol, double rtol);
   void setMaxNumSteps(uint maxn);
   int  initialise();
   int  reinit(realtype starttime);

   int  run(realtype endtime);

   void checkpoint(std::fstream &);
   void restore(std::fstream &);
 };

void check_flag(void *flagvalue, const char *funcname, int opt)
{
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == nullptr)
    {
        std::ostringstream os;
        os << "\nSUNDIALS_ERROR: " << funcname << "() failed - returned NULL pointer\n\n";
        SysErrLog(os.str());
    }

    /* Check if flag < 0 */
    else if (opt == 1)
    {
        errflag = static_cast<int *>(flagvalue);
        if (*errflag < 0)
        {
            std::ostringstream os;
            os << "\nSUNDIALS_ERROR: " << funcname
               << "() failed with flag = " << *errflag
               << "\n\n";
            SysErrLog(os.str());
        }
     }
}

////////////////////////////////////////////////////////////////////////////////

CVodeState::CVodeState(uint N_, uint maxn, realtype atol, realtype rtol) {
    N = N_;
    Nmax_cvode = maxn;

    // Creates serial vectors for y and absolute tolerances
    y_cvode = N_VNew_Serial(N);
    check_flag(y_cvode, "N_VNew_Serial", 0);

    abstol_cvode = N_VNew_Serial(N);
    check_flag(abstol_cvode, "N_VNew_Serial", 0);

    reltol_cvode = rtol;

    for (uint i=0; i<N; ++i)
    {
        Ith(abstol_cvode, i) = atol;
    }

    // Call CVodeCreate to create the solver memory and specify the
    // Backward Differentiation Formula and the use of a Newton iteration
    //cvode_mem_cvode = CVodeCreate(CV_BDF, CV_NEWTON);
    // NO- the above choice eats up memory like you wouldn't believe
    // and causes segmentation faults

    // ADAMS and FUNCTIONAL are a much much better choice
    cvode_mem_cvode = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL) ;

    check_flag(cvode_mem_cvode, "CVodeCreate", 0);

    // Initialise y:
    for (uint i=0; i<N; ++i)
    {
        Ith(y_cvode, i) = 0.0;
    }

    // Call CVodeInit to initialize the integrator memory and specify the
    // user's right hand side function in y'=f(t,y), the initial time T0, and
    // the initial dependent variable vector y.
    // For the reason that CVode basically doesn't expect changes to be made
    // during any given run (such as injection of molecules) at the moment
    // such features will not be supported. To support them will mean
    // creating and freeing memory, copying structures etc and could be quite tricky
    int flag = CVodeInit(cvode_mem_cvode, f_cvode, 0.0, y_cvode);
    check_flag(&flag, "CVodeInit", 1);
}

CVodeState::~CVodeState() {
    N_VDestroy_Serial(y_cvode);
    N_VDestroy_Serial(abstol_cvode);

    /* Free integrator memory */
    CVodeFree(&cvode_mem_cvode);
}

////////////////////////////////////////////////////////////////////////////////

void CVodeState::checkpoint(std::fstream &cp_file) {
    cp_file.write(reinterpret_cast<char*>(&Nmax_cvode), sizeof(uint));
    cp_file.write(reinterpret_cast<char*>(&reltol_cvode), sizeof(realtype));
    {
        auto content = static_cast<N_VectorContent_Serial>(abstol_cvode->content);
        cp_file.write(reinterpret_cast<char *>(content->data), sizeof(realtype) * N);
    }
    {
        auto content = static_cast<N_VectorContent_Serial>(y_cvode->content);
        cp_file.write(reinterpret_cast<char*>(content->data), sizeof(realtype) * N);
    }

}

////////////////////////////////////////////////////////////////////////////////

void CVodeState::restore(std::fstream &cp_file) {
    cp_file.read(reinterpret_cast<char*>(&Nmax_cvode), sizeof(uint));
    cp_file.read(reinterpret_cast<char*>(&reltol_cvode), sizeof(realtype));
    {
        auto content = static_cast<N_VectorContent_Serial>(abstol_cvode->content);
        cp_file.read(reinterpret_cast<char*>(content->data), sizeof(realtype) * N);
    }
    {
        auto content = static_cast<N_VectorContent_Serial>(y_cvode->content);
        cp_file.read(reinterpret_cast<char*>(content->data), sizeof(realtype) * N);
    }

}

////////////////////////////////////////////////////////////////////////////////

void CVodeState::setTolerances(double atol, double rtol) {
    // I suppose they shouldn't be negative
    if (atol < 0.0 or rtol < 0.0)
    {
        std::ostringstream os;
        os << "Neither absolute tolerance nor relative tolerance should ";
        os << "be negative.\n";
        ArgErrLog(os.str());
    }

    reltol_cvode = rtol;

    for (uint i=0; i<N; ++i)
    {
        Ith(abstol_cvode, i) = atol;
    }
}

////////////////////////////////////////////////////////////////////////////////

void CVodeState::setMaxNumSteps(uint maxn) {
    int flag = CVodeSetMaxNumSteps(cvode_mem_cvode, maxn);
    check_flag(&flag, "CVodeSetMaxNumSteps", 1);

    Nmax_cvode = maxn;
}

////////////////////////////////////////////////////////////////////////////////

int CVodeState::initialise() {
    int flag;

    flag = CVodeSetMaxNumSteps(cvode_mem_cvode, Nmax_cvode);
    check_flag(&flag, "CVodeSetMaxNumSteps", 1);

    // Call CVodeSVtolerances to specify the scalar relative tolerance
    // and vector absolute tolerances */
    flag = CVodeSVtolerances(cvode_mem_cvode, reltol_cvode, abstol_cvode);
    check_flag(&flag, "CVodeSVtolerances", 1);

    // Call CVDense to specify the CVDENSE dense linear solver */
    //flag = CVDense(cvode_mem_cvode, pSpecs_tot);
    //check_flag(&flag, "CVDense", 1);

    return flag;
}

int CVodeState::reinit(realtype starttime) {
    int flag = CVodeReInit(cvode_mem_cvode, starttime, y_cvode);
    check_flag(&flag, "CVodeInit", 1);

    return flag;
}

int CVodeState::run(realtype endtime) {
    realtype t;
    return CVode(cvode_mem_cvode, endtime, y_cvode, &t, CV_NORMAL);
}

////////////////////////////////////////////////////////////////////////////////

TetODE::TetODE(steps::model::Model *m, steps::wm::Geom *g, const rng::RNGptr &r,
               int calcMembPot)
: API(m, g, r)
, pEFoption(static_cast<EF_solver>(calcMembPot))
{
    _setup();
}

////////////////////////////////////////////////////////////////////////////////

TetODE::~TetODE()
{
    for (auto const& c: pComps) delete c;
    for (auto const& p: pPatches) delete p;

    for (auto const& t: pTets) {
        delete t;
    }

    for (auto const& t: pTris) {
        delete t;
    }

    delete pCVodeState;

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

////////////////////////////////////////////////////////////////////////////////

std::string TetODE::getSolverName() const
{
    return "tetODE";
}

////////////////////////////////////////////////////////////////////////////////

std::string TetODE::getSolverDesc() const
{
    return "Reaction-diffusion ODE solver in tetrahedral mesh";
}

////////////////////////////////////////////////////////////////////////////////

std::string TetODE::getSolverAuthors() const
{
    return "Iain Hepburn";
}

////////////////////////////////////////////////////////////////////////////////

std::string TetODE::getSolverEmail() const
{
    return "steps.dev@gmail.com";
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::checkpoint(std::string const & file_name)
{
    std::fstream cp_file;

    cp_file.open(file_name.c_str(),
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    statedef().checkpoint(cp_file);

    for (auto &pComp : pComps) {
        pComp->checkpoint(cp_file);
    }

    for (auto const& pPatch: pPatches) {
        pPatch->checkpoint(cp_file);
    }

    for (auto const& pTri: pTris) {
        pTri->checkpoint(cp_file);
    }

    for (auto &pTet : pTets) {
        pTet->checkpoint(cp_file);
    }

    pCVodeState->checkpoint(cp_file);

    if (efflag()) {
        cp_file.write(reinterpret_cast<char*>(&pTemp), sizeof(double));
        cp_file.write(reinterpret_cast<char*>(&pEFDT), sizeof(double));
        pEField->checkpoint(cp_file);
    }

    cp_file.close();
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::restore(std::string const & file_name)
{
    std::fstream cp_file;

    cp_file.open(file_name.c_str(),
                 std::fstream::in | std::fstream::binary);

    cp_file.seekg(0);

    statedef().restore(cp_file);

    for (auto &pComp : pComps) {
        pComp->restore(cp_file);
    }

    for (auto &pPatch : pPatches) {
        pPatch->restore(cp_file);
    }

    for (auto &pTri : pTris) {
        pTri->restore(cp_file);
    }

    for (auto & pTet : pTets) {
        pTet->restore(cp_file);
    }

    pCVodeState->restore(cp_file);

    if (efflag()) {
        cp_file.read(reinterpret_cast<char*>(&pTemp), sizeof(double));
        cp_file.read(reinterpret_cast<char*>(&pEFDT), sizeof(double));
        pEField->restore(cp_file);
    }

    cp_file.close();

    pTolsset = true;
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::reset()
{
    std::ostringstream os;
    os << "reset() not implemented for steps::solver::TetODE solver";
    NotImplErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::setMaxNumSteps(uint maxn)
{
    pCVodeState->setMaxNumSteps(maxn);
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::getTime() const
{
    return statedef().time();
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setup()
{
    // Perform upcast.
    pMesh = dynamic_cast<steps::tetmesh::Tetmesh *>(geom());
    if (!pMesh)
        ArgErrLog("Geometry description to steps::solver::TetODE solver "
                "constructor is not a valid steps::tetmesh::Tetmesh object.");

    const auto ntets = pMesh->countTets();
    const auto ntris = pMesh->countTris();
    const auto ncomps = pMesh->_countComps();
    const auto npatches = pMesh->_countPatches();

    AssertLog(npatches == statedef().countPatches());
    AssertLog(ncomps == statedef().countComps());

    // Now create the actual compartments.
    for (auto const& c : statedef().comps()) {
        const auto compdef_gidx = c->gidx();
        const auto comp_idx = _addComp(c);
        AssertLog(compdef_gidx == comp_idx);
    }
    // Create the actual patches.
    for (auto const& p : statedef().patches()) {
        const auto patchdef_gidx = p->gidx();
        const auto patch_idx = _addPatch(p);
        AssertLog(patchdef_gidx == patch_idx);
    }

    AssertLog(pPatches.size() == npatches);
    AssertLog(pComps.size() == ncomps);

    pTets.assign(ntets, nullptr);
    pTris.assign(ntris, nullptr);

    for (uint p = 0; p < npatches; ++p)
    {
        // Now add the tris for this patch
        steps::wm::Patch * wmpatch = pMesh->_getPatch(p);

        // sanity check
        AssertLog(statedef().getPatchIdx(wmpatch) == p );

        // Perform upcast
        auto tmpatch = dynamic_cast<steps::tetmesh::TmPatch*>(wmpatch);
        if (!tmpatch)
            ArgErrLog("Well-mixed patches not supported in steps::solver::TetODE solver.");

        steps::tetode::Patch *localpatch = pPatches[p];
        std::map<bar_id_t, std::vector<triangle_id_t> > bar2tri;
        for (auto tri: tmpatch->_getAllTriIndices())
        {
            const auto& bars = pMesh->_getTriBars(tri);
            for (int i = 0; i < 3; ++i)
                bar2tri[bars[i]].push_back(tri);
        }

        for (auto tri: tmpatch->_getAllTriIndices())
        {
            AssertLog(pMesh->getTriPatch(tri) == tmpatch);

            double area = pMesh->getTriArea(tri);

            // NB: Tri vertices may not be in consistent order, so use bar interface.
            const auto tri_bars = pMesh->_getTriBars(tri);
            double l[3] = {0, 0, 0};

            for (uint j=0; j<3; ++j) {
                const auto v = pMesh->_getBar(tri_bars[j]);
                l[j] = distance(pMesh->_getVertex(v[0]), pMesh->_getVertex(v[1]));
            }

            std::vector<triangle_id_t> tris(3, UNKNOWN_TRI);
            for (int j = 0; j < 3; ++j)
            {
                const auto& neighb_tris = bar2tri[tri_bars[j]];
                for (auto neighb_tri : neighb_tris) {
                    if (neighb_tri == tri || pMesh->getTriPatch(neighb_tri) != tmpatch)
                        continue;
                    tris[j] = neighb_tri;
                    break;
                }
            }
            point3d baryc = pMesh->_getTriBarycenter(tri);

            double d[3] = {0, 0, 0};
            for (uint j = 0; j < 3; ++j) {
                if (tris[j] == UNKNOWN_TRI) continue;
                d[j] = distance(baryc, pMesh->_getTriBarycenter(tris[j]));
            }

            const auto tri_tets = pMesh->_getTriTetNeighb(tri);
            _addTri(tri, localpatch, area, l[0], l[1], l[2], d[0], d[1], d[2], tri_tets[0], tri_tets[1], tris[0], tris[1], tris[2]);
        }
    }

    for (uint c = 0; c < ncomps; ++c)
    {
        // Now add the tets for this comp
        steps::wm::Comp * wmcomp = pMesh->_getComp(c);

        // sanity check
        AssertLog(statedef().getCompIdx(wmcomp) == c);

        // Perform upcast
        auto *tmcomp = dynamic_cast<steps::tetmesh::TmComp*>(wmcomp);
        if (!tmcomp)
            ArgErrLog("Well-mixed compartments not supported in steps::solver::TetODE solver.");
        
        steps::tetode::Comp *localcomp = pComps[c];
        for (auto tet: tmcomp->_getAllTetIndices())
        {
            AssertLog(pMesh->getTetComp(tet) == tmcomp);

            double vol = pMesh->getTetVol(tet);

            const auto *tris = pMesh->_getTetTriNeighb(tet);

            double a[4] = {0, 0, 0, 0};
            for (uint j = 0; j < 4; ++j) {
                a[j] = pMesh->getTriArea(tris[j]);
            }

            const auto tets = pMesh->_getTetTetNeighb(tet);
            point3d baryc = pMesh->_getTetBarycenter(tet);

            double d[4] = {0, 0, 0, 0};
            for (uint j = 0; j < 4; ++j) {
                if (tets[j] == UNKNOWN_TET) continue;
                d[j] = distance(baryc, pMesh->_getTetBarycenter(tets[j]));
            }

            _addTet(tet, localcomp, vol, a[0], a[1], a[2], a[3], d[0], d[1], d[2], d[3],
                    tets[0], tets[1], tets[2], tets[3]);
        }
    }


    // All tets and tris that belong to some comp or patch have been created
    // locally- now we can connect them locally
    // NOTE: currently if a tetrahedron's neighbour belongs to a different
    // comp they do not talk to each other (see stex::Tet::setNextTet())
    //

    AssertLog(ntets == pTets.size());
    // pTets member size of all tets in geometry, but may not be filled with
    // local tets if they have not been added to a compartment
    for (uint t = 0; t < ntets; ++t)
    {
        if (pTets[t] == nullptr) continue;

        for (uint j = 0; j < 4; ++j) {
            auto tet = pTets[t]->tet(j);
            if (tet != UNKNOWN_TET && pTets[tet.get()] != nullptr) pTets[t]->setNextTet(j, pTets[tet.get()]);
        }
        // Not setting Tet triangles at this point- only want to set
        // for surface triangles
    }
    AssertLog(ntris == pTris.size());

    for (uint t = 0; t < ntris; ++t)
    {
        // Looping over all possible tris, but only some have been added to a patch
        if (pTris[t] == nullptr) continue;

        for (uint j = 0; j < 3; ++j) {
            auto tri = pTris[t]->tri(j);
            if (tri != UNKNOWN_TRI && pTris[tri.get()] != nullptr) pTris[t]->setNextTri(j, pTris[tri.get()]);
        }

        // By convention, triangles in a patch should have an inner tetrahedron defined
        // (neighbouring tets 'flipped' if necessary in Tetmesh)
        // but not necessarily an outer tet
        // 17/3/10- actually this is not the case any more with well-mixed compartments
        //
        auto tetinner = pTris[t]->tet(0);
        auto tetouter = pTris[t]->tet(1);

        // Now correct check, previously didn't allow for tet index == 0
        AssertLog(tetinner != UNKNOWN_TET);
        AssertLog(pTets[tetinner.get()] != nullptr);


        if (pTets[tetinner.get()] != nullptr)
        {
            // A triangle may already have an inner tet defined as a well-mixed
            // volume, but that should not be the case here:
            AssertLog(pTris[t]->iTet() == nullptr);

            pTris[t]->setInnerTet(pTets[tetinner.get()]);
            // Now add this triangle to inner tet's list of neighbours
            for (uint i=0; i <= 4; ++i)
            {
                // include assert for debugging purposes and remove
                // once this is tested
                AssertLog(i < 4);                                                        //////////
                // check if there is already a neighbouring tet or tri
                // In theory if there is a tri to add, the tet should
                // have less than 4 neighbouring tets added because
                // a neighbouring tet(s) is in a different compartment

                //     Also check tris because in some cases a surface tet
                // may have more than 1 neighbouring tri
                // NOTE: The order here will end up being different to the
                // neighbour order at the Tetmesh level

                Tet * tet_in = pTets[tetinner.get()];
                if (tet_in->nextTet(i) != nullptr) continue;

                if (tet_in->nextTri(i) != nullptr) continue;
                tet_in->setNextTri(i, pTris[t]);
                break;
            }
        }

        if (tetouter != UNKNOWN_TET)
        {
            if (pTets[tetouter.get()] != nullptr)
            {
                // A triangle may already have an inner tet defined as a well-mixed
                // volume, but that should not be the case here:
                AssertLog(pTris[t]->oTet() == nullptr);

                pTris[t]->setOuterTet(pTets[tetouter.get()]);
                // Add this triangle to outer tet's list of neighbours
                for (uint i=0; i <= 4; ++i)
                {
                    AssertLog(i < 4);

                    // See above in that tets now store tets from different comps
                    Tet * tet_out = pTets[tetouter.get()];

                    if (tet_out->nextTet(i) != nullptr) continue;

                    if (tet_out->nextTri(i) != nullptr) continue;
                    tet_out->setNextTri(i, pTris[t]);
                    break;
                }
            }
        }
    }


    // Ok, now to set up the (big) reaction structures

    /// cumulative number of species from each compartment and patch
    pSpecs_tot = 0;
    /// cumulative number of reacs from each compartment and patch
    pReacs_tot = 0;

    uint Comps_N = statedef().countComps();
    uint Patches_N = statedef().countPatches();

    /*
    uint Tets_N = 0.0;
    uint Tris_N = 0.0;
    CompPVecCI comp_end = pComps.end();
    PatchPVecCI patch_end = pPatches.end();
    for(CompPVecCI comp = pComps.begin(); comp != comp_end; ++comp) Tets_N+= (*comp)->countTets;
    for (PatchPVec patch = pPatches.begin(); patch!=patch_end; ++patch) Tris_N+=(*patch)->countTris;
    */


    CompPVecCI comp_end = pComps.end();
    for(CompPVecCI comp = pComps.begin(); comp != comp_end; ++comp)
    {
        uint comp_specs = (*comp)->def()->countSpecs();
        pSpecs_tot+=(comp_specs * (*comp)->countTets());

        uint comp_reacs = (*comp)->def()->countReacs();
        uint comp_diffs = (*comp)->def()->countDiffs();
        // This is not enough indices for diffusion if we should
        // require to change local dcsts in the future
        pReacs_tot+=((comp_reacs+comp_diffs) * (*comp)->countTets());
    }

    PatchPVecCI patch_end = pPatches.end();
    for (PatchPVecCI patch = pPatches.begin(); patch!=patch_end; ++patch)
    {
        uint patch_specs = (*patch)->def()->countSpecs();
        pSpecs_tot+= (patch_specs * (*patch)->countTris());

        uint patch_sreacs = (*patch)->def()->countSReacs();
        uint patch_vdepsreacs = (*patch)->def()->countVDepSReacs();
        uint patch_sdiffs = (*patch)->def()->countSurfDiffs();

        pReacs_tot+= ((patch_sreacs+patch_vdepsreacs+patch_sdiffs) * (*patch)->countTris());
    }

    pSpec_matrixsub =  std::vector< std::vector<steps::tetode::structA> >(pSpecs_tot, std::vector<steps::tetode::structA>());

    //pCcst = new double[pReacs_tot];

    /// set row marker to beginning of matrix for first compartment  (previous rowp)
    uint reac_gidx = 0;
    /// set column marker to beginning of matrix for first compartment (previous colp)
    uint spec_gidx = 0;

    for (uint i=0; i< Comps_N; ++i)
    {
        Comp * comp = pComps[i];
        solver::CompDefP cdef = comp->def();

        uint compReacs_N = cdef->countReacs();
        uint compSpecs_N = cdef->countSpecs();

        uint compTets_N = comp->countTets();

        for (uint t =0; t < compTets_N; ++t)
        {
            for(uint j=0; j< compReacs_N; ++j)
            {
                /// set scaled reaction constant
                double reac_kcst = cdef->kcst(j);
                double tet_vol = comp->getTet(t)->vol();
                uint reac_order = cdef->reacdef(j)->order();
                //pCcst[reac_gidx+j] =_ccst(reac_kcst, comp_vol, reac_order);
                double ccst = _ccst(reac_kcst, tet_vol, reac_order);
                for(uint k=0; k< compSpecs_N; ++k)
                {
                    uint * lhs = cdef->reac_lhs_bgn(j);
                    int upd = cdef->reac_upd_bgn(j)[k];

                    if (upd != 0)
                    {
                        //structB btmp = {std::vector<uint>(), std::vector<uint>()};
                        structB btmp = {std::vector<structC>()};

                        for (uint l=0; l < compSpecs_N; ++l)
                        {
                            uint lhs_spec = lhs[l];
                            if (lhs_spec != 0)
                            {
                                structC ctmp = {lhs_spec, spec_gidx +l};
                                //btmp.order.push_back(lhs_spec);
                                //btmp.spec_idx.push_back(spec_gidx +l);
                                btmp.info.push_back(ctmp);
                            }
                        }
                        structA atmp = {ccst,reac_gidx+j, upd, std::vector<steps::tetode::structB>()};
                        atmp.players.push_back(btmp);

                        pSpec_matrixsub[spec_gidx+k].push_back(atmp);
                    }
                }
            }

            /// step up markers for next compartment
            reac_gidx += compReacs_N;
            spec_gidx += compSpecs_N;
        }

        // No diffusion between compartments (yet) so we can set up diffusion here.
        uint compDiffs_N = cdef->countDiffs();
        for (uint t =0; t < compTets_N; ++t)
        {
            for (uint j=0; j<compDiffs_N; ++j)
            {
                for(uint k=0; k< compSpecs_N; ++k)
                {
                    if (cdef->diff_dep(j, k))
                    {
                        Tet * tet_base = comp->getTet(t);
                        AssertLog(tet_base != nullptr);

                        // The tricky part is to get the correct locations in the
                        // (imaginary) matrix
                        for (uint l = 0; l < 4; ++l)
                        {
                            double dcst = cdef->dcst(j);
                            Tet * tet_neighb = tet_base->nextTet(l);
                            if (tet_neighb== nullptr) continue;

                            // If we are here we found a connection- set up the diffusion 'reaction'
                            double dist = tet_base->dist(l);
                            double dccst = (tet_base->area(l) * dcst) / (tet_base->vol() * dist);

                            // Find the locations in the imaginary spec 'matrix'
                            uint spec_base_idx = 0;

                            for (uint m=0; m<i; ++m)
                            {
                                spec_base_idx += (statedef().compdef(m)->countSpecs())*(pComps[m]->countTets());
                            }

                            // Need to convert neighbour to local index:
                            auto tet_neighb_gidx = tet_neighb->idx();
                            auto tet_neighb_lidx = comp->getTet_GtoL(tet_neighb_gidx);

                            // At the moment spec_base_idx points to the start position in the 'matrix' for this comp (start of the tets)
                            auto spec_neighb_idx = spec_base_idx + (compSpecs_N * tet_neighb_lidx.get()) + k;

                            // Update the base index to point to the correct species (the one that is diffusing out)
                            spec_base_idx+=(compSpecs_N*t)+k;

                            // Create the reaction of diffusion out of tet, affects two species
                            structB btmp_out = {std::vector<structC>()};
                            structC ctmp_out = {1, spec_base_idx};
                            btmp_out.info.push_back(ctmp_out);

                            // The diffusion index is not strictly correct because of all the
                            // different directions. In the future, for local dccsts to be changed,
                            // this needs to be changed
                            structA atmp_out = {dccst,reac_gidx + j,-1, std::vector<steps::tetode::structB>()};
                            atmp_out.players.push_back(btmp_out);

                            structB btmp_in = {std::vector<structC>()};
                            structC ctmp_in = {1, spec_base_idx};
                            btmp_in.info.push_back(ctmp_in);

                            structA atmp_in = {dccst,reac_gidx + j,1, std::vector<steps::tetode::structB>()};
                            atmp_in.players.push_back(btmp_in);

                            pSpec_matrixsub[spec_base_idx].push_back(atmp_out);
                            pSpec_matrixsub[spec_neighb_idx].push_back(atmp_in);
                        }
                    }
                    // Can only depend on one species
                    continue;
                }
            }
            // This is not sufficient if we need to change local dccsts in the future
            // we'll need a new index for each of the 4 directions
            reac_gidx += compDiffs_N;
        }
    } // end of loop over compartments

    // This is going to be used later as a base for patch spec idx for setting up surface diffusion
    uint spec_base_idx_patch = spec_gidx;

    for(uint i=0; i< Patches_N; ++i)
    {
        Patch * patch = pPatches[i];
        solver::PatchDefP pdef = patch->def();

        uint patchReacs_N = pdef->countSReacs();
        uint patchSpecs_N_S = pdef->countSpecs();
        uint patchSpecs_N_I = pdef->countSpecs_I();
        uint patchSpecs_N_O = pdef->countSpecs_O();

        uint patchTris_N = patch->countTris();

        uint patchVDepSReacs_N = pdef->countVDepSReacs();

        for (uint t=0; t < patchTris_N; ++t)
        {
            Tri * tri = patch->getTri(t);

            for (uint j=0; j< patchReacs_N; ++j)
            {
                double ccst=0.0;
                if (!pdef->sreacdef(j)->surf_surf())
                {
                    double reac_kcst = pdef->kcst(j);
                    double vol=0.0;
                    if (pdef->sreacdef(j)->inside())
                    {
                        AssertLog(pdef->icompdef() != nullptr);
                        Tet * itet = tri->iTet();
                        AssertLog(itet!= nullptr);
                        vol = itet->vol();
                    }
                    else
                    {
                        AssertLog(pdef->ocompdef() != nullptr);
                        Tet * otet = tri->oTet();
                        AssertLog(otet != nullptr);
                        vol = otet->vol();
                    }
                    uint sreac_order = pdef->sreacdef(j)->order();
                    ccst = _ccst(reac_kcst, vol, sreac_order);
                }
                else
                {
                    // 2D reaction
                    double area = tri->area();
                    double reac_kcst = pdef->sreacdef(j)->kcst();
                    uint sreac_order = pdef->sreacdef(j)->order();
                    ccst = _ccst2D(reac_kcst, area, sreac_order);
                }

                // This structure will hold all species players for the reaction, which can be in any of 3 locations
                // - the surface, inner comp or outer comp
                structB btmp = {std::vector<structC>()};

                // First we need to collect all LHS information,
                // then do another loop to add the reaction for any
                // species whose upd value is non-zero, but all species
                // can appear in 3 locations - the patch, the inner comp
                // and the outer comp
                uint * slhs = pdef->sreac_lhs_S_bgn(j);
                for (uint l=0; l < patchSpecs_N_S; ++l)
                {
                    uint slhs_spec = slhs[l];
                    if (slhs_spec != 0)
                    {
                        // spec_gidx is up to date for this triangle:
                        structC ctmp = {slhs_spec, spec_gidx +l};
                        btmp.info.push_back(ctmp);
                    }
                }

                auto * icompdef = pdef->icompdef();
                if (icompdef != nullptr)
                {
                    // Sanity check
                    AssertLog(icompdef == tri->iTet()->compdef());
                    uint icompidx = icompdef->gidx();
                    Comp * localicomp = pComps[icompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_itetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< icompidx; ++l)
                    {
                        mtx_itetidx += (statedef().compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    auto tet_gidx = tri->iTet()->idx();
                    auto tet_lidx = localicomp->getTet_GtoL(tet_gidx);
                    mtx_itetidx += tet_lidx.get() * icompdef->countSpecs();

                    uint * ilhs = pdef->sreac_lhs_I_bgn(j);
                    for (uint l=0; l<patchSpecs_N_I; ++l)
                    {
                        uint ilhs_spec = ilhs[l];
                        if (ilhs_spec != 0)
                        {
                            structC ctmp = {ilhs_spec, mtx_itetidx+l};
                            btmp.info.push_back(ctmp);
                        }
                    }
                }

                auto * ocompdef = pdef->ocompdef();
                if (ocompdef != nullptr)
                {
                    // Sanity check
                    AssertLog(ocompdef == tri->oTet()->compdef());
                    uint ocompidx = ocompdef->gidx();
                    Comp * localocomp = pComps[ocompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_otetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< ocompidx; ++l)
                    {
                        mtx_otetidx += (statedef().compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    auto tet_gidx = tri->oTet()->idx();
                    auto tet_lidx = localocomp->getTet_GtoL(tet_gidx);
                    mtx_otetidx += tet_lidx.get() * ocompdef->countSpecs();

                    uint * olhs = pdef->sreac_lhs_O_bgn(j);
                    for (uint l=0; l<patchSpecs_N_O; ++l)
                    {
                        uint olhs_spec = olhs[l];
                        if (olhs_spec != 0)
                        {
                            structC ctmp = {olhs_spec, mtx_otetidx+l};
                            btmp.info.push_back(ctmp);
                        }
                    }
                }

                // I can't see any alternative but to do the loops twice- once to fill the 'players' (lhs's),
                // then go round again and add the reaction to every species involved (update not equal to 1)
                for (uint k=0; k< patchSpecs_N_S; ++k)
                {
                    int supd = pdef->sreac_upd_S_bgn(j)[k];
                    if (supd != 0)
                    {
                        structA atmp = {ccst, reac_gidx+j,supd, std::vector<steps::tetode::structB>()};
                        atmp.players.push_back(btmp);
                        // PROBLEM here that I am going to add the same structure to a lot of vectors- need to figure out if I should copy, or
                        // something else fancy like storing it once and using pointers

                        // WHAT I COULD DO is have a big array of structBs somewhere (these are a bit like reactions)
                        // and add to the array as I go (keeping track of indices- actually is that necessary?)
                        // then the structAs simply store the pointer. This could also be useful for diffusion where the two 'reactions' depend on the same 'players'
                        // PRoblem is that we don't know how big it'll be- so use vectors instead? Then how to store pointer- store vector iterator??
                        // Could also do similar for structAs, then this pSpec_matrixsub only stores pointer to vector

                        // Each time a copy of the vector is made as a new object, which is fine. It will mean perhaps a
                        // 2 or 3 fold increase in memory to using pointers, but there should be some gain to efficiency.
                        pSpec_matrixsub[spec_gidx+k].push_back(atmp);
                    }
                }

                if (icompdef != nullptr)
                {
                    uint icompidx = icompdef->gidx();
                    Comp * localicomp = pComps[icompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_itetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< icompidx; ++l)
                    {
                        mtx_itetidx += (statedef().compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    auto tet_gidx = tri->iTet()->idx();
                    auto tet_lidx = localicomp->getTet_GtoL(tet_gidx);
                    mtx_itetidx += tet_lidx.get() * icompdef->countSpecs();

                    for (uint k=0; k<patchSpecs_N_I; ++k)
                    {
                        int upd = pdef->sreac_upd_I_bgn(j)[k];
                        if (upd!=0)
                        {
                            structA atmp = {ccst, reac_gidx+j,upd,  std::vector<steps::tetode::structB>()};
                            atmp.players.push_back(btmp);
                            pSpec_matrixsub[mtx_itetidx+k].push_back(atmp);
                        }
                    }
                }

                if (ocompdef != nullptr)
                {
                    uint ocompidx = ocompdef->gidx();
                    Comp * localocomp = pComps[ocompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_otetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< ocompidx; ++l)
                    {
                        mtx_otetidx += (statedef().compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    auto tet_gidx = tri->oTet()->idx();
                    auto tet_lidx = localocomp->getTet_GtoL(tet_gidx);
                    mtx_otetidx += tet_lidx.get() * ocompdef->countSpecs();

                    for (uint k=0; k<patchSpecs_N_O; ++k)
                    {
                        int upd = pdef->sreac_upd_O_bgn(j)[k];
                        if (upd!=0)
                        {
                            structA atmp = {ccst, reac_gidx+j,upd,  std::vector<steps::tetode::structB>()};
                            atmp.players.push_back(btmp);
                            pSpec_matrixsub[mtx_otetidx+k].push_back(atmp);
                        }
                    }
                }
            } // end of loop over patch surface reactions

            reac_gidx += patchReacs_N;



            // Just initialise all kcsts as 0 initially

            for (uint j=0; j< patchVDepSReacs_N; ++j)
            {
                double ccst=0.0;

                // This structure will hold all species players for the reaction, which can be in any of 3 locations
                // - the surface, inner comp or outer comp
                structB btmp = {std::vector<structC>()};

                // First we need to collect all LHS information,
                // then do another loop to add the reaction for any
                // species whose upd value is non-zero, but all species
                // can appear in 3 locations - the patch, the inner comp
                // and the outer comp
                uint * slhs = pdef->vdepsreac_lhs_S_bgn(j);
                for (uint l=0; l < patchSpecs_N_S; ++l)
                {
                    uint slhs_spec = slhs[l];
                    if (slhs_spec != 0)
                    {
                        // spec_gidx is up to date for this triangle:
                        structC ctmp = {slhs_spec, spec_gidx +l};
                        btmp.info.push_back(ctmp);
                    }
                }

                auto * icompdef = pdef->icompdef();
                if (icompdef != nullptr)
                {
                    // Sanity check
                    AssertLog(icompdef == tri->iTet()->compdef());
                    uint icompidx = icompdef->gidx();
                    Comp * localicomp = pComps[icompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_itetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< icompidx; ++l)
                    {
                        mtx_itetidx += (statedef().compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    auto tet_gidx = tri->iTet()->idx();
                    auto tet_lidx = localicomp->getTet_GtoL(tet_gidx);
                    mtx_itetidx += tet_lidx.get() * icompdef->countSpecs();

                    uint * ilhs = pdef->vdepsreac_lhs_I_bgn(j);
                    for (uint l=0; l<patchSpecs_N_I; ++l)
                    {
                        uint ilhs_spec = ilhs[l];
                        if (ilhs_spec != 0)
                        {
                            structC ctmp = {ilhs_spec, mtx_itetidx+l};
                            btmp.info.push_back(ctmp);
                        }
                    }
                }

                auto * ocompdef = pdef->ocompdef();
                if (ocompdef != nullptr)
                {
                    // Sanity check
                    AssertLog(ocompdef == tri->oTet()->compdef());
                    uint ocompidx = ocompdef->gidx();
                    Comp * localocomp = pComps[ocompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_otetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< ocompidx; ++l)
                    {
                        mtx_otetidx += (statedef().compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    auto tet_gidx = tri->oTet()->idx();
                    auto tet_lidx = localocomp->getTet_GtoL(tet_gidx);
                    mtx_otetidx += tet_lidx.get() * ocompdef->countSpecs();

                    uint * olhs = pdef->vdepsreac_lhs_O_bgn(j);
                    for (uint l=0; l<patchSpecs_N_O; ++l)
                    {
                        uint olhs_spec = olhs[l];
                        if (olhs_spec != 0)
                        {
                            structC ctmp = {olhs_spec, mtx_otetidx+l};
                            btmp.info.push_back(ctmp);
                        }
                    }
                }

                // I can't see any alternative but to do the loops twice- once to fill the 'players' (lhs's),
                // then go round again and add the reaction to every species involved (update not equal to 1)
                for (uint k=0; k< patchSpecs_N_S; ++k)
                {
                    int supd = pdef->vdepsreac_upd_S_bgn(j)[k];
                    if (supd != 0)
                    {
                        structA atmp = {ccst, reac_gidx+j,supd, std::vector<steps::tetode::structB>()};
                        atmp.players.push_back(btmp);
                        // PROBLEM here that I am going to add the same structure to a lot of vectors- need to figure out if I should copy, or
                        // something else fancy like storing it once and using pointers

                        // WHAT I COULD DO is have a big array of structBs somewhere (these are a bit like reactions)
                        // and add to the array as I go (keeping track of indices- actually is that necessary?)
                        // then the structAs simply store the pointer. This could also be useful for diffusion where the two 'reactions' depend on the same 'players'
                        // PRoblem is that we don't know how big it'll be- so use vectors instead? Then how to store pointer- store vector iterator??
                        // Could also do similar for structAs, then this pSpec_matrixsub only stores pointer to vector

                        // Each time a copy of the vector is made as a new object, which is fine. It will mean perhaps a
                        // 2 or 3 fold increase in memory to using pointers, but there should be some gain to efficiency.
                        pSpec_matrixsub[spec_gidx+k].push_back(atmp);
                    }
                }

                if (icompdef != nullptr)
                {
                    uint icompidx = icompdef->gidx();
                    Comp * localicomp = pComps[icompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_itetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< icompidx; ++l)
                    {
                        mtx_itetidx += (statedef().compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    auto tet_gidx = tri->iTet()->idx();
                    auto tet_lidx = localicomp->getTet_GtoL(tet_gidx);
                    mtx_itetidx += tet_lidx.get() * icompdef->countSpecs();

                    for (uint k=0; k<patchSpecs_N_I; ++k)
                    {
                        int upd = pdef->vdepsreac_upd_I_bgn(j)[k];
                        if (upd!=0)
                        {
                            structA atmp = {ccst, reac_gidx+j,upd,  std::vector<steps::tetode::structB>()};
                            atmp.players.push_back(btmp);
                            pSpec_matrixsub[mtx_itetidx+k].push_back(atmp);
                        }
                    }
                }

                if (ocompdef != nullptr)
                {
                    uint ocompidx = ocompdef->gidx();
                    Comp * localocomp = pComps[ocompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_otetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< ocompidx; ++l)
                    {
                        mtx_otetidx += (statedef().compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    auto tet_gidx = tri->oTet()->idx();
                    auto tet_lidx = localocomp->getTet_GtoL(tet_gidx);
                    mtx_otetidx += tet_lidx.get() * ocompdef->countSpecs();

                    for (uint k=0; k<patchSpecs_N_O; ++k)
                    {
                        int upd = pdef->vdepsreac_upd_O_bgn(j)[k];
                        if (upd!=0)
                        {
                            structA atmp = {ccst, reac_gidx+j,upd,  std::vector<steps::tetode::structB>()};
                            atmp.players.push_back(btmp);
                            pSpec_matrixsub[mtx_otetidx+k].push_back(atmp);
                        }
                    }
                }
            } // end of loop over patch vdep surface reactions

            reac_gidx += patchVDepSReacs_N;

            spec_gidx += patchSpecs_N_S;

        } // end of loop over triangles

        // Set up surface diffusion rules
        uint patchSDiffs_N = pdef->countSurfDiffs();
        for (uint t=0; t < patchTris_N; ++t)
        {
            for (uint j=0; j<patchSDiffs_N; ++j)
            {
                for(uint k=0; k< patchSpecs_N_S; ++k)
                {
                    if (pdef->surfdiff_dep(j,k))
                    {
                        Tri * tri_base = patch->getTri(t);
                        AssertLog(tri_base != nullptr);

                        // Find the correct location in the 'matrix'
                        for (uint l = 0; l < 3; ++l)
                        {
                            double dcst = pdef->dcst(j);
                            Tri * tri_neighb = tri_base->nextTri(l);
                            if (tri_neighb == nullptr) continue;

                            double dist = tri_base->dist(l);
                            double dccst = (tri_base->length(l)*dcst)/(tri_base->area()*dist);

                            // Find the right location in the 'matrix'
                            // BUG HERE - needed to start at gidx after comps: uint spec_base_idx = 0; NO
                            uint spec_base_idx  = spec_base_idx_patch;
                            for (uint m=0; m<i; ++m)
                            {
                                spec_base_idx += (statedef().patchdef(m)->countSpecs())*(pPatches[m]->countTris());
                            }

                            // Need to convert neighbour to local index
                            auto tri_neighb_gidx = tri_neighb->idx();
                            auto tri_neighb_lidx = patch->getTri_GtoL(tri_neighb_gidx);

                            // At the moment spec_base_idx points to the start position
                            // in the 'matrix' for this patch (start of the tris)
                            auto spec_neighb_idx = spec_base_idx+(patchSpecs_N_S*tri_neighb_lidx.get())+k;

                            // Update the base index to point to the correct species (the one that is diffusing out)
                            spec_base_idx+=(patchSpecs_N_S*t)+k;

                            // Create the reaction of diffusion out of tri, affects two species
                            structB btmp_out = {std::vector<structC>()};
                            structC ctmp_out = {1, spec_base_idx};
                            btmp_out.info.push_back(ctmp_out);

                            // The diffusion index is not strictly correct because of all the
                            // different directions. In the future, for local dccsts to be changed,
                            // this needs to be changed
                            structA atmp_out = {dccst,reac_gidx + j,-1, std::vector<steps::tetode::structB>()};
                            atmp_out.players.push_back(btmp_out);

                            structB btmp_in = {std::vector<structC>()};
                            structC ctmp_in = {1, spec_base_idx};
                            btmp_in.info.push_back(ctmp_in);

                            structA atmp_in = {dccst,reac_gidx + j,1, std::vector<steps::tetode::structB>()};
                            atmp_in.players.push_back(btmp_in);

                            pSpec_matrixsub[spec_base_idx].push_back(atmp_out);
                            pSpec_matrixsub[spec_neighb_idx].push_back(atmp_in);
                        }
                    }
                    // Can only depend on one species
                    continue;
                }
            }
            // This is not sufficient if we need to change local dccsts in the future
            // we'll need a new index for each of the 4 directions
            reac_gidx += patchSDiffs_N;
        }

    }// end of loop over patches

    // make sure we added what we expected to
    AssertLog(spec_gidx == pSpecs_tot);
    AssertLog(reac_gidx == pReacs_tot);

    ////////// Now to setup the cvode structures ///////////

    pCVodeState = new CVodeState(pSpecs_tot, 10000, 1.0e-3, 1.0e-3);

    if (efflag()) _setupEField();

}

////////////////////////////////////////////////////////////////////////////////

void TetODE::setTemp(double t)
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "\nWARNING: Temperature set in simulation without membrane ";
        os << "potential calculation will be ignored.\n";
        CLOG(INFO, "general_log") << os.str() << std::endl;
    }
    AssertLog(t >= 0.0);
    pTemp = t;
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setupEField()
{
    using steps::math::point3d;
    using namespace steps::solver::efield;

    //// Note to self: for now roughly following flow from original code in sim/controller.py.
    //// code for setting up a mesh was in func_tetmesh constructor and called functions.

    AssertLog(efflag());

    switch (pEFoption) {
    case EF_DEFAULT:
    case EF_DV_BDSYS:
        pEField = make_EField<dVSolverBanded>();
        break;
    default:
        ArgErrLog("Unsupported E-Field solver.");
    }

    // Give temperature a default value of 20c
    pTemp = 293.15;

    uint nmembs = mesh()->_countMembs();

    if  (nmembs != 1)
    {
        std::ostringstream os;
        os << "Membrane potential solver currently supports only one ";
        os << "membrane description object.";
        ArgErrLog(os.str());
    }

    steps::tetmesh::Memb * memb = mesh()->_getMemb(0);
    AssertLog(memb != nullptr);

    // TODO: Decide what checks are needed for the membrane and implement them here

    pEFNTets = memb->countVolTets();
    pEFNTris = memb->countTris();
    pEFNVerts = memb->countVerts();

    pEFTets = new vertex_id_t[neftets() * 4];

    // All the triangles we will count here are the real membrane triangles,
    // virtual triangles will not require a capacitance.
    pEFTris = new vertex_id_t[neftris() * 3];

    pEFVerts = new double[nefverts() * 3];

    const auto nverts = mesh()->countVertices();
    const auto ntris = mesh()->countTris();
    const auto ntets= mesh()->countTets();

    pEFVert_GtoL = new vertex_id_t[nverts];
    for (uint i=0; i < nverts; ++i) {
        pEFVert_GtoL[i] = UNKNOWN_VER;
    }
    pEFTri_GtoL = new triangle_id_t[ntris];
    for (uint i=0; i< ntris; ++i) {
        pEFTri_GtoL[i] = UNKNOWN_TRI;
    }
    pEFTet_GtoL = new tetrahedron_id_t[ntets];
    for (uint i=0; i < ntets; ++i) {
        pEFTet_GtoL[i] = UNKNOWN_TET;
    }

    pEFTri_LtoG = new triangle_id_t[neftris()];

    // Copy the data to local structures.

    auto const& membverts = memb->_getAllVertIndices();
    AssertLog(membverts.size() == nefverts());
    for (uint efv = 0; efv < nefverts(); ++efv)
    {
        const auto vertidx = membverts[efv];
        point3d verttemp = mesh()->_getVertex(vertidx);
        uint efv2 = efv*3;

        // CONVERTING TO MICRONS HERE. EFIELD OBJECT WILL NOT PERFORM THIS CONVERSION
        verttemp *= 1.0e6;
        pEFVerts[efv2] = verttemp[0];
        pEFVerts[efv2+1] = verttemp[1];
        pEFVerts[efv2+2] = verttemp[2];

        pEFVert_GtoL[vertidx.get()] = efv;
    }

    const auto& membtets = memb->_getAllVolTetIndices();
    AssertLog(membtets.size() == neftets());
    for (uint eft=0; eft < neftets(); ++eft)
    {
        auto tetidx = membtets[eft];
        const auto tettemp = mesh()->_getTet(tetidx);
        uint eft2 = eft*4;

        // Convert to indices used by EField object
        auto tv0 =  pEFVert_GtoL[tettemp[0].get()];
        auto tv1 = pEFVert_GtoL[tettemp[1].get()];
        auto tv2 = pEFVert_GtoL[tettemp[2].get()];
        auto tv3 = pEFVert_GtoL[tettemp[3].get()];
        if  (tv0 == UNKNOWN_VER || tv1 == UNKNOWN_VER || tv2 == UNKNOWN_VER || tv3 == UNKNOWN_VER)
        {
            std::ostringstream os;
            os << "Failed to create EField structures.";
            ProgErrLog(os.str());
        }

        pEFTets[eft2] = tv0;
        pEFTets[eft2 + 1] = tv1;
        pEFTets[eft2 + 2] = tv2;
        pEFTets[eft2 + 3] = tv3;

        pEFTet_GtoL[tetidx.get()] = eft;
    }

    auto const& membtris = memb->_getAllTriIndices();
    AssertLog(membtris.size() == neftris());

    pEFTris_vec.resize(neftris());

    for (uint eft = 0; eft < neftris(); ++eft)
    {
        auto triidx = membtris[eft];
        const auto tritemp = mesh()->_getTri(triidx);
        uint eft2 = eft*3;

        // Convert to indices used by EField object
        auto tv0 =  pEFVert_GtoL[tritemp[0].get()];
        auto tv1 = pEFVert_GtoL[tritemp[1].get()];
        auto tv2 = pEFVert_GtoL[tritemp[2].get()];
        if  (tv0 == UNKNOWN_VER || tv1 == UNKNOWN_VER || tv2 == UNKNOWN_VER)
        {
            std::ostringstream os;
            os << "Failed to create EField structures.";
            ProgErrLog(os.str());
        }

        pEFTris[eft2] = tv0;
        pEFTris[eft2 + 1] = tv1;
        pEFTris[eft2 + 2] = tv2;

        pEFTri_GtoL[triidx.get()] = eft;
        pEFTri_LtoG[eft] = triidx;

        // This is added now for quicker iteration during run()
        // Extremely important for larger meshes, orders of magnitude times faster
        pEFTris_vec[eft] = pTris[triidx.get()];
    }

    using namespace steps::solver::efield;

    pEField->initMesh(nefverts(), pEFVerts, neftris(), pEFTris, neftets(), pEFTets, memb->_getOpt_method(), memb->_getOpt_file_name(), memb->_getSearch_percent());
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::setTolerances(double atol, double rtol)
{
    pCVodeState->setTolerances(atol, rtol);
    pTolsset = true;
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_addTet(tetrahedron_id_t tetidx,
                     steps::tetode::Comp *comp, double vol,
                     double a1, double a2, double a3, double a4,
                     double d1, double d2, double d3, double d4,
                     tetrahedron_id_t tet0, tetrahedron_id_t tet1, tetrahedron_id_t tet2, tetrahedron_id_t tet3)
{
    steps::solver::Compdef * compdef  = comp->def();
    auto localtet = new Tet(tetidx, compdef, vol, a1, a2, a3, a4, d1, d2, d3, d4,
                                         tet0, tet1, tet2, tet3);
    AssertLog(tetidx < static_cast<index_t>(pTets.size()));
    AssertLog(pTets[tetidx.get()] == nullptr);
    pTets[tetidx.get()] = localtet;
    comp->addTet(localtet);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_addTri(triangle_id_t triidx,
                     steps::tetode::Patch *patch,
                     double area,
                     double l0,
                     double l1,
                     double l2,
                     double d0,
                     double d1,
                     double d2,
                     tetrahedron_id_t tinner,
                     tetrahedron_id_t touter,
                     triangle_id_t tri0,
                     triangle_id_t tri1,
                     triangle_id_t tri2)
{
    steps::solver::Patchdef * patchdef = patch->def();
    auto tri = new Tri(triidx, patchdef, area, l0, l1, l2, d0, d1, d2,  tinner, touter, tri0, tri1, tri2);
    AssertLog(triidx < static_cast<index_t>(pTris.size()));
    AssertLog(pTris[triidx.get()] == nullptr);
    pTris[triidx.get()] = tri;
    patch->addTri(tri);
}

////////////////////////////////////////////////////////////////////////////////

std::size_t TetODE::_addComp(steps::solver::Compdef * cdef)
{
    auto comp = new Comp(cdef);
    auto compidx = pComps.size();
    pComps.push_back(comp);
    return compidx;
}

////////////////////////////////////////////////////////////////////////////////

std::size_t TetODE::_addPatch(steps::solver::Patchdef * pdef)
{
    /* Comp * icomp = 0;
     Comp * ocomp = 0;
     if (pdef->icompdef()) icomp = pCompMap[pdef->icompdef()];
     if (pdef->ocompdef()) ocomp = pCompMap[pdef->ocompdef()];
     */
    auto patch = new Patch(pdef);
    auto patchidx = pPatches.size();
    pPatches.push_back(patch);
    return patchidx;
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_ccst(double kcst, double vol, uint order)
{
    double vscale = 1.0e3 * vol * AVOGADRO;
    int o1 = static_cast<int>(order) - 1;

    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;
    return kcst * pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_ccst2D(double kcst, double area, uint order)
{
    double vscale = area * AVOGADRO;
    int o1 = static_cast<int>(order) - 1;
    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;
    return kcst * pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////

void TetODE::advance(double adv)
{
    if (adv < 0.0)
        ArgErrLog("Time to advance cannot be negative.");

    double endtime = statedef().time() + adv;
    run(endtime);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::run(double endtime)
{
    if (endtime < statedef().time())
        ArgErrLog("Endtime is before current simulation time.");

    if (endtime == 0.0) return;

    int flag = 0;

    if (not pInitialised)
    {
        if (not pTolsset)
        {
            CLOG(INFO, "general_log") << "Warning: tolerances have not been set and will ";
            CLOG(INFO, "general_log") << "retain default values\n";
        }

        flag = pCVodeState->initialise();
        pInitialised = true;
    }

    // Call CVodeInit to re- initialize the integrator memory and specify the
    // user's right hand side function in y'=f(t,y), the initial time T0, and
    // the initial dependent variable vector y.
    // Re-initialising here allows for injection of molecules, and possibly other
    // additions in the future such as flags (though this would be a little tricky)

    if (pReinit)
    {
        if (efflag())
        {
            // THis flag is set to true initially so this is a good place to set the VDepSReac constants


            TriPVecCI eftri_end = pEFTris_vec.end();
            uint tlidx = 0;
            for (TriPVecCI eft = pEFTris_vec.begin(); eft != eftri_end; ++eft)
            {
                const auto tgidx = pEFTri_LtoG[tlidx];

                double voltage = pEField->getTriV(tlidx);

                Tri * tri = pTris[tgidx.get()];

                steps::solver::Patchdef * pdef = tri->patchdef();

                uint pidx = pdef->gidx();

                uint nvdepsreacs = pdef->countVDepSReacs();

                for (uint vlidx = 0; vlidx < nvdepsreacs; ++vlidx)
                {
                    double kcst = pdef->vdepsreacdef(vlidx)->getVDepK(voltage);
                    // Calculate the reaction constant
                    double ccst=0.0;
                    if (!pdef->vdepsreacdef(vlidx)->surf_surf())
                    {
                        double vol=0.0;
                        if (pdef->vdepsreacdef(vlidx)->inside())
                        {
                            AssertLog(pdef->icompdef() != nullptr);
                            Tet * itet = tri->iTet();
                            AssertLog(itet!=nullptr);
                            vol = itet->vol();
                        }
                        else
                        {
                            AssertLog(pdef->ocompdef() != nullptr);
                            Tet * otet = tri->oTet();
                            AssertLog(otet != nullptr);
                            vol = otet->vol();
                        }
                        uint sreac_order = pdef->vdepsreacdef(vlidx)->order();
                        ccst = _ccst(kcst, vol, sreac_order);
                    }
                    else
                    {
                        // 2D reaction
                        double area = tri->area();
                        uint sreac_order = pdef->vdepsreacdef(vlidx)->order();
                        ccst = _ccst2D(kcst, area, sreac_order);
                    }


                    // First do the easy part- update the surface species
                    uint reac_idx = 0;
                    uint spec_idx = 0;

                    const auto ncomps = pComps.size();
                    for (uint i=0; i< ncomps; ++i)
                    {
                        spec_idx += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
                        reac_idx += (statedef().compdef(i)->countReacs())*(pComps[i]->countTets());
                        reac_idx += (statedef().compdef(i)->countDiffs())*(pComps[i]->countTets());
                    }

                    // Step up to the correct patch:
                    for (uint i=0; i < pidx; ++i)
                    {
                        spec_idx += (statedef().patchdef(i)->countSpecs())*(pPatches[i]->countTris());
                        reac_idx += (statedef().patchdef(i)->countSReacs())*(pPatches[i]->countTris());
                        reac_idx += (statedef().patchdef(i)->countVDepSReacs())*(pPatches[i]->countTris());
                        reac_idx += (statedef().patchdef(i)->countSurfDiffs())*(pPatches[i]->countTris());
                    }

                    uint patchSpecs_N = pdef->countSpecs();
                    uint patchSReacs_N = pdef->countSReacs();

                    uint patchVDepSReacs_N = pdef->countVDepSReacs();

                    // Step up the index to the right triangle
                    Patch * localpatch = pPatches[pidx];
                    auto tri_lpidx = localpatch->getTri_GtoL(tgidx);

                    // Step up indices to the correct triangle
                    spec_idx += (patchSpecs_N * tri_lpidx.get());
                    reac_idx += (patchSReacs_N * tri_lpidx.get());
                    // The following is right because SReacs and VDepSReacs are added within the same loop over Tris:
                    reac_idx += (patchVDepSReacs_N * tri_lpidx.get());

                    // Not necessary because of separate loops:    reac_idx += (patchSDiffs_N*tri_lpidx);

                    // And set the correct reaction index
                    // I think it's right to add one more patch SReac
                    reac_idx += patchSReacs_N;
                    reac_idx+=vlidx;

                    for (uint k=0; k < patchSpecs_N; ++k)
                    {
                        for (auto& r : pSpec_matrixsub[spec_idx+k]) {
                            if (r.r_idx == reac_idx) {
                                r.ccst = ccst;
                            }
                        }
                    }

                    // Now the complicated part, which is to change the constants in the inner
                    // and/or outer tets

                    // DO not do this! spec_idx = 0;

                    if (pdef->vdepsreacdef(vlidx)->reqInside())
                    {
                        Tet * itet = tri->iTet();
                        AssertLog(itet != nullptr);

                        // Fetch the global index of the comp
                        uint icidx = itet->compdef()->gidx();

                        // First step up the species index to the correct comp
                        uint spec_idx_i = 0;
                        for (uint i=0; i< icidx; ++i)
                        {
                            spec_idx_i += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
                        }
                        uint icompSpecs_N = statedef().compdef(icidx)->countSpecs();
                        AssertLog(icompSpecs_N == pdef->countSpecs_I());

                        // Fetch the comps-local index for the tet
                        auto tlidx2 = pComps[icidx]->getTet_GtoL(itet->idx());

                        // Step up the indices to the right tet:
                        spec_idx_i += icompSpecs_N * tlidx2.get();
                        for(uint k=0; k< icompSpecs_N; ++k)
                        {
                            for (auto& r : pSpec_matrixsub[spec_idx_i+k]) {
                                if (r.r_idx == reac_idx) {
                                    r.ccst = ccst;
                                }
                            }
                        }
                    }

                    if (pdef->vdepsreacdef(vlidx)->reqOutside())
                    {
                        Tet * otet = tri->oTet();
                        AssertLog(otet != nullptr);

                        // Fetch the global index of the comp
                        uint ocidx = otet->compdef()->gidx();

                        // First step up the species index to the correct comp
                        uint spec_idx_o = 0;
                        for (uint i=0; i< ocidx; ++i)
                        {
                            spec_idx_o += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
                        }
                        uint ocompSpecs_N = statedef().compdef(ocidx)->countSpecs();
                        AssertLog(ocompSpecs_N == pdef->countSpecs_O());

                        // Fetch the comps-local index for the tet
                        auto tlidx2 = pComps[ocidx]->getTet_GtoL(otet->idx());

                        // Step up the indices to the right tet:
                        spec_idx_o += ocompSpecs_N * tlidx2.get();
                        for(uint k=0; k< ocompSpecs_N; ++k)
                        {
                            for (auto& r : pSpec_matrixsub[spec_idx_o+k]) {
                                if (r.r_idx == reac_idx) {
                                  r.ccst = ccst;
                                }
                            }
                        }
                    }
                }
                tlidx += 1;
            }
        }

        flag = pCVodeState->reinit(statedef().time());

        pReinit = false;
    }

    flag = pCVodeState->run(endtime);

    if (flag != CV_SUCCESS)
    {
          std::ostringstream os;
          os << "\nCVODE iteration failed\n\n";
          SysErrLog(os.str());
    }

    if (efflag())
    {
        double dt = endtime - statedef().time();

        TriPVecCI eftri_end = pEFTris_vec.end();
        uint tlidx = 0;
        for (TriPVecCI eft = pEFTris_vec.begin(); eft != eftri_end; ++eft)
        {
            double v = pEField->getTriV(tlidx);
            //double cur = (*eft)->getOhmicI(v, dt, this);
            double ohmcur = (*eft)->getOhmicI(v, this);

            // The following method will also move the ions
            double ghkcur = (*eft)->getGHKI(v, dt, this);

            pEField->setTriI(tlidx, ohmcur+ghkcur);
            tlidx++;

        }

        pEField->advance(dt); //Now got to figure out how to update the voltage-dependent reactions, must have to be
        // at the top of this function somewhere

        // TODO: Replace this with something that only resets voltage-dependent things
        pReinit = true;
    }

    statedef().setTime(endtime);
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getCompVol(uint cidx) const
{
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    Comp * comp = pComps[cidx];
    AssertLog(comp != nullptr);
    return comp->vol();
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getCompAmount(uint cidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getCompCount(cidx, sidx);
    return (count / AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setCompAmount(uint cidx, uint sidx, double a)
{
    // convert amount in mols to number of molecules
    double a2 = a * AVOGADRO;
    // the following method does all the necessary argument checking
    _setCompCount(cidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getCompConc(uint cidx, uint sidx) const
{
    // the following methods do all the necessary argument checking
    double count = _getCompCount(cidx, sidx);
    double vol = _getCompVol(cidx);

    return count/ (1.0e3 * vol * AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setCompConc(uint cidx, uint sidx, double c)
{
    // the following method does cidx argument checking
    double vol = _getCompVol(cidx);

    double count = c * (1.0e3 * vol * AVOGADRO);

    // the following method does more argument checking
    _setCompCount(cidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool TetODE::_getCompClamped(uint /*cidx*/, uint /*sidx*/) const
{
    std::ostringstream os;
    os << "getCompClamped not implemented for steps::solver::TetODE solver";
    NotImplErrLog(os.str());

    return false;
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setCompClamped(uint /*cidx*/, uint /*sidx*/, bool /*b*/)
{
    std::ostringstream os;
    os << "setCompClamped not implemented for steps::solver::TetODE solver";
    NotImplErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getCompReacK(uint /*cidx*/, uint /*ridx*/) const
{
    std::ostringstream os;
    os << "getCompReacK not implemented for steps::solver::TetODE solver";
    NotImplErrLog(os.str());

    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setCompReacK(uint cidx, uint ridx, double kf)
{
    AssertLog(cidx < statedef().countComps());
    AssertLog(statedef().countComps() == pComps.size());
    Comp * comp = pComps[cidx];
    AssertLog(comp != nullptr);

    for (auto const& tet : comp->tets()) {
        // The following method will check the ridx argument
        _setTetReacK(tet->idx(), ridx, kf);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool TetODE::_getCompReacActive(uint /*cidx*/, uint /*ridx*/) const
{
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setCompReacActive(uint /*cidx*/, uint /*ridx*/, bool /*a*/)
{
    std::ostringstream os;
    os << "setCompReacActive not implemented for steps::solver::TetODE solver";
    NotImplErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getPatchArea(uint pidx) const
{
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    Patch * patch = pPatches[pidx];
    AssertLog(patch != nullptr);
    return patch->area();
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getPatchAmount(uint pidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getPatchCount(pidx, sidx);
    return (count / AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setPatchAmount(uint pidx, uint sidx, double a)
{
    AssertLog(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * AVOGADRO;
    // the following method does all the necessary argument checking
    _setPatchCount(pidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

bool TetODE::_getPatchClamped(uint /*pidx*/, uint /*sidx*/) const
{
    std::ostringstream os;
    os << "getPatchClamped not implemented for steps::solver::TetODE solver";
    NotImplErrLog(os.str());

    return false;
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setPatchClamped(uint /*pidx*/, uint /*sidx*/, bool /*buf*/)
{
    std::ostringstream os;
    os << "setPatchClamped not implemented for steps::solver::TetODE solver";
    NotImplErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getPatchSReacK(uint /*pidx*/, uint /*ridx*/) const
{
    std::ostringstream os;
    os << "getPatchSReacK not implemented for steps::solver::TetODE solver";
    NotImplErrLog(os.str());
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setPatchSReacK(uint pidx, uint ridx, double kf)
{
    AssertLog(pidx < statedef().countPatches());
    AssertLog(statedef().countPatches() == pPatches.size());
    Patch * patch = pPatches[pidx];
    AssertLog(patch != nullptr);

    for (auto const& tri: patch->tris()) {
        // The following method will check the ridx argument
        _setTriSReacK(tri->idx(), ridx, kf);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool TetODE::_getPatchSReacActive(uint /*pidx*/, uint /*ridx*/) const
{
    std::ostringstream os;
    os << "getPatchSReacActive not implemented for steps::solver::TetODE solver";
    NotImplErrLog(os.str());
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setPatchSReacActive(uint /*pidx*/, uint /*ridx*/, bool /*a*/)
{
    std::ostringstream os;
    os << "setPatchSReacActive not implemented for steps::solver::TetODE solver";
    NotImplErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getCompCount(uint cidx, uint sidx) const
{
    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    Compdef * comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    uint slidx = comp->specG2L(sidx);
    if (slidx == LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    uint idx = 0;
    /// step up marker to correct comp
    for (uint i=0; i< cidx; ++i)
    {
        idx += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
    }

    uint comp_nspecs = comp->countSpecs();
    uint ntets = pComps[cidx]->countTets();

    double count = 0.0;

    AssertLog((idx+(((ntets-1)*comp_nspecs)+slidx)) < pSpecs_tot);
    for (uint i=0; i< ntets; ++i)
    {
        count += Ith(pCVodeState->y_cvode, idx+((i*comp_nspecs)+slidx));
    }

    return count;
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setCompCount(uint cidx, uint sidx, double n)
{

    AssertLog(cidx < statedef().countComps());
    AssertLog(sidx < statedef().countSpecs());
    Compdef * comp = statedef().compdef(cidx);
    AssertLog(comp != nullptr);
    uint slidx = comp->specG2L(sidx);
    if (slidx == LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        ArgErrLog(os.str());
    }

    uint idx = 0;
    /// step up marker to correct comp
    for (uint i=0; i< cidx; ++i)
    {
        idx += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
    }

    Comp* localcomp = pComps[cidx];
    uint ntets = localcomp->countTets();
    uint comp_nspecs = comp->countSpecs();

    double comp_vol = localcomp->vol();

    AssertLog((idx+(((ntets-1)*comp_nspecs)+slidx)) < pSpecs_tot);

    for (uint i=0; i< ntets; ++i)
    {
        double tetvol = localcomp->getTet(i)->vol();
        Ith(pCVodeState->y_cvode, idx+((i*comp_nspecs)+slidx)) = n*(tetvol/comp_vol);
    }
    // Reinitialise CVode structures
    pReinit = true;
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getPatchCount(uint pidx, uint sidx) const
{
    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countSpecs());
    Patchdef * patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);

    uint slidx = patch->specG2L(sidx);
    if (slidx == LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    uint idx = 0;
    // Step up to correct species index. Comps comes first
    for (uint i=0; i< pComps.size(); ++i)
    {
        idx += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
    }
    // quick sanity check
    AssertLog(idx < pSpecs_tot);

    // Now step up to the correct species index
    for (uint i=0; i< pidx; ++i)
    {
        idx+= (statedef().patchdef(i)->countSpecs())*(pPatches[i]->countTris());
    }

    Patch * localpatch = pPatches[pidx];
    uint patch_nspecs = patch->countSpecs();
    uint ntris = localpatch->countTris();

    AssertLog(idx+((ntris-1)*patch_nspecs)+slidx < pSpecs_tot);

    double count = 0.0;
    for (uint i = 0; i < ntris; ++i)
    {
        count += Ith(pCVodeState->y_cvode, idx+((i*patch_nspecs)+slidx));
    }

    return count;
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setPatchCount(uint pidx, uint sidx, double n)
{

    AssertLog(pidx < statedef().countPatches());
    AssertLog(sidx < statedef().countSpecs());
    Patchdef * patch = statedef().patchdef(pidx);
    AssertLog(patch != nullptr);

    uint slidx = patch->specG2L(sidx);
    if (slidx == LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        ArgErrLog(os.str());
    }

    uint idx = 0;
    // Step up to correct species index. Comps comes first
    for (uint i=0; i< pComps.size(); ++i)
    {
        idx += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
    }
    // quick sanity check
    AssertLog(idx < pSpecs_tot);

    // Now step up to the correct species index
    for (uint i=0; i< pidx; ++i)
    {
        idx+= (statedef().patchdef(i)->countSpecs())*(pPatches[i]->countTris());
    }

    Patch * localpatch = pPatches[pidx];
    uint patch_nspecs = patch->countSpecs();
    uint ntris = localpatch->countTris();

    AssertLog((idx+((ntris-1)*patch_nspecs+slidx)) < pSpecs_tot);

    double patch_area = localpatch->area();

    for (uint i=0; i<ntris; ++i)
    {
        double triarea = localpatch->getTri(i)->area();

        Ith(pCVodeState->y_cvode, idx+((i*patch_nspecs)+slidx)) = n*(triarea/patch_area);
    }

    // Reinitialise CVode structures
    pReinit = true;
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getTetCount(tetrahedron_id_t tidx, uint sidx) const
{
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(tidx < static_cast<index_t>(pTets.size()));

    if (pTets[tidx.get()] == nullptr)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }

    Tet * tet = pTets[tidx.get()];

    Compdef * comp = tet->compdef();
    uint cidx = comp->gidx();

    uint lsidx = comp->specG2L(sidx);
    if (lsidx == LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    Comp * localcomp = pComps[cidx];
    auto tet_lcidx = localcomp->getTet_GtoL(tidx);

    uint idx = 0;
    /// step up marker to correct comp
    for (uint i=0; i< cidx; ++i)
    {
        idx += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
    }

    AssertLog((idx + (comp->countSpecs() * tet_lcidx.get()) + lsidx) < pSpecs_tot);

    return Ith(pCVodeState->y_cvode, idx + comp->countSpecs() * tet_lcidx.get() + lsidx);
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getTetConc(tetrahedron_id_t tidx, uint sidx) const
{
    // Following does arg checking on tidx and sidx
    double count = _getTetCount(tidx, sidx);
    double vol = pTets[tidx.get()]->vol();

    return count/(1.0e3*vol* AVOGADRO);

}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setTetCount(tetrahedron_id_t tidx, uint sidx, double n)
{

    AssertLog(sidx < statedef().countSpecs());
    AssertLog(tidx < static_cast<index_t>(pTets.size()));

    if (pTets[tidx.get()] == nullptr)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }

    Tet * tet = pTets[tidx.get()];

    Compdef * comp = tet->compdef();
    uint cidx = comp->gidx();

    uint lsidx = comp->specG2L(sidx);
    if (lsidx == LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    Comp * localcomp = pComps[cidx];
    auto tet_lcidx = localcomp->getTet_GtoL(tidx);

    uint idx = 0;
    /// step up marker to correct comp
    for (uint i=0; i< cidx; ++i)
    {
        idx += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
    }

    AssertLog((idx + (comp->countSpecs() * tet_lcidx.get()) + lsidx) < pSpecs_tot);

    Ith(pCVodeState->y_cvode, idx + comp->countSpecs() * tet_lcidx.get() + lsidx) = n;

    // Reinitialise CVode structures
    pReinit = true;
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getTetReacK(tetrahedron_id_t /*tidx*/, uint /*ridx*/) const
{
    std::ostringstream os;
    os << "getTetReacK not implemented for steps::solver::TetODE solver";
    NotImplErrLog(os.str());
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setTetReacK(tetrahedron_id_t tidx, uint ridx, double kf)
{
    AssertLog(tidx < static_cast<index_t>(pTets.size()));
    AssertLog(ridx < statedef().countReacs());

    if (pTets[tidx.get()] == nullptr)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }

    Tet * tet = pTets[tidx.get()];
    Compdef * comp = tet->compdef();

    uint lridx = comp->reacG2L(ridx);
    if (lridx == LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        ArgErrLog(os.str());
    }

    // Fetch the global index of the comp
    uint cidx = tet->compdef()->gidx();

    // Now the tricky part
    // First step up the species and reaction indices to the correct comp
    uint reac_idx = 0;
    uint spec_idx = 0;
    for (uint i=0; i< cidx; ++i)
    {
        spec_idx += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
        // Diffusion rules are counted as 'reacs' too
        reac_idx += (statedef().compdef(i)->countReacs())*(pComps[i]->countTets());
        reac_idx += (statedef().compdef(i)->countDiffs())*(pComps[i]->countTets());
    }

    uint compSpecs_N = comp->countSpecs();
    uint compReacs_N = comp->countReacs();


    auto tlidx = pComps[cidx]->getTet_GtoL(tidx);
    // Step up the indices to the right tet:
    spec_idx += compSpecs_N * tlidx.get();
    reac_idx += compReacs_N * tlidx.get();

    // This not necessary because 1st loop through tets adds reactions, then later loop adds diffusion: reac_idx += (compDiffs_N*tlidx);

    // And finally to the right reaction:
    reac_idx+=lridx;

    // We have the 'local reac index'
    // We need to loop over all species in this tet and see if it is involved in this
    // reaction by looking at all reaction indices in all reactions for each species
    // First find the ccst
    double tet_vol = tet->vol();
    uint reac_order = comp->reacdef(lridx)->order();
    double ccst = _ccst(kf, tet_vol, reac_order);

    for(uint k=0; k< compSpecs_N; ++k) {
        for (auto& r : pSpec_matrixsub[spec_idx+k]) {
            if (r.r_idx == reac_idx) {
              r.ccst = ccst;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setTetConc(tetrahedron_id_t tidx, uint sidx, double c)
{
    AssertLog(tidx < static_cast<index_t>(pTets.size()));
    if (pTets[tidx.get()] == nullptr)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        ArgErrLog(os.str());
    }

    double vol = pTets[tidx.get()]->vol();

    double count = c * (1.0e3 * vol * AVOGADRO);

    // Further (and repeated) arg checking done by next method
    _setTetCount(tidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getTetAmount(tetrahedron_id_t tidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getTetCount(tidx, sidx);
    return (count / AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setTetAmount(tetrahedron_id_t tidx, uint sidx, double a)
{
    // convert amount in mols to number of molecules
    double a2 = a * AVOGADRO;
    // the following method does all the necessary argument checking
    _setTetCount(tidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getTriCount(triangle_id_t tidx, uint sidx) const
{
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(tidx < static_cast<index_t>(pTris.size()));

    if (pTris[tidx.get()] == nullptr)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        ArgErrLog(os.str());
    }

    Tri * tri = pTris[tidx.get()];

    Patchdef * patch = tri->patchdef();
    uint pidx = patch->gidx();

    uint lsidx = patch->specG2L(sidx);
    if (lsidx == LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    Patch * localpatch = pPatches[pidx];
    auto tri_lpidx = localpatch->getTri_GtoL(tidx.get());

    uint idx = 0;
    // Step over compartments' tetrahedrons
    for (uint i=0; i< pComps.size(); ++i)
    {
        idx += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
    }
    // And step up further over previous patches' triangles
    for (uint i=0; i < pidx; ++i)
    {
        idx += (statedef().patchdef(i)->countSpecs())*(pPatches[i]->countTris());
    }

    AssertLog((idx + (patch->countSpecs() * tri_lpidx.get()) + lsidx) < pSpecs_tot);

    return Ith(pCVodeState->y_cvode, idx + (patch->countSpecs() * tri_lpidx.get()) + lsidx);

}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setTriCount(triangle_id_t tidx, uint sidx, double n)
{
    AssertLog(sidx < statedef().countSpecs());
    AssertLog(tidx < static_cast<index_t>(pTris.size()));

    if (pTris[tidx.get()] == nullptr)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        ArgErrLog(os.str());
    }

    Tri * tri = pTris[tidx.get()];

    Patchdef * patch = tri->patchdef();
    uint pidx = patch->gidx();

    uint lsidx = patch->specG2L(sidx);
    if (lsidx == LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    Patch * localpatch = pPatches[pidx];
    auto tri_lpidx = localpatch->getTri_GtoL(tidx.get());

    uint idx = 0;

    const auto ncomps = pComps.size();
    // Step over compartments' tetrahedrons
    for (uint i=0; i< ncomps; ++i)
    {
        idx += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
    }
    // And step up further over previous patches' triangles
    for (uint i=0; i < pidx; ++i)
    {
        idx += (statedef().patchdef(i)->countSpecs())*(pPatches[i]->countTris());
    }

    AssertLog((idx + (patch->countSpecs() * tri_lpidx.get()) + lsidx) < pSpecs_tot);

    Ith(pCVodeState->y_cvode, idx + (patch->countSpecs() * tri_lpidx.get()) + lsidx) = n;

    // Reinitialise CVode structures
    pReinit = true;
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getTriAmount(triangle_id_t tidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getTriCount(tidx, sidx);
    return (count / AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setTriAmount(triangle_id_t tidx, uint sidx, double a)
{
    AssertLog(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * AVOGADRO;
    // the following method does all the necessary argument checking
    _setTriCount(tidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getTriSReacK(triangle_id_t /*tidx*/, uint /*ridx*/) const
{
    std::ostringstream os;
    os << "getTriSReacK not implemented for steps::solver::TetODE solver";
    NotImplErrLog(os.str());
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setTriSReacK(triangle_id_t tidx, uint ridx, double kf)
{
    AssertLog(ridx < statedef().countSReacs());
    AssertLog(tidx < static_cast<index_t>(pTris.size()));

    if (pTris[tidx.get()] == nullptr)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        ArgErrLog(os.str());
    }

    Tri * tri = pTris[tidx.get()];

    Patchdef * patch = tri->patchdef();

    // Fetch the global index of the patch
    uint pidx = patch->gidx();

    uint lsridx = patch->sreacG2L(ridx);

    if (lsridx == LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface Reaction undefined in triangle.\n";
        ArgErrLog(os.str());
    }

    // Calculate the reaction constant
    double ccst=0.0;
    if (!patch->sreacdef(lsridx)->surf_surf())
    {
        double vol=0.0;
        if (patch->sreacdef(lsridx)->inside())
        {
            AssertLog(patch->icompdef() != nullptr);
            Tet * itet = tri->iTet();
            AssertLog(itet!=0);
            vol = itet->vol();
        }
        else
        {
            AssertLog(patch->ocompdef() != nullptr);
            Tet * otet = tri->oTet();
            AssertLog(otet != nullptr);
            vol = otet->vol();
        }
        uint sreac_order = patch->sreacdef(lsridx)->order();
        ccst = _ccst(kf, vol, sreac_order);
    }
    else
    {
        // 2D reaction
        double area = tri->area();
        uint sreac_order = patch->sreacdef(lsridx)->order();
        ccst = _ccst2D(kf, area, sreac_order);
    }

    // First do the easy part- update the surface species
    uint reac_idx = 0;
    uint spec_idx = 0;

    const auto ncomps = pComps.size();
    for (uint i=0; i< ncomps; ++i)
    {
        spec_idx += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
        reac_idx += (statedef().compdef(i)->countReacs())*(pComps[i]->countTets());
        reac_idx += (statedef().compdef(i)->countDiffs())*(pComps[i]->countTets());
    }

    // Step up to the correct patch:
    for (uint i=0; i < pidx; ++i)
    {
        spec_idx += (statedef().patchdef(i)->countSpecs())*(pPatches[i]->countTris());
        reac_idx += (statedef().patchdef(i)->countSReacs())*(pPatches[i]->countTris());
        reac_idx += (statedef().patchdef(i)->countVDepSReacs())*(pPatches[i]->countTris());
        reac_idx += (statedef().patchdef(i)->countSurfDiffs())*(pPatches[i]->countTris());
    }

    uint patchSpecs_N = patch->countSpecs();
    uint patchSReacs_N = patch->countSReacs();
    uint patchVDepSReacs_N = patch->countVDepSReacs();

    // Step up the index to the right triangle
    Patch * localpatch = pPatches[pidx];
    auto tri_lpidx = localpatch->getTri_GtoL(tidx.get());

    // Step up indices to the correct triangle
    spec_idx += patchSpecs_N * tri_lpidx.get();
    reac_idx += patchSReacs_N * tri_lpidx.get();
    // The following is right because SReacs and VDepSReacs are added within the same loop over Tris:
    reac_idx += patchVDepSReacs_N * tri_lpidx.get();

    // This is not necessary because separate loop over tris for sreacs, then sdiffs:     reac_idx += (patchSDiffs_N*tri_lpidx);

    // And set the correct reaction index
    reac_idx+=lsridx;

    for (uint k=0; k < patchSpecs_N; ++k) {
        for (auto& r : pSpec_matrixsub[spec_idx+k]) {
            if (r.r_idx == reac_idx) {
              r.ccst = ccst;
            }
        }
    }

    // Now the complicated part, which is to change the constants in the inner
    // and/or outer tets

    // reset the spec idx
    spec_idx = 0;

    if (patch->sreacdef(lsridx)->reqInside())
    {
        Tet * itet = tri->iTet();
        AssertLog(itet != nullptr);

        // Fetch the global index of the comp
        uint icidx = itet->compdef()->gidx();

        // First step up the species index to the correct comp
        uint spec_idx2 = 0;
        for (uint i=0; i< icidx; ++i)
        {
            spec_idx2 += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
        }
        uint icompSpecs_N = statedef().compdef(icidx)->countSpecs();
        AssertLog(icompSpecs_N == patch->countSpecs_I());

        // Fetch the comps-local index for the tet
        auto tlidx = pComps[icidx]->getTet_GtoL(itet->idx());

        // Step up the indices to the right tet:
        spec_idx2 += icompSpecs_N * tlidx.get();
        for(uint k=0; k< icompSpecs_N; ++k)
        {
            for (auto& r : pSpec_matrixsub[spec_idx2+k]) {
                if (r.r_idx == reac_idx) {
                  r.ccst = ccst;
                }
            }
        }
    }

    if (patch->sreacdef(lsridx)->reqOutside())
    {
        Tet * otet = tri->oTet();
        AssertLog(otet != nullptr);

        // Fetch the global index of the comp
        uint ocidx = otet->compdef()->gidx();

        // First step up the species index to the correct comp
        uint spec_idx2 = 0;
        for (uint i=0; i< ocidx; ++i)
        {
            spec_idx2 += (statedef().compdef(i)->countSpecs())*(pComps[i]->countTets());
        }
        uint ocompSpecs_N = statedef().compdef(ocidx)->countSpecs();
        AssertLog(ocompSpecs_N == patch->countSpecs_O());

        // Fetch the comps-local index for the tet
        auto tlidx = pComps[ocidx]->getTet_GtoL(otet->idx());

        // Step up the indices to the right tet:
        spec_idx2 += ocompSpecs_N * tlidx.get();
        for(uint k=0; k< ocompSpecs_N; ++k)
        {
            for (auto& r : pSpec_matrixsub[spec_idx2+k]) {
                if (r.r_idx == reac_idx) {
                  r.ccst = ccst;
                }
            }
        }

    }
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getTriArea(triangle_id_t tidx) const
{
    AssertLog(tidx < static_cast<index_t>(pTris.size()));

    if (pTris[tidx.get()] == nullptr)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.";
        ArgErrLog(os.str());
    }

    return pTris[tidx.get()]->area();
}
////////////////////////////////////////////////////////////////////////////////

double TetODE::_getTetVol(tetrahedron_id_t tidx) const
{
    AssertLog(tidx < static_cast<index_t>(pTets.size()));
    if (pTets[tidx.get()] == nullptr)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
        ArgErrLog(os.str());
    }
    return pTets[tidx.get()]->vol();
}

////////////////////////////////////////////////////////////////////////////////


double TetODE::_getTetV(tetrahedron_id_t tidx) const
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto loctidx = pEFTet_GtoL[tidx.get()];
    if (loctidx == UNKNOWN_TET)
    {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    // EField object should convert value to base s.i. units
    return pEField->getTetV(loctidx);

}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setTetV(tetrahedron_id_t tidx, double v)
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto loctidx = pEFTet_GtoL[tidx.get()];
    if (loctidx == UNKNOWN_TET)
    {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    // EField object should convert to millivolts
    pEField->setTetV(loctidx, v);
}

////////////////////////////////////////////////////////////////////////////////

bool TetODE::_getTetVClamped(tetrahedron_id_t tidx) const
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto loctidx = pEFTet_GtoL[tidx.get()];
    if (loctidx == UNKNOWN_TET)
    {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    return pEField->getTetVClamped(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setTetVClamped(tetrahedron_id_t tidx, bool cl)
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto loctidx = pEFTet_GtoL[tidx.get()];
    if (loctidx == UNKNOWN_TET)
    {
        std::ostringstream os;
        os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
        ArgErrLog(os.str());
    }

    pEField->setTetVClamped(loctidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getTriV(triangle_id_t tidx) const
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto loctidx = pEFTri_GtoL[tidx.get()];
    if (loctidx == UNKNOWN_TRI)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }
    // EField object should convert value to base s.i. units
    return pEField->getTriV(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setTriV(triangle_id_t tidx, double v)
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto loctidx = pEFTri_GtoL[tidx.get()];
    if (loctidx == UNKNOWN_TRI)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    // EField object should convert to millivolts
    pEField->setTriV(loctidx, v);
}

////////////////////////////////////////////////////////////////////////////////

bool TetODE::_getTriVClamped(triangle_id_t tidx) const
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto loctidx = pEFTri_GtoL[tidx.get()];
    if (loctidx == UNKNOWN_TRI)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    return pEField->getTriVClamped(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setTriVClamped(triangle_id_t tidx, bool cl)
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto loctidx = pEFTri_GtoL[tidx.get()];
    if (loctidx == UNKNOWN_TRI)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    pEField->setTriVClamped(loctidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getTriI(triangle_id_t tidx) const
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }

    auto loctidx = pEFTri_GtoL[tidx.get()];
    if (loctidx == UNKNOWN_TRI)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    // EField object should convert to required units
    return pEField->getTriI(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setVertIClamp(vertex_id_t vidx, double cur)
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto locvidx = pEFVert_GtoL[vidx.get()];
    if (locvidx == UNKNOWN_VER)
    {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }

    // EField object should convert to required units
    pEField->setVertIClamp(locvidx, cur);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setTriIClamp(triangle_id_t tidx, double cur)
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto loctidx = pEFTri_GtoL[tidx.get()];
    if (loctidx == UNKNOWN_TRI)
    {
        std::ostringstream os;
        os << "Triangle index " << tidx << " not assigned to a membrane.";
        ArgErrLog(os.str());
    }

    // EField object should convert to required units
    pEField->setTriIClamp(loctidx, cur);
}

////////////////////////////////////////////////////////////////////////////////

double TetODE::_getVertV(vertex_id_t vidx) const
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto locvidx = pEFVert_GtoL[vidx.get()];
    if (locvidx == UNKNOWN_VER)
    {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }
    // EField object should convert value to base s.i. units
    return pEField->getVertV(locvidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setVertV(vertex_id_t vidx, double v)
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto locvidx = pEFVert_GtoL[vidx.get()];
    if (locvidx == UNKNOWN_VER)
    {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }
    // EField object should convert to millivolts
    pEField->setVertV(locvidx, v);
}

////////////////////////////////////////////////////////////////////////////////

bool TetODE::_getVertVClamped(vertex_id_t vidx) const
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto locvidx = pEFVert_GtoL[vidx.get()];
    if (locvidx == UNKNOWN_VER)
    {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }

    return pEField->getVertVClamped(locvidx);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setVertVClamped(vertex_id_t vidx, bool cl)
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    auto locvidx = pEFVert_GtoL[vidx.get()];
    if (locvidx == UNKNOWN_VER)
    {
        std::ostringstream os;
        os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
        ArgErrLog(os.str());
    }
    // EField object should convert to millivolts
    pEField->setVertVClamped(locvidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setMembRes(uint midx, double ro, double vrev)
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    if (ro <= 0.0)
    {
        std::ostringstream os;
        os << "Resistivity must be greater than zero.";
        ArgErrLog(os.str());
    }
    // EField object should convert to required units
    AssertLog(midx == 0);
    pEField->setSurfaceResistivity(midx, ro, vrev);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setMembPotential(uint midx, double v)
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    // EField object should convert to millivolts
    AssertLog(midx == 0);
    pEField->setMembPotential(midx, v);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setMembCapac(uint midx, double cm)
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    if (cm < 0.0)
    {
        std::ostringstream os;
        os << "Capacitance must be greater than or equal to zero.";
        ArgErrLog(os.str());
    }


    // EField object should convert to required units
    AssertLog(midx == 0);
    pEField->setMembCapac(midx, cm);
}

////////////////////////////////////////////////////////////////////////////////

void TetODE::_setMembVolRes(uint midx, double ro)
{
    if (!efflag())
    {
        std::ostringstream os;
        os << "Method not available: EField calculation not included in simulation.";
        ArgErrLog(os.str());
    }
    if (ro < 0.0)
    {
        std::ostringstream os;
        os << "Resistivity must be greater than or equal to zero.";
        ArgErrLog(os.str());
    }
    // EField object should convert to required units
    AssertLog(midx == 0);
    pEField->setMembVolRes(midx, ro);
}

////////////////////////////////////////////////////////////////////////////////
// END

 }  // namespace tetode
 }  // namespace steps
