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


#include <iostream>
#include <sstream>
#include <cmath>
#include <memory>

#include "steps/tetode/tetode.hpp"

#include "steps/tetode/tet.hpp"
#include "steps/tetode/tri.hpp"
#include "steps/tetode/comp.hpp"
#include "steps/tetode/patch.hpp"

#include "steps/common.h"
#include "steps/math/constants.hpp"
#include "steps/math/point.hpp"

#include "steps/solver/compdef.hpp"
#include "steps/solver/patchdef.hpp"
#include "steps/solver/reacdef.hpp"
#include "steps/solver/sreacdef.hpp"
#include "steps/solver/vdepsreacdef.hpp"

#include "steps/geom/tetmesh.hpp"

#include "steps/error.hpp"

#include "third_party/cvode-2.6.0/src/cvode/cvode.h"                 /* prototypes for CVODE fcts., consts. */
#include "third_party/cvode-2.6.0/src/nvec_ser/nvector_serial.h"      /* serial N_Vector types, fcts., macros */
#include "third_party/cvode-2.6.0/src/cvode/cvode_dense.h"          /* prototype for CVDense */
#include "third_party/cvode-2.6.0/src/sundials/sundials_dense.h"     /* definitions DlsMat DENSE_ELEM */
#include "third_party/cvode-2.6.0/src/sundials/sundials_types.h"     /* definition of type realtype */
#include "third_party/cvode-2.6.0/src/sundials/sundials_nvector.h"

#include "steps/solver/efield/efield.hpp"
#include "steps/solver/efield/dVsolver.hpp"

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
std::vector<std::vector<steps::tetode::structA> >  pSpec_matrixsub = std::vector<std::vector<steps::tetode::structA> > ();

////////////////////////////////////////////////////////////////////////////////

namespace stode = steps::tetode;
namespace ssolver = steps::solver;
namespace smath = steps::math;

using steps::math::point3d;

// CVODE stuff
struct stode::CVodeState {
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

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace tetode {

void check_flag(void *flagvalue, const char *funcname, int opt)
{
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL)
    {
        std::ostringstream os;
        os << "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",  funcname;
        throw steps::SysErr(os.str());
    }

    /* Check if flag < 0 */
    else if (opt == 1)
    {
        errflag = (int *) flagvalue;
        if (*errflag < 0)
        {
            std::ostringstream os;
            os << "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag;
            throw steps::SysErr(os.str());
        }
     }
}

}
}

////////////////////////////////////////////////////////////////////////////////

int f_cvode(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    uint i = 0;
    std::vector< std::vector<steps::tetode::structA> >::const_iterator sp_end = pSpec_matrixsub.end();
    //for (uint i = 0; i < tetode->pSpecs_tot; ++i)
    for (std::vector< std::vector<steps::tetode::structA> >::const_iterator sp= pSpec_matrixsub.begin(); sp!=sp_end; ++sp)
    {
        double dydt = 0.0;
        std::vector<steps::tetode::structA>::const_iterator r_end = (*sp).end();
        for (std::vector<steps::tetode::structA>::const_iterator r = (*sp).begin(); r!=r_end; ++r)
        {
            double dydt_r = (*r).upd*(*r).ccst;
            std::vector<steps::tetode::structB>::const_iterator p_end = (*r).players.end();
            for (std::vector<steps::tetode::structB>::const_iterator p = (*r).players.begin(); p!=p_end; ++p)
            {
                std::vector<steps::tetode::structC>::const_iterator q_end = (*p).info.end();
                for (std::vector<steps::tetode::structC>::const_iterator q = (*p).info.begin(); q!=q_end; ++q)
                {
                    double val = Ith(y,(*q).spec_idx);
                    uint order = (*q).order;
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

    return (0);
}

////////////////////////////////////////////////////////////////////////////////

stode::CVodeState::CVodeState(uint N_, uint maxn, realtype atol, realtype rtol) {
    N = N_;
    Nmax_cvode = maxn;

    // Creates serial vectors for y and absolute tolerances
    y_cvode = N_VNew_Serial(N);
    check_flag((void *)y_cvode, "N_VNew_Serial", 0);

    abstol_cvode = N_VNew_Serial(N);
    check_flag((void *)abstol_cvode, "N_VNew_Serial", 0);

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

    check_flag((void *)cvode_mem_cvode, "CVodeCreate", 0);

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

stode::CVodeState::~CVodeState() {
    N_VDestroy_Serial(y_cvode);
    N_VDestroy_Serial(abstol_cvode);

    /* Free integrator memory */
    CVodeFree(&cvode_mem_cvode);
}

////////////////////////////////////////////////////////////////////////////////

void stode::CVodeState::checkpoint(std::fstream &cp_file) {
    cp_file.write((char*)&Nmax_cvode, sizeof(uint));
    cp_file.write((char*)&reltol_cvode, sizeof(realtype));
    cp_file.write((char*)(((N_VectorContent_Serial)(abstol_cvode->content))->data), sizeof(realtype) * N);
    cp_file.write((char*)(((N_VectorContent_Serial)(y_cvode->content))->data), sizeof(realtype) * N);
}

////////////////////////////////////////////////////////////////////////////////

void stode::CVodeState::restore(std::fstream &cp_file) {
    cp_file.read((char*)&Nmax_cvode, sizeof(uint));
    cp_file.read((char*)&reltol_cvode, sizeof(realtype));
    cp_file.read((char*)(((N_VectorContent_Serial)(abstol_cvode->content))->data), sizeof(realtype) * N);
    cp_file.read((char*)(((N_VectorContent_Serial)(y_cvode->content))->data), sizeof(realtype) * N);
}

////////////////////////////////////////////////////////////////////////////////

void stode::CVodeState::setTolerances(double atol, double rtol) {
    // I suppose they shouldn't be negative
    if (atol < 0.0 or rtol < 0.0)
    {
        std::ostringstream os;
        os << "Neither absolute tolerance nor relative tolerance should ";
        os << "be negative.\n";
        throw(steps::ArgErr(os.str()));
    }

    reltol_cvode = rtol;

    for (uint i=0; i<N; ++i)
    {
        Ith(abstol_cvode, i) = atol;
    }
}

////////////////////////////////////////////////////////////////////////////////

void stode::CVodeState::setMaxNumSteps(uint maxn) {
    int flag = CVodeSetMaxNumSteps(cvode_mem_cvode, maxn);
    check_flag(&flag, "CVodeSetMaxNumSteps", 1);

    Nmax_cvode = maxn;
}

////////////////////////////////////////////////////////////////////////////////

int stode::CVodeState::initialise() {
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

int stode::CVodeState::reinit(realtype starttime) {
    int flag = CVodeReInit(cvode_mem_cvode, starttime, y_cvode);
    check_flag(&flag, "CVodeInit", 1);

    return flag;
}

int stode::CVodeState::run(realtype endtime) {
    realtype t;
    return CVode(cvode_mem_cvode, endtime, y_cvode, &t, CV_NORMAL);
}

////////////////////////////////////////////////////////////////////////////////

stode::TetODE::TetODE(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r,
         int calcMembPot)
: API(m, g, r)
, pMesh(0)
, pComps()
, pPatches()
, pTris()
, pTets()
, pSpecs_tot(0)
, pReacs_tot(0)
, pCVodeState(0)
, pInitialised(false)
, pTolsset(false)
, pReinit(true)
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
    _setup();
}

////////////////////////////////////////////////////////////////////////////////

stode::TetODE::~TetODE(void)
{
    CompPVecCI comp_e = pComps.end();
    for (CompPVecCI c = pComps.begin(); c != comp_e; ++c) delete *c;
    PatchPVecCI patch_e = pPatches.end();
    for (PatchPVecCI p = pPatches.begin(); p != patch_e; ++p) delete *p;

    TetPVecCI tet_e = pTets.end();
    for (TetPVecCI t = pTets.begin(); t != tet_e; ++t)
    {
        if ((*t) != 0) delete (*t);
    }

    TriPVecCI tri_e = pTris.end();
    for (TriPVecCI t = pTris.end(); t != tri_e; ++t)
    {
        if ((*t) != 0) delete (*t);
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

std::string stode::TetODE::getSolverName(void) const
{
    return "tetODE";
}

////////////////////////////////////////////////////////////////////////////////

std::string stode::TetODE::getSolverDesc(void) const
{
    return "Reaction-diffusion ODE solver in tetrahedral mesh";
}

////////////////////////////////////////////////////////////////////////////////

std::string stode::TetODE::getSolverAuthors(void) const
{
    return "Iain Hepburn";
}

////////////////////////////////////////////////////////////////////////////////

std::string stode::TetODE::getSolverEmail(void) const
{
    return "steps.dev@gmail.com";
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::checkpoint(std::string const & file_name)
{
    std::fstream cp_file;

    cp_file.open(file_name.c_str(),
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    statedef()->checkpoint(cp_file);

    for (uint c = 0; c < pComps.size(); c++) {
        pComps[c]->checkpoint(cp_file);
    }

    for (uint p = 0; p < pPatches.size(); p++) {
        pPatches[p]->checkpoint(cp_file);
    }

    for (uint tri = 0; tri < pTris.size(); tri++) {
        pTris[tri]->checkpoint(cp_file);
    }

    for (uint tet = 0; tet < pTets.size(); tet++) {
        pTets[tet]->checkpoint(cp_file);
    }

    pCVodeState->checkpoint(cp_file);

    if (efflag()) {
        cp_file.write((char*)&pTemp, sizeof(double));
        cp_file.write((char*)&pEFDT, sizeof(double));
        pEField->checkpoint(cp_file);
    }

    cp_file.close();
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::restore(std::string const & file_name)
{
    std::fstream cp_file;

    cp_file.open(file_name.c_str(),
                 std::fstream::in | std::fstream::binary);

    cp_file.seekg(0);

    statedef()->restore(cp_file);

    for (uint c = 0; c < pComps.size(); c++) {
        pComps[c]->restore(cp_file);
    }

    for (uint p = 0; p < pPatches.size(); p++) {
        pPatches[p]->restore(cp_file);
    }

    for (uint tri = 0; tri < pTris.size(); tri++) {
        pTris[tri]->restore(cp_file);
    }

    for (uint tet = 0; tet < pTets.size(); tet++) {
        pTets[tet]->restore(cp_file);
    }

    pCVodeState->restore(cp_file);

    if (efflag()) {
        cp_file.read((char*)&pTemp, sizeof(double));
        cp_file.read((char*)&pEFDT, sizeof(double));
        pEField->restore(cp_file);
    }

    cp_file.close();

    pTolsset = true;
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::reset(void)
{
    std::ostringstream os;
    os << "reset() not implemented for steps::solver::TetODE solver";
    throw steps::NotImplErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::setMaxNumSteps(uint maxn)
{
    pCVodeState->setMaxNumSteps(maxn);
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::getTime(void) const
{
    return statedef()->time();
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setup(void)
{
    // Perform upcast.
    pMesh = dynamic_cast<steps::tetmesh::Tetmesh *>(geom());
    if (!pMesh)
        throw steps::ArgErr("Geometry description to steps::solver::TetODE solver "
                "constructor is not a valid steps::tetmesh::Tetmesh object.");

    uint ntets = pMesh->countTets();
    uint ntris = pMesh->countTris();
    uint ncomps = pMesh->_countComps();
    uint npatches = pMesh->_countPatches();

    assert(npatches == statedef()->countPatches());
    assert(ncomps == statedef()->countComps());

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

    assert(pPatches.size() == npatches);
    assert(pComps.size() == ncomps);

    pTets.assign(ntets, nullptr);
    pTris.assign(ntris, nullptr);

    for (uint p = 0; p < npatches; ++p)
    {
        // Now add the tris for this patch
        steps::wm::Patch * wmpatch = pMesh->_getPatch(p);

        // sanity check
        assert (statedef()->getPatchIdx(wmpatch) == p );

        // Perform upcast
        steps::tetmesh::TmPatch *tmpatch = dynamic_cast<steps::tetmesh::TmPatch*>(wmpatch);
        if (!tmpatch)
            throw steps::ArgErr("Well-mixed patches not supported in steps::solver::TetODE solver.");

        steps::tetode::Patch *localpatch = pPatches[p];

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

    for (uint c = 0; c < ncomps; ++c)
    {
        // Now add the tets for this comp
        steps::wm::Comp * wmcomp = pMesh->_getComp(c);

        // sanity check
        assert(statedef()->getCompIdx(wmcomp) == c);

        // Perform upcast
        steps::tetmesh::TmComp *tmcomp = dynamic_cast<steps::tetmesh::TmComp*>(wmcomp);
        if (!tmcomp)
            throw steps::ArgErr("Well-mixed compartments not supported in steps::solver::TetODE solver.");
        
        steps::tetode::Comp *localcomp = pComps[c];
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

        // Now correct check, previously didn't allow for tet index == 0
        assert(tetinner >= 0);
        assert(pTets[tetinner] != 0 );


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

                //     Also check tris because in some cases a surface tet
                // may have more than 1 neighbouring tri
                // NOTE: The order here will end up being different to the
                // neighbour order at the Tetmesh level

                stode::Tet * tet_in = pTets[tetinner];
                if (tet_in->nextTet(i) != 0) continue;

                if (tet_in->nextTri(i) != 0) continue;
                tet_in->setNextTri(i, pTris[t]);
                break;
            }
        }

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
                    stode::Tet * tet_out = pTets[tetouter];

                    if (tet_out->nextTet(i) != 0) continue;

                    if (tet_out->nextTri(i) != 0) continue;
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

    uint Comps_N = statedef()->countComps();
    uint Patches_N = statedef()->countPatches();

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
        ssolver::CompDefP cdef = comp->def();

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
                        assert(tet_base != 0);

                        // The tricky part is to get the correct locations in the
                        // (imaginary) matrix
                        for (uint l = 0; l < 4; ++l)
                        {
                            double dcst = cdef->dcst(j);
                            Tet * tet_neighb = tet_base->nextTet(l);
                            if (tet_neighb==0) continue;

                            // If we are here we found a connection- set up the diffusion 'reaction'
                            double dist = tet_base->dist(l);
                            double dccst = (tet_base->area(l) * dcst) / (tet_base->vol() * dist);

                            // Find the locations in the imaginary spec 'matrix'
                            uint spec_base_idx = 0;

                            for (uint m=0; m<i; ++m)
                            {
                                spec_base_idx += (statedef()->compdef(m)->countSpecs())*(pComps[m]->countTets());
                            }

                            // Need to convert neighbour to local index:
                            uint tet_neighb_gidx = tet_neighb->idx();
                            uint tet_neighb_lidx = comp->getTet_GtoL(tet_neighb_gidx);

                            // At the moment spec_base_idx points to the start position in the 'matrix' for this comp (start of the tets)
                            uint spec_neighb_idx = spec_base_idx+(compSpecs_N*tet_neighb_lidx)+k;

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

    for(uint i=0; i< Patches_N; ++i)
    {
        Patch * patch = pPatches[i];
        ssolver::PatchDefP pdef = patch->def();

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
                if (pdef->sreacdef(j)->surf_surf() == false)
                {
                    double reac_kcst = pdef->kcst(j);
                    double vol=0.0;
                    if (pdef->sreacdef(j)->inside() == true)
                    {
                        assert(pdef->icompdef() != 0);
                        Tet * itet = tri->iTet();
                        assert(itet!=0);
                        vol = itet->vol();
                    }
                    else
                    {
                        assert(pdef->ocompdef() != 0);
                        Tet * otet = tri->oTet();
                        assert(otet != 0);
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

                ssolver::Compdef * icompdef = pdef->icompdef();
                if (icompdef != 0)
                {
                    // Sanity check
                    assert(icompdef == tri->iTet()->compdef());
                    uint icompidx = icompdef->gidx();
                    Comp * localicomp = pComps[icompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_itetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< icompidx; ++l)
                    {
                        mtx_itetidx += (statedef()->compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    uint tet_gidx = tri->iTet()->idx();
                    uint tet_lidx = localicomp->getTet_GtoL(tet_gidx);
                    mtx_itetidx+=(tet_lidx*icompdef->countSpecs());

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

                ssolver::Compdef * ocompdef = pdef->ocompdef();
                if (ocompdef != 0)
                {
                    // Sanity check
                    assert(ocompdef == tri->oTet()->compdef());
                    uint ocompidx = ocompdef->gidx();
                    Comp * localocomp = pComps[ocompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_otetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< ocompidx; ++l)
                    {
                        mtx_otetidx += (statedef()->compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    uint tet_gidx = tri->oTet()->idx();
                    uint tet_lidx = localocomp->getTet_GtoL(tet_gidx);
                    mtx_otetidx+=(tet_lidx*ocompdef->countSpecs());

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

                if (icompdef != 0)
                {
                    uint icompidx = icompdef->gidx();
                    Comp * localicomp = pComps[icompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_itetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< icompidx; ++l)
                    {
                        mtx_itetidx += (statedef()->compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    uint tet_gidx = tri->iTet()->idx();
                    uint tet_lidx = localicomp->getTet_GtoL(tet_gidx);
                    mtx_itetidx+=(tet_lidx*icompdef->countSpecs());

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

                if (ocompdef != 0)
                {
                    uint ocompidx = ocompdef->gidx();
                    Comp * localocomp = pComps[ocompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_otetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< ocompidx; ++l)
                    {
                        mtx_otetidx += (statedef()->compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    uint tet_gidx = tri->oTet()->idx();
                    uint tet_lidx = localocomp->getTet_GtoL(tet_gidx);
                    mtx_otetidx+=(tet_lidx*ocompdef->countSpecs());

                    for (uint k=0; k<patchSpecs_N_O; ++k)
                    {
                        uint * olhs = pdef->sreac_lhs_O_bgn(j);
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

                ssolver::Compdef * icompdef = pdef->icompdef();
                if (icompdef != 0)
                {
                    // Sanity check
                    assert(icompdef == tri->iTet()->compdef());
                    uint icompidx = icompdef->gidx();
                    Comp * localicomp = pComps[icompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_itetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< icompidx; ++l)
                    {
                        mtx_itetidx += (statedef()->compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    uint tet_gidx = tri->iTet()->idx();
                    uint tet_lidx = localicomp->getTet_GtoL(tet_gidx);
                    mtx_itetidx+=(tet_lidx*icompdef->countSpecs());

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

                ssolver::Compdef * ocompdef = pdef->ocompdef();
                if (ocompdef != 0)
                {
                    // Sanity check
                    assert(ocompdef == tri->oTet()->compdef());
                    uint ocompidx = ocompdef->gidx();
                    Comp * localocomp = pComps[ocompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_otetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< ocompidx; ++l)
                    {
                        mtx_otetidx += (statedef()->compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    uint tet_gidx = tri->oTet()->idx();
                    uint tet_lidx = localocomp->getTet_GtoL(tet_gidx);
                    mtx_otetidx+=(tet_lidx*ocompdef->countSpecs());

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

                if (icompdef != 0)
                {
                    uint icompidx = icompdef->gidx();
                    Comp * localicomp = pComps[icompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_itetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< icompidx; ++l)
                    {
                        mtx_itetidx += (statedef()->compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    uint tet_gidx = tri->iTet()->idx();
                    uint tet_lidx = localicomp->getTet_GtoL(tet_gidx);
                    mtx_itetidx+=(tet_lidx*icompdef->countSpecs());

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

                if (ocompdef != 0)
                {
                    uint ocompidx = ocompdef->gidx();
                    Comp * localocomp = pComps[ocompidx];
                    // marker for correct position of inner tetrahedron in matrix
                    uint mtx_otetidx = 0;
                    /// step up marker to correct comp
                    for (uint l=0; l< ocompidx; ++l)
                    {
                        mtx_otetidx += (statedef()->compdef(l)->countSpecs())*(pComps[l]->countTets());
                    }
                    // We need to loop up to the correct tet
                    uint tet_gidx = tri->oTet()->idx();
                    uint tet_lidx = localocomp->getTet_GtoL(tet_gidx);
                    mtx_otetidx+=(tet_lidx*ocompdef->countSpecs());

                    for (uint k=0; k<patchSpecs_N_O; ++k)
                    {
                        uint * olhs = pdef->vdepsreac_lhs_O_bgn(j);
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
                        assert(tri_base != 0);

                        // Find the correct location in the 'matrix'
                        for (uint l = 0; l < 3; ++l)
                        {
                            double dcst = pdef->dcst(j);
                            Tri * tri_neighb = tri_base->nextTri(l);
                            if (tri_neighb == 0) continue;

                            double dist = tri_base->dist(l);
                            double dccst = (tri_base->length(l)*dcst)/(tri_base->area()*dist);

                            // Find the right location in the 'matrix'
                            uint spec_base_idx = 0;

                            for (uint m=0; m<i; ++m)
                            {
                                spec_base_idx += (statedef()->patchdef(m)->countSpecs())*(pPatches[m]->countTris());
                            }

                            // Need to convert neighbour to local index
                            uint tri_neighb_gidx = tri_neighb->idx();
                            uint tri_neighb_lidx = patch->getTri_GtoL(tri_neighb_gidx);

                            // At the moment spec_base_idx points to the start position
                            // in the 'matrix' for this patch (start of the tris)
                            uint spec_neighb_idx = spec_base_idx+(patchSpecs_N_S*tri_neighb_lidx)+k;

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
    assert(spec_gidx == pSpecs_tot);
    assert(reac_gidx == pReacs_tot);

    ////////// Now to setup the cvode structures ///////////

    pCVodeState = new CVodeState(pSpecs_tot, 10000, 1.0e-3, 1.0e-3);

    if (efflag() == true) _setupEField();

}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::setTemp(double t)
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

void stode::TetODE::_setupEField(void)
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

    using namespace steps::solver::efield;

    pEField = make_EField<dVSolverBanded>();
    pEField->initMesh(nefverts(), pEFVerts, neftris(), pEFTris, neftets(), pEFTets, memb->_getOpt_method(), memb->_getOpt_file_name(), memb->_getSearch_percent());
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::setTolerances(double atol, double rtol)
{
    pCVodeState->setTolerances(atol, rtol);
    pTolsset = true;
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_addTet(uint tetidx,
                            steps::tetode::Comp * comp, double vol,
                             double a1, double a2, double a3, double a4,
                             double d1, double d2, double d3, double d4,
                             int tet0, int tet1, int tet2, int tet3)
{
    steps::solver::Compdef * compdef  = comp->def();
    stode::Tet * localtet = new stode::Tet(tetidx, compdef, vol, a1, a2, a3, a4, d1, d2, d3, d4,
                                         tet0, tet1, tet2, tet3);
    assert(localtet != 0);
    assert(tetidx < pTets.size());
    assert(pTets[tetidx] == 0);
    pTets[tetidx] = localtet;
    comp->addTet(localtet);
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_addTri(uint triidx, steps::tetode::Patch * patch, double area,
        double l0, double l1, double l2, double d0, double d1, double d2,
        int tinner, int touter, int tri0, int tri1, int tri2)
{
    steps::solver::Patchdef * patchdef = patch->def();
    stode::Tri * tri = new stode::Tri(triidx, patchdef, area, l0, l1, l2, d0, d1, d2,  tinner, touter, tri0, tri1, tri2);
    assert(tri != 0);
    assert (triidx < pTris.size());
    assert (pTris[triidx] == 0);
    pTris[triidx] = tri;
    patch->addTri(tri);
}

////////////////////////////////////////////////////////////////////////////////

uint stode::TetODE::_addComp(steps::solver::Compdef * cdef)
{
    stode::Comp * comp = new Comp(cdef);
    assert(comp != 0);
    uint compidx = pComps.size();
    pComps.push_back(comp);
    return compidx;
}

////////////////////////////////////////////////////////////////////////////////

uint stode::TetODE::_addPatch(steps::solver::Patchdef * pdef)
{
    /* Comp * icomp = 0;
     Comp * ocomp = 0;
     if (pdef->icompdef()) icomp = pCompMap[pdef->icompdef()];
     if (pdef->ocompdef()) ocomp = pCompMap[pdef->ocompdef()];
     */
    stode::Patch * patch = new Patch(pdef);
    assert(patch != 0);
    uint patchidx = pPatches.size();
    pPatches.push_back(patch);
    return patchidx;
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_ccst(double kcst, double vol, uint order)
{
    double vscale = 1.0e3 * vol * steps::math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;

    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;
    return kcst * pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_ccst2D(double kcst, double area, uint order)
{
    double vscale = area * steps::math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;
    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;
    return kcst * pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////

void stode::TetODE::advance(double adv)
{
    if (adv < 0.0)
        throw steps::ArgErr("Time to advance cannot be negative.");

    double endtime = statedef()->time() + adv;
    run(endtime);
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::run(double endtime)
{
    if (endtime < statedef()->time())
        throw steps::ArgErr("Endtime is before current simulation time.");

    if (endtime == 0.0) return;

    int flag = 0;

    if (not pInitialised)
    {
        if (not pTolsset)
        {
            std::cout << "Warning: tolerances have not been set and will ";
            std::cout << "retain default values\n";
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
        if (efflag() == true)
        {
            // THis flag is set to true initially so this is a good place to set the VDepSReac constants


            TriPVecCI eftri_end = pEFTris_vec.end();
            uint tlidx = 0;
            for (TriPVecCI eft = pEFTris_vec.begin(); eft != eftri_end; ++eft)
            {
                uint tgidx = pEFTri_LtoG[tlidx];

                double voltage = pEField->getTriV(tlidx);

                Tri * tri = pTris[tgidx];

                steps::solver::Patchdef * pdef = tri->patchdef();

                uint pidx = pdef->gidx();

                uint nvdepsreacs = pdef->countVDepSReacs();

                for (uint vlidx = 0; vlidx < nvdepsreacs; ++vlidx)
                {
                    double kcst = pdef->vdepsreacdef(vlidx)->getVDepK(voltage);
                    // Calculate the reaction constant
                    double ccst=0.0;
                    if (pdef->vdepsreacdef(vlidx)->surf_surf() == false)
                    {
                        double vol=0.0;
                        if (pdef->vdepsreacdef(vlidx)->inside() == true)
                        {
                            assert(pdef->icompdef() != 0);
                            Tet * itet = tri->iTet();
                            assert(itet!=0);
                            vol = itet->vol();
                        }
                        else
                        {
                            assert(pdef->ocompdef() != 0);
                            Tet * otet = tri->oTet();
                            assert(otet != 0);
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

                    uint ncomps = pComps.size();
                    for (uint i=0; i< ncomps; ++i)
                    {
                        spec_idx += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
                        reac_idx += (statedef()->compdef(i)->countReacs())*(pComps[i]->countTets());
                        reac_idx += (statedef()->compdef(i)->countDiffs())*(pComps[i]->countTets());
                    }

                    // Step up to the correct patch:
                    for (uint i=0; i < pidx; ++i)
                    {
                        spec_idx += (statedef()->patchdef(i)->countSpecs())*(pPatches[i]->countTris());
                        reac_idx += (statedef()->patchdef(i)->countSReacs())*(pPatches[i]->countTris());
                        reac_idx += (statedef()->patchdef(i)->countVDepSReacs())*(pPatches[i]->countTris());
                        reac_idx += (statedef()->patchdef(i)->countSurfDiffs())*(pPatches[i]->countTris());
                    }

                    uint patchSpecs_N = pdef->countSpecs();
                    uint patchSReacs_N = pdef->countSReacs();

                    uint patchVDepSReacs_N = pdef->countVDepSReacs();

                    // Step up the index to the right triangle
                    Patch * localpatch = pPatches[pidx];
                    uint tri_lpidx = localpatch->getTri_GtoL(tgidx);

                    // Step up indices to the correct triangle
                    spec_idx += (patchSpecs_N*tri_lpidx);
                    reac_idx += (patchSReacs_N*tri_lpidx);
                    // The following is right because SReacs and VDepSReacs are added within the same loop over Tris:
                    reac_idx += (patchVDepSReacs_N*tri_lpidx);

                    // Not necessary because of separate loops:    reac_idx += (patchSDiffs_N*tri_lpidx);

                    // And set the correct reaction index
                    // I think it's right to add one more patch SReac
                    reac_idx += patchSReacs_N;
                    reac_idx+=vlidx;

                    for (uint k=0; k < patchSpecs_N; ++k)
                    {
                        std::vector<stode::structA>::iterator r_end = pSpec_matrixsub[spec_idx+k].end();
                        for (std::vector<stode::structA>::iterator r = pSpec_matrixsub[spec_idx+k].begin(); r!=r_end; ++r)
                        {
                            if ((*r).r_idx == reac_idx) (*r).ccst = ccst;
                        }
                    }

                    // Now the complicated part, which is to change the constants in the inner
                    // and/or outer tets

                    // DO not do this! spec_idx = 0;

                    if (pdef->vdepsreacdef(vlidx)->reqInside())
                    {
                        Tet * itet = tri->iTet();
                        assert(itet != 0);

                        // Fetch the global index of the comp
                        uint icidx = itet->compdef()->gidx();

                        Comp * icomp = pComps[icidx];

                        // First step up the species index to the correct comp
                        uint spec_idx_i = 0;
                        for (uint i=0; i< icidx; ++i)
                        {
                            spec_idx_i += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
                        }
                        uint icompSpecs_N = statedef()->compdef(icidx)->countSpecs();
                        assert(icompSpecs_N == pdef->countSpecs_I());

                        // Fetch the comps-local index for the tet
                        uint tlidx = pComps[icidx]->getTet_GtoL(itet->idx());

                        // Step up the indices to the right tet:
                        spec_idx_i += (icompSpecs_N*tlidx);
                        for(uint k=0; k< icompSpecs_N; ++k)
                        {
                            std::vector<stode::structA>::iterator r_end = pSpec_matrixsub[spec_idx_i+k].end();
                            for (std::vector<stode::structA>::iterator r = pSpec_matrixsub[spec_idx_i+k].begin(); r!=r_end; ++r)
                            {
                                if ((*r).r_idx == reac_idx) (*r).ccst = ccst;
                            }
                        }
                    }

                    if (pdef->vdepsreacdef(vlidx)->reqOutside())
                    {
                        Tet * otet = tri->oTet();
                        assert(otet != 0);

                        // Fetch the global index of the comp
                        uint ocidx = otet->compdef()->gidx();

                        Comp * ocomp = pComps[ocidx];

                        // First step up the species index to the correct comp
                        uint spec_idx_o = 0;
                        for (uint i=0; i< ocidx; ++i)
                        {
                            spec_idx_o += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
                        }
                        uint ocompSpecs_N = statedef()->compdef(ocidx)->countSpecs();
                        assert(ocompSpecs_N == pdef->countSpecs_O());

                        // Fetch the comps-local index for the tet
                        uint tlidx = pComps[ocidx]->getTet_GtoL(otet->idx());

                        // Step up the indices to the right tet:
                        spec_idx_o += (ocompSpecs_N*tlidx);
                        for(uint k=0; k< ocompSpecs_N; ++k)
                        {
                            std::vector<stode::structA>::iterator r_end = pSpec_matrixsub[spec_idx_o+k].end();
                            for (std::vector<stode::structA>::iterator r = pSpec_matrixsub[spec_idx_o+k].begin(); r!=r_end; ++r)
                            {
                                if ((*r).r_idx == reac_idx) (*r).ccst = ccst;
                            }
                        }
                    }
                }
                tlidx += 1;
            }
        }

        flag = pCVodeState->reinit(statedef()->time());

        pReinit = false;
    }

    flag = pCVodeState->run(endtime);

    if (flag != CV_SUCCESS)
    {
          std::ostringstream os;
          os << "\nCVODE iteration failed\n\n";
          throw steps::SysErr(os.str());
    }

    if (efflag() == true)
    {
        double dt = endtime - statedef()->time();

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

    statedef()->setTime(endtime);
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getCompVol(uint cidx) const
{
    assert(cidx < statedef()->countComps());
    assert (statedef()->countComps() == pComps.size());
    stode::Comp * comp = pComps[cidx];
    assert(comp != 0);
    return comp->vol();
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getCompAmount(uint cidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getCompCount(cidx, sidx);
    return (count / smath::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setCompAmount(uint cidx, uint sidx, double a)
{
    // convert amount in mols to number of molecules
    double a2 = a * steps::math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setCompCount(cidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getCompConc(uint cidx, uint sidx) const
{
    // the following methods do all the necessary argument checking
    double count = _getCompCount(cidx, sidx);
    double vol = _getCompVol(cidx);

    return count/ (1.0e3 * vol * steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setCompConc(uint cidx, uint sidx, double c)
{
    // the following method does cidx argument checking
    double vol = _getCompVol(cidx);

    double count = c * (1.0e3 * vol * steps::math::AVOGADRO);

    // the following method does more argument checking
    _setCompCount(cidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool stode::TetODE::_getCompClamped(uint cidx, uint sidx) const
{
    std::ostringstream os;
    os << "getCompClamped not implemented for steps::solver::TetODE solver";
    throw steps::NotImplErr(os.str());

    return false;
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setCompClamped(uint cidx, uint sidx, bool b)
{
    std::ostringstream os;
    os << "setCompClamped not implemented for steps::solver::TetODE solver";
    throw steps::NotImplErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getCompReacK(uint cidx, uint ridx) const
{
    std::ostringstream os;
    os << "getCompReacK not implemented for steps::solver::TetODE solver";
    throw steps::NotImplErr(os.str());

    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setCompReacK(uint cidx, uint ridx, double kf)
{
    assert(cidx < statedef()->countComps());
    assert (statedef()->countComps() == pComps.size());
    stode::Comp * comp = pComps[cidx];
    assert(comp != 0);

    TetPVecCI tet_end = comp->endTet();
    for (TetPVecCI tet = comp->bgnTet(); tet != tet_end; ++tet)
    {
        // The following method will check the ridx argument
        _setTetReacK((*tet)->idx(), ridx, kf);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool stode::TetODE::_getCompReacActive(uint cidx, uint ridx) const
{
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setCompReacActive(uint cidx, uint ridx, bool a)
{
    std::ostringstream os;
    os << "setCompReacActive not implemented for steps::solver::TetODE solver";
    throw steps::NotImplErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getPatchArea(uint pidx) const
{
    assert(pidx < statedef()->countPatches());
    assert (statedef()->countPatches() == pPatches.size());
    stode::Patch * patch = pPatches[pidx];
    assert(patch != 0);
    return patch->area();
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getPatchAmount(uint pidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getPatchCount(pidx, sidx);
    return (count / steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setPatchAmount(uint pidx, uint sidx, double a)
{
    assert(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * steps::math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setPatchCount(pidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

bool stode::TetODE::_getPatchClamped(uint pidx, uint sidx) const
{
    std::ostringstream os;
    os << "getPatchClamped not implemented for steps::solver::TetODE solver";
    throw steps::NotImplErr(os.str());

    return false;
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setPatchClamped(uint pidx, uint sidx, bool buf)
{
    std::ostringstream os;
    os << "setPatchClamped not implemented for steps::solver::TetODE solver";
    throw steps::NotImplErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getPatchSReacK(uint pidx, uint ridx) const
{
    std::ostringstream os;
    os << "getPatchSReacK not implemented for steps::solver::TetODE solver";
    throw steps::NotImplErr(os.str());
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setPatchSReacK(uint pidx, uint ridx, double kf)
{
    assert(pidx < statedef()->countPatches());
    assert (statedef()->countPatches() == pPatches.size());
    stode::Patch * patch = pPatches[pidx];
    assert(patch != 0);

    TriPVecCI tri_end = patch->endTri();
    for (TriPVecCI tri = patch->bgnTri(); tri != tri_end; ++tri)
    {
        // The following method will check the ridx argument
        _setTriSReacK((*tri)->idx(), ridx, kf);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool stode::TetODE::_getPatchSReacActive(uint pidx, uint ridx) const
{
    std::ostringstream os;
    os << "getPatchSReacActive not implemented for steps::solver::TetODE solver";
    throw steps::NotImplErr(os.str());
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setPatchSReacActive(uint pidx, uint ridx, bool a)
{
    std::ostringstream os;
    os << "setPatchSReacActive not implemented for steps::solver::TetODE solver";
    throw steps::NotImplErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getCompCount(uint cidx, uint sidx) const
{
    assert(cidx < statedef()->countComps());
    assert(sidx < statedef()->countSpecs());
    Compdef * comp = statedef()->compdef(cidx);
    assert(comp != 0);
    uint slidx = comp->specG2L(sidx);
    if (slidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    uint idx = 0;
    /// step up marker to correct comp
    for (uint i=0; i< cidx; ++i)
    {
        idx += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
    }

    uint comp_nspecs = comp->countSpecs();
    uint ntets = pComps[cidx]->countTets();

    double count = 0.0;

    assert((idx+(((ntets-1)*comp_nspecs)+slidx)) < pSpecs_tot);
    for (uint i=0; i< ntets; ++i)
    {
        count += Ith(pCVodeState->y_cvode, idx+((i*comp_nspecs)+slidx));
    }

    return count;
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setCompCount(uint cidx, uint sidx, double n)
{

    assert(cidx < statedef()->countComps());
    assert(sidx < statedef()->countSpecs());
    Compdef * comp = statedef()->compdef(cidx);
    assert(comp != 0);
    uint slidx = comp->specG2L(sidx);
    if (slidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in compartment.\n";
        throw steps::ArgErr(os.str());
    }

    uint idx = 0;
    /// step up marker to correct comp
    for (uint i=0; i< cidx; ++i)
    {
        idx += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
    }

    Comp* localcomp = pComps[cidx];
    uint ntets = localcomp->countTets();
    uint comp_nspecs = comp->countSpecs();

    double comp_vol = localcomp->vol();

    assert((idx+(((ntets-1)*comp_nspecs)+slidx)) < pSpecs_tot);

    for (uint i=0; i< ntets; ++i)
    {
        double tetvol = localcomp->getTet(i)->vol();
        Ith(pCVodeState->y_cvode, idx+((i*comp_nspecs)+slidx)) = n*(tetvol/comp_vol);
    }
    // Reinitialise CVode structures
    pReinit = true;
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getPatchCount(uint pidx, uint sidx) const
{
    assert(pidx < statedef()->countPatches());
    assert(sidx < statedef()->countSpecs());
    Patchdef * patch = statedef()->patchdef(pidx);
    assert(patch != 0);

    uint slidx = patch->specG2L(sidx);
    if (slidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }

    uint idx = 0;
    // Step up to correct species index. Comps comes first
    for (uint i=0; i< pComps.size(); ++i)
    {
        idx += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
    }
    // quick sanity check
    assert(idx < pSpecs_tot);

    // Now step up to the correct species index
    for (uint i=0; i< pidx; ++i)
    {
        idx+= (statedef()->patchdef(i)->countSpecs())*(pPatches[i]->countTris());
    }

    stode::Patch * localpatch = pPatches[pidx];
    uint patch_nspecs = patch->countSpecs();
    uint ntris = localpatch->countTris();

    assert(idx+((ntris-1)*patch_nspecs)+slidx < pSpecs_tot);

    double count = 0.0;
    for (uint i = 0; i < ntris; ++i)
    {
        count += Ith(pCVodeState->y_cvode, idx+((i*patch_nspecs)+slidx));
    }

    return count;
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setPatchCount(uint pidx, uint sidx, double n)
{

    assert(pidx < statedef()->countPatches());
    assert(sidx < statedef()->countSpecs());
    Patchdef * patch = statedef()->patchdef(pidx);
    assert(patch != 0);

    uint slidx = patch->specG2L(sidx);
    if (slidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in patch.\n";
        throw steps::ArgErr(os.str());
    }

    uint idx = 0;
    // Step up to correct species index. Comps comes first
    for (uint i=0; i< pComps.size(); ++i)
    {
        idx += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
    }
    // quick sanity check
    assert(idx < pSpecs_tot);

    // Now step up to the correct species index
    for (uint i=0; i< pidx; ++i)
    {
        idx+= (statedef()->patchdef(i)->countSpecs())*(pPatches[i]->countTris());
    }

    stode::Patch * localpatch = pPatches[pidx];
    uint patch_nspecs = patch->countSpecs();
    uint ntris = localpatch->countTris();

    assert((idx+((ntris-1)*patch_nspecs+slidx)) < pSpecs_tot);

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

double stode::TetODE::_getTetCount(uint tidx, uint sidx) const
{
    assert (sidx < statedef()->countSpecs());
    assert (tidx < pTets.size());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    Tet * tet = pTets[tidx];

    Compdef * comp = tet->compdef();
    uint cidx = comp->gidx();

    uint lsidx = comp->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    Comp * localcomp = pComps[cidx];
    uint tet_lcidx = localcomp->getTet_GtoL(tidx);

    uint idx = 0;
    /// step up marker to correct comp
    for (uint i=0; i< cidx; ++i)
    {
        idx += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
    }

    assert((idx + (comp->countSpecs()*tet_lcidx) + lsidx) < pSpecs_tot);

    return Ith(pCVodeState->y_cvode, idx + (comp->countSpecs()*tet_lcidx) + lsidx);
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getTetConc(uint tidx, uint sidx) const
{
    // Following does arg checking on tidx and sidx
    double count = _getTetCount(tidx, sidx);
    double vol = pTets[tidx]->vol();

    return count/(1.0e3*vol* steps::math::AVOGADRO);

}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setTetCount(uint tidx, uint sidx, double n)
{

    assert (sidx < statedef()->countSpecs());
    assert (tidx < pTets.size());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    Tet * tet = pTets[tidx];

    Compdef * comp = tet->compdef();
    uint cidx = comp->gidx();

    uint lsidx = comp->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    Comp * localcomp = pComps[cidx];
    uint tet_lcidx = localcomp->getTet_GtoL(tidx);

    uint idx = 0;
    /// step up marker to correct comp
    for (uint i=0; i< cidx; ++i)
    {
        idx += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
    }

    assert((idx + (comp->countSpecs()*tet_lcidx) + lsidx) < pSpecs_tot);

    Ith(pCVodeState->y_cvode, idx + (comp->countSpecs()*tet_lcidx) + lsidx) = n;

    // Reinitialise CVode structures
    pReinit = true;
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getTetReacK(uint tidx, uint ridx) const
{
    std::ostringstream os;
    os << "getTetReacK not implemented for steps::solver::TetODE solver";
    throw steps::NotImplErr(os.str());
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setTetReacK(uint tidx, uint ridx, double kf)
{
    assert (tidx < pTets.size());
    assert (ridx < statedef()->countReacs());

    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    stode::Tet * tet = pTets[tidx];
    Compdef * comp = tet->compdef();

    uint lridx = comp->reacG2L(ridx);
    if (lridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Reaction undefined in tetrahedron.\n";
        throw steps::ArgErr(os.str());
    }

    // Fetch the global index of the comp
    uint cidx = tet->compdef()->gidx();

    // Now the tricky part
    // First step up the species and reaction indices to the correct comp
    uint reac_idx = 0;
    uint spec_idx = 0;
    for (uint i=0; i< cidx; ++i)
    {
        spec_idx += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
        // Diffusion rules are counted as 'reacs' too
        reac_idx += (statedef()->compdef(i)->countReacs())*(pComps[i]->countTets());
        reac_idx += (statedef()->compdef(i)->countDiffs())*(pComps[i]->countTets());
    }

    uint compSpecs_N = comp->countSpecs();
    uint compReacs_N = comp->countReacs();


    uint tlidx = pComps[cidx]->getTet_GtoL(tidx);
    // Step up the indices to the right tet:
    spec_idx += (compSpecs_N*tlidx);
    reac_idx += (compReacs_N*tlidx);

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

    for(uint k=0; k< compSpecs_N; ++k)
    {
        std::vector<stode::structA>::iterator r_end = pSpec_matrixsub[spec_idx+k].end();
        for (std::vector<stode::structA>::iterator r = pSpec_matrixsub[spec_idx+k].begin(); r!=r_end; ++r)
        {
            if ((*r).r_idx == reac_idx) (*r).ccst = ccst;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setTetConc(uint tidx, uint sidx, double c)
{
    assert (tidx < pTets.size());
    if (pTets[tidx] == 0)
    {
        std::ostringstream os;
        os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
        throw steps::ArgErr(os.str());
    }

    double vol = pTets[tidx]->vol();

    double count = c * (1.0e3 * vol * steps::math::AVOGADRO);

    // Further (and repeated) arg checking done by next method
    _setTetCount(tidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getTetAmount(uint tidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getTetCount(tidx, sidx);
    return (count / smath::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setTetAmount(uint tidx, uint sidx, double a)
{
    // convert amount in mols to number of molecules
    double a2 = a * steps::math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setTetCount(tidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getTriCount(uint tidx, uint sidx) const
{
    assert(sidx < statedef()->countSpecs());
    assert(tidx < pTris.size());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    Tri * tri = pTris[tidx];

    Patchdef * patch = tri->patchdef();
    uint pidx = patch->gidx();

    uint lsidx = patch->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    Patch * localpatch = pPatches[pidx];
    uint tri_lpidx = localpatch->getTri_GtoL(tidx);

    uint idx = 0;
    // Step over compartments' tetrahedrons
    for (uint i=0; i< pComps.size(); ++i)
    {
        idx += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
    }
    // And step up further over previous patches' triangles
    for (uint i=0; i < pidx; ++i)
    {
        idx += (statedef()->patchdef(i)->countSpecs())*(pPatches[i]->countTris());
    }

    assert((idx + (patch->countSpecs()*tri_lpidx) + lsidx) < pSpecs_tot);

    return Ith(pCVodeState->y_cvode, idx + (patch->countSpecs()*tri_lpidx) + lsidx);

}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setTriCount(uint tidx, uint sidx, double n)
{
    assert(sidx < statedef()->countSpecs());
    assert(tidx < pTris.size());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    Tri * tri = pTris[tidx];

    Patchdef * patch = tri->patchdef();
    uint pidx = patch->gidx();

    uint lsidx = patch->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Species undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    Patch * localpatch = pPatches[pidx];
    uint tri_lpidx = localpatch->getTri_GtoL(tidx);

    uint idx = 0;

    uint ncomps = pComps.size();
    // Step over compartments' tetrahedrons
    for (uint i=0; i< ncomps; ++i)
    {
        idx += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
    }
    // And step up further over previous patches' triangles
    for (uint i=0; i < pidx; ++i)
    {
        idx += (statedef()->patchdef(i)->countSpecs())*(pPatches[i]->countTris());
    }

    assert((idx + (patch->countSpecs()*tri_lpidx) + lsidx) < pSpecs_tot);

    Ith(pCVodeState->y_cvode, idx + (patch->countSpecs()*tri_lpidx) + lsidx) = n;

    // Reinitialise CVode structures
    pReinit = true;
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getTriAmount(uint tidx, uint sidx) const
{
    // the following method does all the necessary argument checking
    double count = _getTriCount(tidx, sidx);
    return (count / steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setTriAmount(uint tidx, uint sidx, double a)
{
    assert(a >= 0.0);
    // convert amount in mols to number of molecules
    double a2 = a * steps::math::AVOGADRO;
    // the following method does all the necessary argument checking
    _setTriCount(tidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getTriSReacK(uint tidx, uint ridx) const
{
    std::ostringstream os;
    os << "getTriSReacK not implemented for steps::solver::TetODE solver";
    throw steps::NotImplErr(os.str());
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void stode::TetODE::_setTriSReacK(uint tidx, uint ridx, double kf)
{
    assert(ridx < statedef()->countSReacs());
    assert(tidx < pTris.size());

    if (pTris[tidx] == 0)
    {
        std::ostringstream os;
        os << "Triangle " << tidx << " has not been assigned to a patch.\n";
        throw steps::ArgErr(os.str());
    }

    Tri * tri = pTris[tidx];

    Patchdef * patch = tri->patchdef();

    // Fetch the global index of the patch
    uint pidx = patch->gidx();

    uint lsridx = patch->sreacG2L(ridx);

    if (lsridx == ssolver::LIDX_UNDEFINED)
    {
        std::ostringstream os;
        os << "Surface Reaction undefined in triangle.\n";
        throw steps::ArgErr(os.str());
    }

    // Calculate the reaction constant
    double ccst=0.0;
    if (patch->sreacdef(lsridx)->surf_surf() == false)
    {
        double vol=0.0;
        if (patch->sreacdef(lsridx)->inside() == true)
        {
            assert(patch->icompdef() != 0);
            Tet * itet = tri->iTet();
            assert(itet!=0);
            vol = itet->vol();
        }
        else
        {
            assert(patch->ocompdef() != 0);
            Tet * otet = tri->oTet();
            assert(otet != 0);
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

    uint ncomps = pComps.size();
    for (uint i=0; i< ncomps; ++i)
    {
        spec_idx += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
        reac_idx += (statedef()->compdef(i)->countReacs())*(pComps[i]->countTets());
        reac_idx += (statedef()->compdef(i)->countDiffs())*(pComps[i]->countTets());
    }

    // Step up to the correct patch:
    for (uint i=0; i < pidx; ++i)
    {
        spec_idx += (statedef()->patchdef(i)->countSpecs())*(pPatches[i]->countTris());
        reac_idx += (statedef()->patchdef(i)->countSReacs())*(pPatches[i]->countTris());
        reac_idx += (statedef()->patchdef(i)->countVDepSReacs())*(pPatches[i]->countTris());
        reac_idx += (statedef()->patchdef(i)->countSurfDiffs())*(pPatches[i]->countTris());
    }

    uint patchSpecs_N = patch->countSpecs();
    uint patchSReacs_N = patch->countSReacs();
    uint patchVDepSReacs_N = patch->countVDepSReacs();

    // Step up the index to the right triangle
    Patch * localpatch = pPatches[pidx];
    uint tri_lpidx = localpatch->getTri_GtoL(tidx);

    // Step up indices to the correct triangle
    spec_idx += (patchSpecs_N*tri_lpidx);
    reac_idx += (patchSReacs_N*tri_lpidx);
    // The following is right because SReacs and VDepSReacs are added within the same loop over Tris:
    reac_idx += (patchVDepSReacs_N*tri_lpidx);

    // This is not necessary because separate loop over tris for sreacs, then sdiffs:     reac_idx += (patchSDiffs_N*tri_lpidx);

    // And set the correct reaction index
    reac_idx+=lsridx;

    for (uint k=0; k < patchSpecs_N; ++k)
    {
        std::vector<stode::structA>::iterator r_end = pSpec_matrixsub[spec_idx+k].end();
        for (std::vector<stode::structA>::iterator r = pSpec_matrixsub[spec_idx+k].begin(); r!=r_end; ++r)
        {
            if ((*r).r_idx == reac_idx) (*r).ccst = ccst;
        }
    }

    // Now the complicated part, which is to change the constants in the inner
    // and/or outer tets

    // reset the spec idx
    spec_idx = 0;

    if (patch->sreacdef(lsridx)->reqInside())
    {
        Tet * itet = tri->iTet();
        assert(itet != 0);

        // Fetch the global index of the comp
        uint icidx = itet->compdef()->gidx();

        Comp * icomp = pComps[icidx];

        // First step up the species index to the correct comp
        uint spec_idx = 0;
        for (uint i=0; i< icidx; ++i)
        {
            spec_idx += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
        }
        uint icompSpecs_N = statedef()->compdef(icidx)->countSpecs();
        assert(icompSpecs_N == patch->countSpecs_I());

        // Fetch the comps-local index for the tet
        uint tlidx = pComps[icidx]->getTet_GtoL(itet->idx());

        // Step up the indices to the right tet:
        spec_idx += (icompSpecs_N*tlidx);
        for(uint k=0; k< icompSpecs_N; ++k)
        {
            std::vector<stode::structA>::iterator r_end = pSpec_matrixsub[spec_idx+k].end();
            for (std::vector<stode::structA>::iterator r = pSpec_matrixsub[spec_idx+k].begin(); r!=r_end; ++r)
            {
                if ((*r).r_idx == reac_idx) (*r).ccst = ccst;
            }
        }
    }

    if (patch->sreacdef(lsridx)->reqOutside())
    {
        Tet * otet = tri->oTet();
        assert(otet != 0);

        // Fetch the global index of the comp
        uint ocidx = otet->compdef()->gidx();

        Comp * ocomp = pComps[ocidx];

        // First step up the species index to the correct comp
        uint spec_idx = 0;
        for (uint i=0; i< ocidx; ++i)
        {
            spec_idx += (statedef()->compdef(i)->countSpecs())*(pComps[i]->countTets());
        }
        uint ocompSpecs_N = statedef()->compdef(ocidx)->countSpecs();
        assert(ocompSpecs_N == patch->countSpecs_O());

        // Fetch the comps-local index for the tet
        uint tlidx = pComps[ocidx]->getTet_GtoL(otet->idx());

        // Step up the indices to the right tet:
        spec_idx += (ocompSpecs_N*tlidx);
        for(uint k=0; k< ocompSpecs_N; ++k)
        {
            std::vector<stode::structA>::iterator r_end = pSpec_matrixsub[spec_idx+k].end();
            for (std::vector<stode::structA>::iterator r = pSpec_matrixsub[spec_idx+k].begin(); r!=r_end; ++r)
            {
                if ((*r).r_idx == reac_idx) (*r).ccst = ccst;
            }
        }

    }
}

////////////////////////////////////////////////////////////////////////////////

double stode::TetODE::_getTriArea(uint tidx) const
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

double stode::TetODE::_getTetVol(uint tidx) const
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


double stode::TetODE::_getTetV(uint tidx) const
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

void stode::TetODE::_setTetV(uint tidx, double v)
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

bool stode::TetODE::_getTetVClamped(uint tidx) const
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

void stode::TetODE::_setTetVClamped(uint tidx, bool cl)
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

double stode::TetODE::_getTriV(uint tidx) const
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

void stode::TetODE::_setTriV(uint tidx, double v)
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

bool stode::TetODE::_getTriVClamped(uint tidx) const
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

void stode::TetODE::_setTriVClamped(uint tidx, bool cl)
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
/*
double stode::TetODE::_getTriOhmicI(uint tidx) const
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

double stode::TetODE::_getTriOhmicI(uint tidx, uint ocidx) const
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

double stode::TetODE::_getTriGHKI(uint tidx) const
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

double stode::TetODE::_getTriGHKI(uint tidx, uint ghkidx) const
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
*/
double stode::TetODE::_getTriI(uint tidx) const
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

void stode::TetODE::_setVertIClamp(uint vidx, double cur)
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

void stode::TetODE::_setTriIClamp(uint tidx, double cur)
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

double stode::TetODE::_getVertV(uint vidx) const
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

void stode::TetODE::_setVertV(uint vidx, double v)
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

bool stode::TetODE::_getVertVClamped(uint vidx) const
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

void stode::TetODE::_setVertVClamped(uint vidx, bool cl)
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

void stode::TetODE::_setMembRes(uint midx, double ro, double vrev)
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

void stode::TetODE::_setMembPotential(uint midx, double v)
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

void stode::TetODE::_setMembCapac(uint midx, double cm)
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

void stode::TetODE::_setMembVolRes(uint midx, double ro)
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
// END
