/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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


/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

// STL headers.
#include <string>
#include <sstream>
#include <cassert>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/geom/geom.hpp"
#include "steps/model/model.hpp"
#include "steps/rng/rng.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/specdef.hpp"
#include "steps/solver/compdef.hpp"
#include "steps/solver/patchdef.hpp"
#include "steps/solver/reacdef.hpp"
#include "steps/solver/sreacdef.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/solver/chandef.hpp"
#include "steps/solver/vdeptransdef.hpp"
#include "steps/solver/vdepsreacdef.hpp"
#include "steps/solver/ohmiccurrdef.hpp"
#include "steps/solver/ghkcurrdef.hpp"
#include "steps/solver/diffboundarydef.hpp"
#include "steps/solver/sdiffboundarydef.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

ssolver::Statedef::Statedef(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r)
: pModel(m)
, pGeom(g)
, pRNG(r)
, pTime(0.0)
, pNSteps(0)
, pSpecdefs()
, pChandefs()
, pCompdefs()
, pPatchdefs()
, pReacdefs()
, pSReacdefs()
, pVDepTransdefs()
, pVDepSReacdefs()
, pDiffdefs()
, pOhmicCurrdefs()
, pGHKcurrdefs()
, pDiffBoundarydefs()
, pSDiffBoundarydefs()

{
    AssertLog(pModel != 0);
    AssertLog(pGeom != 0);


    // Create the def objects.
    // NOTE: The order is very important. For example all objects after SpecDef need
    // to know the number of species in the state so SpecDef must be first; CompDef must know what reacs and
    // diffs are in the system, so Reacdef and Diffdef must be created before Compdef.
    // Compdef MUST COME BEFORE Patchef, etc.
    //
    uint nspecs = pModel->_countSpecs();
    AssertLog(nspecs > 0);
    for (uint sidx = 0; sidx < nspecs; ++sidx)
    {
        ssolver::Specdef * specdef = new Specdef(this, sidx,  pModel->_getSpec(sidx));
        AssertLog(specdef != 0);
        pSpecdefs.push_back(specdef);
    }

    uint nchans = pModel->_countChans();
    for (uint chidx = 0; chidx < nchans; ++chidx)
    {
        ssolver::Chandef * chandef = new Chandef(this, chidx, pModel->_getChan(chidx));
        AssertLog(chandef != 0);
        pChandefs.push_back(chandef);
    }

    uint nreacs = pModel->_countReacs();
    for (uint ridx = 0; ridx < nreacs; ++ridx)
    {
        ssolver::Reacdef * reacdef = new Reacdef(this, ridx, pModel->_getReac(ridx));
        AssertLog(reacdef != 0);
        pReacdefs.push_back(reacdef);
    }

    uint nvdiffs = pModel->_countVDiffs();
    for (uint didx = 0; didx < nvdiffs; ++didx)
    {
           ssolver::Diffdef * diffdef = new Diffdef(this, didx, pModel->_getVDiff(didx));
           AssertLog(diffdef != 0);
           pDiffdefs.push_back(diffdef);
    }

    uint nsdiffs = pModel->_countSDiffs();
    for (uint didx = 0; didx < nsdiffs; ++didx)
    {
           ssolver::Diffdef * surfdiffdef = new Diffdef(this, didx, pModel->_getSDiff(didx));
           AssertLog(surfdiffdef != 0);
           pSurfDiffdefs.push_back(surfdiffdef);
    }

    uint nsreacs = pModel->_countSReacs();
    for (uint sridx = 0; sridx < nsreacs; ++sridx)
    {
          ssolver::SReacdef * sreacdef = new SReacdef(this, sridx, pModel->_getSReac(sridx));
           AssertLog(sreacdef != 0);
           pSReacdefs.push_back(sreacdef);
    }

    uint nvdtrans = pModel->_countVDepTrans();
    for (uint vdtidx = 0; vdtidx < nvdtrans; ++vdtidx)
    {
        ssolver::VDepTransdef * vdtdef = new VDepTransdef(this, vdtidx, pModel->_getVDepTrans(vdtidx));
        AssertLog(vdtdef != 0);
        pVDepTransdefs.push_back(vdtdef);
    }

    uint nvdsreacs = pModel->_countVDepSReacs();
    for (uint vdsridx = 0; vdsridx < nvdsreacs; ++vdsridx)
    {
        ssolver::VDepSReacdef * vdsrdef = new VDepSReacdef(this, vdsridx, pModel->_getVDepSReac(vdsridx));
        AssertLog(vdsrdef != 0);
        pVDepSReacdefs.push_back(vdsrdef);
    }

    uint nohmiccurrs = pModel->_countOhmicCurrs();
    for (uint ocidx = 0; ocidx < nohmiccurrs; ++ocidx)
    {
        ssolver::OhmicCurrdef * ocdef = new OhmicCurrdef(this, ocidx, pModel->_getOhmicCurr(ocidx));
        AssertLog(ocdef != 0);
        pOhmicCurrdefs.push_back(ocdef);
    }

    uint nghkcurrs = pModel->_countGHKcurrs();
    for (uint ghkidx = 0; ghkidx < nghkcurrs; ++ghkidx)
    {
        ssolver::GHKcurrdef * ghkdef = new GHKcurrdef(this, ghkidx, pModel->_getGHKcurr(ghkidx));
        AssertLog(ghkdef != 0);
        pGHKcurrdefs.push_back(ghkdef);
    }

    uint ncomps = pGeom->_countComps();
    AssertLog(ncomps >0);
    for (uint cidx = 0; cidx < ncomps; ++cidx)
    {
        ssolver::Compdef * compdef = new Compdef(this, cidx, pGeom->_getComp(cidx));
        AssertLog(compdef != 0);
        pCompdefs.push_back(compdef);
    }

    uint npatches = pGeom->_countPatches();
    for (uint pidx = 0; pidx < npatches; ++pidx)
    {
        ssolver::Patchdef * patchdef = new Patchdef(this, pidx, pGeom->_getPatch(pidx));
        AssertLog(patchdef != 0);
        pPatchdefs.push_back(patchdef);
    }

    if (steps::tetmesh::Tetmesh * tetmesh = dynamic_cast<steps::tetmesh::Tetmesh *>(pGeom))
    {
        uint ndiffbs = tetmesh->_countDiffBoundaries();
        for (uint dbidx = 0; dbidx < ndiffbs; ++dbidx)
        {
            ssolver::DiffBoundarydef * diffboundarydef = new DiffBoundarydef(this, dbidx, tetmesh->_getDiffBoundary(dbidx));
            AssertLog(diffboundarydef != 0);
            pDiffBoundarydefs.push_back(diffboundarydef);
        }

        uint nsdiffbs = tetmesh->_countSDiffBoundaries();
        for (uint sdbidx = 0; sdbidx < nsdiffbs; ++sdbidx)
        {
            ssolver::SDiffBoundarydef * sdiffboundarydef = new SDiffBoundarydef(this, sdbidx, tetmesh->_getSDiffBoundary(sdbidx));
            AssertLog(sdiffboundarydef != 0);
            pSDiffBoundarydefs.push_back(sdiffboundarydef);
        }
    }

    // Now setup all the def objects. This can't be achieved purely with
    // the constructors, e.g.  a patch may need to add species from its
    // surface reactions to inner, outer comp
    // NOTE: Again, order is important.
    //
    for (SpecdefPVecI s = pSpecdefs.begin(); s != pSpecdefs.end(); ++s)
        (*s)->setup();
    for (ChandefPVecI ch = pChandefs.begin(); ch != pChandefs.end(); ++ch)
        (*ch)->setup();
    for (ReacdefPVecI r = pReacdefs.begin(); r != pReacdefs.end(); ++r)
        (*r)->setup();
    for (DiffdefPVecI d = pDiffdefs.begin(); d != pDiffdefs.end(); ++d)
        (*d)->setup();
    for (SurfDiffdefPVecI sd = pSurfDiffdefs.begin(); sd != pSurfDiffdefs.end(); ++sd)
        (*sd)->setup();

    for (SReacdefPVecI sr = pSReacdefs.begin(); sr != pSReacdefs.end(); ++sr)
        (*sr)->setup();
    for (VDepSReacdefPVecI vdsr = pVDepSReacdefs.begin(); vdsr != pVDepSReacdefs.end(); ++vdsr)
        (*vdsr)->setup();
    for (VDepTransdefPVecI vdt = pVDepTransdefs.begin(); vdt != pVDepTransdefs.end(); ++vdt)
        (*vdt)->setup();
    for (OhmicCurrdefPVecI oc = pOhmicCurrdefs.begin(); oc != pOhmicCurrdefs.end(); ++oc)
        (*oc)->setup();
    for (GHKcurrdefPVecI ghk = pGHKcurrdefs.begin(); ghk != pGHKcurrdefs.end(); ++ghk)
        (*ghk)->setup();
    for (CompdefPVecI c = pCompdefs.begin(); c != pCompdefs.end(); ++c)
        (*c)->setup_references();
    for (PatchdefPVecI p = pPatchdefs.begin(); p != pPatchdefs.end(); ++p)
        (*p)->setup_references();

    // Make local indices for species, (surface) reactions, diffusion rules
    // in compartments then patches. Separate method necessary since e.g.
    // Patchdef::setup_references can add species to Compdef
    for (CompdefPVecI c = pCompdefs.begin(); c != pCompdefs.end(); ++c)
        (*c)->setup_indices();
    for (PatchdefPVecI p = pPatchdefs.begin(); p != pPatchdefs.end(); ++p)
        (*p)->setup_indices();

    for (DiffBoundaryDefPVecI db = pDiffBoundarydefs.begin(); db != pDiffBoundarydefs.end(); ++db)
        (*db)->setup();
    for (SDiffBoundaryDefPVecI sdb = pSDiffBoundarydefs.begin(); sdb != pSDiffBoundarydefs.end(); ++sdb)
        (*sdb)->setup();
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Statedef::~Statedef()
{
    CompdefPVecCI c_end = pCompdefs.end();
    for (CompdefPVecCI c = pCompdefs.begin(); c != c_end; ++c) delete *c;

    PatchdefPVecCI p_end = pPatchdefs.end();
    for (PatchdefPVecCI p = pPatchdefs.begin(); p != p_end; ++p) delete *p;

    DiffBoundarydefPVecCI db_end = pDiffBoundarydefs.end();
    for (DiffBoundarydefPVecCI db = pDiffBoundarydefs.begin(); db != db_end; ++db) delete *db;

    SDiffBoundarydefPVecCI sdb_end = pSDiffBoundarydefs.end();
    for (SDiffBoundarydefPVecCI sdb = pSDiffBoundarydefs.begin(); sdb != sdb_end; ++sdb) delete *sdb;

    ReacdefPVecCI r_end = pReacdefs.end();
    for (ReacdefPVecCI r = pReacdefs.begin(); r != r_end; ++r) delete *r;

    SReacdefPVecCI sr_end = pSReacdefs.end();
    for (SReacdefPVecCI sr = pSReacdefs.begin(); sr != sr_end; ++sr) delete *sr;

    SurfDiffdefPVecCI sd_end = pSurfDiffdefs.end();
    for (SurfDiffdefPVecCI sd = pSurfDiffdefs.begin(); sd != sd_end; ++sd) delete *sd;

    DiffdefPVecCI d_end = pDiffdefs.end();
    for (DiffdefPVecCI d = pDiffdefs.begin(); d != d_end; ++d) delete *d;

    ChandefPVecCI ch_end = pChandefs.end();
    for (ChandefPVecCI ch = pChandefs.begin(); ch != ch_end; ++ch) delete *ch;

    VDepTransdefPVecCI vdt_end = pVDepTransdefs.end();
    for (VDepTransdefPVecCI vdt = pVDepTransdefs.begin(); vdt != vdt_end; ++vdt) delete *vdt;

    VDepSReacdefPVecCI vsr_end = pVDepSReacdefs.end();
    for (VDepSReacdefPVecCI vsr = pVDepSReacdefs.begin(); vsr != vsr_end; ++vsr) delete *vsr;

    OhmicCurrdefPVecCI oc_end = pOhmicCurrdefs.end();
    for (OhmicCurrdefPVecCI oc = pOhmicCurrdefs.begin(); oc != oc_end; ++oc) delete *oc;

    GHKcurrdefPVecCI ghk_end = pGHKcurrdefs.end();
    for (GHKcurrdefPVecCI ghk = pGHKcurrdefs.begin(); ghk != ghk_end; ++ghk) delete *ghk;

    SpecdefPVecCI s_end = pSpecdefs.end();
    for (SpecdefPVecCI s = pSpecdefs.begin(); s != s_end; ++s) delete *s;

}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Statedef::checkpoint(std::fstream & cp_file)
{

    SpecdefPVecCI s_end = pSpecdefs.end();
    for (SpecdefPVecCI s = pSpecdefs.begin(); s != s_end; ++s) {
        (*s)->checkpoint(cp_file);
    }

    ChandefPVecCI ch_end = pChandefs.end();
    for (ChandefPVecCI ch = pChandefs.begin(); ch != ch_end; ++ch) {
        (*ch)->checkpoint(cp_file);
    }

    CompdefPVecCI c_end = pCompdefs.end();
    for (CompdefPVecCI c = pCompdefs.begin(); c != c_end; ++c) {
        (*c)->checkpoint(cp_file);
    }

    PatchdefPVecCI p_end = pPatchdefs.end();
    for (PatchdefPVecCI p = pPatchdefs.begin(); p != p_end; ++p) {
        (*p)->checkpoint(cp_file);
    }

    ReacdefPVecCI r_end = pReacdefs.end();
    for (ReacdefPVecCI r = pReacdefs.begin(); r != r_end; ++r) {
        (*r)->checkpoint(cp_file);
    }

    SReacdefPVecCI sr_end = pSReacdefs.end();
    for (SReacdefPVecCI sr = pSReacdefs.begin(); sr != sr_end; ++sr) {
        (*sr)->checkpoint(cp_file);
    }

    DiffdefPVecCI d_end = pDiffdefs.end();
    for (DiffdefPVecCI d = pDiffdefs.begin(); d != d_end; ++d) {
        (*d)->checkpoint(cp_file);
    }

    SurfDiffdefPVecCI sd_end = pSurfDiffdefs.end();
    for (SurfDiffdefPVecCI sd = pSurfDiffdefs.begin(); sd != sd_end; ++sd) {
        (*sd)->checkpoint(cp_file);
    }

    DiffBoundarydefPVecCI db_end = pDiffBoundarydefs.end();
    for (DiffBoundarydefPVecCI db = pDiffBoundarydefs.begin(); db != db_end; ++db) {
        (*db)->checkpoint(cp_file);
    }

    SDiffBoundarydefPVecCI sdb_end = pSDiffBoundarydefs.end();
    for (SDiffBoundarydefPVecCI sdb = pSDiffBoundarydefs.begin(); sdb != sdb_end; ++sdb) {
        (*sdb)->checkpoint(cp_file);
    }

    VDepTransdefPVecCI vdeptrans_end = pVDepTransdefs.end();
    for (VDepTransdefPVecCI vdt = pVDepTransdefs.begin(); vdt != vdeptrans_end; ++vdt) {
        (*vdt)->checkpoint(cp_file);
    }

    VDepSReacdefPVecCI vdsr_end = pVDepSReacdefs.end();
    for (VDepSReacdefPVecCI vdsr = pVDepSReacdefs.begin(); vdsr != vdsr_end; ++vdsr) {
        (*vdsr)->checkpoint(cp_file);
    }

    OhmicCurrdefPVecCI oc_end = pOhmicCurrdefs.end();
    for (OhmicCurrdefPVecCI oc = pOhmicCurrdefs.begin(); oc != oc_end; ++oc) {
        (*oc)->checkpoint(cp_file);
    }

    GHKcurrdefPVecCI ghkc_end = pGHKcurrdefs.end();
    for (GHKcurrdefPVecCI ghkc = pGHKcurrdefs.begin(); ghkc != ghkc_end; ++ghkc) {
        (*ghkc)->checkpoint(cp_file);
    }
    cp_file.write((char*)&pTime, sizeof(double));
    cp_file.write((char*)&pNSteps, sizeof(uint));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Statedef::restore(std::fstream & cp_file)
{

    SpecdefPVecCI s_end = pSpecdefs.end();
    for (SpecdefPVecCI s = pSpecdefs.begin(); s != s_end; ++s) {
        (*s)->restore(cp_file);
    }

    ChandefPVecCI ch_end = pChandefs.end();
    for (ChandefPVecCI ch = pChandefs.begin(); ch != ch_end; ++ch) {
        (*ch)->restore(cp_file);
    }

    CompdefPVecCI c_end = pCompdefs.end();
    for (CompdefPVecCI c = pCompdefs.begin(); c != c_end; ++c) {
        (*c)->restore(cp_file);
    }

    PatchdefPVecCI p_end = pPatchdefs.end();
    for (PatchdefPVecCI p = pPatchdefs.begin(); p != p_end; ++p) {
        (*p)->restore(cp_file);
    }

    ReacdefPVecCI r_end = pReacdefs.end();
    for (ReacdefPVecCI r = pReacdefs.begin(); r != r_end; ++r) {
        (*r)->restore(cp_file);
    }

    SReacdefPVecCI sr_end = pSReacdefs.end();
    for (SReacdefPVecCI sr = pSReacdefs.begin(); sr != sr_end; ++sr) {
        (*sr)->restore(cp_file);
    }

    DiffdefPVecCI d_end = pDiffdefs.end();
    for (DiffdefPVecCI d = pDiffdefs.begin(); d != d_end; ++d) {
        (*d)->restore(cp_file);
    }

    SurfDiffdefPVecCI sd_end = pSurfDiffdefs.end();
    for (SurfDiffdefPVecCI sd = pSurfDiffdefs.begin(); sd != sd_end; ++sd) {
        (*sd)->restore(cp_file);
    }

    DiffBoundarydefPVecCI db_end = pDiffBoundarydefs.end();
    for (DiffBoundarydefPVecCI db = pDiffBoundarydefs.begin(); db != db_end; ++db) {
        (*db)->restore(cp_file);
    }

    SDiffBoundarydefPVecCI sdb_end = pSDiffBoundarydefs.end();
    for (SDiffBoundarydefPVecCI sdb = pSDiffBoundarydefs.begin(); sdb != sdb_end; ++sdb) {
        (*sdb)->restore(cp_file);
    }

    VDepTransdefPVecCI vdeptrans_end = pVDepTransdefs.end();
    for (VDepTransdefPVecCI vdt = pVDepTransdefs.begin(); vdt != vdeptrans_end; ++vdt) {
        (*vdt)->restore(cp_file);
    }

    VDepSReacdefPVecCI vdsr_end = pVDepSReacdefs.end();
    for (VDepSReacdefPVecCI vdsr = pVDepSReacdefs.begin(); vdsr != vdsr_end; ++vdsr) {
        (*vdsr)->restore(cp_file);
    }

    OhmicCurrdefPVecCI oc_end = pOhmicCurrdefs.end();
    for (OhmicCurrdefPVecCI oc = pOhmicCurrdefs.begin(); oc != oc_end; ++oc) {
        (*oc)->restore(cp_file);
    }

    GHKcurrdefPVecCI ghkc_end = pGHKcurrdefs.end();
    for (GHKcurrdefPVecCI ghkc = pGHKcurrdefs.begin(); ghkc != ghkc_end; ++ghkc) {
        (*ghkc)->restore(cp_file);
    }

    cp_file.read((char*)&pTime, sizeof(double));
    cp_file.read((char*)&pNSteps, sizeof(uint));
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Compdef * ssolver::Statedef::compdef(uint gidx) const
{
    AssertLog(gidx < pCompdefs.size());
    return pCompdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getCompIdx(std::string const & c) const
{
    uint maxcidx = pCompdefs.size();
    AssertLog(maxcidx > 0);
    AssertLog(maxcidx == pGeom->_countComps());
    uint cidx = 0;
    while(cidx < maxcidx)
    {
        if (c == pGeom->_getComp(cidx)->getID()) return cidx;
        ++cidx;
    }
    std::ostringstream os;
    os << "Geometry does not contain comp with string identifier '" << c << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getCompIdx(steps::wm::Comp * comp) const
{
    uint maxcidx = pCompdefs.size();
    AssertLog(maxcidx > 0);
    AssertLog(maxcidx == pGeom->_countComps());
    uint cidx = 0;
    while(cidx < maxcidx)
    {
        if (comp == pGeom->_getComp(cidx)) return cidx;
        ++cidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Patchdef * ssolver::Statedef::patchdef(uint gidx) const
{
    AssertLog(gidx < pPatchdefs.size());
    return pPatchdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getPatchIdx(std::string const & p) const
{
    uint maxpidx = pPatchdefs.size();
    AssertLog(maxpidx == pGeom->_countPatches());
    uint pidx = 0;
    while(pidx < maxpidx)
    {
        if (p == pGeom->_getPatch(pidx)->getID()) return pidx;
        ++pidx;
    }
    std::ostringstream os;
    os << "Geometry does not contain patch with string identifier '" << p << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getPatchIdx(steps::wm::Patch * patch) const
{
    uint maxpidx = pPatchdefs.size();
    AssertLog(maxpidx == pGeom->_countPatches());
    uint pidx = 0;
    while(pidx < maxpidx)
    {
        if (patch == pGeom->_getPatch(pidx)) return pidx;
        ++pidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Specdef * ssolver::Statedef::specdef(uint gidx) const
{
    AssertLog(gidx < pSpecdefs.size());
    return pSpecdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getSpecIdx(std::string const & s) const
{
    uint maxsidx = pSpecdefs.size();
    AssertLog(maxsidx > 0);
    AssertLog(maxsidx == pModel->_countSpecs());
    uint sidx = 0;
    while(sidx < maxsidx)
    {
        if (s == pModel->_getSpec(sidx)->getID()) return sidx;
        ++sidx;
    }
    std::ostringstream os;
    os << "Model does not contain species with string identifier '" << s << "'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getSpecIdx(steps::model::Spec * spec) const
{
    uint maxsidx = pSpecdefs.size();
    AssertLog(maxsidx > 0);
    AssertLog(maxsidx == pModel->_countSpecs());
    uint sidx = 0;
    while(sidx < maxsidx)
    {
        if (spec == pModel->_getSpec(sidx)) return sidx;
        ++sidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Reacdef * ssolver::Statedef::reacdef(uint gidx) const
{
    AssertLog(gidx < pReacdefs.size());
    return pReacdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getReacIdx(std::string const & r) const
{
    uint maxridx = pReacdefs.size();
    AssertLog(maxridx == pModel->_countReacs());
    uint ridx = 0;
    while(ridx < maxridx)
    {
        if (r == pModel->_getReac(ridx)->getID()) return ridx;
        ++ridx;
    }
    std::ostringstream os;
    os << "Model does not contain reac with string identifier '" << r <<"'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getReacIdx(steps::model::Reac * reac) const
{
    uint maxridx = pReacdefs.size();
    AssertLog(maxridx == pModel->_countReacs());
    uint ridx = 0;
    while(ridx < maxridx)
    {
        if (reac == pModel->_getReac(ridx)) return ridx;
        ++ridx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::SReacdef * ssolver::Statedef::sreacdef(uint gidx) const
{
    AssertLog(gidx < pSReacdefs.size());
    return pSReacdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getSReacIdx(std::string const & sr) const
{
    uint maxsridx = pSReacdefs.size();
    AssertLog(maxsridx == pModel->_countSReacs());
    uint sridx = 0;
    while(sridx < maxsridx)
    {
        if (sr == pModel->_getSReac(sridx)->getID()) return sridx;
        ++sridx;
    }
    std::ostringstream os;
    os << "Model does not contain sreac with string identifier '" << sr <<"'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getSReacIdx(steps::model::SReac * sreac) const
{
    uint maxsridx = pSReacdefs.size();
    AssertLog(maxsridx == pModel->_countSReacs());
    uint sridx = 0;
    while(sridx < maxsridx)
    {
        if (sreac == pModel->_getSReac(sridx)) return sridx;
        ++sridx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Diffdef * ssolver::Statedef::diffdef(uint gidx) const
{
    AssertLog(gidx < pDiffdefs.size());
    return pDiffdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getDiffIdx(std::string const & d) const
{
    uint maxdidx = pDiffdefs.size();
    AssertLog(maxdidx == pModel->_countVDiffs());
    uint didx = 0;
    while(didx < maxdidx)
    {
        if (d == pModel->_getVDiff(didx)->getID()) return didx;
        ++didx;
    }
    std::ostringstream os;
    os << "Model does not contain diff with string identifier '" << d <<"'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getDiffIdx(steps::model::Diff * diff) const
{
    uint maxdidx = pDiffdefs.size();
    AssertLog(maxdidx == pModel->_countVDiffs());
    uint didx = 0;
    while(didx < maxdidx)
    {
        if (diff == pModel->_getVDiff(didx)) return didx;
        ++didx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Diffdef * ssolver::Statedef::surfdiffdef(uint gidx) const
{
    AssertLog(gidx < pSurfDiffdefs.size());
    return pSurfDiffdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getSurfDiffIdx(std::string const & d) const
{
    uint maxdidx = pSurfDiffdefs.size();
    AssertLog(maxdidx == pModel->_countSDiffs());
    uint didx = 0;
    while(didx < maxdidx)
    {
        if (d == pModel->_getSDiff(didx)->getID()) return didx;
        ++didx;
    }
    std::ostringstream os;
    os << "Model does not contain diff with string identifier '" << d <<"'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getSurfDiffIdx(steps::model::Diff * diff) const
{
    uint maxdidx = pSurfDiffdefs.size();
    AssertLog(maxdidx == pModel->_countSDiffs());
    uint didx = 0;
    while(didx < maxdidx)
    {
        if (diff == pModel->_getSDiff(didx)) return didx;
        ++didx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::OhmicCurrdef * ssolver::Statedef::ohmiccurrdef(uint gidx) const
{
    AssertLog(gidx < pOhmicCurrdefs.size());
    return pOhmicCurrdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getOhmicCurrIdx(std::string const & oc) const
{
    uint maxocidx = pOhmicCurrdefs.size();
    AssertLog(maxocidx == pModel->_countOhmicCurrs());
    uint ocidx = 0;
    while(ocidx < maxocidx)
    {
        if (oc == pModel->_getOhmicCurr(ocidx)->getID()) return ocidx;
        ++ocidx;
    }
    std::ostringstream os;
    os << "Model does not contain ohmic current with string identifier '" << oc <<"'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getOhmicCurrIdx(steps::model::OhmicCurr * ohmiccurr) const
{
    uint maxocidx = pOhmicCurrdefs.size();
    AssertLog(maxocidx == pModel->_countOhmicCurrs());
    uint ocidx = 0;
    while(ocidx < maxocidx)
    {
        if (ohmiccurr == pModel->_getOhmicCurr(ocidx)) return ocidx;
        ++ocidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::VDepTransdef * ssolver::Statedef::vdeptransdef(uint gidx) const
{
    AssertLog(gidx < pVDepTransdefs.size());
    return pVDepTransdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getVDepTransIdx(std::string const & vdt) const
{
    uint maxvdtidx = pVDepTransdefs.size();
    AssertLog(maxvdtidx == pModel->_countVDepTrans());
    uint vdtidx = 0;
    while(vdtidx < maxvdtidx)
    {
        if (vdt == pModel->_getVDepTrans(vdtidx)->getID()) return vdtidx;
        ++vdtidx;
    }
    std::ostringstream os;
    os << "Model does not contain voltage-dependent transition with string identifier '" << vdt <<"'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getVDepTransIdx(steps::model::VDepTrans * vdeptrans) const
{
    uint maxvdtidx = pVDepTransdefs.size();
    AssertLog(maxvdtidx == pModel->_countVDepTrans());
    uint vdtidx = 0;
    while(vdtidx < maxvdtidx)
    {
        if (vdeptrans == pModel->_getVDepTrans(vdtidx)) return vdtidx;
        ++vdtidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::VDepSReacdef * ssolver::Statedef::vdepsreacdef(uint gidx) const
{
    AssertLog(gidx < pVDepSReacdefs.size());
    return pVDepSReacdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getVDepSReacIdx(std::string const & vdsr) const
{
    uint maxvdsridx = pVDepSReacdefs.size();
    AssertLog(maxvdsridx == pModel->_countVDepSReacs());
    uint vdsridx = 0;
    while(vdsridx < maxvdsridx)
    {
        if (vdsr == pModel->_getVDepSReac(vdsridx)->getID()) return vdsridx;
        ++vdsridx;
    }
    std::ostringstream os;
    os << "Model does not contain voltage-dependent reaction with string identifier '" << vdsr <<"'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getVDepSReacIdx(steps::model::VDepSReac * vdepsreac) const
{
    uint maxvdsridx = pVDepSReacdefs.size();
    AssertLog(maxvdsridx == pModel->_countVDepSReacs());
    uint vdsridx = 0;
    while(vdsridx < maxvdsridx)
    {
        if (vdepsreac == pModel->_getVDepSReac(vdsridx)) return vdsridx;
        ++vdsridx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////
\
ssolver::GHKcurrdef * ssolver::Statedef::ghkcurrdef(uint gidx) const
{
    AssertLog(gidx < pGHKcurrdefs.size());
    return pGHKcurrdefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getGHKcurrIdx(std::string const & ghk) const
{
    uint maxghkidx = pGHKcurrdefs.size();
    AssertLog(maxghkidx == pModel->_countGHKcurrs());
    uint ghkidx = 0;
    while(ghkidx < maxghkidx)
    {
        if (ghk == pModel->_getGHKcurr(ghkidx)->getID()) return ghkidx;
        ++ghkidx;
    }
    std::ostringstream os;
    os << "Model does not contain ghk current with string identifier '" << ghk <<"'.";
    ArgErrLog(os.str());
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getGHKcurrIdx(steps::model::GHKcurr * ghkcurr) const
{
    uint maxghkidx = pGHKcurrdefs.size();
    AssertLog(maxghkidx == pModel->_countGHKcurrs());
    uint ghkidx = 0;
    while(ghkidx < maxghkidx)
    {
        if (ghkcurr == pModel->_getGHKcurr(ghkidx)) return ghkidx;
        ++ghkidx;
    }
    // Argument should be valid so we should not get here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getMembIdx(std::string const & m) const
{
    if (steps::tetmesh::Tetmesh * tetmesh = dynamic_cast<steps::tetmesh::Tetmesh *>(pGeom))
    {
        /* This isn't the right place for this check, rather the setup should check for
         * multiple membranes
        uint nmembs = tetmesh->_countMembs();
        if (nmembs != 1)
        {
            std::ostringstream os;
            os << "Only one Membrane may exist in simulation";
            ArgErrLog(os.str());
        }
        */
        steps::tetmesh::Memb * membrane = tetmesh->_getMemb(0);
        if (m == membrane->getID())
        {
            return 0;
        }
        else
        {
            std::ostringstream os;
            os << "Geometry does not contain membrane with string identifier '" <<  m <<"'.";
            ArgErrLog(os.str());
        }

    }
    else
    {
        std::ostringstream os;
        os << "Membrane methods not available with well-mixed geometry";
        ArgErrLog(os.str());
    }

}

////////////////////////////////////////////////////////////////////////////////

ssolver::DiffBoundarydef * ssolver::Statedef::diffboundarydef(uint gidx) const
{
    AssertLog(gidx < pDiffBoundarydefs.size());
    return pDiffBoundarydefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getDiffBoundaryIdx(std::string const & d) const
{
    uint maxdidx = pDiffBoundarydefs.size();
    if (steps::tetmesh::Tetmesh * tetmesh = dynamic_cast<steps::tetmesh::Tetmesh *>(pGeom))
    {
        AssertLog(maxdidx == tetmesh->_countDiffBoundaries());
        uint didx = 0;
        while(didx < maxdidx)
        {
            if (d == tetmesh->_getDiffBoundary(didx)->getID()) return didx;
            ++didx;
        }
        std::ostringstream os;
        os << "Geometry does not contain diff boundary with string identifier '" << d <<"'.";
        ArgErrLog(os.str());
    }
    else
    {
        std::ostringstream os;
        os << "Diffusion boundary methods not available with well-mixed geometry";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getDiffBoundaryIdx(steps::tetmesh::DiffBoundary * diffb) const
{
    uint maxdidx = pDiffBoundarydefs.size();
    if (steps::tetmesh::Tetmesh * tetmesh = dynamic_cast<steps::tetmesh::Tetmesh *>(pGeom))
    {
        AssertLog(maxdidx == tetmesh->_countDiffBoundaries());
        uint didx = 0;
        while(didx < maxdidx)
        {
            if (diffb == tetmesh->_getDiffBoundary(didx)) return didx;
            ++didx;
        }
        // Argument should be valid so we should not get here
        AssertLog(false);
    }
    else
    {
        std::ostringstream os;
        os << "Diffusion boundary methods not available with well-mixed geometry";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

ssolver::SDiffBoundarydef * ssolver::Statedef::sdiffboundarydef(uint gidx) const
{
    AssertLog(gidx < pSDiffBoundarydefs.size());
    return pSDiffBoundarydefs[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getSDiffBoundaryIdx(std::string const & sd) const
{
    uint maxsdidx = pSDiffBoundarydefs.size();
    if (steps::tetmesh::Tetmesh * tetmesh = dynamic_cast<steps::tetmesh::Tetmesh *>(pGeom))
    {
        AssertLog(maxsdidx == tetmesh->_countSDiffBoundaries());
        uint sdidx = 0;
        while(sdidx < maxsdidx)
        {
            if (sd == tetmesh->_getSDiffBoundary(sdidx)->getID()) return sdidx;
            ++sdidx;
        }
        std::ostringstream os;
        os << "Geometry does not contain surface diffusion boundary with string identifier '" << sd <<"'.";
        ArgErrLog(os.str());
    }
    else
    {
        std::ostringstream os;
        os << "Surface Diffusion Boundary methods not available with well-mixed geometry";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Statedef::getSDiffBoundaryIdx(steps::tetmesh::SDiffBoundary * sdiffb) const
{
    uint maxsdidx = pSDiffBoundarydefs.size();
    if (steps::tetmesh::Tetmesh * tetmesh = dynamic_cast<steps::tetmesh::Tetmesh *>(pGeom))
    {
        AssertLog(maxsdidx == tetmesh->_countSDiffBoundaries());
        uint sdidx = 0;
        while(sdidx < maxsdidx)
        {
            if (sdiffb == tetmesh->_getSDiffBoundary(sdidx)) return sdidx;
            ++sdidx;
        }
        // Argument should be valid so we should not get here
        AssertLog(false);
    }
    else
    {
        std::ostringstream os;
        os << "Surface Diffusion Boundary methods not available with well-mixed geometry";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Statedef::setTime(double t)
{
    AssertLog(t >= 0.0);
    pTime = t;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Statedef::incTime(double dt)
{
    AssertLog(dt >= 0.0);
    pTime += dt;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Statedef::incNSteps(uint i)
{
    AssertLog(i != 0);
    pNSteps += i;
}

////////////////////////////////////////////////////////////////////////////////

// END
