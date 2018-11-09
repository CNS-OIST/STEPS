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
#include <cassert>
#include <sstream>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/geom/patch.hpp"
#include "steps/model/ghkcurr.hpp"
#include "steps/model/ohmiccurr.hpp"
#include "steps/model/sreac.hpp"
#include "steps/model/vdepsreac.hpp"
#include "steps/model/vdeptrans.hpp"
#include "steps/solver/compdef.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/solver/ghkcurrdef.hpp"
#include "steps/solver/ohmiccurrdef.hpp"
#include "steps/solver/patchdef.hpp"
#include "steps/solver/sreacdef.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/types.hpp"
#include "steps/solver/vdepsreacdef.hpp"
#include "steps/solver/vdeptransdef.hpp"

// logging
#include "easylogging++.h"

namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

ssolver::Patchdef::Patchdef(Statedef * sd, uint idx, steps::wm::Patch * p)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pArea()
, pPssys()
, pIcomp(nullptr)
, pOcomp(nullptr)
, pInner(nullptr)
, pOuter(nullptr)
, pSetupRefsdone(false)
, pSetupIndsdone(false)
, pSpecsN_I(0)
, pSpecsN_S(0)
, pSpecsN_O(0)
, pSpec_G2L(nullptr)
, pSpec_L2G(nullptr)
, pPoolCount(nullptr)
, pPoolFlags(nullptr)
, pSReacsN(0)
, pSReac_G2L(nullptr)
, pSReac_L2G(nullptr)
, pSReacKcst(nullptr)
, pSReacFlags(nullptr)
, pSReac_DEP_I_Spec(nullptr)
, pSReac_DEP_S_Spec(nullptr)
, pSReac_DEP_O_Spec(nullptr)
, pSReac_LHS_I_Spec(nullptr)
, pSReac_LHS_S_Spec(nullptr)
, pSReac_LHS_O_Spec(nullptr)
, pSReac_UPD_I_Spec(nullptr)
, pSReac_UPD_S_Spec(nullptr)
, pSReac_UPD_O_Spec(nullptr)
, pSurfDiffsN(0)
, pSurfDiff_G2L(nullptr)
, pSurfDiff_L2G(nullptr)
, pSurfDiff_DEP_Spec(nullptr)
, pSurfDiff_LIG(nullptr)
, pOhmicCurrsN(0)
, pOhmicCurr_G2L(nullptr)
, pOhmicCurr_L2G(nullptr)
, pOhmicCurr_DEP_Spec(nullptr)
, pOhmicCurr_CHANSTATE(nullptr)
, pGHKcurrsN(0)
, pGHKcurr_G2L(nullptr)
, pGHKcurr_DEP_Spec(nullptr)
, pGHKcurr_CHANSTATE(nullptr)
, pGHKcurr_ION(nullptr)
, pVDepTransN(0)
, pVDepTrans_G2L(nullptr)
, pVDepTrans_L2G(nullptr)
, pVDepTrans_DEP_Spec(nullptr)
, pVDepTrans_SRCCHANSTATE(nullptr)
, pVDepTrans_DSTCHANSTATE(nullptr)
, pVDepSReacsN(0)
, pVDepSReac_G2L(nullptr)
, pVDepSReac_L2G(nullptr)
, pVDepSReac_DEP_I_Spec(nullptr)
, pVDepSReac_DEP_S_Spec(nullptr)
, pVDepSReac_DEP_O_Spec(nullptr)
, pVDepSReac_LHS_I_Spec(nullptr)
, pVDepSReac_LHS_S_Spec(nullptr)
, pVDepSReac_LHS_O_Spec(nullptr)
, pVDepSReac_UPD_I_Spec(nullptr)
, pVDepSReac_UPD_S_Spec(nullptr)
, pVDepSReac_UPD_O_Spec(nullptr)
{
    AssertLog(pStatedef != 0);
    AssertLog(p != 0);

    pName = p->getID();
    pArea = p->getArea();
    pPssys = p->getSurfsys();
    pIcomp = p->getIComp();
    pOcomp = p->getOComp();

    uint nspecs = pStatedef->countSpecs();
    if (nspecs > 0)
    {
        pSpec_G2L = new uint[nspecs];
        std::fill_n(pSpec_G2L, nspecs, LIDX_UNDEFINED);
    }

    uint nsreacs = pStatedef->countSReacs();
    if (nsreacs > 0)
    {
        pSReac_G2L = new uint[nsreacs];
        std::fill_n(pSReac_G2L, nsreacs, LIDX_UNDEFINED);
    }

    uint nsdiffs = pStatedef->countSurfDiffs();
    if (nsdiffs > 0)
    {
        pSurfDiff_G2L = new uint[nsdiffs];
        std::fill_n(pSurfDiff_G2L, nsdiffs, LIDX_UNDEFINED);
    }

    uint nohmiccurrs = pStatedef->countOhmicCurrs();
    if (nohmiccurrs > 0)
    {
        pOhmicCurr_G2L = new uint[nohmiccurrs];
        std::fill_n(pOhmicCurr_G2L, nohmiccurrs, LIDX_UNDEFINED);
    }

    uint nghkcurrs = pStatedef->countGHKcurrs();
    if (nghkcurrs > 0)
    {
        pGHKcurr_G2L = new uint[nghkcurrs];
        std::fill_n(pGHKcurr_G2L, nghkcurrs, LIDX_UNDEFINED);
    }

    uint nvdeptrans = pStatedef->countVDepTrans();
    if (nvdeptrans > 0)
    {
        pVDepTrans_G2L = new uint[nvdeptrans];
        std::fill_n(pVDepTrans_G2L, nvdeptrans, LIDX_UNDEFINED);
    }

    uint nvdepsreacs = pStatedef->countVDepSReacs();
    if (nvdepsreacs > 0)
    {
        pVDepSReac_G2L = new uint[nvdepsreacs];
        std::fill_n(pVDepSReac_G2L, nvdepsreacs, LIDX_UNDEFINED);
    }
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Patchdef::~Patchdef()
{
    if (pStatedef->countSpecs() > 0) delete[] pSpec_G2L;
    if (pStatedef->countSReacs() > 0) { delete[] pSReac_G2L;
}
    if (pStatedef->countSurfDiffs() > 0) { delete[] pSurfDiff_G2L;
}
    if (pStatedef->countOhmicCurrs() > 0) { delete[] pOhmicCurr_G2L;
}
    if (pStatedef->countGHKcurrs() > 0) { delete[] pGHKcurr_G2L;
}
    if (pStatedef->countVDepTrans() > 0) { delete[] pVDepTrans_G2L;
}
    if (pStatedef->countVDepSReacs() > 0) { delete[] pVDepSReac_G2L;
}

    if (pSpecsN_S != 0) { delete[] pSpec_L2G;
}

    if (pSReacsN != 0)
    {
        delete[] pSReac_L2G;
        delete[] pSReac_DEP_S_Spec;
        delete[] pSReac_LHS_S_Spec;
        delete[] pSReac_UPD_S_Spec;
        delete[] pSReac_DEP_I_Spec;
        delete[] pSReac_LHS_I_Spec;
        delete[] pSReac_UPD_I_Spec;
        if (pOuter != nullptr)
        {
            delete[] pSReac_DEP_O_Spec;
            delete[] pSReac_LHS_O_Spec;
            delete[] pSReac_UPD_O_Spec;
        }
    }

    if (pVDepSReacsN != 0)
    {
        delete[] pVDepSReac_L2G;
        delete[] pVDepSReac_DEP_S_Spec;
        delete[] pVDepSReac_LHS_S_Spec;
        delete[] pVDepSReac_UPD_S_Spec;
        delete[] pVDepSReac_DEP_I_Spec;
        delete[] pVDepSReac_LHS_I_Spec;
        delete[] pVDepSReac_UPD_I_Spec;
        if (pOuter != nullptr)
        {
            delete[] pVDepSReac_DEP_O_Spec;
            delete[] pVDepSReac_LHS_O_Spec;
            delete[] pVDepSReac_UPD_O_Spec;
        }
    }

    if (pOhmicCurrsN != 0)
    {
        delete[] pOhmicCurr_L2G;
        delete[] pOhmicCurr_DEP_Spec;
        delete[] pOhmicCurr_CHANSTATE;
    }

    if (pGHKcurrsN != 0)
    {
        delete[] pGHKcurr_L2G;
        delete[] pGHKcurr_DEP_Spec;
        delete[] pGHKcurr_CHANSTATE;
        delete[] pGHKcurr_ION;
    }

    if (pVDepTransN != 0)
    {
        delete[] pVDepTrans_L2G;
        delete[] pVDepTrans_DEP_Spec;
        delete[] pVDepTrans_SRCCHANSTATE;
        delete[] pVDepTrans_DSTCHANSTATE;
    }

    if (pSpecsN_S != 0)
    {
        delete[] pPoolCount;
        delete[] pPoolFlags;
    }

    if (pSurfDiffsN != 0)
    {
        delete[] pSurfDiff_L2G;
        delete[] pSurfDiff_DEP_Spec;
        delete[] pSurfDiff_LIG;
        delete[] pSurfDiffDcst;

    }

    if (pSReacsN != 0)
    {
        delete[] pSReacFlags;
        delete[] pSReacKcst;
    }
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)pPoolCount, sizeof(double) * pSpecsN_S);
    cp_file.write((char*)pPoolFlags, sizeof(uint) * pSpecsN_S);
    cp_file.write((char*)pSReacKcst, sizeof(double) * pSReacsN);
    cp_file.write((char*)pSReacFlags, sizeof(uint) * pSReacsN);
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::restore(std::fstream & cp_file)
{
    cp_file.read((char*)pPoolCount, sizeof(double) * pSpecsN_S);
    cp_file.read((char*)pPoolFlags, sizeof(uint) * pSpecsN_S);
    cp_file.read((char*)pSReacKcst, sizeof(double) * pSReacsN);
    cp_file.read((char*)pSReacFlags, sizeof(uint) * pSReacsN);
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::setup_references()
{
    AssertLog(pSetupRefsdone == false);
    AssertLog(pSetupIndsdone == false);

    // first find the inner and outer comps of this patch
    AssertLog(pIcomp != 0);
    uint icompidx = pStatedef->getCompIdx(pIcomp);                ///// bit long-winded, add new method to statedef??
    pInner = pStatedef->compdef(icompidx);
    if (pOcomp != nullptr)
    {
        uint ocompidx = pStatedef->getCompIdx(pOcomp);
        pOuter = pStatedef->compdef(ocompidx);
    }

    const uint ngspecs = pStatedef->countSpecs();
    const uint ngsreacs = pStatedef->countSReacs();
    const uint ngvdepsreacs = pStatedef->countVDepSReacs();
    const uint ngohmiccurrs = pStatedef->countOhmicCurrs();
    const uint ngghkcurrs = pStatedef->countGHKcurrs();
    const uint ngvdeptrans = pStatedef->countVDepTrans();
    const uint ngsdiffs = pStatedef->countSurfDiffs();


    if (ngspecs == 0) AssertLog(pSpec_G2L == 0);
    if (ngsreacs == 0) AssertLog(pSReac_G2L == 0);
    if (ngvdepsreacs == 0) AssertLog(pVDepSReac_G2L == 0);
    if (ngohmiccurrs == 0) AssertLog(pOhmicCurr_G2L == 0);
    if (ngghkcurrs == 0) AssertLog(pGHKcurr_G2L == 0);
    if (ngvdeptrans == 0) AssertLog(pVDepTrans_G2L == 0);
    if (ngsdiffs == 0) AssertLog(pSurfDiff_G2L == 0);


    // set up local sreac indices
    std::set<std::string>::const_iterator s_end = pPssys.end();
    for(std::set<std::string>::const_iterator s = pPssys.begin();
        s != s_end; ++s)
    {
        std::map<std::string, steps::model::SReac *> ssreacs = pStatedef->model()->getSurfsys(*s)->_getAllSReacs();
        if (ngsreacs == 0) AssertLog(ssreacs.empty() == true);
        std::map<std::string, steps::model::SReac*>::const_iterator sr_end = ssreacs.end();
        for(std::map<std::string, steps::model::SReac *>::const_iterator sr = ssreacs.begin(); sr != sr_end; ++sr)
        {
            uint gidx = pStatedef->getSReacIdx((sr->second));
            AssertLog(gidx < ngsreacs);
            if(pSReac_G2L[gidx] != LIDX_UNDEFINED) continue;
            pSReac_G2L[gidx] = pSReacsN++;
        }

           std::map<std::string, steps::model::Diff *> sdiffs = pStatedef->model()->getSurfsys(*s)->_getAllDiffs();
        if (ngsdiffs == 0) AssertLog(sdiffs.empty() == true);
           std::map<std::string, steps::model::Diff*>::const_iterator sd_end = sdiffs.end();
           for (std::map<std::string, steps::model::Diff*>::const_iterator sd = sdiffs.begin(); sd != sd_end; ++sd)
           {
               uint gidx = pStatedef->getSurfDiffIdx((sd->second));
               AssertLog(gidx < ngsdiffs);
               if (pSurfDiff_G2L[gidx] != LIDX_UNDEFINED) continue;
               pSurfDiff_G2L[gidx] = pSurfDiffsN++;
           }

        std::map<std::string, steps::model::VDepSReac *> vdssreacs = pStatedef->model()->getSurfsys(*s)->_getAllVDepSReacs();
        if (ngvdepsreacs == 0) AssertLog(vdssreacs.empty() == true);
        std::map<std::string, steps::model::VDepSReac*>::const_iterator vdsr_end = vdssreacs.end();
        for(std::map<std::string, steps::model::VDepSReac *>::const_iterator vdsr = vdssreacs.begin(); vdsr != vdsr_end; ++vdsr)
        {
            uint gidx = pStatedef->getVDepSReacIdx((vdsr->second));
            AssertLog(gidx < ngvdepsreacs);
            if(pVDepSReac_G2L[gidx] != LIDX_UNDEFINED) continue;
            pVDepSReac_G2L[gidx] = pVDepSReacsN++;
        }

        std::map<std::string, steps::model::OhmicCurr *> ocs = pStatedef->model()->getSurfsys(*s)->_getAllOhmicCurrs();
        if (ngohmiccurrs == 0) AssertLog(ocs.empty() == true);
        std::map<std::string, steps::model::OhmicCurr *>::const_iterator oc_end = ocs.end();
        for(std::map<std::string, steps::model::OhmicCurr *>::const_iterator oc = ocs.begin(); oc != oc_end; ++oc)
        {
            uint gidx = pStatedef->getOhmicCurrIdx((oc->second));
            AssertLog(gidx < ngohmiccurrs);
            if(pOhmicCurr_G2L[gidx] != LIDX_UNDEFINED) continue;
            pOhmicCurr_G2L[gidx] = pOhmicCurrsN++;
        }

        std::map<std::string, steps::model::GHKcurr *> ghks = pStatedef->model()->getSurfsys(*s)->_getAllGHKcurrs();
        if (ngghkcurrs == 0) AssertLog(ghks.empty() == true);
        std::map<std::string, steps::model::GHKcurr *>::const_iterator ghk_end = ghks.end();
        for(std::map<std::string, steps::model::GHKcurr *>::const_iterator ghk = ghks.begin(); ghk != ghk_end; ++ghk)
        {
            uint gidx = pStatedef->getGHKcurrIdx((ghk->second));
            AssertLog(gidx < ngghkcurrs);
            if(pGHKcurr_G2L[gidx] != LIDX_UNDEFINED) continue;
            pGHKcurr_G2L[gidx] = pGHKcurrsN++;
        }

        std::map<std::string, steps::model::VDepTrans *> vdts = pStatedef->model()->getSurfsys(*s)->_getAllVDepTrans();
        if (ngvdeptrans == 0) AssertLog(vdts.empty() == true);
        std::map<std::string, steps::model::VDepTrans *>::const_iterator vdt_end = vdts.end();
        for(std::map<std::string, steps::model::VDepTrans *>::const_iterator vdt = vdts.begin(); vdt != vdt_end; ++vdt)
        {
            uint gidx = pStatedef->getVDepTransIdx((vdt->second));
            AssertLog(gidx < ngvdeptrans);
            if(pVDepTrans_G2L[gidx] != LIDX_UNDEFINED) continue;
            pVDepTrans_G2L[gidx] = pVDepTransN++;
        }
    }


    // Now add all species that appear in all surface reactions, ohmic currents,
    // ghk currents and voltage-dependent transitions/reactions that can occur
    // on this patch: to the patch, inner or outer compartment.
    for(uint sr = 0; sr < ngsreacs; ++sr)
    {
        if(pSReac_G2L[sr] == LIDX_UNDEFINED) { continue;
}
        SReacdef * srdef = pStatedef->sreacdef(sr);
        AssertLog(srdef != 0);
        for (uint s = 0; s < ngspecs; ++s)
        {
            if (srdef->reqspec_S(s) == true)
            {
                AssertLog(pStatedef->specdef(s) != 0);
                if (pSpec_G2L[s] == LIDX_UNDEFINED) { pSpec_G2L[s] = pSpecsN_S++;
}
            }
            if (srdef->reqspec_I(s) == true)
            {
                AssertLog(pInner != 0);
                pInner->addSpec(s);
            }
            if (srdef->reqspec_O(s) == true)
            {
                if (pOuter == nullptr)
                {
                    std::ostringstream os;
                    os << "Can't add surface reaction '" << srdef->name() << "' to patch '";
                    os << name() << "'. Outer compartment not defined for this patch.";
                    ArgErrLog(os.str());
                }
                pOuter->addSpec(s);
            }
        }
    }

    for (uint sd = 0; sd < ngsdiffs; ++sd)
    {
        if (pSurfDiff_G2L[sd] == LIDX_UNDEFINED) { continue;
}
        Diffdef * sddef = pStatedef->surfdiffdef(sd);
        AssertLog(sddef != 0);
        for (uint s = 0; s < ngspecs; ++s)
        {
            if (sddef->reqspec(s) == true)
            {
                if (pSpec_G2L[s] == LIDX_UNDEFINED) { pSpec_G2L[s] = pSpecsN_S++;
}
            }
        }
    }

    for(uint vdsr = 0; vdsr < ngvdepsreacs; ++vdsr)
    {
        if(pVDepSReac_G2L[vdsr] == LIDX_UNDEFINED) { continue;
}
        VDepSReacdef * vdsrdef = pStatedef->vdepsreacdef(vdsr);
        AssertLog(vdsrdef != 0);
        for (uint s = 0; s < ngspecs; ++s)
        {
            if (vdsrdef->reqspec_S(s) == true)
            {
                AssertLog(pStatedef->specdef(s) != 0);
                if (pSpec_G2L[s] == LIDX_UNDEFINED) { pSpec_G2L[s] = pSpecsN_S++;
}
            }
            if (vdsrdef->reqspec_I(s) == true)
            {
                AssertLog(pInner != 0);
                pInner->addSpec(s);
            }
            if (vdsrdef->reqspec_O(s) == true)
            {
                if (pOuter == nullptr)
                {
                    std::ostringstream os;
                    os << "Can't add voltage-dependent reaction '" << vdsrdef->name() << "' to patch '";
                    os << name() << "'. Outer compartment not defined for this patch.";
                    ArgErrLog(os.str());
                }
                pOuter->addSpec(s);
            }
        }
    }

    for (uint oc = 0; oc < ngohmiccurrs; ++oc)
    {
        if(pOhmicCurr_G2L[oc] == LIDX_UNDEFINED) { continue;
}
        OhmicCurrdef * ocdef = pStatedef->ohmiccurrdef(oc);
        AssertLog(ocdef != 0);
        uint added = 0;
        for (uint s = 0; s < ngspecs; ++s)
        {
            // Add the channel state
            if (ocdef->req(s) == true)
            {
                AssertLog(pStatedef->specdef(s) != 0);
                if (pSpec_G2L[s] == LIDX_UNDEFINED) { pSpec_G2L[s] = pSpecsN_S++;
}
                added +=1;
            }
        }
        // Only one channel state should be added per ohmic current
        AssertLog(added == 1);
    }

    for (uint ghk = 0; ghk < ngghkcurrs; ++ghk)
    {
        if (pGHKcurr_G2L[ghk] == LIDX_UNDEFINED) { continue;
}
        GHKcurrdef * ghkdef = pStatedef->ghkcurrdef(ghk);
        AssertLog(ghkdef != 0);
        uint added = 0;
        for (uint s = 0; s < ngspecs; ++s)
        {
            // Add the channel state
            if (ghkdef->req(s) == true)
            {
                AssertLog(pStatedef->specdef(s) != 0);
                // Only add the channel state, not the volume ion species (that affects the GHK rate)
                if (ghkdef->req_v(s) == false)
                {
                    if (pSpec_G2L[s] == LIDX_UNDEFINED) { pSpec_G2L[s] = pSpecsN_S++;
}
                    added += 1;
                }
            }
            // Add the volume ion species to the inner and outer compartment.
            if (ghkdef->req_v(s) == true)
            {
                AssertLog(pInner != 0);
                pInner->addSpec(s);
                if (pOuter == nullptr)
                {
                    if (ghkdef->voconc() < 0.0)
                    {
                        std::ostringstream os;
                        os << "Can't add GHK current '" << ghkdef->name() << "' to patch '";
                        os << name() << "'. Outer compartment not defined for this patch ";
                        os << "and no virtual concentration has been defined.";
                        ArgErrLog(os.str());
                    }
                }
                else if (ghkdef->voconc() < 0.0) { pOuter->addSpec(s);
}
            }
        }
        // Only one channel state should be added per ghk current
        AssertLog(added == 1);
    }

    for (uint vdt = 0; vdt < ngvdeptrans; ++vdt)
    {
        if (pVDepTrans_G2L[vdt] == LIDX_UNDEFINED) { continue;
}
        VDepTransdef * vdtdef = pStatedef->vdeptransdef(vdt);
        for (uint s = 0; s < ngspecs; ++s)
        {
            if (vdtdef->req(s) == true)
            {
                AssertLog(pStatedef->specdef(s) != 0);
                if (pSpec_G2L[s] == LIDX_UNDEFINED) { pSpec_G2L[s] = pSpecsN_S++;
}
            }
        }
    }

    pSetupRefsdone = true;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::setup_indices()
{
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == false);

    // 1 -- DEAL WITH PATCH SPECIES
    // (Only if any species have been added to the patch)
    //   -> Setup local indices for all species.
    //
    // 2 -- COPY #SPECIES FOR INNER AND OUTER COMPS
    //   (These are required a lot during simulation, so it's sound
    //   to have them ready here to avoid an extra level of pointer
    //   lookup.)
    //
    // 3 -- DEAL WITH PATCH SREAC'S
    // (Only if any surface reactions have been added to the patch)
    //   -> Setup local indices for all surface reactions.
    //   -> The SReac objects have LHS, DEP and UPD vectors expressed in
    //      global species indices. Transform this to local indices:
    //      -> Pre-create pSReac_DEP, _LHS and _UPD vectors
    //          -> Always for surface (_S)
    //          -> For inner comp (_I)
    //          -> For outer comp (_O) if outer comp has been defined
    //      -> Loop over the SReacDef objects added to this patch:
    //          -> Fill out the newly created vectors by appropriate
    //             copying of the vectors defined the SReacDef object.
    //          -> While doing this, check whether everything can be
    //             resolved.
    //
    // 4 -- DEAL WITH PATCH VOLTAGE-DEPENDENT REACTIONS
    // (Only if any vdep reactions have been added to the patch)
    //   -> Setup local indices for all vdep reactions.
    //   -> The VDepSReac objects have LHS, DEP and UPD vectors expressed in
    //      global species indices. Transform this to local indices:
    //      -> Pre-create pVDepSReac_DEP, _LHS and _UPD vectors
    //          -> Always for surface (_S)
    //          -> For inner comp (_I)
    //          -> For outer comp (_O) if outer comp has been defined
    //      -> Loop over the VDepSReacDef objects added to this patch:
    //          -> Fill out the newly created vectors by appropriate
    //             copying of the vectors defined the VDepSReacDef object.
    //          -> While doing this, check whether everything can be
    //             resolved.
    // 5 -- DEAL WITH OHMIC CURRENTS
    // 6 -- DEAL WITH GHK CURRENTS
    // 7 -- DEAL WITH V-DEPENDENT TRANSITIONS

    // 1 -- DEAL WITH PATCH SPECIES
    uint ngspecs = pStatedef->countSpecs();
    if (pSpecsN_S != 0)
    {
        pSpec_L2G = new uint[pSpecsN_S];
        for (uint i = 0; i < ngspecs; ++i)
        {
            uint lidx = pSpec_G2L[i];
            if (lidx == LIDX_UNDEFINED) { continue;
}
            pSpec_L2G[lidx] = i;
        }
    }

    // 2 -- COPY #SPECS FOR INNER AND OUTER COMPS
    if (pInner != nullptr) { pSpecsN_I = pInner->countSpecs();
}
    if (pOuter != nullptr) { pSpecsN_O = pOuter->countSpecs();
}

    // 3 -- DEAL WITH PATCH SREAC'S
    if (pSReacsN != 0)
    {
        // Set up local indices.
        pSReac_L2G = new uint[pSReacsN];
        uint ngsreacs = pStatedef->countSReacs();
        for (uint i = 0; i < ngsreacs; ++i)
        {
            uint lidx = pSReac_G2L[i];
            if (lidx == LIDX_UNDEFINED) { continue;
}
            pSReac_L2G[lidx] = i;
        }

        // Create _DEP, _LHS and _UPD vectors.
        uint arrsize_i = 0;
        uint arrsize_s = pSpecsN_S * pSReacsN;
        uint arrsize_o = 0;
        pSReac_DEP_S_Spec = new int[arrsize_s];
        pSReac_LHS_S_Spec = new uint[arrsize_s];
        pSReac_UPD_S_Spec = new int[arrsize_s];
        std::fill_n(pSReac_DEP_S_Spec, arrsize_s, 0);
        std::fill_n(pSReac_LHS_S_Spec, arrsize_s, 0);
        std::fill_n(pSReac_UPD_S_Spec, arrsize_s, 0);

        AssertLog(pInner != 0); // Inner comp should exist
        {
            arrsize_i = pSpecsN_I * pSReacsN;
            pSReac_DEP_I_Spec = new int[arrsize_i];
            pSReac_LHS_I_Spec = new uint[arrsize_i];
            pSReac_UPD_I_Spec = new int[arrsize_i];
            std::fill_n(pSReac_DEP_I_Spec, arrsize_i, 0);
            std::fill_n(pSReac_LHS_I_Spec, arrsize_i, 0);
            std::fill_n(pSReac_UPD_I_Spec, arrsize_i, 0);
        }
        if (pOuter != nullptr) // Only create if outer comp exists.
        {
            arrsize_o = pSpecsN_O * pSReacsN;
            pSReac_DEP_O_Spec = new int[arrsize_o];
            pSReac_LHS_O_Spec = new uint[arrsize_o];
            pSReac_UPD_O_Spec = new int[arrsize_o];
            std::fill_n(pSReac_DEP_O_Spec, arrsize_o, 0);
            std::fill_n(pSReac_LHS_O_Spec, arrsize_o, 0);
            std::fill_n(pSReac_UPD_O_Spec, arrsize_o, 0);
        }

        // Fill the vectors with all kinds of useful information.
        for (uint ri = 0; ri < pSReacsN; ++ri)
        {
            SReacdef * srdef = sreacdef(ri);

            // Handle surface stuff.
            for (uint si = 0; si < ngspecs; ++si)
            {
                if (srdef->reqspec_S(si) == false) { continue;
}

                // TODO: turn into error check?
                uint sil = pSpec_G2L[si];
                AssertLog(sil != LIDX_UNDEFINED);

                uint aridx = _IDX_SReac_S_Spec(ri, sil);
                pSReac_DEP_S_Spec[aridx] = srdef->dep_S(si);
                pSReac_LHS_S_Spec[aridx] = srdef->lhs_S(si);
                pSReac_UPD_S_Spec[aridx] = srdef->upd_S(si);
            }

            // Handle the inside comp stuff.
            if (srdef->reqInside() == true)
            {
                // TODO: turn into real error check?
                AssertLog(pInner != 0);

                for (uint si = 0; si < ngspecs; ++si)
                {
                    if (srdef->reqspec_I(si) == false) { continue;
}

                    // TODO: turn into error check?
                    uint sil = specG2L_I(si);
                    AssertLog(sil != LIDX_UNDEFINED);

                    uint aridx = _IDX_SReac_I_Spec(ri, sil);
                    pSReac_DEP_I_Spec[aridx] = srdef->dep_I(si);
                    pSReac_LHS_I_Spec[aridx] = srdef->lhs_I(si);
                    pSReac_UPD_I_Spec[aridx] = srdef->upd_I(si);
                }
            }

            // Handle the outside comp stuff.
            if (srdef->reqOutside() == true)
            {
                // TODO: turn into real error check?
                AssertLog(pOuter != 0);

                for (uint si = 0; si < ngspecs; ++si)
                {
                    if (srdef->reqspec_O(si) == false) { continue;
}

                    // TODO: turn into error check?
                    uint sil = specG2L_O(si);
                    AssertLog(sil != LIDX_UNDEFINED);

                    uint aridx = _IDX_SReac_O_Spec(ri, sil);
                    pSReac_DEP_O_Spec[aridx] = srdef->dep_O(si);
                    pSReac_LHS_O_Spec[aridx] = srdef->lhs_O(si);
                    pSReac_UPD_O_Spec[aridx] = srdef->upd_O(si);
                }
            }
        }
    }

    // 3.5 -- DEAL WITH PATCH SURFACE-DIFFUSION
    if (pSurfDiffsN != 0)
    {
        pSurfDiff_L2G = new uint[pSurfDiffsN];
        uint ngsdiffs = pStatedef->countSurfDiffs();

        for (uint i = 0; i < ngsdiffs; ++i)
        {
            uint lidx = pSurfDiff_G2L[i];
            if (lidx == LIDX_UNDEFINED) { continue;
}
            pSurfDiff_L2G[lidx] = i;
        }

        uint arrsize = pSpecsN_S * pSurfDiffsN;
        pSurfDiff_DEP_Spec = new uint[arrsize];
        std::fill_n(pSurfDiff_DEP_Spec, arrsize, 0);
        pSurfDiff_LIG = new uint[pSurfDiffsN];
        for (uint di = 0; di < pSurfDiffsN; ++di)
        {
            Diffdef * sddef = surfdiffdef(di);
            pSurfDiff_LIG[di] = pSpec_G2L[sddef->lig()];
            for (uint si = 0; si < ngspecs; ++si)
            {
                if (sddef->reqspec(si) == false) { continue;
}
                uint sil = pSpec_G2L[si];
                AssertLog(sil != LIDX_UNDEFINED);
                uint aridx = _IDX_SurfDiff_Spec(di, sil);
                pSurfDiff_DEP_Spec[aridx] = sddef->dep(si);
            }
        }
    }

    // 4 -- DEAL WITH PATCH VOLTAGE-DEPENDENT SURFACE REACTIONS
    if (pVDepSReacsN != 0)
    {
        // Set up local indices.
        pVDepSReac_L2G = new uint[pVDepSReacsN];
        uint ngvdsreacs = pStatedef->countVDepSReacs();
        for (uint i = 0; i < ngvdsreacs; ++i)
        {
            uint lidx = pVDepSReac_G2L[i];
            if (lidx == LIDX_UNDEFINED) { continue;
}
            pVDepSReac_L2G[lidx] = i;
        }

        // Create _DEP, _LHS and _UPD vectors.
        uint arrsize_i = 0;
        uint arrsize_s = pSpecsN_S * pVDepSReacsN;
        uint arrsize_o = 0;
        pVDepSReac_DEP_S_Spec = new int[arrsize_s];
        pVDepSReac_LHS_S_Spec = new uint[arrsize_s];
        pVDepSReac_UPD_S_Spec = new int[arrsize_s];
        std::fill_n(pVDepSReac_DEP_S_Spec, arrsize_s, 0);
        std::fill_n(pVDepSReac_LHS_S_Spec, arrsize_s, 0);
        std::fill_n(pVDepSReac_UPD_S_Spec, arrsize_s, 0);

        AssertLog(pInner != 0); // Inner comp should exist
        {
            arrsize_i = pSpecsN_I * pVDepSReacsN;
            pVDepSReac_DEP_I_Spec = new int[arrsize_i];
            pVDepSReac_LHS_I_Spec = new uint[arrsize_i];
            pVDepSReac_UPD_I_Spec = new int[arrsize_i];
            std::fill_n(pVDepSReac_DEP_I_Spec, arrsize_i, 0);
            std::fill_n(pVDepSReac_LHS_I_Spec, arrsize_i, 0);
            std::fill_n(pVDepSReac_UPD_I_Spec, arrsize_i, 0);
        }
        if (pOuter != nullptr) // Only create if outer comp exists.
        {
            arrsize_o = pSpecsN_O * pVDepSReacsN;
            pVDepSReac_DEP_O_Spec = new int[arrsize_o];
            pVDepSReac_LHS_O_Spec = new uint[arrsize_o];
            pVDepSReac_UPD_O_Spec = new int[arrsize_o];
            std::fill_n(pVDepSReac_DEP_O_Spec, arrsize_o, 0);
            std::fill_n(pVDepSReac_LHS_O_Spec, arrsize_o, 0);
            std::fill_n(pVDepSReac_UPD_O_Spec, arrsize_o, 0);
        }

        // Fill the vectors with all kinds of useful information.
        for (uint ri = 0; ri < pVDepSReacsN; ++ri)
        {
            VDepSReacdef * vdsrdef = vdepsreacdef(ri);

            // Handle surface stuff.
            for (uint si = 0; si < ngspecs; ++si)
            {
                if (vdsrdef->reqspec_S(si) == false) { continue;
}

                // TODO: turn into error check?
                uint sil = pSpec_G2L[si];
                AssertLog(sil != LIDX_UNDEFINED);

                uint aridx = _IDX_VDepSReac_S_Spec(ri, sil);
                pVDepSReac_DEP_S_Spec[aridx] = vdsrdef->dep_S(si);
                pVDepSReac_LHS_S_Spec[aridx] = vdsrdef->lhs_S(si);
                pVDepSReac_UPD_S_Spec[aridx] = vdsrdef->upd_S(si);
            }

            // Handle the inside comp stuff.
            if (vdsrdef->reqInside() == true)
            {
                // TODO: turn into real error check?
                AssertLog(pInner != 0);

                for (uint si = 0; si < ngspecs; ++si)
                {
                    if (vdsrdef->reqspec_I(si) == false) { continue;
}

                    // TODO: turn into error check?
                    uint sil = specG2L_I(si);
                    AssertLog(sil != LIDX_UNDEFINED);

                    uint aridx = _IDX_VDepSReac_I_Spec(ri, sil);
                    pVDepSReac_DEP_I_Spec[aridx] = vdsrdef->dep_I(si);
                    pVDepSReac_LHS_I_Spec[aridx] = vdsrdef->lhs_I(si);
                    pVDepSReac_UPD_I_Spec[aridx] = vdsrdef->upd_I(si);
                }
            }

            // Handle the outside comp stuff.
            if (vdsrdef->reqOutside() == true)
            {
                // TODO: turn into real error check?
                AssertLog(pOuter != 0);

                for (uint si = 0; si < ngspecs; ++si)
                {
                    if (vdsrdef->reqspec_O(si) == false) { continue;
}

                    // TODO: turn into error check?
                    uint sil = specG2L_O(si);
                    AssertLog(sil != LIDX_UNDEFINED);

                    uint aridx = _IDX_VDepSReac_O_Spec(ri, sil);
                    pVDepSReac_DEP_O_Spec[aridx] = vdsrdef->dep_O(si);
                    pVDepSReac_LHS_O_Spec[aridx] = vdsrdef->lhs_O(si);
                    pVDepSReac_UPD_O_Spec[aridx] = vdsrdef->upd_O(si);
                }
            }
        }
    }
    // 5 -- DEAL WITH OHMIC CURRENTS
    if (pOhmicCurrsN != 0)
    {
        // Set up local indices.
        pOhmicCurr_L2G= new uint[countOhmicCurrs()];
        uint ngohmiccurrs = pStatedef->countOhmicCurrs();
        for (uint i = 0; i < ngohmiccurrs; ++i)
        {
            uint lidx = ohmiccurrG2L(i);
            if (lidx == LIDX_UNDEFINED) { continue;
}
            pOhmicCurr_L2G[lidx] = i;
        }

        // Create local _DEP and _CHANSTATE vectors
        uint arrsize1 = countSpecs() * countOhmicCurrs();
        uint arrsize2 = countOhmicCurrs();
        pOhmicCurr_DEP_Spec = new int[arrsize1];
        pOhmicCurr_CHANSTATE = new uint[arrsize2];
        std::fill_n(pOhmicCurr_DEP_Spec, arrsize1, 0);
        std::fill_n(pOhmicCurr_CHANSTATE, arrsize2, 0);

        // Fill the vectors with useful information
        for (uint ri = 0; ri < countOhmicCurrs(); ++ri)
        {
            OhmicCurrdef * ocdef = ohmiccurrdef(ri);
            for (uint si = 0; si < ngspecs; ++si)
            {
                if (ocdef->req(si) == false) { continue;
}
                // TODO: turn into error check?
                uint slidx = specG2L(si);
                AssertLog(slidx != LIDX_UNDEFINED);

                uint aridx = _IDX_OhmicCurr_Spec(ri, slidx);
                pOhmicCurr_DEP_Spec[aridx] = ocdef->dep(si);
            }
            pOhmicCurr_CHANSTATE[ri] = specG2L(ocdef->chanstate());
        }
    }

    // 6 -- DEAL WITH GHK CURRENTS
    if (pGHKcurrsN != 0)
    {
        // Set up local indices.
        pGHKcurr_L2G = new uint[countGHKcurrs()];
        uint ngghkcurrs = pStatedef->countGHKcurrs();
        for (uint i = 0; i < ngghkcurrs; ++i)
        {
            uint lidx = ghkcurrG2L(i);
            if (lidx == LIDX_UNDEFINED) { continue;
}
            pGHKcurr_L2G[lidx] = i;
        }
        // Create local _DEP and _CHANSTATE vectors.
        uint arrsize1 = countSpecs() * countGHKcurrs();
        uint arrsize2 = countGHKcurrs();
        pGHKcurr_DEP_Spec = new int[arrsize1];
        pGHKcurr_CHANSTATE = new uint[arrsize2];
        pGHKcurr_ION = new uint[arrsize2];
        std::fill_n(pGHKcurr_DEP_Spec, arrsize1, 0);
        std::fill_n(pGHKcurr_CHANSTATE, arrsize2, 0);
        std::fill_n(pGHKcurr_ION, arrsize2, 0);

        // Fill the vectors with useful information
        for (uint ri = 0; ri < countGHKcurrs(); ++ri)
        {
            GHKcurrdef * ghkdef = ghkcurrdef(ri);
            for (uint si = 0; si < ngspecs; ++si)
            {
                if (ghkdef->req(si) == false) { continue;
}
                // TODO: turn into error check?
                uint slidx = specG2L(si);
                // If not the volume ion species local index should be defined.
                if (ghkdef->req_v(si) == false) AssertLog(slidx != LIDX_UNDEFINED);

                uint aridx = _IDX_GHKcurr_Spec(ri, slidx);
                // DEP information can be rate (value:2) not just stoichiometry (value:1)
                pGHKcurr_DEP_Spec[aridx] = ghkdef->dep(si);
            }
            pGHKcurr_CHANSTATE[ri] = specG2L(ghkdef->chanstate());
            pGHKcurr_ION[ri] = specG2L(ghkdef->ion());
        }
    }

    // 7 -- DEAL WITH V-DEPENDENT TRANSITIONS
    if (pVDepTransN != 0)
    {
        // Set up local indices
        pVDepTrans_L2G = new uint[countVDepTrans()];
        uint ngvdeptrans = pStatedef->countVDepTrans();
        for (uint i = 0; i < ngvdeptrans; ++i)
        {
            uint lidx = vdeptransG2L(i);
            if (lidx == LIDX_UNDEFINED) { continue;
}
            pVDepTrans_L2G[lidx] = i;
        }
        // Create _DEP and _CHANSTATE vectors
        uint arrsize1 = countSpecs() * countVDepTrans();
        uint arrsize2 = countVDepTrans();
        pVDepTrans_DEP_Spec = new int[arrsize1];
        pVDepTrans_SRCCHANSTATE = new uint[arrsize2];
        pVDepTrans_DSTCHANSTATE = new uint[arrsize2];
        std::fill_n(pVDepTrans_DEP_Spec, arrsize1, 0);
        std::fill_n(pVDepTrans_SRCCHANSTATE, arrsize2, 0);
        std::fill_n(pVDepTrans_DSTCHANSTATE, arrsize2, 0);

        // Fill the vectors with all kinds of useful information.
        for (uint ri = 0; ri < countVDepTrans(); ++ri)
        {
            VDepTransdef * vdtdef = vdeptransdef(ri);
            AssertLog(vdtdef != 0);

            for (uint si = 0; si < ngspecs; ++si)
            {
                if (vdtdef->req(si) == false) { continue;
}

                // TODO: turn into error check?
                uint slidx = specG2L(si);
                AssertLog(slidx != LIDX_UNDEFINED);

                uint aridx = _IDX_VDepTrans_Spec(ri, slidx);
                pVDepTrans_DEP_Spec[aridx] = vdtdef->dep(si);
            }
            pVDepTrans_SRCCHANSTATE[ri] = specG2L(vdtdef->srcchanstate());
            pVDepTrans_DSTCHANSTATE[ri] = specG2L(vdtdef->dstchanstate());
        }
    }


    // Initialise the pools and flags members to zeros.
    if (pSpecsN_S != 0)
    {
        pPoolCount = new double[pSpecsN_S];
        pPoolFlags = new uint[pSpecsN_S];
        std::fill_n(pPoolCount, pSpecsN_S, 0.0);
        std::fill_n(pPoolFlags, pSpecsN_S, 0);
    }
    if (pSReacsN != 0)
    {
        pSReacFlags = new uint[pSReacsN];
        std::fill_n(pSReacFlags, pSReacsN, 0);

        // Finally initialise Kcsts to user-supplied values
        pSReacKcst = new double[pSReacsN];
        for (uint i = 0; i < pSReacsN; ++i)
        {
            // sreacdef() returns global Reacdef by local index
            ssolver::SReacdef * sreac = sreacdef(i);
            pSReacKcst[i] = sreac->kcst();
        }
    }

    if (pSurfDiffsN != 0)
    {
        pSurfDiffDcst = new double[pSurfDiffsN];

        for (uint i = 0; i <pSurfDiffsN; ++i)
        {
            // sdiffdef() returns global SDiffdef by local index
            ssolver::Diffdef * sdiff = surfdiffdef(i);
            pSurfDiffDcst[i] = sdiff->dcst();
        }
    }

    pSetupIndsdone = true;
}

////////////////////////////////////////////////////////////////////////////////

double ssolver::Patchdef::area() const
{
    return pArea;
}

////////////////////////////////////////////////////////////////////////////////

std::string const ssolver::Patchdef::name() const
{
    return pName;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::setArea(double a)
{
    AssertLog(a > 0.0);
    pArea = a;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::reset()
{
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    std::fill_n(pPoolCount, pSpecsN_S, 0.0);
    std::fill_n(pPoolFlags, pSpecsN_S, 0);
    std::fill_n(pSReacFlags, pSReacsN, 0);

    for (uint i = 0; i < pSReacsN; ++i)
    {
        // sreacdef() returns global Reacdef by local index
        ssolver::SReacdef * sreac = sreacdef(i);
        pSReacKcst[i] = sreac->kcst();
    }
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::setCount(uint slidx, double count)
{
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    AssertLog(slidx < pSpecsN_S);
    AssertLog(count >= 0.0);
    pPoolCount[slidx] = count;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::setClamped(uint slidx, bool clamp)
{
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    AssertLog(slidx < pSpecsN_S);
    if (clamp == true) { pPoolFlags[slidx] |= CLAMPED;
    } else { pPoolFlags[slidx] &= ~CLAMPED;
}
}

////////////////////////////////////////////////////////////////////////////////

ssolver::SReacdef * ssolver::Patchdef::sreacdef(uint lidx) const
{
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < pSReacsN);
    return pStatedef->sreacdef(pSReac_L2G[lidx]);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Diffdef * ssolver::Patchdef::surfdiffdef(uint dlidx) const
{
    AssertLog(pSetupRefsdone == true);
    AssertLog(dlidx < pSurfDiffsN);
    return pStatedef->surfdiffdef(pSurfDiff_L2G[dlidx]);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::VDepSReacdef * ssolver::Patchdef::vdepsreacdef(uint lidx) const
{
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < pVDepSReacsN);
    return pStatedef->vdepsreacdef(pVDepSReac_L2G[lidx]);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::OhmicCurrdef * ssolver::Patchdef::ohmiccurrdef(uint lidx) const
{
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < pOhmicCurrsN);
    return pStatedef->ohmiccurrdef(pOhmicCurr_L2G[lidx]);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::GHKcurrdef * ssolver::Patchdef::ghkcurrdef(uint lidx) const
{
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < pGHKcurrsN);
    return pStatedef->ghkcurrdef(pGHKcurr_L2G[lidx]);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::VDepTransdef * ssolver::Patchdef::vdeptransdef(uint lidx) const
{
    AssertLog(pSetupRefsdone == true);
    AssertLog(lidx < pVDepTransN);
    return pStatedef->vdeptransdef(pVDepTrans_L2G[lidx]);
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Patchdef::sreac_dep_I(uint srlidx, uint splidx) const
{
    return pSReac_DEP_I_Spec[splidx + (srlidx * pSpecsN_I)];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Patchdef::sreac_dep_S(uint srlidx, uint splidx) const
{
    return pSReac_DEP_S_Spec[splidx + (srlidx * pSpecsN_S)];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Patchdef::sreac_dep_O(uint srlidx, uint splidx) const
{
    return pSReac_DEP_O_Spec[splidx + (srlidx * pSpecsN_O)];
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::sreac_lhs_I_bgn(uint lidx) const
{
    return pSReac_LHS_I_Spec + (lidx * pSpecsN_I);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::sreac_lhs_I_end(uint lidx) const
{
    return pSReac_LHS_I_Spec + ((lidx + 1) * pSpecsN_I);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::sreac_lhs_S_bgn(uint lidx) const
{
    return pSReac_LHS_S_Spec + (lidx * pSpecsN_S);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::sreac_lhs_S_end(uint lidx) const
{
    return pSReac_LHS_S_Spec + ((lidx + 1) * pSpecsN_S);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::sreac_lhs_O_bgn(uint lidx) const
{
    return pSReac_LHS_O_Spec + (lidx * pSpecsN_O);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::sreac_lhs_O_end(uint lidx) const
{
    return pSReac_LHS_O_Spec + ((lidx + 1) * pSpecsN_O);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::sreac_upd_I_bgn(uint lidx) const
{
    return pSReac_UPD_I_Spec + (lidx * pSpecsN_I);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::sreac_upd_I_end(uint lidx) const
{
    return pSReac_UPD_I_Spec + ((lidx + 1) * pSpecsN_I);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::sreac_upd_S_bgn(uint lidx) const
{
    return pSReac_UPD_S_Spec + (lidx * pSpecsN_S);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::sreac_upd_S_end(uint lidx) const
{
    return pSReac_UPD_S_Spec + ((lidx + 1) * pSpecsN_S);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::sreac_upd_O_bgn(uint lidx) const
{
    return pSReac_UPD_O_Spec + (lidx * pSpecsN_O);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::sreac_upd_O_end(uint lidx) const
{
    return pSReac_UPD_O_Spec + ((lidx + 1) * pSpecsN_O);
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Patchdef::vdepsreac_dep_I(uint vdsrlidx, uint splidx) const
{
    return pVDepSReac_DEP_I_Spec[splidx + (vdsrlidx * pSpecsN_I)];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Patchdef::vdepsreac_dep_S(uint vdsrlidx, uint splidx) const
{
    return pVDepSReac_DEP_S_Spec[splidx + (vdsrlidx * pSpecsN_S)];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Patchdef::vdepsreac_dep_O(uint vdsrlidx, uint splidx) const
{
    return pVDepSReac_DEP_O_Spec[splidx + (vdsrlidx * pSpecsN_O)];
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::vdepsreac_lhs_I_bgn(uint lidx) const
{
    return pVDepSReac_LHS_I_Spec + (lidx * pSpecsN_I);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::vdepsreac_lhs_I_end(uint lidx) const
{
    return pVDepSReac_LHS_I_Spec + ((lidx + 1) * pSpecsN_I);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::vdepsreac_lhs_S_bgn(uint lidx) const
{
    return pVDepSReac_LHS_S_Spec + (lidx * pSpecsN_S);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::vdepsreac_lhs_S_end(uint lidx) const
{
    return pVDepSReac_LHS_S_Spec + ((lidx + 1) * pSpecsN_S);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::vdepsreac_lhs_O_bgn(uint lidx) const
{
    return pVDepSReac_LHS_O_Spec + (lidx * pSpecsN_O);
}

////////////////////////////////////////////////////////////////////////////////

uint * ssolver::Patchdef::vdepsreac_lhs_O_end(uint lidx) const
{
    return pVDepSReac_LHS_O_Spec + ((lidx + 1) * pSpecsN_O);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::vdepsreac_upd_I_bgn(uint lidx) const
{
    return pVDepSReac_UPD_I_Spec + (lidx * pSpecsN_I);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::vdepsreac_upd_I_end(uint lidx) const
{
    return pVDepSReac_UPD_I_Spec + ((lidx + 1) * pSpecsN_I);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::vdepsreac_upd_S_bgn(uint lidx) const
{
    return pVDepSReac_UPD_S_Spec + (lidx * pSpecsN_S);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::vdepsreac_upd_S_end(uint lidx) const
{
    return pVDepSReac_UPD_S_Spec + ((lidx + 1) * pSpecsN_S);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::vdepsreac_upd_O_bgn(uint lidx) const
{
    return pVDepSReac_UPD_O_Spec + (lidx * pSpecsN_O);
}

////////////////////////////////////////////////////////////////////////////////

int * ssolver::Patchdef::vdepsreac_upd_O_end(uint lidx) const
{
    return pVDepSReac_UPD_O_Spec + ((lidx + 1) * pSpecsN_O);
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::setKcst(uint srlidx, double kcst)
{
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    AssertLog(srlidx < pSReacsN);
    AssertLog(kcst >= 0.0);
    pSReacKcst[srlidx] = kcst;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Patchdef::setActive(uint srlidx, bool active)
{
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    AssertLog(srlidx < pSReacsN);
    if (active == true) { pSReacFlags[srlidx] &= ~INACTIVATED;
    } else { pSReacFlags[srlidx] |= INACTIVATED;
}
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Patchdef::specG2L_I(uint gidx) const
{
    if (pInner == nullptr) { return LIDX_UNDEFINED;
}
    return pInner->specG2L(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Patchdef::specG2L_O(uint gidx) const
{
    if (pOuter == nullptr) { return LIDX_UNDEFINED;
}
    return pOuter->specG2L(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Patchdef::ohmiccurr_dep_S(uint oclidx, uint splidx) const
{
    return pOhmicCurr_DEP_Spec[splidx + (oclidx*countOhmicCurrs())];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Patchdef::ohmiccurr_chanstate(uint oclidx) const
{
    return pOhmicCurr_CHANSTATE[oclidx];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Patchdef::ghkcurr_dep_S(uint ghklidx, uint splidx) const
{
    return pGHKcurr_DEP_Spec[splidx + (ghklidx*countGHKcurrs())];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Patchdef::ghkcurr_chanstate(uint ghklidx) const
{
    return pGHKcurr_CHANSTATE[ghklidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Patchdef::ghkcurr_ion(uint ghklidx) const
{
    return pGHKcurr_ION[ghklidx];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Patchdef::vdeptrans_dep_S(uint vdtlidx, uint splidx) const
{
    return pVDepTrans_DEP_Spec[splidx + (vdtlidx * countVDepTrans())];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Patchdef::vdeptrans_srcchanstate(uint vdtlidx) const
{
    return pVDepTrans_SRCCHANSTATE[vdtlidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Patchdef::vdeptrans_dstchanstate(uint vdtlidx) const
{
    return pVDepTrans_DSTCHANSTATE[vdtlidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Patchdef::surfdiff_dep(uint dlidx, uint slidx) const
{
    return pSurfDiff_DEP_Spec[slidx + (dlidx * pSpecsN_S)];
}
////////////////////////////////////////////////////////////////////////////////

// END
