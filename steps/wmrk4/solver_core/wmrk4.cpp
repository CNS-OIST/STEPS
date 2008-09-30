////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// This file is part of STEPS.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
//
// 
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/constants.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/reacdef.hpp>
#include <steps/sim/shared/patchdef.hpp>
#include <steps/sim/shared/sreacdef.hpp>
#include <steps/wmrk4/solver_core/comp.hpp>
#include <steps/wmrk4/solver_core/patch.hpp>
#include <steps/wmrk4/solver_core/wmrk4.hpp>


Wmrk4::Wmrk4(State * s)
: pState(s)
, pReacMtx(0)
, pUpdMtx(0)
, pSpecs_tot(0)
, pReacs_tot(0)
, pCcst()
, pVals()
, pSFlags()
, pRFlags()
, pNewVals()
, pDyDx()
, pDyDxlhs(0)
, pdt(0.0)										
, yt()
, dyt()
, dym()
{
	uint nspecstot = 0;
	uint nreacstot = 0;
	for(uint i=0; i< state()->countComps(); ++i)
	{
		nspecstot += state()->comp(i)->def()->countSpecs();
		nreacstot += state()->comp(i)->def()->countReacs();
	}
	for(uint i=0; i< state()->countPatches(); ++i)
	{
		nspecstot += state()->patch(i)->def()->countSpecs(); 
		nreacstot += state()->patch(i)->def()->countSReacs(); 
	}
	
	pReacMtx = new uint * [nreacstot];
	for (uint i=0; i< nreacstot; ++i)
	{
		pReacMtx[i]= new uint[nspecstot];
		std::fill_n(pReacMtx[i], nspecstot, 0);
	}
	
	pUpdMtx = new int * [nreacstot];
	for (uint i=0; i< nreacstot; ++i)
	{
		pUpdMtx[i] = new int[nspecstot];
		std::fill_n(pUpdMtx[i], nspecstot, 0);
	}
	
	pDyDxlhs = new double * [nspecstot];
	for (uint i=0; i< nspecstot; ++i)
	{
		pDyDxlhs[i] = new double[nreacstot];
		std::fill_n(pDyDxlhs[i], nreacstot, 0.0);
	}
	
}

///////////////////////////////////////////////////////////////////////////////

Wmrk4::~Wmrk4(void) 
{   
    for(uint i=0; i< pReacs_tot; ++i) delete[] pReacMtx[i];
	delete[] pReacMtx;
	for(uint i=0; i< pReacs_tot; ++i) delete[] pUpdMtx[i];
	delete[] pUpdMtx;
	for (uint i=0; i< pSpecs_tot; ++i) delete[] pDyDxlhs[i];
	delete[] pDyDxlhs;
}

///////////////////////////////////////////////////////////////////////////////

double Wmrk4::comp_ccst(double kcst, double vol, uint order)
{
    double vscale = 1.0e3 * vol * steps::math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;
    if (o1 < 0) o1 = 0;
    return kcst * pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

void Wmrk4::setup(void)

{   
	/// cumulative number of species from each compartment and patch
    pSpecs_tot = 0;	
	/// cumulative number of reacs from each compartment and patch
	pReacs_tot = 0;
	
	uint Comps_N = state()->countComps();
	uint Patches_N = state()->countPatches();
	
	for(uint i=0; i< Comps_N; ++i)
	{
		pSpecs_tot += state()->comp(i)->def()->countSpecs();
		pReacs_tot += state()->comp(i)->def()->countReacs();
	}
	for(uint i=0; i< Patches_N; ++i)
	{
		pSpecs_tot += state()->patch(i)->def()->countSpecs(); 
		pReacs_tot += state()->patch(i)->def()->countSReacs(); 
	}
	assert(pSpecs_tot > 0);
	assert(pReacs_tot > 0);
	
	/// initialise vectors 
	for(uint i=0; i< pReacs_tot; ++i) 
	{
		pRFlags.push_back(0);
	} 															
	for(uint i=0; i< pSpecs_tot; ++i)		
	{
		pVals.push_back(0.0);
		pSFlags.push_back(0);
		pNewVals.push_back(0.0);
		pDyDx.push_back(0.0);
		yt.push_back(0.0);
		dyt.push_back(0.0);
		dym.push_back(0.0);
	}
	
	/// fill the reaction matrix 
	/// loop over compartments, 
	/// then comp reacs and copy compdef LHS values to correct index
	
	/// set row marker to beginning of matrix for first compartment
	uint rowp = 0; 
	/// set column marker to beginning of matrix for first compartment
	uint colp = 0;
	
	for (uint i=0; i< Comps_N; ++i)  
	{   
		uint compReacs_N = state()->comp(i)->def()->countReacs();
		uint compSpecs_N = state()->comp(i)->def()->countSpecs();
		
        for(uint j=0; j< compReacs_N; ++j)
		{
			for(uint k=0; k< compSpecs_N; ++k)
			{
				uint lhs = state()->comp(i)->def()->reac_lhs_bgn(j)[k];
				int upd = state()->comp(i)->def()->reac_upd_bgn(j)[k];
				pReacMtx[rowp + j][colp + k] = lhs;
				pUpdMtx[rowp + j][colp + k] = upd;
			}
			/// set scaled reaction constant
			double reac_kcst = state()->comp(i)->def()->reac(j)->kcst();
			double comp_vol = state()->comp(i)->def()->vol();
			uint reac_order = state()->comp(i)->def()->reac(j)->order();
			pCcst.push_back(comp_ccst(reac_kcst, comp_vol, reac_order));
		}	
		/// step up markers for next compartment
		rowp += compReacs_N;
		colp += compSpecs_N;
	}
	
	/// now loop over patches, 
	/// then sreacs, filling for correct inner and outer compartments
	for(uint i=0; i< Patches_N; ++i)
	{
		uint patchReacs_N = state()->patch(i)->def()->countSReacs();
		uint patchSpecs_N_S = state()->patch(i)->def()->countSpecs();		
	    uint patchSpecs_N_I = state()->patch(i)->def()->countSpecs_I();		
		uint patchSpecs_N_O = state()->patch(i)->def()->countSpecs_O();      
		
		for (uint j=0; j< patchReacs_N; ++j)
		{
			for(uint k=0; k< patchSpecs_N_S; ++k)
			{   uint slhs = state()->patch(i)->def()->sreac_lhs_S_bgn(j)[k];
				int supd = state()->patch(i)->def()->sreac_upd_S_bgn(j)[k];
				pReacMtx[rowp + j][colp + k] = slhs;
				pUpdMtx[rowp + j][colp + k] = supd;
			}
			
			/// fill for inner and outer compartments involved in sreac j
			/// only perform if inner comp exists, similarly for outer comp
			if (state()->patch(i)->def()->icompdef() != 0) 
			{		
				/// fetch global index of inner compartment
				uint icompidx = state()->patch(i)->def()->icompdef()->gidx();
				// marker for correct position of inner compartment in matrix
				uint mtx_icompidx = 0;  
				/// step up marker to correct comp
				for (uint l=0; l< icompidx; ++l) 
				{
					mtx_icompidx += state()->comp(l)->def()->countSpecs();
				}
				for(uint k=0; k< patchSpecs_N_I; ++k)
				{
					uint ilhs = state()->patch(i)->def()->sreac_lhs_I_bgn(j)[k];
					int iupd = state()->patch(i)->def()->sreac_upd_I_bgn(j)[k];
					pReacMtx[rowp + j][mtx_icompidx + k] = ilhs;
					pUpdMtx[rowp + j][mtx_icompidx + k] = iupd;
				}
			}
			if (state()->patch(i)->def()->ocompdef() != 0)
			{
				uint ocompidx = state()->patch(i)->def()->ocompdef()->gidx();
				uint mtx_ocompidx =0;  
				for (uint l=0; l< ocompidx; ++l) 
				{
					mtx_ocompidx += state()->comp(l)->def()->countSpecs();
				}
				for(uint k=0; k< patchSpecs_N_O; ++k)
				{
					uint olhs = state()->patch(i)->def()->sreac_lhs_O_bgn(j)[k];
					int oupd = state()->patch(i)->def()->sreac_upd_O_bgn(j)[k];
					pReacMtx[rowp + j][mtx_ocompidx + k] = olhs;
					pUpdMtx[rowp + j][mtx_ocompidx + k] = oupd;
				}
			}
			/// set scaled reaction constant
			/// depends on volume of lhs reaction compartment
			double vol;
			if (state()->patch(i)->def()->sreac(j)->inside() == true)
			{
				assert(state()->patch(i)->iComp() != 0);
				vol = state()->patch(i)->iComp()->vol();
			}
			else
			{
				assert(state()->patch(i)->oComp() != 0);
				vol = state()->patch(i)->oComp()->vol();
			}
			double sreac_kcst = state()->patch(i)->def()->sreac(j)->kcst();
			uint sreac_order = state()->patch(i)->def()->sreac(j)->order();
			pCcst.push_back(comp_ccst(sreac_kcst, vol, sreac_order));
		} 
		/// move markers to next point in matrix
		rowp += patchReacs_N;									
		colp += patchSpecs_N_S;
	}
	
	assert (pCcst.size() == pReacs_tot); 
}

////////////////////////////////////////////////////////////////////////////////

void Wmrk4::refill(void)
{
	uint Comps_N = state()->countComps();
	uint Patches_N = state()->countPatches();
	assert (Comps_N > 0);
	
	uint c_marker = 0;
	uint r_marker = 0;
	
    for(uint i=0; i< Comps_N; ++i)
	{
		uint comp_Specs_N = state()->comp(i)->def()->countSpecs();
		uint comp_Reacs_N = state()->comp(i)->def()->countReacs();
		Comp * comp = state()->comp(i);
		assert(comp != 0);
		for (uint j=0; j < comp_Specs_N; ++j)
		{
			pVals[c_marker + j] = comp->pools()[j];
			pSFlags[c_marker + j] = state()->comp(i)->flags()[j];
		}
		for(uint k=0; k< comp_Reacs_N; ++k)
		{
			pRFlags[r_marker + k] = state()->comp(i)->rflags()[k];
		}
		c_marker += comp_Specs_N;
		r_marker += comp_Reacs_N;
	}
	
	for(uint i=0; i< Patches_N; ++i)
	{
		uint patch_Specs_N = state()->patch(i)->def()->countSpecs();
		uint patch_Reacs_N = state()->patch(i)->def()->countSReacs();
		for (uint j=0; j< patch_Specs_N; ++j)
		{
			pVals[c_marker +j] = state()->patch(i)->pools()[j];
			pSFlags[c_marker + j] = state()->patch(i)->flags()[j];
		}
		for (uint k=0; k< patch_Reacs_N; ++k)
		{
			pRFlags[r_marker + k] = state()->patch(i)->srflags()[k];
		}
		c_marker += patch_Specs_N;
		r_marker += patch_Reacs_N;
	}
	
	assert(c_marker == pVals.size());
	assert(pVals.size() == pSFlags.size());
	assert(r_marker == pRFlags.size());
	assert(pReacs_tot == pRFlags.size());
	assert(pSFlags.size() == pSpecs_tot);
}

////////////////////////////////////////////////////////////////////////////////

void Wmrk4::setderivs(dVec& vals, dVec& dydx)
{
	for (uint n=0; n< pSpecs_tot; ++n)
	{
		for(uint r=0; r< pReacs_tot; ++r)
		{	
			/// check reaction flags
			if (pRFlags[r] & State::INACTIVE_REACFLAG)						
			{
				pDyDxlhs[n][r] = 0.0;
				continue;
			}
			
			if (pUpdMtx[r][n] == 0)
			{
				pDyDxlhs[n][r] = 0.0;
				continue;
			}
			
			pDyDxlhs[n][r] = pUpdMtx[r][n] * pCcst[r];
			for (uint i=0; i < pSpecs_tot; ++i)
			{
				double dydx_lhs_temp = 1.0;
				double val = vals[i];
				switch (pReacMtx[r][i])
				{
					case 3: dydx_lhs_temp *= val;
					case 2: dydx_lhs_temp *= val;
					case 1: dydx_lhs_temp *= val;
					case 0: break;
					/// allow maximum 3 molecules of one species in reaction
					default: assert(0);
				}
				pDyDxlhs[n][r] *= dydx_lhs_temp;
			}
		}
	}
	
	for (uint n=0; n< pSpecs_tot; ++n)
	{
		dydx[n] = 0.0;
		for (uint r=0; r< pReacs_tot; ++r)
		{
			dydx[n] += pDyDxlhs[n][r];
		}
	}	
}

////////////////////////////////////////////////////////////////////////////////

void Wmrk4::rk4(void)
{
	double dt_2 = pdt/2.0;
	double dt_6 = pdt/6.0;
	
	for(uint i=0; i< pSpecs_tot; ++i) yt[i] = pVals[i] + (dt_2 * pDyDx[i]);
	setderivs(yt, dyt);
	for(uint i =0; i< pSpecs_tot; ++i) yt[i]= pVals[i] + (dt_2 * dyt[i]);
	setderivs(yt, dym);
	for (uint i=0; i< pSpecs_tot; ++i)
	{
		yt[i] = pVals[i] + (pdt * dym[i]);
		dym[i] += dyt[i];
	}
	setderivs(yt, dyt);
	for (uint i=0; i< pSpecs_tot; ++i) 
	{
		pNewVals[i]=pVals[i]+dt_6*(pDyDx[i]+dyt[i]+(2.0*dym[i]));
	}
}

////////////////////////////////////////////////////////////////////////////////

void Wmrk4::rksteps(double t1, double t2)
{
	assert(t1 < t2);
	double t = t1;
	pdt = state()->getDT();			
	assert (pdt > 0.0);
	assert (pdt <= (t2-t1));
	
	/// step up until over maximum time
	while(t < t2)
	{
		setderivs(pVals, pDyDx);
		rk4();						
		update();
		t += pdt;
	}
}

////////////////////////////////////////////////////////////////////////////////

void Wmrk4::update(void)
{
	/// update local values vector with computed counts
	for (uint i=0; i< pSpecs_tot; ++i)
	{
		/// check clamped flag and only update if not clamped
		if (pSFlags[i] & State::CLAMPED_POOLFLAG)
		{
			continue;
		}
		else
		{
			pVals[i] = pNewVals[i];
		}
	}
	
	/// update pools with computed values
	uint Comps_N = state()->countComps();
	uint Patches_N = state()->countPatches();
	uint c_marker = 0;
	
	for(uint i=0; i< Comps_N; ++i)
	{
		uint comp_Specs_N = state()->comp(i)->def()->countSpecs();
		for (uint j=0; j< comp_Specs_N; ++j)
		{	
			double c = pVals[c_marker + j];
		    state()->comp(i)->pools()[j] = c;
		}
		c_marker += comp_Specs_N;
	}
	
	for(uint i=0; i< Patches_N; ++i)
	{
		uint patch_Specs_N = state()->patch(i)->def()->countSpecs();
		for (uint j=0; j< patch_Specs_N; ++j)
		{
			double c = pVals[c_marker + j];
			state()->patch(i)->pools()[j] = c;
		}
		c_marker += patch_Specs_N;
	}
}

////////////////////////////////////////////////////////////////////////////////

// END
