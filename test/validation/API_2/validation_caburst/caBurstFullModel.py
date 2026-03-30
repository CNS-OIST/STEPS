####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2026 Okinawa Institute of Science and Technology, Japan.
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
###

########################################################################
import steps.interface

from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *

import numpy as np
import sys, os
from math import pow
import time

sys.path.append('..')
from . import constants_withampa_yunliang as par


class SimulationError(Exception):
    pass


def run(seed, mesh_path, steps_version, niter):
    seed, steps_version = int(seed), int(steps_version)

    if steps_version not in [3, 4]:
        raise SimulationError(f"Steps number: {steps_version} is not 3 or 4")

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    ######Glutamate transient#######
    # Reference (Rudolph et al. 2011)
    # Units (mM)
    with open("validation_caburst/Glut_Pulse_MVR.dat", "r") as f:
        Glut = list(map(float, f.read().split()))

    ########################### BIOCHEMICAL MODEL ###############################
    Vrange = (-150e-3, 100e-3, 0.01e-3)

    mdl = Model()
    r = ReactionManager()
    with mdl:
        # Species
        Ca = Species.Create(valence=2)
        Pump, CaPump, iCBsf, iCBsCa, iCBCaf, iCBCaCa, CBsf, CBsCa, CBCaf, CBCaCa, PV, PVMg, PVCa, Mg = Species.Create()

        # Channels
        CaP_m0, CaP_m1, CaP_m2, CaP_m3 = SubUnitState.Create()
        CaPchan = Channel.Create([CaP_m0, CaP_m1, CaP_m2, CaP_m3])

        BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4 = SubUnitState.Create()
        BKchan = Channel.Create([BK_C0, BK_C1, BK_C2, BK_C3, BK_C4, BK_O0, BK_O1, BK_O2, BK_O3, BK_O4])

        SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2 = SubUnitState.Create()
        SKchan = Channel.Create([SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2])

        AMPA_C, AMPA_C1, AMPA_C2, AMPA_D1, AMPA_D2, AMPA_O = SubUnitState.Create()
        AMPA = Channel.Create([AMPA_C, AMPA_C1, AMPA_C2, AMPA_D1, AMPA_D2, AMPA_O])

        Leak = SubUnitState.Create()
        L = Channel.Create([Leak])

        # Vol/surface systems
        vsys = VolumeSystem.Create()
        ssys = SurfaceSystem.Create()

        with vsys:
            diff_Ca = Diffusion.Create(Ca, par.DCST)
            diff_CBsf = Diffusion.Create(CBsf, par.DCB)
            diff_CBsCa = Diffusion.Create(CBsCa, par.DCB)
            diff_CBCaf = Diffusion.Create(CBCaf, par.DCB)
            diff_CBCaCa = Diffusion.Create(CBCaCa, par.DCB)
            diff_PV = Diffusion.Create(PV, par.DPV)
            diff_PVCa = Diffusion.Create(PVCa, par.DPV)
            diff_PVMg = Diffusion.Create(PVMg, par.DPV)

            # iCBsf fast and slow
            (iCBsf + Ca < r[1] > iCBsCa) + Ca < r[2] > iCBCaCa
            (iCBsf + Ca < r[3] > iCBCaf) + Ca < r[4] > iCBCaCa
            r[1].K = par.iCBsf1_f_kcst, par.iCBsf1_b_kcst
            r[2].K = par.iCBsCa_f_kcst, par.iCBsCa_b_kcst
            r[3].K = par.iCBsf2_f_kcst, par.iCBsf2_b_kcst
            r[4].K = par.iCBCaf_f_kcst, par.iCBCaf_b_kcst

            # CBsf fast and slow
            (CBsf + Ca < r[1] > CBsCa) + Ca < r[2] > CBCaCa
            (CBsf + Ca < r[3] > CBCaf) + Ca < r[4] > CBCaCa
            r[1].K = par.CBsf1_f_kcst, par.CBsf1_b_kcst
            r[2].K = par.CBsCa_f_kcst, par.CBsCa_b_kcst
            r[3].K = par.CBsf2_f_kcst, par.CBsf2_b_kcst
            r[4].K = par.CBCaf_f_kcst, par.CBCaf_b_kcst

            # PVCa
            PV + Ca < r[1] > PVCa
            r[1].K = par.PVca_f_kcst, par.PVca_b_kcst

            # PVMg
            PV + Mg < r[1] > PVMg
            r[1].K = par.PVmg_f_kcst, par.PVmg_b_kcst

        with ssys:
            # Ca Pump
            Pump.s + Ca.i < r[1] > CaPump.s > r[2] > Pump.s
            r[1].K = par.P_f_kcst, par.P_b_kcst
            r[2].K = par.P_k_kcst

            # CaP channel
            with CaPchan[...]:
                CaP_m0.s < r[1] > CaP_m1.s < r[2] > CaP_m2.s < r[3] > CaP_m3.s
                r[1].K = 3 * VDepRate(par.alpha_cap, vrange=Vrange), 1 * VDepRate(par.beta_cap, vrange=Vrange)
                r[2].K = 2 * VDepRate(par.alpha_cap, vrange=Vrange), 2 * VDepRate(par.beta_cap, vrange=Vrange)
                r[3].K = 1 * VDepRate(par.alpha_cap, vrange=Vrange), 3 * VDepRate(par.beta_cap, vrange=Vrange)
            OC_CaP = GHKCurr.Create(CaPchan[CaP_m3], Ca, par.CaP_P, computeflux=True, virtual_oconc=par.Ca_oconc)

            # BK channel
            with BKchan[...]:
                (((BK_C0.s + Ca.i < r[1] > BK_C1.s) \
                  + Ca.i < r[2] > BK_C2.s) \
                 + Ca.i < r[3] > BK_C3.s) \
                + Ca.i < r[4] > BK_C4.s
                r[1].K = par.c_01, par.c_10
                r[2].K = par.c_12, par.c_21
                r[3].K = par.c_23, par.c_32
                r[4].K = par.c_34, par.c_43

                (((BK_O0.s + Ca.i < r[1] > BK_O1.s) \
                  + Ca.i < r[2] > BK_O2.s) \
                 + Ca.i < r[3] > BK_O3.s) \
                + Ca.i < r[4] > BK_O4.s
                r[1].K = par.o_01, par.o_10
                r[2].K = par.o_12, par.o_21
                r[3].K = par.o_23, par.o_32
                r[4].K = par.o_34, par.o_43

                BK_C0.s < r[1] > BK_O0.s
                BK_C1.s < r[2] > BK_O1.s
                BK_C2.s < r[3] > BK_O2.s
                BK_C3.s < r[4] > BK_O3.s
                BK_C4.s < r[5] > BK_O4.s
                r[1].K = VDepRate(par.f_0, vrange=Vrange), VDepRate(par.b_0, vrange=Vrange)
                r[2].K = VDepRate(par.f_1, vrange=Vrange), VDepRate(par.b_1, vrange=Vrange)
                r[3].K = VDepRate(par.f_2, vrange=Vrange), VDepRate(par.b_2, vrange=Vrange)
                r[4].K = VDepRate(par.f_3, vrange=Vrange), VDepRate(par.b_3, vrange=Vrange)
                r[5].K = VDepRate(par.f_4, vrange=Vrange), VDepRate(par.b_4, vrange=Vrange)
            OC_BK0 = OhmicCurr.Create(BKchan[BK_O0], par.BK_G, par.BK_rev)
            OC_BK1 = OhmicCurr.Create(BKchan[BK_O1], par.BK_G, par.BK_rev)
            OC_BK2 = OhmicCurr.Create(BKchan[BK_O2], par.BK_G, par.BK_rev)
            OC_BK3 = OhmicCurr.Create(BKchan[BK_O3], par.BK_G, par.BK_rev)
            OC_BK4 = OhmicCurr.Create(BKchan[BK_O4], par.BK_G, par.BK_rev)

            # SK channel
            with SKchan[...]:
                ((SK_C1.s + Ca.i < r[1] > SK_C2.s) \
                 + Ca.i < r[2] > SK_C3.s) \
                + Ca.i < r[3] > SK_C4.s
                r[1].K = par.dirc2_t, par.invc1_t
                r[2].K = par.dirc3_t, par.invc2_t
                r[3].K = par.dirc4_t, par.invc3_t

                SK_C3.s < r[1] > SK_O1.s
                SK_C4.s < r[2] > SK_O2.s
                r[1].K = par.diro1_t, par.invo1_t
                r[2].K = par.diro2_t, par.invo2_t
            OC1_SK = OhmicCurr.Create(SKchan[SK_O1], par.SK_G, par.SK_rev)
            OC2_SK = OhmicCurr.Create(SKchan[SK_O2], par.SK_G, par.SK_rev)

            # AMPA channel
            with AMPA[...]:
                AMPA_C.s < r['AMPACC1'] > AMPA_C1.s < r['AMPAC1C2'] > AMPA_C2.s < r[3] > AMPA_O.s
                r['AMPACC1'].K = 0, par.ru1
                r['AMPAC1C2'].K = 0, par.ru2
                r[3].K = par.ro, par.rc

                AMPA_C1.s < r[1] > AMPA_D1.s
                AMPA_C2.s < r[2] > AMPA_D2.s
                r[1].K = par.rd, par.rr
                r[2].K = par.rd, par.rr
            OC_AMPAR1 = OhmicCurr.Create(AMPA[AMPA_O], par.AMPA_G, par.AMPA_rev)

            # Leak current channel
            OC_L = OhmicCurr.Create(L[Leak], par.L_G, par.L_rev)

    ########## MESH & COMPARTMENTALIZATION #################
    mesh = DistMesh(mesh_path, 1e-6) if steps_version == 4 else TetMesh.LoadGmsh(mesh_path, 1e-6)

    with mesh:
        z_cut = -4.0e-6 # pyramid mesh is just length 10 on z-axis -5 to 5, size tapers from 1 to 0.2

        smooth_tris = TriList(tri for tri in mesh.surface if tri.center.z <= z_cut)
        spiny_tris = TriList(tri for tri in mesh.surface if tri.center.z > z_cut)

        if steps_version == 4:
            cyto = Compartment.Create(vsys, tetLst=mesh.tets, conductivity=1 / par.Ra)
            smooth = Patch.Create(smooth_tris, cyto, None, ssys)
            spiny = Patch.Create(spiny_tris, cyto, None, ssys)
            memb_spiny = Membrane.Create([spiny], capacitance=par.memb_capac_spiny)
            memb_smooth = Membrane.Create([smooth], capacitance=par.memb_capac_proximal)
        else:
            cyto = Compartment.Create(mesh.tets, vsys)
            smooth = Patch.Create(smooth_tris, cyto, None, ssys)
            spiny = Patch.Create(spiny_tris, cyto, None, ssys)

            memb = Membrane.Create([spiny, smooth])
        
        record_points = [
            [0,0,-4.5e-6],
            [0,0,-1.5e-6],
            [0.1e-6,0,1.5e-6],
            [0,0,4.5e-6],
        ]
        
        record_tets = TetList(mesh.tets[point] for point in record_points)

    # # # # # # # # # # # # # # # # # # # # # # # # SIMULATION  # # # # # # # # # # # # # # # # # # # # # #

    rng = RNG('mt19937', 512, seed)

    if steps_version == 4:
        sim = Simulation('DistTetOpSplit', mdl, mesh, rng, searchMethod=NextEventSearchMethod.GIBSON_BRUCK)
    else:
        part = GmshPartition(mesh)
        sim = Simulation('TetOpSplit', mdl, mesh, rng, MPI.EF_DV_PETSC, part)

    # set temperature for ghk reactions
    sim.Temp = 273.15 + par.TEMPERATURE

    if MPI.rank == 0:
        print("solver created.", flush=True)

    rs = ResultSelector(sim)

    Pots = rs.TETS(record_tets).V

    sim.toSave(Pots, dt=par.TIMECONVERTER * 10)

    for i in range (0, niter):    
        sim.newRun()

        smooth_area = smooth.Area
        spiny_area = spiny.Area

        # Total pump is 1e-15 mol/cm2 ---> 1e-11 mol/m2
        # pumpnbs per unit area (im m2) is Total pump times AVOGADRO's NUMBER (1e-11 mol/m2 * 6.022e23 /mol )
        pumpnbs = 6.022141e12 * smooth_area

        sim.smooth.Pump.Count = round(pumpnbs)
        sim.smooth.CaPump.Count = 0

        sim.smooth.CaPchan[CaP_m0].Count = round(par.CaP_ro * smooth_area * par.CaP_m0_p)
        sim.smooth.CaPchan[CaP_m1].Count = round(par.CaP_ro * smooth_area * par.CaP_m1_p)
        sim.smooth.CaPchan[CaP_m2].Count = round(par.CaP_ro * smooth_area * par.CaP_m2_p)
        sim.smooth.CaPchan[CaP_m3].Count = round(par.CaP_ro * smooth_area * par.CaP_m3_p)

        sim.smooth.BKchan[BK_C0].Count = round(par.BK_ro * smooth_area * par.BK_C0_p)
        sim.smooth.BKchan[BK_C1].Count = round(par.BK_ro * smooth_area * par.BK_C1_p)
        sim.smooth.BKchan[BK_C2].Count = round(par.BK_ro * smooth_area * par.BK_C2_p)
        sim.smooth.BKchan[BK_C3].Count = round(par.BK_ro * smooth_area * par.BK_C3_p)
        sim.smooth.BKchan[BK_C4].Count = round(par.BK_ro * smooth_area * par.BK_C4_p)

        sim.smooth.BKchan[BK_O0].Count = round(par.BK_ro * smooth_area * par.BK_O0_p)
        sim.smooth.BKchan[BK_O1].Count = round(par.BK_ro * smooth_area * par.BK_O1_p)
        sim.smooth.BKchan[BK_O2].Count = round(par.BK_ro * smooth_area * par.BK_O2_p)
        sim.smooth.BKchan[BK_O3].Count = round(par.BK_ro * smooth_area * par.BK_O3_p)
        sim.smooth.BKchan[BK_O4].Count = round(par.BK_ro * smooth_area * par.BK_O4_p)

        sim.smooth.SKchan[SK_C1].Count = round(par.SK_ro * smooth_area * par.SK_C1_p)
        sim.smooth.SKchan[SK_C2].Count = round(par.SK_ro * smooth_area * par.SK_C2_p)
        sim.smooth.SKchan[SK_C3].Count = round(par.SK_ro * smooth_area * par.SK_C3_p)
        sim.smooth.SKchan[SK_C4].Count = round(par.SK_ro * smooth_area * par.SK_C4_p)

        sim.smooth.SKchan[SK_O1].Count = round(par.SK_ro * smooth_area * par.SK_O1_p)
        sim.smooth.SKchan[SK_O2].Count = round(par.SK_ro * smooth_area * par.SK_O2_p)

        sim.smooth.AMPA[AMPA_C].Count = round(par.AMPA_receptors)
        sim.smooth.AMPA[AMPA_C1].Count = 0
        sim.smooth.AMPA[AMPA_C2].Count = 0
        sim.smooth.AMPA[AMPA_O].Count = 0
        sim.smooth.AMPA[AMPA_D1].Count = 0
        sim.smooth.AMPA[AMPA_D2].Count = 0

        sim.smooth.L[Leak].Count = round(par.L_ro_proximal * smooth_area)
        

        # Total pump is 1e-15 mol/cm2 ---> 1e-11 mol/m2
        # pumpnbs per unit area (im m2) is Total pump times AVOGADRO's NUMBER (1e-11 mol/m2 * 6.022e23 /mol )
        pumpnbs = 6.022141e12 * spiny_area

        sim.spiny.Pump.Count = round(pumpnbs)
        sim.spiny.CaPump.Count = 0

        sim.spiny.CaPchan[CaP_m0].Count = round(par.CaP_ro * spiny_area * par.CaP_m0_p)
        sim.spiny.CaPchan[CaP_m1].Count = round(par.CaP_ro * spiny_area * par.CaP_m1_p)
        sim.spiny.CaPchan[CaP_m2].Count = round(par.CaP_ro * spiny_area * par.CaP_m2_p)
        sim.spiny.CaPchan[CaP_m3].Count = round(par.CaP_ro * spiny_area * par.CaP_m3_p)

        sim.spiny.BKchan[BK_C0].Count = round(par.BK_ro * spiny_area * par.BK_C0_p)
        sim.spiny.BKchan[BK_C1].Count = round(par.BK_ro * spiny_area * par.BK_C1_p)
        sim.spiny.BKchan[BK_C2].Count = round(par.BK_ro * spiny_area * par.BK_C2_p)
        sim.spiny.BKchan[BK_C3].Count = round(par.BK_ro * spiny_area * par.BK_C3_p)
        sim.spiny.BKchan[BK_C4].Count = round(par.BK_ro * spiny_area * par.BK_C4_p)

        sim.spiny.BKchan[BK_O0].Count = round(par.BK_ro * spiny_area * par.BK_O0_p)
        sim.spiny.BKchan[BK_O1].Count = round(par.BK_ro * spiny_area * par.BK_O1_p)
        sim.spiny.BKchan[BK_O2].Count = round(par.BK_ro * spiny_area * par.BK_O2_p)
        sim.spiny.BKchan[BK_O3].Count = round(par.BK_ro * spiny_area * par.BK_O3_p)
        sim.spiny.BKchan[BK_O4].Count = round(par.BK_ro * spiny_area * par.BK_O4_p)

        sim.spiny.SKchan[SK_C1].Count = round(par.SK_ro * spiny_area * par.SK_C1_p)
        sim.spiny.SKchan[SK_C2].Count = round(par.SK_ro * spiny_area * par.SK_C2_p)
        sim.spiny.SKchan[SK_C3].Count = round(par.SK_ro * spiny_area * par.SK_C3_p)
        sim.spiny.SKchan[SK_C4].Count = round(par.SK_ro * spiny_area * par.SK_C4_p)

        sim.spiny.SKchan[SK_O1].Count = round(par.SK_ro * spiny_area * par.SK_O1_p)
        sim.spiny.SKchan[SK_O2].Count = round(par.SK_ro * spiny_area * par.SK_O2_p)

        sim.spiny.AMPA[AMPA_C].Count = 0
        sim.spiny.AMPA[AMPA_C1].Count = 0
        sim.spiny.AMPA[AMPA_C2].Count = 0
        sim.spiny.AMPA[AMPA_O].Count = 0
        sim.spiny.AMPA[AMPA_D1].Count = 0
        sim.spiny.AMPA[AMPA_D2].Count = 0

        sim.spiny.L[Leak].Count = round(par.L_ro_spiny * spiny_area)


        sim.cyto.Ca.Conc = par.Ca_iconc
        
        sim.cyto.Mg.Conc = par.Mg_conc
        
        sim.cyto.iCBsf.Conc = par.iCBsf_conc
        sim.cyto.iCBCaf.Conc = par.iCBCaf_conc
        sim.cyto.iCBsCa.Conc = par.iCBsCa_conc
        sim.cyto.iCBCaCa.Conc = par.iCBCaCa_conc
        
        sim.cyto.CBsf.Conc = par.CBsf_conc
        sim.cyto.CBCaf.Conc = par.CBCaf_conc
        sim.cyto.CBsCa.Conc = par.CBsCa_conc
        sim.cyto.CBCaCa.Conc = par.CBCaCa_conc

        sim.cyto.PV.Conc = par.PV_conc
        sim.cyto.PVCa.Conc = par.PVCa_conc
        sim.cyto.PVMg.Conc = par.PVMg_conc


        if steps_version != 4:
            sim.ALL(Membrane).VolRes = par.Ra
            sim.TRIS(smooth.tris).Capac = par.memb_capac_proximal
            sim.TRIS(spiny.tris).Capac = par.memb_capac_spiny

            
        sim.EfieldDT = par.EF_DT

        sim.ALL(Membrane).Potential = par.init_pot


        for l in range(par.NTIMEPOINTS):

            sim.smooth.AMPACC1['fwd'].K = 1.0e-3 * par.rb * Glut[l + 2000]
            sim.smooth.AMPAC1C2['fwd'].K = 1.0e-3 * par.rb * Glut[l + 2000]

            sim.run(par.TIMECONVERTER * l)


    return Pots.data, Pots.time[0]


