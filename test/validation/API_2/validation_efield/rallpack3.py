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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import math
import pathlib
import sys

import numpy as np
import steps.interface
from steps.saving import *
from steps.sim import *
from steps.utils import *

from steps.geom import *
from steps.model import *
from steps.rng import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

sim_parameters = {
    ## Rallpack3:
    "K_G": 360,  # Potassium conductance, Siemens/m^2
    "Na_G": 1200,  # Sodium conductance, Siemens/m^2
    "L_G": 0.25,  # Leak conductance, Siemens/m^2
    "K_rev": -77e-3,  # Potassium reversal potential, V
    "Na_rev": 50e-3,  # Sodium reversal potential, V
    "leak_rev": -65.0e-3,  # Leak reveral potential, V
    "K_ro": 18.0e12
    * 5,  # this should avoid the missing spike issue# Potassium channel density
    "Na_ro": 60.0e12
    * 5,  # this should avoid the missing spike issue# Sodium channel density
    "Ra": 1.0,  # Ohm.m
    "C_M": 0.01,  # F.m^-2
    "Iinj": 0.1e-9,  # The current injection in amps
    "K_FACS": [
        0.216750577045,
        0.40366011853,
        0.281904943772,
        0.0874997924409,
        0.0101845682113,
    ],  # A table of potassium density factors at -65mV, found in getpops. n0, n1, n2, n3, n4
    "NA_FACS": [
        0.343079175644,
        0.0575250437508,
        0.00321512825945,
        5.98988373918e-05,
        0.506380603793,
        0.0849062503811,
        0.00474548939393,
        8.84099403236e-05,
    ],  # A table of sodium density factors. m0h1, m1h1, m2h1, m3h1, m0h0, m1h0, m2h0, m3h0
    "surfarea_cyl": 1.0 * math.pi * 1000 * 1e-12,
    "vol_cyl": math.pi * 0.5 * 0.5 * 1000 * 1e-18,
    "K_i_conc": 155e-3,
    "K_o_conc": 4e-3,
    "Na_i_conc": 50e-3,
    "Na_o_conc": 440e-3,
    "K_DCST": 1e-9,
    "Na_DCST": 1e-9,
    # STEPS
    "SIM_END": 0.01,  # simulation stop time s
    "EF_DT": 1e-5,  # E-field evaluation time step s
    "SAVE_DT": 0.000005,  # smaller dt to check how the ks test deals with discretizations
    "solver": None,
    "currType": "Ohmic",
    "ion_diff": False,
}


def build_geometry(mesh_path, param):
    axonTag = "Axon"
    ecsTag = "ECS"
    if param["solver"].endswith("DistTetOpSplit"):
        mesh = DistMesh(mesh_path, 1e-6)
        axon_tets = mesh.tetGroups[axonTag]
        ecs_tets = mesh.tetGroups[ecsTag]
    else:
        mesh = TetMesh.LoadGmsh(mesh_path, 1e-6)
        axon_tets = mesh.tetGroups[(0, axonTag)]
        ecs_tets = mesh.tetGroups[(0, ecsTag)]

    with mesh:

        axon = Compartment.Create(axon_tets, "axonvsys")
        ecs = Compartment.Create(ecs_tets, "ecsvsys")

        zmin_tris = TriList(
            [tri for tri in axon.surface if tri.center.z == mesh.bbox.min.z]
        )
        zmax_tris = TriList(
            [tri for tri in axon.surface if tri.center.z == mesh.bbox.max.z]
        )
        memb = Patch.Create(
            axon.surface & ecs.surface - zmax_tris - zmin_tris, axon, ecs, "ssys"
        )
        z_min = Patch.Create(zmin_tris, axon, None)
        z_max = Patch.Create(zmax_tris, axon, None)

        corr_fac_vol = axon.Vol / sim_parameters["vol_cyl"]
        corr_fac_area = memb.Area / sim_parameters["surfarea_cyl"]

        if param["solver"].endswith("DistTetOpSplit"):
            axon.Conductivity = 1 / (sim_parameters["Ra"] * corr_fac_vol)
            membrane = Membrane.Create(
                [memb], capacitance=sim_parameters["C_M"] / corr_fac_area
            )
        elif param["solver"] == "Tetexact":
            membrane = Membrane.Create([memb, z_min, z_max])
        else:
            membrane = Membrane.Create([memb])

    return mesh


def build_model(mesh, param):
    mdl = Model()
    r = ReactionManager()
    with mdl:
        axonvsys, ecsvsys = VolumeSystem.Create()
        ssys = SurfaceSystem.Create()

        # K channel
        n0, n1, n2, n3, n4 = SubUnitState.Create()
        VGKC = Channel.Create([n0, n1, n2, n3, n4])

        # Na channel
        h0, h1, m0, m1, m2, m3 = SubUnitState.Create()
        mNaSU, hNaSU = SubUnit.Create(
            [m0, m1, m2, m3],
            [h0, h1],
        )
        VGNaC = Channel.Create([mNaSU, hNaSU])

        # Leak
        leaksus = SubUnitState.Create()
        Leak = Channel.Create([leaksus])

        if param["currType"] == "GHK":
            K_ion = Species.Create(valence=1)
            Na_ion = Species.Create(valence=1)

        # Gating kinetics
        a_n = VDepRate(
            lambda V: 1e3
            * (
                (
                    0.01
                    * (10 - (V * 1e3 + 65))
                    / (math.exp((10 - (V * 1e3 + 65)) / 10) - 1)
                )
            )
        )
        b_n = VDepRate(lambda V: 1e3 * ((0.125 * math.exp(-(V * 1e3 + 65) / 80))))
        a_m = VDepRate(
            lambda V: 1e3
            * (
                (
                    0.1
                    * (25 - (V * 1e3 + 65))
                    / (math.exp((25 - (V * 1e3 + 65)) / 10) - 1)
                )
            )
        )
        b_m = VDepRate(lambda V: 1e3 * ((4 * math.exp(-(V * 1e3 + 65) / 18))))
        a_h = VDepRate(lambda V: 1e3 * ((0.07 * math.exp(-(V * 1e3 + 65) / 20))))
        b_h = VDepRate(
            lambda V: 1e3 * ((1 / (math.exp((30 - (V * 1e3 + 65)) / 10) + 1)))
        )

        with ssys:
            with VGKC[...]:
                n0.s < r[1] > n1.s < r[2] > n2.s < r[3] > n3.s < r[4] > n4.s
                r[1].K = 4 * a_n, b_n
                r[2].K = 3 * a_n, 2 * b_n
                r[3].K = 2 * a_n, 3 * b_n
                r[4].K = a_n, 4 * b_n

            with VGNaC[...]:
                h1.s < r[1] > h0.s
                r[1].K = a_h, b_h

                m0.s < r[1] > m1.s < r[2] > m2.s < r[3] > m3.s
                r[1].K = 3 * a_m, b_m
                r[2].K = 2 * a_m, 2 * b_m
                r[3].K = a_m, 3 * b_m

            if param["currType"] == "Ohmic":
                VGKC_I = OhmicCurr.Create(
                    VGKC[n4], param["K_G"] / param["K_ro"], param["K_rev"]
                )
                VGNaC_I = OhmicCurr.Create(
                    VGNaC[m3, h0], param["Na_G"] / param["Na_ro"], param["Na_rev"]
                )
            elif param["currType"] == "GHK":
                K_Pinfo = GHKCurr.PInfo(
                    g=param["K_G"] / param["K_ro"],
                    V=param["K_rev"],
                    T=293.15,
                    oconc=param["K_o_conc"],
                    iconc=param["K_i_conc"],
                )
                VGKC_I = GHKCurr.Create(
                    VGKC[n4],
                    K_ion,
                    K_Pinfo,
                    computeflux=True,
                )
                Na_Pinfo = GHKCurr.PInfo(
                    g=param["Na_G"] / param["Na_ro"],
                    V=param["Na_rev"],
                    T=293.15,
                    oconc=param["Na_o_conc"],
                    iconc=param["Na_i_conc"],
                )
                VGNaC_I = GHKCurr.Create(
                    VGNaC[m3, h0],
                    Na_ion,
                    Na_Pinfo,
                    computeflux=True,
                )

            else:
                raise NotImplementedError()

            # Set the single-channel conductance:
            L_G_tot = param["L_G"] * param["surfarea_cyl"]
            g_leak_sc = L_G_tot / len(mesh.memb.tris)
            OC_L = OhmicCurr.Create(Leak[leaksus], g_leak_sc, param["leak_rev"])

        if param["ion_diff"]:
            with axonvsys:
                Diffusion(K_ion, param["K_DCST"])
                Diffusion(Na_ion, param["Na_DCST"])
            with ecsvsys:
                Diffusion(K_ion, param["K_DCST"])
                Diffusion(Na_ion, param["Na_DCST"])

    return mdl


def create_sim(model, mesh, param):
    rng = RNG("mt19937", 512, param["seed"])

    if param["solver"] == "Tetexact":
        sim = Simulation("Tetexact", model, mesh, rng, calcMembPot=MPI.EF_DV_BDSYS)
    elif param["solver"] == "TetODE":
        sim = Simulation("TetODE", model, mesh, rng, calcMembPot=MPI.EF_DV_BDSYS)
    elif param["solver"] == "TetOpSplit":
        part = LinearMeshPartition(mesh, 1, 1, MPI.nhosts)
        sim = Simulation("TetOpSplit", model, mesh, rng, MPI.EF_DV_PETSC, part)
    elif param["solver"] == "TetVesicle":
        sim = Simulation("TetVesicle", model, mesh, rng, MPI.EF_DV_PETSC)
    elif param["solver"] == "DistTetOpSplit":
        sim = Simulation("DistTetOpSplit", model, mesh, rng, isEfield=True)
    elif param["solver"] == "RLeapDistTetOpSplit":
        sim = Simulation(
            "DistTetOpSplit",
            model,
            mesh,
            rng,
            isEfield=True,
            SSAMethod=SSAMethod.RLEAPING,
            searchMethod=NextEventSearchMethod.RLEAPING,
        )
    else:
        raise ValueError(f"SSA solver {param['solver']} not available")
    if param["solver"] == "TetODE":
        sim.setTolerances(1e-6, 1e-6)
    else:
        sim.EfieldDT = param["EF_DT"]

    return sim


def init_sim(sim, param):
    mdl, mesh = sim.model, sim.geom

    sim.TRIS(sim.geom.memb.tris).Leak[mdl.leaksus].Count = 1

    Na_ro, K_ro, NA_FACS, K_FACS = (
        param["Na_ro"],
        param["K_ro"],
        param["NA_FACS"],
        param["K_FACS"],
    )
    surfarea_cyl = param["surfarea_cyl"]

    sim.memb.VGNaC[mdl.m0, mdl.h1].Count = Na_ro * surfarea_cyl * NA_FACS[0]
    sim.memb.VGNaC[mdl.m1, mdl.h1].Count = Na_ro * surfarea_cyl * NA_FACS[1]
    sim.memb.VGNaC[mdl.m2, mdl.h1].Count = Na_ro * surfarea_cyl * NA_FACS[2]
    sim.memb.VGNaC[mdl.m3, mdl.h1].Count = Na_ro * surfarea_cyl * NA_FACS[3]
    sim.memb.VGNaC[mdl.m0, mdl.h0].Count = Na_ro * surfarea_cyl * NA_FACS[4]
    sim.memb.VGNaC[mdl.m1, mdl.h0].Count = Na_ro * surfarea_cyl * NA_FACS[5]
    sim.memb.VGNaC[mdl.m2, mdl.h0].Count = Na_ro * surfarea_cyl * NA_FACS[6]
    sim.memb.VGNaC[mdl.m3, mdl.h0].Count = Na_ro * surfarea_cyl * NA_FACS[7]
    sim.memb.VGKC[mdl.n0].Count = K_ro * surfarea_cyl * K_FACS[0]
    sim.memb.VGKC[mdl.n1].Count = K_ro * surfarea_cyl * K_FACS[1]
    sim.memb.VGKC[mdl.n2].Count = K_ro * surfarea_cyl * K_FACS[2]
    sim.memb.VGKC[mdl.n3].Count = K_ro * surfarea_cyl * K_FACS[3]
    sim.memb.VGKC[mdl.n4].Count = K_ro * surfarea_cyl * K_FACS[4]

    sim.membrane.Potential = -65e-3

    minzverts = mesh.z_min.tris.verts
    sim.VERTS(minzverts).IClamp = param["Iinj"] / len(minzverts)

    if not param["solver"].endswith("DistTetOpSplit"):
        corr_fac_vol = mesh.axon.Vol / param["vol_cyl"]
        corr_fac_area = mesh.memb.Area / param["surfarea_cyl"]
        sim.membrane.Capac = param["C_M"] / corr_fac_area
        sim.membrane.VolRes = param["Ra"] * corr_fac_vol

    if param["currType"] == "GHK":
        sim.axon.K_ion.Conc = param["K_i_conc"]
        sim.axon.Na_ion.Conc = param["Na_i_conc"]
        sim.ecs.K_ion.Conc = param["K_o_conc"]
        sim.ecs.Na_ion.Conc = param["Na_o_conc"]


def run_comparison(mesh_file, params):
    params = {**sim_parameters, **params}

    mesh = build_geometry(mesh_file, params)
    model = build_model(mesh, params)

    sim = create_sim(model, mesh, params)

    rs = ResultSelector(sim)

    # record potential at the two extremes along (z) axis
    Vzmin = rs.SUM(rs.VERTS(mesh.z_min.tris.verts).V) / len(mesh.z_min.tris.verts)
    Vzmax = rs.SUM(rs.VERTS(mesh.z_max.tris.verts).V) / len(mesh.z_max.tris.verts)

    sim.toSave(Vzmin, Vzmax, dt=params["SAVE_DT"])

    if params["saving_mode"] == "xdmf":
        pots = rs.TETS(mesh.tets).V
        cyto_specs = rs.TETS(mesh.tets).ALL(Species).Conc
        sim.toSave(pots, cyto_specs, dt=params["SAVE_DT"])

        hdf = XDMFHandler("Rallpack3", mode="w")
        sim.toDB(hdf, "Rallpack3")

    sim.newRun()

    init_sim(sim, params)

    sim.run(params["SIM_END"])

    if params["saving_mode"] == "xdmf":
        del hdf

    # Reference data
    v0path = pathlib.Path(params["datadir"]).joinpath(params["v0data"])
    v1path = pathlib.Path(params["datadir"]).joinpath(params["v1data"])
    with v0path.open("r") as f:
        V0ref = np.array(
            [list(map(float, line.split())) for line in f.read().split("\n")[2:-1]]
        )
        V0ref = np.interp(Vzmin.time[0] * 1e3, V0ref[:, 0], V0ref[:, 1])
    with v1path.open("r") as f:
        V1ref = np.array(
            [list(map(float, line.split())) for line in f.read().split("\n")[2:-1]]
        )
        V1ref = np.interp(Vzmin.time[0] * 1e3, V1ref[:, 0], V1ref[:, 1])

    return (
        Vzmin.time[...] * 1e3,
        Vzmin.data[...] * 1e3,
        Vzmax.data[...] * 1e3,
        V0ref,
        V1ref,
    )
