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
#
# Rallpack1 model
# Author Iain Hepburn

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.utils import *

import numpy as np
import numpy.linalg as la
import sys

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

sim_parameters = {
## Rallpack1:
    'R_A'  :    1.0,          # axial resistivity Ω·m
    'R_M'  :    4.0,          # membrane resistivity Ω·m²
    'C_M'  :    0.01,         # membrane capacity F/m²
    'E_M'  :   -0.065,        # p.d. across membrane V
    'Iinj' :    0.1e-9,       # injection current A
    'diameter': 1.0e-6,       # cylinder diameter m
    'length':   1.0e-3,       # cylinder length m
# STEPS
    'sim_end':  0.25,         # simulation stop time s
    'EF_dt':    1.0e-5,        # E-field evaluation time step s
    'solver': None,
    'currType': 'Ohmic',
}

# Build steps geom object, and report lower and upper surface vertixes
# for later use.

def build_geometry(mesh_path, param):
    if param['solver'] == 'DistTetOpSplit':
        mesh = DistMesh(mesh_path)
    else:
        mesh = TetMesh.Load(mesh_path)

    with mesh:

        zmin_tris = TriList(tri for tri in mesh.surface if all(v.z <= mesh.bbox.min.z for v in tri.verts))
        zmax_tris = TriList(tri for tri in mesh.surface if all(v.z >= mesh.bbox.max.z for v in tri.verts))

        if param['comp_mode'] == 'single':
            cyto = Compartment.Create(mesh.tets)

            memb_tris = mesh.surface - (zmin_tris | zmax_tris)
            memb = Patch.Create(memb_tris, cyto, None, 'ssys')
            membrane = Membrane.Create([memb])
        elif param['comp_mode'].startswith('multi'):
            # Multiple compartments with different conductivity
            fact = 0.5
            if param['comp_mode'].endswith('seq'):
                # Compartments are sequential (one after the other on the z axis)
                # So the voltage at the end should be equivalent when the resistivity of each compartment is scaled such that:
                # 2R =        R       +        R        # uniform resistivity case
                #    = R * (1 + fact) + R * (1 - fact)  # non-uniform resistivity case
                #
                # We should only check that the end point voltage is as expected, not the start
                tets = TetList([tet for tet in mesh.tets if tet.center.z < mesh.bbox.center.z])
                cyto1 = Compartment.Create(tets, conductivity=1 / ((1 + fact) * param['R_A']))
                cyto2 = Compartment.Create(mesh.tets - tets, conductivity=1 / ((1 - fact) * param['R_A']))
            elif param['comp_mode'].endswith('par'):
                # Compartments are parallel (one on top of the other on the z axis)
                # So the voltage should be equivalent when the conductivity of each compartment is scaled such that:
                # 2G =        G       +        G        # uniform conductivity case
                #    = G * (1 + fact) + G * (1 - fact)  # non-uniform conductivity case
                # 2R = 1 / ((1 + fact) / R + (1 - fact) / R)
                #
                # We can check both start and end point voltages
                tets = TetList([tet for tet in mesh.tets if tet.center.x < mesh.bbox.center.x])
                cyto1 = Compartment.Create(tets, conductivity=(1 + fact) / param['R_A'])
                cyto2 = Compartment.Create(mesh.tets - tets, conductivity=(1 - fact) / param['R_A'])
            elif param['comp_mode'].endswith('rand'):
                # Randomly distribute tetrahedrons to compartments
                # It is unclear how the conductivity should be set in that case.
                # So we only check that it does not blow up.
                tets = TetList([tet for tet in mesh.tets if hash(tet.idx) % 2 == 0])
                cyto1 = Compartment.Create(tets, conductivity=(1 + fact) / param['R_A'])
                cyto2 = Compartment.Create(mesh.tets - tets, conductivity=(1 - fact) / param['R_A'])
            elif param['comp_mode'].endswith('maxcomps'):
                # Create as many compartments as possible to be part of the conduction volume
                memb_patches = []
                tet2comp = {}
                *tris, last_tri = mesh.surface - (zmin_tris | zmax_tris)
                for tri in tris:
                    tet = tri.tetNeighbs[0]
                    if tet not in tet2comp:
                        comp = Compartment(tet.toList(), conductivity= 1 / param['R_A'])
                        tet2comp[tet] = comp
                    memb_patches.append(Patch(tri.toList(), tet2comp[tet], None, 'ssys'))
                last_comp = Compartment(mesh.tets - TetList(tet2comp.keys()), conductivity= 1 / param['R_A'])
                memb_patches.append(Patch(last_tri.toList(), last_comp, None, 'ssys'))
            else:
                raise NotImplementedError()

            if param['comp_mode'].endswith('maxcomps'):
                membrane = Membrane.Create(memb_patches)
            else:
                memb1 = Patch.Create((mesh.surface & cyto1.surface) - (zmin_tris | zmax_tris), cyto1, None, 'ssys')
                memb2 = Patch.Create((mesh.surface & cyto2.surface) - (zmin_tris | zmax_tris), cyto2, None, 'ssys')

                membrane = Membrane.Create([memb1, memb2])
        else:
            raise NotImplementedError()

        v_zmin = zmin_tris.verts

        # Find the Z-axially closest and furthest vertices
        distFunc = lambda v: la.norm(v[:2])
        v_zmin_sample = VertList([min(zmin_tris.verts, key=distFunc), max(zmin_tris.verts, key=distFunc)])
        v_zmax_sample = VertList([min(zmax_tris.verts, key=distFunc), max(zmax_tris.verts, key=distFunc)])

    return mesh, v_zmin, v_zmin_sample, v_zmax_sample


def build_model(mesh, param):
    mdl = Model()
    r = ReactionManager()
    with mdl:
        ssys = SurfaceSystem.Create()

        if param['currType'] == 'Ohmic':
            Leak = SubUnitState.Create()
            L = Channel.Create([Leak])
        elif param['currType'] == 'SReacCharge':
            Leak = Species.Create()
        else:
            raise NotImplementedError()

        # membrane conductance
        area_cylinder = np.pi * param['diameter'] * param['length']
        L_G_tot = area_cylinder / param['R_M']
        g_leak_sc = L_G_tot / sum(len(memb.tris) for memb in mesh.ALL(Patch))
        with ssys:
            if param['currType'] == 'Ohmic':
                OC_L = OhmicCurr.Create(L[Leak], g_leak_sc, param['E_M'])
            elif param['currType'] == 'SReacCharge':
                Leak.s > r[1] > Leak.s
                r[1].K = VDepRate(lambda V: (1.0/1.60217663e-19) * g_leak_sc*max(0, V-param['E_M']), vrange=(-150e-3, 150e-3, 0.01e-3))
                r[1].Charge = -1
            else:
                raise NotImplementedError()

    return mdl


def create_sim(model, mesh, param):
    # Create the solver objects
    if param['solver'] == 'TetODE':
        sim = Simulation('TetODE', model, mesh, calcMembPot=MPI.EF_DV_BDSYS)
        sim.setTolerances(1.0e-6, 1e-6)
    else:
        rng = RNG('mt19937', 512, param['seed'])
        if param['solver'] == 'Tetexact':
            sim = Simulation('Tetexact', model, mesh, rng, calcMembPot=MPI.EF_DV_BDSYS)
        elif param['solver'] == 'TetOpSplit':
            part = LinearMeshPartition(mesh, 1, 1, MPI.nhosts)
            sim = Simulation('TetOpSplit', model, mesh, rng, MPI.EF_DV_PETSC, part)
        elif param['solver'] == 'TetVesicle':
            sim = Simulation('TetVesicle', model, mesh, rng, MPI.EF_DV_PETSC)
        elif param['solver'] == 'DistTetOpSplit':
            sim = Simulation('DistTetOpSplit', model, mesh, rng, isEfield=True)
        else :
            raise ValueError(f"SSA solver {param['solver']} not available")
        sim.EfieldDT = param['EF_dt']
        
    print(f"Running Rallpack1 test with {param['solver']}")

    return sim


def init_sim(sim, model, mesh, v_zmin, param):
    SetVerbosity(0)

    # Set initial conditions

    for memb in mesh.ALL(Patch):
        if param['currType'] == 'Ohmic':
            sim.TRIS(memb.tris).L[model.Leak].Count = 1
        elif param['currType'] == 'SReacCharge':
            sim.TRIS(memb.tris).Leak.Count = 1
        else:
            raise NotImplementedError()

    sim.ALL(Membrane).Potential = param['E_M']
    # Correction factor for deviation between mesh and model cylinder:
    area_cylinder = np.pi * param['diameter'] * param['length']
    area_mesh_factor = sum(memb.Area for memb in mesh.ALL(Patch)) / area_cylinder
    if param['comp_mode'] == 'single':
        sim.ALL(Membrane).VolRes = param['R_A']
    sim.ALL(Membrane).Capac = param['C_M']/area_mesh_factor

    sim.VERTS(v_zmin).IClamp = param['Iinj']/len(v_zmin)


# Returns RMS error, table containing computed end-point voltages
# and reference voltage data.

def run_comparison(mesh_file, v0_datafile, v1_datafile, params, verbose=False):
    params = {**sim_parameters, **params}

    # sample at same interval as rallpack1 reference data
    sim_dt = 5.0e-5

    def snarf(fname):
        with open(fname, 'r') as f:
            return [tuple(map(float, line.split())) for line in f]

    vref_0um = np.array([v for (t,v) in snarf(v0_datafile)])
    vref_1000um = np.array([v for (t,v) in snarf(v1_datafile)])

    geom, v_zmin, zmin_sample, zmax_sample = build_geometry(mesh_file, params)
    model = build_model(geom, params)

    # grab sample vertices
    n_zmin_sample = len(zmin_sample)
    vertices =  zmin_sample + zmax_sample

    sim = create_sim(model, geom, params)

    rs = ResultSelector(sim)

    result = rs.VERTS(vertices).V

    sim.toSave(result, dt=sim_dt)

    if params['saving_mode'] == 'xdmf':
        pots = rs.TETS(sim.geom.tets).V
        sim.toSave(pots, dt=sim_dt)

        hdf = XDMFHandler('Rallpack1')
        sim.toDB(hdf)

    # Run simulation and sample potential every dt until t_end
    sim.newRun()

    init_sim(sim, model, geom, v_zmin, params)

    sim.run(params['sim_end'])

    if params['saving_mode'] == 'xdmf':
        del hdf

    result = result.data[0]

    vmean_0um = np.mean(result[1:,0:n_zmin_sample], axis=1)
    vmean_1000um = np.mean(result[1:,n_zmin_sample:], axis=1)
    npt = min(len(vmean_0um),len(vref_0um))

    data = np.zeros((5,npt))
    data[0,:] = np.linspace(0, stop=npt*sim_dt, num=npt, endpoint=False)
    data[1,:] = vmean_0um[0:npt]
    data[2,:] = vref_0um[0:npt]
    data[3,:] = vmean_1000um[0:npt]
    data[4,:] = vref_1000um[0:npt]

    # rms difference
    err_0um = data[2,:] - data[1,:] 
    rms_err_0um = la.norm(err_0um)/np.sqrt(npt)

    err_1000um = data[4,:] - data[3,:]
    rms_err_1000um = la.norm(err_1000um)/np.sqrt(npt)

    return data, rms_err_0um, rms_err_1000um

