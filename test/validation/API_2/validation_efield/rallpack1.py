###___license_placeholder___###
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
    'SSA_solver': 'Tetexact'
}

# Build steps geom object, and report lower and upper surface vertixes
# for later use.

def build_geometry(mesh_path):
    mesh = TetMesh.Load(mesh_path)
    with mesh:

        cyto = Compartment.Create(mesh.tets)

        zmin_tris = TriList(tri for tri in mesh.surface if all(v.z <= mesh.bbox.min.z for v in tri.verts))
        zmax_tris = TriList(tri for tri in mesh.surface if all(v.z >= mesh.bbox.max.z for v in tri.verts))

        memb_tris = mesh.surface - (zmin_tris | zmax_tris)
        memb = Patch.Create(memb_tris, cyto)
        membrane = Membrane.Create([memb])

        v_zmin = ROI.Create(zmin_tris.verts)

        # Find the Z-axially closest and furthest vertices
        distFunc = lambda v: la.norm(v[:2])
        v_zmin_sample = ROI.Create(VertList([min(zmin_tris.verts, key=distFunc), max(zmin_tris.verts, key=distFunc)]))
        v_zmax_sample = ROI.Create(VertList([min(zmax_tris.verts, key=distFunc), max(zmax_tris.verts, key=distFunc)]))

    return mesh


def build_model(mesh, param):
    mdl = Model()
    with mdl:
        ssys = SurfaceSystem.Create()
        mesh.memb.addSystem('ssys')

        Leak = SubUnitState.Create()
        L = Channel.Create([Leak])

        # membrane conductance
        area_cylinder = np.pi * param['diameter'] * param['length']
        L_G_tot = area_cylinder / param['R_M']
        g_leak_sc = L_G_tot / len(mesh.memb.tris)
        with ssys:
            OC_L = OhmicCurr.Create(L[Leak], g_leak_sc, param['E_M'])

    return mdl


def create_sim(model, mesh, seed, param):
    # Create the solver objects
    if param['SSA_solver'] == 'TetODE':
        sim = Simulation('TetODE', model, mesh, calcMembPot=MPI.EF_DV_BDSYS)
        sim.setTolerances(1.0e-6, 1e-6)
    elif param['SSA_solver'] == 'Tetexact':
        rng = RNG('mt19937', 512, seed)
        sim = Simulation('Tetexact', model, mesh, rng, calcMembPot=MPI.EF_DV_BDSYS)
        sim.EfieldDT = param['EF_dt']
    else :
        raise ValueError('SSA solver ' + param['SSA_solver'] + 'not available')
        
    print('Running Rallpack1 test with ' + param['SSA_solver'])

    return sim


def init_sim(sim, model, mesh, param):
    SetVerbosity(0)

    # Correction factor for deviation between mesh and model cylinder:
    area_cylinder = np.pi * param['diameter'] * param['length']
    area_mesh_factor = sim.memb.Area / area_cylinder

    # Set initial conditions

    sim.TRIS(mesh.memb.tris).L[model.Leak].Count = 1

    sim.membrane.Potential = param['E_M']
    sim.membrane.VolRes = param['R_A']
    sim.membrane.Capac = param['C_M']/area_mesh_factor

    v_zmin = mesh.v_zmin.verts
    sim.VERTS(v_zmin).IClamp = param['Iinj']/len(v_zmin)


# Returns RMS error, table containing computed end-point voltages
# and reference voltage data.

def run_comparison(seed, mesh_file, v0_datafile, v1_datafile, verbose=False):
    # sample at same interval as rallpack1 reference data
    sim_dt = 5.0e-5

    def snarf(fname):
        F = open(fname, 'r')
        for line in F: yield tuple([float(x) for x in line.split()])
        F.close()

    vref_0um = np.array([v for (t,v) in snarf(v0_datafile)])
    vref_1000um = np.array([v for (t,v) in snarf(v1_datafile)])

    geom = build_geometry(mesh_file)
    model = build_model(geom, sim_parameters)

    # grab sample vertices
    zmin_sample = geom.v_zmin_sample.verts
    n_zmin_sample = len(zmin_sample)
    zmax_sample = geom.v_zmax_sample.verts
    vertices =  zmin_sample + zmax_sample

    sim = create_sim(model, geom, seed, sim_parameters)

    rs = ResultSelector(sim)

    result = rs.VERTS(vertices).V

    sim.toSave(result, dt=sim_dt)

    # Run simulation and sample potential every dt until t_end
    sim.newRun()

    init_sim(sim, model, geom, sim_parameters)

    sim.run(sim_parameters['sim_end'])

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

