# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# -*- coding: utf-8 -*-
#
# Rallpack1 model
# Author Iain Hepburn

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


import steps.quiet
from steps.geom import UNKNOWN_TET
import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.utilities.meshio as meshio
import steps.solver as ssolver

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
    'sim_end':  0.25,           # simulation stop time s
    'EF_dt':    1.0e-5,         # E-field evaluation time step s
    'EF_solver': 'EF_DV_BDSYS', # E-field lin algebra solver
    'SSA_solver': 'Tetexact'
}

def ROIset(x):
    return set_uint(list(x))

def boundary_tets(mesh):
    return (i for i in range(mesh.ntets)
            if any([nb == UNKNOWN_TET for nb in mesh.getTestNeighb(i)]))

def boundary_tris(mesh):
    def btris(tet):
        return [mesh.getTetTriNeighb(tet)[face] for face in range(4)
                if mesh.getTetTetNeighb(tet)[face] == UNKNOWN_TET]

    return (tri for tet in range(mesh.ntets) for tri in btris(tet))

def zminmax_tris(mesh):
    minz = mesh.getBoundMin()[2]
    maxz = mesh.getBoundMax()[2]
    zeps = (maxz - minz)/mesh.ntets

    minz_check = minz + zeps
    maxz_check = maxz - zeps

    zmin_tris = []
    zmin_vset = set()
    zmax_tris = []
    zmax_vset = set()

    for tri in boundary_tris(mesh):
        vertices = mesh.getTri(tri)
        ptverts = [mesh.getVertex(vi) for vi in vertices]

        if all([v[2] <= minz_check for v in ptverts]):
            zmin_tris.append(tri)
            for v in vertices: zmin_vset.add(v)
        elif all([v[2] >= maxz_check for v in ptverts]):
            zmax_tris.append(tri)
            for v in vertices: zmax_vset.add(v)

    return (zmin_tris,zmin_vset,zmax_tris,zmax_vset)

# Find the Z-axially closest and furthest vertices

def radial_extrema(geom, vset):
    r2min = r2max = None
    vmin  = vmax  = None

    for v in vset:
        x = geom.getVertex(v)
        s = x[0]*x[0] + x[1]*x[1]
        if  r2min == None or r2min > s:
            r2min = s
            vmin = v
        if  r2max == None or r2max < s:
            r2max = s
            vmax = v

    return (vmin,vmax)


# Build steps geom object, and report lower and upper surface vertixes
# for later use.

def build_geometry(mesh_path):
    mesh = meshio.loadMesh(mesh_path)[0]

    cyto = sgeom.TmComp('cyto', mesh, range(mesh.ntets))

    (zmin_tris,zmin_vset,zmax_tris,zmax_vset) = zminmax_tris(mesh)
    memb_tris = list(mesh.getSurfTris())
    for tri in zmin_tris: memb_tris.remove(tri)
    for tri in zmax_tris: memb_tris.remove(tri)
    memb = sgeom.TmPatch('memb', mesh, memb_tris, cyto)
    membrane = sgeom.Memb('membrane', mesh, [memb] )

    #mesh.addROI('v_zmin',sgeom.ELEM_VERTEX,ROIset(zmin_vset))
    #mesh.addROI('v_zmin_sample',sgeom.ELEM_VERTEX,ROIset(radial_extrema(mesh,zmin_vset)))
    #mesh.addROI('v_zmax_sample',sgeom.ELEM_VERTEX,ROIset(radial_extrema(mesh,zmax_vset)))
    mesh.addROI('v_zmin',sgeom.ELEM_VERTEX,zmin_vset)
    mesh.addROI('v_zmin_sample',sgeom.ELEM_VERTEX,radial_extrema(mesh,zmin_vset))
    mesh.addROI('v_zmax_sample',sgeom.ELEM_VERTEX,radial_extrema(mesh,zmax_vset))
    return mesh


def build_model(mesh, param):
    mdl = smodel.Model()
    memb = sgeom.castToTmPatch(mesh.getPatch('memb'))

    ssys = smodel.Surfsys('ssys', mdl)
    memb.addSurfsys('ssys')

    L = smodel.Chan('L', mdl)
    Leak = smodel.ChanState('Leak', mdl, L)

    # membrane conductance
    area_cylinder = np.pi * param['diameter'] * param['length']
    L_G_tot = area_cylinder / param['R_M']
    g_leak_sc = L_G_tot / len(memb.tris)
    OC_L = smodel.OhmicCurr('OC_L', ssys, chanstate = Leak, erev = param['E_M'], g = g_leak_sc)

    return mdl


def init_sim(model, mesh, seed, param):
    # previous setup
    # rng = srng.create('r123', 512)
    # rng.initialize(seed)
    # sim = ssolver.Tetexact(model, mesh, rng, True)

    # Create the solver objects
    EFSolver = getattr(ssolver, sim_parameters['EF_solver'])
    if param['SSA_solver'] == 'TetODE':
        sim = ssolver.TetODE(model, mesh, calcMembPot=EFSolver)
        sim.setTolerances(1.0e-6, 1e-6)
    elif param['SSA_solver'] == 'Tetexact':
        rng = srng.create('mt19937', 512)
        rng.initialize(seed)
        sim = ssolver.Tetexact(model, mesh, rng, calcMembPot=EFSolver)
        sim.reset()
        sim.setEfieldDT(param['EF_dt'])
    else :
        raise ValueError('SSA solver ' + param['SSA_solver'] + 'not available')
        
    print('Running Rallpack1 test with ' + param['SSA_solver'])

    # Correction factor for deviation between mesh and model cylinder:
    area_cylinder = np.pi * param['diameter'] * param['length']
    area_mesh_factor = sim.getPatchArea('memb') / area_cylinder

    # Set initial conditions

    memb = sgeom.castToTmPatch(mesh.getPatch('memb'))
    for t in memb.tris: sim.setTriCount(t, 'Leak', 1)

    sim.setMembPotential('membrane', param['E_M'])
    sim.setMembVolRes('membrane', param['R_A'])
    sim.setMembCapac('membrane', param['C_M']/area_mesh_factor)

    v_zmin = mesh.getROIData('v_zmin')
    I = param['Iinj']/len(v_zmin)
    for v in v_zmin: sim.setVertIClamp(v, I)

    return sim

# Run simulation and sample potential every dt until t_end
def run_sim(sim, dt, t_end, vertices, verbose=False):
    N = int(np.ceil(t_end/dt))+1
    result = np.zeros((N, len(vertices)))
    nvert = len(vertices)

    for l in range(N):
        if verbose and not l%100:  print(str(l)+" out of "+str(N))    
        sim.run(l*dt)
        result[l,:] = [sim.getVertV(v) for v in vertices]

    return result


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
    sim = init_sim(model, geom, seed, sim_parameters)

    # grab sample vertices
    zmin_sample = geom.getROIData('v_zmin_sample')
    n_zmin_sample = len(zmin_sample)
    zmax_sample = geom.getROIData('v_zmax_sample')
    n_zmax_sample = len(zmax_sample)
    vertices =  zmin_sample + zmax_sample

    result = run_sim(sim, sim_dt, sim_parameters['sim_end'], vertices, verbose=verbose)

    vmean_0um = np.mean(result[:,0:n_zmin_sample], axis=1)
    vmean_1000um = np.mean(result[:,n_zmin_sample:], axis=1)
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

