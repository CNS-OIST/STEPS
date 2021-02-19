########################################################################

# Stochastic degradation-diffusion process.

# AIMS: to verify if zero order reactions are inserted correctly to the SSA
# after bugfix https://github.com/CNS-OIST/HBP_STEPS/pull/146
# Note: 
# 1. This is not validated against analytical result, but 
# only check if the zero order reaction can be found through SSA search.
# Future improvement is needed.
# 2. sim.reset() should be avoided in this test
  
########################################################################
from __future__ import print_function, absolute_import

import steps.model as smod
import steps.geom as sgeom
import steps.rng as srng
import steps.mpi
import steps.mpi.solver as ssolv

import math
import time
import steps.utilities.meshio as meshio
import steps.utilities.geom_decompose as gd
########################################################################

def get_model():
    mdl  = smod.Model()
    A = smod.Spec('A', mdl)
    volsys = smod.Volsys('vsys',mdl)
    R1 = smod.Reac('R1', volsys, lhs = [], rhs = [A])
    R1.setKcst(1e-3)
    return mdl

def get_geom():
    if __name__== "__main__":
        mesh_loc = "meshes/1x1x1.inp"
    else:
        mesh_loc = "validation_rd/meshes/1x1x1.inp"
    mesh = meshio.importAbaqus(mesh_loc, 1e-6)[0]
    comp = sgeom.TmComp("comp", mesh, range(mesh.ntets))
    comp.addVolsys("vsys")
    return mesh

def get_solver(mdl, geom):
    r = srng.create('r123', 1000)
    tet_hosts = gd.binTetsByAxis(geom, steps.mpi.nhosts)
    solver = ssolv.TetOpSplit(mdl, geom, r, ssolv.EF_NONE, tet_hosts)
    return solver

def validate(sim):
    sim.run(1)
    count = sim.getCompCount("comp", "A")
    assert(count > 0)
    if __name__== "__main__" and steps.mpi.rank == 0: print("Count: ", count)

def test_zeroorder_TetOpSplit():
    "Zero order reaction (TetOpSplit)"
    if __name__== "__main__" and steps.mpi.rank == 0: print("Zero order reaction (TetOpSplit)")
    mdl = get_model()
    geom = get_geom()
    sim = get_solver(mdl, geom)
    validate(sim)

if __name__== "__main__":
  test_zeroorder_TetOpSplit()
