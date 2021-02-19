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
import steps.solver as ssolv

import math
import time
import steps.utilities.meshio as meshio

########################################################################

def get_model():
    mdl  = smod.Model()
    A = smod.Spec('A', mdl)
    volsys = smod.Volsys('vsys',mdl)
    R1 = smod.Reac('R1', volsys, lhs = [], rhs = [A])
    R1.setKcst(1e-3)
    return mdl

def get_geom(geom_type):
    if geom_type == "well_mixed":
        geom = sgeom.Geom()
        comp = sgeom.Comp("comp", geom, 1e-18)
        comp.addVolsys("vsys")
        return geom
    elif geom_type == "tetmesh":
        mesh_loc = None
        if __name__== "__main__":
            mesh_loc = "meshes/1x1x1.inp"
        else:
            mesh_loc = "validation_rd/meshes/1x1x1.inp"
        mesh = meshio.importAbaqus(mesh_loc, 1e-6)[0]
        comp = sgeom.TmComp("comp", mesh, range(mesh.ntets))
        comp.addVolsys("vsys")
        return mesh
    else:
        assert(False)

def get_solver(solver_type, mdl, geom):
    r = srng.create('r123', 1000)
    if solver_type == "Wmrk4":
        solver = ssolv.Wmrk4(mdl, geom)
        solver.setDT(1e-5)
        return solver
    elif solver_type == "Wmdirect":
        solver = ssolv.Wmdirect(mdl, geom, r)
        return solver
    if solver_type == "Wmrssa":
        solver = ssolv.Wmrssa(mdl, geom, r)
        return solver
    elif solver_type == "Tetexact":
        solver = ssolv.Tetexact(mdl, geom, r)
        return solver
    else:
        assert(False)

def validate(sim):
    sim.run(1)
    count = sim.getCompCount("comp", "A")
    assert(count > 0)
    if __name__== "__main__": print("Count: ", count)

def test_zeroorder_wmrk4():
    "Zero order reaction (Wmrk4)"
    if __name__== "__main__": print("Zero order reaction (Wmrk4)")
    mdl = get_model()
    geom = get_geom("well_mixed")
    sim = get_solver("Wmrk4", mdl, geom)
    validate(sim)

def test_zeroorder_wmdirect():
    "Zero order reaction (Wmdirect)"
    if __name__== "__main__": print("Zero order reaction (Wmdirect)")
    mdl = get_model()
    geom = get_geom("well_mixed")
    sim = get_solver("Wmdirect", mdl, geom)
    validate(sim)

def test_zeroorder_wmrssa():
    "Zero order reaction (Wmrssa)"
    if __name__== "__main__": print("Zero order reaction (Wmrssa)")
    mdl = get_model()
    geom = get_geom("well_mixed")
    sim = get_solver("Wmrssa", mdl, geom)
    validate(sim)

def test_zeroorder_tetexact():
    "Zero order reaction (Tetexact)"
    if __name__== "__main__": print("Zero order reaction (Tetexact)")
    mdl = get_model()
    geom = get_geom("tetmesh")
    sim = get_solver("Tetexact", mdl, geom)
    validate(sim)

if __name__== "__main__":
  test_zeroorder_wmrk4()
  test_zeroorder_wmdirect()
  test_zeroorder_wmrssa()
  test_zeroorder_tetexact()

