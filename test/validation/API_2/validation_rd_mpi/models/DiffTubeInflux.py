import numpy as np
import steps.interface
from steps.saving import *
from steps.sim import *

from steps.geom import *
from steps.model import *
from steps.rng import *

from . import DiffPoint

########################################################################
# 1D diffusion in a finite tube with a constant and equal influx of the
# same species of molecule at both ends
########################################################################

set_parameters = DiffPoint.set_parameters

########################################################################

setup_model = DiffPoint.setup_model


def setup_model(mdl, args, params):
    r = ReactionManager()

    DiffPoint.setup_model(mdl, args, params)
    with mdl:
        SA = Species.Create()
        with mdl.vsys:
            SA > r[1] > SA + mdl.X
            r[1].K = args.KCST

    return mdl


########################################################################


def setup_mesh(mesh, mdl, args, params):
    DiffPoint.setup_mesh(mesh, mdl, args, params)
    with mesh:
        zmin, zmax = mesh.bbox.min.z, mesh.bbox.max.z
        params["borderTets"] = TetList(
            tet for tet in mesh.tets if sum(v.z in (zmin, zmax) for v in tet.verts) > 2
        )


########################################################################


def setup_sim(sim, args, params, getConcFunc=None):
    dists, tets = zip(
        *sorted(
            [(tet.center.z, tet) for tet in sim.geom.tets],
            key=lambda x: x[0],
        )
    )
    tets = TetList(tets, mesh=sim.geom)

    rs = ResultSelector(sim)
    if getConcFunc is None:
        Concs = rs.TETS(tets).X.Conc
    else:
        Concs = getConcFunc(rs, tets)

    Concs.metaData["distance"] = dists
    Concs.description = "Concentrations"

    Concs.metaData["volume"] = [tet.Vol for tet in tets]

    sim.toSave(Concs, dt=params["DT"])


########################################################################


def setup_group(group, sim, args, params):
    DiffPoint.setup_group(group, sim, args, params)

    mesh = sim.geom
    borderTris = TriList(
        [
            tri
            for tri in mesh.surface
            if all(v.z in [mesh.bbox.min.z, mesh.bbox.max.z] for v in tri.verts)
        ],
        mesh=mesh,
    )
    borderArea = borderTris.Area / 2
    borderNbTets = len(params["borderTets"]) / 2
    length = mesh.bbox.max.z - mesh.bbox.min.z

    if group is not None:
        group.staticData["BorderArea"] = borderArea
        group.staticData["BorderNbTets"] = borderNbTets
        group.staticData["Length"] = length


########################################################################


def initialize_sim(sim, args, params):
    sim.TETS(params["borderTets"]).SA.Count = 1
