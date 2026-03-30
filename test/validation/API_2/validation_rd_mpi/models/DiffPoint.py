import numpy as np
import steps.interface
from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *

########################################################################
# Inject species in the center of the mesh and let it diffuse
########################################################################


def set_parameters(args, params):
    if not args.noscaletime:
        params["timefact"] = args.DCST / 1e-11
        params["DT"] = args.DT / params["timefact"]
        params["ENDT"] = args.ENDT / params["timefact"]


########################################################################


def setup_model(mdl, args, params):
    r = ReactionManager()

    with mdl:
        X = Species.Create()
        vsys = VolumeSystem.Create()
        with vsys:
            Diffusion(X, args.DCST)

    return mdl


########################################################################


def setup_mesh(mesh, mdl, args, params):
    with mesh:
        comp = Compartment.Create(mesh.tets, mdl.vsys)
        params["injTet"] = mesh.tets[mesh.bbox.center]


########################################################################


def setup_sim(sim, args, params):
    injTet = params["injTet"]

    dists, tets = zip(
        *sorted(
            [(np.linalg.norm(injTet.center - tet.center), tet) for tet in sim.geom.tets],
            key=lambda x: x[0],
        )
    )
    tets = TetList(tets, mesh=sim.geom)

    rs = ResultSelector(sim)
    Concs = rs.TETS(tets).X.Conc

    Concs.metaData["distance"] = dists
    Concs.description = "Concentrations"

    Concs.metaData["volume"] = [tet.Vol for tet in tets]

    sim.toSave(Concs, dt=params["DT"])


########################################################################


def setup_group(group, sim, args, params):
    pass


########################################################################


def initialize_sim(sim, args, params):
    sim.TET(params["injTet"]).X.Count = args.INJ
