import numpy as np
import steps.interface
from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *

from . import DiffPoint, DiffTubeInflux

########################################################################
# Two separated reactants diffuse towards each other in the center of a
# cylinder and annihilate.
########################################################################

KCST = 1e11  # annihilation reaction rate


########################################################################

set_parameters = DiffTubeInflux.set_parameters

########################################################################


def setup_model(mdl, args, params):
    r = ReactionManager()

    DiffPoint.setup_model(mdl, args, params)
    with mdl:
        Y = Species.Create()
        with mdl.vsys:
            Diffusion(Y, args.DCST)

            mdl.X + Y > r[1] > None
            r[1].K = KCST


########################################################################


def setup_mesh(mesh, mdl, args, params):
    DiffPoint.setup_mesh(mesh, mdl, args, params)
    with mesh:
        centerz = mesh.bbox.center.z
        params["half_tets1"] = TetList(tet for tet in mesh.tets if tet.center.z < centerz)
        params["half_tets2"] = mesh.tets - params["half_tets1"]
        params["half_tets1_vol"] = np.array([tet.Vol for tet in params["half_tets1"]])
        params["half_tets2_vol"] = np.array([tet.Vol for tet in params["half_tets2"]])


########################################################################


def setup_sim(sim, args, params):
    DiffTubeInflux.setup_sim(
        sim,
        args,
        params,
        getConcFunc=lambda rs, tets: rs.JOIN(
            rs.TET(tet).X.Conc + rs.TET(tet).Y.Conc for tet in tets
        ),
    )


########################################################################


def setup_group(group, sim, args, params):
    DiffPoint.setup_group(group, sim, args, params)

    mesh = sim.geom
    length = mesh.bbox.max.z - mesh.bbox.min.z
    if group is not None:
        group.staticData["Length"] = length


########################################################################


def distribute(weights, qte, method="uniform"):
    total = np.sum(weights)
    match method:
        case "binomial":
            res = []
            for w in weights:
                res.append(np.random.binomial(qte, min(1, w / total)))
                total -= w
                qte -= res[-1]
            return np.array(res)
        case "uniform":
            probs = weights / total
            res = np.floor(qte * probs).astype(int)
            missing = qte - np.sum(res)
            inds = np.random.choice(range(len(weights)), missing, replace=False, p=probs)
            for i in inds:
                res[i] += 1
            return res
        case _:
            raise NotImplementedError()


def initialize_sim(sim, args, params):
    sim.TETS(params["half_tets1"]).X.Count = distribute(params["half_tets1_vol"], int(args.INJ))
    sim.TETS(params["half_tets2"]).Y.Count = distribute(params["half_tets2_vol"], int(args.INJ))
