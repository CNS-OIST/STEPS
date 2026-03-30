import numpy as np
import steps.interface
from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *

from . import DiffPoint

########################################################################
# Inject species in the whole mesh and let them degrade through a first
# order reaction. Diffusion is present but does not play much of a role.
########################################################################


def setup_parser(parser):
    DiffPoint.setup_parser(parser)
    parser.add_argument("--KCST", type=float, default=100)
    parser.set_defaults(ENDT=0.1, DT=0.01)


########################################################################


def set_parameters(args, params):
    if not args.noscaletime:
        params["timefact"] = args.KCST / 5e1
        params["DT"] = args.DT / params["timefact"]
        params["ENDT"] = args.ENDT / params["timefact"]


########################################################################


def setup_model(mdl, args, params):
    r = ReactionManager()

    DiffPoint.setup_model(mdl, args, params)
    with mdl:
        with mdl.vsys:
            mdl.X > r[1] > None
            r[1].K = args.KCST


########################################################################

setup_mesh = DiffPoint.setup_mesh

########################################################################

setup_sim = DiffPoint.setup_sim

########################################################################

setup_group = DiffPoint.setup_group

########################################################################


def initialize_sim(sim, args, params):
    sim.comp.X.Count = args.INJ
