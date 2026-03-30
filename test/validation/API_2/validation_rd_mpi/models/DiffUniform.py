import numpy as np
import steps.interface
from steps.saving import *
from steps.sim import *

from steps.geom import *
from steps.model import *
from steps.rng import *

from . import DiffPoint

########################################################################
# Inject species in the whole mesh and let them diffuse
########################################################################

set_parameters = DiffPoint.set_parameters

########################################################################

setup_model = DiffPoint.setup_model

########################################################################

setup_mesh = DiffPoint.setup_mesh

########################################################################

setup_sim = DiffPoint.setup_sim

########################################################################

setup_group = DiffPoint.setup_group

########################################################################


def initialize_sim(sim, args, params):
    sim.comp.X.Count = args.INJ
