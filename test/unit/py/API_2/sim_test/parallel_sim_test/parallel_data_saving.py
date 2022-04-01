####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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

""" Unit tests for parallel data saving."""

import unittest

from steps import interface

from steps.geom import *
from steps.rng import *
from steps.sim import *

import steps.API_1.mpi.solver as ssolver
import steps.API_1.mpi as smpi
import steps.API_1.rng as srng
import steps.API_1.utilities.geom_decompose as gd

from ..test_dataSaving import TetSimDataSaving

class TetOpSplitSimDataSaving(TetSimDataSaving):
    """Test data access, setting, and saving with tetmeshes."""

    def _get_API2_Sim(self, nmdl, ngeom):
        nrng = RNG('mt19937', 512, self.seed)
        part = LinearMeshPartition(ngeom, 1, MPI.nhosts, 1)
        nsim = Simulation('TetOpSplit', nmdl, ngeom, nrng, MPI.EF_DEFAULT, part)
        nsim.EfieldDT = self.efielddt
        return nsim

    def _get_API1_Sim(self, omdl, ogeom):
        orng=srng.create('mt19937',512)
        orng.initialize(self.seed)
        tet_hosts = gd.linearPartition(ogeom, [1, smpi.nhosts, 1])
        tri_hosts = gd.partitionTris(ogeom, tet_hosts, (self.ptris + self.ptris2).indices)
        if smpi.rank == 0:
            gd.validatePartition(ogeom, tet_hosts)
        osim = ssolver.TetOpSplit(omdl, ogeom, orng, ssolver.EF_DEFAULT, tet_hosts, tri_hosts)
        osim.setEfieldDT(self.efielddt)
        return osim

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TetOpSplitSimDataSaving, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
