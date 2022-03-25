#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

import unittest
import steps.quiet
import steps.geom
import steps.rng
import steps.model
import steps.solver
from functools import partial
from . import distribution_common as D
            
class TestSetConc(unittest.TestCase):
    def setUp(self):
        self.rng = steps.rng.create('r123', 512)
        self.rng.initialize(12345)

    # Convert counts to mol/L

    def count_to_conc(self, count, volume):
        avogadro_constant = 6.022140857e23 # 1/mol
        return count/(volume*1000*avogadro_constant)

    # Check set/getConc pair is sane

    def test_getsetconc(self):
        mesh = D.make_prism_mesh(9, 7.0, 'geometric', 1.3)
        total_vol = mesh.getMeshVolume() # m^3
        molarity = self.count_to_conc(100, total_vol)

        model = D.make_model(mesh)
        solver = steps.solver.Tetexact(model, mesh, self.rng)
        
        solver.setCompConc('interior','A', molarity)
        m2 = solver.getCompConc('interior','A')

        self.assertTrue(abs(molarity-m2)<=1)

    # Produce x, mu pair from calling setCompConc() on
    # a given mesh.

    def set_conc_kernel(self, count, T, mesh):
        total_vol = mesh.getMeshVolume() # m^3
        molarity = self.count_to_conc(count, total_vol)

        model = D.make_model(mesh)
        solver = steps.solver.Tetexact(model, mesh, self.rng)

        n = mesh.countTets()

        # expected means
        mu = [ count*mesh.getTetVol(i)/total_vol for i in range(n)]

        # empirical means
        x = [0]*n
        
        for s in range(T):
            solver.setCompConc('interior', 'A', molarity)
            for t in range(n):
                x[t] += (solver.getTetCount(t,'A')-x[t])/(1+s)

        return x, mu


    def test_distribution(self):
        # even volumes

        D.run_distribution_check(self.set_conc_kernel, alpha=0.005, verbose=False)

        # uneven volumes

        D.run_distribution_check(self.set_conc_kernel, alpha=0.005, ratio=0.1, verbose=False)
        D.run_distribution_check(self.set_conc_kernel, alpha=0.005, ratio=10, verbose=False)

def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(TestSetConc, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
