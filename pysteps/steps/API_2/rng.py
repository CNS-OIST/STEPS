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

from steps import stepslib

__all__ = [
    'RNG',
]


class RNG:
    """Random number generator class

    :param algoStr: Algorithm used to generate random numbers (see below).
    :type algoStr: str
    :param buffSz: Pre-allocated buffer size
    :type buffSz: int
    :param seed: Seed for the random number generator
    :type seed: int

    Available algorithms:
        - ``'mt19937'`` (Mersenne Twister, based on the original mt19937.c)
        - ``'r123'``

    Method and attributes are the same as in :py:class:`steps.API_1.rng.RNG`.
    """

    def __init__(self, algoStr, buffSz, seed):
        self.stepsrng = stepslib._py_rng_create(algoStr, buffSz)
        self.stepsrng.initialize(seed)

    def __getattr__(self, name):
        return getattr(self.stepsrng, name)

    def __call__(self, *args, **kwargs):
        return self.stepsrng.__call__(*args, **kwargs)
