####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

import pickle
import sys

if __name__ == '__main__':
    # Not using argparse here because this module is not meant to be ran by users
    paramPath, sitepackages = sys.argv[sys.argv.index('--') + 1:]

    sys.path.append(sitepackages)

    from stepsblender.blenderloader import HDF5BlenderLoader

    with open(paramPath, 'rb') as f:
        parameters = pickle.load(f)

    ld = HDF5BlenderLoader(parameters=parameters)
