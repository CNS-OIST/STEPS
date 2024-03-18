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

# Generate data saving files for the currently installed STEPS version
# This is used to then check that files created with older versions can still be read
# in newer versions.

import steps.interface

from sim_test import test_dataSaving as tds

import os
import shutil

FILEDIR = os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':
    allTests = [
        ('testFileSaving', ['dat']),
        ('testFileSavingSQLite', ['db']),
        ('testFileSavingHDF5', ['h5']),
    ]
    version_dir = os.path.join(FILEDIR, 'saved_files', str(steps.__version__))
    os.makedirs(version_dir, exist_ok=True)
    for test, extensions in allTests:
        rsTests = tds.SimDataSaving()
        rsTests.setUp()
        # Run test
        getattr(rsTests, test)()
        i = 0
        for path in rsTests.createdFiles:
            if os.path.isfile(path):
                *_, ext = path.split('.')
                if ext in extensions:
                    shutil.copyfile(path, os.path.join(version_dir, f'{test}_{i}.{ext}'))
                    i += 1
        rsTests.tearDown()
