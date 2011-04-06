# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2011 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Example: Diffusion boundary
# http://steps.sourceforge.net/manual/diffusion_boundary.html

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as solvmod
import steps.utilities.meshio as meshio

import numpy
from pylab import *
	
########################################################################

def gen_model():
    
    # Create the model container object
    mdl = smodel.Model()
    
    # Create the chemical species
    X = smodel.Spec('X', mdl)
    Y = smodel.Spec('Y', mdl)
    
    # Create separate volume systems for compartments A and B
    vsysA = smodel.Volsys('vsysA', mdl)
    vsysB = smodel.Volsys('vsysB', mdl)
    
    # Describe diffusion of molecules in compartments A and B
    # NOTE: diffusion is not defined for species X in compartment B
    diff_X_A = smodel.Diff('diff_X_A', vsysA, X, dcst = 0.1e-9)
    diff_X_B = smodel.Diff('diff_X_B', vsysB, X, dcst = 0.1e-9)
    diff_Y_A = smodel.Diff('diff_Y_A', vsysA, Y, dcst = 0.1e-9)    
    diff_Y_B = smodel.Diff('diff_Y_B', vsysB, Y, dcst = 0.1e-9) 
    
    # Return the container object
    return mdl
	
########################################################################

def gen_geom():
    mesh = meshio.loadMesh('meshes/cyl_len10_diam1')[0]
    ntets = mesh.countTets()
    
    tets_compA = []
    tets_compB = []
    tris_compA = set()
    tris_compB = set()
    z_max = mesh.getBoundMax()[2]
    z_min = mesh.getBoundMin()[2]
    z_mid = z_min+(z_max-z_min)/2.0
    
    for t in range(ntets):
        # Fetch the z coordinate of the barycenter
        barycz = mesh.getTetBarycenter(t)[2]
        # Fetch the triangle indices of the tetrahedron, a tuple of length 4:
        tris = mesh.getTetTriNeighb(t)
        if barycz < z_mid: 
            tets_compA.append(t)
            tris_compA.add(tris[0])
            tris_compA.add(tris[1])
            tris_compA.add(tris[2])
            tris_compA.add(tris[3])
        else: 
            tets_compB.append(t)
            tris_compB.add(tris[0])
            tris_compB.add(tris[1])
            tris_compB.add(tris[2])
            tris_compB.add(tris[3])
    
    # Create the mesh compartments
    compA = sgeom.TmComp('compA', mesh, tets_compA)
    compB = sgeom.TmComp('compB', mesh, tets_compB)
    compA.addVolsys('vsysA')
    compB.addVolsys('vsysB')
    
    # Make the diff boundary tris with the intersection and convert to list
    tris_DB = tris_compA.intersection(tris_compB)
    tris_DB = list(tris_DB)

    # Create the diffusion boundary between compA and compB
    diffb = sgeom.DiffBoundary('diffb', mesh, tris_DB)
    
    return mesh, tets_compA, tets_compB

########################################################################

mdl = gen_model()
mesh, tets_compA, tets_compB = gen_geom()

rng = srng.create_mt19937(256) 
rng.initialize(432) 

sim = solvmod.Tetexact(mdl, mesh, rng)
sim.reset()

tpnts = numpy.arange(0.0, 0.101, 0.001)
ntpnts = tpnts.shape[0]

# And fetch the total number of tets to make the data structures
ntets = mesh.countTets()
    
# Create the data structures: time points x tetrahedrons
resX = numpy.zeros((ntpnts, ntets))
resY = numpy.zeros((ntpnts, ntets))

tetx = mesh.findTetByPoint([0, 0, -4.99e-6])
tety = mesh.findTetByPoint([0,0, 4.99e-6])

sim.setTetCount(tetx , 'X', 1000)
sim.setTetCount(tety, 'Y', 500)

sim.setDiffBoundaryDiffusionActive('diffb', 'Y', True)

for i in range(ntpnts):
    sim.run(tpnts[i])
    for k in range(ntets):
        resY[i, k] = sim.getTetCount(k, 'Y')
        resX[i, k] = sim.getTetCount(k, 'X')

########################################################################

def plot_binned(t_idx, bin_n = 100):
    if (t_idx > tpnts.size):
        print "Time index is out of range."
        return
    
    z_tets = numpy.zeros(ntets)  
    zbound_min = mesh.getBoundMin()[2]
      
    # Now find the distance of the centre of the tets to the Z lower face
    for i in range(ntets):
        baryc = mesh.getTetBarycenter(i)
        z = baryc[2] - zbound_min
        # Convert to microns
        z_tets[i] = z*1.0e6
    
    # Find the maximum and minimum z of all tetrahedrons 
    z_max = z_tets.max()
    z_min = z_tets.min()
    
    # Set up the bin structures, recording the individual bin volumes
    z_seg = (z_max-z_min)/bin_n
    bin_mins = numpy.zeros(bin_n+1)
    z_tets_binned = numpy.zeros(bin_n)
    bin_vols = numpy.zeros(bin_n)    
    
    # Now sort the counts into bins for species 'X' 
    z = z_min
    for b in range(bin_n + 1):
        bin_mins[b] = z
        if (b!=bin_n): z_tets_binned[b] = z +z_seg/2.0
        z+=z_seg
    bin_counts = [None]*bin_n
    for i in range(bin_n): bin_counts[i] = []
    for i in range((resX[t_idx].size)):
        i_z = z_tets[i]
        for b in range(bin_n):
            if(i_z>=bin_mins[b] and i_z<bin_mins[b+1]):
                bin_counts[b].append(resX[t_idx][i])
                bin_vols[b]+=sim.getTetVol(i)
                break
    
    # Convert to concentration in arbitrary units              
    bin_concs = numpy.zeros(bin_n)
    for c in range(bin_n): 
        for d in range(bin_counts[c].__len__()):
            bin_concs[c] += bin_counts[c][d]
        bin_concs[c]/=(bin_vols[c]*1.0e18)
    
    t = tpnts[t_idx]
    
    scatter(z_tets_binned, bin_concs, label = 'X', color = 'blue')
    
    # Repeat the process for species 'Y'- separate from 'X' for clarity:
    z = z_min
    for b in range(bin_n + 1):
        bin_mins[b] = z
        if (b!=bin_n): z_tets_binned[b] = z +z_seg/2.0
        z+=z_seg
    bin_counts = [None]*bin_n
    for i in range(bin_n): bin_counts[i] = []
    for i in range((resY[t_idx].size)):
        i_z = z_tets[i]
        for b in range(bin_n):
            if(i_z>=bin_mins[b] and i_z<bin_mins[b+1]):
                bin_counts[b].append(resY[t_idx][i])
                break
                
    bin_concs = numpy.zeros(bin_n)
    for c in range(bin_n): 
        for d in range(bin_counts[c].__len__()):
            bin_concs[c] += bin_counts[c][d]
        bin_concs[c]/=(bin_vols[c]*1.0e18)
    
    scatter(z_tets_binned, bin_concs, label = 'Y', color = 'red')
    xlabel('Z axis (microns)', fontsize=16)
    ylabel('Bin concentration (N/m^3)', fontsize=16)
    ylim(0)    
    xlim(0, 10)
    legend(numpoints=1)
    show()

########################################################################

plot_binned(100, 50)

