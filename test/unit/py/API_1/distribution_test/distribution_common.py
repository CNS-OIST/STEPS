#!/usr/bin/env python
# -*- coding: utf-8 -*-

####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

import math
import steps.quiet
import steps.geom
import steps.model

import numpy as np

# Approximate erf and inverf from Winitzki (2008), "A handy approximation
# for the error function and its inverse".
# <URL: https://sites.google.com/site/winitzki/sergei-winitzkis-files/erf-approx.pdf>

def approx_inverf(x):
    c1 = 6.80272108844  # 1/a, a=0.147
    c2 = 4.3307467508   # 2/π·c₁
    log_omx2 = np.log(1-x*x)
    t1 = c1*log_omx2
    t2 = c2+0.5*log_omx2
    r = np.sqrt(np.sqrt(t2*t2-t1)-t2)
    if x<0:
        return -r
    else:
        return r
    
def approx_erf(x):
    a  = 0.147
    c3 = 1.27323954474 # 4/π

    x2 = x*x
    r = np.sqrt(1-np.exp(-x2*(c3+a*x2)/(1+a*x2)))
    if x<0:
        return -r
    else:
        return r

def p_twotailed(z):
    return 1-approx_erf(abs(z)*0.707106781187)

# 'Monadic' min and argmin: None values ignored,
# None returned if no non-None arguments.

def m_argmin(xs):
    nns = [x for x in enumerate(xs) if x[1]!=None]
    if not nns: return None
    return min(nns, key = lambda i: i[1])[0]

def m_min(xs, **kwargs):
    try: return min((x for x in xs if x != None), **kwargs)
    except ValueError: return None

# List utilities

def partial_sums(xs):
    s = 0
    for x in xs:
        s += x
        yield s

def flatten(xss):
    return [x for xs in xss for x in xs]

# Make two species, A and B, that live
# in volumes and surfaces respectively.
#
#
# Add compartment for full intetior volume
# and patch for full mesh surface that will
# hold species A and B respectively.

def make_model(mesh):
    model = steps.model.Model()
    vsys = steps.model.Volsys('vsys', model)
    ssys = steps.model.Surfsys('ssys', model)

    A = steps.model.Spec('A', model)
    steps.model.Diff('diff_A', vsys, A)
    
    B = steps.model.Spec('B', model)
    steps.model.Diff('diff_B', ssys, B)

    #mesh.delComp('interior')
    comp = steps.geom.TmComp('interior', mesh, range(mesh.countTets()))
    comp.addVolsys('vsys')

    #mesh.delPatch('surface')
    patch = steps.geom.TmPatch('surface', mesh, mesh.getSurfTris(), icomp = comp)
    patch.addSurfsys('ssys')

    return model


# Construct a mesh consisting of stacked triangular prisms in the z-axis,
# at heights given by zlevels. Each prism is decomposed into three tetrahedra.

def prism_mesh_by_z(zlevels = [0,1], scale = 1):
    nprism = len(zlevels)-1

    # make vertices
    verts = [(x,y,z) for z in zlevels for y in range(0,2) for x in range(0,2-y)]
                
    # shear and scale ...
    sqrt3o2 = np.sqrt(3.0)/2
    verts = [(x*scale+0.5*y*scale, sqrt3o2*scale*y, scale*z) for (x,y,z) in verts]
   
    # vertex at (x,y,z), x+y <= 1.
    def xyz_to_vidx(x, y, z):
        return z*3 + x + 2*y

    def vidx_to_xyz(v):
        z = v/3
        x = v%3
        if x==2:
            x = 0
            y = 1

        return (x,y,z)

    # tets in jth simple triangular prism
    def tet_verts(j):
        vs = [[0,1,2,3], [1,2,3,5], [1,3,4,5]]
        return [[v+3*j for v in tetv] for tetv in vs]

    tets = [tet for level in range(0,len(zlevels)-1)
                for tet in tet_verts(level)]

    return (verts,tets)

# Make a prism with n/3 divisions (with 3 tetrahedra each) with given total
# volume and the specified ratio between the volumes of first and last
# segments. Segments in between have volume as specified by distribution.

def make_prism_mesh(n, total_volume, distribution, ratio = 1):
    n //= 3

    if distribution=='constant':
        zstep = [1]*n
    elif distribution=='linear':
        a = 2.0/(n-1)*(ratio-1.0)/(ratio+1.0)
        zstep = [1+a*(j-(n-1)/2.0) for j in range(0,n)]
    elif distribution=='geometric':
        a = np.power(ratio,1/(n-1.0))
        zstep = [np.power(a,j) for j in range(0,n)]
    else:
        raise ValueError('unrecognized volume distribution scheme')

    zlevels = [0] + [z for z in partial_sums(zstep)]
    scale = np.power(4*total_volume/(zlevels[-1]*np.sqrt(3)), 1/3.0)
    (vx,tetx) = prism_mesh_by_z(zlevels, scale)

    return steps.geom.Tetmesh(flatten(vx), flatten(tetx), [])


# Compare observed means x against expected means mu, where the observed
# means result from T distribution trials, and the distribution is
# considered to be a process of rounding-down to a count followed by
# with- or without-replacement sampling.
#
# Report critical value for rejection of null hypothesis (that the
# distribution does follow such a sampling procedure), after applying a
# Bonferonni correction. Also returns list of p-values and  z-scores for
# each element of x. 

def run_stats(x, mu, T, sample_type):
    if sample_type not in ['wor', 'wr']:
        raise ValueError("sample type must be 'wr' or 'wor'")

    n = len(x)
    if n != len(mu):
        raise ValueError("observed and expected vector lengths differ")

    # If c is the count in a tetrahedron after distribution, and µ is the
    # fractional count that should be met on average, then the mean of T trials
    # of c should be distributed as 1/T·Binomial(T,{µ}) where {µ} = µ-⌊µ⌋ is
    # the fractional part of µ, if the sampling is done without replacement.
    # With replacement, this becomes 1/T·Binomial(Tm,{µ}/m) where m = n - Σ⌊µ⌋ is
    # the population size the sample is drawn from.  These have standard
    # deviations of √({μ}(1-{μ})/T) and √({μ}(1-{μ}/m)/T) respectively.
    #
    # Assuming a normal approximation to the binomial distribution for large T,
    # we can reject the hypothesis that the samples come from a with- or
    # without-replacement sampling process by computing the two-sided p-value
    # associated with the Z-score, and applying the Bonferroni correction.
    # (Owing to the linear dependence in the sample results, this will be an
    # over-correction.)
    #
    # The normal approximation should not be used for Binomial(N,p) if Np or
    # N(1-p)  < 5. In this cases, report 'None' for p-values.

    if sample_type=='wor':
        r = 1
    else:
        r = sum(mu) - sum([int(y) for y in mu])

    expected_sd = [np.sqrt(r*psi*(1-psi)/T) for y in mu for psi in [y-int(y)]]

    # use min of {µ} and 1-{µ} in mask test to avoid 
    pmask = [T*psi > 5 and T*(1-psi) > 5 for y in mu for psi in [y-int(y)]]
    b = sum([1 for k in pmask if k])

    z = [0]*n
    p = [0]*n
    for i in range(n):
        deviation = abs(x[i]-mu[i])
        z[i] = deviation/expected_sd[i]
        p[i] = p_twotailed(z[i]) if pmask[i] else None

    pmin = m_min(p)
    if pmin!=None: pmin *= b

    return pmin, p, z
        
    
def run_distribution_check(kernel, alpha, ratio=1, min_trials=500, max_trials=10000, verbose=False):
    # Test equal count distributions

    for n in [30, 1200]:
        for count in [8*n+n/5, 0.99*n, 4.1*n]:
            if ratio==1:
                mesh = make_prism_mesh(n, 1.0, 'constant')
            else:
                mesh = make_prism_mesh(n, 1.0, 'geometric', ratio)

            T = min(max_trials, max(min_trials, 500000/n))
            x, mu = kernel(count, T, mesh)

            p_crit, p, z = run_stats(x, mu, T, 'wr')
            i = m_argmin(p) 

            if verbose:
                print('# bins={0}, count={1}, trials={2}\nx,mu,z'.format(n, count, T))
                for a, b, c in zip(x, mu, z):
                    print('{0},{1},{2}'.format(a, b, c))

            message = 'Failure with {0} bins, {1} molecules, over {2} trials:\n'.format(n, count, T)

            if p_crit==None:
                message += 'Too few trials to determine representative z-score.'
            elif p_crit<alpha:
                message += 'bin {0}: x={1}, mu={2}, p-value={3}, z-score={4}'.format(i, x[i], mu[i], p[i], z[i])

            assert p_crit!=None and p_crit>alpha, message

