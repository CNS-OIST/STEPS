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

import numpy

from steps import stepslib

from . import model as nmodel
from . import geom as ngeom
from . import utils as nutils

RATE_ZSCORE_THRESH = 3
RATE_NBREACS_THRESH = 10


def Check(sim, printmsgs=True):
    """
    Check that the model is consistent, return errors for major issues and warnings for
    potential issues.
    """

    def getLocations(reac, elem):
        if elem.loc is None:
            return reac._parents[nmodel.VolumeSystem].locations
        else:
            locations = reac._parents[nmodel.SurfaceSystem].locations
            if elem.loc == nmodel.Location.SURF:
                return [l for l in locations if isinstance(l, ngeom.Patch)]
            elif elem.loc == nmodel.Location.IN:
                return [l.innerComp for l in locations if isinstance(l, ngeom.Patch)]
            elif elem.loc == nmodel.Location.OUT:
                return [l.outerComp for l in locations if isinstance(l, ngeom.Patch)]
            else:
                raise NotImplementedError()

    errors, warnings = [], []
    LHSObjs, RHSObjs = set(), set()
    reacRates = []
    for reac in sim.model._getChildrenOfType(nmodel.Reaction):
        if reac._added:
            for lhs, rhs, rm in reac[nmodel.Reaction._FwdSpecifier]._LRP:
                tmpLhs, tmpRhs = set(), set()
                for e in lhs:
                    tmpLhs |= set((e._elem, loc) for loc in getLocations(reac, e))
                for e in rhs:
                    tmpRhs |= set((e._elem, loc) for loc in getLocations(reac, e))
                LHSObjs |= tmpLhs - tmpRhs
                RHSObjs |= tmpRhs - tmpLhs
            if reac._bidir:
                for lhs, rhs, rm in reac[nmodel.Reaction._BkwSpecifier]._LRP:
                    tmpLhs, tmpRhs = set(), set()
                    for e in lhs:
                        tmpLhs |= set((e._elem, loc) for loc in getLocations(reac, e))
                    for e in rhs:
                        tmpRhs |= set((e._elem, loc) for loc in getLocations(reac, e))
                    LHSObjs |= tmpLhs - tmpRhs
                    RHSObjs |= tmpRhs - tmpLhs
            # TODO Later release: Checking reaction with real complexes

            # Check reaction rates
            for stepsReac in reac['fwd']:
                stepsReac = stepsReac._stepsReac
                if isinstance(stepsReac, (stepslib._py_Reac, stepslib._py_SReac)):
                    reacRates.append((reac, stepsReac, stepsReac.getKcst()))
            if reac._bidir:
                for stepsReac in reac['bkw']:
                    stepsReac = stepsReac._stepsReac
                    if isinstance(stepsReac, (stepslib._py_Reac, stepslib._py_SReac)):
                        reacRates.append((reac, stepsReac, stepsReac.getKcst()))
        else:
            errors.append(f'Reaction {reac} was not added to STEPS.')

    # Inspect objects that appear only on the RHS of reactions
    for obj, loc in RHSObjs - LHSObjs:
        if isinstance(obj, nmodel.ComplexState):
            errors.append(f'Complex state {obj} in {loc} is only ever present on the RHS of reactions.')
        else:
            warnings.append(f'{obj} in {loc} is only ever present on the RHS of reactions.')

    # Check that rates are non-zero and check for outliers
    ratesByOrder = {}
    for reac, sr, rate in reacRates:
        order = sr.getOrder()
        ratesByOrder.setdefault(order, []).append((reac, sr, rate))

    for order, reacs in ratesByOrder.items():
        rates = numpy.array([rate for _, _, rate in reacs if rate > 0])
        if len(rates) > 0:
            gmeanRate = numpy.exp(numpy.mean(numpy.log(rates)))
            gstdDevRate = numpy.exp(numpy.std(numpy.log(rates)))
        for reac, sr, rate in reacs:
            if rate == 0:
                warnings.append(
                    f'The rate of steps reaction {sr.getID()} declared in reaction {reac} is set to 0.'
                )
            # Only check zscore if there is some threshold amount of reactions
            elif len(reacs) >= RATE_NBREACS_THRESH and gstdDevRate > 1:
                gzscore = numpy.log(rate / gmeanRate) / numpy.log(gstdDevRate)
                if abs(gzscore) > RATE_ZSCORE_THRESH:
                    warnings.append(
                        f'The rate of steps reaction {sr.getID()} declared in reaction '
                        f'{reac} is very different ({rate}, gzscore={gzscore}) from '
                        f'other reactions of the same order (gmean={gmeanRate}, '
                        f'gstddev={gstdDevRate}).'
                    )

    # Check model volume system constraints
    for reac, vsys, ssys, loc in sim.model.volSysConstraints:
        for patch in ssys.locations:
            if isinstance(patch, ngeom.Patch):
                comp = patch.innerComp if loc == nmodel.Location.IN else patch.outerComp
                if (vsys, loc) not in comp.systems:
                    errors.append(
                        f'Reac {reac} implied that volume system {vsys} was always '
                        f'associated to the location {loc} of surface system {ssys}. '
                        f'It is not the case for {patch}.'
                    )

    ok = len(errors) == 0
    if printmsgs:
        nutils._print('Model checking:', 1)
        if len(errors) == len(warnings) == 0:
            nutils._print('No errors were found', 1)
        if len(errors) > 0:
            nutils._print('Errors:', 1)
        for e in errors:
            nutils._print(e, 1)
        if len(warnings) > 0:
            nutils._print('Warnings:', 1)
        for w in warnings:
            nutils._print(w, 1)
        return ok
    else:
        return ok, errors, warnings
