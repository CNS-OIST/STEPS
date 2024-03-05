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

import numpy
import numbers
from scipy.stats import chi2

from steps import stepslib

from . import model as nmodel
from . import geom as ngeom
from . import utils as nutils

RATE_ZSCORE_THRESH = 3
RATE_NBREACS_THRESH = 10

ACCEPTABLE_JUMP_FRACTION = 0.05
MESH_SMALL_TET_RATIO = 0.1
MESH_SMALL_TET_VOL = 20e-9 ** 3 / (numpy.sqrt(2) * 6)
VESICLE_TET_FULL_OVERLAP_RATIO = 0.1


def Check(sim, printmsgs=True):
    """
    Check that the model is consistent, return errors for major issues and warnings for
    potential issues.
    """

    def getLocations(reac, elem, geom):
        if elem.loc is None:
            return reac._parents[nmodel.VolumeSystem].locations
        else:
            locations = reac._parents[nmodel.SurfaceSystem].locations
            if elem.loc == nmodel.Location.SURF:
                return [l for l in locations if isinstance(l, ngeom.Patch)]
            elif elem.loc == nmodel.Location.IN:
                compLocs = []
                for l in locations:
                    if isinstance(l, ngeom.Patch):
                        compLocs.append(l.innerComp)
                    elif isinstance(l, nmodel.Raft):
                        compLocs += [
                            patch.innerComp for patch in geom.ALL(ngeom.Patch) if patch.innerComp is not None
                        ]
                return compLocs
            elif elem.loc == nmodel.Location.OUT:
                compLocs = []
                for l in locations:
                    if isinstance(l, ngeom.Patch):
                        compLocs.append(l.outerComp)
                    elif isinstance(l, nmodel.Raft):
                        compLocs += [
                            patch.outerComp for patch in geom.ALL(ngeom.Patch) if patch.outerComp is not None
                        ]
                    elif isinstance(l, nmodel.Vesicle):
                        compLocs += geom.ALL(ngeom.Compartment)
                return compLocs
            elif elem.loc == nmodel.Location.VESSURF:
                return [l for l in locations if isinstance(l, nmodel.Vesicle)]
            elif elem.loc == nmodel.Location.RAFTSURF:
                return [l for l in locations if isinstance(l, nmodel.Raft)]
            else:
                raise NotImplementedError()

    def getReactionsData(pred, errors, warnings):
        LHSObjs, RHSObjs = set(), set()
        reacRates = []
        for reac in sim.model._getChildrenOfType(nmodel.Reaction):
            if reac._added:
                if pred(reac):
                    for lhs, rhs, rm in reac[nmodel.Reaction._FwdSpecifier]._LRP:
                        tmpLhs, tmpRhs = set(), set()
                        for e in lhs:
                            tmpLhs |= set((e._elem, loc)
                                          for loc in getLocations(reac, e, sim.geom))
                        for e in rhs:
                            tmpRhs |= set((e._elem, loc)
                                          for loc in getLocations(reac, e, sim.geom))
                        LHSObjs |= tmpLhs - tmpRhs
                        RHSObjs |= tmpRhs - tmpLhs
                    if reac._bidir:
                        for lhs, rhs, rm in reac[nmodel.Reaction._BkwSpecifier]._LRP:
                            tmpLhs, tmpRhs = set(), set()
                            for e in lhs:
                                tmpLhs |= set((e._elem, loc)
                                              for loc in getLocations(reac, e, sim.geom))
                            for e in rhs:
                                tmpRhs |= set((e._elem, loc)
                                              for loc in getLocations(reac, e, sim.geom))
                            LHSObjs |= tmpLhs - tmpRhs
                            RHSObjs |= tmpRhs - tmpLhs
                    # TODO Later release: Checking reaction with real complexes

                    # Check reaction rates
                    for stepsReac in reac['fwd']:
                        stepsReac = stepsReac._stepsReac
                        if isinstance(stepsReac, (stepslib._py_Reac, stepslib._py_SReac)):
                            reacRates.append(
                                (reac, stepsReac, stepsReac.getKcst()))
                    if reac._bidir:
                        for stepsReac in reac['bkw']:
                            stepsReac = stepsReac._stepsReac
                            if isinstance(stepsReac, (stepslib._py_Reac, stepslib._py_SReac)):
                                reacRates.append(
                                    (reac, stepsReac, stepsReac.getKcst()))
            else:
                errors.append(f'Reaction {reac} was not added to STEPS.')

        return LHSObjs, RHSObjs, reacRates

    errors, warnings = [], []
    LHSObjs, RHSObjs, reacRates = getReactionsData(
        lambda r: True, errors, warnings)
    # Inspect objects that appear only on the RHS of reactions, considering all reactions
    for obj, loc in RHSObjs - LHSObjs:
        if isinstance(obj, nmodel.ComplexState):
            errors.append(
                f'Complex state {obj} in {loc} is only ever present on the RHS of reactions.')
        else:
            warnings.append(
                f'{obj} in {loc} is only ever present on the RHS of reactions.')

    # Inspect objects that appear only on the RHS of reactions, considering only conventional reactions
    LHSObjs_noves, RHSObjs_noves, _ = getReactionsData(
        lambda r: not (r._isVesicleSurfReac() or r._isRaftSurfReac()), [], [])
    for obj, loc in (RHSObjs_noves - LHSObjs_noves) - (RHSObjs - LHSObjs):
        # If the object is only present on the RHS of conventional reactions, but is present on the LHS
        # of vesicle or raft reactions, display different warnings
        if isinstance(obj, nmodel.ComplexState):
            errors.append(
                f'Complex state {obj} in {loc} is only ever present on the RHS of conventional '
                f'reactions. It is present on the LHS of vesicle or raft reactions, but rafts or '
                f'vesicles of the corresponding types might not be added to the simulation.'
            )
        else:
            warnings.append(
                f'{obj} in {loc} is only ever present on the RHS of conventional reactions. It is '
                f'present on the LHS of vesicle or raft reactions, but rafts or vesicles of the '
                f'corresponding types might not be added to the simulation.'
            )

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

    if isinstance(sim.geom, ngeom._BaseTetMesh):
        # Check that the tetrahedrons are not too small
        smallTetRatio = sum(1 for tet in sim.geom.tets if tet.Vol <
                            MESH_SMALL_TET_VOL) / len(sim.geom.tets)
        if smallTetRatio > MESH_SMALL_TET_RATIO:
            warnings.append(
                f'{smallTetRatio*100:.2f}% of the tetrahedrons in the mesh are below a size of 20nm, '
                f'check that the scaling of the mesh is correct and, if so, consider using a coarser '
                f'mesh.'
            )

        # Check that the sizes and diffusion constants of vesicles are coherent with the size of the mesh
        meshDiag = numpy.linalg.norm(sim.geom.bbox.max - sim.geom.bbox.min)
        smallestDim = numpy.min(sim.geom.bbox.max - sim.geom.bbox.min)
        if sim._solverStr == 'TetVesicle':
            vesDt = sim.solver.getVesicleDT()

            # Check that vesicle dt is not too high for second-order vesicle reactions with immobile reactants
            # First get mobile species
            vsys2SpecDcsts = {}
            allSpecDcsts = {}
            for vsys in sim.model.ALL(nmodel.VolumeSystem):
                specDcsts = vsys2SpecDcsts.setdefault(vsys, {})
                for diff in vsys.ALL(nmodel.Diffusion):
                    if isinstance(diff.Dcst, numbers.Number):
                        specDcsts.setdefault(diff._elem, 0)
                        specDcsts[diff._elem] += diff.Dcst
                    elif isinstance(diff.Dcst, nmodel.CompDepRate):
                        cmplx = [diff._elem] if isinstance(
                            diff._elem, nmodel.ComplexState) else diff._elem
                        for state in cmplx:
                            specDcsts.setdefault(state, 0)
                            specDcsts[state] += diff.Dcst(state)
                for spec, D in specDcsts.items():
                    allSpecDcsts.setdefault(spec, 0)
                    allSpecDcsts[spec] += D

            mobileSpecs = set(
                spec for spec, dcst in allSpecDcsts.items() if dcst > 0)
            vesSurfSyss2reacs = {}
            vesSurfSyss2specs = {}
            for vssys in sim.model.ALL(nmodel.VesicleSurfaceSystem):
                for reac in vssys.ALL(nmodel.Reaction):
                    for spec in reac._getAllElems(nmodel.Location.OUT, sides='L'):
                        if not isinstance(spec, nmodel.Complex) and spec not in mobileSpecs:
                            vesSurfSyss2reacs.setdefault(vssys, set()).add(reac)
                            vesSurfSyss2specs.setdefault(vssys, set()).add(spec)
            tetRadius = (sim.geom.Vol / len(sim.geom.tets)
                         * 3 / 4 / numpy.pi) ** (1 / 3)
            for vssys, reacs in vesSurfSyss2reacs.items():
                for ves in vssys.locations:
                    if isinstance(ves, nmodel.Vesicle):
                        vesDtThresh = tetRadius ** 2 / (15 * ves.Dcst)
                        if vesDtThresh < vesDt:
                            reacsStr = '\n'.join(map(str, reacs))
                            warnings.append(
                                f'Vesicle {ves} is involved in vesicle surface reactions that involve '
                                f'immobile volume reactants {vesSurfSyss2specs[vssys]}. With the current '
                                f'vesicle Dt ({vesDt} s), the rate of these reactions would be '
                                f'underestimated. It is recommended to set the vesicle Dt to a value '
                                f'below {vesDtThresh:.3g} s. The affected reactions are:\n{reacsStr}'
                            )

            # Check that vesicle surface diffusion is not too high
            message = (
                '{} has a vesicle surface diffusion constant ({} m^2.s^-1, declared in {}) '
                'that, with the current vesicle Dt ({} s), would result in being randomly '
                'positioned on the vesicle surface after each vesicle step. If its value is '
                'correct, the diffusion constant could be reduced to lower computational load.'
            )
            for vssys in sim.model.ALL(nmodel.VesicleSurfaceSystem):
                for diff in vssys.ALL(nmodel.Diffusion):
                    for ves in vssys.locations:
                        if isinstance(ves, nmodel.Vesicle):
                            avg_jump = numpy.sqrt(4 * diff.Dcst * vesDt)
                            if avg_jump > numpy.pi * ves.Diameter:
                                warnings.append(message.format(
                                    diff._elem, diff.Dcst, diff, vesDt))
            for vesbind in sim.model.ALL(nmodel.VesicleBind):
                for ves, linkspec in zip(vesbind.Vesicles, vesbind._linkProducts):
                    if numpy.sqrt(4 * linkspec.Dcst * vesDt) > numpy.pi * ves.Diameter:
                        warnings.append(message.format(
                            linkspec, linkspec.Dcst, vesbind, vesDt))

            # Check that the sizes of rafts are coherent with the size of the mesh
            for raft in sim.model.ALL(nmodel.Raft):
                if raft.Diameter > meshDiag:
                    errors.append(
                        f'Raft {raft} has a diameter ({raft.Diameter} m) higher than the diagonal of '
                        f'the mesh bounding box ({meshDiag} m).'
                    )

            if len(list(sim.model.ALL(nmodel.Vesicle))) > 0:
                smallestVesDiam = min(ves.Diameter for ves in sim.model.ALL(nmodel.Vesicle))
                for ves in sim.model.ALL(nmodel.Vesicle):
                    if ves.Diameter > smallestDim:
                        errors.append(
                            f'Vesicle {ves} has a diameter ({ves.Diameter} m) higher than the smallest '
                            f'side of the mesh bounding box ({smallestDim} m).'
                        )
                    threshDist = smallestVesDiam + ves.Diameter
                    frac = 1 - chi2.cdf(3 * threshDist ** 2 /
                                        (6 * ves.Dcst * vesDt), 3)
                    if frac > ACCEPTABLE_JUMP_FRACTION:
                        recDt = (3 / chi2.ppf(1 - ACCEPTABLE_JUMP_FRACTION, 3)
                                * threshDist ** 2) / (6 * ves.Dcst)
                        warnings.append(
                            f'Vesicle {ves} has a diffusion constant ({ves.Dcst} m^2.s^-1) for which, with '
                            f'the current vesicle Dt ({vesDt} s), {frac * 100:.2f}% of the vesicle jumps could '
                            f'skip over the smallest vesicle (diameter {smallestVesDiam} m). This can happen '
                            f'for less than {ACCEPTABLE_JUMP_FRACTION * 100} percent of the jumps '
                            f'if the vesicle Dt is set to {recDt:.4e} s.'
                        )

                # Check whether some species could stay in tetrahedrons even if fully overlapped by vesicles
                # This can only happen if the diameter of the biggest vesicle is bigger than the diameter of
                # the circumsphere of the tetrahedron
                # We compute the volume of the regular tetrahedron that fits in the biggest vesicle
                maxVesDiam = max(ves.Diameter for ves in sim.model.ALL(nmodel.Vesicle))
                volThresh = maxVesDiam ** 3 / (9 * numpy.sqrt(3))
                overlap_ratio = sum(1 for tet in sim.geom.tets if tet.Vol < volThresh) / len(sim.geom.tets)
                if overlap_ratio > VESICLE_TET_FULL_OVERLAP_RATIO:
                    for comp in sim.geom.ALL(ngeom.Compartment):
                        allElems = set(s.getID() for s in comp.stepsComp.getAllSpecs(sim.model.stepsModel))
                        diffElems = set()
                        for vsys in comp.systems:
                            diffElems |= set(
                                s._getStepsObjects()[0].getID() for s in vsys2SpecDcsts.get(vsys, {})
                            )
                        noDiffElems = allElems - diffElems
                        warnings.append(
                            f'The following elements do not diffuse in compartment {comp.name} and could '
                            f'thus stay present in tetrahedrons that are fully overlapped by a vesicle: '
                            f'{", ".join(map(str, noDiffElems))}. Note that depending on the shapes of'
                            f'tetrahedrons, this might not happen for this specific mesh.'
                        )

                # Check whether non-diffusing species can react in second order or higher reactions
                for vsys in sim.model.ALL(nmodel.VolumeSystem):
                    for reac in vsys.ALL(nmodel.Reaction):
                        specifiers = [nmodel.Reaction._FwdSpecifier]
                        if reac._bidir:
                            specifiers.append(nmodel.Reaction._BkwSpecifier)
                        for specif in specifiers:
                            for lhs, rhs, rm in reac[specif]._LRP:
                                nonDiffElems = set(e._elem for e in lhs) - set(vsys2SpecDcsts.get(vsys, {}))
                                if len(nonDiffElems) > 1:
                                    warnings.append(
                                        f'Reaction {reac} involves at least two reactants that are '
                                        f'non-diffusing ({", ".join(map(str, nonDiffElems))}). In presence '
                                        f'of vesicles, the reaction rate will be overestimated in the '
                                        f'parts of the mesh that are overlapped by vesicles. It is '
                                        f'advised to either add diffusion rules for these reactants '
                                        f'or to not use vesicles in compartments {vsys.locations}.'
                                    )

    ok = len(errors) == 0
    if printmsgs:
        nutils._print('Model checking:', 1, tpe=nutils.MessagesTypes.HEADER)
        if len(errors) == len(warnings) == 0:
            nutils._print('No errors were found', 1,
                          tpe=nutils.MessagesTypes.SUCCESS)
        else:
            nutils._print(
                'Errors and warnings can be turned off by providing check=False when creating the '
                'Simulation. Checks can be manually run with steps.simcheck.Check(sim, True).', 1
            )
        if len(errors) > 0:
            nutils._print(
                'Errors:', 1, tpe=nutils.MessagesTypes.ERROR + nutils.MessagesTypes.BOLD)
        for e in errors:
            nutils._print(e, 1, tpe=nutils.MessagesTypes.ERROR, indent=1)
        if len(warnings) > 0:
            nutils._print(
                'Warnings:', 1, tpe=nutils.MessagesTypes.WARNING + nutils.MessagesTypes.BOLD)
        for w in warnings:
            nutils._print(w, 1, tpe=nutils.MessagesTypes.WARNING, indent=1)
        return ok
    else:
        return ok, errors, warnings
