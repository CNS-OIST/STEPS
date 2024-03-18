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

from steps import stepslib
from enum import Enum

# -----------------------------------------------------
# Steps classes are named here, and possibly extended
# -----------------------------------------------------


class Model(stepslib._py_Model)     :
    """
    Top-level container for the objects in a kinetic model.
    """
    pass


class Volsys(stepslib._py_Volsys)   :
    """
    A container that groups reactions and diffusion rules involving Species
    located in a volume.
    """
    pass


class Surfsys(stepslib._py_Surfsys) :
    """
    A container that groups reactions, diffusion rules and other transport
    mechanisms involving a Species located in a membrane.
    """
    pass


class Spec(stepslib._py_Spec)       :
    """
    A chemical species which can be a reactant and/or product in reaction
    stoichiometry and/or associated with a diffusion rule, or other transport
    mechanisms.
    """
    pass


class Reac(stepslib._py_Reac)       :
    """
    A kinetic reaction rule in a volume.

    The reaction rule is specified by:

    * Species on the left hand side of the reaction: the reactants (lhs).
    * Species on the right hand side of the reaction: the products (rhs).
    * Rate constant for the reaction, supplied in s.i. units (kcst).
    """
    pass


class SReac(stepslib._py_SReac)     :
    """
    A reaction rule where at least one reactant or product is embedded in a surface.

    In a surface reaction, the species can be classified as:

    * Reactants on the left hand side of the reaction:
        * Species in the 'outer' compartment (olhs)
        * or Species in the 'inner' compartment (ilhs)
        * Species on the surface (slhs)
    * Products on right hand side of the reaction:
        * Species in the 'outer' compartment (orhs)
        * Species in the 'inner' compartment (irhs)
        * Species on the surface (srhs)
    * The reaction rate is defined by kcst, supplied in s.i. units.
    """
    pass


class Diff(stepslib._py_Diff)       :
    """
    A diffusion rule for a chemical species in a volume or on a surface.

    A diffusion rule is described by:

    * Species to which the diffusion rule applies (lig).
    * Diffusion constant (dcst) specified in s.i. units.
    """
    pass


class Chan(stepslib._py_Chan)       :
    """
    A channel object which groups a set of channel states. Whilst not involved
    directly in reactions, currents, transitions etc it provides necessary
    grouping for channel states.
    """
    pass


class ChanState(stepslib._py_ChanState):
    """
    A channel state object which may be involved in channel events such as
    currents, voltage-dependent transitions and any other steps.model.Spec events,
    such as surface reactions and diffusion.
    """
    pass


class OhmicCurr(stepslib._py_OhmicCurr):
    """
    An ohmic current object that describes an ohmic current through a channel in
    a particular conducting state.

    An ohmic current object is a simple current based on a single, fixed value
    for single-channel conductance and a constant reversal potential. An ohmic current
    does not result in movement of ions between compartments, but simply
    contributes a continous current to the EField solver for every conducting channel
    on a membrane surface that contains the ohmic current.
    """
    pass


class GHKcurr(stepslib._py_GHKcurr) :
    """
    A current object based on the Goldman-Hodgkin-Katz flux equation,
    that describes a current through a channel in a particular permeable state.

    Each GHK current in the simulation is solved within the SSA with a rate determined
    from the simulation state, i.e. membrane potential, 'outer' and 'inner'
    concentration of the ion and temperature (which is fixed), and constant parameters
    ion valence and permeability.
    For a permeabilty to be found the user must supply
    information from a channel measurement:

    * the single-channel conductance
    * the potential
    * the temperature (this may be different from the STEPS simulation temperature)
    * the 'outer' concentration of the ion
    * the 'inner' concentration of the ion

    This information must be set with the member function setPInfo()

    A GHK current involves, optionally, a real transfer of ions between compartments
    separated by a membrane in the STEPS simulation. It is important that
    these ions implement diffusion if compartments are not well-mixed.
    """
    pass


class VDepSReac(stepslib._py_VDepSReac):
    """
    A voltage-dependent reaction involving any species or channel state.

    A voltage-dependent reaction is similar to a surface reaction, except the
    reaction parameter is voltage-dependent.

    A function object must be passed to the constructor returning the reaction
    'constant' as a function of voltage. This function is used to fill a table of
    reaction 'constants' at a range of voltage values (default: -150mV to 100mV,
    in steps of 0.1mV) These default values may be altered with optional argument
    'vrange'.
    """

    def __init__(self, id, surfsys, **kwargs):
        """
        Construction::

            vdepsreac = steps.model.VDepSReac(id, surfsys,
            ilhs = [ ], olhs = [ ], slhs = [ ],
            irhs = [ ], orhs = [ ], srhs = [ ],
            k = <function>,
            vrange = [-150.0e-3, 100.0e-3, 1.0e-4] )

        Construct a voltage-dependent reaction object with identifier string id
        and assign surfsys as the parent surface system. A list of
        left hand reactants are assigned with ilhs, olhs and slhs
        (default for each is an empty list). A list of right hand side
        products are assigned with irhs, orhs and srhs (default for each
        is an empty list).
        A function that returns the kinetic reaction 'constant' in ordinary, Molar
        units at any voltage (in volts) is supplied with argument k.
        A 'voltage range' over which to calculate the reaction rate is
        optionally provided by the argument vrange (in Volts) as a list:
        [minimum voltage, maximum voltage, voltage step]

        Arguments:
            * string id
            * steps.model.Surfsys surfsys
            * list(steps.model.Spec) ilhs (default = [ ])
            * list(steps.model.Spec) olhs (default = [ ])
            * list(steps.model.Spec) slhs (default = [ ])
            * list(steps.model.Spec) irhs (default = [ ])
            * list(steps.model.Spec) orhs (default = [ ])
            * list(steps.model.Spec) srhs (default = [ ])
            * function k
            * list vrange (default = [-150.0e-3, 100.0e-3, 1.0e-4])
        """
        # Get optional arguments that must be passed without change
        all_args = {hs: val for hs, val in kwargs.items() if hs in {'olhs', 'ilhs',
                                                                    'slhs', 'irhs', 'srhs', 'orhs'}}
        # Process vrange
        rate_f = kwargs['k']
        vtable = _VoltageTable(rate_f, kwargs.get('vrange'))
        # Construct
        all_args.update(vtable.as_dict())
        super(self.__class__, self).__init__(id, surfsys, **all_args)

class Vesicle(stepslib._py_Vesicle):
    """
    A vesicle object, described by a spherical diameter and diffusion rate,
    models the many behaviours and interactions of these complex biological
    entities. A vesicle may be formed by endocytosis or by direct user input,
    may exit by exocytosis or user input, may transport species on its surface or
    luminally, and may undergo many interactions with its environment. These
    interactions are individually described in each corresponding class.
    """
    pass

class Raft(stepslib._py_Raft):
    """
    The membrane analogy of vesicles, rafts effectively group surface species
    within a defined radius and exist in and may diffuse within patches. Rafts
    and the species they contain undergo special interactions such as raft
    endocytosis, raft generation and raft dissolution.
    """
    pass

class LinkSpec(stepslib._py_LinkSpec):
    """
    A link species is a special kind of chemical species that links two vesicles
    (by forming a bond with another link species).
    A link species is formed by a vesicle binding event and exists within a
    specified upper and lower bound of length for the duration of its existence.
    A link species may diffuse on a vesicle surface, but only within this length
    bound. Link species are destroyed by a vesicle unbinding event.
    """
    pass

class VesSurfsys(stepslib._py_VesSurfsys):
    """
    A container that groups certain vesicle interactions such as VesSDiff,
    VesSReac and Exocytosis.
    """
    pass

class VesSReac(stepslib._py_VesSReac):
    """
    A reaction rule involving reactants and products in a vesicle surface, or
    products inside the vesicle.

    In a vesicle surface reaction, the species can be classified as:

    * Reactants on the left hand side of the reaction:
        * Species in the compartment (olhs)
        * Species on a patch (slhs)
        * Species on the vesicle surface (vlhs)
        * Link species on the vesicle surface (llhs)

    * Products on the right hand side of the reaction:
        * Species inside the vesicle lumen (irhs)
        * Species in the compartment (orhs)
        * Species on a patch (srhs)
        * Species on the vesicle surface (vrhs)
        * Link species on the vesicle surface (lrhs)

    * Dependencies for the reaction:
        * Species on the vesicle surface (vdeps)

    * The reaction rate is defined by kcst, supplied in s.i. units.

    * Optional immobility may be defined, which must be 1 (immobilizing), 0 (no
          effect), or -1 (mobilizing).

    * Optional max_distance may be defined (units m). The reaction will not occur
        with reactants on the vesicle surface that are beyond this distance to
        patch triangle barycenters. This is intended for modeling docking reactions
        at realistic distances.

    """
    pass

class VesSDiff(stepslib._py_VesSDiff):
    """
    A diffusion rule for a chemical species on a vesicle surface.

    A vesicle surface diffusion rule is described by:

    * Species to which the diffusion rule applies (lig).
    * Diffusion constant (dcst) specified in s.i. units.
    """
    pass

class VesBind(stepslib._py_VesBind):
    """
    A vesicle binding event binds two vesicles by creating link species between them.
    The types of vesicle are arbitrarily termed 1 and 2, and may be the same.
    A reactant on vesicle surface of vesicle 1 binds with a reactant on vesicle
    2, and link species are created, one on vesicle 1 and one on vesicle 2.
    Internally these two link species will be associated with each other and
    exist within a defined length bound.

    * The reactants for the vesicle binding reaction are described as:
        * The reactant species (r1) on the surface of vesicle 1 (ves1)
        * The reactant species (r2) on the surface of vesicle 2 (ves2)

    * The products of the vesicle binding reaction are described as:
        * Link species on vesicle 1 surface (l1)
        * Link species on vesicle 2 surface (l2)

    * The optional dependencies are:
        * Species on vesicle 1 surface (vdeps1)
        * Species on vesicle 2 surface (vdeps2)
        * Link species on vesicle 1 surface (ldeps1)
        * Link species on vesicle 2 surface (ldeps2)

    * The rate constant is defined by kcst, supplied in s.i. units of /(M.s)

    * Optional immobility may be defined, which must be 1 (immobilizing), 0 (no
          effect), or -1 (mobilizing).
    """
    pass

class VesUnbind(stepslib._py_VesUnbind):
    """
    A vesicle unbinding event unbinds two vesicles bound by two link species.
    The types of vesicle are arbitrarily termed 1 and 2, and may be the same.
    Upon application of this reaction, the link species on vesicle 1 becomes a
    species on the vesicle 1 surface and the link species on vesicle 2 becomes a
    species on the vesicle 2 surface.

    * The two link species and the complex formed between then are effectively
        the reactant for this reaction and are described as:
        * The link species on the vesicle 1 surface (l1)
        * The link species on the vesicle 2 surface (l2)

    * The products of the vesicle unbinding reaction are described as:
        * Species (p1) on the surface of vesicle 1 (ves1)
        * Species (p2) on the surface of vesicle 2 (ves2)

    * The rate constant is defined by kcst, supplied in s.i. units of /s.

    * Optional immobility may be defined, which must be 1 (immobilizing), 0 (no
          effect), or -1 (mobilizing).
    """
    pass

class Exocytosis(stepslib._py_Exocytosis):
    """
    An exocytosis event models the process of vesicle exocytosis. By default the vesicle
    is destroyed, species in the vesicle surface are deposited in the patch at
    the location at which exocytosis occurs, and species inside the vesicle lumen
    are deposited in the opposite compartment.

    This behavior can be modified by optional arguments to the exocytosis
    constructor:

    * raft : if this argument is given the exocytosis event will create a raft
    on the patch at the location that exocytosis occurs. Any species on the
    vesicle surface will be deposited in the raft instead of in the patch.

    * kiss_and_run : models a kiss-and-run exocytosis event. This does not result
    in collapse of the vesicle, instead the vesicle is maintained in position
    after the release of vesicle lumen contents into the opposite compartment.
    The vesicle surface species are maintained on the vesicle.

    * kiss_and_run_spec_changes : to aid kiss-and-run modeling, this argument can
    be used to specify any species changes that take place on the vesicle surface
    upon a kiss-and-run exocytosis event. This can be useful for modeling
    maturation of a complex for which undocking may be dependent.

    * kiss_and_run_partial_release : to aid kiss-and-run modeling, this argument can
    be used to specify partial release of any inner species that are released
    upon a kiss-and-run exocytosis event. 

    In addition there are the following optional arguments:

    * deps is a list of species dependencies on the vesicle surface.
    * kcst is the rate constant, supplied in s.i. units of /s.
    """
    pass

class Endocytosis(stepslib._py_Endocytosis):
    """
    An endocytosis event models the process of vesicle endocytosis by creating a
    new vesicle within a compartment. A vesicle will be created at a given rate
    when an optional species signature is met within an endocytic zone of a
    patch that is defined by solver function calls. The species signature is
    defined by argument deps and the rate constant kcst is supplied in s.i. units
    of /s.

    The vesicle will be created in either the 'inner' or 'outer' compartment to
    the patch. This is specified by defining which type of vesicle will be
    created by argument irhs or orhs respectively, only one of which must be
    defined.

    """
    pass

class Raftsys(stepslib._py_Raftsys):
    """
    A container that groups certain raft interactions such as RaftEndocytosis,
    RaftSReac and RaftDis.
    """
    pass

class RaftSReac(stepslib._py_RaftSReac):
    """
    A reaction rule involving reactants and products on a raft surface.

    In a raft surface reaction, the species can be classified as:

    * Reactants on the left hand side of the reaction:
        * Species in the outer compartment (olhs)
        * Species in the inner compartment (ilhs)
        * Species on a patch (slhs)
        * Species on the raft surface (rlhs)

    * Products on the right hand side of the reaction:
        * Species in the outer compartment (orhs)
        * Species in the inner compartment (irhs)
        * Species on a patch (srhs)
        * Species on the raft surface (rrhs)

    * Dependencies for the reaction:
        * Species on the raft surface (rdeps)

    * Anti-dependencies for the reaction (if these species are present the
        reaction will not happen):
        * Species on the raft surface (anti_rdeps)

    * The reaction rate is defined by kcst, supplied in s.i. units.

    * Optional immobility may be defined, which must be 1 (immobilizing), 0 (no
          effect), or -1 (mobilizing).
    """
    pass

class RaftEndocytosis(stepslib._py_RaftEndocytosis):
    """
    A raft endocytosis event models the process of vesicle endocytosis by
    creating a new vesicle within a compartment from a raft. A vesicle will be
    created at a given rate when an optional species signature is met within the
    raft. The argument deps defines the species signature in the raft, and kcst
    is the rate constant in s.i. units of /s.

    The vesicle will be created in either the 'inner' or 'outer' compartment to
    the patch in which the raft resides. This is specified by defining which
    type of vesicle will be created by argument irhs or orhs respectively, only
    one of which must be defined.

    """
    pass

class RaftGen(stepslib._py_RaftGen):
    """
    Generate a raft at a given rate when a defined species signature is met
    within a patch triangle. The argument deps defines the species signature,
    raft is the type of raft that is created and kcst is the rate constant
    in s.i. units of /s.

    """
    pass

class RaftDis(stepslib._py_RaftDis):
    """
    A raft dissolution event results in removal of a raft at a given rate when
    the population of species within a raft are at or below a given signature.
    Any remaining species in the raft are inserted into patch triangles. The
    arguments antideps defines the species signature and kcst is the rate constant
    in s.i. units of /s.
    """
    pass


# -----------------------------------------------------
# Aux: Voltage table creator
# -----------------------------------------------------
class _VoltageTable:
    # These parameters are the hard-coded default voltage bounds
    vmin = -150.0e-3
    vmax = 100.0e-3
    dv = 1.0e-4

    def __init__(self, k_func, vrange=None):
        self.k_func = k_func
        if vrange:
            assert len(vrange) == 3
            self.vmin = vrange[0]
            self.vmax = vrange[1]
            self.dv = vrange[2]

        self.tablesize = int((self.vmax - self.vmin) / self.dv) + 1
        klist = [0.0] * self.tablesize

        #Calc table
        v = self.vmin
        for i in range(self.tablesize):
            klist[i] = k_func(v)
            v += self.dv

        #save klist as ktab
        self.ktab = klist

    #----------------
    def as_dict(self):
        """Return all important properties as a dictionary"""
        indexes = ['vmin', 'vmax', 'dv', 'ktab', 'tablesize']
        return {idx: getattr(self, idx)
                for idx in indexes}


# -----------------------------------------------------
# Aux: Immobilization enum
# -----------------------------------------------------

try:
    Immobilization = stepslib._py_Immobilization
except AttributeError:
    # Older versions of cython do not export this enum
    class Immobilization(Enum):
        IMMOBILIZING = 0
        MOBILIZING   = 1
        NO_EFFECT    = 2

# Aliases
IMMOBILIZING = Immobilization.IMMOBILIZING
MOBILIZING = Immobilization.MOBILIZING
NO_EFFECT = Immobilization.NO_EFFECT
