# cython:language_level=3str
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

include "cysteps.pyx"

# ======================================================================================================================
# Python bindings to namespace steps::mpi
# ======================================================================================================================
cimport steps_mpi
from steps_mpi cimport TetOpSplitP
from steps_mpi cimport TetVesicleRDEF
from steps_mpi cimport TetVesicleVesRaft

cdef extern from "solver/fwd.hpp":
    cdef index_t _UNDEFINED_VESICLE "steps::solver::vesicle_individual_id::unknown_value()"
    cdef index_t _UNDEFINED_RAFT "steps::solver::raft_individual_id::unknown_value()"

UNDEFINED_VESICLE = _UNDEFINED_VESICLE
UNDEFINED_RAFT = _UNDEFINED_RAFT

def mpiInit():
    """
    Initialise the MPI solver. NOTE: handled automatically, should not be called by user.
    """
    steps_mpi.mpiInit()

def getRank():
    """
    Get the MPI rank. NOTE: handled automatically, should not be called by user.
    """
    return steps_mpi.getRank()

def getNHosts():
    """
    Get the number of hosts. NOTE: handled automatically, should not be called by user.
    """
    return steps_mpi.getNHosts()

def mpiFinish():
    """
    Finalise the MPI solver. NOTE: handled automatically, should not be called by user.
    """
    steps_mpi.mpiFinish()

def mpiAbort():
    """
    Aborts the MPI solver. NOTE: handled automatically, should not be called by user.
    """
    steps_mpi.mpiAbort()

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TetOpSplitP(_py_TetAPI):
    """Bindings for MPI TetOpSplitP"""
# ----------------------------------------------------------------------------------------------------------------------
    cdef TetOpSplitP *ptrx(self):
        return <TetOpSplitP*> self._ptr

    def __init__(self, _py_Model model, _py_Geom geom, _py_RNG rng, int calcMembPot=0, std.vector[int] tet_hosts = [], dict tri_hosts = {}, std.vector[int] wm_hosts = []):
        """
        Construction::

            sim = steps.solver.TetOpSplit(model, geom, rng, tet_hosts=[], tri_hosts={}, wm_hosts=[], calcMembPot=0)

        Create a spatial stochastic solver based on operator splitting: reaction events are partitioned and diffusion is approximated.
        If voltage is to be simulated, argument calcMembPot specifies the solver. E.g. calcMembPot=steps.solver.EF_DV_PETSC will utilise the PETSc library. calcMembPot=0 means that voltage will not be simulated.

        Arguments:
        steps.model.Model model
        steps.geom.Geom geom
        steps.rng.RNG rng
        list<int> tet_hosts (default=[])
        dict<index_t, int> tri_hosts (default={})
        list<int> wm_hosts (default=[])
        int calcMemPot (default=0)

        """
        cdef std.map[steps.triangle_global_id, int] _tri_hosts
        for key, elem in tri_hosts.items():
            _tri_hosts[steps.triangle_global_id(key)] = elem
        # We constructed a map. Now call constructor
        if model == None:
            raise TypeError('The Model object is empty.')
        if geom == None:
            raise TypeError('The Geom object is empty.')
        if rng == None:
            raise TypeError('The RNG object is empty.')
        self._ptr = new TetOpSplitP(model.ptr(), geom.ptr(), rng.ptr(), calcMembPot, tet_hosts, _tri_hosts, wm_hosts)

    def getSolverName(self, ):
        """
        Returns a string of the solver's name.

        Syntax::

            getSolverName()

        Arguments:
        None

        Return:
        string
        """
        return self.ptrx().getSolverName()

    def getSolverDesc(self, ):
        """
        Returns a string giving a short description of the solver.

        Syntax::

            getSolverDesc()

        Arguments:
        None

        Return:
        string
        """
        return self.ptrx().getSolverDesc()

    def getSolverAuthors(self, ):
        """
        Returns a string of the solver authors names.

        Syntax::

            getSolverAuthors()

        Arguments:
        None

        Return:
        string
        """
        return self.ptrx().getSolverAuthors()

    def getSolverEmail(self, ):
        """
        Returns a string giving the author's email address.

        Syntax::

            getSolverEmail()

        Arguments:
        None

        Return:
        string

        """
        return self.ptrx().getSolverEmail()

    def reset(self, ):
        """
        Reset the simulation to the state the solver was initialised to.

        Syntax::

            reset()

        Arguments:
        None

        Return:
        None
        """
        self.ptrx().reset()

    def run(self, double endtime):
        """
        Advance the simulation until endtime (given in seconds) is reached.
        The endtime must be larger or equal to the current simulation time.

        Syntax::

            run(endtime)

        Arguments:
        float endtime

        Return:
        None
        """
        self.ptrx().run(endtime)

    def advance(self, double adv):
        """
        Advance the simulation for adv seconds.

        Syntax::

            advance(adv)

        Arguments:
        float adv

        Return:
        None

        """
        self.ptrx().advance(adv)

    def step(self, ):
        """
        Advance the simulation for one 'step'. In stochastic solvers this is one
        'realization' of the Gillespie SSA (one reaction 'event').
        In numerical solvers (currently Wmrk4) this is one time-step, with the
        stepsize defined with the setDT method.

        Syntax::

            step()

        Arguments:
        None

        Return:
        None

        """
        self.ptrx().step()

    def checkpoint(self, str file_name):
        """
        Checkpoint data to a file.

        Syntax::

            checkpoint(file_name)

        Arguments:
        string file_name

        Return:
        None

        """
        fullName = file_name + str(f'_{getRank()}')
        self.ptrx().checkpoint(to_std_string(fullName))

    def restore(self, str file_name):
        """
        Restore data from a file.

        Syntax::

            restore(file_name)

        Arguments:
        string file_name

        Return:
        None
        """
        fullName = file_name + str(f'_{getRank()}')
        self.ptrx().restore(to_std_string(fullName))

    def setEfieldDT(self, double efdt):
        """
        Set the stepsize for membrane potential solver (default 1us).
        This is the time for each voltage calculation step. The SSA will
        run until passing this stepsize, so in fact each membrane potential
        time step will vary slightly around the dt so as to be aligned with the SSA.

        Syntax::

            setEFieldDT(dt)

        Arguments:
        float dt

        Return:
        None

        """
        self.ptrx().setEfieldDT(efdt)

    def getEfieldDT(self, ):
        """
        Get the stepsize for the membrane potential solver.

        Syntax::

            getEFieldDT(dt)

        Arguments:
        None

        Return:
        float

        """
        return self.ptrx().getEfieldDT()

    def setTemp(self, double t):
        """
        Set the simulation temperature. Currently, this will only
        influence the GHK flux rate, so will only influence simulations
        including membrane potential calculation.

        Syntax::

        	setTemp(temp)

        Arguments:
        float temp

        Return:
        None

        """
        self.ptrx().setTemp(t)

    def getTemp(self, ):
        """
        Return the simulation temperature.

        Syntax::

        	getTemp()

        Arguments:
        None

        Return:
        float

        """
        return self.ptrx().getTemp()

    def saveMembOpt(self, str opt_file_name):
        """
        Save the membrane optimisation.

        Syntax::

            saveMembOpt(opt_file_name)

        Arguments:
        string opt_file_name

        Return:
        None
        """
        self.ptrx().saveMembOpt(to_std_string(opt_file_name))

    def getTime(self, ):
        """
        Returns the current simulation time in seconds.

        Syntax::

            getTime()

        Arguments:
        None

        Return:
        float
        """
        return self.ptrx().getTime()

    def getA0(self, ):
        """
        Returns the total propensity of the current simulation state
        (the total propensity multiplied by an infinitesimally small
        time dt gives the probability that a reaction will occur in that dt).
        For Tetexact this includes the propensity from the extension of the SSA
        for diffusive flux between tetrahedral elements in the mesh.

        Syntax::

            getA0()

        Arguments:
        None

        Return:
        float

        """
        return self.ptrx().getA0()

    def getNSteps(self, ):
        """
        Return the number of 'realizations' of the SSA, the number of reaction
        (and diffusion) events in stochastic solvers.

        Syntax::

            getNSteps()

        Arguments:
        None

        Return:
        int

        """
        return self.ptrx().getNSteps()

    def setTime(self, double time):
        """
        Set the current simulation time.

        Syntax::

            setTime(time)

        Arguments:
        float time

        Return:
        None

        """
        self.ptrx().setTime(time)

    def setNSteps(self, uint nsteps):
        """
        Set the number of 'realizations' of the SSA, the number of reaction
        (and diffusion) events in stochastic solvers.

        Syntax::

            setNSteps(nsteps)

        Arguments:
        uint nsteps

        Return:
        None

        """
        self.ptrx().setNSteps(nsteps)


    def getBatchTetSpecCounts(self, std.vector[index_t] tets, str s):
        """
        Get the counts of a species s in a list of tetrahedrons.

        Syntax::

            getBatchTetSpecCounts(tets, s)

        Arguments:
        list<index_t> tets
        string s

        Return:
        list<double>

        """
        return self.ptrx().getBatchTetSpecCounts(tets, to_std_string(s))

    def getBatchTriSpecCounts(self, std.vector[index_t] tris, str s):
        """
        Get the counts of a species s in a list of triangles.

        Syntax::

            getBatchTriSpecCounts(tris, s)

        Arguments:
        list<index_t> tris
        string s

        Return:
        list<double>

        """
        return self.ptrx().getBatchTriSpecCounts(tris, to_std_string(s))

    def setBatchTetSpecConcs(self, std.vector[index_t] tets, str s, std.vector[double] concs):
        """
        Set the concentration of a species s in a list of tetrahedrons individually.

        Syntax::

            setBatchTetConcs(tets, s, concs)

        Arguments:
        list<index_t> tets
        string s
        list<double> concs

        Return:
        None

        """
        self.ptrx().setBatchTetSpecConcs(tets, to_std_string(s), concs)

    def getBatchTetSpecConcs(self, std.vector[index_t] tets, str s):
        """
        DEPRECATED
        Get the individual concentration of a species s in a list of tetrahedrons.

        Syntax::

            getBatchTetSpecConcs(tets, s)

        Arguments:
        list<index_t> tets
        string s

        Return:
        list<double>

        """
        return self.ptrx().getBatchTetSpecConcs(tets, to_std_string(s))

    # ---------------------------------------------------------------------------------
    # NUMPY section - we accept numpy arrays and generically typed memory-views
    # ---------------------------------------------------------------------------------
    def getBatchTetSpecCountsNP(self, index_t[:] index_array, str s, double[:] counts):
        """
        Get the counts of a species s in a list of tetrahedrons.

        Syntax::
            getBatchTetSpecCountsNP(indices, s, counts)

        Arguments:
        numpy.array<index_t> indices
        string s
        numpy.array<double, length = len(indices)> counts

        Return:
        None

        """
        self.ptrx().getBatchTetSpecCountsNP(&index_array[0], index_array.shape[0], to_std_string(s), &counts[0], counts.shape[0])

    def getBatchTetSpecConcsNP(self, index_t[:] index_array, str s, double[:] concs):
        """
        DEPRECATED
        Get the individual concentrations of a species s in a list of tetrahedrons.

        Syntax::
            getBatchTetSpecConcsNP(indices, s, concs)

        Arguments:
        numpy.array<index_t> indices
        string s
        numpy.array<double, length = len(indices)> concs

        Return:
        None

        """
        self.ptrx().getBatchTetSpecConcsNP(&index_array[0], index_array.shape[0], to_std_string(s), &concs[0], concs.shape[0])

    def getBatchTriSpecCountsNP(self, index_t[:] index_array, str s, double[:] counts):
        """
        Get the counts of a species s in a list of triangles.

        Syntax::
            getBatchTriSpecCountsNP(indices, s, counts)

        Arguments:
        numpy.array<index_t> indices
        string s
        numpy.array<double, length = len(indices)> counts

        Return:
            None

        """
        self.ptrx().getBatchTriSpecCountsNP(&index_array[0], index_array.shape[0], to_std_string(s), &counts[0], counts.shape[0])

    def setBatchTetSpecConcsNP(self, index_t[:] index_array, str s, double[:] concs):
        """
        Set the concetration of a species s in a list of tetrahedrons.

        Syntax::
            setBatchTetSpecConcsNP(indices, s, concs)

        Arguments:
        numpy.array<index_t> indices
        string s
        numpy.array<double, length = len(indices)>

        Return:
        None

        """
        self.ptrx().setBatchTetSpecConcsNP(&index_array[0], index_array.shape[0], to_std_string(s), &concs[0], concs.shape[0])

    def sumBatchTetCountsNP(self, index_t[:] tet_array, str s):
        """
        Return the accumulated sum of species s in a batch of tetrahedrons.

        This function requires NumPy array as input, and called globally in all processes.

        Syntax::

            sumBatchTetCountsNP(tet_array, s)

        Arguments:
        numpy.array<index_t> tet_array
        string s

        Return:
        float
        """
        return self.ptrx().sumBatchTetCountsNP(&tet_array[0], tet_array.shape[0], to_std_string(s))

    def sumBatchTriCountsNP(self, index_t[:] tri_array, str s):
        """
        Return the accumulated sum of species s in a batch of triangles.

        This function requires NumPy array as input, and called globally in all processes.

        Syntax::

            sumBatchTriCountsNP(tri_array, s)

        Arguments:
        numpy.array<index_t> tri_array
        string s

        Return:
        float
        """
        return self.ptrx().sumBatchTriCountsNP(&tri_array[0], tri_array.shape[0], to_std_string(s))

    def sumBatchTriGHKIsNP(self, index_t[:] tri_array, str ghk):
        """
        Return the accumulated sum of GHK currents in a batch of triangles.

        This function requires NumPy array as input, and called globally in all processes.

        Syntax::

            sumBatchTriGHKIsNP(tri_array, ghk)

        Arguments:
        numpy.array<index_t> tri_array
        string ghk

        Return:
        float
        """
        return self.ptrx().sumBatchTriGHKIsNP(&tri_array[0], tri_array.shape[0], to_std_string(ghk))

    def sumBatchTriOhmicIsNP(self, index_t[:] tri_array, str ghk):
        """
        Return the accumulated sum of Ohmic currents in a batch of triangles.

        This function requires NumPy array as input, and called globally in all processes.

        Syntax::

            sumBatchTriOhmicIsNP(tri_array, oc)

        Arguments:
        numpy.array<index_t> tri_array
        string oc

        Return:
        float
        """
        return self.ptrx().sumBatchTriOhmicIsNP(&tri_array[0], tri_array.shape[0], to_std_string(ghk))

    def getBatchTriOhmicIsNP(self, index_t[:] index_array, str oc, double[:] counts):
        """
        Get the Ohmic currents in a list of triangles.

        Syntax::
            getBatchTriOhmicIsNP(indices, oc, counts)

        Arguments:
        numpy.array<index_t> indices
        string oc
        numpy.array<double, length = len(indices)> counts

        Return:
            None

        """
        self.ptrx().getBatchTriOhmicIsNP(&index_array[0], index_array.shape[0], to_std_string(oc), &counts[0], counts.shape[0])

    def getBatchTriGHKIsNP(self, index_t[:] index_array, str ghk, double[:] counts):
        """
        Get the GHK currents in a list of triangles.

        Syntax::
            getBatchTriGHKIsNP(indices, ghk, counts)

        Arguments:
        numpy.array<index_t> indices
        string ghk
        numpy.array<double, length = len(indices)> counts

        Return:
            None

        """
        self.ptrx().getBatchTriGHKIsNP(&index_array[0], index_array.shape[0], to_std_string(ghk), &counts[0], counts.shape[0])

    def getBatchTriVsNP(self, index_t[:] index_array, double[:] counts):
        """
        Get the Voltages in a list of triangles.

        Syntax::
            getBatchTriVsNP(indices, counts)

        Arguments:
        numpy.array<index_t> indices
        numpy.array<double, length = len(indices)> counts

        Return:
            None

        """
        self.ptrx().getBatchTriVsNP(&index_array[0], index_array.shape[0], &counts[0], counts.shape[0])

    def getBatchTetVsNP(self, index_t[:] index_array, double[:] counts):
        """
        Get the Voltages in a list of tetrahedrons.

        Syntax::
            getBatchTetVsNP(indices, counts)

        Arguments:
        numpy.array<index_t> indices
        numpy.array<double, length = len(indices)> counts

        Return:
            None

        """
        self.ptrx().getBatchTetVsNP(&index_array[0], index_array.shape[0], &counts[0], counts.shape[0])

    def getBatchTriBatchOhmicIsNP(self, index_t[:] index_array, ocs, double[:] counts):
        """
        Get the values of a list of Ohmic currents in a list of triangles,
        store in a flatten 2d array.
        The value of current ocs[j] of triangle index_array[i] is stored in counts[i * len(ocs) + j]

        Syntax::
            getBatchTriBatchOhmicIsNP(indices, ocs, counts)

        Arguments:
        numpy.array<index_t> indices
        std.vector[string] ocs
        numpy.array<double, length = len(indices) * len(ocs)> counts

        Return:
            None

        """
        cdef std.vector[string] std_ocs = to_vec_std_strings(ocs)
        self.ptrx().getBatchTriBatchOhmicIsNP(&index_array[0], index_array.shape[0], std_ocs, &counts[0], counts.shape[0])

    def getBatchTriBatchGHKIsNP(self, index_t[:] index_array, ghks, double[:] counts):
        """
        Get the values of a list of GHK currents in a list of triangles,
        store in a flatten 2d array.
        The value of current ghks[j] of triangle index_array[i] is stored in counts[i * len(ghks) + j]

        Syntax::
            getBatchTriBatchGHKIsNP(indices, ghks, counts)

        Arguments:
        numpy.array<index_t> indices
        std.vector[string] ghks
        numpy.array<double, length = len(indices) * len(ghks)> counts

        Return:
            None

        """
        cdef std.vector[string] std_ghks = to_vec_std_strings(ghks)
        self.ptrx().getBatchTriBatchGHKIsNP(&index_array[0], index_array.shape[0], std_ghks, &counts[0], counts.shape[0])

    # ---------------------------------------------------------------------------------
    # ROI section
    # ---------------------------------------------------------------------------------

    def getROITetSpecCounts(self, str roi, str s):
        """
        Get the counts of a species s in tetrehedrons of a ROI.

        Syntax::

            getROITetSpecCounts(roi, s)

        Arguments:
        string roi
        string s

        Return:
        list<float>

        """
        return self.ptrx().getROITetSpecCounts(to_std_string(roi), to_std_string(s))

    def getROITriSpecCounts(self, str roi, str s):
        """
        Get the counts of a species s in triangles of a ROI.

        Syntax::

            getROITriSpecCounts(roi, s)

        Arguments:
        string roi
        string s

        Return:
        list<float>

        """
        return self.ptrx().getROITriSpecCounts(to_std_string(roi), to_std_string(s))

    def getROITetSpecCountsNP(self, str roi, str s, double[:] counts):
        """
        Get the counts of a species s in tetrehedrons of a ROI.

        Syntax::
            getROITetSpecCountsNP(roi, s, counts)

        Arguments:
        string roi
        string s
        numpy.array<float, length = len(indices)>

        Return:
            None

        """
        self.ptrx().getROITetSpecCountsNP(to_std_string(roi), to_std_string(s), &counts[0], counts.shape[0])

    def getROITriSpecCountsNP(self, str roi, str s, double[:] counts):
        """
        Get the counts of a species s in triangles of a ROI.

        Syntax::
            getROITriSpecCountsNP(roi, s, counts)

        Arguments:
        string roi
        string s
        numpy.array<float, length = len(indices)>

        Return:
            None

        """
        self.ptrx().getROITriSpecCountsNP(to_std_string(roi), to_std_string(s), &counts[0], counts.shape[0])

    def getROIVol(self, str roi):
        """
        Get the volume of a ROI.

        Syntax::
            getROIVol(roi)

        Arguments:
        string roi

        Return:
        float

        """
        return self.ptrx().getROIVol(to_std_string(roi))

    def getROIArea(self, str roi):
        """
        Get the area of a ROI.

        Syntax::
            getROIArea(roi)

        Arguments:
        string roi

        Return:
        float

        """
        return self.ptrx().getROIArea(to_std_string(roi))

    def getROISpecCount(self, str roi, str s):
        """
        Get the count of a species in a ROI.

        Syntax::
            getROISpecCount(roi, s)

        Arguments:
        string roi
        string s

        Return:
        float

        """
        return self.ptrx().getROISpecCount(to_std_string(roi), to_std_string(s))

    def setROISpecCount(self, str roi, str s, double count):
        """
        Set the count of a species in a ROI.

        Syntax::
            setROISpecCount(roi, s, count)

        Arguments:
        string roi
        string s
        float count

        Return:
        None

        """
        self.ptrx().setROISpecCount(to_std_string(roi), to_std_string(s), count)

    def getROISpecAmount(self, str roi, str s):
        """
        Get the amount of a species in a ROI.

        Syntax::
            getROISpecAmount(roi, s, count)

        Arguments:
        string roi
        string s

        Return:
        float

        """
        return self.ptrx().getROISpecAmount(to_std_string(roi), to_std_string(s))

    def setROISpecAmount(self, str roi, str s, double amount):
        """
        Set the amount of a species in a ROI.

        Syntax::
            setROISpecAmount(roi, s, amount)

        Arguments:
        string roi
        string s
        float amount

        Return:
        None

        """
        return self.ptrx().setROISpecAmount(to_std_string(roi), to_std_string(s), amount)

    def getROISpecConc(self, str roi, str s):
        """
        Get the concentration of a species in a ROI.

        Syntax::
            getROISpecConc(roi, s)

        Arguments:
        string roi
        string s

        Return:
        float

        """
        return self.ptrx().getROISpecConc(to_std_string(roi), to_std_string(s))

    def setROISpecConc(self, str roi, str s, double conc):
        """
        Set the concentration of a species in a ROI.

        Syntax::
            setROISpecConc(roi, s, conc)

        Arguments:
        string roi
        string s
        float conc

        Return:
        None

        """
        self.ptrx().setROISpecConc(to_std_string(roi), to_std_string(s), conc)

    def setROISpecClamped(self, str roi, str s, bool clamped):
        """
        Set a species in a ROI to be clamped or not. The count of species spec in the ROI is clamped if
        clamped is True, not clamped if clamped is False.

        Syntax::
            setROISpecClamped(roi, s, clamped)

        Arguments:
        string roi
        string s
        bool clamped

        Return:
        None

        """
        self.ptrx().setROISpecClamped(to_std_string(roi), to_std_string(s), clamped)

    def setROIReacK(self, str roi, str r, double kf):
        """
        Sets the macroscopic reaction constant of reaction with identifier string r
        in a ROI with identifier string roi to kf. The unit of the reaction constant
        depends on the order of the reaction.

        Note: The default value still comes from the steps.model description, so
        calling reset() will return the reaction constant to that value.

        Syntax::
            setROIReacK(roi, r, kf)

        Arguments:
        string roi
        string r
        float kf

        Return:
        None

        """
        self.ptrx().setROIReacK(to_std_string(roi), to_std_string(r), kf)

    def setROISReacK(self, str roi, str sr, double kf):
        """
        Sets the macroscopic reaction constant of surface reaction with identifier string sr
        in a ROI with identifier string roi to kf. The unit of the reaction constant
        depends on the order of the reaction.

        Note: The default value still comes from the steps.model description, so
        calling reset() will return the reaction constant to that value.

        Syntax::
            setROISReacK(roi, sr, kf)

        Arguments:
        string roi
        string sr
        float kf

        Return:
        None

        """
        self.ptrx().setROISReacK(to_std_string(roi), to_std_string(sr), kf)

    def setROIDiffD(self, str roi, str diff, double dcst):
        """
        Sets the macroscopic diffusion constant of diffusion with identifier string diff
        in a ROI with identifier string roi to dcst.

        Note: The default value still comes from the steps.model description, so
        calling reset() will return the diffusion constant to that value.

        Syntax::
            setROIDiffD(roi, diff, dcst)

        Arguments:
        string roi
        string diff
        float dcst

        Return:
            None

        """
        self.ptrx().setROIDiffD(to_std_string(roi), to_std_string(diff), dcst)

    def setROIReacActive(self, str roi, str r, bool a):
        """
        Set reaction r in a ROI to be active or not.

        Syntax::
            setROIReacActive(roi, r, a)

        Arguments:
        string roi
        string r
        bool a

        Return:
        None

        """
        self.ptrx().setROIReacActive(to_std_string(roi), to_std_string(r), a)

    def setROISReacActive(self, str roi, str sr, bool a):
        """
        Set surface reaction sr in a ROI to be active or not.

        Syntax::
            setROISReacActive(roi, sr, a)

        Arguments:
        string roi
        string sr
        bool a

        Return:
        None

        """
        self.ptrx().setROISReacActive(to_std_string(roi), to_std_string(sr), a)

    def setROIDiffActive(self, str roi, str d, bool act):
        """
        Set diffusion d in a ROI to be active or not.

        Syntax::
            setROIDiffActive(roi, sr, a)

        Arguments:
        string roi
        string sr
        bool a

        Return:
        None

        """
        self.ptrx().setROIDiffActive(to_std_string(roi), to_std_string(d), act)

    def setROIVDepSReacActive(self, str roi, str vsr, bool a):
        """
        Set voltage dependent surface reaction vsr in a ROI to be active or not.

        Syntax::
            setROIVDepSReacActive(roi, vsr, a)

        Arguments:
        string roi
        string vsr
        bool a

        Return:
        None

        """
        self.ptrx().setROIVDepSReacActive(to_std_string(roi), to_std_string(vsr), a)

    def getROIReacExtent(self, str roi, str r):
        """
        Return the extent of reaction with identifier string reac in ROI with
        identifier string roi, that is the number of times the reaction has occurred up
        to the current simulation time.

        Syntax::
            getROIReacExtent(roi, reac)

        Arguments:
        string roi
        string r

        Return:
        index_t

        """
        return self.ptrx().getROIReacExtent(to_std_string(roi), to_std_string(r))

    def resetROIReacExtent(self, str roi, str r):
        """
        Reset the extent of reaction with identifier string reac in ROI with
        identifier string roi, that is the number of times the reaction has occurred up
        to the current simulation time, to 0.

        Syntax::
            resetROIReacExtent(roi, reac)

        Arguments:
        string roi
        string r

        Return:
        None

        """
        self.ptrx().resetROIReacExtent(to_std_string(roi), to_std_string(r))

    def getROISReacExtent(self, str roi, str sr):
        """
        Return the extent of surface reaction with identifier string sr in ROI with
        identifier string roi, that is the number of times the reaction has occurred up
        to the current simulation time.

        Syntax::
            getROISReacExtent(roi, sr)

        Arguments:
        string roi
        string sr

        Return:
        index_t

        """
        return self.ptrx().getROISReacExtent(to_std_string(roi), to_std_string(sr))

    def resetROISReacExtent(self, str roi, str sr):
        """
        Reset the extent of surface reaction with identifier string reac in ROI with
        identifier string roi, that is the number of times the reaction has occurred up
        to the current simulation time, to 0.

        Syntax::
            resetROISReacExtent(roi, reac)

        Arguments:
        string roi
        string sr

        Return:
        None

        """
        self.ptrx().resetROISReacExtent(to_std_string(roi), to_std_string(sr))

    def getROIDiffExtent(self, str roi, str d):
        """
        Return the extent of diffusion with identifier string diff in ROI with
        identifier string roi, that is the number of times the diffusion has occurred up
        to the current simulation time.

        Syntax::
            getROIDiffExtent(roi, diff)

        Arguments:
        string roi
        string diff

        Return:
        index_t

        """
        return self.ptrx().getROIDiffExtent(to_std_string(roi), to_std_string(d))

    def resetROIDiffExtent(self, str roi, str s):
        """
        Reset the extent of diffusion with identifier string diff in ROI with
        identifier string roi, that is the number of times the diffusion has occurred up
        to the current simulation time, to 0.

        Syntax::
            resetROIDiffExtent(roi, diff)

        Arguments:
        string roi
        string diff

        Return:
        None

        """
        self.ptrx().resetROIDiffExtent(to_std_string(roi), to_std_string(s))


    # ------------------------------------------------------------------------------------------------------------

    def setDiffApplyThreshold(self, int threshold):
        """
        Set the threshold for using binomial distribution for molecule diffusion instead of
        single molecule diffusion.

        If the number of molecules in a tetrahedron await for diffusion is higher than this
        threshold, the solver will use binomial function to distribute these molecules to
        each neighboring tetrahedron. Otherwise the molecules will diffuse one by one.

        The default threshold is 10.

        Syntax::

            setDiffApplyThreshold(threshold)

        Arguments:
        int threshold

        Return:
        None
        """
        self.ptrx().setDiffApplyThreshold(threshold)

    def getReacExtent(self, bool local=False):
        """
        Return the number of reaction events that have happened in the simulation.

        if all processes call this function, it will return the accumulated
        result across all processes. It can also be called in individual process with
        the local argument set to true, in which case it returns the local result of this process.

        By default it is called globally and return the accumulated result.

        Syntax::

            getReacExtent(local)

        Arguments:
        bool local (default = False)

        Return:
        index_t
        """
        return self.ptrx().getReacExtent(local)

    def getDiffExtent(self, bool local=False):
        """
        Return the number of diffusion events that have happened in the simulation.

        if all processes call this function, it will return the accumulated
        result accross all processes. It can also be called in individual process with
        the local argument set to true, in which case it returns the local result of this process.

        By default it is called globally and return the accumlated result.

        Syntax::

            getDiffExtent(local)

        Arguments:
        bool local (default = False)

        Return:
        index_t
        """
        return self.ptrx().getDiffExtent(local)

    def getNIteration(self, ):
        """
        Return the number of Operator-Splitting iterations that have happened in the simulation.

        See (Hepburn et al, 2016) and (Chen et al, 2017) for more detail.

        This function can be called locally.

        Syntax::

            getNIteration()

        Arguments:
        None

        Return:
        float
        """
        return self.ptrx().getNIteration()


    def getUpdPeriod(self, ):
        """
        Return the update period tau of the Operator-Splitting solution.
        See (Hepburn et al, 2016) and (Chen et al, 2017) for more detail.

        Syntax::

            getUpdPeriod()

        Arguments:
        None

        Return:
        float
        """
        return self.ptrx().getUpdPeriod()

    def getCompTime(self, ):
        """
        Return the accumulated computation time of the process.

        To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
        line in the src/CmakeLists.txt file. This function is always called and return result locally.

        See (Chen, 2017) for more detail.

        Syntax::

            getCompTime()

        Arguments:
        None

        Return:
        float
        """
        return self.ptrx().getCompTime()

    def getSyncTime(self, ):
        """
        Return the accumulated synchronization time of the process.

        To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
        line in the src/CmakeLists.txt file. This function is always called and return result locally.

        See (Chen, 2017) for more detail.

        Syntax::

            getSyncTime()

        Arguments:
        None

        Return:
        float
        """
        return self.ptrx().getSyncTime()

    def getIdleTime(self, ):
        """
        Return the accumulated idle time of the process.

        To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
        line in the src/CmakeLists.txt file. This function is always called and return result locally.

        See (Chen, 2017) for more detail.

        Syntax::

            getSyncTime()

        Arguments:
        None

        Return:
        float
        """
        return self.ptrx().getIdleTime()

    def getEFieldTime(self, ):
        """
        Return the accumulated EField run time of the process.

        To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
        line in the src/CmakeLists.txt file. This function is always called and return result locally.


        Syntax::

            getEFieldTime()

        Arguments:
        None

        Return:
        float
        """
        return self.ptrx().getEFieldTime()

    def getRDTime(self, ):
        """
        Return the accumulated reaction-diffusion run time of the process.

        To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
        line in the src/CmakeLists.txt file. This function is always called and return result locally.

        Syntax::

            getRDTime()

        Arguments:
        None

        Return:
        float
        """
        return self.ptrx().getRDTime()

    def getDataExchangeTime(self, ):
        """
        Return the accumulated data exchanging time between RD and EField solvers of the process.

        To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
        line in the src/CmakeLists.txt file. This function is always called and return result locally.

        Syntax::

            getDataExchangeTime()

        Arguments:
        None

        Return:
        float
        """
        return self.ptrx().getDataExchangeTime()

    def repartitionAndReset(self, std.vector[int] tet_hosts=(), dict tri_hosts=None, std.vector[int] wm_hosts=()):
        """
        Repartition and reset the simulation.

        Note: It is possible to repartitioning the mesh using any subset of the processes,
        in which case processes with no assigned subvolumes will mostly idle for the working processes.

        Syntax::

            repartitionAndReset(tet_hosts, tri_hosts, wm_hosts)

        Arguments:
        list tet_hosts
        dict tri_hosts (default = {})
        dict wm_hosts (default = {})

        Return:
            None
        """
        if tri_hosts is None: tri_hosts = {}
        cdef std.map[unsigned int, int] _tri_hosts = tri_hosts
        self.ptrx().repartitionAndReset(tet_hosts, _tri_hosts, wm_hosts)

    @staticmethod
    cdef _py_TetOpSplitP from_ptr(TetOpSplitP *ptr):
        cdef _py_TetOpSplitP obj = _py_TetOpSplitP.__new__(_py_TetOpSplitP )
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_TetOpSplitP from_ref(const TetOpSplitP &ref):
        _py_TetOpSplitP.from_ptr(<TetOpSplitP*>&ref)

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TetVesicleRDEF(_py_TetAPI):
    """Bindings for MPI TetVesicleRDEF"""
# ----------------------------------------------------------------------------------------------------------------------
    cdef TetVesicleRDEF *ptrx(self):
        return <TetVesicleRDEF*> self._ptr

    def __init__(self, _py_Model model, _py_Geom geom, _py_RNG rng, int calcMembPot=0):
        """
        Construction::

            sim = steps.mpi.solver.TetVesicleRDEF(model, geom, rng, calcMembPot=0)

        Create a spatial stochastic solver based on operator-splitting, which also supports vesicles, 'rafts' and related phenomena such as exocytosis and endocytosis.
        If voltage is to be simulated, argument calcMembPot specifies the solver. E.g. calcMembPot=steps.solver.EF_DV_PETSC will utilise the PETSc library. calcMembPot=0 means that voltage will not be simulated.

        Arguments:
        steps.model.Model model
        steps.geom.Geom geom
        steps.rng.RNG rng
        int calcMemPot (default=0)

        """
        if model == None:
            raise TypeError('The Model object is empty.')
        if geom == None:
            raise TypeError('The Geom object is empty.')
        if rng == None:
            raise TypeError('The RNG object is empty.')
        self._ptr = new TetVesicleRDEF(model.ptr(), geom.ptr(), rng.ptr(), calcMembPot)

    def getSolverName(self, ):
        """
        Return the solver's name.


        :rtype: str
        """
        return from_std_string(self.ptrx().getSolverName())

    def getSolverDesc(self, ):
        """
        Return the solver's description.


        :rtype: str
        """
        return from_std_string(self.ptrx().getSolverDesc())

    def getSolverAuthors(self, ):
        """
        Return the solver's author.


        :rtype: str
        """
        return from_std_string(self.ptrx().getSolverAuthors())

    def getSolverEmail(self, ):
        """
        Return the solver's email.


        :rtype: str
        """
        return from_std_string(self.ptrx().getSolverEmail())

    def reset(self, ):
        """
        Reset the solver.


        :rtype: None
        """
        self.ptrx().reset()

    def run(self, double endtime):
        """
        Advance the simulation until endtime (given in seconds) is reached.
        The endtime must be larger or equal to the current simulation time.

        :param endtime: Time to end the solver.
        :type endtime: float


        :rtype: None
        """
        self.ptrx().run(endtime)

    def advance(self, double adv):
        """
        Advance the solver a given amount of time.

        :param adv: Time to advance the solver (in seconds)
        :type adv: float


        :rtype: None
        """
        self.ptrx().advance(adv)

    def step(self, ):
        """
        Advance the simulation for one 'step'. In stochastic solvers this is
        one 'realization' of the Gillespie SSA (one reaction 'event'). In
        numerical solvers (currently Wmrk4) this is one time-step, with the
        stepsize defined with the setDT method.


        :rtype: None
        """
        self.ptrx().step()

    def checkpoint(self, str file_name):
        """
        checkpoint simulator state to a file

        :param file_name: file name prefix (will be suffixed with the rank of each process)
        :type file_name: str


        :rtype: None
        """
        fullName = file_name + str(f'_{getRank()}')
        self.ptrx().checkpoint(to_std_string(fullName))

    def restore(self, str file_name):
        """
        restore simulator state from a file

        :param file_name: file name prefix (does not include the rank suffix)
        :type file_name: str


        :rtype: None
        """
        fullName = file_name + str(f'_{getRank()}')
        self.ptrx().restore(to_std_string(fullName))

    def setEfieldDT(self, double efdt):
        """
        Set the stepsize for membrane potential solver (default 1e-6s).
        This is the time for each voltage calculation step. The SSA will
        run until passing this stepsize, so in fact each membrane potential
        time step will vary slightly around the dt so as to be aligned with the
        SSA.

        :param efdt:
        :type efdt: float


        :rtype: None
        """
        self.ptrx().setEfieldDT(efdt)

    def setTemp(self, double temp):
        """
        Set the simulation temperature. Currently, this will only
        influence the GHK flux rate, so will only influence simulations
        including membrane potential calculation.

        :param temp: EField temperature (in Kelvin)
        :type temp: float


        :rtype: None
        """
        self.ptrx().setTemp(temp)

    def getTemp(self, ):
        """
        Return the simulation temperature (in Kelvin)


        :rtype: float
        """
        return self.ptrx().getTemp()

    def getTime(self, ):
        """
        Returns the current simulation time in seconds.


        :rtype: float
        """
        return self.ptrx().getTime()

    def getA0(self, ):
        """
        Returns the total propensity of the current simulation state
        (the total propensity multiplied by an infinitesimally small
        time dt gives the probability that a reaction will occur in that dt).
        For Tetexact this includes the propensity from the extension of the SSA
        for diffusive flux between tetrahedral elements in the mesh.


        :rtype: float
        """
        return self.ptrx().getA0()

    def getNSteps(self, ):
        """
        Return the number of 'realizations' of the SSA, the number of reaction
        (and diffusion) events in stochastic solvers.


        :rtype: int
        """
        return self.ptrx().getNSteps()

    def setTime(self, double time):
        """
        Set the current simulation time.

        Syntax::

            setTime(time)

        Arguments:
        float time

        Return:
        None

        """
        self.ptrx().setTime(time)

    def setVesicleDT(self, double dt):
        """
        Set the default vesicle dt (Note: the actual vesicle dt
        used in the simulation can be lower than this number depending
        on simulation conditions).

        :param dt: Vesicle dt (in seconds)
        :type dt: float


        :rtype: None
        """
        self.ptrx().setVesicleDT(dt)

    def getVesicleDT(self, ):
        """
        Return the present vesicle dt (in seconds)


        :rtype: float
        """
        return self.ptrx().getVesicleDT()

    def createPath(self, str path):
        """
        Create a path

        :param path: Name of the path
        :type path: str


        :rtype: None
        """
        self.ptrx().createPath(to_std_string(path))

    def addPathPoint(self, str path, uint point_idx, std.vector[double] position):
        """
        Add a point to a path

        :param path: Name of the path
        :type path: str
        :param point_idx: An index for the point, positive integer
        :type point_idx: int
        :param position: Position of the point in cartesian coordinates
        :type position: List[float]


        :rtype: None
        """
        self.ptrx().addPathPoint(to_std_string(path), point_idx, position)

    def addPathBranch(self, str path, uint sourcepoint_idx, std.map[uint, double] destpoints_indxs):
        """
        Create branching in the path from a point

        :param path: Name of the path
        :type path: str
        :param sourcepoint_idx: An index for the source point, positive integer
        :type sourcepoint_idx: int
        :param destpoints_indxs:
        :type destpoints_indxs: Dict[int, float]


        :rtype: None
        """
        self.ptrx().addPathBranch(to_std_string(path), sourcepoint_idx, destpoints_indxs)

    def getBatchTetSpecCounts(self, std.vector[steps.index_t] tets, str s):
        """
        Get species counts of a list of tetrahedrons

        :param tets:
        :type tets: List[steps.index_t]
        :param s:
        :type s: str


        :rtype: List[float]
        """
        return [_val1 for _val1 in self.ptrx().getBatchTetSpecCounts(tets, to_std_string(s))]

    def getBatchTriSpecCounts(self, std.vector[steps.index_t] tris, str s):
        """
        Get species counts of a list of triangles

        :param tris:
        :type tris: List[steps.index_t]
        :param s:
        :type s: str


        :rtype: List[float]
        """
        return [_val1 for _val1 in self.ptrx().getBatchTriSpecCounts(tris, to_std_string(s))]

    def getBatchTetSpecCountsNP(self, index_t[:] indices, str s, double[:] counts):
        """
        Get species counts of a list of tetrahedrons

        :param indices:
        :type indices: numpy.array<index_t>
        :param s:
        :type s: str
        :param counts:
        :type counts: numpy.array<float>


        :rtype: None
        """
        self.ptrx().getBatchTetSpecCountsNP(&indices[0], indices.shape[0], to_std_string(s), &counts[0], counts.shape[0])

    def getBatchTriSpecCountsNP(self, index_t[:] indices, str s, double[:] counts):
        """
        Get species counts of a list of triangles

        :param indices:
        :type indices: numpy.array<index_t>
        :param s:
        :type s: str
        :param counts:
        :type counts: numpy.array<float>


        :rtype: None
        """
        self.ptrx().getBatchTriSpecCountsNP(&indices[0], indices.shape[0], to_std_string(s), &counts[0], counts.shape[0])

    def setROITetSpecClamped(self, std.vector[steps.index_t] triangles, str s, bool b):
        """

        :param triangles:
        :type triangles: List[steps.index_t]
        :param s:
        :type s: str
        :param b:
        :type b: bool


        :rtype: None
        """
        cdef std.vector[steps.tetrahedron_global_id] std_triangles
        std_triangles.reserve(len(triangles))
        for _val in triangles:
            std_triangles.push_back(steps.tetrahedron_global_id(_val))
        self.ptrx().setROITetSpecClamped(std_triangles, to_std_string(s), b)

    def setROITriSpecClamped(self, std.vector[steps.index_t] triangles, str s, bool b):
        """

        :param triangles:
        :type triangles: List[steps.index_t]
        :param s:
        :type s: str
        :param b:
        :type b: bool


        :rtype: None
        """
        cdef std.vector[steps.triangle_global_id] std_triangles
        std_triangles.reserve(len(triangles))
        for _val in triangles:
            std_triangles.push_back(steps.triangle_global_id(_val))
        self.ptrx().setROITriSpecClamped(std_triangles, to_std_string(s), b)

    def getROITetSpecCount(self, std.vector[steps.index_t] triangles, str s):
        """

        :param triangles:
        :type triangles: List[steps.index_t]
        :param s:
        :type s: str


        :rtype: float
        """
        cdef std.vector[steps.tetrahedron_global_id] std_triangles
        std_triangles.reserve(len(triangles))
        for _val in triangles:
            std_triangles.push_back(steps.tetrahedron_global_id(_val))
        return self.ptrx().getROITetSpecCount(std_triangles, to_std_string(s))

    def getROITriSpecCount(self, std.vector[steps.index_t] triangles, str s):
        """

        :param triangles:
        :type triangles: List[steps.index_t]
        :param s:
        :type s: str


        :rtype: float
        """
        cdef std.vector[steps.triangle_global_id] std_triangles
        std_triangles.reserve(len(triangles))
        for _val in triangles:
            std_triangles.push_back(steps.triangle_global_id(_val))
        return self.ptrx().getROITriSpecCount(std_triangles, to_std_string(s))

    def setROITetSpecCount(self, std.vector[steps.index_t] triangles, str s, double count):
        """

        :param triangles:
        :type triangles: List[steps.index_t]
        :param s:
        :type s: str
        :param count:
        :type count: float


        :rtype: None
        """
        cdef std.vector[steps.tetrahedron_global_id] std_triangles
        std_triangles.reserve(len(triangles))
        for _val in triangles:
            std_triangles.push_back(steps.tetrahedron_global_id(_val))
        self.ptrx().setROITetSpecCount(std_triangles, to_std_string(s), count)

    def setROITriSpecCount(self, std.vector[steps.index_t] triangles, str s, double count):
        """

        :param triangles:
        :type triangles: List[steps.index_t]
        :param s:
        :type s: str
        :param count:
        :type count: float


        :rtype: None
        """
        cdef std.vector[steps.triangle_global_id] std_triangles
        std_triangles.reserve(len(triangles))
        for _val in triangles:
            std_triangles.push_back(steps.triangle_global_id(_val))
        self.ptrx().setROITriSpecCount(std_triangles, to_std_string(s), count)

    def getROITetSpecCounts(self, str roi, str s):
        """
        Get species counts of a list of tetrahedrons

        :param roi:
        :type roi: str
        :param s:
        :type s: str


        :rtype: List[float]
        """
        return [_val1 for _val1 in self.ptrx().getROITetSpecCounts(to_std_string(roi), to_std_string(s))]

    def getROITriSpecCounts(self, str roi, str s):
        """
        Get species counts of a list of triangles

        :param roi:
        :type roi: str
        :param s:
        :type s: str


        :rtype: List[float]
        """
        return [_val1 for _val1 in self.ptrx().getROITriSpecCounts(to_std_string(roi), to_std_string(s))]

    def getROITetSpecCountsNP(self, str roi, str s, double[:] counts):
        """
        Get species counts of a list of tetrahedrons

        :param roi:
        :type roi: str
        :param s:
        :type s: str
        :param counts:
        :type counts: numpy.array<float>


        :rtype: None
        """
        self.ptrx().getROITetSpecCountsNP(to_std_string(roi), to_std_string(s), &counts[0], counts.shape[0])

    def getROITriSpecCountsNP(self, str roi, str s, double[:] counts):
        """
        Get species counts of a list of triangles

        :param roi:
        :type roi: str
        :param s:
        :type s: str
        :param counts:
        :type counts: numpy.array<float>


        :rtype: None
        """
        self.ptrx().getROITriSpecCountsNP(to_std_string(roi), to_std_string(s), &counts[0], counts.shape[0])

    def getROIVol(self, str roi):
        """
        Get the volume of a ROI.

        :param roi:
        :type roi: str


        :rtype: float
        """
        return self.ptrx().getROIVol(to_std_string(roi))

    def getROIArea(self, str roi):
        """
        Get the area of a ROI.

        :param roi:
        :type roi: str


        :rtype: float
        """
        return self.ptrx().getROIArea(to_std_string(roi))

    def getROISpecCount(self, str roi, str s):
        """
        Get the count of a species in a ROI.

        :param roi:
        :type roi: str
        :param s:
        :type s: str


        :rtype: float
        """
        return self.ptrx().getROISpecCount(to_std_string(roi), to_std_string(s))

    def setROISpecCount(self, str roi, str s, double count):
        """
        Set the count of a species in a ROI.

        :param roi:
        :type roi: str
        :param s:
        :type s: str
        :param count:
        :type count: float


        :rtype: None
        """
        self.ptrx().setROISpecCount(to_std_string(roi), to_std_string(s), count)

    def getROISpecAmount(self, str roi, str s):
        """
        Get the amount of a species in a ROI.

        :param roi:
        :type roi: str
        :param s:
        :type s: str


        :rtype: float
        """
        return self.ptrx().getROISpecAmount(to_std_string(roi), to_std_string(s))

    def setROISpecAmount(self, str roi, str s, double amount):
        """
        Set the amount of a species in a ROI.

        :param roi:
        :type roi: str
        :param s:
        :type s: str
        :param amount:
        :type amount: float


        :rtype: None
        """
        self.ptrx().setROISpecAmount(to_std_string(roi), to_std_string(s), amount)

    def getROISpecConc(self, str roi, str s):
        """
        Get the concentration of a species in a ROI.

        :param roi:
        :type roi: str
        :param s:
        :type s: str


        :rtype: float
        """
        return self.ptrx().getROISpecConc(to_std_string(roi), to_std_string(s))

    def setROISpecConc(self, str roi, str s, double conc):
        """
        Set the concentration of a species in a ROI.

        :param roi:
        :type roi: str
        :param s:
        :type s: str
        :param conc:
        :type conc: float


        :rtype: None
        """
        self.ptrx().setROISpecConc(to_std_string(roi), to_std_string(s), conc)

    def setROISpecClamped(self, str roi, str s, bool clamped):
        """
        Set a species in a ROI to be clamped or not. The count of species spec in the ROI is clamped if
        clamped is True, not clamped if clamped is False.

        Syntax::
            setROISpecClamped(roi, s, clamped)

        Arguments:
        string roi
        string s
        bool clamped

        Return:
        None

        """
        self.ptrx().setROISpecClamped(to_std_string(roi), to_std_string(s), clamped)

    def setROIReacK(self, str roi, str r, double kf):
        """
        Sets the macroscopic reaction constant of reaction with identifier
        string r in a ROI with identifier string roi to kf. The unit of the
        reaction constant depends on the order of the reaction.
        Note: The default value still comes from the steps.model description,
        so calling reset() will return the reaction constant to that value.

        :param roi:
        :type roi: str
        :param r:
        :type r: str
        :param kf:
        :type kf: float


        :rtype: None
        """
        self.ptrx().setROIReacK(to_std_string(roi), to_std_string(r), kf)

    def setROISReacK(self, str roi, str sr, double kf):
        """
        Sets the macroscopic reaction constant of surface reaction with
        identifier string sr in a ROI with identifier string roi to kf. The
        unit of the reaction constant depends on the order of the reaction.
        Note: The default value still comes from the steps.model description,
        so calling reset() will return the reaction constant to that value.

        :param roi:
        :type roi: str
        :param sr:
        :type sr: str
        :param kf:
        :type kf: float


        :rtype: None
        """
        self.ptrx().setROISReacK(to_std_string(roi), to_std_string(sr), kf)

    def setROIDiffD(self, str roi, str diff, double dcst):
        """
        Sets the macroscopic diffusion constant of diffusion with identifier string diff
        in a ROI with identifier string roi to dcst.

        Note: The default value still comes from the steps.model description, so
        calling reset() will return the diffusion constant to that value.

        Syntax::
            setROIDiffD(roi, diff, dcst)

        Arguments:
        string roi
        string diff
        float dcst

        Return:
            None

        """
        self.ptrx().setROIDiffD(to_std_string(roi), to_std_string(diff), dcst)

    def setROIReacActive(self, str roi, str r, bool a):
        """
        Set reaction r in a ROI to be active or not.

        :param roi:
        :type roi: str
        :param r:
        :type r: str
        :param a:
        :type a: bool


        :rtype: None
        """
        self.ptrx().setROIReacActive(to_std_string(roi), to_std_string(r), a)

    def setROISReacActive(self, str roi, str sr, bool a):
        """
        Set surface reaction sr in a ROI to be active or not.

        :param roi:
        :type roi: str
        :param sr:
        :type sr: str
        :param a:
        :type a: bool


        :rtype: None
        """
        self.ptrx().setROISReacActive(to_std_string(roi), to_std_string(sr), a)

    def setROIDiffActive(self, str roi, str d, bool act):
        """
        Set diffusion d in a ROI to be active or not.

        :param roi:
        :type roi: str
        :param d:
        :type d: str
        :param act:
        :type act: bool


        :rtype: None
        """
        self.ptrx().setROIDiffActive(to_std_string(roi), to_std_string(d), act)

    def setROIVDepSReacActive(self, str roi, str vsr, bool a):
        """
        Set voltage dependent surface reaction vsr in a ROI to be active or
        not.

        :param roi:
        :type roi: str
        :param vsr:
        :type vsr: str
        :param a:
        :type a: bool


        :rtype: None
        """
        self.ptrx().setROIVDepSReacActive(to_std_string(roi), to_std_string(vsr), a)

    def getROIReacExtent(self, str roi, str r):
        """
        Return the extent of reaction with identifier string reac in ROI with
        identifier string roi, that is the number of times the reaction has
        occurred up to the current simulation time.

        :param roi:
        :type roi: str
        :param r:
        :type r: str


        :rtype: unsigned long long
        """
        return self.ptrx().getROIReacExtent(to_std_string(roi), to_std_string(r))

    def resetROIReacExtent(self, str roi, str r):
        """
        Reset the extent of reaction with identifier string reac in ROI with
        identifier string roi, that is the number of times the reaction has
        occurred up to the current simulation time, to 0.

        :param roi:
        :type roi: str
        :param r:
        :type r: str


        :rtype: None
        """
        self.ptrx().resetROIReacExtent(to_std_string(roi), to_std_string(r))

    def getROISReacExtent(self, str roi, str sr):
        """
        Return the extent of surface reaction with identifier string sr in ROI
        with identifier string roi, that is the number of times the reaction
        has occurred up to the current simulation time.

        :param roi:
        :type roi: str
        :param sr:
        :type sr: str


        :rtype: unsigned long long
        """
        return self.ptrx().getROISReacExtent(to_std_string(roi), to_std_string(sr))

    def resetROISReacExtent(self, str roi, str sr):
        """
        Reset the extent of surface reaction with identifier string reac in ROI
        with identifier string roi, that is the number of times the reaction
        has occurred up to the current simulation time, to 0.

        :param roi:
        :type roi: str
        :param sr:
        :type sr: str


        :rtype: None
        """
        self.ptrx().resetROISReacExtent(to_std_string(roi), to_std_string(sr))

    def getROIDiffExtent(self, str roi, str d):
        """
        Return the extent of diffusion with identifier string diff in ROI with
        identifier string roi, that is the number of times the diffusion has
        occurred up to the current simulation time.

        :param roi:
        :type roi: str
        :param d:
        :type diff: str


        :rtype: unsigned long long
        """
        return self.ptrx().getROIDiffExtent(to_std_string(roi), to_std_string(d))

    def resetROIDiffExtent(self, str roi, str s):
        """
        Reset the extent of diffusion with identifier string diff in ROI with
        identifier string roi, that is the number of times the diffusion has
        occurred up to the current simulation time, to 0.

        :param roi:
        :type roi: str
        :param s:
        :type s: str


        :rtype: None
        """
        self.ptrx().resetROIDiffExtent(to_std_string(roi), to_std_string(s))

    def setDiffApplyThreshold(self, int threshold):
        """
        Set the threshold for using binomial distribution for molecule diffusion instead of
        single molecule diffusion.
        If the number of molecules in a tetrahedron await for diffusion is higher than this
        threshold, the solver will use binomial function to distribute these molecules to
        each neighboring tetrahedron. Otherwise the molecules will diffuse one by one.
        The default threshold is 10.

        :param threshold:
        :type threshold: int


        :rtype: None
        """
        self.ptrx().setDiffApplyThreshold(threshold)

    def getReacExtent(self, bool local):
        """
        Return the number of reaction events that have happened in the simulation.
        if all processes call this function, it will return the accumulated
        result across all processes. It can also be called in individual process with
        the local argument set to true, in which case it returns the local result of this process.
        By default it is called globally and return the accumulated result.

        :param local:
        :type local: bool


        :rtype: unsigned long long
        """
        return self.ptrx().getReacExtent(local)

    def getDiffExtent(self, bool local):
        """
        Return the number of diffusion events that have happened in the simulation.
        if all processes call this function, it will return the accumulated
        result accross all processes. It can also be called in individual process with
        the local argument set to true, in which case it returns the local result of this process.
        By default it is called globally and return the accumlated result.

        :param local:
        :type local: bool


        :rtype: unsigned long long
        """
        return self.ptrx().getDiffExtent(local)

    def getNIteration(self, ):
        """
        Return the number of Operator-Splitting iterations that have happened in the simulation.
        See (Hepburn et al, 2016) and (Chen et al, 2017) for more detail.
        This function can be called locally.


        :rtype: float
        """
        return self.ptrx().getNIteration()

    def getUpdPeriod(self, ):
        """
        Return the update period tau of the Operator-Splitting solution.
        See (Hepburn et al, 2016) and (Chen et al, 2017) for more detail.


        :rtype: float
        """
        return self.ptrx().getUpdPeriod()

    def getCompTime(self, ):
        """
        Return the accumulated computation time of the process.
        To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
        line in the src/CmakeLists.txt file. This function is always called and return result
        locally.
        See (Chen, 2017) for more detail.


        :rtype: float
        """
        return self.ptrx().getCompTime()

    def getSyncTime(self, ):
        """
        Return the accumulated synchronization time of the process.
        To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
        line in the src/CmakeLists.txt file. This function is always called and return result
        locally.
        See (Chen, 2017) for more detail.


        :rtype: float
        """
        return self.ptrx().getSyncTime()

    def getIdleTime(self, ):
        """
        Return the accumulated idle time of the process.
        To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
        line in the src/CmakeLists.txt file. This function is always called and return result
        locally.
        See (Chen, 2017) for more detail.


        :rtype: float
        """
        return self.ptrx().getIdleTime()

    def getEFieldTime(self, ):
        """
        Return the accumulated EField run time of the process.
        To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
        line in the src/CmakeLists.txt file. This function is always called and return result
        locally.


        :rtype: float
        """
        return self.ptrx().getEFieldTime()

    def getRDTime(self, ):
        """
        Return the accumulated reaction-diffusion run time of the process.
        To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
        line in the src/CmakeLists.txt file. This function is always called and return result
        locally.


        :rtype: float
        """
        return self.ptrx().getRDTime()

    def getDataExchangeTime(self, ):
        """
        Return the accumulated data exchanging time between RD and EField solvers of the process.
        To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
        line in the src/CmakeLists.txt file. This function is always called and return result
        locally.


        :rtype: float
        """
        return self.ptrx().getDataExchangeTime()

    def getEfieldDT(self, ):
        """
        Return the stepsize for membrane potential solver (in seconds)


        :rtype: float
        """
        return self.ptrx().getEfieldDT()

    def getCompVol(self, str c):
        """
        Returns the volume of compartment c (in m^3).

        :param c: Name of the compartment.
        :type c: str


        :rtype: float
        """
        return self.ptrx().getCompVol(to_std_string(c))

    def getCompSpecCount(self, str c, str s):
        """
        Returns the number of molecules of species s in compartment c.
        NOTE: in a mesh-based simulation, the total count is computed as
        the sum of the counts in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getCompSpecCount(to_std_string(c), to_std_string(s))

    def setCompSpecCount(self, str c, str s, double n):
        """
        Sets the number of molecules of species s in compartment c.
        NOTE: in a mesh-based simulation, the total amount is equally divided
        over all tetrahedrons in the compartment (i.e. a uniform distribution).

        :param c: Name of the compartment.
        :type c: str
        :param s: Name of the species.
        :type s: str
        :param n: Number of molecules of the species.
        :type n: float


        :rtype: None
        """
        self.ptrx().setCompSpecCount(to_std_string(c), to_std_string(s), n)

    def getCompSpecAmount(self, str c, str s):
        """
        Returns the amount (in mols) of species s in compartment c.
        NOTE: in a mesh-based simulation, the total amount is computed as
        the sum of the amounts in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getCompSpecAmount(to_std_string(c), to_std_string(s))

    def setCompSpecAmount(self, str c, str s, double a):
        """
        Set the amount (in mols) of species s in compartment c.
        NOTE: in a mesh-based simulation, the total amount is equally divided
        over all tetrahedrons in the compartment (i.e. a uniform distribution).

        :param c: Name of the compartment.
        :type c: str
        :param s: Name of the species.
        :type s: str
        :param a: Amount of the species.
        :type a: float


        :rtype: None
        """
        self.ptrx().setCompSpecAmount(to_std_string(c), to_std_string(s), a)

    def getCompSpecConc(self, str c, str s):
        """
        Returns the concentration (in molar units) of species s in compartment
        c.
        NOTE: in a mesh-based simulation, the overall concentration in a
        compartment is computed by taking the volume-weighted sum of the
        concentration in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getCompSpecConc(to_std_string(c), to_std_string(s))

    def setCompSpecConc(self, str c, str s, double conc):
        """
        Sets the concentration (in molar units) of species s in compartment c.
        NOTE: in a mesh-based simulation, this method changes the
        concentration to the same value in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param s: Name of the species.
        :type s: str
        :param conc: Concentration of the species.
        :type conc: float


        :rtype: None
        """
        self.ptrx().setCompSpecConc(to_std_string(c), to_std_string(s), conc)

    def getCompSpecClamped(self, str c, str s):
        """
        Returns whether the concentration of species s in compartment c
        remains constant over time (unless changed explicitly).
        NOTE: in a mesh-based simulation, this method will only return true
        only if the species has been clamped in all tetrahedrons of the
        compartment. \param c Name of the compartment. \param s Name of the
        species.

        :param c:
        :type c: str
        :param s:
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getCompSpecClamped(to_std_string(c), to_std_string(s))

    def setCompSpecClamped(self, str c, str s, bool b):
        """
        Turns clamping of species s in compartment c on or off.
        NOTE: in a mesh based simulation, this method turns clamping on/off
        in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param s: Name of the species.
        :type s: str
        :param b: Flag to trun clamping of species on / off.
        :type b: bool


        :rtype: None
        """
        self.ptrx().setCompSpecClamped(to_std_string(c), to_std_string(s), b)

    def getCompReacK(self, str c, str r):
        """
        Returns the macroscopic reaction constant of reaction r in
        compartment c.
        Note: in a mesh-based simulation, the value is computed as the
        volume-weighted sum of the reaction constants in all tetrahedrons of
        the compartment.

        :param c:
        :type c: str
        :param r:
        :type r: str


        :rtype: float
        """
        return self.ptrx().getCompReacK(to_std_string(c), to_std_string(r))

    def setCompReacK(self, str c, str r, double kf):
        """
        Sets the macroscopic reaction constant of reaction r in compartment c
        (units vary according to the order of the reaction).
        NOTE: in a mesh-based simulation, this method changes the reaction
        constant equally in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of te reaction.
        :type r: str
        :param kf: Reaction constant.
        :type kf: float


        :rtype: None
        """
        self.ptrx().setCompReacK(to_std_string(c), to_std_string(r), kf)

    def getCompReacActive(self, str c, str r):
        """
        Returns whether reaction r in compartment c is active or not
        NOTE: in a mesh-based simulation, this method returns false only when
        the reaction has been inactivated in all tetrahedrons.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: bool
        """
        return self.ptrx().getCompReacActive(to_std_string(c), to_std_string(r))

    def setCompReacActive(self, str c, str r, bool a):
        """
        Activate or inactivate reaction r in compartment c.
        NOTE: in a mesh-based simulation, activation/inactivation of a reaction
        turns it on/off in all tetrahedrons.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of the reaction.
        :type r: str
        :param a: Flag to activate or deactivate the reaction.
        :type a: bool


        :rtype: None
        """
        self.ptrx().setCompReacActive(to_std_string(c), to_std_string(r), a)

    def getCompDiffD(self, str c, str d):
        """
        Returns the diffusion constant of diffusion rule d in compartment c.

        :param c: Name of the compartment.
        :type c: str
        :param d: Name of the diffusion.
        :type d: str


        :rtype: float
        """
        return self.ptrx().getCompDiffD(to_std_string(c), to_std_string(d))

    def setCompDiffD(self, str c, str d, double dcst):
        """
        Set the diffusion constant of diffusion rule d in compartment c.

        :param c: Name of the compartment.
        :type c: str
        :param d: Name of the diffusion.
        :type d: str
        :param dcst: Rate constant of the diffusion.
        :type dcst: float


        :rtype: None
        """
        self.ptrx().setCompDiffD(to_std_string(c), to_std_string(d), dcst)

    def getCompDiffActive(self, str c, str d):
        """
        Returns whether diffusion rule d in compartment c is active or not.

        :param c: Name of the compartment.
        :type c: str
        :param d: Name of the diffusion.
        :type d: str


        :rtype: bool
        """
        return self.ptrx().getCompDiffActive(to_std_string(c), to_std_string(d))

    def setCompDiffActive(self, str c, str d, bool act):
        """
        Activate or deactivate diffusion rule d in compartment c.

        :param c: Name of the compartment.
        :type c: str
        :param d: Name of the diffusion.
        :type d: str
        :param act: Flag to activate or deactivate the diffusion.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setCompDiffActive(to_std_string(c), to_std_string(d), act)

    def getCompReacC(self, str c, str r):
        """
        Returns c_mu, the mesoscopic reaction constant of reaction r in
        compartment c.
        NOTE: in a mesh-based simulation, the mesoscopic reaction constant is
        computed as the sum of the mesoscopic constants in all tetrahedrons of
        the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getCompReacC(to_std_string(c), to_std_string(r))

    def getCompReacH(self, str c, str r):
        """
        Returns h_mu, the distinct number of ways in which reaction r can
        occur in compartment c, by computing the product of its reactants.
        NOTE: in a mesh-based simulation, it returns the sum of the h_mu's
        over all tetrahedrons of the compartment. This can become a very large
        value.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getCompReacH(to_std_string(c), to_std_string(r))

    def getCompReacA(self, str c, str r):
        """
        Returns the propensity, a_mu, of reaction r in compartment c.
        The propensity value gives the probability per unit time that this
        reaction will occur in the current state.
        NOTE: in a mesh-based simulation, a_mu is computed as the sum of the
        a_mu in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getCompReacA(to_std_string(c), to_std_string(r))

    def getCompReacExtent(self, str c, str r):
        """
        Returns the extent of reaction r in compartment c.
        NOTE: in a mesh-based simulation, returns the sum of the reaction
        extents in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: unsigned long long
        """
        return self.ptrx().getCompReacExtent(to_std_string(c), to_std_string(r))

    def resetCompReacExtent(self, str c, str r):
        """
        Resets the extent of reaction r in compartment c to zero.
        NOTE: in a mesh-based simulation, resets the extents of the reaction
        in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: None
        """
        self.ptrx().resetCompReacExtent(to_std_string(c), to_std_string(r))

    def getCompVesicleCount(self, str c, str v):
        """
        Returns the number of vesicles v in compartment c

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str


        :rtype: int
        """
        return self.ptrx().getCompVesicleCount(to_std_string(c), to_std_string(v))

    def setCompVesicleCount(self, str c, str v, uint n):
        """
        Sets the number of vesicles v in compartment c

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str
        :param n: Number of vesicles
        :type n: int


        :rtype: None
        """
        self.ptrx().setCompVesicleCount(to_std_string(c), to_std_string(v), n)

    def addCompVesicle(self, str c, str v):
        """
        Adds one vesicle v to compartment c

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str


        :rtype: steps.index_t
        """
        return self.ptrx().addCompVesicle(to_std_string(c), to_std_string(v)).get()

    def deleteSingleVesicle(self, str v, steps.index_t ves_unique_index):
        """
        Deletes individual vesicle of type v with unique index ves_unique_index

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t


        :rtype: None
        """
        self.ptrx().deleteSingleVesicle(to_std_string(v), steps.vesicle_individual_id(ves_unique_index))

    def getSingleVesicleSurfaceLinkSpecCount(self, str v, steps.index_t ves_unique_index, str ls):
        """
        Returns the count of 'link' species ls on vesicle of type v with unique index
        vesicle_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param ls: Name of the link species
        :type ls: str


        :rtype: int
        """
        return self.ptrx().getSingleVesicleSurfaceLinkSpecCount(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(ls))

    def getSingleVesicleSurfaceLinkSpecIndices(self, str v, steps.index_t ves_unique_index, str ls):
        """
        Returns the indices of 'link' species ls on vesicle of type v with unique index
        vesicle_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param ls: Name of the link species
        :type ls: str


        :rtype: List[steps.index_t]
        """
        return [_val1.get() for _val1 in self.ptrx().getSingleVesicleSurfaceLinkSpecIndices(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(ls))]

    def getSingleVesicleSurfaceSpecIndices(self, str v, steps.index_t ves_unique_index, str s):
        """
        Returns the indices of point species s on vesicle of type v with unique index
        vesicle_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the point species
        :type s: str


        :rtype: List[steps.index_t]
        """
        return [_val1.get() for _val1 in self.ptrx().getSingleVesicleSurfaceSpecIndices(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s))]

    def getAllVesicleIndices(self):
        """
        Returns all the unique indices of all vesicles in the simulation


        :rtype: List[steps.index_t]
        """
        return [_val1.get() for _val1 in self.ptrx().getAllVesicleIndices()]

    def getAllVesicleIndicesOnPath(self):
        """
        Returns all the unique indices of all vesicles in the simulation currently on a path


        :rtype: List[steps.index_t]
        """
        return [_val1.get() for _val1 in self.ptrx().getAllVesicleIndicesOnPath()]

    def getCompVesicleIndices(self, str c, str v):
        """
        Returns all the unique indices of vesicles v currently in compartment c

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str


        :rtype: List[steps.index_t]
        """
        return [_val1.get() for _val1 in self.ptrx().getCompVesicleIndices(to_std_string(c), to_std_string(v))]

    def getSingleVesicleCompartment(self, str v, steps.index_t ves_unique_index):
        """
        Get the compartment of vesicle of type v and unique index ves_unique_index

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t


        :rtype: str
        """
        return from_std_string(self.ptrx().getSingleVesicleCompartment(to_std_string(v), steps.vesicle_individual_id(ves_unique_index)))

    def getSingleVesicleCompartment(self, str v, steps.index_t ves_unique_index):
        """
        Get the compartment of vesicle of type v and unique index ves_unique_index

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t


        :rtype: str
        """
        return from_std_string(self.ptrx().getSingleVesicleCompartment(to_std_string(v), steps.vesicle_individual_id(ves_unique_index)))

    def getSingleVesiclePos(self, str v, steps.index_t ves_unique_index):
        """
        Returns the position of vesicle of type v with unique index
        ves_unique_index(in cartesian coordinates)

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t


        :rtype: List[float]
        """
        return [_val1 for _val1 in self.ptrx().getSingleVesiclePos(to_std_string(v), steps.vesicle_individual_id(ves_unique_index))]

    def setSingleVesiclePos(self, str v, steps.index_t ves_unique_index, std.vector[double] pos, bool force=False):
        """
        Set the position of vesicle of type v with unique index
        ves_unique_index to pos (cartesian coordinates)

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param pos: Position of the vesicle in cartesian corrdinates
        :type pos: List[float]
        :param force:
        :type force: bool

        :rtype: None
        """
        self.ptrx().setSingleVesiclePos(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), pos, force)

    # TODO Remove this method eventually and use setSingleVesiclePos instead
    def setCompSingleVesiclePos(self, str c, str v, steps.index_t ves_unique_index, std.vector[double] pos, bool force=False):
        """
        Set the position of vesicle of type v with unique index
        ves_unique_index in compartment c to pos (cartesian coordinates)

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param pos: Position of the vesicle in cartesian corrdinates
        :type pos: List[float]
        :param force:
        :type force: bool


        :rtype: None
        """
        self.ptrx().setSingleVesiclePos(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), pos, force)

    def getCompVesicleSurfaceSpecCountDict(self, str c, str v, str s):
        """
        Returns the surface count of species s on vesicles v in compartment c.
        Return is a map, vesicle_unique_index : count of s

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str
        :param s: Name of the species
        :type s: str


        :rtype: Dict[steps.index_t, int]
        """
        return { _pair1.first.get(): _pair1.second for _pair1 in self.ptrx().getCompVesicleSurfaceSpecCountDict(to_std_string(c), to_std_string(v), to_std_string(s)) }

    def getCompVesicleSurfaceSpecCount(self, str c, str v, str s):
        """
        Returns the summed surface count of species s on vesicles v in compartment c.

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str
        :param s: Name of the species
        :type s: str


        :rtype: int
        """
        return self.ptrx().getCompVesicleSurfaceSpecCount(to_std_string(c), to_std_string(v), to_std_string(s))

    def getCompVesicleInnerSpecCount(self, str c, str v, str s):
        """
        Returns the summed inner count of species s on vesicles v in compartment c.

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str
        :param s: Name of the species
        :type s: str


        :rtype: int
        """
        return self.ptrx().getCompVesicleInnerSpecCount(to_std_string(c), to_std_string(v), to_std_string(s))

    def getSingleVesicleSurfaceSpecCount(self, str v, steps.index_t ves_unique_index, str s):
        """
        Get the surface count of species s on vesicle of type v and unique
        index ves_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str


        :rtype: int
        """
        return self.ptrx().getSingleVesicleSurfaceSpecCount(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s))

    def getSingleVesicleInnerSpecCount(self, str v, steps.index_t ves_unique_index, str s):
        """
        Get the inner count of species s on vesicle of type v and unique
        index ves_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str


        :rtype: int
        """
        return self.ptrx().getSingleVesicleInnerSpecCount(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s))

    def setSingleVesicleSurfaceSpecCount(self, str v, steps.index_t ves_unique_index, str s, uint count):
        """
        Set the surface count of species s on vesicle of type v and unique
        index ves_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str
        :param count: The number of the surface species
        :type count: int


        :rtype: None
        """
        self.ptrx().setSingleVesicleSurfaceSpecCount(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s), count)

    def getSingleVesicleSurfaceSpecPos(self, str v, steps.index_t ves_unique_index, str s):
        """
        Get the cartesian coordinates of species s on vesicle of type v and
        unique index ves_unique_index. Position is absolute,

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str


        :rtype: List[List[float]]
        """
        return [[_val2 for _val2 in _val1] for _val1 in self.ptrx().getSingleVesicleSurfaceSpecPos(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s))]

    def setSingleVesicleInnerSpecCount(self, str v, steps.index_t ves_unique_index, str s, uint count):
        """
        Set the count of the inner species, that is inside the vesicle volume,
        on vesicle of type v and unique index ves_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str
        :param count: The number of the inner species
        :type count: int


        :rtype: None
        """
        self.ptrx().setSingleVesicleInnerSpecCount(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s), count)

    def getSingleVesicleSurfaceSpecPosSpherical(self, str v, steps.index_t ves_unique_index, str s):
        """
        Get the spherical coordinates of species s on vesicle of type v and
        unique index ves_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str

        :rtype: List[List[float]]
        """
        return [list(pos) for pos in self.ptrx().getSingleVesicleSurfaceSpecPosSpherical(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s))]

    def setSingleVesicleSurfaceSpecPosSpherical(self, str v, steps.index_t ves_unique_index, str s, list pos_spherical):
        """
        Set the spherical coordinates of species s on vesicle of type v and
        unique index ves_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str
        :param pos_spherical: Position of the molecule in spherical coordinates (relative to the vesicle)
        :type pos_spherical: List[float]


        :rtype: None
        """
        cdef std.vector[double] pos
        cdef std.vector[std.vector[double]] allPos
        if len(pos_spherical) == 2 and all(isinstance(val, (float, int)) for val in pos_spherical):
            for val in pos_spherical:
                pos.push_back(val)
            for i in range(self.getSingleVesicleSurfaceSpecCount(v, ves_unique_index, s)):
                allPos.push_back(pos)
        else:
            for ppos in pos_spherical:
                pos.clear()
                for val in ppos:
                    pos.push_back(val)
                allPos.push_back(pos)
        self.ptrx().setSingleVesicleSurfaceSpecPosSpherical(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s), allPos)

    def getSingleSpecPosSpherical(self, str s, steps.index_t spec_unique_index):
        """
        Get the spherical coordinates of point species s with unique index spec_unique_index.

        :param s: Name of the species
        :type s: str
        :param spec_unique_index: Unique index of the individual point species
        :type spec_unique_index: steps.index_t

        :rtype: List[float]
        """
        return self.ptrx().getSingleSpecPosSpherical(to_std_string(s), steps.pointspec_individual_id(spec_unique_index))

    def getCompVesicleSurfaceLinkSpecCountDict(self, str c, str v, str ls):
        """
        Returns the count of 'link' species ls on vesicles v in compartment c.
        Return is a map, vesicle_unique_index : count of ls

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str
        :param ls: Name of the link species
        :type ls: str


        :rtype: Dict[steps.index_t, int]
        """
        return { _pair1.first.get(): _pair1.second for _pair1 in self.ptrx().getCompVesicleSurfaceLinkSpecCountDict(to_std_string(c), to_std_string(v), to_std_string(ls)) }

    def getCompVesicleSurfaceLinkSpecCount(self, str c, str v, str ls):
        """
        Returns the summed count of 'link' species ls on vesicles v in compartment c.

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str
        :param ls: Name of the link species
        :type ls: str


        :rtype: int
        """
        return self.ptrx().getCompVesicleSurfaceLinkSpecCount(to_std_string(c), to_std_string(v), to_std_string(ls))

    def getSingleVesicleSurfaceLinkSpecPos(self, str v, steps.index_t ves_unique_index, str ls):
        """
        Returns the positions of 'link' species ls on vesicle of type v and
        unique index ves_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param ls: Name of the link species
        :type ls: str


        :rtype: List[List[float]]
        """
        return [[_val2 for _val2 in _val1] for _val1 in self.ptrx().getSingleVesicleSurfaceLinkSpecPos(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(ls))]

    def getSingleLinkSpecPos(self, steps.index_t ls_unique_index):
        """
        Returns the position of the link species with unique index ls_unique_id.

        :param ls_unique_index: Unique index of the individual link species
        :type ls_unique_index: steps.index_t

        :rtype: List[float]
        """
        return [_val1 for _val1 in self.ptrx().getSingleLinkSpecPos(steps.linkspec_individual_id(ls_unique_index))]

    def getSingleLinkSpecLinkedTo(self, steps.index_t ls_unique_index):
        """
        Returns the unique index of the link species linked to the link species
        with unique index ls_unique_id.

        :param ls_unique_index: Unique index of the individual link species
        :type ls_unique_index: steps.index_t

        :rtype: steps.index_t
        """
        return self.ptrx().getSingleLinkSpecLinkedTo(steps.linkspec_individual_id(ls_unique_index)).get()

    def getSingleLinkSpecVes(self, steps.index_t ls_unique_index):
        """
        Returns the unique index of the vesicle that contains the link species
        with unique index ls_unique_id or UNDEFINED_VESICLE if the link species no longer exists.

        :param ls_unique_index: Unique index of the individual link species
        :type ls_unique_index: steps.index_t

        :rtype: steps.index_t
        """
        return self.ptrx().getSingleLinkSpecVes(steps.linkspec_individual_id(ls_unique_index)).get()

    def getSingleVesicleImmobility(self, str v, steps.index_t ves_unique_index):
        """
        Get the 'immobility' of vesicle of type v and unique index
        ves_unique_index. All non-zero numbers mean vesicle is
        immobile.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t


        :rtype: int
        """
        return self.ptrx().getSingleVesicleImmobility(to_std_string(v), steps.vesicle_individual_id(ves_unique_index))

    def getSingleVesicleOverlapTets(self, str v, steps.index_t ves_unique_index):
        """
        Get the indexes of tetrahedrons that overlap vesicle of type v and
        unique index ves_unique_index

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t


        :rtype: List[steps.index_t]
        """
        return [_val1.get() for _val1 in self.ptrx().getSingleVesicleOverlapTets(to_std_string(v), steps.vesicle_individual_id(ves_unique_index))]

    def setTetVesicleDcst(self, steps.index_t idx, str v, double dcst):
        """
        Set the diffusion rate per tetrahedron of vesicles of type v. Vesicles
        will use this diffusion rate when vesicle centre is in this tet.

        :param idx: Tetrahrdron index
        :type idx: steps.index_t
        :param v: Name of the vesicle.
        :type v: str
        :param dcst: Diffusion coefficient
        :type dcst: float


        :rtype: None
        """
        self.ptrx().setTetVesicleDcst(steps.tetrahedron_global_id(idx), to_std_string(v), dcst)

    def setVesicleSurfaceLinkSpecSDiffD(self, str v, str ls, double dcst):
        """
        Set the diffusion rate of 'link' species ls on vesicles of type v.

        :param v: Name of the vesicle.
        :type v: str
        :param ls: Name of the link species
        :type ls: str
        :param dcst: Diffusion coefficient
        :type dcst: float


        :rtype: None
        """
        self.ptrx().setVesicleSurfaceLinkSpecSDiffD(to_std_string(v), to_std_string(ls), dcst)

    def setVesSReacK(self, str vsr, double kf):
        """

        :param vsr: Name of the vesicle surface reaction
        :type vsr: str
        :param kf: Rate
        :type kf: float


        :rtype: None
        """
        self.ptrx().setVesSReacK(to_std_string(vsr), kf)

    def getVesSReacExtent(self, str vsr):
        """
        Get the reaction extent of vesicle surface reaction. Not
        compartment-specific.

        :param vsr: Name of the vesicle surface reaction
        :type vsr: str


        :rtype: int
        """
        return self.ptrx().getVesSReacExtent(to_std_string(vsr))

    def setExocytosisK(self, str exo, double kf):
        """
        Set rate of exocytosis events. Not compartment-specific.

        :param exo: Name of the vesicle exocytosis
        :type exo: str
        :param kf: Rate
        :type kf: float


        :rtype: None
        """
        self.ptrx().setExocytosisK(to_std_string(exo), kf)

    def getExocytosisExtent(self, str exo):
        """
        Get the extent of vesicle exocytosis. Not compartment-specific.

        :param exo: Name of the vesicle exocytosis
        :type exo: str


        :rtype: int
        """
        return self.ptrx().getExocytosisExtent(to_std_string(exo))

    def getExocytosisEvents(self, str exo):
        """
        Get the vesicle exocytosis events that happened since last call. Not compartment-specific.

        Returns a list of tuple with each tuple containing:
        (time of event, individual index of vesicle, triangle index where the exocytosis happened, individual index of raft or None)

        :param exo: Name of the vesicle exocytosis
        :type exo: str

        :rtype: List[Tuple[float, steps.index_t, steps.index_t, steps.index_t]]
        """
        return [(event.time, event.vidx.get(), event.tidx.get(), event.ridx.get() if event.ridx.valid() else None) for event in self.ptrx().getExocytosisEvents(to_std_string(exo))]

    def getRaftEndocytosisExtent(self, str rendo):
        """
        Get the extent of raft endocytosis. Not patch-specific.

        :param rendo: Name of the raft endocytosis
        :type rendo: str


        :rtype: int
        """
        return self.ptrx().getRaftEndocytosisExtent(to_std_string(rendo))

    def getRaftEndocytosisEvents(self, str rendo):
        """
        Get the raft endocytosis events that happened since last call. Not patch-specific.

        Returns a list of tuple with each tuple containing:
        (time of event, individual index of raft, triangle index where the endocytosis happened, individual index of vesicle)

        :param rendo: Name of the raft endocytosis
        :type rendo: str

        :rtype: List[Tuple[float, steps.index_t, steps.index_t, steps.index_t]]
        """
        return [(event.time, event.ridx.get(), event.tidx.get(), event.vidx.get()) for event in self.ptrx().getRaftEndocytosisEvents(to_std_string(rendo))]

    def setRaftEndocytosisK(self, str rendo, double kcst):
        """
        Set the rate of raft endocytosis. Not patch-specific.

        :param rendo: Name of the raft endocytosis
        :type rendo: str
        :param kcst: Rate of the raft endocytosis
        :type kcst: double

        """
        self.ptrx().setRaftEndocytosisK(to_std_string(rendo), kcst)

    def addVesicleDiffusionGroup(self, str v, list comps):
        """
        Add a 'diffusion group' for vesicles of type v. Vesicles will diffuse
        freely amongst this group (if they border each other).

        :param v: Name of the vesicle.
        :type v: str
        :param comps: List of compartment names.
        :type comps: List[str]


        :rtype: None
        """
        cdef std.vector[std.string] std_comps
        std_comps.reserve(len(comps))
        for _elem in comps:
            std_comps.push_back(to_std_string(_elem))
        self.ptrx().addVesicleDiffusionGroup(to_std_string(v), std_comps)

    def addPathVesicle(self, str path, str ves, double speed, dict spec_deps={}, stoch_stepsize=1e-9):
        """
        Add a vesicle to this path. This means a vesicle of this type can
        interact with this path upon overlapping it

        :param path: Name of the path
        :type path: str
        :param ves: Name of the vesicle
        :type ves: str
        :param speed: Speed of the vesicle on this path in m/s
        :type speed: float
        :param spec_deps: Optional species dependencies, in map of species names to number of the species required
        :type spec_deps: Dict[str, int]
        :param stoch_stepsize: Stochastic step length. This may be a single float value where a single-exponential will be applied. If a list of length 2, a double-exponential will be applied by the proportionality specified in the 2nd element.
        :type stoch_stepsize: Float or list[float]

        :rtype: None
        """
        cdef std.map[std.string, uint] std_spec_deps
        for _elem1, _elem2 in spec_deps.items():
            std_spec_deps[to_std_string(_elem1)] = _elem2

        cdef std.vector[double] std_stoch_stepsize
        if not hasattr(stoch_stepsize, '__iter__'):
            stoch_stepsize = [stoch_stepsize]
        std_stoch_stepsize.reserve(len(stoch_stepsize))
        for _elem in stoch_stepsize:
            std_stoch_stepsize.push_back(_elem)
        self.ptrx().addPathVesicle(to_std_string(path), to_std_string(ves), speed, std_spec_deps, std_stoch_stepsize)

    def getTetVol(self, steps.index_t idx):
        """
        Returns the volume of a tetrahedron (in m^3).

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTetVol(steps.tetrahedron_global_id(idx))

    def getTetReducedVol(self, steps.index_t idx):
        """
        Returns the reduced volume of a tetrahedron (in m^3).

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTetReducedVol(steps.tetrahedron_global_id(idx))

    def setTetVol(self, steps.index_t idx, double vol):
        """
        Set the volume of a tetrahedron (in m^3).

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param vol: Volume of the tetrahedron.
        :type vol: float


        :rtype: None
        """
        self.ptrx().setTetVol(steps.tetrahedron_global_id(idx), vol)

    def getTetSpecDefined(self, steps.index_t idx, str s):
        """
        Returns whether species s is defined in a tetrahedron

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getTetSpecDefined(steps.tetrahedron_global_id(idx), to_std_string(s))

    def getTetSpecCount(self, steps.index_t idx, str s):
        """
        Returns the number of molecules of species s in a tetrahedron

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getTetSpecCount(steps.tetrahedron_global_id(idx), to_std_string(s))

    def setTetSpecCount(self, steps.index_t idx, str s, double n):
        """
        Sets the number of molecules of species s in a tetrahedron

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str
        :param n: Number of molecules of the species.
        :type n: float


        :rtype: None
        """
        self.ptrx().setTetSpecCount(steps.tetrahedron_global_id(idx), to_std_string(s), n)

    def getTetSpecAmount(self, steps.index_t idx, str s):
        """
        Returns the amount (in mols) of species s in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getTetSpecAmount(steps.tetrahedron_global_id(idx), to_std_string(s))

    def setTetSpecAmount(self, steps.index_t idx, str s, double m):
        """
        Sets the amount (in mols) of species s in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str
        :param m: Amount of the species.
        :type m: float


        :rtype: None
        """
        self.ptrx().setTetSpecAmount(steps.tetrahedron_global_id(idx), to_std_string(s), m)

    def getTetSpecConc(self, steps.index_t idx, str s):
        """
        Returns the concentration (in molar units) of species s in a
        tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getTetSpecConc(steps.tetrahedron_global_id(idx), to_std_string(s))

    def setTetSpecConc(self, steps.index_t idx, str s, double c):
        """
        Sets the concentration (in molar units) of species s in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str
        :param c: Concentration of the species.
        :type c: float


        :rtype: None
        """
        self.ptrx().setTetSpecConc(steps.tetrahedron_global_id(idx), to_std_string(s), c)

    def getTetSpecClamped(self, steps.index_t idx, str s):
        """
        Returns whether the concentration of species s in a tetrahedron
        remains constant over time (unless changed explicitly).

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getTetSpecClamped(steps.tetrahedron_global_id(idx), to_std_string(s))

    def setTetSpecClamped(self, steps.index_t idx, str s, bool buf):
        """
        Sets clamping of species s in a tetrahedron on or off.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str
        :param buf: Flag to turn the clamping of species on or off.
        :type buf: bool


        :rtype: None
        """
        self.ptrx().setTetSpecClamped(steps.tetrahedron_global_id(idx), to_std_string(s), buf)

    def getTetReacK(self, steps.index_t idx, str r):
        """
        Returns the macroscopic reaction constant of reaction r in a
        tetrahedron (units vary with order of reaction).

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTetReacK(steps.tetrahedron_global_id(idx), to_std_string(r))

    def setTetReacK(self, steps.index_t idx, str r, double kf):
        """
        Sets the macroscopic reaction constant of reaction r in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param r: Name of the reaction.
        :type r: str
        :param kf: Rate constant of the reaction.
        :type kf: float


        :rtype: None
        """
        self.ptrx().setTetReacK(steps.tetrahedron_global_id(idx), to_std_string(r), kf)

    def getTetReacActive(self, steps.index_t idx, str r):
        """
        Returns whether reaction r in a tetrahedron is active or not

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param r: Name of the reaction.
        :type r: str


        :rtype: bool
        """
        return self.ptrx().getTetReacActive(steps.tetrahedron_global_id(idx), to_std_string(r))

    def setTetReacActive(self, steps.index_t idx, str r, bool act):
        """
        Activates/deactivates reaction r in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param r: Name of the reaction.
        :type r: str
        :param act: Flag to activate or deactivate the reaction.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setTetReacActive(steps.tetrahedron_global_id(idx), to_std_string(r), act)

    def getTetDiffD(self, steps.index_t idx, str d, steps.index_t direction_tet=UNKNOWN_TET):
        """
        Returns the diffusion constant of diffusion rule d in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param d: Name of the deffusion.
        :type d: str
        :param direction_tet: Tetrahedron index which specifies diffusion direction.
        :type direction_tet: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTetDiffD(steps.tetrahedron_global_id(idx), to_std_string(d), steps.tetrahedron_global_id(direction_tet))

    def setTetDiffD(self, steps.index_t idx, str d, double dk, steps.index_t direction_tet=UNKNOWN_TET):
        """
        Sets the diffusion constant of diffusion rule d in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param d: Name of the diffusion.
        :type d: str
        :param dk: Rate constant of the diffusion.
        :type dk: float
        :param direction_tet: Tetrahedron index which the diffusion towards.
        :type direction_tet: steps.index_t


        :rtype: None
        """
        self.ptrx().setTetDiffD(steps.tetrahedron_global_id(idx), to_std_string(d), dk, steps.tetrahedron_global_id(direction_tet))

    def getTetDiffActive(self, steps.index_t idx, str d):
        """
        Returns whether diffusion rule d in a tetrahedron is active or not.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param d: Name of the diffusion.
        :type d: str


        :rtype: bool
        """
        return self.ptrx().getTetDiffActive(steps.tetrahedron_global_id(idx), to_std_string(d))

    def setTetDiffActive(self, steps.index_t idx, str d, bool act):
        """
        Activates/deactivates diffusion rule d in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param d: Name of the diffusion.
        :type d: str
        :param act: Flag to activate / deactivate the diffusion.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setTetDiffActive(steps.tetrahedron_global_id(idx), to_std_string(d), act)

    def getTetReacC(self, steps.index_t idx, str r):
        """
        Returns c_mu, the mesoscopic reaction constant of reaction r in
        a tetrahedron

        :param idx: Index of the diffusion.
        :type idx: steps.index_t
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTetReacC(steps.tetrahedron_global_id(idx), to_std_string(r))

    def getTetReacH(self, steps.index_t idx, str r):
        """
        Returns h_mu, the distinct number of ways in which reaction r can
        occur in a tetrahedron, by computing the product of its reactants.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param r:
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTetReacH(steps.tetrahedron_global_id(idx), to_std_string(r))

    def getTetReacA(self, steps.index_t idx, str r):
        """
        Returns the propensity, a_mu, of reaction r in a tetrahedron.
        The propensity value gives the probability per unit time that this
        reaction will occur in the current state.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTetReacA(steps.tetrahedron_global_id(idx), to_std_string(r))

    def getTetDiffA(self, steps.index_t idx, str d):
        """
        Returns the propensity, a_mu of diffusion rule d in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param d: Name of the diffusion.
        :type d: str


        :rtype: float
        """
        return self.ptrx().getTetDiffA(steps.tetrahedron_global_id(idx), to_std_string(d))

    def getTetV(self, steps.index_t idx):
        """
        Returns the potential of tetrahedron in Volts.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTetV(steps.tetrahedron_global_id(idx))

    def setTetV(self, steps.index_t idx, double v):
        """
        Set the potential of tetrahedron.

        :param idx: Index of the tetrahedron
        :type idx: steps.index_t
        :param v:
        :type v: float


        :rtype: None
        """
        self.ptrx().setTetV(steps.tetrahedron_global_id(idx), v)

    def getTetVClamped(self, steps.index_t idx):
        """
        Returns whether the potential of tetrahedron is clamped over time
        (unless changed explicitly)

        :param idx: Index of the tetrahedron
        :type idx: steps.index_t


        :rtype: bool
        """
        return self.ptrx().getTetVClamped(steps.tetrahedron_global_id(idx))

    def setTetVClamped(self, steps.index_t idx, bool cl):
        """
        Sets voltage clamp in tetrahedron.

        :param idx: Index of the tetrahedron
        :type idx: steps.index_t
        :param cl: Flag to turn the clamping on or off.
        :type cl: bool


        :rtype: None
        """
        self.ptrx().setTetVClamped(steps.tetrahedron_global_id(idx), cl)

    def getPatchArea(self, str p):
        """
        Returns the area of patch p (in m^2)

        :param p: Name of the patch.
        :type p: str


        :rtype: float
        """
        return self.ptrx().getPatchArea(to_std_string(p))

    def getPatchSpecCount(self, str p, str s):
        """
        Returns the number of molecules of species s in patch p.
        NOTE: in a mesh-based simulation, the total count is computed as
        the sum of the counts in all triangles of the patch.

        :param p: Name of the path.
        :type p: str
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getPatchSpecCount(to_std_string(p), to_std_string(s))

    def setPatchSpecCount(self, str p, str s, double n):
        """
        Sets the number of molecules of species s in patch p.
        NOTE: in a mesh-based simulation, the total amount is equally divided
        over all triangles in the patch.

        :param p: Name of the patch.
        :type p: str
        :param s: Name of the species.
        :type s: str
        :param n: Number of molecules of species.
        :type n: float


        :rtype: None
        """
        self.ptrx().setPatchSpecCount(to_std_string(p), to_std_string(s), n)

    def getPatchSpecAmount(self, str p, str s):
        """
        Returns the amount (in mols) of species s in patch p.
        NOTE: in a mesh-based simulation, the total amount is computed as
        the sum of the amounts in all triangles of the patch.

        :param p: Name of the patch.
        :type p: str
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getPatchSpecAmount(to_std_string(p), to_std_string(s))

    def setPatchSpecAmount(self, str p, str s, double a):
        """
        Sets the amount (in mols) of species s in patch p.
        NOTE: in a mesh-based simulation, the total amount is equally divided
        over all triangles in the patch.

        :param p: Name of the patch.
        :type p: str
        :param s: Name of the species.
        :type s: str
        :param a: Amount of the species.
        :type a: float


        :rtype: None
        """
        self.ptrx().setPatchSpecAmount(to_std_string(p), to_std_string(s), a)

    def getPatchSpecClamped(self, str p, str s):
        """
        Returns whether the count of species s in patch p remains constant.
        over time (unless changed explicitly).
        NOTE: this method will only return true if the species has been
        clamped in all triangles in the patch.

        :param p: Name of the patch.
        :type p: str
        :param s: Name of the species.
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getPatchSpecClamped(to_std_string(p), to_std_string(s))

    def setPatchSpecClamped(self, str p, str s, bool buf):
        """
        Turns clamping of species in patch p on or off.
        NOTE: in a mesh-based simulation, this method turns clamping on/off
        in all triangles in the patch.

        :param p: Name of the patch.
        :type p: str
        :param s: Name of the species.
        :type s: str
        :param buf: Flag to turn clamping of species on /off.
        :type buf: bool


        :rtype: None
        """
        self.ptrx().setPatchSpecClamped(to_std_string(p), to_std_string(s), buf)

    def getPatchSReacK(self, str p, str r):
        """
        Returns the macroscopic reaction constant of surface reaction r
        in patch p.
        NOTE: in a mesh-based simulation, the value is computed as the
        area-weighted sum of the reaction constants in all triangles of
        the patch.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getPatchSReacK(to_std_string(p), to_std_string(r))

    def setPatchSReacK(self, str p, str r, double kf):
        """
        Sets the macroscopic reaction constant of surface reaction r
        in patch p.
        NOTE: in a mesh-based simulation this method changes the reaction
        constant equally in all triangles of the patch.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str
        :param kf: Rate constant of the reaction.
        :type kf: float


        :rtype: None
        """
        self.ptrx().setPatchSReacK(to_std_string(p), to_std_string(r), kf)

    def getPatchSReacActive(self, str p, str r):
        """
        Returns whether surface reaction r in patch p is active or not.
        NOTE: in a mesh-based simulation, only returns false when the
        reaction has been inactivated in all triangles.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: bool
        """
        return self.ptrx().getPatchSReacActive(to_std_string(p), to_std_string(r))

    def setPatchSReacActive(self, str p, str r, bool a):
        """
        Activate or inactivate surface reaction r in patch p.
        NOTE: in a mesh-based simulation, activation/inactivation of a
        surface reaction turns it on/off in all triangles.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str
        :param a: Flag to activate / deactivate the reaction.
        :type a: bool


        :rtype: None
        """
        self.ptrx().setPatchSReacActive(to_std_string(p), to_std_string(r), a)

    def getPatchSReacC(self, str p, str r):
        """
        Returns c_mu, the mesoscopic reaction constant of surface reaction r
        in patch p.
        NOTE: in a mesh_based simulation, the mesoscopic reaction constant
        is computed as the sum of the mesoscopic reaction constants from all
        triangles of the patch.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reacton.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getPatchSReacC(to_std_string(p), to_std_string(r))

    def getPatchSReacH(self, str p, str r):
        """
        Returns h_mu, the distinct number of ways in which a surface reaction
        r can occur in patch p, by computing the product of its reactants.
        NOTE: in a mesh-based simulation, it returns the sum of the h_mu's
        over all triangles triangles of the patch.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getPatchSReacH(to_std_string(p), to_std_string(r))

    def getPatchSReacA(self, str p, str r):
        """
        Returns the propensity, a_mu of surface reaction r in patch p.
        This propensity value gives the probability per unit time that this
        surface reaction will occur in the current state.
        NOTE: in a mesh-based simulation, a_mu is computed as the sum of the
        a_mu in all triangles in the patch.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getPatchSReacA(to_std_string(p), to_std_string(r))

    def getPatchSReacExtent(self, str p, str r):
        """
        Returns the extent of surface reaction r in patch p.
        NOTE: in a mesh-based simulation, returns the sum of the extents in
        all triangles of the patch.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: unsigned long long
        """
        return self.ptrx().getPatchSReacExtent(to_std_string(p), to_std_string(r))

    def resetPatchSReacExtent(self, str p, str r):
        """
        Resets the extent of surface reaction r in patch p to zero.
        NOTE: in a mesh-based simulation, resets the extents of the
        surface reaction in all triangles of the patch.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: None
        """
        self.ptrx().resetPatchSReacExtent(to_std_string(p), to_std_string(r))

    def getPatchVDepSReacActive(self, str p, str vsr):
        """
        Returns whether voltage-dependent surface reaction vsr in patch p is
        active or not.
        NOTE: only returns false when the voltage-dependent surface
        reaction has been inactivated in all triangles.

        :param p: Name of the patch.
        :type p: str
        :param vsr: Name of the voltage-dependent surface reaction.
        :type vsr: str


        :rtype: bool
        """
        return self.ptrx().getPatchVDepSReacActive(to_std_string(p), to_std_string(vsr))

    def setPatchVDepSReacActive(self, str p, str vsr, bool a):
        """
        Activate or inactivate voltage-dependent surface reaction vsr in patch
        p.
        NOTE: activation/inactivation of a voltage-dependent
        surface reaction turns it on/off in all triangles.

        :param p: Name of the patch.
        :type p: str
        :param vsr: Name of the voltage-dependent surface reaction.
        :type vsr: str
        :param a: Flag to activate / deactivate the reaction.
        :type a: bool


        :rtype: None
        """
        self.ptrx().setPatchVDepSReacActive(to_std_string(p), to_std_string(vsr), a)

    def getPatchRaftCount(self, str p, str r):
        """
        Returns the number of rafts r in patch p

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the raft.
        :type r: str


        :rtype: int
        """
        return self.ptrx().getPatchRaftCount(to_std_string(p), to_std_string(r))

    def setPatchRaftCount(self, str p, str r, uint n):
        """
        Set the number of rafts r in patch p

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the raft.
        :type r: str
        :param n: Number of rafts
        :type n: int


        :rtype: None
        """
        self.ptrx().setPatchRaftCount(to_std_string(p), to_std_string(r), n)

    def getSingleRaftPos(self, str r, steps.index_t raft_unique_index):
        """
        Returns the position of raft of type r with unique index
        raft_unique_index (in cartesian coordinates)

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t


        :rtype: List[float]
        """
        return [_val1 for _val1 in self.ptrx().getSingleRaftPos(to_std_string(r), steps.raft_individual_id(raft_unique_index))]

    def getPatchRaftSpecCountDict(self, str p, str r, str s):
        """
        Returns the count of species s in raft r in patch p. Return is a map,
        raft_unique_index : count of s

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the raft.
        :type r: str
        :param s: Name of the species
        :type s: str


        :rtype: Dict[steps.index_t, int]
        """
        return { _pair1.first.get(): _pair1.second for _pair1 in self.ptrx().getPatchRaftSpecCountDict(to_std_string(p), to_std_string(r), to_std_string(s)) }

    def getPatchRaftSpecCount(self, str p, str r, str s):
        """
        Returns the summed count of species s in raft r in patch p.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the raft.
        :type r: str
        :param s: Name of the species
        :type s: str


        :rtype: int
        """
        return self.ptrx().getPatchRaftSpecCount(to_std_string(p), to_std_string(r), to_std_string(s))

    def getSingleRaftSpecCount(self, str r, steps.index_t raft_unique_index, str s):
        """
        Get the count of species s on raft of type r and unique index
        raft_unique_index.

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str


        :rtype: int
        """
        return self.ptrx().getSingleRaftSpecCount(to_std_string(r), steps.raft_individual_id(raft_unique_index), to_std_string(s))

    def setSingleRaftSpecCount(self, str r, steps.index_t raft_unique_index, str s, uint count):
        """
        Set the count of species s on raft of type r and unique index
        raft_unique_index.

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str
        :param count: The number of the species
        :type count: int

        :rtype: None
        """
        self.ptrx().setSingleRaftSpecCount(to_std_string(r), steps.raft_individual_id(raft_unique_index), to_std_string(s), count)

    def getSingleRaftImmobility(self, str r, steps.index_t raft_unique_index):
        """
        Get the 'immobility' of raft of type r and unique index
        raft_unique_index. All non-zero numbers mean raft is
        immobile.

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t


        :rtype: int
        """
        return self.ptrx().getSingleRaftImmobility(to_std_string(r), steps.raft_individual_id(raft_unique_index))

    def getSingleRaftRaftEndocytosisK(self, str r, steps.index_t raft_unique_index, str rendo):
        """
        Get the raft endocytosis rate of a raft

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t
        :param rendo: Name of the raft endocytosis
        :type rendo: str


        :rtype: float
        """
        return self.ptrx().getSingleRaftRaftEndocytosisK(to_std_string(r), steps.raft_individual_id(raft_unique_index), to_std_string(rendo))

    def setSingleRaftRaftEndocytosisK(self, str r, steps.index_t raft_unique_index, str rendo, double k):
        """
        Set the raft endocytosis rate of a raft

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t
        :param rendo: Name of the raft endocytosis
        :type rendo: str
        :param k: Rate.
        :type k: float


        :rtype: None
        """
        self.ptrx().setSingleRaftRaftEndocytosisK(to_std_string(r), steps.raft_individual_id(raft_unique_index), to_std_string(rendo), k)


    def setSingleRaftSReacActive(self, str r, steps.index_t raft_unique_index, str rsreac, bool active):
        """
        Activate or de-activate a raft surface reaction on a raft

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t
        :param rsreac: Name of the raft surface reaction
        :type rsreac: str
        :param active: Whether the raftsreac should be active
        :type active: bool


        :rtype: None
        """
        self.ptrx().setSingleRaftSReacActive(to_std_string(r), steps.raft_individual_id(raft_unique_index), to_std_string(rsreac), active)


    def getSingleRaftSReacActive(self, str r, steps.index_t raft_unique_index, str rsreac):
        """
        Get whether raft surface reaction on a raft is active or not

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t
        :param rsreac: Name of the raft surface reaction
        :type rsreac: str


        :rtype: bool
        """
        return self.ptrx().getSingleRaftSReacActive(to_std_string(r), steps.raft_individual_id(raft_unique_index), to_std_string(rsreac))


    def setPatchEndocyticZoneEndocytosisActive(self, str patch, str zone, str endo, bool active):
        """
        Activate or de-activate endocytosis reaction in an endocytic zone.

        :param patch:
        :type patch: str
        :param zone: Name of the zone
        :type zone: str
        :param endo: Name of the endocytosis
        :type endo: str
        :param active: Whether the endocytosis should be active
        :type active: bool


        :rtype: None
        """
        self.ptrx().setPatchEndocyticZoneEndocytosisActive(to_std_string(patch), to_std_string(zone), to_std_string(endo), active)

    def setPatchEndocyticZoneEndocytosisK(self, str patch, str zone, str endo, double k):
        """
        Change endocytosis rate in an endocytic zone.

        :param patch:
        :type patch: str
        :param zone: Name of the zone
        :type zone: str
        :param endo: Name of the endocytosis
        :type endo: str
        :param k: Rate of the endocytosis in this endocytic zone
        :type k: double


        :rtype: None
        """
        self.ptrx().setPatchEndocyticZoneEndocytosisK(to_std_string(patch), to_std_string(zone), to_std_string(endo), k)

    def getPatchEndocyticZoneEndocytosisExtent(self, str patch, str zone, str endo):
        """
        Returns the endocytosis extent in an endocytic zone

        :param patch:
        :type patch: str
        :param zone: Name of the zone
        :type zone: str
        :param endo: Name of the endocytosis.
        :type endo: str


        :rtype: int
        """
        return self.ptrx().getPatchEndocyticZoneEndocytosisExtent(to_std_string(patch), to_std_string(zone), to_std_string(endo))

    def getPatchEndocyticZoneEndocytosisEvents(self, str patch, str zone, str endo):
        """
        Get the vesicle endocytosis events in an endocytic zone that happened since last call.

        :param patch: Name of the patch
        :type patch: str
        :param zone: Name of the zone
        :type zone: str
        :param endo: Name of the endocytosis.
        :type endo: str


        :rtype: List[Tuple[float, steps.index_t, steps.index_t]]
        """
        return [(event.time, event.tidx.get(), event.vidx.get()) for event in self.ptrx().getPatchEndocyticZoneEndocytosisEvents(to_std_string(patch), to_std_string(zone), to_std_string(endo))]

    def getPatchRaftIndices(self, str p, str r):
        """
        Returns all the unique indices of rafts r currently in patch p

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the raft.
        :type r: str


        :rtype: List[steps.index_t]
        """
        return [_val1.get() for _val1 in self.ptrx().getPatchRaftIndices(to_std_string(p), to_std_string(r))]

    def getSingleRaftPatch(self, str r, steps.index_t raft_unique_index):
        """
        Get the patch of raft of type r and unique index raft_unique_index

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft
        :type raft_unique_index: steps.index_t


        :rtype: str
        """
        return from_std_string(self.ptrx().getSingleRaftPatch(to_std_string(r), steps.raft_individual_id(raft_unique_index)))

    def setDiffBoundarySpecDiffusionActive(self, str db, str s, bool act):
        """
        Activate or inactivate diffusion across a diffusion boundary for a
        species.

        :param db: Name of the diffusion boundary.
        :type db: str
        :param s: Name of the species.
        :type s: str
        :param act: Bool to activate (true) or inactivate (false) diffusion.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setDiffBoundarySpecDiffusionActive(to_std_string(db), to_std_string(s), act)

    def getDiffBoundarySpecDiffusionActive(self, str db, str s):
        """
        Returns whether diffusion is active across a diffusion boundary for a
        species.

        :param db: Name of the diffusion boundary.
        :type db: str
        :param s: Name of the species.
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getDiffBoundarySpecDiffusionActive(to_std_string(db), to_std_string(s))

    def setDiffBoundarySpecDcst(self, str db, str s, double dcst, str direction_comp=''):
        """
        Set the diffusion constant across a diffusion boundary.

        :param db: Name of the diffusion boundary.
        :type db: str
        :param s: Name of the species.
        :type s: str
        :param dcst: diffusion constant.
        :type dcst: float
        :param direction_comp: ID of the compartment which the diffusion towards to.
        :type direction_comp: str


        :rtype: None
        """
        self.ptrx().setDiffBoundarySpecDcst(to_std_string(db), to_std_string(s), dcst, to_std_string(direction_comp))

    def setSDiffBoundarySpecDiffusionActive(self, str sdb, str s, bool act):
        """
        Activate or inactivate diffusion across a surface diffusion boundary
        for a species.

        :param sdb: Name of the surface diffusion boundary.
        :type sdb: str
        :param s: Name of the species.
        :type s: str
        :param act: Bool to activate (true) or inactivate (false) diffusion.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setSDiffBoundarySpecDiffusionActive(to_std_string(sdb), to_std_string(s), act)

    def getSDiffBoundarySpecDiffusionActive(self, str sdb, str s):
        """
        Returns whether diffusion is active across a surface diffusion boundary
        for a species.

        :param sdb: Name of the surface diffusion boundary.
        :type sdb: str
        :param s: Name of the species.
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getSDiffBoundarySpecDiffusionActive(to_std_string(sdb), to_std_string(s))

    def setSDiffBoundarySpecDcst(self, str sdb, str s, double dcst, str direction_patch=''):
        """
        Set the diffusion constant across a surface diffusion boundary.

        :param sdb: Name of the surface diffusion boundary.
        :type sdb: str
        :param s: Name of the species.
        :type s: str
        :param dcst: diffusion constant.
        :type dcst: float
        :param direction_patch: ID of the patch which the diffusion is towards to.
        :type direction_patch: str


        :rtype: None
        """
        self.ptrx().setSDiffBoundarySpecDcst(to_std_string(sdb), to_std_string(s), dcst, to_std_string(direction_patch))

    def getTriArea(self, steps.index_t idx):
        """
        Returns the area of the triangle (in m^2).

        :param idx: Index of the triangle.
        :type idx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTriArea(steps.triangle_global_id(idx))

    def setTriArea(self, steps.index_t idx, double area):
        """
        Set the area (in m^2) of the triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param area: Area of teh triangle.
        :type area: float


        :rtype: None
        """
        self.ptrx().setTriArea(steps.triangle_global_id(idx), area)

    def getTriSpecDefined(self, steps.index_t idx, str s):
        """
        Returns whether species s is defined in a triangle

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getTriSpecDefined(steps.triangle_global_id(idx), to_std_string(s))

    def getTriSpecCount(self, steps.index_t idx, str s):
        """
        Returns the number of molecules of species s in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getTriSpecCount(steps.triangle_global_id(idx), to_std_string(s))

    def setTriSpecCount(self, steps.index_t idx, str s, double n):
        """
        Sets the number of molecules of species s in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str
        :param n: Number of molecules of the species.
        :type n: float


        :rtype: None
        """
        self.ptrx().setTriSpecCount(steps.triangle_global_id(idx), to_std_string(s), n)

    def getTriSpecAmount(self, steps.index_t idx, str s):
        """
        Returns the amount (in mols) of species s in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getTriSpecAmount(steps.triangle_global_id(idx), to_std_string(s))

    def setTriSpecAmount(self, steps.index_t idx, str s, double m):
        """
        Sets the amount (in mols) of species s in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str
        :param m: Amount of the species.
        :type m: float


        :rtype: None
        """
        self.ptrx().setTriSpecAmount(steps.triangle_global_id(idx), to_std_string(s), m)

    def getTriSpecClamped(self, steps.index_t idx, str s):
        """
        Returns whether the number of molecules of species s in a triangle
        remains constant over time (unless changed explicitly)

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param s: name of the species.
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getTriSpecClamped(steps.triangle_global_id(idx), to_std_string(s))

    def setTriSpecClamped(self, steps.index_t idx, str s, bool buf):
        """
        Sets clamping of species s in a triangle on or off.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param s: name of the species.
        :type s: str
        :param buf: Flag to set clamping of species on /off.
        :type buf: bool


        :rtype: None
        """
        self.ptrx().setTriSpecClamped(steps.triangle_global_id(idx), to_std_string(s), buf)

    def getTriSReacK(self, steps.index_t idx, str r):
        """
        Returns the macroscopic reaction constant of surface reaction r
        in a triangle (units vary with order of reaction).

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTriSReacK(steps.triangle_global_id(idx), to_std_string(r))

    def setTriSReacK(self, steps.index_t idx, str r, double kf):
        """
        Sets the macroscopic reaction constant of surface reaction r in
        a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: name of the reaction.
        :type r: str
        :param kf: Rate constant of the reaction.
        :type kf: float


        :rtype: None
        """
        self.ptrx().setTriSReacK(steps.triangle_global_id(idx), to_std_string(r), kf)

    def getTriSReacActive(self, steps.index_t idx, str r):
        """
        Returns whether surface reaction r in a triangle is active or not.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: name of the reaction.
        :type r: str


        :rtype: bool
        """
        return self.ptrx().getTriSReacActive(steps.triangle_global_id(idx), to_std_string(r))

    def setTriSReacActive(self, steps.index_t idx, str r, bool act):
        """
        Activates/inactivates surface reaction r in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: name of the reaction.
        :type r: str
        :param act: Flag to activate / deactivate the reaction.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setTriSReacActive(steps.triangle_global_id(idx), to_std_string(r), act)

    def getTriSReacC(self, steps.index_t idx, str r):
        """
        Returns c_mu, the mesoscopic reaction constant of surface reaction r
        in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTriSReacC(steps.triangle_global_id(idx), to_std_string(r))

    def getTriSReacH(self, steps.index_t idx, str r):
        """
        Returns h_mu, the distinct number of ways in which surface reaction r
        can occur in a triangle, by computing the product of it's reactants.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTriSReacH(steps.triangle_global_id(idx), to_std_string(r))

    def getTriSReacA(self, steps.index_t idx, str r):
        """
        Returns the propensity, a_mu, of surface reaction r in a triangle.
        The propensity value gives the probability per unit time that this
        surface reaction will occur in the current state.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTriSReacA(steps.triangle_global_id(idx), to_std_string(r))

    def getTriDiffD(self, steps.index_t idx, str d, uint direction_tri=UNKNOWN_TRI):
        """
        outdate function

        :param idx:
        :type idx: steps.index_t
        :param d:
        :type d: str
        :param direction_tri:
        :type direction_tri: int


        :rtype: float
        """
        return self.ptrx().getTriDiffD(steps.triangle_global_id(idx), to_std_string(d), direction_tri)

    def getTriSDiffD(self, steps.index_t idx, str d, steps.index_t direction_tri=UNKNOWN_TRI):
        """
        Returns the diffusion constant of diffusion rule d in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param d: Name of the diffusion.
        :type d: str
        :param direction_tri:
        :type direction_tri: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTriSDiffD(steps.triangle_global_id(idx), to_std_string(d), steps.triangle_global_id(direction_tri))

    def setTriDiffD(self, steps.index_t idx, str d, double dk, steps.index_t direction_tri=UNKNOWN_TRI):
        """
        outdated function

        :param idx:
        :type idx: steps.index_t
        :param d:
        :type d: str
        :param dk:
        :type dk: float
        :param direction_tri:
        :type direction_tri: steps.index_t


        :rtype: None
        """
        self.ptrx().setTriDiffD(steps.triangle_global_id(idx), to_std_string(d), dk, steps.triangle_global_id(direction_tri))

    def setTriSDiffD(self, steps.index_t idx, str d, double dk, steps.index_t direction_tri=UNKNOWN_TRI):
        """
        Sets the diffusion constant of diffusion rule d on a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param d: Name of the diffusion.
        :type d: str
        :param dk: Rate constant of the diffusion.
        :type dk: float
        :param direction_tri: Triangle index which the diffusion towards
        :type direction_tri: steps.index_t


        :rtype: None
        """
        self.ptrx().setTriSDiffD(steps.triangle_global_id(idx), to_std_string(d), dk, steps.triangle_global_id(direction_tri))

    def getTriRaftCount(self, steps.index_t idx, str r):
        """
        Returns the number of rafts in a triangle

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: Name of the raft.
        :type r: str


        :rtype: int
        """
        return self.ptrx().getTriRaftCount(steps.triangle_global_id(idx), to_std_string(r))

    def setTriRaftCount(self, steps.index_t idx, str r, uint n):
        """
        Returns the number of rafts in a triangle

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: Name of the raft.
        :type r: str
        :param n: Number of rafts
        :type n: int


        :rtype: None
        """
        self.ptrx().setTriRaftCount(steps.triangle_global_id(idx), to_std_string(r), n)

    def addTriRaft(self, steps.index_t idx, str r):
        """
        Add one raft to a triangle

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: Name of the raft.
        :type r: str


        :rtype: steps.index_t
        """
        return self.ptrx().addTriRaft(steps.triangle_global_id(idx), to_std_string(r)).get()

    def getTriExocytosisActive(self, steps.index_t idx, str er):
        """
        Returns whether exocytotic reaction in a triangle is active or not.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param er: name of the exocytotic reaction.
        :type er: str


        :rtype: bool
        """
        return self.ptrx().getTriExocytosisActive(steps.triangle_global_id(idx), to_std_string(er))

    def setTriExocytosisActive(self, steps.index_t idx, str er, bool act):
        """
        Activates/inactivates exocytotic reaction in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param er: name of the exocytotic reaction.
        :type er: str
        :param act: Flag to activate / deactivate the exocytotic reaction.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setTriExocytosisActive(steps.triangle_global_id(idx), to_std_string(er), act)

    def getTriV(self, steps.index_t idx):
        """
        Returns the potential of triangle in Volts.

        :param idx: Index of the triangle.
        :type idx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTriV(steps.triangle_global_id(idx))

    def setTriV(self, steps.index_t idx, double v):
        """
        Set the potential of triangle.

        :param idx: Index of the triangle
        :type idx: steps.index_t
        :param v:
        :type v: float


        :rtype: None
        """
        self.ptrx().setTriV(steps.triangle_global_id(idx), v)

    def getTriVClamped(self, steps.index_t idx):
        """
        Returns whether the potential of triangle is clamped over time
        (unless changed explicitly).

        :param idx: Index of the triangle
        :type idx: steps.index_t


        :rtype: bool
        """
        return self.ptrx().getTriVClamped(steps.triangle_global_id(idx))

    def setTriVClamped(self, steps.index_t idx, bool cl):
        """
        Sets voltage clamp in triangle.

        :param idx: Index of the triangle
        :type idx: steps.index_t
        :param cl: Flag to turn the clamping on or off.
        :type cl: bool


        :rtype: None
        """
        self.ptrx().setTriVClamped(steps.triangle_global_id(idx), cl)

    def getTriOhmicI(self, steps.index_t idx, str oc):
        """
        Returns the ohmic current of triangle in amperes.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param oc: name of the ohmic current
        :type oc: str


        :rtype: float
        """
        return self.ptrx().getTriOhmicI(steps.triangle_global_id(idx), to_std_string(oc))

    def getTriGHKI(self, steps.index_t idx, str ghk):
        """
        Returns the GHK current of triangle in amperes.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param ghk: name of the ghk current
        :type ghk: str


        :rtype: float
        """
        return self.ptrx().getTriGHKI(steps.triangle_global_id(idx), to_std_string(ghk))

    def getTriI(self, steps.index_t idx):
        """
        Returns the current of a triangle in amperes from the last EField
        calculation step.

        :param idx: Index of the triangle
        :type idx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTriI(steps.triangle_global_id(idx))

    def getTriIClamp(self, steps.index_t idx):
        """
        Gets current injection to triangle.

        :param idx: Index of the triangle
        :type idx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTriIClamp(steps.triangle_global_id(idx))

    def setTriIClamp(self, steps.index_t idx, double i):
        """
        Sets current injection to triangle.
        Will be assumed to be constant for one EField DT

        :param idx: Index of the triangle
        :type idx: steps.index_t
        :param i:
        :type i: float


        :rtype: None
        """
        self.ptrx().setTriIClamp(steps.triangle_global_id(idx), i)

    def getTriVDepSReacActive(self, steps.index_t idx, str vsr):
        """
        Returns whether voltage-dependent surface reaction vsr in a triangle is
        active or not.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param vsr: name of the voltage-dependent surface reaction.
        :type vsr: str


        :rtype: bool
        """
        return self.ptrx().getTriVDepSReacActive(steps.triangle_global_id(idx), to_std_string(vsr))

    def setTriVDepSReacActive(self, steps.index_t idx, str vsr, bool act):
        """
        Activates/inactivates voltage-dependent surface reaction vsr in a
        triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param vsr: name of the voltage-dependent surface reaction.
        :type vsr: str
        :param act: Flag to activate / deactivate the reaction.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setTriVDepSReacActive(steps.triangle_global_id(idx), to_std_string(vsr), act)

    def setTriCapac(self, steps.index_t idx, double cm):
        """
        Set the specific capacitance of a triangle surface element.

        :param idx: Index of the triangle surface element
        :type idx: steps.index_t
        :param cm: Specific membrane capacitance (farad / m^2)
        :type cm: float


        :rtype: None
        """
        self.ptrx().setTriCapac(steps.triangle_global_id(idx), cm)

    def getVertV(self, steps.index_t vidx):
        """
        Returns the potential of vertex in Volts.

        :param vidx: Index of the vertex.
        :type vidx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getVertV(steps.vertex_id_t(vidx))

    def setVertV(self, steps.index_t vidx, double v):
        """
        Set the potential of vertex.

        :param vidx: Index of the vertex
        :type vidx: steps.index_t
        :param v:
        :type v: float


        :rtype: None
        """
        self.ptrx().setVertV(steps.vertex_id_t(vidx), v)

    def getVertVClamped(self, steps.index_t vidx):
        """
        Returns whether the potential of vertex is clamped over time
        (unless changed explicitly).

        :param vidx: Index of the vertex
        :type vidx: steps.index_t


        :rtype: bool
        """
        return self.ptrx().getVertVClamped(steps.vertex_id_t(vidx))

    def setVertVClamped(self, steps.index_t vidx, bool cl):
        """
        Sets voltage clamp in vertex.

        :param vidx: Index of the vertex
        :type vidx: steps.index_t
        :param cl: Flag to turn the clamping on or off.
        :type cl: bool


        :rtype: None
        """
        self.ptrx().setVertVClamped(steps.vertex_id_t(vidx), cl)

    def getVertIClamp(self, steps.index_t vidx):
        """
        Gets current injection to vertex.

        :param vidx: Index of the vertex
        :type vidx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getVertIClamp(steps.vertex_id_t(vidx))

    def setVertIClamp(self, steps.index_t vidx, double i):
        """
        Sets current injection to vertex.
        Will be assumed to be constant for one EField DT

        :param vidx: Index of the vertex
        :type vidx: steps.index_t
        :param i:
        :type i: float


        :rtype: None
        """
        self.ptrx().setVertIClamp(steps.vertex_id_t(vidx), i)

    def setMembPotential(self, str m, double v):
        """
        Set the electric potential of the membrane, including all nodes
        in the conduction volume.

        :param m: Name of the membrane
        :type m: str
        :param v: Potential (volts)
        :type v: float


        :rtype: None
        """
        self.ptrx().setMembPotential(to_std_string(m), v)

    def setMembCapac(self, str m, double cm):
        """
        Set the specific membrane capacitance of the membrane

        :param m: Name of the membrane
        :type m: str
        :param cm: Specific membrane capacitance (farad / m^2)
        :type cm: float


        :rtype: None
        """
        self.ptrx().setMembCapac(to_std_string(m), cm)

    def setMembVolRes(self, str m, double ro):
        """
        Set the bulk electrical resistivity of the section of the mesh
        representing the volume conductor

        :param m: Name of the membrane
        :type m: str
        :param ro: Electrical resistivity (ohm.m)
        :type ro: float


        :rtype: None
        """
        self.ptrx().setMembVolRes(to_std_string(m), ro)

    def setMembRes(self, str m, double ro, double vrev):
        """
        Sets the surface electrical resistivity ro (in ohm.m^2) of the membrane with string identifier memb. Reversal potential vrev is required in Volts.

        :param m: Name of the membrane
        :type m: str
        :param ro: membrane resistivity (ohm.m^2)
        :type ro: float
        :param vrev: Reversal potential (Volts)
        :type vrev: float


        :rtype: None
        """
        self.ptrx().setMembRes(to_std_string(m), ro, vrev)

    def setOutputSync(self, bool enable_sync, int output_rank):
        """
        Set if the outputs are synced across all ranks.
        If not, set the output rank.

        :param enable_sync:
        :type enable_sync: bool
        :param output_rank:
        :type output_rank: int


        :rtype: None
        """
        self.ptrx().setOutputSync(enable_sync, output_rank)

    def getOutputSyncStatus(self, ):
        """
        Get if the outputs are synced across all ranks.


        :rtype: bool
        """
        return self.ptrx().getOutputSyncStatus()

    def getOutputSyncRank(self, ):
        """
        Get the rank id for data output.


        :rtype: int
        """
        return self.ptrx().getOutputSyncRank()

    @staticmethod
    cdef _py_TetVesicleRDEF from_ptr(TetVesicleRDEF *ptr):
        cdef _py_TetVesicleRDEF obj = _py_TetVesicleRDEF.__new__(_py_TetVesicleRDEF )
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_TetVesicleRDEF from_ref(const TetVesicleRDEF &ref):
        _py_TetVesicleRDEF.from_ptr(<TetVesicleRDEF*>&ref)

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TetVesicleVesRaft(_py_TetAPI):
    """Bindings for MPI TetVesicleVesRaft"""
# ----------------------------------------------------------------------------------------------------------------------
    cdef TetVesicleVesRaft *ptrx(self):
        return <TetVesicleVesRaft*> self._ptr

    def __init__(self, _py_Model model, _py_Geom geom, _py_RNG rng, int calcMembPot=0):
        """
        Construction::

            sim = steps.mpi.solver.TetVesicle(model, geom, rng, calcMembPot=0)

        Create a spatial stochastic solver based on operator-splitting, which also supports vesicles, 'rafts' and related phenomena such as exocytosis and endocytosis.
        If voltage is to be simulated, argument calcMembPot specifies the solver. E.g. calcMembPot=steps.solver.EF_DV_PETSC will utilise the PETSc library. calcMembPot=0 means that voltage will not be simulated.

        Arguments:
        steps.model.Model model
        steps.geom.Geom geom
        steps.rng.RNG rng
        int calcMemPot (default=0)

        """
        # We constructed a map. Now call constructor
        if model == None:
            raise TypeError('The Model object is empty.')
        if geom == None:
            raise TypeError('The Geom object is empty.')
        if rng == None:
            raise TypeError('The RNG object is empty.')
        self._ptr = new TetVesicleVesRaft(model.ptr(), geom.ptr(), rng.ptr(), calcMembPot)

    def getSolverName(self, ):
        """
        Return the solver's name.


        :rtype: str
        """
        return from_std_string(self.ptrx().getSolverName())

    def getSolverDesc(self, ):
        """
        Return the solver's description.


        :rtype: str
        """
        return from_std_string(self.ptrx().getSolverDesc())

    def getSolverAuthors(self, ):
        """
        Return the solver's author.


        :rtype: str
        """
        return from_std_string(self.ptrx().getSolverAuthors())

    def getSolverEmail(self, ):
        """
        Return the solver's email.


        :rtype: str
        """
        return from_std_string(self.ptrx().getSolverEmail())

    def reset(self, ):
        """
        Reset the solver.


        :rtype: None
        """
        self.ptrx().reset()

    def run(self, double endtime):
        """
        Advance the simulation until endtime (given in seconds) is reached.
        The endtime must be larger or equal to the current simulation time.

        :param endtime: Time to end the solver.
        :type endtime: float


        :rtype: None
        """
        self.ptrx().run(endtime)

    def advance(self, double adv):
        """
        Advance the solver a given amount of time.

        :param adv: Time to advance the solver (in seconds)
        :type adv: float


        :rtype: None
        """
        self.ptrx().advance(adv)

    def step(self, ):
        """
        Advance the simulation for one 'step'. In stochastic solvers this is
        one 'realization' of the Gillespie SSA (one reaction 'event'). In
        numerical solvers (currently Wmrk4) this is one time-step, with the
        stepsize defined with the setDT method.


        :rtype: None
        """
        self.ptrx().step()

    def checkpoint(self, str file_name):
        """
        checkpoint simulator state to a file

        :param file_name: file name prefix (will be suffixed with the rank of each process)
        :type file_name: str


        :rtype: None
        """
        fullName = file_name + str(f'_{getRank()}')
        self.ptrx().checkpoint(to_std_string(fullName))

    def restore(self, str file_name):
        """
        restore simulator state from a file

        :param file_name: file name prefix (does not include the rank suffix)
        :type file_name: str


        :rtype: None
        """
        fullName = file_name + str(f'_{getRank()}')
        self.ptrx().restore(to_std_string(fullName))

    def getTime(self, ):
        """
        Returns the current simulation time in seconds.


        :rtype: float
        """
        return self.ptrx().getTime()

    def getA0(self, ):
        """
        Returns the total propensity of the current simulation state
        (the total propensity multiplied by an infinitesimally small
        time dt gives the probability that a reaction will occur in that dt).
        For Tetexact this includes the propensity from the extension of the SSA
        for diffusive flux between tetrahedral elements in the mesh.


        :rtype: float
        """
        return self.ptrx().getA0()

    def getNSteps(self, ):
        """
        Return the number of 'realizations' of the SSA, the number of reaction
        (and diffusion) events in stochastic solvers.


        :rtype: int
        """
        return self.ptrx().getNSteps()

    def setTime(self, double time):
        """
        Set the current simulation time.

        Syntax::

            setTime(time)

        Arguments:
        float time

        Return:
        None

        """
        self.ptrx().setTime(time)

    def setVesicleDT(self, double dt):
        """
        Set the default vesicle dt (Note: the actual vesicle dt
        used in the simulation can be lower than this number depending
        on simulation conditions).

        :param dt: Vesicle dt (in seconds)
        :type dt: float


        :rtype: None
        """
        self.ptrx().setVesicleDT(dt)

    def getVesicleDT(self, ):
        """
        Return the present vesicle dt (in seconds)


        :rtype: float
        """
        return self.ptrx().getVesicleDT()

    def createPath(self, str path):
        """
        Create a path

        :param path: Name of the path
        :type path: str


        :rtype: None
        """
        self.ptrx().createPath(to_std_string(path))

    def addPathPoint(self, str path, uint point_idx, std.vector[double] position):
        """
        Add a point to a path

        :param path: Name of the path
        :type path: str
        :param point_idx: An index for the point, positive integer
        :type point_idx: int
        :param position: Position of the point in cartesian coordinates
        :type position: List[float]


        :rtype: None
        """
        self.ptrx().addPathPoint(to_std_string(path), point_idx, position)

    def addPathBranch(self, str path, uint sourcepoint_idx, std.map[uint, double] destpoints_indxs):
        """
        Create branching in the path from a point

        :param path: Name of the path
        :type path: str
        :param sourcepoint_idx: An index for the source point, positive integer
        :type sourcepoint_idx: int
        :param destpoints_indxs:
        :type destpoints_indxs: Dict[int, float]


        :rtype: None
        """
        self.ptrx().addPathBranch(to_std_string(path), sourcepoint_idx, destpoints_indxs)

    def _getAllPaths(self):
        """
        Return a structure describing all the paths added so far
        Only rank 0 should call this method

        :rtype: Dict[str, Dict[int, Tuple[List[float], Dict[int, float]]]]
        """
        return {from_std_string(_p1.first): {_p2.first: [[v for v in _p2.second.first], {_p3.first: _p3.second for _p3 in _p2.second.second}] for _p2 in _p1.second} for _p1 in self.ptrx().getAllPaths()}

    def getBatchTetSpecCounts(self, std.vector[steps.index_t] tets, str s):
        """
        Get species counts of a list of tetrahedrons

        :param tets:
        :type tets: List[steps.index_t]
        :param s:
        :type s: str


        :rtype: List[float]
        """
        return [_val1 for _val1 in self.ptrx().getBatchTetSpecCounts(tets, to_std_string(s))]

    def getBatchTriSpecCounts(self, std.vector[steps.index_t] tris, str s):
        """
        Get species counts of a list of triangles

        :param tris:
        :type tris: List[steps.index_t]
        :param s:
        :type s: str


        :rtype: List[float]
        """
        return [_val1 for _val1 in self.ptrx().getBatchTriSpecCounts(tris, to_std_string(s))]

    def getBatchTetSpecCountsNP(self, index_t[:] indices, str s, double[:] counts):
        """
        Get species counts of a list of tetrahedrons

        :param indices:
        :type indices: numpy.array<index_t>
        :param s:
        :type s: str
        :param counts:
        :type counts: numpy.array<float>


        :rtype: None
        """
        self.ptrx().getBatchTetSpecCountsNP(&indices[0], indices.shape[0], to_std_string(s), &counts[0], counts.shape[0])

    def getBatchTriSpecCountsNP(self, index_t[:] indices, str s, double[:] counts):
        """
        Get species counts of a list of triangles

        :param indices:
        :type indices: numpy.array<index_t>
        :param s:
        :type s: str
        :param counts:
        :type counts: numpy.array<float>


        :rtype: None
        """
        self.ptrx().getBatchTriSpecCountsNP(&indices[0], indices.shape[0], to_std_string(s), &counts[0], counts.shape[0])

    def setROITetSpecClamped(self, std.vector[steps.index_t] triangles, str s, bool b):
        """

        :param triangles:
        :type triangles: List[steps.index_t]
        :param s:
        :type s: str
        :param b:
        :type b: bool


        :rtype: None
        """
        cdef std.vector[steps.tetrahedron_global_id] std_triangles
        std_triangles.reserve(len(triangles))
        for _val in triangles:
            std_triangles.push_back(steps.tetrahedron_global_id(_val))
        self.ptrx().setROITetSpecClamped(std_triangles, to_std_string(s), b)

    def setROITriSpecClamped(self, std.vector[steps.index_t] triangles, str s, bool b):
        """

        :param triangles:
        :type triangles: List[steps.index_t]
        :param s:
        :type s: str
        :param b:
        :type b: bool


        :rtype: None
        """
        cdef std.vector[steps.triangle_global_id] std_triangles
        std_triangles.reserve(len(triangles))
        for _val in triangles:
            std_triangles.push_back(steps.triangle_global_id(_val))
        self.ptrx().setROITriSpecClamped(std_triangles, to_std_string(s), b)

    def getROITetSpecCount(self, std.vector[steps.index_t] triangles, str s):
        """

        :param triangles:
        :type triangles: List[steps.index_t]
        :param s:
        :type s: str


        :rtype: float
        """
        cdef std.vector[steps.tetrahedron_global_id] std_triangles
        std_triangles.reserve(len(triangles))
        for _val in triangles:
            std_triangles.push_back(steps.tetrahedron_global_id(_val))
        return self.ptrx().getROITetSpecCount(std_triangles, to_std_string(s))

    def getROITriSpecCount(self, std.vector[steps.index_t] triangles, str s):
        """

        :param triangles:
        :type triangles: List[steps.index_t]
        :param s:
        :type s: str


        :rtype: float
        """
        cdef std.vector[steps.triangle_global_id] std_triangles
        std_triangles.reserve(len(triangles))
        for _val in triangles:
            std_triangles.push_back(steps.triangle_global_id(_val))
        return self.ptrx().getROITriSpecCount(std_triangles, to_std_string(s))

    def setROITetSpecCount(self, std.vector[steps.index_t] triangles, str s, double count):
        """

        :param triangles:
        :type triangles: List[steps.index_t]
        :param s:
        :type s: str
        :param count:
        :type count: float


        :rtype: None
        """
        cdef std.vector[steps.tetrahedron_global_id] std_triangles
        std_triangles.reserve(len(triangles))
        for _val in triangles:
            std_triangles.push_back(steps.tetrahedron_global_id(_val))
        self.ptrx().setROITetSpecCount(std_triangles, to_std_string(s), count)

    def setROITriSpecCount(self, std.vector[steps.index_t] triangles, str s, double count):
        """

        :param triangles:
        :type triangles: List[steps.index_t]
        :param s:
        :type s: str
        :param count:
        :type count: float


        :rtype: None
        """
        cdef std.vector[steps.triangle_global_id] std_triangles
        std_triangles.reserve(len(triangles))
        for _val in triangles:
            std_triangles.push_back(steps.triangle_global_id(_val))
        self.ptrx().setROITriSpecCount(std_triangles, to_std_string(s), count)

    def getROITetSpecCounts(self, str roi, str s):
        """
        Get species counts of a list of tetrahedrons

        :param roi:
        :type roi: str
        :param s:
        :type s: str


        :rtype: List[float]
        """
        return [_val1 for _val1 in self.ptrx().getROITetSpecCounts(to_std_string(roi), to_std_string(s))]

    def getROITriSpecCounts(self, str roi, str s):
        """
        Get species counts of a list of triangles

        :param roi:
        :type roi: str
        :param s:
        :type s: str


        :rtype: List[float]
        """
        return [_val1 for _val1 in self.ptrx().getROITriSpecCounts(to_std_string(roi), to_std_string(s))]

    def getROITetSpecCountsNP(self, str roi, str s, double[:] counts):
        """
        Get species counts of a list of tetrahedrons

        :param roi:
        :type roi: str
        :param s:
        :type s: str
        :param counts:
        :type counts: numpy.array<float>


        :rtype: None
        """
        self.ptrx().getROITetSpecCountsNP(to_std_string(roi), to_std_string(s), &counts[0], counts.shape[0])

    def getROITriSpecCountsNP(self, str roi, str s, double[:] counts):
        """
        Get species counts of a list of triangles

        :param roi:
        :type roi: str
        :param s:
        :type s: str
        :param counts:
        :type counts: numpy.array<float>


        :rtype: None
        """
        self.ptrx().getROITriSpecCountsNP(to_std_string(roi), to_std_string(s), &counts[0], counts.shape[0])

    def getROIVol(self, str roi):
        """
        Get the volume of a ROI.

        :param roi:
        :type roi: str


        :rtype: float
        """
        return self.ptrx().getROIVol(to_std_string(roi))

    def getROIArea(self, str roi):
        """
        Get the area of a ROI.

        :param roi:
        :type roi: str


        :rtype: float
        """
        return self.ptrx().getROIArea(to_std_string(roi))

    def getROISpecCount(self, str roi, str s):
        """
        Get the count of a species in a ROI.

        :param roi:
        :type roi: str
        :param s:
        :type s: str


        :rtype: float
        """
        return self.ptrx().getROISpecCount(to_std_string(roi), to_std_string(s))

    def setROISpecCount(self, str roi, str s, double count):
        """
        Set the count of a species in a ROI.

        :param roi:
        :type roi: str
        :param s:
        :type s: str
        :param count:
        :type count: float


        :rtype: None
        """
        self.ptrx().setROISpecCount(to_std_string(roi), to_std_string(s), count)

    def getROISpecAmount(self, str roi, str s):
        """
        Get the amount of a species in a ROI.

        :param roi:
        :type roi: str
        :param s:
        :type s: str


        :rtype: float
        """
        return self.ptrx().getROISpecAmount(to_std_string(roi), to_std_string(s))

    def setROISpecAmount(self, str roi, str s, double amount):
        """
        Set the amount of a species in a ROI.

        :param roi:
        :type roi: str
        :param s:
        :type s: str
        :param amount:
        :type amount: float


        :rtype: None
        """
        self.ptrx().setROISpecAmount(to_std_string(roi), to_std_string(s), amount)

    def getROISpecConc(self, str roi, str s):
        """
        Get the concentration of a species in a ROI.

        :param roi:
        :type roi: str
        :param s:
        :type s: str


        :rtype: float
        """
        return self.ptrx().getROISpecConc(to_std_string(roi), to_std_string(s))

    def setROISpecConc(self, str roi, str s, double conc):
        """
        Set the concentration of a species in a ROI.

        :param roi:
        :type roi: str
        :param s:
        :type s: str
        :param conc:
        :type conc: float


        :rtype: None
        """
        self.ptrx().setROISpecConc(to_std_string(roi), to_std_string(s), conc)

    def setROISpecClamped(self, str roi, str s, bool clamped):
        """
        Set a species in a ROI to be clamped or not. The count of species spec in the ROI is clamped if
        clamped is True, not clamped if clamped is False.

        Syntax::
            setROISpecClamped(roi, s, clamped)

        Arguments:
        string roi
        string s
        bool clamped

        Return:
        None

        """
        self.ptrx().setROISpecClamped(to_std_string(roi), to_std_string(s), clamped)

    def setROIReacK(self, str roi, str r, double kf):
        """
        Sets the macroscopic reaction constant of reaction with identifier
        string r in a ROI with identifier string roi to kf. The unit of the
        reaction constant depends on the order of the reaction.
        Note: The default value still comes from the steps.model description,
        so calling reset() will return the reaction constant to that value.

        :param roi:
        :type roi: str
        :param r:
        :type r: str
        :param kf:
        :type kf: float


        :rtype: None
        """
        self.ptrx().setROIReacK(to_std_string(roi), to_std_string(r), kf)

    def setROISReacK(self, str roi, str sr, double kf):
        """
        Sets the macroscopic reaction constant of surface reaction with
        identifier string sr in a ROI with identifier string roi to kf. The
        unit of the reaction constant depends on the order of the reaction.
        Note: The default value still comes from the steps.model description,
        so calling reset() will return the reaction constant to that value.

        :param roi:
        :type roi: str
        :param sr:
        :type sr: str
        :param kf:
        :type kf: float


        :rtype: None
        """
        self.ptrx().setROISReacK(to_std_string(roi), to_std_string(sr), kf)

    def setROIDiffD(self, str roi, str diff, double dcst):
        """
        Sets the macroscopic diffusion constant of diffusion with identifier string diff
        in a ROI with identifier string roi to dcst.

        Note: The default value still comes from the steps.model description, so
        calling reset() will return the diffusion constant to that value.

        Syntax::
            setROIDiffD(roi, diff, dcst)

        Arguments:
        string roi
        string diff
        float dcst

        Return:
            None

        """
        self.ptrx().setROIDiffD(to_std_string(roi), to_std_string(diff), dcst)

    def setROIReacActive(self, str roi, str r, bool a):
        """
        Set reaction r in a ROI to be active or not.

        :param roi:
        :type roi: str
        :param r:
        :type r: str
        :param a:
        :type a: bool


        :rtype: None
        """
        self.ptrx().setROIReacActive(to_std_string(roi), to_std_string(r), a)

    def setROISReacActive(self, str roi, str sr, bool a):
        """
        Set surface reaction sr in a ROI to be active or not.

        :param roi:
        :type roi: str
        :param sr:
        :type sr: str
        :param a:
        :type a: bool


        :rtype: None
        """
        self.ptrx().setROISReacActive(to_std_string(roi), to_std_string(sr), a)

    def setROIDiffActive(self, str roi, str d, bool act):
        """
        Set diffusion d in a ROI to be active or not.

        :param roi:
        :type roi: str
        :param d:
        :type d: str
        :param act:
        :type act: bool


        :rtype: None
        """
        self.ptrx().setROIDiffActive(to_std_string(roi), to_std_string(d), act)

    def setROIVDepSReacActive(self, str roi, str vsr, bool a):
        """
        Set voltage dependent surface reaction vsr in a ROI to be active or
        not.

        :param roi:
        :type roi: str
        :param vsr:
        :type vsr: str
        :param a:
        :type a: bool


        :rtype: None
        """
        self.ptrx().setROIVDepSReacActive(to_std_string(roi), to_std_string(vsr), a)

    def getROIReacExtent(self, str roi, str r):
        """
        Return the extent of reaction with identifier string reac in ROI with
        identifier string roi, that is the number of times the reaction has
        occurred up to the current simulation time.

        :param roi:
        :type roi: str
        :param r:
        :type r: str


        :rtype: unsigned long long
        """
        return self.ptrx().getROIReacExtent(to_std_string(roi), to_std_string(r))

    def resetROIReacExtent(self, str roi, str r):
        """
        Reset the extent of reaction with identifier string reac in ROI with
        identifier string roi, that is the number of times the reaction has
        occurred up to the current simulation time, to 0.

        :param roi:
        :type roi: str
        :param r:
        :type r: str


        :rtype: None
        """
        self.ptrx().resetROIReacExtent(to_std_string(roi), to_std_string(r))

    def getROISReacExtent(self, str roi, str sr):
        """
        Return the extent of surface reaction with identifier string sr in ROI
        with identifier string roi, that is the number of times the reaction
        has occurred up to the current simulation time.

        :param roi:
        :type roi: str
        :param sr:
        :type sr: str


        :rtype: unsigned long long
        """
        return self.ptrx().getROISReacExtent(to_std_string(roi), to_std_string(sr))

    def resetROISReacExtent(self, str roi, str sr):
        """
        Reset the extent of surface reaction with identifier string reac in ROI
        with identifier string roi, that is the number of times the reaction
        has occurred up to the current simulation time, to 0.

        :param roi:
        :type roi: str
        :param sr:
        :type sr: str


        :rtype: None
        """
        self.ptrx().resetROISReacExtent(to_std_string(roi), to_std_string(sr))

    def getROIDiffExtent(self, str roi, str d):
        """
        Return the extent of diffusion with identifier string diff in ROI with
        identifier string roi, that is the number of times the diffusion has
        occurred up to the current simulation time.

        :param roi:
        :type roi: str
        :param d:
        :type d: str


        :rtype: unsigned long long
        """
        return self.ptrx().getROIDiffExtent(to_std_string(roi), to_std_string(d))

    def resetROIDiffExtent(self, str roi, str s):
        """
        Reset the extent of diffusion with identifier string diff in ROI with
        identifier string roi, that is the number of times the diffusion has
        occurred up to the current simulation time, to 0.

        :param roi:
        :type roi: str
        :param s:
        :type s: str


        :rtype: None
        """
        self.ptrx().resetROIDiffExtent(to_std_string(roi), to_std_string(s))


    def setEfieldDT(self, double efdt):
        """
        Set the stepsize for membrane potential solver (default 1e-6s).
        This is the time for each voltage calculation step. The SSA will
        run until passing this stepsize, so in fact each membrane potential
        time step will vary slightly around the dt so as to be aligned with the
        SSA.

        :param efdt:
        :type efdt: float


        :rtype: None
        """
        self.ptrx().setEfieldDT(efdt)

    def setTemp(self, double temp):
        """
        Set the simulation temperature. Currently, this will only
        influence the GHK flux rate, so will only influence simulations
        including membrane potential calculation.

        :param temp: EField temperature (in Kelvin)
        :type temp: float


        :rtype: None
        """
        self.ptrx().setTemp(temp)

    def getEfieldDT(self, ):
        """
        Return the stepsize for membrane potential solver (in seconds)


        :rtype: float
        """
        return self.ptrx().getEfieldDT()

    def getTemp(self, ):
        """
        Return the simulation temperature (in Kelvin)


        :rtype: float
        """
        return self.ptrx().getTemp()

    def getCompVol(self, str c):
        """
        Returns the volume of compartment c (in m^3).

        :param c: Name of the compartment.
        :type c: str


        :rtype: float
        """
        return self.ptrx().getCompVol(to_std_string(c))

    def getCompSpecCount(self, str c, str s):
        """
        Returns the number of molecules of species s in compartment c.
        NOTE: in a mesh-based simulation, the total count is computed as
        the sum of the counts in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getCompSpecCount(to_std_string(c), to_std_string(s))

    def setCompSpecCount(self, str c, str s, double n):
        """
        Sets the number of molecules of species s in compartment c.
        NOTE: in a mesh-based simulation, the total amount is equally divided
        over all tetrahedrons in the compartment (i.e. a uniform distribution).

        :param c: Name of the compartment.
        :type c: str
        :param s: Name of the species.
        :type s: str
        :param n: Number of molecules of the species.
        :type n: float


        :rtype: None
        """
        self.ptrx().setCompSpecCount(to_std_string(c), to_std_string(s), n)

    def getCompSpecAmount(self, str c, str s):
        """
        Returns the amount (in mols) of species s in compartment c.
        NOTE: in a mesh-based simulation, the total amount is computed as
        the sum of the amounts in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getCompSpecAmount(to_std_string(c), to_std_string(s))

    def setCompSpecAmount(self, str c, str s, double a):
        """
        Set the amount (in mols) of species s in compartment c.
        NOTE: in a mesh-based simulation, the total amount is equally divided
        over all tetrahedrons in the compartment (i.e. a uniform distribution).

        :param c: Name of the compartment.
        :type c: str
        :param s: Name of the species.
        :type s: str
        :param a: Amount of the species.
        :type a: float


        :rtype: None
        """
        self.ptrx().setCompSpecAmount(to_std_string(c), to_std_string(s), a)

    def getCompSpecConc(self, str c, str s):
        """
        Returns the concentration (in molar units) of species s in compartment
        c.
        NOTE: in a mesh-based simulation, the overall concentration in a
        compartment is computed by taking the volume-weighted sum of the
        concentration in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getCompSpecConc(to_std_string(c), to_std_string(s))

    def setCompSpecConc(self, str c, str s, double conc):
        """
        Sets the concentration (in molar units) of species s in compartment c.
        NOTE: in a mesh-based simulation, this method changes the
        concentration to the same value in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param s: Name of the species.
        :type s: str
        :param conc: Concentration of the species.
        :type conc: float


        :rtype: None
        """
        self.ptrx().setCompSpecConc(to_std_string(c), to_std_string(s), conc)

    def getCompSpecClamped(self, str c, str s):
        """
        Returns whether the concentration of species s in compartment c
        remains constant over time (unless changed explicitly).
        NOTE: in a mesh-based simulation, this method will only return true
        only if the species has been clamped in all tetrahedrons of the
        compartment. \param c Name of the compartment. \param s Name of the
        species.

        :param c:
        :type c: str
        :param s:
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getCompSpecClamped(to_std_string(c), to_std_string(s))

    def setCompSpecClamped(self, str c, str s, bool b):
        """
        Turns clamping of species s in compartment c on or off.
        NOTE: in a mesh based simulation, this method turns clamping on/off
        in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param s: Name of the species.
        :type s: str
        :param b: Flag to trun clamping of species on / off.
        :type b: bool


        :rtype: None
        """
        self.ptrx().setCompSpecClamped(to_std_string(c), to_std_string(s), b)

    def getCompReacK(self, str c, str r):
        """
        Returns the macroscopic reaction constant of reaction r in
        compartment c.
        Note: in a mesh-based simulation, the value is computed as the
        volume-weighted sum of the reaction constants in all tetrahedrons of
        the compartment.

        :param c:
        :type c: str
        :param r:
        :type r: str


        :rtype: float
        """
        return self.ptrx().getCompReacK(to_std_string(c), to_std_string(r))

    def setCompReacK(self, str c, str r, double kf):
        """
        Sets the macroscopic reaction constant of reaction r in compartment c
        (units vary according to the order of the reaction).
        NOTE: in a mesh-based simulation, this method changes the reaction
        constant equally in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of te reaction.
        :type r: str
        :param kf: Reaction constant.
        :type kf: float


        :rtype: None
        """
        self.ptrx().setCompReacK(to_std_string(c), to_std_string(r), kf)

    def getCompReacActive(self, str c, str r):
        """
        Returns whether reaction r in compartment c is active or not
        NOTE: in a mesh-based simulation, this method returns false only when
        the reaction has been inactivated in all tetrahedrons.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: bool
        """
        return self.ptrx().getCompReacActive(to_std_string(c), to_std_string(r))

    def setCompReacActive(self, str c, str r, bool a):
        """
        Activate or inactivate reaction r in compartment c.
        NOTE: in a mesh-based simulation, activation/inactivation of a reaction
        turns it on/off in all tetrahedrons.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of the reaction.
        :type r: str
        :param a: Flag to activate or deactivate the reaction.
        :type a: bool


        :rtype: None
        """
        self.ptrx().setCompReacActive(to_std_string(c), to_std_string(r), a)

    def getCompDiffD(self, str c, str d):
        """
        Returns the diffusion constant of diffusion rule d in compartment c.

        :param c: Name of the compartment.
        :type c: str
        :param d: Name of the diffusion.
        :type d: str


        :rtype: float
        """
        return self.ptrx().getCompDiffD(to_std_string(c), to_std_string(d))

    def setCompDiffD(self, str c, str d, double dcst):
        """
        Set the diffusion constant of diffusion rule d in compartment c.

        :param c: Name of the compartment.
        :type c: str
        :param d: Name of the diffusion.
        :type d: str
        :param dcst: Rate constant of the diffusion.
        :type dcst: float


        :rtype: None
        """
        self.ptrx().setCompDiffD(to_std_string(c), to_std_string(d), dcst)

    def getCompDiffActive(self, str c, str d):
        """
        Returns whether diffusion rule d in compartment c is active or not.

        :param c: Name of the compartment.
        :type c: str
        :param d: Name of the diffusion.
        :type d: str


        :rtype: bool
        """
        return self.ptrx().getCompDiffActive(to_std_string(c), to_std_string(d))

    def setCompDiffActive(self, str c, str d, bool act):
        """
        Activate or deactivate diffusion rule d in compartment c.

        :param c: Name of the compartment.
        :type c: str
        :param d: Name of the diffusion.
        :type d: str
        :param act: Flag to activate or deactivate the diffusion.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setCompDiffActive(to_std_string(c), to_std_string(d), act)

    def getCompReacC(self, str c, str r):
        """
        Returns c_mu, the mesoscopic reaction constant of reaction r in
        compartment c.
        NOTE: in a mesh-based simulation, the mesoscopic reaction constant is
        computed as the sum of the mesoscopic constants in all tetrahedrons of
        the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getCompReacC(to_std_string(c), to_std_string(r))

    def getCompReacH(self, str c, str r):
        """
        Returns h_mu, the distinct number of ways in which reaction r can
        occur in compartment c, by computing the product of its reactants.
        NOTE: in a mesh-based simulation, it returns the sum of the h_mu's
        over all tetrahedrons of the compartment. This can become a very large
        value.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getCompReacH(to_std_string(c), to_std_string(r))

    def getCompReacA(self, str c, str r):
        """
        Returns the propensity, a_mu, of reaction r in compartment c.
        The propensity value gives the probability per unit time that this
        reaction will occur in the current state.
        NOTE: in a mesh-based simulation, a_mu is computed as the sum of the
        a_mu in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getCompReacA(to_std_string(c), to_std_string(r))

    def getCompReacExtent(self, str c, str r):
        """
        Returns the extent of reaction r in compartment c.
        NOTE: in a mesh-based simulation, returns the sum of the reaction
        extents in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: unsigned long long
        """
        return self.ptrx().getCompReacExtent(to_std_string(c), to_std_string(r))

    def resetCompReacExtent(self, str c, str r):
        """
        Resets the extent of reaction r in compartment c to zero.
        NOTE: in a mesh-based simulation, resets the extents of the reaction
        in all tetrahedrons of the compartment.

        :param c: Name of the compartment.
        :type c: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: None
        """
        self.ptrx().resetCompReacExtent(to_std_string(c), to_std_string(r))

    def getCompVesicleCount(self, str c, str v):
        """
        Returns the number of vesicles v in compartment c

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str


        :rtype: int
        """
        return self.ptrx().getCompVesicleCount(to_std_string(c), to_std_string(v))

    def setCompVesicleCount(self, str c, str v, uint n):
        """
        Sets the number of vesicles v in compartment c

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str
        :param n: Number of vesicles
        :type n: int


        :rtype: None
        """
        self.ptrx().setCompVesicleCount(to_std_string(c), to_std_string(v), n)

    def addCompVesicle(self, str c, str v):
        """
        Adds one vesicle v to compartment c

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str


        :rtype: steps.index_t
        """
        return self.ptrx().addCompVesicle(to_std_string(c), to_std_string(v)).get()

    def deleteSingleVesicle(self, str v, steps.index_t ves_unique_index):
        """
        Deletes individual vesicle of type v with unique index ves_unique_index

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t


        :rtype: None
        """
        self.ptrx().deleteSingleVesicle(to_std_string(v), steps.vesicle_individual_id(ves_unique_index))

    def getSingleVesicleSurfaceLinkSpecCount(self, str v, steps.index_t ves_unique_index, str ls):
        """
        Returns the count of 'link' species ls on vesicle of type v with unique index
        vesicle_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param ls: Name of the link species
        :type ls: str


        :rtype: int
        """
        return self.ptrx().getSingleVesicleSurfaceLinkSpecCount(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(ls))

    def getSingleVesicleSurfaceLinkSpecIndices(self, str v, steps.index_t ves_unique_index, str ls):
        """
        Returns the indices of 'link' species ls on vesicle of type v with unique index
        vesicle_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param ls: Name of the link species
        :type ls: str


        :rtype: List[steps.index_t]
        """
        return [_val1.get() for _val1 in self.ptrx().getSingleVesicleSurfaceLinkSpecIndices(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(ls))]

    def getSingleVesicleSurfaceSpecIndices(self, str v, steps.index_t ves_unique_index, str s):
        """
        Returns the indices of point species s on vesicle of type v with unique index
        vesicle_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the point species
        :type s: str


        :rtype: List[steps.index_t]
        """
        return [_val1.get() for _val1 in self.ptrx().getSingleVesicleSurfaceSpecIndices(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s))]

    def getAllVesicleIndices(self):
        """
        Returns all the unique indices of all vesicles in the simulation


        :rtype: List[steps.index_t]
        """
        return [_val1.get() for _val1 in self.ptrx().getAllVesicleIndices()]

    def getAllVesicleIndicesOnPath(self):
        """
        Returns all the unique indices of all vesicles in the simulation currently on a path


        :rtype: List[steps.index_t]
        """
        return [_val1.get() for _val1 in self.ptrx().getAllVesicleIndicesOnPath()]

    def getCompVesicleIndices(self, str c, str v):
        """
        Returns all the unique indices of vesicles v currently in compartment c

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str


        :rtype: List[steps.index_t]
        """
        return [_val1.get() for _val1 in self.ptrx().getCompVesicleIndices(to_std_string(c), to_std_string(v))]

    def getSingleVesicleCompartment(self, str v, steps.index_t ves_unique_index):
        """
        Get the compartment of vesicle of type v and unique index ves_unique_index

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t


        :rtype: str
        """
        return from_std_string(self.ptrx().getSingleVesicleCompartment(to_std_string(v), steps.vesicle_individual_id(ves_unique_index)))

    def getSingleVesiclePos(self, str v, steps.index_t ves_unique_index):
        """
        Returns the position of vesicle of type v with unique index
        ves_unique_index  (in cartesian coordinates)

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t


        :rtype: List[float]
        """
        return [_val1 for _val1 in self.ptrx().getSingleVesiclePos(to_std_string(v), steps.vesicle_individual_id(ves_unique_index))]

    def setSingleVesiclePos(self, str v, steps.index_t ves_unique_index, std.vector[double] pos, bool force=False):
        """
        Set the position of vesicle of type v with unique index
        ves_unique_index to pos (cartesian coordinates)

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param pos: Position of the vesicle in cartesian corrdinates
        :type pos: List[float]
        :param force:
        :type force: bool

        :rtype: None
        """
        self.ptrx().setSingleVesiclePos(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), pos, force)

    # TODO Remove this method eventually and use setSingleVesiclePos instead
    def setCompSingleVesiclePos(self, str c, str v, steps.index_t ves_unique_index, std.vector[double] pos, bool force=False):
        """
        Set the position of vesicle of type v with unique index
        ves_unique_index in compartment c to pos (cartesian coordinates)

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param pos: Position of the vesicle in cartesian corrdinates
        :type pos: List[float]
        :param force:
        :type force: bool


        :rtype: None
        """
        self.ptrx().setSingleVesiclePos(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), pos, force)

    def getCompVesicleSurfaceSpecCountDict(self, str c, str v, str s):
        """
        Returns the surface count of species s on vesicles v in compartment c.
        Return is a map, vesicle_unique_index : count of s

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str
        :param s: Name of the species
        :type s: str


        :rtype: Dict[steps.index_t, int]
        """
        return { _pair1.first.get(): _pair1.second for _pair1 in self.ptrx().getCompVesicleSurfaceSpecCountDict(to_std_string(c), to_std_string(v), to_std_string(s)) }

    def getCompVesicleSurfaceSpecCount(self, str c, str v, str s):
        """
        Returns the summed surface count of species s on vesicles v in compartment c.

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str
        :param s: Name of the species
        :type s: str


        :rtype: int
        """
        return self.ptrx().getCompVesicleSurfaceSpecCount(to_std_string(c), to_std_string(v), to_std_string(s))

    def getCompVesicleInnerSpecCount(self, str c, str v, str s):
        """
        Returns the summed inner count of species s on vesicles v in compartment c.

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str
        :param s: Name of the species
        :type s: str


        :rtype: int
        """
        return self.ptrx().getCompVesicleInnerSpecCount(to_std_string(c), to_std_string(v), to_std_string(s))

    def getSingleVesicleSurfaceSpecCount(self, str v, steps.index_t ves_unique_index, str s):
        """
        Get the surface count of species s on vesicle of type v and unique
        index ves_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str


        :rtype: int
        """
        return self.ptrx().getSingleVesicleSurfaceSpecCount(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s))

    def getSingleVesicleInnerSpecCount(self, str v, steps.index_t ves_unique_index, str s):
        """
        Get the inner count of species s on vesicle of type v and unique
        index ves_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str


        :rtype: int
        """
        return self.ptrx().getSingleVesicleInnerSpecCount(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s))

    def setSingleVesicleSurfaceSpecCount(self, str v, steps.index_t ves_unique_index, str s, uint count):
        """
        Set the surface count of species s on vesicle of type v and unique
        index ves_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str
        :param count: The number of the surface species
        :type count: int


        :rtype: None
        """
        self.ptrx().setSingleVesicleSurfaceSpecCount(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s), count)

    def getSingleVesicleSurfaceSpecPos(self, str v, steps.index_t ves_unique_index, str s):
        """
        Get the cartesian coordinates of species s on vesicle of type v and
        unique index ves_unique_index. Position is absolute,

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str


        :rtype: List[List[float]]
        """
        return [[_val2 for _val2 in _val1] for _val1 in self.ptrx().getSingleVesicleSurfaceSpecPos(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s))]

    def setSingleVesicleInnerSpecCount(self, str v, steps.index_t ves_unique_index, str s, uint count):
        """
        Set the count of the inner species, that is inside the vesicle volume,
        on vesicle of type v and unique index ves_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str
        :param count: The number of the inner species
        :type count: int


        :rtype: None
        """
        self.ptrx().setSingleVesicleInnerSpecCount(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s), count)

    def getSingleVesicleSurfaceSpecPosSpherical(self, str v, steps.index_t ves_unique_index, str s):
        """
        Get the spherical coordinates of species s on vesicle of type v and
        unique index ves_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str

        :rtype: List[List[float]]
        """
        return [list(pos) for pos in self.ptrx().getSingleVesicleSurfaceSpecPosSpherical(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s))]

    def setSingleVesicleSurfaceSpecPosSpherical(self, str v, steps.index_t ves_unique_index, str s, list pos_spherical):
        """
        Set the spherical coordinates of species s on vesicle of type v and
        unique index ves_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str
        :param pos_spherical: Position of the molecule in spherical coordinates (relative to the vesicle)
        :type pos_spherical: List[float]


        :rtype: None
        """
        cdef std.vector[double] pos
        cdef std.vector[std.vector[double]] allPos
        if len(pos_spherical) == 2 and all(isinstance(val, (float, int)) for val in pos_spherical):
            for val in pos_spherical:
                pos.push_back(val)
            for i in range(self.getSingleVesicleSurfaceSpecCount(v, ves_unique_index, s)):
                allPos.push_back(pos)
        else:
            for ppos in pos_spherical:
                pos.clear()
                for val in ppos:
                    pos.push_back(val)
                allPos.push_back(pos)
        self.ptrx().setSingleVesicleSurfaceSpecPosSpherical(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(s), allPos)

    def getSingleSpecPosSpherical(self, str s, steps.index_t spec_unique_index):
        """
        Get the spherical coordinates of point species s with unique index spec_unique_index.

        :param s: Name of the species
        :type s: str
        :param spec_unique_index: Unique index of the individual point species
        :type spec_unique_index: steps.index_t

        :rtype: List[float]
        """
        return self.ptrx().getSingleSpecPosSpherical(to_std_string(s), steps.pointspec_individual_id(spec_unique_index))

    def getCompVesicleSurfaceLinkSpecCountDict(self, str c, str v, str ls):
        """
        Returns the count of 'link' species ls on vesicles v in compartment c.
        Return is a map, vesicle_unique_index : count of ls

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str
        :param ls: Name of the link species
        :type ls: str


        :rtype: Dict[steps.index_t, int]
        """
        return { _pair1.first.get(): _pair1.second for _pair1 in self.ptrx().getCompVesicleSurfaceLinkSpecCountDict(to_std_string(c), to_std_string(v), to_std_string(ls)) }

    def getCompVesicleSurfaceLinkSpecCount(self, str c, str v, str ls):
        """
        Returns the summed count of 'link' species ls on vesicles v in compartment c.

        :param c: Name of the compartment.
        :type c: str
        :param v: Name of the vesicle.
        :type v: str
        :param ls: Name of the link species
        :type ls: str


        :rtype: int
        """
        return self.ptrx().getCompVesicleSurfaceLinkSpecCount(to_std_string(c), to_std_string(v), to_std_string(ls))

    def getSingleVesicleSurfaceLinkSpecPos(self, str v, steps.index_t ves_unique_index, str ls):
        """
        Returns the positions of 'link' species ls on vesicle of type v and
        unique index ves_unique_index.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t
        :param ls: Name of the link species
        :type ls: str


        :rtype: List[List[float]]
        """
        return [[_val2 for _val2 in _val1] for _val1 in self.ptrx().getSingleVesicleSurfaceLinkSpecPos(to_std_string(v), steps.vesicle_individual_id(ves_unique_index), to_std_string(ls))]

    def getSingleLinkSpecPos(self, steps.index_t ls_unique_index):
        """
        Returns the position of the link species with unique index ls_unique_id.

        :param ls_unique_index: Unique index of the individual link species
        :type ls_unique_index: steps.index_t

        :rtype: List[float]
        """
        return [_val1 for _val1 in self.ptrx().getSingleLinkSpecPos(steps.linkspec_individual_id(ls_unique_index))]

    def getSingleLinkSpecLinkedTo(self, steps.index_t ls_unique_index):
        """
        Returns the unique index of the link species linked to the link species
        with unique index ls_unique_id.

        :param ls_unique_index: Unique index of the individual link species
        :type ls_unique_index: steps.index_t

        :rtype: steps.index_t
        """
        return self.ptrx().getSingleLinkSpecLinkedTo(steps.linkspec_individual_id(ls_unique_index)).get()

    def getSingleLinkSpecVes(self, steps.index_t ls_unique_index):
        """
        Returns the unique index of the vesicle that contains the link species
        with unique index ls_unique_id or UNDEFINED_VESICLE if the link species no longer exists.

        :param ls_unique_index: Unique index of the individual link species
        :type ls_unique_index: steps.index_t

        :rtype: steps.index_t
        """
        return self.ptrx().getSingleLinkSpecVes(steps.linkspec_individual_id(ls_unique_index)).get()

    def getSingleVesicleImmobility(self, str v, steps.index_t ves_unique_index):
        """
        Get the 'immobility' of vesicle of type v and unique index
        ves_unique_index. All non-zero numbers mean vesicle is
        immobile.

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t


        :rtype: int
        """
        return self.ptrx().getSingleVesicleImmobility(to_std_string(v), steps.vesicle_individual_id(ves_unique_index))

    def getSingleVesicleOverlapTets(self, str v, steps.index_t ves_unique_index):
        """
        Get the indexes of tetrahedrons that overlap vesicle of type v and
        unique index ves_unique_index

        :param v: Name of the vesicle.
        :type v: str
        :param ves_unique_index: Unique index of the individual vesicle
        :type ves_unique_index: steps.index_t


        :rtype: List[steps.index_t]
        """
        return [_val1.get() for _val1 in self.ptrx().getSingleVesicleOverlapTets(to_std_string(v), steps.vesicle_individual_id(ves_unique_index))]

    def setTetVesicleDcst(self, steps.index_t idx, str v, double dcst):
        """
        Set the diffusion rate per tetrahedron of vesicles of type v. Vesicles
        will use this diffusion rate when vesicle centre is in this tet.

        :param idx: Tetrahrdron index
        :type idx: steps.index_t
        :param v: Name of the vesicle.
        :type v: str
        :param dcst: Diffusion coefficient
        :type dcst: float


        :rtype: None
        """
        self.ptrx().setTetVesicleDcst(steps.tetrahedron_global_id(idx), to_std_string(v), dcst)

    def setVesicleSurfaceLinkSpecSDiffD(self, str v, str ls, double dcst):
        """
        Set the diffusion rate of 'link' species ls on vesicles of type v.

        :param v: Name of the vesicle.
        :type v: str
        :param ls: Name of the link species
        :type ls: str
        :param dcst: Diffusion coefficient
        :type dcst: float


        :rtype: None
        """
        self.ptrx().setVesicleSurfaceLinkSpecSDiffD(to_std_string(v), to_std_string(ls), dcst)

    def setVesSReacK(self, str vsr, double kf):
        """

        :param vsr: Name of the vesicle surface reaction
        :type vsr: str
        :param kf: Rate
        :type kf: float


        :rtype: None
        """
        self.ptrx().setVesSReacK(to_std_string(vsr), kf)

    def getVesSReacExtent(self, str vsr):
        """
        Get the reaction extent of vesicle surface reaction. Not
        compartment-specific.

        :param vsr: Name of the vesicle surface reaction
        :type vsr: str


        :rtype: int
        """
        return self.ptrx().getVesSReacExtent(to_std_string(vsr))

    def setExocytosisK(self, str exo, double kf):
        """
        Set rate of exocytosis events. Not compartment-specific.

        :param exo: Name of the vesicle exocytosis
        :type exo: str
        :param kf: Rate
        :type kf: float


        :rtype: None
        """
        self.ptrx().setExocytosisK(to_std_string(exo), kf)

    def getExocytosisExtent(self, str exo):
        """
        Get the extent of vesicle exocytosis. Not compartment-specific.

        :param exo: Name of the vesicle exocytosis
        :type exo: str


        :rtype: int
        """
        return self.ptrx().getExocytosisExtent(to_std_string(exo))

    def getExocytosisEvents(self, str exo):
        """
        Get the vesicle exocytosis events that happened since last call. Not compartment-specific.

        Returns a list of tuple with each tuple containing:
        (time of event, individual index of vesicle, triangle index where the exocytosis happened, individual index of raft or None)

        :param exo: Name of the vesicle exocytosis
        :type exo: str

        :rtype: List[Tuple[float, steps.index_t, steps.index_t, steps.index_t]]
        """
        return [(event.time, event.vidx.get(), event.tidx.get(), event.ridx.get() if event.ridx.valid() else None) for event in self.ptrx().getExocytosisEvents(to_std_string(exo))]

    def getRaftEndocytosisExtent(self, str rendo):
        """
        Get the extent of raft endocytosis. Not patch-specific.

        :param rendo: Name of the raft endocytosis
        :type rendo: str


        :rtype: int
        """
        return self.ptrx().getRaftEndocytosisExtent(to_std_string(rendo))

    def getRaftEndocytosisEvents(self, str rendo):
        """
        Get the raft endocytosis events that happened since last call. Not patch-specific.

        Returns a list of tuple with each tuple containing:
        (time of event, individual index of raft, triangle index where the endocytosis happened, individual index of vesicle)

        :param rendo: Name of the raft endocytosis
        :type rendo: str

        :rtype: List[Tuple[float, steps.index_t, steps.index_t, steps.index_t]]
        """
        return [(event.time, event.ridx.get(), event.tidx.get(), event.vidx.get()) for event in self.ptrx().getRaftEndocytosisEvents(to_std_string(rendo))]

    def setRaftEndocytosisK(self, str rendo, double kcst):
        """
        Set the rate of raft endocytosis. Not patch-specific.

        :param rendo: Name of the raft endocytosis
        :type rendo: str
        :param kcst: Rate of the raft endocytosis
        :type kcst: double

        """
        self.ptrx().setRaftEndocytosisK(to_std_string(rendo), kcst)

    def addVesicleDiffusionGroup(self, str v, list comps):
        """
        Add a 'diffusion group' for vesicles of type v. Vesicles will diffuse
        freely amongst this group (if they border each other).

        :param v: Name of the vesicle.
        :type v: str
        :param comps: List of compartment names.
        :type comps: List[str]


        :rtype: None
        """
        cdef std.vector[std.string] std_comps
        std_comps.reserve(len(comps))
        for _elem in comps:
            std_comps.push_back(to_std_string(_elem))
        self.ptrx().addVesicleDiffusionGroup(to_std_string(v), std_comps)

    def addPathVesicle(self, str path, str ves, double speed, dict spec_deps={}, stoch_stepsize=1e-9):
        """
        Add a vesicle to this path. This means a vesicle of this type can
        interact with this path upon overlapping it

        :param path: Name of the path
        :type path: str
        :param ves: Name of the vesicle
        :type ves: str
        :param speed: Speed of the vesicle on this path in m/s
        :type speed: float
        :param spec_deps: Optional species dependencies, in map of species names to number of the species required
        :type spec_deps: Dict[str, int]
        :param stoch_stepsize: Stochastic step length. This may be a single float value where a single-exponential will be applied. If a list of length 2, a double-exponential will be applied by the proportionality specified in the 2nd element.
        :type stoch_stepsize: Float or list[float]

        :rtype: None
        """
        cdef std.map[std.string, uint] std_spec_deps
        for _elem1, _elem2 in spec_deps.items():
            std_spec_deps[to_std_string(_elem1)] = _elem2

        cdef std.vector[double] std_stoch_stepsize
        if not hasattr(stoch_stepsize, '__iter__'):
            stoch_stepsize = [stoch_stepsize]
        std_stoch_stepsize.reserve(len(stoch_stepsize))
        for _elem in stoch_stepsize:
            std_stoch_stepsize.push_back(_elem)
        self.ptrx().addPathVesicle(to_std_string(path), to_std_string(ves), speed, std_spec_deps, std_stoch_stepsize)

    def getTetVol(self, steps.index_t idx):
        """
        Returns the volume of a tetrahedron (in m^3).

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTetVol(steps.tetrahedron_global_id(idx))

    def getTetReducedVol(self, steps.index_t idx):
        """
        Returns the reduced volume of a tetrahedron (in m^3).

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTetReducedVol(steps.tetrahedron_global_id(idx))

    def setTetVol(self, steps.index_t idx, double vol):
        """
        Set the volume of a tetrahedron (in m^3).

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param vol: Volume of the tetrahedron.
        :type vol: float


        :rtype: None
        """
        self.ptrx().setTetVol(steps.tetrahedron_global_id(idx), vol)

    def getTetSpecDefined(self, steps.index_t idx, str s):
        """
        Returns whether species s is defined in a tetrahedron

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getTetSpecDefined(steps.tetrahedron_global_id(idx), to_std_string(s))

    def getTetSpecCount(self, steps.index_t idx, str s):
        """
        Returns the number of molecules of species s in a tetrahedron

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getTetSpecCount(steps.tetrahedron_global_id(idx), to_std_string(s))

    def setTetSpecCount(self, steps.index_t idx, str s, double n):
        """
        Sets the number of molecules of species s in a tetrahedron

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str
        :param n: Number of molecules of the species.
        :type n: float


        :rtype: None
        """
        self.ptrx().setTetSpecCount(steps.tetrahedron_global_id(idx), to_std_string(s), n)

    def getTetSpecAmount(self, steps.index_t idx, str s):
        """
        Returns the amount (in mols) of species s in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getTetSpecAmount(steps.tetrahedron_global_id(idx), to_std_string(s))

    def setTetSpecAmount(self, steps.index_t idx, str s, double m):
        """
        Sets the amount (in mols) of species s in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str
        :param m: Amount of the species.
        :type m: float


        :rtype: None
        """
        self.ptrx().setTetSpecAmount(steps.tetrahedron_global_id(idx), to_std_string(s), m)

    def getTetSpecConc(self, steps.index_t idx, str s):
        """
        Returns the concentration (in molar units) of species s in a
        tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getTetSpecConc(steps.tetrahedron_global_id(idx), to_std_string(s))

    def setTetSpecConc(self, steps.index_t idx, str s, double c):
        """
        Sets the concentration (in molar units) of species s in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str
        :param c: Concentration of the species.
        :type c: float


        :rtype: None
        """
        self.ptrx().setTetSpecConc(steps.tetrahedron_global_id(idx), to_std_string(s), c)

    def getTetSpecClamped(self, steps.index_t idx, str s):
        """
        Returns whether the concentration of species s in a tetrahedron
        remains constant over time (unless changed explicitly).

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getTetSpecClamped(steps.tetrahedron_global_id(idx), to_std_string(s))

    def setTetSpecClamped(self, steps.index_t idx, str s, bool buf):
        """
        Sets clamping of species s in a tetrahedron on or off.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str
        :param buf: Flag to turn the clamping of species on or off.
        :type buf: bool


        :rtype: None
        """
        self.ptrx().setTetSpecClamped(steps.tetrahedron_global_id(idx), to_std_string(s), buf)

    def getTetReacK(self, steps.index_t idx, str r):
        """
        Returns the macroscopic reaction constant of reaction r in a
        tetrahedron (units vary with order of reaction).

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTetReacK(steps.tetrahedron_global_id(idx), to_std_string(r))

    def setTetReacK(self, steps.index_t idx, str r, double kf):
        """
        Sets the macroscopic reaction constant of reaction r in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param r: Name of the reaction.
        :type r: str
        :param kf: Rate constant of the reaction.
        :type kf: float


        :rtype: None
        """
        self.ptrx().setTetReacK(steps.tetrahedron_global_id(idx), to_std_string(r), kf)

    def getTetReacActive(self, steps.index_t idx, str r):
        """
        Returns whether reaction r in a tetrahedron is active or not

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param r: Name of the reaction.
        :type r: str


        :rtype: bool
        """
        return self.ptrx().getTetReacActive(steps.tetrahedron_global_id(idx), to_std_string(r))

    def setTetReacActive(self, steps.index_t idx, str r, bool act):
        """
        Activates/deactivates reaction r in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param r: Name of the reaction.
        :type r: str
        :param act: Flag to activate or deactivate the reaction.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setTetReacActive(steps.tetrahedron_global_id(idx), to_std_string(r), act)

    def getTetDiffD(self, steps.index_t idx, str d, steps.index_t direction_tet=UNKNOWN_TET):
        """
        Returns the diffusion constant of diffusion rule d in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param d: Name of the deffusion.
        :type d: str
        :param direction_tet: Tetrahedron index which specifies diffusion direction.
        :type direction_tet: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTetDiffD(steps.tetrahedron_global_id(idx), to_std_string(d), steps.tetrahedron_global_id(direction_tet))

    def setTetDiffD(self, steps.index_t idx, str d, double dk, steps.index_t direction_tet=UNKNOWN_TET):
        """
        Sets the diffusion constant of diffusion rule d in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param d: Name of the diffusion.
        :type d: str
        :param dk: Rate constant of the diffusion.
        :type dk: float
        :param direction_tet: Tetrahedron index which the diffusion towards.
        :type direction_tet: steps.index_t


        :rtype: None
        """
        self.ptrx().setTetDiffD(steps.tetrahedron_global_id(idx), to_std_string(d), dk, steps.tetrahedron_global_id(direction_tet))

    def getTetDiffActive(self, steps.index_t idx, str d):
        """
        Returns whether diffusion rule d in a tetrahedron is active or not.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param d: Name of the diffusion.
        :type d: str


        :rtype: bool
        """
        return self.ptrx().getTetDiffActive(steps.tetrahedron_global_id(idx), to_std_string(d))

    def setTetDiffActive(self, steps.index_t idx, str d, bool act):
        """
        Activates/deactivates diffusion rule d in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param d: Name of the diffusion.
        :type d: str
        :param act: Flag to activate / deactivate the diffusion.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setTetDiffActive(steps.tetrahedron_global_id(idx), to_std_string(d), act)

    def getTetReacC(self, steps.index_t idx, str r):
        """
        Returns c_mu, the mesoscopic reaction constant of reaction r in
        a tetrahedron

        :param idx: Index of the diffusion.
        :type idx: steps.index_t
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTetReacC(steps.tetrahedron_global_id(idx), to_std_string(r))

    def getTetReacH(self, steps.index_t idx, str r):
        """
        Returns h_mu, the distinct number of ways in which reaction r can
        occur in a tetrahedron, by computing the product of its reactants.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param r:
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTetReacH(steps.tetrahedron_global_id(idx), to_std_string(r))

    def getTetReacA(self, steps.index_t idx, str r):
        """
        Returns the propensity, a_mu, of reaction r in a tetrahedron.
        The propensity value gives the probability per unit time that this
        reaction will occur in the current state.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTetReacA(steps.tetrahedron_global_id(idx), to_std_string(r))

    def getTetDiffA(self, steps.index_t idx, str d):
        """
        Returns the propensity, a_mu of diffusion rule d in a tetrahedron.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t
        :param d: Name of the diffusion.
        :type d: str


        :rtype: float
        """
        return self.ptrx().getTetDiffA(steps.tetrahedron_global_id(idx), to_std_string(d))

    def getTetV(self, steps.index_t idx):
        """
        Returns the potential of tetrahedron in Volts.

        :param idx: Index of the tetrahedron.
        :type idx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTetV(steps.tetrahedron_global_id(idx))

    def setTetV(self, steps.index_t idx, double v):
        """
        Set the potential of tetrahedron.

        :param idx: Index of the tetrahedron
        :type idx: steps.index_t
        :param v:
        :type v: float


        :rtype: None
        """
        self.ptrx().setTetV(steps.tetrahedron_global_id(idx), v)

    def getTetVClamped(self, steps.index_t idx):
        """
        Returns whether the potential of tetrahedron is clamped over time
        (unless changed explicitly)

        :param idx: Index of the tetrahedron
        :type idx: steps.index_t


        :rtype: bool
        """
        return self.ptrx().getTetVClamped(steps.tetrahedron_global_id(idx))

    def setTetVClamped(self, steps.index_t idx, bool cl):
        """
        Sets voltage clamp in tetrahedron.

        :param idx: Index of the tetrahedron
        :type idx: steps.index_t
        :param cl: Flag to turn the clamping on or off.
        :type cl: bool


        :rtype: None
        """
        self.ptrx().setTetVClamped(steps.tetrahedron_global_id(idx), cl)

    def getPatchArea(self, str p):
        """
        Returns the area of patch p (in m^2)

        :param p: Name of the patch.
        :type p: str


        :rtype: float
        """
        return self.ptrx().getPatchArea(to_std_string(p))

    def getPatchSpecCount(self, str p, str s):
        """
        Returns the number of molecules of species s in patch p.
        NOTE: in a mesh-based simulation, the total count is computed as
        the sum of the counts in all triangles of the patch.

        :param p: Name of the path.
        :type p: str
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getPatchSpecCount(to_std_string(p), to_std_string(s))

    def setPatchSpecCount(self, str p, str s, double n):
        """
        Sets the number of molecules of species s in patch p.
        NOTE: in a mesh-based simulation, the total amount is equally divided
        over all triangles in the patch.

        :param p: Name of the patch.
        :type p: str
        :param s: Name of the species.
        :type s: str
        :param n: Number of molecules of species.
        :type n: float


        :rtype: None
        """
        self.ptrx().setPatchSpecCount(to_std_string(p), to_std_string(s), n)

    def getPatchSpecAmount(self, str p, str s):
        """
        Returns the amount (in mols) of species s in patch p.
        NOTE: in a mesh-based simulation, the total amount is computed as
        the sum of the amounts in all triangles of the patch.

        :param p: Name of the patch.
        :type p: str
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getPatchSpecAmount(to_std_string(p), to_std_string(s))

    def setPatchSpecAmount(self, str p, str s, double a):
        """
        Sets the amount (in mols) of species s in patch p.
        NOTE: in a mesh-based simulation, the total amount is equally divided
        over all triangles in the patch.

        :param p: Name of the patch.
        :type p: str
        :param s: Name of the species.
        :type s: str
        :param a: Amount of the species.
        :type a: float


        :rtype: None
        """
        self.ptrx().setPatchSpecAmount(to_std_string(p), to_std_string(s), a)

    def getPatchSpecClamped(self, str p, str s):
        """
        Returns whether the count of species s in patch p remains constant.
        over time (unless changed explicitly).
        NOTE: this method will only return true if the species has been
        clamped in all triangles in the patch.

        :param p: Name of the patch.
        :type p: str
        :param s: Name of the species.
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getPatchSpecClamped(to_std_string(p), to_std_string(s))

    def setPatchSpecClamped(self, str p, str s, bool buf):
        """
        Turns clamping of species in patch p on or off.
        NOTE: in a mesh-based simulation, this method turns clamping on/off
        in all triangles in the patch.

        :param p: Name of the patch.
        :type p: str
        :param s: Name of the species.
        :type s: str
        :param buf: Flag to turn clamping of species on /off.
        :type buf: bool


        :rtype: None
        """
        self.ptrx().setPatchSpecClamped(to_std_string(p), to_std_string(s), buf)

    def getPatchSReacK(self, str p, str r):
        """
        Returns the macroscopic reaction constant of surface reaction r
        in patch p.
        NOTE: in a mesh-based simulation, the value is computed as the
        area-weighted sum of the reaction constants in all triangles of
        the patch.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getPatchSReacK(to_std_string(p), to_std_string(r))

    def setPatchSReacK(self, str p, str r, double kf):
        """
        Sets the macroscopic reaction constant of surface reaction r
        in patch p.
        NOTE: in a mesh-based simulation this method changes the reaction
        constant equally in all triangles of the patch.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str
        :param kf: Rate constant of the reaction.
        :type kf: float


        :rtype: None
        """
        self.ptrx().setPatchSReacK(to_std_string(p), to_std_string(r), kf)

    def getPatchSReacActive(self, str p, str r):
        """
        Returns whether surface reaction r in patch p is active or not.
        NOTE: in a mesh-based simulation, only returns false when the
        reaction has been inactivated in all triangles.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: bool
        """
        return self.ptrx().getPatchSReacActive(to_std_string(p), to_std_string(r))

    def setPatchSReacActive(self, str p, str r, bool a):
        """
        Activate or inactivate surface reaction r in patch p.
        NOTE: in a mesh-based simulation, activation/inactivation of a
        surface reaction turns it on/off in all triangles.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str
        :param a: Flag to activate / deactivate the reaction.
        :type a: bool


        :rtype: None
        """
        self.ptrx().setPatchSReacActive(to_std_string(p), to_std_string(r), a)

    def getPatchSReacC(self, str p, str r):
        """
        Returns c_mu, the mesoscopic reaction constant of surface reaction r
        in patch p.
        NOTE: in a mesh_based simulation, the mesoscopic reaction constant
        is computed as the sum of the mesoscopic reaction constants from all
        triangles of the patch.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reacton.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getPatchSReacC(to_std_string(p), to_std_string(r))

    def getPatchSReacH(self, str p, str r):
        """
        Returns h_mu, the distinct number of ways in which a surface reaction
        r can occur in patch p, by computing the product of its reactants.
        NOTE: in a mesh-based simulation, it returns the sum of the h_mu's
        over all triangles triangles of the patch.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getPatchSReacH(to_std_string(p), to_std_string(r))

    def getPatchSReacA(self, str p, str r):
        """
        Returns the propensity, a_mu of surface reaction r in patch p.
        This propensity value gives the probability per unit time that this
        surface reaction will occur in the current state.
        NOTE: in a mesh-based simulation, a_mu is computed as the sum of the
        a_mu in all triangles in the patch.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getPatchSReacA(to_std_string(p), to_std_string(r))

    def getPatchSReacExtent(self, str p, str r):
        """
        Returns the extent of surface reaction r in patch p.
        NOTE: in a mesh-based simulation, returns the sum of the extents in
        all triangles of the patch.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: unsigned long long
        """
        return self.ptrx().getPatchSReacExtent(to_std_string(p), to_std_string(r))

    def resetPatchSReacExtent(self, str p, str r):
        """
        Resets the extent of surface reaction r in patch p to zero.
        NOTE: in a mesh-based simulation, resets the extents of the
        surface reaction in all triangles of the patch.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the reaction.
        :type r: str


        :rtype: None
        """
        self.ptrx().resetPatchSReacExtent(to_std_string(p), to_std_string(r))

    def getPatchVDepSReacActive(self, str p, str vsr):
        """
        Returns whether voltage-dependent surface reaction vsr in patch p is
        active or not.
        NOTE: only returns false when the voltage-dependent surface
        reaction has been inactivated in all triangles.

        :param p: Name of the patch.
        :type p: str
        :param vsr: Name of the voltage-dependent surface reaction.
        :type vsr: str


        :rtype: bool
        """
        return self.ptrx().getPatchVDepSReacActive(to_std_string(p), to_std_string(vsr))

    def setPatchVDepSReacActive(self, str p, str vsr, bool a):
        """
        Activate or inactivate voltage-dependent surface reaction vsr in patch
        p.
        NOTE: activation/inactivation of a voltage-dependent
        surface reaction turns it on/off in all triangles.

        :param p: Name of the patch.
        :type p: str
        :param vsr: Name of the voltage-dependent surface reaction.
        :type vsr: str
        :param a: Flag to activate / deactivate the reaction.
        :type a: bool


        :rtype: None
        """
        self.ptrx().setPatchVDepSReacActive(to_std_string(p), to_std_string(vsr), a)

    def getPatchRaftCount(self, str p, str r):
        """
        Returns the number of rafts r in patch p

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the raft.
        :type r: str


        :rtype: int
        """
        return self.ptrx().getPatchRaftCount(to_std_string(p), to_std_string(r))

    def setPatchRaftCount(self, str p, str r, uint n):
        """
        Set the number of rafts r in patch p

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the raft.
        :type r: str
        :param n: Number of rafts
        :type n: int


        :rtype: None
        """
        self.ptrx().setPatchRaftCount(to_std_string(p), to_std_string(r), n)

    def getSingleRaftPos(self, str r, steps.index_t raft_unique_index):
        """
        Returns the position of raft of type r with unique index
        raft_unique_index (in cartesian coordinates)

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t


        :rtype: List[float]
        """
        return [_val1 for _val1 in self.ptrx().getSingleRaftPos(to_std_string(r), steps.raft_individual_id(raft_unique_index))]

    def getPatchRaftSpecCountDict(self, str p, str r, str s):
        """
        Returns the count of species s in raft r in patch p. Return is a map,
        raft_unique_index : count of s

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the raft.
        :type r: str
        :param s: Name of the species
        :type s: str


        :rtype: Dict[steps.index_t, int]
        """
        return { _pair1.first.get(): _pair1.second for _pair1 in self.ptrx().getPatchRaftSpecCountDict(to_std_string(p), to_std_string(r), to_std_string(s)) }

    def getPatchRaftSpecCount(self, str p, str r, str s):
        """
        Returns the summed count of species s in raft r in patch p.

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the raft.
        :type r: str
        :param s: Name of the species
        :type s: str


        :rtype: int
        """
        return self.ptrx().getPatchRaftSpecCount(to_std_string(p), to_std_string(r), to_std_string(s))

    def getSingleRaftSpecCount(self, str r, steps.index_t raft_unique_index, str s):
        """
        Get the count of species s on raft of type r and unique index
        raft_unique_index.

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str


        :rtype: int
        """
        return self.ptrx().getSingleRaftSpecCount(to_std_string(r), steps.raft_individual_id(raft_unique_index), to_std_string(s))

    def setSingleRaftSpecCount(self, str r, steps.index_t raft_unique_index, str s, uint count):
        """
        Set the count of species s on raft of type r and unique index
        raft_unique_index.

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t
        :param s: Name of the species
        :type s: str
        :param count: The number of the species
        :type count: int


        :rtype: None
        """
        self.ptrx().setSingleRaftSpecCount(to_std_string(r), steps.raft_individual_id(raft_unique_index), to_std_string(s), count)

    def getSingleRaftImmobility(self, str r, steps.index_t raft_unique_index):
        """
        Get the 'immobility' of raft of type r and unique index
        raft_unique_index. All non-zero numbers mean raft is
        immobile.

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t


        :rtype: int
        """
        return self.ptrx().getSingleRaftImmobility(to_std_string(r), steps.raft_individual_id(raft_unique_index))

    def getSingleRaftRaftEndocytosisK(self, str r, steps.index_t raft_unique_index, str rendo):
        """
        Get the raft endocytosis rate of a raft

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t
        :param rendo: Name of the raft endocytosis
        :type rendo: str


        :rtype: float
        """
        return self.ptrx().getSingleRaftRaftEndocytosisK(to_std_string(r), steps.raft_individual_id(raft_unique_index), to_std_string(rendo))

    def setSingleRaftRaftEndocytosisK(self, str r, steps.index_t raft_unique_index, str rendo, double k):
        """
        Set the raft endocytosis rate of a raft

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t
        :param rendo: Name of the raft endocytosis
        :type rendo: str
        :param k: Rate.
        :type k: float


        :rtype: None
        """
        self.ptrx().setSingleRaftRaftEndocytosisK(to_std_string(r), steps.raft_individual_id(raft_unique_index), to_std_string(rendo), k)


    def setSingleRaftSReacActive(self, str r, steps.index_t raft_unique_index, str rsreac, bool active):
        """
        Activate or de-activate a raft surface reaction on a raft

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t
        :param rsreac: Name of the raft surface reaction
        :type rsreac: str
        :param active: Whether the raftsreac should be active
        :type active: bool


        :rtype: None
        """
        self.ptrx().setSingleRaftSReacActive(to_std_string(r), steps.raft_individual_id(raft_unique_index), to_std_string(rsreac), active)


    def getSingleRaftSReacActive(self, str r, steps.index_t raft_unique_index, str rsreac):
        """
        Get whether raft surface reaction on a raft is active or not

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft.
        :type raft_unique_index: steps.index_t
        :param rsreac: Name of the raft surface reaction
        :type rsreac: str


        :rtype: bool
        """
        return self.ptrx().getSingleRaftSReacActive(to_std_string(r), steps.raft_individual_id(raft_unique_index), to_std_string(rsreac))


    def setPatchEndocyticZoneEndocytosisActive(self, str patch, str zone, str endo, bool active):
        """
        Activate or de-activate endocytosis reaction in an endocytic zone.

        :param patch:
        :type patch: str
        :param zone: Name of the zone
        :type zone: str
        :param endo: Name of the endocytosis
        :type endo: str
        :param active: Whether the endocytosis should be active
        :type active: bool


        :rtype: None
        """
        self.ptrx().setPatchEndocyticZoneEndocytosisActive(to_std_string(patch), to_std_string(zone), to_std_string(endo), active)

    def setPatchEndocyticZoneEndocytosisK(self, str patch, str zone, str endo, double k):
        """
        Change endocytosis rate in an endocytic zone.

        :param patch:
        :type patch: str
        :param zone: Name of the zone
        :type zone: str
        :param endo: Name of the endocytosis
        :type endo: str
        :param k: Rate of the endocytosis in this endocytic zone
        :type k: double


        :rtype: None
        """
        self.ptrx().setPatchEndocyticZoneEndocytosisK(to_std_string(patch), to_std_string(zone), to_std_string(endo), k)

    def getPatchEndocyticZoneEndocytosisExtent(self, str patch, str zone, str endo):
        """
        Returns the endocytosis extent in an endocytic zone

        :param patch: Name of the patch
        :type patch: str
        :param zone: Name of the zone
        :type zone: str
        :param endo: Name of the endocytosis.
        :type endo: str


        :rtype: int
        """
        return self.ptrx().getPatchEndocyticZoneEndocytosisExtent(to_std_string(patch), to_std_string(zone), to_std_string(endo))

    def getPatchEndocyticZoneEndocytosisEvents(self, str patch, str zone, str endo):
        """
        Get the vesicle endocytosis events in an endocytic zone that happened since last call.

        :param patch: Name of the patch
        :type patch: str
        :param zone: Name of the zone
        :type zone: str
        :param endo: Name of the endocytosis.
        :type endo: str


        :rtype: List[Tuple[float, steps.index_t, steps.index_t]]
        """
        return [(event.time, event.tidx.get(), event.vidx.get()) for event in self.ptrx().getPatchEndocyticZoneEndocytosisEvents(to_std_string(patch), to_std_string(zone), to_std_string(endo))]

    def getPatchRaftIndices(self, str p, str r):
        """
        Returns all the unique indices of rafts r currently in patch p

        :param p: Name of the patch.
        :type p: str
        :param r: Name of the raft.
        :type r: str


        :rtype: List[steps.index_t]
        """
        return [_val1.get() for _val1 in self.ptrx().getPatchRaftIndices(to_std_string(p), to_std_string(r))]

    def getSingleRaftPatch(self, str r, steps.index_t raft_unique_index):
        """
        Get the patch of raft of type r and unique index raft_unique_index

        :param r: Name of the raft.
        :type r: str
        :param raft_unique_index: Unique index of the individual raft
        :type raft_unique_index: steps.index_t


        :rtype: str
        """
        return from_std_string(self.ptrx().getSingleRaftPatch(to_std_string(r), steps.raft_individual_id(raft_unique_index)))

    def setDiffBoundarySpecDiffusionActive(self, str db, str s, bool act):
        """
        Activate or inactivate diffusion across a diffusion boundary for a
        species.

        :param db: Name of the diffusion boundary.
        :type db: str
        :param s: Name of the species.
        :type s: str
        :param act: Bool to activate (true) or inactivate (false) diffusion.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setDiffBoundarySpecDiffusionActive(to_std_string(db), to_std_string(s), act)

    def getDiffBoundarySpecDiffusionActive(self, str db, str s):
        """
        Returns whether diffusion is active across a diffusion boundary for a
        species.

        :param db: Name of the diffusion boundary.
        :type db: str
        :param s: Name of the species.
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getDiffBoundarySpecDiffusionActive(to_std_string(db), to_std_string(s))

    def setDiffBoundarySpecDcst(self, str db, str s, double dcst, str direction_comp=''):
        """
        Set the diffusion constant across a diffusion boundary.

        :param db: Name of the diffusion boundary.
        :type db: str
        :param s: Name of the species.
        :type s: str
        :param dcst: diffusion constant.
        :type dcst: float
        :param direction_comp: ID of the compartment which the diffusion towards to.
        :type direction_comp: str


        :rtype: None
        """
        self.ptrx().setDiffBoundarySpecDcst(to_std_string(db), to_std_string(s), dcst, to_std_string(direction_comp))

    def setSDiffBoundarySpecDiffusionActive(self, str sdb, str s, bool act):
        """
        Activate or inactivate diffusion across a surface diffusion boundary
        for a species.

        :param sdb: Name of the surface diffusion boundary.
        :type sdb: str
        :param s: Name of the species.
        :type s: str
        :param act: Bool to activate (true) or inactivate (false) diffusion.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setSDiffBoundarySpecDiffusionActive(to_std_string(sdb), to_std_string(s), act)

    def getSDiffBoundarySpecDiffusionActive(self, str sdb, str s):
        """
        Returns whether diffusion is active across a surface diffusion boundary
        for a species.

        :param sdb: Name of the surface diffusion boundary.
        :type sdb: str
        :param s: Name of the species.
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getSDiffBoundarySpecDiffusionActive(to_std_string(sdb), to_std_string(s))

    def setSDiffBoundarySpecDcst(self, str sdb, str s, double dcst, str direction_patch=''):
        """
        Set the diffusion constant across a surface diffusion boundary.

        :param sdb: Name of the surface diffusion boundary.
        :type sdb: str
        :param s: Name of the species.
        :type s: str
        :param dcst: diffusion constant.
        :type dcst: float
        :param direction_patch: ID of the patch which the diffusion is towards to.
        :type direction_patch: str


        :rtype: None
        """
        self.ptrx().setSDiffBoundarySpecDcst(to_std_string(sdb), to_std_string(s), dcst, to_std_string(direction_patch))

    def getTriArea(self, steps.index_t idx):
        """
        Returns the area of the triangle (in m^2).

        :param idx: Index of the triangle.
        :type idx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTriArea(steps.triangle_global_id(idx))

    def setTriArea(self, steps.index_t idx, double area):
        """
        Set the area (in m^2) of the triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param area: Area of teh triangle.
        :type area: float


        :rtype: None
        """
        self.ptrx().setTriArea(steps.triangle_global_id(idx), area)

    def getTriSpecDefined(self, steps.index_t idx, str s):
        """
        Returns whether species s is defined in a triangle

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getTriSpecDefined(steps.triangle_global_id(idx), to_std_string(s))

    def getTriSpecCount(self, steps.index_t idx, str s):
        """
        Returns the number of molecules of species s in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getTriSpecCount(steps.triangle_global_id(idx), to_std_string(s))

    def setTriSpecCount(self, steps.index_t idx, str s, double n):
        """
        Sets the number of molecules of species s in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str
        :param n: Number of molecules of the species.
        :type n: float


        :rtype: None
        """
        self.ptrx().setTriSpecCount(steps.triangle_global_id(idx), to_std_string(s), n)

    def getTriSpecAmount(self, steps.index_t idx, str s):
        """
        Returns the amount (in mols) of species s in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str


        :rtype: float
        """
        return self.ptrx().getTriSpecAmount(steps.triangle_global_id(idx), to_std_string(s))

    def setTriSpecAmount(self, steps.index_t idx, str s, double m):
        """
        Sets the amount (in mols) of species s in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param s: Name of the species.
        :type s: str
        :param m: Amount of the species.
        :type m: float


        :rtype: None
        """
        self.ptrx().setTriSpecAmount(steps.triangle_global_id(idx), to_std_string(s), m)

    def getTriSpecClamped(self, steps.index_t idx, str s):
        """
        Returns whether the number of molecules of species s in a triangle
        remains constant over time (unless changed explicitly)

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param s: name of the species.
        :type s: str


        :rtype: bool
        """
        return self.ptrx().getTriSpecClamped(steps.triangle_global_id(idx), to_std_string(s))

    def setTriSpecClamped(self, steps.index_t idx, str s, bool buf):
        """
        Sets clamping of species s in a triangle on or off.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param s: name of the species.
        :type s: str
        :param buf: Flag to set clamping of species on /off.
        :type buf: bool


        :rtype: None
        """
        self.ptrx().setTriSpecClamped(steps.triangle_global_id(idx), to_std_string(s), buf)

    def getTriSReacK(self, steps.index_t idx, str r):
        """
        Returns the macroscopic reaction constant of surface reaction r
        in a triangle (units vary with order of reaction).

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTriSReacK(steps.triangle_global_id(idx), to_std_string(r))

    def setTriSReacK(self, steps.index_t idx, str r, double kf):
        """
        Sets the macroscopic reaction constant of surface reaction r in
        a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: name of the reaction.
        :type r: str
        :param kf: Rate constant of the reaction.
        :type kf: float


        :rtype: None
        """
        self.ptrx().setTriSReacK(steps.triangle_global_id(idx), to_std_string(r), kf)

    def getTriSReacActive(self, steps.index_t idx, str r):
        """
        Returns whether surface reaction r in a triangle is active or not.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: name of the reaction.
        :type r: str


        :rtype: bool
        """
        return self.ptrx().getTriSReacActive(steps.triangle_global_id(idx), to_std_string(r))

    def setTriSReacActive(self, steps.index_t idx, str r, bool act):
        """
        Activates/inactivates surface reaction r in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: name of the reaction.
        :type r: str
        :param act: Flag to activate / deactivate the reaction.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setTriSReacActive(steps.triangle_global_id(idx), to_std_string(r), act)

    def getTriSReacC(self, steps.index_t idx, str r):
        """
        Returns c_mu, the mesoscopic reaction constant of surface reaction r
        in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTriSReacC(steps.triangle_global_id(idx), to_std_string(r))

    def getTriSReacH(self, steps.index_t idx, str r):
        """
        Returns h_mu, the distinct number of ways in which surface reaction r
        can occur in a triangle, by computing the product of it's reactants.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTriSReacH(steps.triangle_global_id(idx), to_std_string(r))

    def getTriSReacA(self, steps.index_t idx, str r):
        """
        Returns the propensity, a_mu, of surface reaction r in a triangle.
        The propensity value gives the probability per unit time that this
        surface reaction will occur in the current state.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: name of the reaction.
        :type r: str


        :rtype: float
        """
        return self.ptrx().getTriSReacA(steps.triangle_global_id(idx), to_std_string(r))

    def getTriDiffD(self, steps.index_t idx, str d, uint direction_tri=UNKNOWN_TRI):
        """
        outdate function

        :param idx:
        :type idx: steps.index_t
        :param d:
        :type d: str
        :param direction_tri:
        :type direction_tri: int


        :rtype: float
        """
        return self.ptrx().getTriDiffD(steps.triangle_global_id(idx), to_std_string(d), direction_tri)

    def getTriSDiffD(self, steps.index_t idx, str d, steps.index_t direction_tri=UNKNOWN_TRI):
        """
        Returns the diffusion constant of diffusion rule d in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param d: Name of the diffusion.
        :type d: str
        :param direction_tri:
        :type direction_tri: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTriSDiffD(steps.triangle_global_id(idx), to_std_string(d), steps.triangle_global_id(direction_tri))

    def setTriDiffD(self, steps.index_t idx, str d, double dk, steps.index_t direction_tri=UNKNOWN_TRI):
        """
        outdated function

        :param idx:
        :type idx: steps.index_t
        :param d:
        :type d: str
        :param dk:
        :type dk: float
        :param direction_tri:
        :type direction_tri: steps.index_t


        :rtype: None
        """
        self.ptrx().setTriDiffD(steps.triangle_global_id(idx), to_std_string(d), dk, steps.triangle_global_id(direction_tri))

    def setTriSDiffD(self, steps.index_t idx, str d, double dk, steps.index_t direction_tri=UNKNOWN_TRI):
        """
        Sets the diffusion constant of diffusion rule d on a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param d: Name of the diffusion.
        :type d: str
        :param dk: Rate constant of the diffusion.
        :type dk: float
        :param direction_tri: Triangle index which the diffusion towards
        :type direction_tri: steps.index_t


        :rtype: None
        """
        self.ptrx().setTriSDiffD(steps.triangle_global_id(idx), to_std_string(d), dk, steps.triangle_global_id(direction_tri))

    def getTriRaftCount(self, steps.index_t idx, str r):
        """
        Returns the number of rafts in a triangle

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: Name of the raft.
        :type r: str


        :rtype: int
        """
        return self.ptrx().getTriRaftCount(steps.triangle_global_id(idx), to_std_string(r))

    def setTriRaftCount(self, steps.index_t idx, str r, uint n):
        """
        Returns the number of rafts in a triangle

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: Name of the raft.
        :type r: str
        :param n: Number of rafts
        :type n: int


        :rtype: None
        """
        self.ptrx().setTriRaftCount(steps.triangle_global_id(idx), to_std_string(r), n)

    def addTriRaft(self, steps.index_t idx, str r):
        """
        Add one raft to a triangle

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param r: Name of the raft.
        :type r: str


        :rtype: steps.index_t
        """
        return self.ptrx().addTriRaft(steps.triangle_global_id(idx), to_std_string(r)).get()

    def getTriExocytosisActive(self, steps.index_t idx, str er):
        """
        Returns whether exocytotic reaction in a triangle is active or not.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param er: name of the exocytotic reaction.
        :type er: str


        :rtype: bool
        """
        return self.ptrx().getTriExocytosisActive(steps.triangle_global_id(idx), to_std_string(er))

    def setTriExocytosisActive(self, steps.index_t idx, str er, bool act):
        """
        Activates/inactivates exocytotic reaction in a triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param er: name of the exocytotic reaction.
        :type er: str
        :param act: Flag to activate / deactivate the exocytotic reaction.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setTriExocytosisActive(steps.triangle_global_id(idx), to_std_string(er), act)

    def getTriV(self, steps.index_t idx):
        """
        Returns the potential of triangle in Volts.

        :param idx: Index of the triangle.
        :type idx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTriV(steps.triangle_global_id(idx))

    def setTriV(self, steps.index_t idx, double v):
        """
        Set the potential of triangle.

        :param idx: Index of the triangle
        :type idx: steps.index_t
        :param v:
        :type v: float


        :rtype: None
        """
        self.ptrx().setTriV(steps.triangle_global_id(idx), v)

    def getTriVClamped(self, steps.index_t idx):
        """
        Returns whether the potential of triangle is clamped over time
        (unless changed explicitly).

        :param idx: Index of the triangle
        :type idx: steps.index_t


        :rtype: bool
        """
        return self.ptrx().getTriVClamped(steps.triangle_global_id(idx))

    def setTriVClamped(self, steps.index_t idx, bool cl):
        """
        Sets voltage clamp in triangle.

        :param idx: Index of the triangle
        :type idx: steps.index_t
        :param cl: Flag to turn the clamping on or off.
        :type cl: bool


        :rtype: None
        """
        self.ptrx().setTriVClamped(steps.triangle_global_id(idx), cl)

    def getTriOhmicI(self, steps.index_t idx, str oc):
        """
        Returns the ohmic current of triangle in amperes.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param oc: name of the ohmic current
        :type oc: str


        :rtype: float
        """
        return self.ptrx().getTriOhmicI(steps.triangle_global_id(idx), to_std_string(oc))

    def getTriGHKI(self, steps.index_t idx, str ghk):
        """
        Returns the GHK current of triangle in amperes.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param ghk: name of the ghk current
        :type ghk: str


        :rtype: float
        """
        return self.ptrx().getTriGHKI(steps.triangle_global_id(idx), to_std_string(ghk))

    def getTriI(self, steps.index_t idx):
        """
        Returns the current of a triangle in amperes from the last EField
        calculation step.

        :param idx: Index of the triangle
        :type idx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTriI(steps.triangle_global_id(idx))

    def getTriIClamp(self, steps.index_t idx):
        """
        Gets current injection to triangle.

        :param idx: Index of the triangle
        :type idx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getTriIClamp(steps.triangle_global_id(idx))

    def setTriIClamp(self, steps.index_t idx, double i):
        """
        Sets current injection to triangle.
        Will be assumed to be constant for one EField DT

        :param idx: Index of the triangle
        :type idx: steps.index_t
        :param i:
        :type i: float


        :rtype: None
        """
        self.ptrx().setTriIClamp(steps.triangle_global_id(idx), i)

    def getTriVDepSReacActive(self, steps.index_t idx, str vsr):
        """
        Returns whether voltage-dependent surface reaction vsr in a triangle is
        active or not.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param vsr: name of the voltage-dependent surface reaction.
        :type vsr: str


        :rtype: bool
        """
        return self.ptrx().getTriVDepSReacActive(steps.triangle_global_id(idx), to_std_string(vsr))

    def setTriVDepSReacActive(self, steps.index_t idx, str vsr, bool act):
        """
        Activates/inactivates voltage-dependent surface reaction vsr in a
        triangle.

        :param idx: Index of the triangle.
        :type idx: steps.index_t
        :param vsr: name of the voltage-dependent surface reaction.
        :type vsr: str
        :param act: Flag to activate / deactivate the reaction.
        :type act: bool


        :rtype: None
        """
        self.ptrx().setTriVDepSReacActive(steps.triangle_global_id(idx), to_std_string(vsr), act)

    def setTriCapac(self, steps.index_t idx, double cm):
        """
        Set the specific capacitance of a triangle surface element.

        :param idx: Index of the triangle surface element
        :type idx: steps.index_t
        :param cm: Specific membrane capacitance (farad / m^2)
        :type cm: float


        :rtype: None
        """
        self.ptrx().setTriCapac(steps.triangle_global_id(idx), cm)

    def getVertV(self, steps.index_t vidx):
        """
        Returns the potential of vertex in Volts.

        :param vidx: Index of the vertex.
        :type vidx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getVertV(steps.vertex_id_t(vidx))

    def setVertV(self, steps.index_t vidx, double v):
        """
        Set the potential of vertex.

        :param vidx: Index of the vertex
        :type vidx: steps.index_t
        :param v:
        :type v: float


        :rtype: None
        """
        self.ptrx().setVertV(steps.vertex_id_t(vidx), v)

    def getVertVClamped(self, steps.index_t vidx):
        """
        Returns whether the potential of vertex is clamped over time
        (unless changed explicitly).

        :param vidx: Index of the vertex
        :type vidx: steps.index_t


        :rtype: bool
        """
        return self.ptrx().getVertVClamped(steps.vertex_id_t(vidx))

    def setVertVClamped(self, steps.index_t vidx, bool cl):
        """
        Sets voltage clamp in vertex.

        :param vidx: Index of the vertex
        :type vidx: steps.index_t
        :param cl: Flag to turn the clamping on or off.
        :type cl: bool


        :rtype: None
        """
        self.ptrx().setVertVClamped(steps.vertex_id_t(vidx), cl)

    def getVertIClamp(self, steps.index_t vidx):
        """
        Gets current injection to vertex.

        :param vidx: Index of the vertex
        :type vidx: steps.index_t


        :rtype: float
        """
        return self.ptrx().getVertIClamp(steps.vertex_id_t(vidx))

    def setVertIClamp(self, steps.index_t vidx, double i):
        """
        Sets current injection to vertex.
        Will be assumed to be constant for one EField DT

        :param vidx: Index of the vertex
        :type vidx: steps.index_t
        :param i:
        :type i: float


        :rtype: None
        """
        self.ptrx().setVertIClamp(steps.vertex_id_t(vidx), i)

    def setMembPotential(self, str m, double v):
        """
        Set the electric potential of the membrane, including all nodes
        in the conduction volume.

        :param m: Name of the membrane
        :type m: str
        :param v: Potential (volts)
        :type v: float


        :rtype: None
        """
        self.ptrx().setMembPotential(to_std_string(m), v)

    def setMembCapac(self, str m, double cm):
        """
        Set the specific membrane capacitance of the membrane

        :param m: Name of the membrane
        :type m: str
        :param cm: Specific membrane capacitance (farad / m^2)
        :type cm: float


        :rtype: None
        """
        self.ptrx().setMembCapac(to_std_string(m), cm)

    def setMembVolRes(self, str m, double ro):
        """
        Set the bulk electrical resistivity of the section of the mesh
        representing the volume conductor

        :param m: Name of the membrane
        :type m: str
        :param ro: Electrical resistivity (ohm.m)
        :type ro: float


        :rtype: None
        """
        self.ptrx().setMembVolRes(to_std_string(m), ro)

    def setMembRes(self, str m, double ro, double vrev):
        """
        Sets the surface electrical resistivity ro (in ohm.m^2) of the membrane with string identifier memb. Reversal potential vrev is required in Volts.

        :param m: Name of the membrane
        :type m: str
        :param ro: membrane resistivity (ohm.m^2)
        :type ro: float
        :param vrev: Reversal potential (Volts)
        :type vrev: float


        :rtype: None
        """
        self.ptrx().setMembRes(to_std_string(m), ro, vrev)

    def setOutputSync(self, bool enable_sync, int output_rank):
        """
        Set if the outputs are synced across all ranks.
        If not, set the output rank.

        :param enable_sync:
        :type enable_sync: bool
        :param output_rank:
        :type output_rank: int


        :rtype: None
        """
        self.ptrx().setOutputSync(enable_sync, output_rank)

    def getOutputSyncStatus(self, ):
        """
        Get if the outputs are synced across all ranks.


        :rtype: bool
        """
        return self.ptrx().getOutputSyncStatus()

    def getOutputSyncRank(self, ):
        """
        Get the rank id for data output.


        :rtype: int
        """
        return self.ptrx().getOutputSyncRank()

    def _setMaxWalkDistFact(self, double fact):
        """
        Set the factor for maximum walking distance in checkPos algorithm

        :param fact:
        :type fact: float

        :rtype: None
        """
        self.ptrx().setMaxWalkDistFact(fact)

    def _setMinNbTetVisited(self, int n):
        """
        Set the minimum number of tetrahedrons explored in checkPos algorithm

        :param fact:
        :type fact: int


        :rtype: None
        """
        self.ptrx().setMinNbTetVisited(n)


    @staticmethod
    cdef _py_TetVesicleVesRaft from_ptr(TetVesicleVesRaft *ptr):
        cdef _py_TetVesicleVesRaft obj = _py_TetVesicleVesRaft.__new__(_py_TetVesicleVesRaft )
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_TetVesicleVesRaft from_ref(const TetVesicleVesRaft &ref):
        _py_TetVesicleVesRaft.from_ptr(<TetVesicleVesRaft*>&ref)
