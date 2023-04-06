# cython: language_level=2
###___license_placeholder___###

include "cysteps.pyx"

# ======================================================================================================================
# Python bindings to namespace steps::mpi
# ======================================================================================================================
cimport steps_mpi
from steps_mpi cimport TetOpSplitP

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

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TetOpSplitP(_py_TetAPI):
    """Bindings for MPI TetOpSplitP"""
# ----------------------------------------------------------------------------------------------------------------------
    cdef TetOpSplitP *ptrx(self):
        return <TetOpSplitP*> self._ptr

    def __init__(self, _py_Model model, _py_Geom geom, _py_RNG rng, int calcMembPot=0, std.vector[uint] tet_hosts = [], dict tri_hosts = {}, std.vector[uint] wm_hosts = []):
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
        cdef std.map[steps.triangle_id_t, uint] _tri_hosts
        for key, elem in tri_hosts.items():
            _tri_hosts[steps.triangle_id_t(key)] = elem
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
        self.ptrx().checkpoint(to_std_string(file_name))

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
        self.ptrx().restore(to_std_string(file_name))

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


    def getBatchTetCounts(self, std.vector[index_t] tets, str spec):
        """
        Get the counts of a species spec in a list of tetrahedrons.

        Syntax::

            getBatchTetCounts(tets, spec)

        Arguments:
        list<index_t> tets
        string spec

        Return:
        list<double>

        """
        return self.ptrx().getBatchTetCounts(tets, to_std_string(spec))

    def getBatchTriCounts(self, std.vector[index_t] tris, str spec):
        """
        Get the counts of a species spec in a list of triangles.

        Syntax::

            getBatchTriCounts(tris, spec)

        Arguments:
        list<index_t> tris
        string spec

        Return:
        list<double>

        """
        return self.ptrx().getBatchTriCounts(tris, to_std_string(spec))

    def setBatchTetConcs(self, std.vector[index_t] tets, str spec, std.vector[double] concs):
        """
        Set the concentration of a species spec in a list of tetrahedrons individually.

        Syntax::

            setBatchTetConcs(tets, spec, concs)

        Arguments:
        list<index_t> tets
        string spec
        list<double> concs

        Return:
        None

        """
        self.ptrx().setBatchTetConcs(tets, to_std_string(spec), concs)

    def getBatchTetConcs(self, std.vector[index_t] tets, str spec):
        """
        Get the individual concentration of a species spec in a list of tetrahedrons.

        Syntax::

            getBatchTetConcs(tets, spec)

        Arguments:
        list<index_t> tets
        string spec

        Return:
        list<double>

        """
        return self.ptrx().getBatchTetConcs(tets, to_std_string(spec))

    # ---------------------------------------------------------------------------------
    # NUMPY section - we accept numpy arrays and generically typed memory-views
    # ---------------------------------------------------------------------------------
    def getBatchTetCountsNP(self, index_t[:] index_array, str spec, double[:] counts):
        """
        Get the counts of a species spec in a list of tetrahedrons.

        Syntax::
            getBatchTetCountsNP(indices, spec, counts)

        Arguments:
        numpy.array<index_t> indices
        string spec
        numpy.array<double, length = len(indices)> counts

        Return:
        None

        """
        self.ptrx().getBatchTetCountsNP(&index_array[0], index_array.shape[0], to_std_string(spec), &counts[0], counts.shape[0])

    def getBatchTetConcsNP(self, index_t[:] index_array, str spec, double[:] concs):
        """
        Get the individual concentrations of a species spec in a list of tetrahedrons.

        Syntax::
            getBatchTetConcsNP(indices, spec, concs)

        Arguments:
        numpy.array<index_t> indices
        string spec
        numpy.array<double, length = len(indices)> concs

        Return:
        None

        """
        self.ptrx().getBatchTetConcsNP(&index_array[0], index_array.shape[0], to_std_string(spec), &concs[0], concs.shape[0])

    def getBatchTriCountsNP(self, index_t[:] index_array, str spec, double[:] counts):
        """
        Get the counts of a species spec in a list of triangles.

        Syntax::
            getBatchTriCountsNP(indices, spec, counts)

        Arguments:
        numpy.array<index_t> indices
        string spec
        numpy.array<double, length = len(indices)> counts

        Return:
            None

        """
        self.ptrx().getBatchTriCountsNP(&index_array[0], index_array.shape[0], to_std_string(spec), &counts[0], counts.shape[0])

    def setBatchTetConcsNP(self, index_t[:] index_array, str spec, double[:] concs):
        """
        Set the concetration of a species spec in a list of tetrahedrons.

        Syntax::
            setBatchTetConcsNP(indices, spec, concs)

        Arguments:
        numpy.array<index_t> indices
        string spec
        numpy.array<double, length = len(indices)> concs

        Return:
        None

        """
        self.ptrx().setBatchTetConcsNP(&index_array[0], index_array.shape[0], to_std_string(spec), &concs[0], concs.shape[0])


    def getBatchTetConcsNP(self, index_t[:] index_array, str spec, double[:] concs):
        """
        Get the concetration of a species spec in a list of tetrahedrons.

        Syntax::
            getBatchTetConcsNP(indices, spec, concs)

        Arguments:
        numpy.array<index_t> indices
        string spec
        numpy.array<double, length = len(indices)>

        Return:
        None

        """
        self.ptrx().getBatchTetConcsNP(&index_array[0], index_array.shape[0], to_std_string(spec), &concs[0], concs.shape[0])

    def sumBatchTetCountsNP(self, index_t[:] tet_array, str spec):
        """
        Return the accumulated sum of species spec in a batch of tetrahedrons.

        This function requires NumPy array as input, and called globally in all processes.

        Syntax::

            sumBatchTetCountsNP(tet_array, spec)

        Arguments:
        numpy.array<index_t> tet_array
        string spec

        Return:
        float
        """
        return self.ptrx().sumBatchTetCountsNP(&tet_array[0], tet_array.shape[0], to_std_string(spec))

    def sumBatchTriCountsNP(self, index_t[:] tri_array, str spec):
        """
        Return the accumulated sum of species spec in a batch of triangles.

        This function requires NumPy array as input, and called globally in all processes.

        Syntax::

            sumBatchTriCountsNP(tri_array, spec)

        Arguments:
        numpy.array<index_t> tri_array
        string spec

        Return:
        float
        """
        return self.ptrx().sumBatchTriCountsNP(&tri_array[0], tri_array.shape[0], to_std_string(spec))

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

    def getROITetCounts(self, str roi, str spec):
        """
        Get the counts of a species spec in tetrehedrons of a ROI.

        Syntax::

            getROITetCounts(roi, spec)

        Arguments:
        string roi
        string spec

        Return:
        list<float>

        """
        return self.ptrx().getROITetCounts(to_std_string(roi), to_std_string(spec))

    def getROITriCounts(self, str roi, str spec):
        """
        Get the counts of a species spec in triangles of a ROI.

        Syntax::

            getROITriCounts(roi, spec)

        Arguments:
        string roi
        string spec

        Return:
        list<float>

        """
        return self.ptrx().getROITriCounts(to_std_string(roi), to_std_string(spec))

    def getROITetCountsNP(self, str roi, str spec, double[:] counts):
        """
        Get the counts of a species spec in tetrehedrons of a ROI.

        Syntax::
            getROITetCountsNP(roi, spec, counts)

        Arguments:
        string roi
        string spec
        numpy.array<float, length = len(indices)>

        Return:
            None

        """
        self.ptrx().getROITetCountsNP(to_std_string(roi), to_std_string(spec), &counts[0], counts.shape[0])

    def getROITriCountsNP(self, str roi, str spec, double[:] counts):
        """
        Get the counts of a species spec in triangles of a ROI.

        Syntax::
            getROITriCountsNP(roi, spec, counts)

        Arguments:
        string roi
        string spec
        numpy.array<float, length = len(indices)>

        Return:
            None

        """
        self.ptrx().getROITriCountsNP(to_std_string(roi), to_std_string(spec), &counts[0], counts.shape[0])

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

    def getROICount(self, str roi, str spec):
        """
        Get the count of a species in a ROI.

        Syntax::
            getROICount(roi, spec)

        Arguments:
        string roi
        string spec

        Return:
        float

        """
        return self.ptrx().getROICount(to_std_string(roi), to_std_string(spec))

    def setROICount(self, str roi, str spec, double count):
        """
        Set the count of a species in a ROI.

        Syntax::
            setROICount(roi, spec, count)

        Arguments:
        string roi
        string spec
        float count

        Return:
        None

        """
        self.ptrx().setROICount(to_std_string(roi), to_std_string(spec), count)

    def getROIAmount(self, str roi, str spec):
        """
        Get the amount of a species in a ROI.

        Syntax::
            getROIAmount(roi, spec, count)

        Arguments:
        string roi
        string spec

        Return:
        float

        """
        return self.ptrx().getROIAmount(to_std_string(roi), to_std_string(spec))

    def setROIAmount(self, str roi, str spec, double amount):
        """
        Set the amount of a species in a ROI.

        Syntax::
            setROIAmount(roi, spec, amount)

        Arguments:
        string roi
        string spec
        float amount

        Return:
        None

        """
        return self.ptrx().setROIAmount(to_std_string(roi), to_std_string(spec), amount)

    def getROIConc(self, str roi, str spec):
        """
        Get the concentration of a species in a ROI.

        Syntax::
            getROIConc(roi, spec)

        Arguments:
        string roi
        string spec

        Return:
        float

        """
        return self.ptrx().getROIConc(to_std_string(roi), to_std_string(spec))

    def setROIConc(self, str roi, str spec, double conc):
        """
        Set the concentration of a species in a ROI.

        Syntax::
            setROIConc(roi, spec, conc)

        Arguments:
        string roi
        string spec
        float conc

        Return:
        None

        """
        self.ptrx().setROIConc(to_std_string(roi), to_std_string(spec), conc)

    def setROIClamped(self, str roi, str spec, bool clamped):
        """
        Set a species in a ROI to be clamped or not. The count of species spec in the ROI is clamped if
        clamped is True, not clamped if clamped is False.

        Syntax::
            setROIClamped(roi, spec, clamped)

        Arguments:
        string roi
        string spec
        bool clamped

        Return:
        None

        """
        self.ptrx().setROIClamped(to_std_string(roi), to_std_string(spec), clamped)

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

    def setROIDiffActive(self, str roi, str diff, bool act):
        """
        Set diffusion diff in a ROI to be active or not.

        Syntax::
            setROIDiffActive(roi, diff, a)

        Arguments:
        string roi
        string diff
        bool a

        Return:
        None

        """
        self.ptrx().setROIDiffActive(to_std_string(roi), to_std_string(diff), act)

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
        Return the extent of reaction with identifier string r in ROI with
        identifier string roi, that is the number of times the reaction has occurred up
        to the current simulation time.

        Syntax::
            getROIReacExtent(roi, r)

        Arguments:
        string roi
        string r

        Return:
        index_t

        """
        return self.ptrx().getROIReacExtent(to_std_string(roi), to_std_string(r))

    def resetROIReacExtent(self, str roi, str r):
        """
        Reset the extent of reaction with identifier string r in ROI with
        identifier string roi, that is the number of times the reaction has occurred up
        to the current simulation time, to 0.

        Syntax::
            resetROIReacExtent(roi, r)

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
        Reset the extent of surface reaction with identifier string r in ROI with
        identifier string roi, that is the number of times the reaction has occurred up
        to the current simulation time, to 0.

        Syntax::
            resetROISReacExtent(roi, r)

        Arguments:
        string roi
        string sr

        Return:
        None

        """
        self.ptrx().resetROISReacExtent(to_std_string(roi), to_std_string(sr))

    def getROIDiffExtent(self, str roi, str diff):
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
        return self.ptrx().getROIDiffExtent(to_std_string(roi), to_std_string(diff))

    def resetROIDiffExtent(self, str roi, str spec):
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
        self.ptrx().resetROIDiffExtent(to_std_string(roi), to_std_string(spec))


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

    def repartitionAndReset(self, std.vector[uint] tet_hosts=(), dict tri_hosts=None, std.vector[uint] wm_hosts=()):
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
        cdef std.map[uint, uint] _tri_hosts = tri_hosts
        self.ptrx().repartitionAndReset(tet_hosts, _tri_hosts, wm_hosts)


    @staticmethod
    cdef _py_TetOpSplitP from_ptr(TetOpSplitP *ptr):
        cdef _py_TetOpSplitP obj = _py_TetOpSplitP.__new__(_py_TetOpSplitP )
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_TetOpSplitP from_ref(const TetOpSplitP &ref):
        _py_TetOpSplitP.from_ptr(<TetOpSplitP*>&ref)
