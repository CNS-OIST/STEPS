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
    Finalsie the MPI solver. NOTE: handled automatically, should not be called by user.
    """
    steps_mpi.mpiFinish()

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TetOpSplitP(_py_API):
    """Bindings for MPI TetOpSplitP"""
# ----------------------------------------------------------------------------------------------------------------------
    cdef TetOpSplitP *ptrx(self):
        return <TetOpSplitP*> self._ptr

    def __init__(self, _py_Model model, _py_Geom geom, _py_RNG rng, int calcMembPot=0, std.vector[uint] tet_hosts = [], dict tri_hosts = {}, std.vector[uint] wm_hosts = []):
        """        
        Construction::
        
            sim = steps.solver.TetOpSplit(model, geom, rng, tet_hosts=[], tri_hosts={}, wm_hosts=[], calcMembPot=0)
        
        Create a spatial stochastic solver based on operator splitting, that is that reaction events are partitioned and diffusion is approximated. 
        If voltage is to be simulated, argument calcMembPot specifies the solver e.g. calcMembPot=steps.solver.EF_DV_PETSC will utilise the PETSc library. calcMembPot=0 means voltage will not be simulated. 
        
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
        Advance the simulation for secs seconds. 

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
        if not isinstance(file_name, bytes):
            file_name = file_name.encode()

        self.ptrx().checkpoint(file_name)

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
        if not isinstance(file_name, bytes):
            file_name = file_name.encode()

        self.ptrx().restore(file_name)

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
        if not isinstance(opt_file_name, bytes):
            opt_file_name = opt_file_name.encode()

        self.ptrx().saveMembOpt(opt_file_name)

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


    def getBatchTetCounts(self, std.vector[index_t] tets, str s):
        """
        Get the counts of a species s in a list of tetrahedrons.

        Syntax::

            getBatchTetCounts(tets, s)

        Arguments:
        list<index_t> tets
        string s

        Return:
        list<double>

        """
        if not isinstance(s, bytes):
            s = s.encode()

        return self.ptrx().getBatchTetCounts(tets, s)

    def getBatchTriCounts(self, std.vector[index_t] tris, str s):
        """
        Get the counts of a species s in a list of triangles.

        Syntax::

            getBatchTriCounts(tris, s)

        Arguments:
        list<index_t> tris
        string s

        Return:
        list<double>

        """
        if not isinstance(s, bytes):
            s = s.encode()

        return self.ptrx().getBatchTriCounts(tris, s)

    def setBatchTetConcs(self, std.vector[index_t] tets, str s, std.vector[double] concs):
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
        t = s
        if not isinstance(s, bytes):
            t = s.encode()

        self.ptrx().setBatchTetConcs(tets, t, concs)

    def getBatchTetConcs(self, std.vector[index_t] tets, str s):
        """
        Get the individual concentration of a species s in a list of tetrahedrons.

        Syntax::

            getBatchTetConcs(tets, s)

        Arguments:
        list<index_t> tets
        string s

        Return:
        list<double>

        """
        t = s
        if not isinstance(s, bytes):
            t = s.encode()

        return self.ptrx().getBatchTetConcs(tets, t)

    # ---------------------------------------------------------------------------------
    # NUMPY section - we accept numpy arrays and generically typed memory-views
    # ---------------------------------------------------------------------------------
    def getBatchTetCountsNP(self, index_t[:] index_array, str s, double[:] counts):
        """
        Get the counts of a species s in a list of tetrahedrons.

        Syntax::
            getBatchTetCountsNP(indices, s, counts)

        Arguments:
        numpy.array<index_t> indices
        string s
        numpy.array<double, length = len(indices)>

        Return:
        None

        """
        if not isinstance(s, bytes):
            s = s.encode()

        self.ptrx().getBatchTetCountsNP(&index_array[0], index_array.shape[0], s, &counts[0], counts.shape[0])

    def getBatchTriCountsNP(self, index_t[:] index_array, str s, double[:] counts):
        """
        Get the counts of a species s in a list of triangles.

        Syntax::
            getBatchTriCountsNP(indices, s, counts)

        Arguments:
        numpy.array<index_t> indices
        string s
        numpy.array<double, length = len(indices)>

        Return:
            None

        """
        if not isinstance(s, bytes):
            s = s.encode()

        self.ptrx().getBatchTriCountsNP(&index_array[0], index_array.shape[0], s, &counts[0], counts.shape[0])

    def sumBatchTetCountsNP(self, uint[:] tet_array, str s):
        """
        Return the accumulated sum of species s in a batch of tetrahedrons.
        
        This function requires NumPy array as input, and called globally in all processes.
        
        Syntax::
            
            sumBatchTetCountsNP(tet_list, s)
            
        Arguments:
        numpy.array<uint> tet_list
        string s
        
        Return:
        float
        """
        if not isinstance(s, bytes):
            s = s.encode()

        return self.ptrx().sumBatchTetCountsNP(&tet_array[0], tet_array.shape[0], s)

    def sumBatchTriCountsNP(self, uint[:] tri_array, str s):
        """
        Return the accumulated sum of species s in a batch of triangles.
        
        This function requires NumPy array as input, and called globally in all processes.
        
        Syntax::
            
            sumBatchTriCountsNP(tri_list, s)
            
        Arguments:
        numpy.array<uint> tri_list
        string s
        
        Return:
        float
        """
        if not isinstance(s, bytes):
            s = s.encode()

        return self.ptrx().sumBatchTriCountsNP(&tri_array[0], tri_array.shape[0], s)

    def sumBatchTriGHKIsNP(self, uint[:] tri_array, str ghk):
        """
        Return the accumulated sum of GHK currents in a batch of triangles.
        
        This function requires NumPy array as input, and called globally in all processes.
        
        Syntax::
            
            sumBatchTriGHKIsNP(tri_list, ghk)
            
        Arguments:
        numpy.array<uint> tri_list
        string ghk
        
        Return:
        float
        """
        if not isinstance(ghk, bytes):
            ghk = ghk.encode()

        return self.ptrx().sumBatchTriGHKIsNP(&tri_array[0], tri_array.shape[0], ghk)

    def sumBatchTriOhmicIsNP(self, uint[:] tri_array, str ghk):
        """
        Return the accumulated sum of Ohmic currents in a batch of triangles.
        
        This function requires NumPy array as input, and called globally in all processes.
        
        Syntax::
            
            sumBatchTriOhmicIsNP(tri_list, oc)
            
        Arguments:
        numpy.array<uint> tri_list
        string oc
        
        Return:
        float
        """
        if not isinstance(ghk, bytes):
            ghk = ghk.encode()

        return self.ptrx().sumBatchTriOhmicIsNP(&tri_array[0], tri_array.shape[0], ghk)

    # ---------------------------------------------------------------------------------
    # ROI section
    # ---------------------------------------------------------------------------------

    def getROITetCounts(self, str ROI_id, str s):
        """
        Get the counts of a species s in tetrehedrons of a ROI.

        Syntax::

            getROITetCounts(ROI_id, s)

        Arguments:
        string ROI_id
        string s

        Return:
        list<float>

        """
        return self.ptrx().getROITetCounts(to_std_string(ROI_id), to_std_string(s))

    def getROITriCounts(self, str ROI_id, str s):
        """
        Get the counts of a species s in triangles of a ROI.

        Syntax::

            getROITriCounts(ROI_id, s)

        Arguments:
        string ROI_id
        string s

        Return:
        list<float>

        """
        return self.ptrx().getROITriCounts(to_std_string(ROI_id), to_std_string(s))

    def getROITetCountsNP(self, str ROI_id, str s, double[:] counts):
        """
        Get the counts of a species s in tetrehedrons of a ROI.

        Syntax::
            getROITetCountsNP(ROI_id, s, counts)

        Arguments:
        string ROI_id
        string s
        numpy.array<float, length = len(indices)>

        Return:
            None

        """
        self.ptrx().getROITetCountsNP(to_std_string(ROI_id), to_std_string(s), &counts[0], counts.shape[0])

    def getROITriCountsNP(self, str ROI_id, str s, double[:] counts):
        """
        Get the counts of a species s in triangles of a ROI.

        Syntax::
            getROITriCountsNP(ROI_id, s, counts)

        Arguments:
        string ROI_id
        string s
        numpy.array<float, length = len(indices)>

        Return:
            None

        """
        self.ptrx().getROITriCountsNP(to_std_string(ROI_id), to_std_string(s), &counts[0], counts.shape[0])

    def getROIVol(self, str ROI_id):
        """
        Get the volume of a ROI.

        Syntax::
            getROIVol(ROI_id)

        Arguments:
        string ROI_id

        Return:
        float

        """
        return self.ptrx().getROIVol(to_std_string(ROI_id))

    def getROIArea(self, str ROI_id):
        """
        Get the area of a ROI.

        Syntax::
            getROIArea(ROI_id)

        Arguments:
        string ROI_id

        Return:
        float

        """
        return self.ptrx().getROIArea(to_std_string(ROI_id))

    def getROICount(self, str ROI_id, str s):
        """
        Get the count of a species in a ROI.

        Syntax::
            getROICount(ROI_id, s)

        Arguments:
        string ROI_id
        string s

        Return:
        float

        """
        return self.ptrx().getROICount(to_std_string(ROI_id), to_std_string(s))

    def setROICount(self, str ROI_id, str s, double count):
        """
        Set the count of a species in a ROI.

        Syntax::
            setROICount(ROI_id, s, count)

        Arguments:
        string ROI_id
        string s
        float count

        Return:
        None

        """
        self.ptrx().setROICount(to_std_string(ROI_id), to_std_string(s), count)

    def getROIAmount(self, str ROI_id, str s):
        """
        Get the amount of a species in a ROI.

        Syntax::
            getROIAmount(ROI_id, s, count)

        Arguments:
        string ROI_id
        string s

        Return:
        float

        """
        return self.ptrx().getROIAmount(to_std_string(ROI_id), to_std_string(s))

    def getROIConc(self, str ROI_id, str s):
        """
        Get the concentration of a species in a ROI.

        Syntax::
            getROIConc(ROI_id, s)

        Arguments:
        string ROI_id
        string s

        Return:
        float

        """
        return self.ptrx().getROIConc(to_std_string(ROI_id), to_std_string(s))

    def setROIConc(self, str ROI_id, str s, double conc):
        """
        Set the concentration of a species in a ROI.

        Syntax::
            setROIConc(ROI_id, s, conc)

        Arguments:
        string ROI_id
        string s
        float conc

        Return:
        None

        """   
        self.ptrx().setROIConc(to_std_string(ROI_id), to_std_string(s), conc)

    def setROIClamped(self, str ROI_id, str s, bool b):
        """
        Set a species in a ROI to be clamped or not. The count of species s in the ROI is clamped if
        b is True, not clamped if b is False.

        Syntax::
            setROIClamped(ROI_id, s, b)

        Arguments:
        string ROI_id
        string s
        bool b

        Return:
        None

        """
        self.ptrx().setROIClamped(to_std_string(ROI_id), to_std_string(s), b)

    def setROIReacK(self, str ROI_id, str r, double kf):
        """
        Sets the macroscopic reaction constant of reaction with identifier string r
        in a ROI with identifier string ROI_id to kf. The unit of the reaction constant
        depends on the order of the reaction.

        Note: The default value still comes from the steps.model description, so 
        calling reset() will return the reaction constant to that value.

        Syntax::
            setROIReacK(ROI_id, r, kf)

        Arguments:
        string ROI_id
        string r
        float kf

        Return:
        None

        """
        self.ptrx().setROIReacK(to_std_string(ROI_id), to_std_string(r), kf)

    def setROISReacK(self, str ROI_id, str sr, double kf):
        """
        Sets the macroscopic reaction constant of surface reaction with identifier string sr
        in a ROI with identifier string ROI_id to kf. The unit of the reaction constant
        depends on the order of the reaction.

        Note: The default value still comes from the steps.model description, so 
        calling reset() will return the reaction constant to that value.

        Syntax::
            setROISReacK(ROI_id, sr, kf)

        Arguments:
        string ROI_id
        string sr
        float kf

        Return:
        None

        """
        self.ptrx().setROISReacK(to_std_string(ROI_id), to_std_string(sr), kf)

    def setROIDiffD(self, str ROI_id, str d, double dk):
        """
        Sets the macroscopic diffusion constant of diffusion with identifier string d
        in a ROI with identifier string ROI_id to dk.

        Note: The default value still comes from the steps.model description, so 
        calling reset() will return the diffusion constant to that value.

        Syntax::
            setROIDiffD(ROI_id, d, dk)

        Arguments:
        string ROI_id
        string d
        float dk

        Return:
            None

        """
        self.ptrx().setROIDiffD(to_std_string(ROI_id), to_std_string(d), dk)

    def setROIReacActive(self, str ROI_id, str r, bool a):
        """
        Set reaction r in a ROI to be active or not.

        Syntax::
            setROIReacActive(ROI_id, r, a)

        Arguments:
        string ROI_id
        string r
        bool a

        Return:
        None

        """
        self.ptrx().setROIReacActive(to_std_string(ROI_id), to_std_string(r), a)

    def setROISReacActive(self, str ROI_id, str sr, bool a):
        """
        Set surface reaction sr in a ROI to be active or not.

        Syntax::
            setROISReacActive(ROI_id, sr, a)

        Arguments:
        string ROI_id
        string sr
        bool a

        Return:
        None

        """
        self.ptrx().setROISReacActive(to_std_string(ROI_id), to_std_string(sr), a)

    def setROIDiffActive(self, str ROI_id, str d, bool act):
        """
        Set diffusion d in a ROI to be active or not.

        Syntax::
            setROIDiffActive(ROI_id, sr, a)

        Arguments:
        string ROI_id
        string sr
        bool a

        Return:
        None

        """
        self.ptrx().setROIDiffActive(to_std_string(ROI_id), to_std_string(d), act)

    def setROIVDepSReacActive(self, str ROI_id, str vsr, bool a):
        """
        Set voltage dependent surface reaction vsr in a ROI to be active or not.

        Syntax::
            setROIVDepSReacActive(ROI_id, vsr, a)

        Arguments:
        string ROI_id
        string vsr
        bool a

        Return:
        None

        """
        self.ptrx().setROIVDepSReacActive(to_std_string(ROI_id), to_std_string(vsr), a)

    def getROIReacExtent(self, str ROI_id, str r):
        """
        Return the extent of reaction with identifier string r in ROI with
        identifier string ROI_id, that is the number of times the reaction has occurred up
        to the current simulation time.

        Syntax::
            getROIReacExtent(ROI_id, r)

        Arguments:
        string ROI_id
        string r

        Return:
        index_t

        """
        return self.ptrx().getROIReacExtent(to_std_string(ROI_id), to_std_string(r))

    def resetROIReacExtent(self, str ROI_id, str r):
        """
        Reset the extent of reaction with identifier string r in ROI with
        identifier string ROI_id, that is the number of times the reaction has occurred up
        to the current simulation time, to 0.

        Syntax::
            resetROIReacExtent(ROI_id, r)

        Arguments:
        string ROI_id
        string r

        Return:
        None

        """
        self.ptrx().resetROIReacExtent(to_std_string(ROI_id), to_std_string(r))

    def getROISReacExtent(self, str ROI_id, str sr):
        """
        Return the extent of surface reaction with identifier string sr in ROI with
        identifier string ROI_id, that is the number of times the reaction has occurred up
        to the current simulation time.

        Syntax::
            getROISReacExtent(ROI_id, sr)

        Arguments:
        string ROI_id
        string sr

        Return:
        index_t

        """
        return self.ptrx().getROISReacExtent(to_std_string(ROI_id), to_std_string(sr))

    def resetROISReacExtent(self, str ROI_id, str sr):
        """
        Reset the extent of surface reaction with identifier string r in ROI with
        identifier string ROI_id, that is the number of times the reaction has occurred up
        to the current simulation time, to 0.

        Syntax::
            resetROISReacExtent(ROI_id, r)

        Arguments:
        string ROI_id
        string sr

        Return:
        None

        """
        self.ptrx().resetROISReacExtent(to_std_string(ROI_id), to_std_string(sr))

    def getROIDiffExtent(self, str ROI_id, str d):
        """
        Return the extent of diffusion with identifier string d in ROI with
        identifier string ROI_id, that is the number of times the diffusion has occurred up
        to the current simulation time.

        Syntax::
            getROIDiffExtent(ROI_id, d)

        Arguments:
        string ROI_id
        string d

        Return:
        index_t

        """
        return self.ptrx().getROIDiffExtent(to_std_string(ROI_id), to_std_string(d))

    def resetROIDiffExtent(self, str ROI_id, str s):
        """
        Reset the extent of diffusion with identifier string d in ROI with
        identifier string ROI_id, that is the number of times the diffusion has occurred up
        to the current simulation time, to 0.

        Syntax::
            resetROIDiffExtent(ROI_id, d)

        Arguments:
        string ROI_id
        string d

        Return:
        None

        """
        self.ptrx().resetROIDiffExtent(to_std_string(ROI_id), to_std_string(s))


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
