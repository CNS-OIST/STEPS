###___license_placeholder___###

from steps_wmrk4 cimport *
from steps_wmdirect cimport *
from steps_wmrssa cimport *
from steps_tetexact cimport *
from steps_tetode cimport *
from steps_solver cimport *
from steps cimport index_t

# ======================================================================================================================
# Python bindings to namespace steps::wmrk4
# ======================================================================================================================

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Wmrk4(_py_API):
    "Python wrapper class for Wmrk4"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Wmrk4] _autodealoc
    cdef Wmrk4 *ptrx(self):
        return <Wmrk4*> self._ptr

    def __init__(self, _py_Model m, _py_Geom g, _py_RNG r = None):
        """        
        Construction::
        
            sim = steps.solver.Wmrk4(model, geom)
        
        Create a non-spatial deterministic solver based on the Runge-Kutta fourth order method.
        
        Arguments:
        steps.model.Model model
        steps.geom.Geom geom
        """

        if m == None:
            raise TypeError('The Model object is empty.')
        if g == None:
            raise TypeError('The Geom object is empty.')

        self._ptr = new Wmrk4(m.ptr(), g.ptr(), r.ptr() if r else shared_ptr[RNG]())
        _py_API.__init__(self, m, g, r)

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
        Typically, this resets all concentrations of all chemical species in 
        all elements (whether compartments and patches in a well-mixed solver 
        or tetrahedrons and triangles in a mesh-based solver) to zero, 
        resets the simulation time to zero and resets reaction (and diffusion) 
        rates to the default values described in the steps.model objects. 
        All reaction (and diffusion) rules are reset to active and all 
        compartment volumes and patch areas are reset to default values 
        described in steps.geom objects (for well-mixed solvers). 
        Usually, this method should be called before starting each simulation iteration.

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

    def setDT(self, double dt):
        """
        Set the stepsize for numerical solvers. Superceded by setRk4DT, but
        Included for backwards compatability.
         
        Syntax::
         
            setDT(dt)
         
        Arguments:
        float dt
         
        Return:
        None

        """
        self.ptrx().setDT(dt)

    def setRk4DT(self, double dt):
        """
        Set the stepsize for numerical solvers. Must be called before running a 
        simulation with these solvers (currently Wmrk4) since there is no default 
        stepsize. The deterministic solver Wmrk4 implements a fixed stepsize 
        (i.e. not adaptive), although the stepsize can be altered at any point 
        during the simulation with this method.

        Syntax::
            
            setRk4tDT(dt)
            
        Arguments:
        float dt

        Return:
        None

        """
        self.ptrx().setRk4DT(dt)

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

    @staticmethod
    cdef _py_Wmrk4 from_ptr(Wmrk4 *ptr):
        cdef _py_Wmrk4 obj = _py_Wmrk4.__new__(_py_Wmrk4)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Wmrk4 from_ref(const Wmrk4 &ref):
        return _py_Wmrk4.from_ptr(<Wmrk4*>&ref)


# ======================================================================================================================
# Python bindings to namespace steps::wmdirect
# ======================================================================================================================

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Wmdirect(_py_API):
    "Python wrapper class for Wmdirect"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Wmdirect] _autodealoc
    cdef Wmdirect *ptrd(self):
        return <Wmdirect*> self._ptr

    def __init__(self, _py_Model m, _py_Geom g, _py_RNG r):
        """        
        Construction::
        
            sim = steps.solver.Wmdirect(model, geom, rng)
        
        Create a non-spatial stochastic solver based on Gillespie's SSA.
        
        Arguments:
        steps.model.Model model
        steps.geom.Geom geom
        steps.rng.RNG rng
        """

        if m == None:
            raise TypeError('The Model object is empty.')
        if g == None:
            raise TypeError('The Geom object is empty.')
        if r == None:
            raise TypeError('The RNG object is empty.')

        self._ptr = new Wmdirect( m.ptr(), g.ptr(), r.ptr() )
        #super(self.__class__, self).__init__(m,g,r)
        _py_API.__init__(self, m, g, r)

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
        self.ptrd().checkpoint(to_std_string(file_name))

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
        self.ptrd().restore(to_std_string(file_name))

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
        return self.ptrd().getSolverName()

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
        return self.ptrd().getSolverDesc()

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
        return self.ptrd().getSolverAuthors()

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
        return self.ptrd().getSolverEmail()

    def reset(self, ):
        """
        Reset the simulation to the state the solver was initialised to. 
        Typically, this resets all concentrations of all chemical species in 
        all elements (whether compartments and patches in a well-mixed solver 
        or tetrahedrons and triangles in a mesh-based solver) to zero, 
        resets the simulation time to zero and resets reaction (and diffusion) 
        rates to the default values described in the steps.model objects. 
        All reaction (and diffusion) rules are reset to active and all 
        compartment volumes and patch areas are reset to default values 
        described in steps.geom objects (for well-mixed solvers). 
        Usually, this method should be called before starting each simulation iteration.

        Syntax::
            
            reset()
            
        Arguments:
        None

        Return:
        None

        """
        self.ptrd().reset()

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
        self.ptrd().run(endtime)

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
        self.ptrd().advance(adv)

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
        self.ptrd().step()

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
        return self.ptrd().getTime()

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
        return self.ptrd().getA0()

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
        return self.ptrd().getNSteps()

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
        self.ptrd().setTime(time)

    def setNSteps(self, uint nsteps):
        """
        Set the number of 'realizations' of the SSA, the number of reaction 
        (and diffusion) events in stochastic solvers.

        Syntax::
            
            setNSteps(nsteps)
            
        Arguments:
        int nsteps

        Return:
        None

        """
        self.ptrd().setNSteps(nsteps)

    # def addKProc(self, steps.wmdirect.KProc* kp):
    #     return _py_void.from_ref(self.ptr().addKProc(kp.ptr()))


    @staticmethod
    cdef _py_Wmdirect from_ptr(Wmdirect *ptr):
        cdef _py_Wmdirect obj = _py_Wmdirect.__new__(_py_Wmdirect)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Wmdirect from_ref(const Wmdirect &ref):
        return _py_Wmdirect.from_ptr(<Wmdirect*>&ref)


# ======================================================================================================================
# Python bindings to namespace steps::wmrssa
# ======================================================================================================================

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Wmrssa(_py_API):
    "Python wrapper class for Wmrssa"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Wmrssa] _autodealoc
    cdef Wmrssa *ptrd(self):
        return <Wmrssa*> self._ptr

    def __init__(self, _py_Model m, _py_Geom g, _py_RNG r):
        """        
        Construction::
        
            sim = steps.solver.Wmrssa(model, geom, rng)
        
        Create a non-spatial stochastic solver based on Gillespie's SSA.
        
        Arguments:
        steps.model.Model model
        steps.geom.Geom geom
        steps.rng.RNG rng
        """
        if m == None:
            raise TypeError('The Model object is empty.')
        if g == None:
            raise TypeError('The Geom object is empty.')
        if r == None:
            raise TypeError('The RNG object is empty.')
        self._ptr = new Wmrssa( m.ptr(), g.ptr(), r.ptr() )
        #super(self.__class__, self).__init__(m,g,r)
        _py_API.__init__(self, m, g, r)

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
        self.ptrd().checkpoint(to_std_string(file_name))

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
        self.ptrd().restore(to_std_string(file_name))

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
        return self.ptrd().getSolverName()

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
        return self.ptrd().getSolverDesc()

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
        return self.ptrd().getSolverAuthors()

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
        return self.ptrd().getSolverEmail()

    def reset(self, ):
        """
        Reset the simulation to the state the solver was initialised to. 
        Typically, this resets all concentrations of all chemical species in 
        all elements (whether compartments and patches in a well-mixed solver 
        or tetrahedrons and triangles in a mesh-based solver) to zero, 
        resets the simulation time to zero and resets reaction (and diffusion) 
        rates to the default values described in the steps.model objects. 
        All reaction (and diffusion) rules are reset to active and all 
        compartment volumes and patch areas are reset to default values 
        described in steps.geom objects (for well-mixed solvers). 
        Usually, this method should be called before starting each simulation iteration.

        Syntax::
            
            reset()
            
        Arguments:
        None

        Return:
        None

        """
        self.ptrd().reset()

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
        self.ptrd().run(endtime)

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
        self.ptrd().advance(adv)

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
        self.ptrd().step()

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
        return self.ptrd().getTime()

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
        return self.ptrd().getNSteps()

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
        self.ptrd().setTime(time)

    def setNSteps(self, uint nsteps):
        """
        Set the number of 'realizations' of the SSA, the number of reaction 
        (and diffusion) events in stochastic solvers.

        Syntax::
            
            setNSteps(nsteps)
            
        Arguments:
        int nsteps

        Return:
        None

        """
        self.ptrd().setNSteps(nsteps)

    # def addKProc(self, steps.Wmrssa.KProc* kp):
    #     return _py_void.from_ref(self.ptr().addKProc(kp.ptr()))


    @staticmethod
    cdef _py_Wmrssa from_ptr(Wmrssa *ptr):
        cdef _py_Wmrssa obj = _py_Wmrssa.__new__(_py_Wmrssa)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Wmrssa from_ref(const Wmrssa &ref):
        return _py_Wmrssa.from_ptr(<Wmrssa*>&ref)


# ======================================================================================================================
# Python bindings to namespace steps::tetexact
# ======================================================================================================================

# ----------------------------------------------------------------------------------------------------------------------
#cdef class _py_TEKProc(_py__base):
#    "Python wrapper class for TEKProc (steps::tetexact::KProc)"
# ----------------------------------------------------------------------------------------------------------------------
#    cdef TEKProc * ptr(self):
#        return <TEKProc* > self._ptr
#
#    @staticmethod
#    cdef _py_TEKProc from_ptr(TEKProc *ptr):
#        cdef _py_TEKProc obj = _py_TEKProc.__new__(_py_TEKProc)
#        obj._ptr = ptr
#       return obj


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Tetexact(_py_TetAPI):
    "Python wrapper class for Tetexact"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Tetexact] _autodealoc
    cdef Tetexact *ptrx(self):
        return <Tetexact*> self._ptr

    def __init__(self, _py_Model m, _py_Geom g, _py_RNG r, int calcMembPot=0):
        """        
        Construction::
        
            sim = steps.solver.Tetexact(model, geom, rng, calcMembPot = 0)
        
        Create a spatial stochastic solver based on Gillespie's SSA, extended with diffusion across elements in a tetrahedral mesh.
        If voltage is to be simulated, argument calcMemPot=1 will set to the default solver. calcMembPot=0 means voltage will not be simulated. 
        
        Arguments:
        steps.model.Model model
        steps.geom.Geom geom
        steps.rng.RNG rng
        int calcMemPot (default=0)
        
        """
        if m == None:
            raise TypeError('The Model object is empty.')
        if g == None:
            raise TypeError('The Geom object is empty.')
        if r == None:
            raise TypeError('The RNG object is empty.')
        self._ptr = new Tetexact(m.ptr(), g.ptr(), r.ptr(), calcMembPot)
        _py_API.__init__(self, m, g, r)

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
        Typically, this resets all concentrations of all chemical species in 
        all elements (whether compartments and patches in a well-mixed solver 
        or tetrahedrons and triangles in a mesh-based solver) to zero, 
        resets the simulation time to zero and resets reaction (and diffusion) 
        rates to the default values described in the steps.model objects. 
        All reaction (and diffusion) rules are reset to active and all 
        compartment volumes and patch areas are reset to default values 
        described in steps.geom objects (for well-mixed solvers). 
        Usually, this method should be called before starting each simulation iteration.

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
        Saves the vertex optimization in the Efield structure.
                     
        Syntax::
                     
            saveMembOpt()
                     
        Arguments:
        string filename
                     
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
        int nsteps

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
        return self.ptrx().getBatchTetCounts(tets, to_std_string(s))

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
        return self.ptrx().getBatchTriCounts(tris, to_std_string(s))

    def getBatchTetCountsNP(self, index_t[:] indices, str s, double[:] counts):
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
        self.ptrx().getBatchTetCountsNP(&indices[0], indices.shape[0], to_std_string(s), &counts[0], counts.shape[0])

    def getBatchTriCountsNP(self, index_t[:] indices, str s, double[:] counts):
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
        self.ptrx().getBatchTriCountsNP(&indices[0], indices.shape[0], to_std_string(s), &counts[0], counts.shape[0])


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
            getROIAmount(ROI_id, s)

        Arguments:
        string ROI_id
        string s

        Return:
        float

        """
        return self.ptrx().getROIAmount(to_std_string(ROI_id), to_std_string(s))

    def setROIAmount(self, str ROI_id, str s, double amount):
        """
        Set the amount of a species in a ROI.

        Syntax::
            setROIAmount(ROI_id, s, amount)

        Arguments:
        string ROI_id
        string s
        float amount

        Return:
        None

        """
        return self.ptrx().setROIAmount(to_std_string(ROI_id), to_std_string(s), amount)

    def getROIConc(self, str ROI_id, str s):
        """
        Get the concentration of a species in a ROI.

        Syntax::
            getROIConc(ROI_id, s, count)

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


    @staticmethod
    cdef _py_Tetexact from_ptr(Tetexact *ptr):
        cdef _py_Tetexact obj = _py_Tetexact.__new__(_py_Tetexact)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Tetexact from_ref(const Tetexact &ref):
        return _py_Tetexact.from_ptr(<Tetexact*>&ref)


# =======================================================_py_wmrk4===============================================================
# Python bindings to namespace steps::tetode
# ======================================================================================================================

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TetODE(_py_TetAPI):
    "Python wrapper class for TetODE"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[TetODE] _autodealoc
    cdef TetODE *ptrx(self):
        return <TetODE*> self._ptr

    def __init__(self, _py_Model m, _py_Geom g, _py_RNG r=None, int calcMembPot=0):
        """        
        Construction::
        
            sim = steps.solver.TetODE(model, geom, rng=None, calcMembPot = 0)
        
        Create a spatial determinstic solver based on the CVODE library.
        If voltage is to be simulated, argument calcMemPot=1 will set to the default solver. calcMembPot=0 means voltage will not be simulated. 
        
        Arguments:
        steps.model.Model model
        steps.geom.Geom geom
        steps.rng.RNG rng (default=None)
        int calcMemPot (default=0)
        
        """
        if m == None:
            raise TypeError('The Model object is empty.')
        if g == None:
            raise TypeError('The Geom object is empty.')

        self._ptr = new TetODE(m.ptr(), g.ptr(), r.ptr() if r else shared_ptr[RNG](), calcMembPot)
        _py_API.__init__(self, m, g, r)

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
        return self.ptrx().setTemp(t)

    def reset(self, ):
        """
        Reset the simulation to the state the solver was initialised to. 
        Typically, this resets all concentrations of all chemical species in 
        all elements (whether compartments and patches in a well-mixed solver 
        or tetrahedrons and triangles in a mesh-based solver) to zero, 
        resets the simulation time to zero and resets reaction (and diffusion) 
        rates to the default values described in the steps.model objects. 
        All reaction (and diffusion) rules are reset to active and all 
        compartment volumes and patch areas are reset to default values 
        described in steps.geom objects (for well-mixed solvers). 
        Usually, this method should be called before starting each simulation iteration.

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


    def setTolerances(self, double atol, double rtol):
        """
        Set the absolute tolerance and the relative tolerance for CVODE.
                     
        Syntax::
                     
            setTolerance(atol, rtol)
                     
        Arguments:
        float atol
        float rtol
                     
        Return:
        None

        """
        self.ptrx().setTolerances(atol, rtol)

    def setMaxNumSteps(self, uint maxn):
        """
        Sets the maximum number of steps in CVODE per call to run().
        Default is 10000 if this function is not called.
                     
        Syntax::
                     
            setMaxNumSteps()
                     
        Arguments:
        int maxn
                     
        Return:
        None

        """
        self.ptrx().setMaxNumSteps(maxn)


    @staticmethod
    cdef _py_TetODE from_ptr(TetODE *ptr):
        cdef _py_TetODE obj = _py_TetODE.__new__(_py_TetODE)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_TetODE from_ref(const TetODE &ref):
        return _py_TetODE.from_ptr(<TetODE*>&ref)


# ======================================================================================================================
# Python bindings to namespace steps::solver
# ======================================================================================================================
cimport steps_solver

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_API(_py__base):
    "Python wrapper class for API"
# ----------------------------------------------------------------------------------------------------------------------
    cdef _py_Model model
    cdef _py_Geom geom

    cdef API *ptr(self):
        return <API*> self._ptr

    # ---- VIRTUAL - doesnt call original constructor ------
    def __init__(self, _py_Model m, _py_Geom g, _py_RNG r):
        self.model = m
        self.geom = g

    def getCompVol(self, str c):
        """
        Returns the volume of compartment with identifier string comp (in m^3).

        Syntax::
            
            getCompVol(comp)
            
        Arguments:
        string comp

        Return:
        float

        """
        return self.ptr().getCompVol(to_std_string(c))

    def setCompVol(self, str c, double vol):
        """
        Set the volume of compartment with identifier string comp (in m^3).

        Syntax::
            
            setCompVol(comp, vol)
            
        Arguments:
        string comp
        float vol

        Return:
        None

        """
        self.ptr().setCompVol(to_std_string(c), vol)

    def getCompCount(self, str c, str s):
        """
        Returns the number of molecules of a species with identifier string spec 
        in compartment with identifier string comp.

        In a mesh-based simulation this is the combined count from 
        all tetrahedral elements in the compartment.

        Syntax::
            
            getCompCount(comp, spec)
            
        Arguments:
        string comp
        string spec

        Return:
        float

        """
        return self.ptr().getCompCount(to_std_string(c), to_std_string(s))

    def setCompCount(self, str c, str s, double n):
        """
        Set the number of molecules of a species with identifier string spec 
        in compartment with identifier string comp.

        In a mesh-based simulation this is the combined count from 
        all tetrahedral elements in the compartment.

        Syntax::
            
            setCompCount(comp, spec, nspec)
            
        Arguments:
        string comp
        string spec
        int nspec

        Return:
        None

        """
        self.ptr().setCompCount(to_std_string(c), to_std_string(s), n)

    def getCompAmount(self, str c, str s):
        """
        Returns the amount (in mols) of species with identifier string spec in compartment 
        with identifier string comp.

        In a mesh-based simulation this is the combined amount from all 
        tetrahedral elements in the compartment.

        Syntax::
            
            getCompAmount(comp, spec)
            
        Arguments:
        string comp
        string spec

        Return:
        float

        """
        return self.ptr().getCompAmount(to_std_string(c), to_std_string(s))

    def setCompAmount(self, str c, str s, double a):
        """
        Set the amount (in mols) of species with identifier string spec in compartment 
        with identifier string comp.

        In a mesh-based simulation this is the combined amount from all 
        tetrahedral elements in the compartment.

        Syntax::
            
            setCompAmount(comp, spec, amount)
            
        Arguments:
        string comp
        string spec
        float amount

        Return:
        None

        """
        self.ptr().setCompAmount(to_std_string(c), to_std_string(s), a)

    def getCompConc(self, str c, str s):
        """
        Returns the concentration (in Molar units) of species with identifier string spec 
        in compartment with identifier string comp.

        Note: in a mesh-based simulation this is calculated from the combined 
        number of molecules from all tetrahedral elements in the compartment and the total 
        volume of the tetrahedrons.

        Syntax::
            
            getCompConc(comp, spec)
            
        Arguments:
        string comp
        string spec

        Return:
        float

        """
        return self.ptr().getCompConc(to_std_string(c), to_std_string(s))

    def setCompConc(self, str c, str s, double conc):
        """
        Sets the concentration (in Molar units) of species with identifier string spec 
        in compartment with identifier string comp to conc. In a discrete solver the 
        continuous concentration is converted to a discrete number of 
        molecules.

        Note: in a mesh-based simulation the molecules are divided as 
        equally as possible over all tetrahedral elements in the compartment (i.e. a 
        uniform distribution).

        Syntax::

            setCompConc(comp, spec, conc)
            
        Arguments:
        string comp
        string spec
        float conc

        Return:
        None

        """
        self.ptr().setCompConc(to_std_string(c), to_std_string(s), conc)

    def getCompClamped(self, str c, str s):
        """
        Returns True if species with identifier string spec in compartment with identifier 
        string comp is clamped, which means the concentration remains the same 
        regardless of reactions that consume or produce molecules of this species. 
        Returns False if not.

        Note: in a mesh-based simulation it returns True only if the species 
        is clamped in all tetrahedral elements of the compartment.

        Syntax::
            
            getCompClamped(comp, spec)
            
        Arguments:
        string comp
        string spec

        Return:
        bool

        """
        return self.ptr().getCompClamped(to_std_string(c), to_std_string(s))

    def setCompClamped(self, str c, str s, bool b):
        """
        Sets whether the concentration of species with identifier string spec in compartment 
        with identifier string comp is clamped (clamped = True) or not (clamped = False). 
        If a species is clamped the concentration stays the same regardless of reactions 
        that consume or produce molecules of the species.

        Note: in a mesh-based simulation this will set the species to be 
        clamped or not in all tetrahedral elements of the compartment.

        Syntax::
            
            setCompClamped(comp, spec, clamped)
            
        Arguments:
        string comp
        string spec
        bool clamped

        Return:
        None

        """
        self.ptr().setCompClamped(to_std_string(c), to_std_string(s), b)

    def getCompReacK(self, str c, str r):
        """
        Returns the macroscopic reaction constant of reaction with identifier string reac 
        in compartment with identifier string comp. The unit of the reaction constant depends 
        on the order of the reaction.

        Note: In a mesh-based simulation the value for the compartment is 
        returned, although individual tetrahedral elements may have different values 
        (set with setTetReacK).

        Syntax::
            
            getCompReacK(comp, reac)
            
        Arguments:
        string comp
        string reac

        Return:
        float

        """
        return self.ptr().getCompReacK(to_std_string(c), to_std_string(r))

    def setCompReacK(self, str c, str r, double kf):
        """
        Sets the macroscopic reaction constant of reaction with identifier string reac 
        in compartment with identifier string comp to kf. The unit of the reaction constant 
        depends on the order of the reaction.

        Note: In a mesh-based simulation this method sets the reaction 
        constant in all tetrahedral elements of the compartment to kf

        Note: The default value still comes from the steps.model description, so 
        calling reset() will return the reaction constant to that value.

        Syntax::
            
            setCompReacK(comp, reac, kf)
            
        Arguments:
        string comp
        string reac
        float kf

        Return:
        None

        """
        self.ptr().setCompReacK(to_std_string(c), to_std_string(r), kf)

    def getCompReacActive(self, str c, str r):
        """
        Returns whether a reaction with identifier string reac in compartment with identifier 
        string comp is active (True) or not (False). If it's not active this means that a 
        reaction will never occur regardless of whether the reactants are present in 
        sufficient numbers or not. 

        Note: In a mesh-based simulation this method will return True only 
        if the reaction is active in all tetrahedral elements in the compartment. 

        Syntax::
            
            getCompReacActive(comp, reac)
            
        Arguments:
        string comp
        string reac

        Return:
        bool

        """
        return self.ptr().getCompReacActive(to_std_string(c), to_std_string(r))

    def setCompReacActive(self, str c, str r, bool a):
        """
        Activate (active = True) or deactivate (active = False) a reaction with identifier 
        string reac in compartment with identifier string comp. If a reaction is not active 
        this means that a reaction will never occur regardless of whether the reactants are 
        present in sufficient numbers or not.

        Note: In a mesh-based simulation this will activate/deactivate the 
        reaction in all tetrahedral elements in the compartment. 

        Syntax::
            
            setCompReacActive(comp, reac, active)
            
        Arguments:
        string comp
        string reac
        bool active

        Return:
        None

        """
        self.ptr().setCompReacActive(to_std_string(c), to_std_string(r), a)

    def getCompDiffD(self, str c, str d):
        """
        Returns the diffusion constant of diffusion rule with identifier string diff 
        in compartment with identifier string comp. This constant is in units m^2/s.

        Note: In a mesh-based solver the value for the compartment is 
        returned, although individual or groups of tetrahedral elements may have different 
        values (set with setTetDiffD). 

        Syntax::
            
            getCompDiffD(comp, diff)
            
        Arguments:
        string comp
        string diff

        Return:
        float

        """
        return self.ptr().getCompDiffD(to_std_string(c), to_std_string(d))

    def setCompDiffD(self, str c, str d, double dcst):
        """
        Sets the diffusion constant of diffusion rule with identifier string diff 
        in compartment with identifier string comp to dcst (in m^2/s).

        Note: This method will set the diffusion constant in all tetrahedral elements 
        in the compartment.

        Note: The default value still comes from the steps.model description, 
        so calling reset() will return the diffusion constants to that value. 

        Syntax::
            
            setCompDiffD(comp, diff, dcst)
            
        Arguments:
        string comp
        string diff
        float dcst

        Return:
            None

        """
        self.ptr().setCompDiffD(to_std_string(c), to_std_string(d), dcst)

    def getCompDiffActive(self, str c, str d):
        """
        Returns whether a diffusion rule with identifier string diff in compartment with 
        identifier string comp is active (True) or not (False). If diffusion of a species 
        is inactive this means the molecules will remain in place and has the same effect 
        as a diffusion constant of zero. 

        Syntax::
            
            getCompDiffActive(comp, diff)
            
        Arguments:
        string comp
        string diff

        Return:
        bool

        """
        return self.ptr().getCompDiffActive(to_std_string(c), to_std_string(d))

    def setCompDiffActive(self, str c, str d, bool act):
        """
        Activate (active = True) or deactivate (active = False) a diffusion rule with 
        identifier string diff in compartment with identifier string comp. If diffusion 
        of a species is inactive this means the molecules will remain in place and is 
        effectively the same as setting the diffusion constant to zero

        Syntax::
            
            setCompDiffActive(comp, diff, active)
            
        Arguments:
        string comp
        string diff
        bool active

        Return:
        None

        """
        self.ptr().setCompDiffActive(to_std_string(c), to_std_string(d), act)

    def getCompReacC(self, str c, str r):
        """
        Returns the 'stochastic reaction constant' (or 'specific probability rate constant') 
        of reaction with identifier string reac in compartment with identifier string comp.

        The 'stochastic reaction constant' multiplied by infinitesimal time interval dt 
        gives the average probability that one reaction channel of this reaction type 
        will react accordingly in dt.

        Note: in a mesh-based simulation (i.e. Tetexact), the stochastic reaction constant 
        is computed as the weighted mean of the stochastic reaction constants in all 
        tetrahedral elements of the compartment.

        Syntax::
            
            getCompReacC(comp, reac)
            
        Arguments:
        string comp
        string reac

        Return:
        float

        """
        return self.ptr().getCompReacC(to_std_string(c), to_std_string(r))

    def getCompReacH(self, str c, str r):
        """
        Returns h_mu, the distinct number of ways in which reaction with identifier string 
        reac can occur in compartment with identifier string comp, by computing the product 
        of its reactants. Note: in a mesh-based simulation (i.e. Tetexact), returns the sum 
        of the h_mu's over all tetrahedral elements in the compartment. 

        Syntax::
            
            getCompReacH(comp, reac)
            
        Arguments:
        string comp
        string reac

        Return:
        float

        """
        return self.ptr().getCompReacH(to_std_string(c), to_std_string(r))

    def getCompReacA(self, str c, str r):
        """
        Returns the propensity of reaction with identifier string reac in compartment 
        with identifier string comp. 

        The propensity of a reaction is a function of state and is defined as the 
        function whose product with infinitesimal time dt gives the probability 
        that the reaction will occur in the next dt. It is the 'stochastic reaction 
        constant' multiplied by 'h_mu'. 

        Note: in a mesh-based simulation (i.e. Tetexact), the propensity of a reaction 
        in a compartment is computed as the sum of the propensities in all tetrahedral 
        elements of the compartment. 

        Syntax::
            
            getCompReacA(comp, reac)
            
        Arguments:
        string comp
        string reac

        Return:
        float

        """
        return self.ptr().getCompReacA(to_std_string(c), to_std_string(r))

    def getCompReacExtent(self, str c, str r):
        """
        Return the extent of reaction with identifier string reac in compartment with 
        identifier string comp, that is the number of times the reaction has occurred up 
        to the current simulation time. 

        Note: in a mesh-based simulation (i.e. Tetexact), returns the sum of the reaction 
        extents in all tetrahedral elements of the compartment.

        Syntax::
            
            getCompReacExtent(comp, reac)
            
        Arguments:
        string comp
        string reac

        Return:
        index_t

        """
        return self.ptr().getCompReacExtent(to_std_string(c), to_std_string(r))

    def resetCompReacExtent(self, str c, str r):
        """
        Resets the extent of reaction with identifier string reac in compartment with 
        identifier string comp to zero. 

        Note: in a mesh-based simulation (i.e. Tetexact), 
        resets the extents of the reaction in all tetrahedral elements of the compartment.

        Syntax::
            
            resetCompReacExtent(comp, reac)
            
        Arguments:
        string comp
        string reac

        Return:
        None

        """
        self.ptr().resetCompReacExtent(to_std_string(c), to_std_string(r))

    def getPatchArea(self, str p):
        """
        Returns the area of patch with identifier string pat (in m^2).

        Syntax::
            
            getPatchArea(pat)
            
        Arguments:
        string pat

        Return:
        float

        """
        return self.ptr().getPatchArea(to_std_string(p))

    def setPatchArea(self, str p, double area):
        """
        Sets the area of patch with identifier string pat to area a (in m^2).

        Syntax::
            
            setPatchArea(pat, area)
            
        Arguments:
        string pat
        float area

        Return:
        None

        """
        self.ptr().setPatchArea(to_std_string(p), area)

    def getPatchCount(self, str p, str s):
        """
        Returns the number of molecules of species with identifier string spec in patch 
        with identifier string pat.Note: in a mesh-based simulation this 
        is the combined count from all triangular elements in the patch. 

        Syntax::
            
            getPatchCount(pat, spec)
            
        Arguments:
        string pat
        string spec

        Return:
        float

        """
        return self.ptr().getPatchCount(to_std_string(p), to_std_string(s))

    def setPatchCount(self, str p, str s, double n):
        """
        Sets the number of molecules of species with identifier string spec in patch 
        with identifier string pat to n. Note: in a mesh-based simulation the molecules 
        are divided as equally as possible over all triangular elements in 
        the patch (i.e. a uniform distribution). 

        Syntax::
            
            setPatchCount(pat, spec, n)
            
        Arguments:
        string pat
        string spec
        int n

        Return:
        float

        """
        self.ptr().setPatchCount(to_std_string(p), to_std_string(s), n)

    def getPatchAmount(self, str p, str s):
        """
        Returns the amount (in mols) of species with identifier string spec in patch 
        with identifier string pat.

        Note: in a mesh-based simulation this is the combined amount 
        from all triangular elements in the patch. 

        Syntax::
            
            getPatchAmount(pat, spec)
            
        Arguments:
        string pat
        string spec

        Return:
        float

        """
        return self.ptr().getPatchAmount(to_std_string(p), to_std_string(s))

    def setPatchAmount(self, str p, str s, double a):
        """
        Sets the amount (in mols) of species with identifier string spec in patch with 
        identifier string pat to a. In a discrete solver, such as Wmdirect and Tetexact, 
        this continuous value is converted internally into a discrete number of molecules 
        by multiplication with Avogadro's number. 

        Note: in a mesh-based simulation the molecules are divided as 
        equally as possible over all triangular elements in the patch (i.e. a uniform 
        distribution).

        Syntax::
            
            setPatchAmount(pat, spec, a)
            
        Arguments:
        string pat
        string spec
        float a

        Return:
        None

        """
        self.ptr().setPatchAmount(to_std_string(p), to_std_string(s), a)

    def getPatchClamped(self, str p, str s):
        """
        Sets the amount (in mols) of species with identifier string spec in patch with 
        identifier string pat to a. In a discrete solver, such as Wmdirect and Tetexact, 
        this continuous value is converted internally into a discrete number of molecules 
        by multiplication with Avogadro's number. 

        Note: in a mesh-based simulation the molecules are divided as equally 
        as possible over all triangular elements in the patch (i.e. a uniform distribution).

        Syntax::
            
            getPatchClamped(pat, spec)
            
        Arguments:
        string pat
        string spec

        Return:
        bool

        """
        return self.ptr().getPatchClamped(to_std_string(p), to_std_string(s))

    def setPatchClamped(self, str p, str s, bool buf):
        """
        Sets whether the species with identifier string spec in patch with identifier 
        string pat is clamped (clamped = True) or not (clamped = False). If a species 
        is clamped the number of molecules stays the same regardless of surface reactions 
        that consume or produce molecules of the species.

        Note: in a mesh-based simulation this will set the species to be clamped in all 
        triangular elements of the patch.

        Syntax::
            
            setPatchClamped(pat, spec, clamped)
            
        Arguments:
        string pat
        string spec
        bool clamped

        Return:
        None

        """
        self.ptr().setPatchClamped(to_std_string(p), to_std_string(s), buf)

    def getPatchSReacK(self, str p, str r):
        """
        Returns the macroscopic reaction constant of surface reaction with identifier 
        string sreac in patch with identifier string pat. The unit of the reaction constant 
        depends on the order of the reaction.

        Note: In a mesh-based solver the value for the patch is returned, 
        although individual triangle elements may have different values 
        (set with setTriSReacK).

        Syntax::
            
            getPatchSReacK(pat, reac)
            
        Arguments:
        string pat
        string reac

        Return:
        float

        """
        return self.ptr().getPatchSReacK(to_std_string(p), to_std_string(r))

    def setPatchSReacK(self, str p, str r, double kf):
        """
        Sets the macroscopic reaction constant of surface reaction with identifier 
        string sreac in patch with identifier string pat to kf. The unit of the reaction 
        constant depends on the order of the reaction. 

        Note: In a mesh-based simulation this method sets the surface 
        reaction constant in all triangular elements of the patch to kf.

        Note: The default value still comes from the steps.model description, so calling 
        reset() will return the surface reaction constant to that value.

        Syntax::
            
            setPatchSReacK(pat, reac, kf)
            
        Arguments:
        string pat
        string reac
        float kf

        Return:
        None

        """
        self.ptr().setPatchSReacK(to_std_string(p), to_std_string(r), kf)

    def getPatchSReacActive(self, str p, str r):
        """
        Returns whether a surface reaction with identifier string sreac in patch with 
        identifier string pat is active (True) or not (False). If it's not active this means 
        that a surface reaction will never occur regardless of whether the reactants are 
        present in sufficient numbers or not. 

        Note: In a mesh-based simulation this method will return True only 
        if the surface reaction is active in all triangular elements in the patch.

        Syntax::
            
            getPatchSReacActive(pat, reac)
            
        Arguments:
        string pat
        string reac

        Return:
        bool

        """
        return self.ptr().getPatchSReacActive(to_std_string(p), to_std_string(r))

    def setPatchSReacActive(self, str p, str r, bool a):
        """
        Activate (active = True) or deactivate (active = False) a surface reaction with 
        identifier string sreac in patch with identifier string pat. If a surface reaction 
        is not active this means that a reaction will never occur regardless of whether the 
        reactants are present in sufficient numbers or not.

        Note: In a mesh-based simulation this will activate/ deactivate the 
        reaction in all triangular elements in the patch.

        Syntax::
            
            setPatchSReacActive(pat, reac, active)
            
        Arguments:
        string pat
        string reac
        bool active

        Return:
        None

        """
        self.ptr().setPatchSReacActive(to_std_string(p), to_std_string(r), a)

    def getPatchSReacC(self, str p, str r):
        """
        Returns the 'stochastic reaction constant' (or 'specific probability rate constant') 
        of surface reaction with identifier string sreac in patch with identifier string pat.

        Note: in a mesh-based simulation (i.e. Tetexact), the stochastic reaction constant is 
        computed as the weighted mean of the stochastic reaction constants in all triangular 
        elements of the patch.

        Syntax::
            
            getPatchSReacC(pat, reac)
            
        Arguments:
        string pat
        string reac

        Return:
        float

        """
        return self.ptr().getPatchSReacC(to_std_string(p), to_std_string(r))

    def getPatchSReacH(self, str p, str r):
        """
        Returns h_mu, the distinct number of ways in which surface reaction with identifier 
        string sreac can occur in patch with identifier string pat, by computing the product 
        of its reactants. Note: in a mesh-based simulation (i.e. Tetexact), returns the sum 
        of the h_mu's over all triangular elements in the patch. 

        Syntax::
            
            getPatchSReacH(pat, reac)
            
        Arguments:
        string pat
        string reac

        Return:
        float

        """
        return self.ptr().getPatchSReacH(to_std_string(p), to_std_string(r))

    def getPatchSReacA(self, str p, str r):
        """
        Returns the propensity of surface reaction with identifier string sreac in patch 
        with identifier string pat. Note: in a mesh-based simulation (i.e. Tetexact), 
        the propensity of a surface reaction in a patch is computed as the sum of the 
        propensities in all triangular elements of the patch.

        Syntax::
            
            getPatchSReacA(pat, reac)
            
        Arguments:
        string pat
        string reac

        Return:
        float

        """
        return self.ptr().getPatchSReacA(to_std_string(p), to_std_string(r))

    def getPatchSReacExtent(self, str p, str r):
        """
        Returns the extent of surface reaction with identifier string sreac in patch 
        with identifier string pat, that is the number of times the surface reaction 
        has occurred up to the current simulation time. 

        Note: in a mesh-based simulation (i.e. Tetexact), returns the sum of the reaction 
        extents in all triangular elements of the patch.

        Syntax::
            
            getPatchSReacExtent(pat,reac)
            
        Arguments:
        string pat
        string reac

        Return:
        index_t

        """
        return self.ptr().getPatchSReacExtent(to_std_string(p), to_std_string(r))

    def resetPatchSReacExtent(self, str p, str r):
        """
        Resets the extent of reaction with identifier string sreac in patch with identifier 
        string pat to zero. 

        Note: in a mesh-based simulation (i.e. Tetexact), resets the extents of the reaction 
        in all triangular elements of the patch.

        Syntax::
            
            resetPatchSReacExtent(pat, reac)
            
        Arguments:
        string pat
        string reac

        Return:
        None

        """
        self.ptr().resetPatchSReacExtent(to_std_string(p), to_std_string(r))

    def getPatchVDepSReacActive(self, str p, str vsr):
        """
        Returns whether a voltage-dependent surface reaction with identifier string vsreac in patch with 
        identifier string pat is active (True) or not (False). If it's not active this means 
        that the voltage-dependent surface reaction will never occur regardless of whether the reactants are 
        present in sufficient numbers or not. 

        Note: In a mesh-based simulation this method will return True only 
        if the voltage-dependent surface reaction is active in all triangular elements in the patch.

        Syntax::

            getPatchVDepSReacActive(pat, vsreac)

        Arguments:
        string pat
        string vsreac

        Return:
        bool

        """
        return self.ptr().getPatchVDepSReacActive(to_std_string(p), to_std_string(vsr))

    def setPatchVDepSReacActive(self, str p, str vsr, bool a):
        """
        Activate (active = True) or deactivate (active = False) a voltage-dependent surface reaction with 
        identifier string vsreac in patch with identifier string pat. If a voltage-dependent surface reaction 
        is not active this means that a reaction will never occur regardless of whether the 
        reactants are present in sufficient numbers or not.

        Note: In a mesh-based simulation this will activate/ deactivate the 
        reaction in all triangular elements in the patch.

        Syntax::

            setPatchVDepSReacActive(pat, vsreac, active)

        Arguments:
        string pat
        string vsreac
        bool active

        Return:
        None

        """
        self.ptr().setPatchVDepSReacActive(to_std_string(p), to_std_string(vsr), a)

    def getNComps(self, ):
        """
        Return the number of compartments in the solver.

        Syntax::

            getNComps()

        Arguments:
        None

        Return:
        uint

        """
        return self.ptr().getNComps()

    def getNPatches(self, ):
        """
        Return the number of patches in the solver.

        Syntax::

            getNPatches()

        Arguments:
        None

        Return:
        uint

        """
        return self.ptr().getNPatches()

    def getCompName(self, uint c_idx):
        """
        Return the name of a compartment by its index in the solver.
            
        Syntax::
            
            getCompName(c_idx)
            
        Arguments:
        uint c_idx
            
        Return:
        string

        """
        return self.ptr().getCompName(c_idx)

    def getPatchName(self, uint p_idx):
        """
        Return the name of a patch by its index in the solver.

        Syntax::

            getPatchName(p_idx)

        Arguments:
        uint p_idx

        Return:
        string

        """
        return self.ptr().getPatchName(p_idx)

    def getNCompSpecs(self, uint c_idx):
        """
        Get number of species in a compartment.

        Syntax::

            getNCompSpecs(c_idx)

        Arguments:
        uint c_idx

        Return:
        uint

        """
        return self.ptr().getNCompSpecs(c_idx)

    def getNPatchSpecs(self, uint p_idx):
        """
        Get number of species in a patch.

        Syntax::

            getNPatchSpecs(p_idx)

        Arguments:
        uint p_idx

        Return:
        uint

        """
        return self.ptr().getNPatchSpecs(p_idx)

    def getCompSpecName(self, uint c_idx, uint s_idx):
        """
        Get the name of a species in a compartment.

        Syntax::

            getCompSpecName(c_idx, s_idx)

        Arguments:
        uint c_idx
        uint s_idx

        Return:
        string

        """
        return self.ptr().getCompSpecName(c_idx, s_idx)

    def getPatchSpecName(self, uint p_idx, uint s_idx):
        """
        Get the name of a species in a patch.

        Syntax::

            getPatchSpecName(p_idx, s_idx)

        Arguments:
        uint p_idx
        uint s_idx

        Return:
        string

        """
        return self.ptr().getPatchSpecName(p_idx, s_idx)


    @staticmethod
    cdef _py_API from_ptr(API *ptr):
        cdef _py_API obj = _py_API.__new__(_py_API )
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_API from_ref(const API &ref):
        return _py_API.from_ptr(<API*>&ref)

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TetAPI(_py_API):
    "Python wrapper class for tetrahedron solvers API"
# ----------------------------------------------------------------------------------------------------------------------

    #Constants
    EF_NONE      = steps_solver.EF_NONE
    EF_DEFAULT   = steps_solver.EF_DEFAULT
    EF_DV_BDSYS  = steps_solver.EF_DV_BDSYS
    EF_DV_PETSC  = steps_solver.EF_DV_PETSC

    # ---- VIRTUAL - doesnt call original constructor ------
    def __init__(self, _py_Model m, _py_Geom g, _py_RNG r):
        _py_API.__init__(m, g, r)

    def getTetVol(self, index_t tidx):
        """
        Returns the volume (in m^3) of the tetrahedral element with index idx.

        Syntax::
            
            getTetVol(idx)
            
        Arguments:
        index_t idx

        Return:
        float

        """
        return self.ptr().getTetVol(tidx)

    def getTetSpecDefined(self, index_t tidx, str s):
        """
        Returns whether species with identifier string spec is defined
        in the tetrahedral element with index idx.

        Syntax::
            
            getTetSpecDefined(idx, spec)
            
        Arguments:
        index_t idx
        string spec

        Return:
        bool

        """
        return self.ptr().getTetSpecDefined(tidx, to_std_string(s))

    def getTetCount(self, index_t tidx, str s):
        """
        Returns the number of molecules of species with identifier string spec 
        in the tetrahedral element with index idx.

        Syntax::
            
            getTetCount(idx, spec)
            
        Arguments:
        index_t idx
        string spec

        Return:
        int

        """
        return self.ptr().getTetCount(tidx, to_std_string(s))

    def setTetCount(self, index_t tidx, str s, double n):
        """
        Sets the number of molecules of species with identifier string spec in 
        tetrahedral element with index idx to n.

        Syntax::
            
            setTetCount(idx, spec, n)
            
        Arguments:
        index_t idx
        string spec
        int n

        Return:
        None

        """
        self.ptr().setTetCount(tidx, to_std_string(s), n)

    def getTetAmount(self, index_t tidx, str s):
        """
        Returns the amount (in mols) of species with identifier string spec in 
        tetrahedral element with index idx.

        Syntax::
            
            getTetAmount(idx, spec)
            
        Arguments:
        index_t idx
        string spec

        Return:
        float

        """
        return self.ptr().getTetAmount(tidx, to_std_string(s))

    def setTetAmount(self, index_t tidx, str s, double m):
        """
        Sets the amount (in mols) of species with identifier string spec in tetrahedral 
        element with index idx to a. This continuous value must be converted internally 
        to a discrete number of molecules by multiplication with Avogadro's 
        number. 

        Due to the small volumes of tetrahedral elements the difference 
        between 'rounding up' and 'rounding down' can be a significant difference in 
        concentration.

        Syntax::
            
            setTetAmount(idx, spec, a)
            
        Arguments:
        index_t idx
        string spec
        float a

        Return:
        None

        """
        self.ptr().setTetAmount(tidx, to_std_string(s), m)

    def getTetConc(self, index_t tidx, str s):
        """
        Returns the concentration (in Molar units) of species with identifier 
        string spec in a tetrahedral element with index idx.

        Syntax::
            
            getTetConc(idx, spec)
            
        Arguments:
        index_t idx
        string spec

        Return:
        float

        """
        return self.ptr().getTetConc(tidx, to_std_string(s))

    def setTetConc(self, index_t tidx, str s, double c):
        """
        Sets the concentration (in Molar units) of species with identifier string spec 
        in a tetrahedral element with index idx to conc.This continuous value must be 
        converted internally to a discrete number of molecules. 

        Due to the small volumes of tetrahedral elements the difference between 'rounding 
        up' and 'rounding down' can be a large difference in concentration.

        Syntax::
            
            setTetConc(idx, spec, conc)
            
        Arguments:
        index_t idx
        string spec
        float conc

        Return:
        None

        """
        self.ptr().setTetConc(tidx, to_std_string(s), c)

    def getTetClamped(self, index_t tidx, str s):
        """
        Returns True if concentration of species with identifier string spec in tetrahedral 
        element with index idx is clamped, which means the concentration stays the 
        same regardless of reactions that consume or produce molecules of this species or 
        diffusion of this species into or out of the tetrahedral element. Returns False if 
        not.

        Syntax::
            
            getTetClamped(idx, spec)
            
        Arguments:
        index_t idx
        string spec

        Return:
        bool

        """
        return self.ptr().getTetClamped(tidx, to_std_string(s))

    def setTetClamped(self, index_t tidx, str s, bool buf):
        """
        Sets whether the concentration of species spec in tetrahedral element with 
        index idx is clamped (clamped = True) or not (clamped = False). 
        If a species is clamped the concentration stays the same regardless 
        of reactions that consume or produce molecules of the species or 
        diffusion of the species into or out of the tetrahedral element.

        Syntax::
            
            setTetClamped(idx, spec, clamped)
            
        Arguments:
        index_t idx
        string spec
        bool clamped

        Return:
        None

        """
        self.ptr().setTetClamped(tidx, to_std_string(s), buf)

    def getTetReacK(self, index_t tidx, str r):
        """
        Returns the macroscopic reaction constant of reaction with identifier string reac 
        in tetrahedral element with index idx. The unit of the reaction constant depends 
        on the order of the reaction.

        Syntax::
            
            getTetReacK(idx, reac)
            
        Arguments:
        index_t idx
        string reac

        Return:
        float

        """
        return self.ptr().getTetReacK(tidx, to_std_string(r))

    def setTetReacK(self, index_t tidx, str r, double kf):
        """
        Sets the macroscopic reaction constant of reaction with identifier string reac 
        in tetrahedral element with index idx to kf. The units of the reaction constant 
        depends on the order of the reaction.

        Syntax::
            
            setTetReacK(idx, reac, kf)
            
        Arguments:
        index_t idx
        string reac
        float kf

        Return:
        None

        """
        self.ptr().setTetReacK(tidx, to_std_string(r), kf)

    def getTetReacActive(self, index_t tidx, str r):
        """
        Returns whether reaction with identifier string reac in tetrahedral element 
        with index idx is active (True) or not (False). If it's not active this means 
        that the reaction will never occur regardless of whether reactants are present 
        in sufficient numbers or not.

        Syntax::
            
            getTetReacActive(idx, reac)
            
        Arguments:
        index_t idx
        string reac

        Return:
        bool

        """
        return self.ptr().getTetReacActive(tidx, to_std_string(r))

    def setTetReacActive(self, index_t tidx, str r, bool act):
        """
        Activate (active = True) or deactivate (active = False) a reaction with identifier 
        string reac in tetrahedral element with index idx. If it's not active this means 
        that the reaction will never occur regardless of whether reactants are present 
        in sufficient numbers or not.

        Syntax::
            
            setTetReacActive(idx, reac, active)
            
        Arguments:
        index_t idx
        string reac
        bool active

        Return:
        None

        """
        self.ptr().setTetReacActive(tidx, to_std_string(r), act)

    def getTetDiffD(self, index_t tidx, str d, index_t direction_tet=UNKNOWN_TET):
        """
        Returns the diffusion constant of diffusion rule with identifier string diff 
        in tetrahedral element with index idx. This constant is in units m^2/s. If direction_tet is specified, return the diffusion constant towards that direction.

        Syntax::
            
            getTetDiffD(idx, diff, direction_tet = UNKNOWN_TET)
            
        Arguments:
        index_t idx
        string diff
        direction_tet
            
        Return:
        float

        """
        return self.ptr().getTetDiffD(tidx, to_std_string(d), direction_tet)

    def setTetDiffD(self, index_t tidx, str d, double dk, index_t direction_tet=UNKNOWN_TET):
        """
        Sets the diffusion constant of diffusion rule with identifier string diff in 
        tetrahedral element with index idx to dcst (in m^2/s). Specify direction_tet to set the constant only towards a given tetrahedron direction.
        Syntax::
            
            setTetDiffD(idx, diff, dcst, direction_tet = UNKNOWN_TET)
            
        Arguments:
        index_t idx
        string diff
        float dcst
        index_t direction_tet

        Return:
        None

        """
        self.ptr().setTetDiffD(tidx, to_std_string(d), dk, direction_tet)

    def getTetDiffActive(self, index_t tidx, str d):
        """
        Returns whether diffusion with identifier string diff in tetrahedral element 
        with index idx is active (True) or not (False). If diffusion of a species 
        is inactive this means the molecules will never diffuse out of the tetrahedron 
        and has the same effect as a diffusion constant of zero.

        Syntax::
            
            getTetDiffActive(idx, diff)
            
        Arguments:
        index_t idx
        string diff

        Return:
        bool

        """
        return self.ptr().getTetDiffActive(tidx, to_std_string(d))

    def setTetDiffActive(self, index_t tidx, str d, bool act):
        """
        Activate (active = True) or deactivate (active = False) diffusion rule with 
        identifier string diff in tetrahedral element with index idx. If diffusion of 
        a species is inactive this means the molecules will never diffuse out of the 
        tetrahedron and has the same effect as a diffusion constant of zero. 

        Syntax::
            
            setTetDiffActive(idx, diff, active)
            
        Arguments:
        index_t idx
        string diff
        bool active

        Return:
        None

        """
        self.ptr().setTetDiffActive(tidx, to_std_string(d), act)

    def getTetReacC(self, index_t tidx, str r):
        """
        Returns the 'stochastic reaction constant' (or 'specific probability rate constant') 
        of reaction with identifier string reac in tetrahedral element with index idx.

        Syntax::
            
            getTetReacC(idx, reac)
            
        Arguments:
        index_t idx
        string reac

        Return:
        float

        """
        return self.ptr().getTetReacC(tidx, to_std_string(r))

    def getTetReacH(self, index_t tidx, str r):
        """
        Returns h_mu, the distinct number of ways in which reaction with identifier string 
        reac can occur in tetrahedral element with index idx, by computing the product of 
        its reactants.

        Syntax::
            
            getTetReacH(idx, reac)
            
        Arguments:
        index_t idx
        string reac

        Return:
        float

        """
        return self.ptr().getTetReacH(tidx, to_std_string(r))

    def getTetReacA(self, index_t tidx, str r):
        """
        Returns the propensity of reaction with identifier string reac in tetrahedral 
        element with index idx.

        Syntax::
            
            getTetReacA(idx, reac)
            
        Arguments:
        index_t idx
        string reac

        Return:
        float

        """
        return self.ptr().getTetReacA(tidx, to_std_string(r))

    def getTetDiffA(self, index_t tidx, str d):
        """
        Returns the propensity of diffusion rule with identifier string diff in 
        tetrahedral element with index idx. 

        Syntax::
            
            getTetDiffA(idx, reac)
            
        Arguments:
        index_t idx
        string reac

        Return:
        float

        """
        return self.ptr().getTetDiffA(tidx, to_std_string(d))

    def getTetV(self, index_t tidx):
        """
        Returns the potential (in volts) of tetrahedral element with index idx, taken at the barycenter.
        			
        Syntax::
        			 
        	getTetV(idx)
        			 
        Arguments:
        index_t idx
        		
        Return:
        float

        """
        return self.ptr().getTetV(tidx)

    def setTetV(self, index_t tidx, double v):
        """
        Set the potential (in volts) of tetrahedral element with index idx.
        			
        Syntax::
        			 
        	setTetV(idx, v)
        			 
        Arguments:
        index_t idx
        float v
        			 
        Return:
        None

        """
        self.ptr().setTetV(tidx, v)

    def getTetVClamped(self, index_t tidx):
        """
        Returns true if the potential of tetrahedral element with index idx is clamped
        to some voltage.
        			
        Syntax::
        			 
        	getTetVClamped(idx)
        			 
        Arguments:
        index_t idx
        		
        Return:
        bool

        """
        return self.ptr().getTetVClamped(tidx)

    def setTetVClamped(self, index_t tidx, bool cl):
        """
        Sets whether the potential of tetrahedral element with index idx is clamped
        (clamped = True) or not (clamped = False).
        			
        Syntax::
        			 
        	setTetVClamped(idx, clamped)
        			 
        Arguments:
        index_t idx
        bool clamped
        			 
        Return:
        None

        """
        self.ptr().setTetVClamped(tidx, cl)

    def setDiffBoundaryDiffusionActive(self, str db, str s, bool act):
        """
        Activates or inactivates diffusion across a diffusion boundary for a species.
                     
        Syntax::
                     
            setDiffBoundaryDiffusionActive(diffb, spec, act)
                     
        Arguments:
        string diffb
        string spec
        bool act
                     
        Return:
        None

        """
        self.ptr().setDiffBoundaryDiffusionActive(to_std_string(db), to_std_string(s), act)

    def getDiffBoundaryDiffusionActive(self, str db, str s):
        """
        Returns whether diffusion is active across a diffusion boundary for a species.
                     
        Syntax::
                     
            getDiffBoundaryDiffusionActive(diffb, spec)
                     
        Arguments:
        string diffb
        string spec
                     
        Return:
        bool

        """
        return self.ptr().getDiffBoundaryDiffusionActive(to_std_string(db), to_std_string(s))

    def setDiffBoundaryDcst(self, str db, str s, double dcst, str direction_comp=""):
        """
        Set the diffusion constant of tetrahedrons across a diffusion boundary. If direction_comp is provided, only set dcsts of diffusion towards it (Directional dcsts of diffusions in tetrahedrons in the other compartment of the diffusion boundary towards tetrahedons in the direction compartment).
                     
        Syntax::
                     
            setDiffBoundaryDcst(diffb, spec, dcst, direction_comp = '')
                     
        Arguments:
        string diffb
        string spec
        float dcst
        string direction_comp
                     
        Return:
        None

        """
        self.ptr().setDiffBoundaryDcst(to_std_string(db), to_std_string(s), dcst, to_std_string(direction_comp))

    def setSDiffBoundaryDiffusionActive(self, str sdb, str s, bool act):
        """
        Activates or inactivates diffusion across a surface diffusion boundary for a species.
                     
        Syntax::
                     
            setSDiffBoundaryDiffusionActive(sdiffb, spec, act)
                     
        Arguments:
        string sdiffb
        string spec
        bool act
                     
        Return:
        None

        """
        self.ptr().setSDiffBoundaryDiffusionActive(to_std_string(sdb), to_std_string(s), act)

    def getSDiffBoundaryDiffusionActive(self, str sdb, str s):
        """
        Returns whether diffusion is active across a surface diffusion boundary for a species.
                     
        Syntax::
                     
            getSDiffBoundaryDiffusionActive(sdiffb, spec)
                     
        Arguments:
        string sdiffb
        string spec
                     
        Return:
        bool

        """    
        return self.ptr().getSDiffBoundaryDiffusionActive(to_std_string(sdb), to_std_string(s))

    def setSDiffBoundaryDcst(self, str sdb, str s, double dcst, str direction_patch=""):
        """
        Set the diffusion constant of triangles across a surface diffusion boundary. 
        If direction_patch is provided, only set dcsts of diffusion towards it 
        (Directional dcsts of diffusions in triangles in the other patches of the diffusion boundary 
        towards triangles in the direction patch).
                     
        Syntax::
                     
            setSDiffBoundaryDcst(sdiffb, spec, dcst, direction_patch = '')
                     
        Arguments:
        string sdiffb
        string spec
        float dcst
        string direction_patch
                     
        Return:
        None

        """
        self.ptr().setSDiffBoundaryDcst(to_std_string(sdb), to_std_string(s), dcst, to_std_string(direction_patch))

    def getTriArea(self, index_t tidx):
        """
        Returns the area (in m^2) of the triangular element with index idx.

        Syntax::
            
            getTriArea(idx)
            
        Arguments:
        index_t idx

        Return:
        float

        """
        return self.ptr().getTriArea(tidx)

    def getTriSpecDefined(self, index_t tidx, str s):
        """
        Returns whether species with identifier string spec is defined
        in the triangle element with index idx.

        Syntax::
            
            getTriSpecDefined(idx, spec)
            
        Arguments:
        index_t idx
        string spec

        Return:
        bool

        """
        return self.ptr().getTriSpecDefined(tidx, to_std_string(s))

    def getTriCount(self, index_t tidx, str s):
        """
        Returns the number of molecules of species with identifier string spec 
        in the triangular element with index idx.

        Syntax::
            
            getTriCount(idx, spec)
            
        Arguments:
        index_t idx
        string spec

        Return:
        float

        """
        return self.ptr().getTriCount(tidx, to_std_string(s))

    def setTriCount(self, index_t tidx, str s, double n):
        """
        Sets the number of molecules of species with identifier string spec in 
        triangular element with index idx to n. 

        Syntax::
            
            setTriCount(idx, spec, n)
            
        Arguments:
        index_t idx
        string spec
        int n

        Return:
        None

        """
        self.ptr().setTriCount(tidx, to_std_string(s), n)

    def getTriAmount(self, index_t tidx, str s):
        """
        Returns the amount (in mols) of species with identifier string spec in triangular 
        element with index idx.  

        Syntax::
            
            getTriAmount(idx, spec)
            
        Arguments:
        index_t idx
        string spec

        Return:
        float

        """
        return self.ptr().getTriAmount(tidx, to_std_string(s))

    def setTriAmount(self, index_t tidx, str s, double m):
        """
        Sets the amount (in mols) of species with identifier string spec in triangular 
        element with index idx to a. This continuous value must be converted internally 
        to a discrete number of molecules by multiplication with Avogadro's number. 

        Syntax::
            
            setTriAmount(idx, spec, a)
            
        Arguments:
        index_t idx
        string spec
        float a

        Return:
        None

        """
        self.ptr().setTriAmount(tidx, to_std_string(s), m)

    def getTriClamped(self, index_t tidx, str s):
        """
        Returns True if the species with identifier string spec in triangular element 
        with index idx is clamped, which means the number of molecules stays 
        the same regardless of reactions that consume or produce molecules of this species. 
        Returns False if not.

        Syntax::
            
            getTriClamped(idx, spec)
            
        Arguments:
        index_t idx
        string spec

        Return:
        bool

        """
        return self.ptr().getTriClamped(tidx, to_std_string(s))

    def setTriClamped(self, index_t tidx, str s, bool buf):
        """
        Sets whether the concentration of species spec in triangular element with index idx 
        is clamped (clamped = True) or not (clamped = False). If a species is clamped the 
        concentration stays the same regardless of reactions that consume or produce 
        molecules of the species. 

        Syntax::
            
            setTriClamped(idx, spec, clamped)
            
        Arguments:
        index_t idx
        string spec
        bool clamped

        Return:
        None

        """
        self.ptr().setTriClamped(tidx, to_std_string(s), buf)

    def getTriSReacK(self, index_t tidx, str r):
        """
        Returns the macroscopic reaction constant of surface reaction with identifier 
        string sreac in triangular element with index idx. The units of the reaction 
        constant depends on the order of the reaction. 

        Syntax::
            
            getTriSReacK(idx, reac)
            
        Arguments:
        index_t idx
        string reac

        Return:
        float

        """
        return self.ptr().getTriSReacK(tidx, to_std_string(r))

    def setTriSReacK(self, index_t tidx, str r, double kf):
        """
        Sets the macroscopic reaction constant of surface reaction with identifier 
        string sreac in triangular element with index idx to kf. The units of the 
        reaction constant depends on the order of the reaction.

        Syntax::
            
            setTriSReacK(idx, reac, kf)
            
        Arguments:
        index_t idx
        string reac
        float kf

        Return:
        None

        """
        self.ptr().setTriSReacK(tidx, to_std_string(r), kf)

    def getTriSReacActive(self, index_t tidx, str r):
        """
        Returns whether surface reaction with identifier string sreac in triangular 
        element with index idx is active (True) or not (False). If it's not active 
        this means that the surface reaction will never occur regardless of whether 
        reactants are present in sufficient numbers or not. 

        Syntax::
            
            getTriSReacActive(idx, reac)
            
        Arguments:
        index_t idx
        string reac

        Return:
        bool

        """
        return self.ptr().getTriSReacActive(tidx, to_std_string(r))

    def setTriSReacActive(self, index_t tidx, str r, bool act):
        """
        Activate (active = True) or deactivate (active = False) a surface reaction 
        with identifier string sreac in triangular element with index idx. If it's 
        not active this means that the surface reaction will never occur regardless 
        of whether reactants are present in sufficient numbers or not.  

        Syntax::
            
            setTriSReacActive(idx, reac, active)
            
        Arguments:
        index_t idx
        string reac
        bool active

        Return:
        None

        """
        self.ptr().setTriSReacActive(tidx, to_std_string(r), act)

    def getTriSReacC(self, index_t tidx, str r):
        """
        Returns the 'stochastic reaction constant' (or 'specific probability rate constant') 
        of surface reaction with identifier string sreac in triangular element with index idx.  

        Syntax::
            
            getTriSReacC(idx, reac)
            
        Arguments:
        index_t idx
        string reac

        Return:
        float

        """
        return self.ptr().getTriSReacC(tidx, to_std_string(r))

    def getTriSReacH(self, index_t tidx, str r):
        """
        Returns h_mu, the distinct number of ways in which surface reaction with identifier 
        string sreac can occur in triangular element with index idx, by computing the product 
        of its reactants. 

        Syntax::
            
            getTriSReacH(idx, reac)
            
        Arguments:
        index_t idx
        string reac

        Return:
        float

        """
        return self.ptr().getTriSReacH(tidx, to_std_string(r))

    def getTriSReacA(self, index_t tidx, str r):
        """
        Returns the propensity of surface reaction with identifier string sreac 
        in triangular element with index idx. 

        Syntax::
            
            getTriSReacA(idx, reac)
            
        Arguments:
        index_t idx
        string reac

        Return:
        float

        """
        return self.ptr().getTriSReacA(tidx, to_std_string(r))

    def getTriSDiffD(self, index_t tidx, str d, index_t direction_tri=UNKNOWN_TRI):
        """
        Returns the diffusion constant of diffusion rule with identifier string diff 
        in triangle element with index idx. If direction_tri is specified, return the diffusion constant towards that direction.

        Syntax::
            
            getTriDiffD(idx, diff, direction_tri = UNKNOWN_TRI)
            
        Arguments:
        index_t tidx
        string diff
        index_t direction_tri
            
        Return:
        float

        """
        return self.ptr().getTriSDiffD(tidx, to_std_string(d), direction_tri)

    def setTriSDiffD(self, index_t tidx, str d, double dk, index_t direction_tri=UNKNOWN_TRI):
        """
        Sets the diffusion constant of diffusion rule with identifier string diff in 
        triangle element with index idx to dcst. Specify direction_tri to set the constant only towards a given triangle direction.
        Syntax::
            
            setTriSDiffD(idx, diff, dcst, direction_tri = UNKNOWN_TRI)
            
        Arguments:
        index_t idx
        string diff
        float dcst
        index_t direction_tri

        Return:
        None

        """
        self.ptr().setTriSDiffD(tidx, to_std_string(d), dk, direction_tri)

    def getTriV(self, index_t tidx):
        """
        Returns the potential (in volts) of triangle element with index idx, taken at the barycenter.
        			
        Syntax::
        			 
        	getTriV(idx)
        			 
        Arguments:
        index_t idx
        		
        Return:
        float

        """
        return self.ptr().getTriV(tidx)

    def setTriV(self, index_t tidx, double v):
        """
        Set the potential (in volts) of triangle element with index idx.
        			
        Syntax::
        			 
        	setTriV(idx, v)
        			 
        Arguments:
        index_t idx
        float v
        			 
        Return:
        None

        """
        self.ptr().setTriV(tidx, v)

    def getTriVClamped(self, index_t tidx):
        """
        Returns true if the potential of triangle element with index idx is clamped
        to some voltage.
        			
        Syntax::
        			 
        	getTriVClamped(idx)
        			 
        Arguments:
        index_t idx
        		
        Return:
        bool

        """
        return self.ptr().getTriVClamped(tidx)

    def setTriVClamped(self, index_t tidx, bool cl):
        """
        Sets whether the potential of triangle element with index idx is clamped
        (clamped = True) or not (clamped = False).
        			
        Syntax::
        			 
        	setTriVClamped(idx, clamped)
        			 
        Arguments:
        index_t idx
        bool clamped
        			 
        Return:
        None

        """
        self.ptr().setTriVClamped(tidx, cl)

    def getTriOhmicI(self, index_t tidx, str oc=''):
        """
        Returns the ohmic current of triangle element with index idx, in amps.
        			
        Syntax::
        			 
        	getTriOhmicI(idx, oc)
        			 
        Arguments:
        index_t idx
        string oc (default = '') 
        		
        Return:
        float

        """
        if oc == '':
            return self.ptr().getTriOhmicI(tidx)
        return self.ptr().getTriOhmicI(tidx, to_std_string(oc))

    def getTriGHKI(self, index_t tidx, str ghk=''):
        """
        Returns the GHK current of triangle element with index idx, in amps.
                     
        Syntax::
        			 
            getTriGHKI(idx)
        			 
        Arguments:
        index_t idx
        string ghk (default = '') 
                     
        Return:
        float

        """
        if ghk == '':
            return self.ptr().getTriGHKI(tidx)
        return self.ptr().getTriGHKI(tidx, to_std_string(ghk))

    def getTriI(self, index_t tidx):
        """
        Returns the current of triangle element with index idx, in amps, 
        at the last EField calculation step.
        			 
        Syntax::
        			 
            getTriI(idx)
        			 
        Arguments:
        index_t idx
        			 
        Return:
        float

        """
        return self.ptr().getTriI(tidx)

    def getTriIClamp(self, index_t tidx):
        """
        Get current clamp to triangle element with index tidx.
        NOTE: Convention is maintained that a positive current clamp is depolarizing, a negative current clamp is hyperpolarizing.
        			
        Syntax::
        			 
        	getTriIClamp(tidx)
        			 
        Arguments:
        index_t tidx
        			 
        Return:
        float

        """
        return self.ptr().getTriIClamp(tidx)

    def setTriIClamp(self, index_t tidx, double i):
        """
        Set current clamp to triangle element with index idx to current i (amps).
        NOTE: Convention is maintained that a positive current clamp is depolarizing, a negative current clamp is hyperpolarizing.
        			
        Syntax::
        			 
        	setTriIClamp(idx, i)
        			 
        Arguments:
        index_t idx
        float i
        			 
        Return:
        None

        """
        self.ptr().setTriIClamp(tidx, i)

    def getTriVDepSReacActive(self, index_t tidx, str vsr):
        """
        Returns whether voltage-dependent surface reaction with identifier string vsreac in triangular 
        element with index idx is active (True) or not (False). If it's not active 
        this means that the voltage-dependent surface reaction will never occur regardless of whether 
        reactants are present in sufficient numbers or not. 

        Syntax::

            getTriVDepSReacActive(idx, reac)

        Arguments:
        index_t idx
        string vsreac

        Return:
        bool

        """
        return self.ptr().getTriVDepSReacActive(tidx, to_std_string(vsr))

    def setTriVDepSReacActive(self, index_t tidx, str vsr, bool act):
        """
        Activate (active = True) or deactivate (active = False) a voltage-dependent surface reaction 
        with identifier string vsreac in triangular element with index idx. If it's 
        not active this means that the voltage-dependent surface reaction will never occur regardless 
        of whether reactants are present in sufficient numbers or not.  

        Syntax::

            setTriVDepSReacActive(idx, vsreac, active)

        Arguments:
        index_t idx
        string vsreac
        bool active

        Return:
        None

        """
        self.ptr().setTriVDepSReacActive(tidx, to_std_string(vsr), act)

    def setTriCapac(self, index_t tidx, double cm):
        """
        Sets the specific membrane capacitance (in farad / m^2) of tri with index tidx.
        			 
        Syntax::
        			 
            setTriCapac(tidx, cm)
        			 
        Arguments:
        index_t tidx
        float cm
        			 
        Return:
        None

        """    
        self.ptr().setTriCapac(tidx, cm)

    def getVertV(self, index_t vidx):
        """
        Returns the potential (in volts) of vertex element with index idx.
        			
        Syntax::
        			 
        	getVertV(idx)
        			 
        Arguments:
        index_t idx
        		
        Return:
        float

        """
        return self.ptr().getVertV(vidx)

    def setVertV(self, index_t vidx, double v):
        """
        Set the potential (in volts) of vertex element with index idx.
        			
        Syntax::
        			 
        	setVertV(idx, v)
        			 
        Arguments:
        index_t idx
        float v
        			 
        Return:
        None

        """
        self.ptr().setVertV(vidx, v)

    def getVertVClamped(self, index_t vidx):
        """
        Returns true if the potential of vertex element with index idx is clamped
        to some voltage.
        			
        Syntax::
        			 
        	getVertVClamped(idx)
        			 
        Arguments:
        index_t idx
        		
        Return:
        bool

        """
        return self.ptr().getVertVClamped(vidx)

    def setVertVClamped(self, index_t vidx, bool cl):
        """
        Sets whether the potential of vertex element with index idx is clamped
        (clamped = True) or not (clamped = False).
        			
        Syntax::
        			 
        	setVertVClamped(idx, clamped)
        			 
        Arguments:
        index_t idx
        bool clamped
        			 
        Return:
        None

        """
        self.ptr().setVertVClamped(vidx, cl)

    def getVertIClamp(self, index_t vidx):
        """
        Get current clamp to vertex element with index vidx (Amps). 
        NOTE: Convention is maintained that a positive current clamp is depolarizing, a negative current clamp is hyperpolarizing.

        Syntax::

            getVertIClamp(vidx)
        			 
        Arguments:
        index_t vidx
        			 
        Return:
        float

        """
        return self.ptr().getVertIClamp(vidx)

    def setVertIClamp(self, index_t vidx, double i):
        """
        Set current clamp to vertex element with index idx to current i (Amps). 
        NOTE: Convention is maintained that a positive current clamp is depolarizing, a negative current clamp is hyperpolarizing.

        Syntax::

            setVertIClamp(idx, i)
        			 
        Arguments:
        index_t idx
        float i
        			 
        Return:
        None

        """
        self.ptr().setVertIClamp(vidx, i)

    def setMembPotential(self, str m, double v):
        """
        Sets the potential (in volts) of membrane with string identifier memb.
        NOTE: This method will set the potential of all nodes in the volume conductor
        to the same value.

        Syntax::
        			
        	setMembPotential(memb, v)
        		
        Arguments:
        string memb
        float v

        Return:
        None

        """
        self.ptr().setMembPotential(to_std_string(m), v)

    def setMembCapac(self, str m, double cm):
        """
        Sets the specific membrane capacitance (in farad / m^2) of membrane with string identifier memb.
        			 
        Syntax::
        			 
            setMembCapac(memb, cm)
        			 
        Arguments:
        string memb
        float cm
        			 
        Return:
        None

        """
        self.ptr().setMembCapac(to_std_string(m), cm)

    def setMembVolRes(self, str m, double ro):
        """
        Sets the bulk electrical resistivity (in ohm.m) of the volume conductor 
        assocaited with membrane with string identifier memb.
        			 
        Syntax::
        			 
            setMembVolRes(memb, ro)
        			 
        Arguments:
        string memb
        float ro
        			 
        Return:
        None

        """
        self.ptr().setMembVolRes(to_std_string(m), ro)

    def setMembRes(self, str m, double ro, double vrev):
        """
        Sets the  electrical resistivity (in ohm/m^2) of the membrane with string identifier memb. Reversal potential vrev is required in Volts.
        			 
        Syntax::
        			 
            setMembRes(memb, ro, vrev)
        			 
        Arguments:
        string memb
        float ro
        float vrev
        			 
        Return:
        None

        """
        self.ptr().setMembRes(to_std_string(m), ro, vrev)

    @staticmethod
    cdef _py_TetAPI from_ptr(API *ptr):
        cdef _py_TetAPI obj = _py_TetAPI.__new__(_py_TetAPI )
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_TetAPI from_ref(const API &ref):
        return _py_TetAPI.from_ptr(<API*>&ref)
